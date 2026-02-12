---
name: scientific-systems-biology
description: |
  システム生物学解析スキル。動的モデリング（ODE / SBML）・
  代謝フラックス解析（FBA / pFBA）・遺伝子制御ネットワーク推定（GRN）・
  シグナル伝達経路モデリング・パラメータ推定・感度解析・
  BioModels/Reactome/KEGG/BiGG 統合パイプライン。
---

# Scientific Systems Biology

システム生物学の定量的モデリングパイプラインを提供する。
ODE ベースの動的モデル、フラックスバランス解析（FBA）、
遺伝子制御ネットワーク（GRN）推定、パラメータ推定・感度解析を扱い、
BioModels・Reactome・KEGG・BiGG の統合的活用を支援する。

## When to Use

- 生物学的パスウェイの動的モデリング（ODE）が必要なとき
- 代謝ネットワークのフラックスバランス解析（FBA）を行うとき
- 遺伝子制御ネットワーク（GRN）を推定するとき
- BioModels / SBML モデルの取得・シミュレーションを行うとき
- モデルパラメータの推定・感度解析を行うとき

---

## Quick Start

## 1. SBML モデルシミュレーション

```python
import numpy as np
import pandas as pd

def simulate_sbml_model(sbml_file, duration=100, n_points=1000):
    """
    SBML モデルのシミュレーション。

    SBML (Systems Biology Markup Language):
      生物学的モデルの標準交換フォーマット。BioModels DB に 1,000+ モデル収録。

    手順:
      1. SBML → RoadRunner ロード
      2. 初期条件設定
      3. 時間発展シミュレーション
      4. 結果抽出・可視化

    対応モデル:
      - ODE ベース（決定論的）
      - Stochastic（Gillespie SSA）
      - Hybrid
    """
    import roadrunner

    rr = roadrunner.RoadRunner(sbml_file)
    result = rr.simulate(0, duration, n_points)

    df = pd.DataFrame(result, columns=result.colnames)
    species = [c for c in df.columns if c != "time"]

    print(f"  SBML: {len(species)} species simulated over t=[0, {duration}]")
    print(f"  Species: {', '.join(species[:5])}{'...' if len(species) > 5 else ''}")
    return df, rr


def steady_state_analysis(rr):
    """
    定常状態解析。

    定常状態: dx/dt = f(x, p) = 0
    ヤコビアン J の固有値 → 安定性判定:
      - Re(λᵢ) < 0 ∀i: 安定平衡点
      - ∃i: Re(λᵢ) > 0: 不安定
    """
    rr.steadyState()
    species_ids = rr.getFloatingSpeciesIds()
    ss_values = rr.getFloatingSpeciesConcentrations()

    # ヤコビアン
    jac = rr.getFullJacobian()
    eigenvalues = np.linalg.eigvals(jac)
    stable = all(np.real(eigenvalues) < 0)

    ss_dict = dict(zip(species_ids, ss_values))
    ss_dict["stable"] = stable
    ss_dict["eigenvalues"] = eigenvalues.tolist()

    print(f"  Steady state: {'Stable' if stable else 'Unstable'}")
    return ss_dict
```

## 2. フラックスバランス解析（FBA）

```python
def flux_balance_analysis(model_path, objective="biomass", method="fba"):
    """
    代謝フラックスバランス解析。

    FBA 定式化:
      max  c^T · v          (目的関数、通常 biomass)
      s.t. S · v = 0        (定常状態制約)
           vₘᵢₙ ≤ v ≤ vₘₐₓ  (フラックス範囲制約)

    method:
      - "fba": 標準 FBA — LP で最適フラックス分布を求める
      - "pfba": Parsimonious FBA — 最小総フラックスで最適化
      - "fva": Flux Variability Analysis — 各反応の許容フラックス範囲
      - "loopless": ループフリー FBA

    入力: SBML / JSON / YAML 形式のゲノムスケール代謝モデル（GEM）
    BiGG Models DB: 100+ 生物種の GEM を収録
    """
    import cobra

    model = cobra.io.read_sbml_model(model_path)
    print(f"  Model: {model.id} — {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

    if method == "fba":
        solution = model.optimize()
    elif method == "pfba":
        solution = cobra.flux_analysis.pfba(model)
    elif method == "fva":
        fva_result = cobra.flux_analysis.flux_variability_analysis(
            model, fraction_of_optimum=0.9)
        return fva_result

    # 結果
    objective_value = solution.objective_value
    fluxes = solution.fluxes

    # Essential genes (single gene knockouts)
    essential = []
    for gene in model.genes:
        with model:
            gene.knock_out()
            ko_sol = model.optimize()
            if ko_sol.objective_value < 0.01 * objective_value:
                essential.append(gene.id)

    print(f"  FBA: objective={objective_value:.4f}, "
          f"{len(essential)} essential genes")

    result = {
        "objective_value": objective_value,
        "n_active_reactions": (fluxes.abs() > 1e-6).sum(),
        "n_essential_genes": len(essential),
        "essential_genes": essential,
    }
    return result, fluxes
```

## 3. 遺伝子制御ネットワーク推定（GRN）

```python
def infer_grn(expression_matrix, method="genie3", n_top=1000):
    """
    遺伝子制御ネットワーク（GRN）推定。

    method:
      - "genie3": GENIE3 — Random Forest ベース
        各遺伝子 gⱼ を他の全遺伝子で回帰し、
        特徴量重要度を制御関係の重みとする。
      - "scenic": SCENIC — cis-regulatory 解析統合
      - "granger": Granger 因果性 — 時系列データ向け

    GENIE3 原理:
      For each target gene gⱼ:
        Train RF: gⱼ = f(g₁, ..., gⱼ₋₁, gⱼ₊₁, ..., gₚ)
        Weight wᵢⱼ = importance of gᵢ for predicting gⱼ
    """
    from arboreto.algo import genie3

    # GENIE3
    if method == "genie3":
        network = genie3(expression_matrix.values,
                          gene_names=expression_matrix.columns.tolist())
        network = network.sort_values("importance", ascending=False).head(n_top)

    # ネットワーク構築
    import networkx as nx
    G = nx.DiGraph()
    for _, row in network.iterrows():
        G.add_edge(row["TF"], row["target"], weight=row["importance"])

    # ハブ TF（高出次数）
    out_degrees = sorted(G.out_degree(), key=lambda x: x[1], reverse=True)
    top_tfs = out_degrees[:10]

    print(f"  GRN: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"  Top TFs: {', '.join([tf for tf, _ in top_tfs[:5]])}")
    return G, network
```

## 4. パラメータ推定・感度解析

```python
from scipy.optimize import differential_evolution
from SALib.sample import saltelli
from SALib.analyze import sobol

def parameter_estimation(model_func, data, param_bounds, method="de"):
    """
    ODE モデルパラメータ推定。

    method:
      - "de": Differential Evolution（グローバル最適化）
      - "mcmc": MCMC — pymc/emcee（事後分布推定）

    目的関数: Σ (y_data - y_model)² / σ² → minimize
    """
    def objective(params):
        y_pred = model_func(params)
        residuals = (data["y_obs"] - y_pred) ** 2
        return np.sum(residuals / data.get("sigma", 1) ** 2)

    if method == "de":
        result = differential_evolution(objective, bounds=param_bounds,
                                         seed=42, maxiter=1000, tol=1e-8)
        return {
            "params": result.x,
            "cost": result.fun,
            "success": result.success,
            "message": result.message,
        }


def global_sensitivity_analysis(model_func, param_names, param_bounds,
                                  n_samples=1024):
    """
    Sobol グローバル感度解析。

    指標:
      - S1: 一次感度指標（主効果）
      - ST: 全次感度指標（主効果＋交互作用全て）
      - S2: 二次感度指標（ペアワイズ交互作用）

    S1 + 交互作用 = ST
    ΣS1 < 1 の場合、交互作用効果が存在する。
    """
    problem = {
        "num_vars": len(param_names),
        "names": param_names,
        "bounds": param_bounds,
    }

    param_values = saltelli.sample(problem, n_samples)
    Y = np.array([model_func(p) for p in param_values])

    Si = sobol.analyze(problem, Y)

    sa_df = pd.DataFrame({
        "parameter": param_names,
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })

    print(f"  Sensitivity: top parameter = {sa_df.loc[sa_df['ST'].idxmax(), 'parameter']} "
          f"(ST={sa_df['ST'].max():.3f})")
    return sa_df, Si
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/simulation_timecourse.csv` | CSV |
| `results/fba_fluxes.csv` | CSV |
| `results/grn_network.csv` | CSV |
| `results/sensitivity_analysis.csv` | CSV |
| `results/parameter_estimates.json` | JSON |
| `figures/timecourse_plot.png` | PNG |
| `figures/flux_map.png` | PNG |
| `figures/grn_graph.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| BioModels | `biomodels_search` | SBML モデル検索 |
| BioModels | `BioModels_get_model` | モデル詳細取得 |
| BioModels | `BioModels_download_model` | モデルダウンロード |
| Reactome | `Reactome_get_pathway` | パスウェイ情報取得 |
| Reactome | `Reactome_get_pathway_reactions` | 反応一覧取得 |
| Reactome | `Reactome_map_uniprot_to_pathways` | UniProt→パスウェイ |
| BiGG | `BiGG_search` | 代謝モデル検索 |
| BiGG | `BiGG_get_model` | GEM モデル取得 |
| BiGG | `BiGG_get_reaction` | 反応詳細取得 |
| KEGG | `kegg_get_pathway_info` | KEGG パスウェイ |
| KEGG | `kegg_get_gene_info` | KEGG 遺伝子情報 |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-network-analysis](../scientific-network-analysis/SKILL.md) | GRN ネットワーク解析 |
| [scientific-multi-omics](../scientific-multi-omics/SKILL.md) | マルチオミクスデータ統合 |
| [scientific-bayesian-statistics](../scientific-bayesian-statistics/SKILL.md) | ベイズパラメータ推定 |
| [scientific-doe](../scientific-doe/SKILL.md) | 実験設計・感度解析 |
| [scientific-metabolomics](../scientific-metabolomics/SKILL.md) | 代謝フラックス-メタボローム統合 |

#### 依存パッケージ

- cobra (cobrapy), roadrunner (libroadrunner), arboreto, SALib, scipy, networkx
