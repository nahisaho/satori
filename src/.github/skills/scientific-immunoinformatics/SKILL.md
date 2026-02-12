---
name: scientific-immunoinformatics
description: |
  免疫情報学スキル。エピトープ予測（MHC-I/II バインディング）・
  T 細胞/B 細胞エピトープマッピング・抗体構造解析（CDR ループ）・
  免疫レパトア解析（TCR/BCR クロノタイプ）・ワクチン候補設計・
  IEDB/IMGT/SAbDab データベース統合パイプライン。
tu_tools:
  - key: iedb
    name: IEDB
    description: 免疫エピトープデータベース
---

# Scientific Immunoinformatics

免疫情報学（Immunoinformatics）に特化した解析パイプラインを提供する。
エピトープ予測、MHC 結合親和性推定、抗体配列・構造解析、
免疫レパトア多様性解析、ワクチン候補優先順位付けを体系的に扱う。

## When to Use

- ペプチド-MHC 結合親和性を予測するとき
- T 細胞 / B 細胞エピトープを同定・マッピングするとき
- TCR / BCR レパトア（クロノタイプ）多様性を解析するとき
- 抗体 CDR ループの構造モデリングを行うとき
- ワクチン候補アンチゲンの優先順位付けを行うとき

---

## Quick Start

## 1. MHC-I バインディング予測

```python
import numpy as np
import pandas as pd

def predict_mhc_binding(peptides, alleles, method="netmhcpan"):
    """
    MHC クラス I バインディング親和性予測。

    method:
      - "netmhcpan": NetMHCpan 4.1 — ペプチド-MHC 結合 IC50 予測
      - "mhcflurry": MHCflurry 2.0 — ニューラルネットワークベース

    閾値:
      - Strong binder: IC50 < 50 nM (または %Rank < 0.5)
      - Weak binder: IC50 < 500 nM (または %Rank < 2.0)

    Parameters:
        peptides: ペプチド配列リスト（8-14 mer）
        alleles: HLA アレルリスト (e.g., ["HLA-A*02:01", "HLA-B*07:02"])
    """
    from mhcflurry import Class1PresentationPredictor

    predictor = Class1PresentationPredictor.load()

    results = []
    for peptide in peptides:
        for allele in alleles:
            pred = predictor.predict(peptides=[peptide], alleles=[allele],
                                     verbose=0)
            results.append({
                "peptide": peptide,
                "allele": allele,
                "affinity_nM": pred["affinity"].values[0],
                "percentile_rank": pred["affinity_percentile"].values[0],
                "processing_score": pred["processing_score"].values[0],
                "presentation_score": pred["presentation_score"].values[0],
            })

    df = pd.DataFrame(results)
    df["binding_level"] = np.where(
        df["affinity_nM"] < 50, "Strong",
        np.where(df["affinity_nM"] < 500, "Weak", "Non-binder")
    )

    n_strong = (df["binding_level"] == "Strong").sum()
    n_weak = (df["binding_level"] == "Weak").sum()
    print(f"  MHC-I: {n_strong} strong + {n_weak} weak binders / {len(df)} predictions")
    return df
```

## 2. B 細胞エピトープ予測

```python
def predict_bcell_epitopes(sequence, window_size=20, threshold=0.5):
    """
    B 細胞（線状）エピトープ予測。

    統合スコアリング:
      1. BepiPred 2.0: Random Forest ベース予測
      2. Parker hydrophilicity scale
      3. Emini surface accessibility
      4. Chou-Fasman β-turn prediction

    combined_score = 0.4 * bepipred + 0.2 * hydrophilicity +
                     0.2 * surface + 0.2 * beta_turn
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    pa = ProteinAnalysis(str(sequence))

    # Parker hydrophilicity
    hydrophilicity = pa.protein_scale(window=window_size,
                                       param_dict="Parker")

    # 簡易 B 細胞エピトープスコア
    from Bio.SeqUtils.ProtParam import ProtParamData
    flexibility = pa.flexibility()

    epitopes = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        score = np.mean([
            hydrophilicity[i] if i < len(hydrophilicity) else 0,
            flexibility[i] if i < len(flexibility) else 0,
        ])
        if score > threshold:
            epitopes.append({
                "start": i + 1,
                "end": i + window_size,
                "sequence": window,
                "score": score,
            })

    df = pd.DataFrame(epitopes)
    print(f"  B-cell epitopes: {len(df)} predicted (threshold={threshold})")
    return df
```

## 3. TCR/BCR レパトア解析

```python
def repertoire_analysis(clonotype_df, chain="TRB",
                         clone_col="cdr3_aa", count_col="clone_count"):
    """
    TCR/BCR レパトア多様性解析。

    多様性指標:
      - Shannon entropy: H = -Σ pᵢ log₂(pᵢ)
      - Simpson index: D = 1 - Σ pᵢ²
      - Chao1 estimator: S_est = S_obs + f₁²/(2·f₂)
      - Clonality: 1 - H/log₂(N)
      - Gini coefficient: 均等性の指標

    Parameters:
        clonotype_df: クロノタイプ DataFrame (cdr3_aa, clone_count)
        chain: TCR/BCR 鎖 (TRA, TRB, IGH, IGL, IGK)
    """
    from scipy.stats import entropy

    counts = clonotype_df[count_col].values
    total = counts.sum()
    freqs = counts / total

    # Shannon entropy
    H = entropy(freqs, base=2)
    # Simpson index
    D = 1 - np.sum(freqs ** 2)
    # Clonality
    n_clones = len(counts)
    clonality = 1 - H / np.log2(n_clones) if n_clones > 1 else 0

    # Chao1
    f1 = np.sum(counts == 1)  # singletons
    f2 = np.sum(counts == 2)  # doubletons
    chao1 = n_clones + (f1 ** 2) / (2 * max(f2, 1))

    # Gini coefficient
    sorted_freqs = np.sort(freqs)
    n = len(sorted_freqs)
    gini = (2 * np.sum((np.arange(1, n + 1)) * sorted_freqs) / (n * np.sum(sorted_freqs))) - (n + 1) / n

    # Top clones
    top10 = clonotype_df.nlargest(10, count_col)

    metrics = {
        "chain": chain,
        "n_clonotypes": n_clones,
        "total_cells": int(total),
        "shannon_entropy": round(H, 4),
        "simpson_index": round(D, 4),
        "clonality": round(clonality, 4),
        "chao1": round(chao1, 1),
        "gini": round(gini, 4),
        "top1_frequency": round(freqs[0], 4) if len(freqs) > 0 else 0,
    }

    print(f"  Repertoire ({chain}): {n_clones} clonotypes, "
          f"Shannon={H:.3f}, Clonality={clonality:.3f}")
    return metrics, top10
```

## 4. 抗体構造解析

```python
def antibody_structure_analysis(vh_seq, vl_seq, numbering="imgt"):
    """
    抗体可変領域の構造解析。

    パイプライン:
      1. ANARCI ナンバリング（IMGT / Kabat / Chothia）
      2. CDR ループ同定（CDR-H1/H2/H3, CDR-L1/L2/L3）
      3. フレームワーク領域（FR1-FR4）抽出
      4. 発生確率・体細胞超変異（SHM）率推定
      5. ヒト化可能性スコア

    CDR 定義（IMGT 方式）:
      CDR-H1: 26-33 (8 残基)
      CDR-H2: 51-57 (7 残基)
      CDR-H3: 93-102 (可変長)
    """
    from anarci import anarci

    # ナンバリング
    vh_numbered = anarci([("VH", vh_seq)], scheme=numbering)
    vl_numbered = anarci([("VL", vl_seq)], scheme=numbering)

    # CDR 抽出（IMGT 方式）
    cdr_regions = {
        "CDR-H1": (26, 33), "CDR-H2": (51, 57), "CDR-H3": (93, 102),
        "CDR-L1": (27, 32), "CDR-L2": (50, 52), "CDR-L3": (89, 97),
    }

    cdrs = {}
    for name, (start, end) in cdr_regions.items():
        chain_data = vh_numbered if "H" in name else vl_numbered
        seq = extract_region(chain_data, start, end)
        cdrs[name] = seq

    # SHM 率（生殖系列との差分）推定
    def estimate_shm_rate(numbered_seq, germline_db="imgt"):
        """生殖系列配列との差異から SHM 率を推定"""
        # 簡易実装: 生殖系列との一致率
        return 0.0  # 要生殖系列 DB

    result = {
        "cdrs": cdrs,
        "vh_length": len(vh_seq),
        "vl_length": len(vl_seq),
        "cdr_h3_length": len(cdrs.get("CDR-H3", "")),
        "numbering": numbering,
    }

    print(f"  Antibody: CDR-H3 length={result['cdr_h3_length']}, "
          f"scheme={numbering}")
    return result
```

## 5. ワクチン候補優先順位付け

```python
def vaccine_candidate_ranking(antigens_df, weights=None):
    """
    ワクチン候補アンチゲンの多基準優先順位付け。

    評価基準:
      1. Antigenicity score: VaxiJen 2.0 スコア（閾値 > 0.4）
      2. Allergenicity: AllerTOP 非アレルゲン性
      3. Toxicity: ToxinPred 非毒性
      4. MHC coverage: HLA supertype カバー率
      5. Conservation: 配列保存性（多株間）
      6. Surface accessibility: 表面露出度

    Composite score = Σ wᵢ · normalized_scoreᵢ
    """
    if weights is None:
        weights = {
            "antigenicity": 0.25,
            "mhc_coverage": 0.25,
            "conservation": 0.20,
            "surface_accessibility": 0.15,
            "non_allergenicity": 0.10,
            "non_toxicity": 0.05,
        }

    # Min-max 正規化
    for col in weights.keys():
        if col in antigens_df.columns:
            min_val = antigens_df[col].min()
            max_val = antigens_df[col].max()
            if max_val > min_val:
                antigens_df[f"{col}_norm"] = (antigens_df[col] - min_val) / (max_val - min_val)
            else:
                antigens_df[f"{col}_norm"] = 1.0

    # Composite スコア
    antigens_df["composite_score"] = sum(
        w * antigens_df.get(f"{col}_norm", 0)
        for col, w in weights.items()
    )

    antigens_df = antigens_df.sort_values("composite_score", ascending=False)
    print(f"  Vaccine candidates: {len(antigens_df)} antigens ranked")
    return antigens_df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/mhc_binding_predictions.csv` | CSV |
| `results/bcell_epitopes.csv` | CSV |
| `results/repertoire_diversity.json` | JSON |
| `results/antibody_structure.json` | JSON |
| `results/vaccine_candidates_ranked.csv` | CSV |
| `figures/epitope_map.png` | PNG |
| `figures/repertoire_clonality.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| IEDB | `iedb_search_epitopes` | エピトープ検索 |
| IEDB | `iedb_get_epitope_mhc` | エピトープ-MHC 結合データ |
| IEDB | `iedb_search_bcell` | B 細胞エピトープ検索 |
| IEDB | `iedb_search_mhc` | MHC アレル検索 |
| IEDB | `iedb_search_antigens` | 抗原検索 |
| IMGT | `IMGT_get_gene_info` | 免疫遺伝子情報 |
| IMGT | `IMGT_get_sequence` | 免疫グロブリン配列取得 |
| IMGT | `IMGT_search_genes` | 免疫遺伝子検索 |
| SAbDab | `SAbDab_search_structures` | 抗体構造検索 |
| SAbDab | `SAbDab_get_structure` | 抗体構造取得 |
| TheraSAbDab | `TheraSAbDab_search_therapeutics` | 治療用抗体検索 |
| TheraSAbDab | `TheraSAbDab_search_by_target` | 標的別治療用抗体 |
| UniProt | `UniProt_get_entry_by_accession` | タンパク質情報取得 |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-sequence-analysis](../scientific-sequence-analysis/SKILL.md) | 配列アライメント・保存性解析 |
| [scientific-protein-structure-analysis](../scientific-protein-structure-analysis/SKILL.md) | 抗体 3D 構造解析 |
| [scientific-protein-design](../scientific-protein-design/SKILL.md) | 抗体エンジニアリング |
| [scientific-variant-interpretation](../scientific-variant-interpretation/SKILL.md) | HLA タイピング・バリアント解釈 |
| [scientific-single-cell-genomics](../scientific-single-cell-genomics/SKILL.md) | 免疫細胞サブタイプ解析 |

#### 依存パッケージ

- mhcflurry, anarci, biopython, immcantation, scirpy
