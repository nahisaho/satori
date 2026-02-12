---
name: scientific-metabolic-flux
description: |
  代謝フラックス解析スキル。13C/15N 安定同位体トレーサー
  データを用いた代謝フラックス推定・EMU モデリング・
  フラックスバランス制約統合パイプライン。
---

# Scientific Metabolic Flux

13C/15N 安定同位体トレーサー実験データを用いた代謝フラックス
推定・EMU (Elementary Metabolite Unit) フレームワーク・
フラックスバランス解析 (FBA) 制約統合パイプラインを提供する。

## When to Use

- 13C 安定同位体トレーサー実験データを解析するとき
- EMU/アイソトポマーモデルを構築するとき
- MID (Mass Isotopomer Distribution) データをフィッティングするとき
- 経路別の代謝フラックスを定量するとき
- FBA 制約とトレーサーデータを統合するとき

---

## Quick Start

## 1. MID (Mass Isotopomer Distribution) データ処理

```python
import numpy as np
import pandas as pd
from scipy.optimize import minimize


def load_mid_data(mid_file, sep="\t"):
    """
    MID データ読み込み・正規化。

    Parameters:
        mid_file: str — MID データファイルパス
            (TSV: metabolite, M+0, M+1, M+2, ...)
        sep: str — 区切り文字
    """
    df = pd.read_csv(mid_file, sep=sep,
                     index_col="metabolite")

    mid_cols = [c for c in df.columns
                if c.startswith("M+")]

    for idx in df.index:
        row_sum = df.loc[idx, mid_cols].sum()
        if row_sum > 0:
            df.loc[idx, mid_cols] /= row_sum

    print(f"MID data: {len(df)} metabolites, "
          f"{len(mid_cols)} isotopomers")
    return df[mid_cols]


def natural_abundance_correction(mid_df, n_carbons):
    """
    天然同位体存在量補正。

    Parameters:
        mid_df: DataFrame — 正規化済み MID データ
        n_carbons: dict — 代謝物名→炭素数マッピング
    """
    C13_NAT = 0.011  # 13C 天然存在比

    corrected = mid_df.copy()
    for met in corrected.index:
        n_c = n_carbons.get(met, 6)
        n_iso = min(corrected.shape[1], n_c + 1)
        raw = corrected.loc[met].values[:n_iso]

        # 補正行列 (簡易)
        corr_matrix = np.zeros((n_iso, n_iso))
        for i in range(n_iso):
            for j in range(i, n_iso):
                from math import comb
                k = j - i
                remain = n_c - i
                if k <= remain:
                    corr_matrix[i, j] = (
                        comb(remain, k)
                        * C13_NAT ** k
                        * (1 - C13_NAT) ** (remain - k)
                    )

        try:
            corrected_vals = np.linalg.solve(
                corr_matrix[:n_iso, :n_iso], raw)
            corrected_vals = np.maximum(corrected_vals, 0)
            corrected_vals /= corrected_vals.sum()
            corrected.loc[met, corrected.columns[:n_iso]] = (
                corrected_vals)
        except np.linalg.LinAlgError:
            pass

    print(f"NA correction: {len(corrected)} metabolites")
    return corrected
```

## 2. EMU フラックスモデル

```python
def build_emu_model(reactions, atom_transitions):
    """
    EMU (Elementary Metabolite Unit) モデル構築。

    Parameters:
        reactions: list[dict] — 反応定義
            [{id, substrates, products, reversible}]
        atom_transitions: dict — 原子遷移マッピング
            {reaction_id: [(from_met, from_atoms,
                            to_met, to_atoms)]}
    """
    emu_network = {}

    for rxn in reactions:
        rxn_id = rxn["id"]
        transitions = atom_transitions.get(rxn_id, [])

        for from_met, f_atoms, to_met, t_atoms in (
            transitions
        ):
            emu_size = len(t_atoms)
            emu_key = (to_met, tuple(sorted(t_atoms)))

            if emu_key not in emu_network:
                emu_network[emu_key] = []

            emu_network[emu_key].append({
                "reaction": rxn_id,
                "precursor": from_met,
                "precursor_atoms": f_atoms,
                "reversible": rxn.get(
                    "reversible", False),
            })

    print(f"EMU model: {len(emu_network)} EMUs, "
          f"{len(reactions)} reactions")
    return emu_network


def simulate_mid(fluxes, emu_model, substrate_labeling,
                   metabolite):
    """
    フラックスからの MID シミュレーション。

    Parameters:
        fluxes: dict — {reaction_id: flux_value}
        emu_model: dict — EMU ネットワーク
        substrate_labeling: dict — 基質ラベリングパターン
            {metabolite: [M+0 fraction, M+1, ...]}
        metabolite: str — シミュレーション対象代謝物
    """
    relevant_emus = {
        k: v for k, v in emu_model.items()
        if k[0] == metabolite
    }

    if not relevant_emus:
        return np.array([1.0])

    max_size = max(len(k[1]) for k in relevant_emus)
    mid = np.zeros(max_size + 1)
    mid[0] = 1.0  # デフォルト: 未標識

    for emu_key, precursors in relevant_emus.items():
        emu_size = len(emu_key[1])
        for prec in precursors:
            rxn_flux = fluxes.get(prec["reaction"], 0)
            prec_label = substrate_labeling.get(
                prec["precursor"],
                [1.0] + [0.0] * emu_size)

            for i, frac in enumerate(
                prec_label[:emu_size + 1]
            ):
                if i <= max_size:
                    mid[i] += rxn_flux * frac

    mid_sum = mid.sum()
    if mid_sum > 0:
        mid /= mid_sum

    return mid
```

## 3. フラックス推定

```python
def estimate_fluxes(observed_mids, emu_model,
                      substrate_labeling,
                      initial_fluxes,
                      metabolites):
    """
    最小二乗法によるフラックス推定。

    Parameters:
        observed_mids: dict — {metabolite: np.array}
            観測 MID データ
        emu_model: dict — EMU ネットワーク
        substrate_labeling: dict — 基質ラベリング
        initial_fluxes: dict — 初期フラックス推定値
        metabolites: list — 対象代謝物リスト
    """
    flux_names = list(initial_fluxes.keys())
    x0 = [initial_fluxes[f] for f in flux_names]

    def objective(x):
        fluxes = dict(zip(flux_names, x))
        residual = 0.0
        for met in metabolites:
            if met not in observed_mids:
                continue
            obs = observed_mids[met]
            sim = simulate_mid(
                fluxes, emu_model,
                substrate_labeling, met)
            n = min(len(obs), len(sim))
            residual += np.sum(
                (obs[:n] - sim[:n]) ** 2)
        return residual

    bounds = [(0, None) for _ in flux_names]
    result = minimize(objective, x0, method="L-BFGS-B",
                      bounds=bounds)

    estimated = dict(zip(flux_names, result.x))
    print(f"Flux estimation: SSR={result.fun:.6f}, "
          f"converged={result.success}")
    return estimated, result
```

## 4. 代謝フラックス統合パイプライン

```python
def metabolic_flux_pipeline(mid_file, reactions,
                              atom_transitions,
                              substrate_labeling,
                              n_carbons,
                              output_dir="results"):
    """
    代謝フラックス統合パイプライン。

    Parameters:
        mid_file: str — MID データファイル
        reactions: list — 反応定義リスト
        atom_transitions: dict — 原子遷移マッピング
        substrate_labeling: dict — 基質ラベリング
        n_carbons: dict — 代謝物→炭素数
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) MID 読み込み・補正
    mid_raw = load_mid_data(mid_file)
    mid_corr = natural_abundance_correction(
        mid_raw, n_carbons)
    mid_corr.to_csv(output_dir / "mid_corrected.csv")

    # 2) EMU モデル構築
    emu_model = build_emu_model(
        reactions, atom_transitions)

    # 3) フラックス推定
    observed = {met: mid_corr.loc[met].values
                for met in mid_corr.index}
    init_fluxes = {r["id"]: 1.0 for r in reactions}
    fluxes, opt_result = estimate_fluxes(
        observed, emu_model, substrate_labeling,
        init_fluxes, list(observed.keys()))

    flux_df = pd.DataFrame([
        {"reaction": k, "flux": v}
        for k, v in fluxes.items()
    ])
    flux_df.to_csv(output_dir / "fluxes.csv",
                   index=False)

    print(f"Metabolic flux pipeline → {output_dir}")
    return {"fluxes": fluxes, "mid_corrected": mid_corr}
```

---

## パイプライン統合

```
metabolic-modeling → metabolic-flux → systems-biology
  (FBA/COBRA)        (13C MFA)        (統合解析)
       │                  │                ↓
flux-balance-analysis ───┘    pathway-enrichment
  (制約ベース)                (パスウェイ集積)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/mid_corrected.csv` | 補正済み MID | → metabolic-modeling |
| `results/fluxes.csv` | 推定フラックス | → systems-biology |
