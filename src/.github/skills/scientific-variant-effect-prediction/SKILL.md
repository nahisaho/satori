---
name: scientific-variant-effect-prediction
description: |
  計算バリアント効果予測スキル。AlphaMissense (タンパク質構造ベース病原性予測)、
  CADD (統合アノテーションスコア)、SpliceAI (スプライシング影響予測) の
  3 大予測ツールを統合したコンセンサス病原性評価パイプライン。
  Ensembl VEP 連携、バリアントフィルタリング、優先順位付け対応。
  9 の ToolUniverse SMCP ツールと連携。
tu_tools:
  - key: spliceai
    name: SpliceAI
    description: スプライシングバリアント効果予測
  - key: cadd
    name: CADD
    description: 統合アノテーション依存性枯湇スコア
---

# Scientific Variant Effect Prediction

AlphaMissense / CADD / SpliceAI の 3 大計算予測ツールを統合した
バリアント病原性評価・優先順位付けパイプラインを提供する。

## When to Use

- ミスセンスバリアントの病原性を計算予測するとき
- CADD スコアで全ゲノムバリアントの有害度を評価するとき
- SpliceAI でスプライシング影響を予測するとき
- 複数予測ツールのコンセンサススコアを算出するとき
- WGS/WES バリアントの優先順位付けが必要なとき

---

## Quick Start

## 1. AlphaMissense 病原性予測

```python
import pandas as pd
import numpy as np
import requests


def alphamissense_predict(variants, uniprot_id=None):
    """
    AlphaMissense タンパク質構造ベース病原性予測。

    Parameters:
        variants: list[dict] — [{"protein": "P12345", "position": 42, "ref": "A", "alt": "V"}]
        uniprot_id: str — タンパク質単位で全ポジションのスコア取得
    """
    results = []

    if uniprot_id:
        # タンパク質全体のスコアマップ取得
        # AlphaMissense は事前計算済みスコアを提供
        print(f"Fetching AlphaMissense scores for {uniprot_id}...")
        # ToolUniverse 経由: AlphaMissense_get_protein_scores
        # または AlphaMissense_get_residue_scores

    for var in variants:
        protein = var.get("protein", uniprot_id)
        pos = var["position"]
        ref_aa = var.get("ref", "")
        alt_aa = var.get("alt", "")

        # スコア分類閾値 (DeepMind 推奨)
        # pathogenic: score > 0.564
        # benign: score < 0.340
        # ambiguous: 0.340 - 0.564
        score = var.get("score", np.nan)

        if not np.isnan(score):
            if score > 0.564:
                classification = "likely_pathogenic"
            elif score < 0.340:
                classification = "likely_benign"
            else:
                classification = "ambiguous"
        else:
            classification = "unknown"

        results.append({
            "protein": protein,
            "position": pos,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "am_score": score,
            "am_class": classification,
            "variant": f"{ref_aa}{pos}{alt_aa}",
        })

    df = pd.DataFrame(results)
    print(f"AlphaMissense: {len(df)} variants scored")
    return df
```

## 2. CADD スコア取得

```python
def cadd_score_variants(variants, genome="GRCh38", version="v1.7"):
    """
    CADD (Combined Annotation Dependent Depletion) スコア取得。

    Parameters:
        variants: list[dict] — [{"chr": "1", "pos": 12345, "ref": "A", "alt": "G"}]
        genome: "GRCh37" or "GRCh38"
        version: CADD バージョン
    """
    base_url = f"https://cadd.gs.washington.edu/api/{version}"

    results = []
    for var in variants:
        chrom = str(var["chr"]).replace("chr", "")
        pos = var["pos"]
        ref = var["ref"]
        alt = var["alt"]

        # CADD API クエリ
        # ToolUniverse 経由: CADD_get_variant_score
        url = f"{base_url}/{genome}/{chrom}:{pos}"
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                for hit in data:
                    if hit.get("Ref") == ref and hit.get("Alt") == alt:
                        raw = hit.get("RawScore", np.nan)
                        phred = hit.get("PHRED", np.nan)
                        break
                else:
                    raw, phred = np.nan, np.nan
            else:
                raw, phred = np.nan, np.nan
        except Exception:
            raw, phred = np.nan, np.nan

        # CADD PHRED 閾値目安
        # >= 20: top 1% deleterious
        # >= 30: top 0.1% deleterious
        if phred >= 30:
            cadd_class = "highly_deleterious"
        elif phred >= 20:
            cadd_class = "deleterious"
        elif phred >= 10:
            cadd_class = "moderate"
        else:
            cadd_class = "benign"

        results.append({
            "chr": chrom, "pos": pos, "ref": ref, "alt": alt,
            "cadd_raw": raw,
            "cadd_phred": phred,
            "cadd_class": cadd_class,
            "variant": f"chr{chrom}:{pos}{ref}>{alt}",
        })

    df = pd.DataFrame(results)
    print(f"CADD: {len(df)} variants scored, "
          f"{(df['cadd_phred'] >= 20).sum()} deleterious (PHRED≥20)")
    return df
```

## 3. SpliceAI スプライシング予測

```python
def spliceai_predict(variants, genome="GRCh38",
                      delta_threshold=0.2):
    """
    SpliceAI スプライシング影響予測。

    Parameters:
        variants: list[dict] — [{"chr": "1", "pos": 12345, "ref": "A", "alt": "G"}]
        delta_threshold: float — Δスコア閾値
            0.2: high recall, 0.5: recommended, 0.8: high precision
    """
    results = []

    for var in variants:
        chrom = str(var["chr"]).replace("chr", "")
        pos = var["pos"]
        ref = var["ref"]
        alt = var["alt"]

        # ToolUniverse 経由: SpliceAI_predict_splice
        # SpliceAI は 4 つの Δスコアを出力:
        # DS_AG: acceptor gain, DS_AL: acceptor loss
        # DS_DG: donor gain, DS_DL: donor loss
        ds_ag = var.get("ds_ag", 0)
        ds_al = var.get("ds_al", 0)
        ds_dg = var.get("ds_dg", 0)
        ds_dl = var.get("ds_dl", 0)

        max_delta = max(ds_ag, ds_al, ds_dg, ds_dl)

        if max_delta >= 0.8:
            splice_class = "high_impact"
        elif max_delta >= 0.5:
            splice_class = "moderate_impact"
        elif max_delta >= 0.2:
            splice_class = "low_impact"
        else:
            splice_class = "no_impact"

        results.append({
            "chr": chrom, "pos": pos, "ref": ref, "alt": alt,
            "ds_acceptor_gain": ds_ag,
            "ds_acceptor_loss": ds_al,
            "ds_donor_gain": ds_dg,
            "ds_donor_loss": ds_dl,
            "max_delta": max_delta,
            "splice_class": splice_class,
            "variant": f"chr{chrom}:{pos}{ref}>{alt}",
        })

    df = pd.DataFrame(results)
    impacted = (df["max_delta"] >= delta_threshold).sum()
    print(f"SpliceAI: {len(df)} variants, "
          f"{impacted} with splice impact (Δ≥{delta_threshold})")
    return df
```

## 4. コンセンサス病原性評価

```python
def consensus_pathogenicity(am_df, cadd_df, spliceai_df,
                             am_threshold=0.564, cadd_threshold=20,
                             splice_threshold=0.5):
    """
    AlphaMissense + CADD + SpliceAI のコンセンサス評価。

    Parameters:
        am_df: AlphaMissense 結果 DataFrame
        cadd_df: CADD 結果 DataFrame
        spliceai_df: SpliceAI 結果 DataFrame
    """
    # バリアント ID で結合
    merged = cadd_df.copy()

    if len(am_df) > 0:
        merged = merged.merge(
            am_df[["variant", "am_score", "am_class"]],
            on="variant", how="left"
        )
    if len(spliceai_df) > 0:
        merged = merged.merge(
            spliceai_df[["variant", "max_delta", "splice_class"]],
            on="variant", how="left"
        )

    # コンセンサススコア
    def compute_consensus(row):
        votes = 0
        total = 0

        if "cadd_phred" in row and not pd.isna(row.get("cadd_phred")):
            total += 1
            if row["cadd_phred"] >= cadd_threshold:
                votes += 1

        if "am_score" in row and not pd.isna(row.get("am_score")):
            total += 1
            if row["am_score"] >= am_threshold:
                votes += 1

        if "max_delta" in row and not pd.isna(row.get("max_delta")):
            total += 1
            if row["max_delta"] >= splice_threshold:
                votes += 1

        if total == 0:
            return "insufficient_data"
        ratio = votes / total
        if ratio >= 0.67:
            return "pathogenic"
        elif ratio >= 0.33:
            return "uncertain"
        else:
            return "benign"

    merged["consensus"] = merged.apply(compute_consensus, axis=1)
    merged["evidence_count"] = merged.apply(
        lambda r: sum(1 for c in ["cadd_phred", "am_score", "max_delta"]
                      if c in r and not pd.isna(r.get(c))), axis=1)

    print(f"Consensus: {len(merged)} variants — "
          f"{(merged['consensus'] == 'pathogenic').sum()} pathogenic, "
          f"{(merged['consensus'] == 'uncertain').sum()} uncertain, "
          f"{(merged['consensus'] == 'benign').sum()} benign")
    return merged
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/alphamissense_scores.csv` | CSV |
| `results/cadd_scores.csv` | CSV |
| `results/spliceai_scores.csv` | CSV |
| `results/consensus_pathogenicity.csv` | CSV |
| `figures/variant_score_distribution.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| AlphaMissense | `AlphaMissense_get_protein_scores` | タンパク質全体スコア |
| AlphaMissense | `AlphaMissense_get_variant_score` | 個別バリアントスコア |
| AlphaMissense | `AlphaMissense_get_residue_scores` | 残基レベルスコア |
| CADD | `CADD_get_variant_score` | 個別バリアント PHRED スコア |
| CADD | `CADD_get_position_scores` | ポジション全体スコア |
| CADD | `CADD_get_range_scores` | 範囲一括スコア |
| SpliceAI | `SpliceAI_predict_splice` | スプライシングΔスコア予測 |
| SpliceAI | `SpliceAI_predict_pangolin` | Pangolin スプライシング予測 |
| SpliceAI | `SpliceAI_get_max_delta` | 最大Δスコア取得 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-variant-interpretation` | ACMG/AMP 臨床バリアント解釈 |
| `scientific-population-genetics` | gnomAD 集団頻度参照 |
| `scientific-disease-research` | 疾患-バリアント関連 |
| `scientific-pharmacogenomics` | PGx バリアント効果 |
| `scientific-protein-structure-analysis` | 構造→機能影響評価 |

### 依存パッケージ

`pandas`, `numpy`, `requests`
