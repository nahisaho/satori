---
name: scientific-nci60-screening
description: |
  NCI-60 がん細胞株薬剤応答スキル。CellMiner API 薬剤感受性・
  NCI-60 GI50/LC50 データ・DepMap cancer dependency 統合・
  薬剤-分子マーカー相関・細胞株パネル比較解析。
---

# Scientific NCI-60 Screening

CellMiner / NCI-60 / DepMap を活用したがん細胞株薬剤応答
パイプラインを提供する。高スループットスクリーニングデータの
統合解析、薬剤感受性マーカー同定、パネル比較。

## When to Use

- NCI-60 細胞株パネルの薬剤応答 (GI50) を解析するとき
- CellMiner から化合物活性データを取得するとき
- 薬剤感受性と分子マーカー (変異/発現) の相関を調べるとき
- DepMap CRISPR/RNAi 依存性データを併用するとき
- 細胞株間の薬剤応答パターンを比較するとき
- 新規化合物のスクリーニング結果を NCI-60 と比較するとき

---

## Quick Start

## 1. CellMiner データ取得

```python
import requests
import pandas as pd
import numpy as np
from io import StringIO

CELLMINER_BASE = "https://discover.nci.nih.gov/cellminer/api"


def cellminer_drug_activity(nsc_id=None, drug_name=None):
    """
    CellMiner — NCI-60 薬剤活性データ取得。

    Parameters:
        nsc_id: str — NSC ID (例: "740")
        drug_name: str — 薬剤名 (例: "Paclitaxel")
    """
    if nsc_id:
        url = f"{CELLMINER_BASE}/compound/{nsc_id}/activity"
    elif drug_name:
        url = f"{CELLMINER_BASE}/compound/search"
        params = {"name": drug_name}
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        compounds = resp.json()
        if not compounds:
            print(f"Drug not found: {drug_name}")
            return pd.DataFrame()
        nsc_id = compounds[0].get("nsc", "")
        url = f"{CELLMINER_BASE}/compound/{nsc_id}/activity"

    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for cell_line, values in data.get("activity", {}).items():
        results.append({
            "cell_line": cell_line,
            "tissue": values.get("tissue", ""),
            "gi50_log": values.get("gi50", None),
            "tgi_log": values.get("tgi", None),
            "lc50_log": values.get("lc50", None),
        })

    df = pd.DataFrame(results)
    print(f"CellMiner: NSC {nsc_id} → {len(df)} cell lines")
    return df
```

## 2. NCI-60 バルクデータ取得

```python
def nci60_bulk_download(data_type="drug_activity"):
    """
    NCI-60 バルクデータセット取得。

    Parameters:
        data_type: str — "drug_activity", "gene_expression",
                         "mutation", "copy_number"
    """
    urls = {
        "drug_activity": "https://discover.nci.nih.gov/cellminer/download/DTP_NCI60_ZSCORE.csv",
        "gene_expression": "https://discover.nci.nih.gov/cellminer/download/GeneExpr_RMA.csv",
        "mutation": "https://discover.nci.nih.gov/cellminer/download/Exome_Mutation.csv",
    }

    url = urls.get(data_type)
    if not url:
        raise ValueError(f"Unknown data type: {data_type}")

    resp = requests.get(url, timeout=120)
    resp.raise_for_status()

    df = pd.read_csv(StringIO(resp.text))
    print(f"NCI-60 bulk: {data_type} → {df.shape}")
    return df
```

## 3. 薬剤-分子マーカー相関

```python
from scipy import stats


def drug_marker_correlation(drug_activity, molecular_data,
                             marker_type="expression", top_n=50):
    """
    薬剤感受性と分子マーカーの相関解析。

    Parameters:
        drug_activity: pd.DataFrame — GI50 データ (cell_line, gi50)
        molecular_data: pd.DataFrame — 分子データ (cell_line, gene, value)
        marker_type: str — "expression", "mutation", "copy_number"
        top_n: int — 上位相関遺伝子数
    """
    # 細胞株一致
    common_lines = set(drug_activity["cell_line"]) & set(molecular_data["cell_line"])
    drug_sub = drug_activity[drug_activity["cell_line"].isin(common_lines)]
    mol_sub = molecular_data[molecular_data["cell_line"].isin(common_lines)]

    # 遺伝子ごとの相関
    correlations = []
    genes = mol_sub["gene"].unique() if "gene" in mol_sub.columns else mol_sub.columns[1:]

    drug_values = drug_sub.set_index("cell_line")["gi50_log"]

    for gene in genes:
        if "gene" in mol_sub.columns:
            gene_data = mol_sub[mol_sub["gene"] == gene].set_index("cell_line")["value"]
        else:
            gene_data = mol_sub.set_index("cell_line")[gene]

        common = drug_values.index.intersection(gene_data.index)
        if len(common) < 10:
            continue

        r, p = stats.pearsonr(drug_values[common], gene_data[common])
        correlations.append({
            "gene": gene,
            "pearson_r": r,
            "p_value": p,
            "n_samples": len(common),
        })

    corr_df = pd.DataFrame(correlations)
    corr_df["adj_p"] = corr_df["p_value"] * len(corr_df)  # Bonferroni
    corr_df = corr_df.sort_values("p_value")

    print(f"Drug-marker correlation: {len(corr_df)} genes tested, "
          f"top |r| = {corr_df['pearson_r'].abs().max():.3f}")
    return corr_df.head(top_n)
```

## 4. 組織別薬剤応答パターン

```python
def tissue_response_pattern(drug_activity, min_lines=3):
    """
    組織別の薬剤応答パターン解析。

    Parameters:
        drug_activity: pd.DataFrame — GI50 データ
        min_lines: int — 最小細胞株数
    """
    tissue_stats = drug_activity.groupby("tissue").agg(
        n_lines=("gi50_log", "count"),
        mean_gi50=("gi50_log", "mean"),
        std_gi50=("gi50_log", "std"),
        min_gi50=("gi50_log", "min"),
        max_gi50=("gi50_log", "max"),
    ).reset_index()

    tissue_stats = tissue_stats[tissue_stats["n_lines"] >= min_lines]
    tissue_stats = tissue_stats.sort_values("mean_gi50")

    # 感受性/耐性スコア
    overall_mean = drug_activity["gi50_log"].mean()
    tissue_stats["sensitivity_z"] = (
        (tissue_stats["mean_gi50"] - overall_mean)
        / drug_activity["gi50_log"].std()
    )

    print(f"Tissue patterns: {len(tissue_stats)} tissues")
    for _, row in tissue_stats.iterrows():
        label = "Sensitive" if row["sensitivity_z"] < -0.5 else (
            "Resistant" if row["sensitivity_z"] > 0.5 else "Neutral"
        )
        print(f"  {row['tissue']}: GI50={row['mean_gi50']:.2f} ({label})")
    return tissue_stats
```

## 5. DepMap 統合スクリーニング

```python
DEPMAP_BASE = "https://depmap.org/portal/api"


def depmap_gene_dependency(gene_symbol, dataset="Chronos_Combined"):
    """
    DepMap — CRISPR/RNAi 遺伝子依存性取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル
        dataset: str — データセット名
    """
    url = f"{DEPMAP_BASE}/download/custom"
    params = {
        "gene": gene_symbol,
        "dataset": dataset,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for entry in data.get("data", []):
        results.append({
            "cell_line": entry.get("cell_line_name", ""),
            "lineage": entry.get("lineage", ""),
            "dependency_score": entry.get("score", None),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        n_dependent = (df["dependency_score"] < -0.5).sum()
        print(f"DepMap {gene_symbol}: {len(df)} lines, "
              f"{n_dependent} dependent (score < -0.5)")
    return df
```

## 6. NCI-60 統合スクリーニングパイプライン

```python
def nci60_screening_pipeline(drug_name=None, nsc_id=None,
                              target_gene=None, output_dir="results"):
    """
    NCI-60 + DepMap 統合スクリーニングパイプライン。

    Parameters:
        drug_name: str — 薬剤名
        nsc_id: str — NSC ID
        target_gene: str — 標的遺伝子 (DepMap 連携)
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) NCI-60 薬剤活性
    drug_data = cellminer_drug_activity(nsc_id=nsc_id, drug_name=drug_name)
    drug_data.to_csv(output_dir / "drug_activity.csv", index=False)

    # 2) 組織別パターン
    tissue_patterns = tissue_response_pattern(drug_data)
    tissue_patterns.to_csv(output_dir / "tissue_patterns.csv", index=False)

    # 3) 発現相関
    expr_data = nci60_bulk_download("gene_expression")
    correlations = drug_marker_correlation(drug_data, expr_data)
    correlations.to_csv(output_dir / "marker_correlations.csv", index=False)

    # 4) DepMap 連携 (標的遺伝子あれば)
    if target_gene:
        depmap_data = depmap_gene_dependency(target_gene)
        depmap_data.to_csv(output_dir / "depmap_dependency.csv", index=False)

    print(f"Pipeline complete: {output_dir}")
    return {
        "drug_activity": drug_data,
        "tissue_patterns": tissue_patterns,
        "correlations": correlations,
    }
```

---

## パイプライン統合

```
compound-screening → nci60-screening → precision-oncology
  (ZINC/VS)           (NCI-60/DepMap)    (MTB レポート)
       │                     │                 ↓
drug-target-profiling ──────┘          cancer-genomics
  (ChEMBL/DGIdb)        │              (COSMIC/DepMap)
                         ↓
                   cell-line-resources
                   (Cellosaurus)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/drug_activity.csv` | NCI-60 GI50 データ | → precision-oncology |
| `results/tissue_patterns.csv` | 組織別応答パターン | → cancer-genomics |
| `results/marker_correlations.csv` | 薬剤-マーカー相関 | → drug-target-profiling |
| `results/depmap_dependency.csv` | DepMap 依存性スコア | → cell-line-resources |
