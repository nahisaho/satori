---
name: scientific-geo-expression
description: |
  GEO (Gene Expression Omnibus) 発現プロファイルスキル。GEO REST
  API データセット検索・サンプル情報・発現マトリクス取得・バルク
  RNA-seq/マイクロアレイ差次的発現解析。ToolUniverse 連携: geo。
tu_tools:
  - key: geo
    name: GEO (Gene Expression Omnibus)
    description: GEO データセット・サンプル情報・発現データ検索
---

# Scientific GEO Expression

GEO REST API を活用したトランスクリプトーム発現プロファイル
解析パイプラインを提供する。

## When to Use

- GEO データセット (GDS/GSE) を検索・ダウンロードするとき
- マイクロアレイ/RNA-seq 発現マトリクスを取得するとき
- 条件間差次的発現解析 (DEG) を実行するとき
- 複数 GEO データセットを横断比較するとき
- GEO メタデータから実験条件を構造化するとき
- 再解析パイプラインで GEO データを再利用するとき

---

## Quick Start

## 1. GEO データセット検索

```python
import requests
import pandas as pd
import GEOparse
from io import StringIO

GEO_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def geo_search_datasets(query, organism="Homo sapiens",
                         study_type=None, limit=20):
    """
    GEO — データセット検索 (E-utilities)。

    Parameters:
        query: str — 検索クエリ (例: "breast cancer RNA-seq")
        organism: str — 生物種
        study_type: str — 研究タイプ ("Expression profiling by array" etc.)
        limit: int — 最大結果数
    """
    search_term = f"{query} AND {organism}[Organism]"
    if study_type:
        search_term += f' AND "{study_type}"[Study Type]'

    # ESearch
    url = f"{GEO_BASE}/esearch.fcgi"
    params = {
        "db": "gds",
        "term": search_term,
        "retmax": limit,
        "retmode": "json",
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    ids = resp.json().get("esearchresult", {}).get("idlist", [])

    if not ids:
        print("No GEO datasets found")
        return pd.DataFrame()

    # ESummary
    url = f"{GEO_BASE}/esummary.fcgi"
    params = {"db": "gds", "id": ",".join(ids), "retmode": "json"}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    summaries = resp.json().get("result", {})

    results = []
    for gds_id in ids:
        info = summaries.get(gds_id, {})
        results.append({
            "accession": info.get("accession", ""),
            "title": info.get("title", ""),
            "summary": info.get("summary", "")[:200],
            "organism": info.get("taxon", ""),
            "platform": info.get("gpl", ""),
            "sample_count": info.get("n_samples", 0),
            "series_type": info.get("gdstype", ""),
            "pub_date": info.get("pdat", ""),
        })

    df = pd.DataFrame(results)
    print(f"GEO search: {len(df)} datasets")
    return df
```

## 2. GEO 発現マトリクス取得

```python
def geo_get_expression_matrix(gse_id, log2_transform=True):
    """
    GEO — GSE 発現マトリクス取得 (GEOparse)。

    Parameters:
        gse_id: str — GSE アクセッション (例: "GSE12345")
        log2_transform: bool — log2 変換を適用するか
    """
    import numpy as np

    gse = GEOparse.get_GEO(geo=gse_id, destdir="/tmp", silent=True)

    # サンプルメタデータ
    samples = []
    for gsm_name, gsm in gse.gsms.items():
        meta = gsm.metadata
        samples.append({
            "sample_id": gsm_name,
            "title": meta.get("title", [""])[0],
            "source": meta.get("source_name_ch1", [""])[0],
            "characteristics": "; ".join(
                meta.get("characteristics_ch1", [])),
            "platform": meta.get("platform_id", [""])[0],
        })
    sample_df = pd.DataFrame(samples)

    # 発現マトリクス
    pivoted = gse.pivot_samples("VALUE")
    if pivoted.empty:
        print(f"No expression data in {gse_id}")
        return sample_df, pd.DataFrame()

    if log2_transform:
        pivoted = np.log2(pivoted.astype(float) + 1)

    print(f"GEO {gse_id}: {pivoted.shape[0]} probes × "
          f"{pivoted.shape[1]} samples")
    return sample_df, pivoted
```

## 3. 差次的発現解析

```python
from scipy import stats
import numpy as np


def geo_differential_expression(expr_matrix, group_a_samples,
                                  group_b_samples, method="ttest",
                                  fdr_threshold=0.05, lfc_threshold=1.0):
    """
    GEO — 差次的発現解析。

    Parameters:
        expr_matrix: pd.DataFrame — 発現マトリクス (genes × samples)
        group_a_samples: list[str] — グループ A サンプル ID
        group_b_samples: list[str] — グループ B サンプル ID
        method: str — "ttest" or "wilcoxon"
        fdr_threshold: float — FDR 閾値
        lfc_threshold: float — log2FC 閾値
    """
    a_data = expr_matrix[group_a_samples]
    b_data = expr_matrix[group_b_samples]

    results = []
    for gene in expr_matrix.index:
        a_vals = a_data.loc[gene].dropna().values
        b_vals = b_data.loc[gene].dropna().values

        if len(a_vals) < 2 or len(b_vals) < 2:
            continue

        lfc = b_vals.mean() - a_vals.mean()

        if method == "ttest":
            stat, pval = stats.ttest_ind(a_vals, b_vals)
        else:
            stat, pval = stats.mannwhitneyu(a_vals, b_vals,
                                             alternative="two-sided")

        results.append({
            "gene": gene,
            "log2fc": lfc,
            "mean_a": a_vals.mean(),
            "mean_b": b_vals.mean(),
            "statistic": stat,
            "p_value": pval,
        })

    df = pd.DataFrame(results)

    # FDR correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    _, df["fdr"], _, _ = multipletests(df["p_value"], method="fdr_bh")

    # DEG フィルタ
    df["is_deg"] = (df["fdr"] < fdr_threshold) & (df["log2fc"].abs() > lfc_threshold)
    n_deg = df["is_deg"].sum()
    n_up = ((df["is_deg"]) & (df["log2fc"] > 0)).sum()
    n_down = ((df["is_deg"]) & (df["log2fc"] < 0)).sum()

    print(f"DEG: {n_deg} genes (↑{n_up} / ↓{n_down}), "
          f"FDR<{fdr_threshold}, |LFC|>{lfc_threshold}")
    return df.sort_values("p_value")
```

## 4. GEO 発現プロファイリングパイプライン

```python
def geo_expression_pipeline(gse_id, group_col="condition",
                              group_a="control", group_b="treatment",
                              output_dir="results"):
    """
    GEO 発現プロファイリング統合パイプライン。

    Parameters:
        gse_id: str — GSE アクセッション
        group_col: str — グループ化カラム
        group_a: str — コントロールグループ
        group_b: str — 処理グループ
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) データ取得
    sample_df, expr = geo_get_expression_matrix(gse_id)
    sample_df.to_csv(output_dir / "samples.csv", index=False)

    # 2) グループ分割
    a_samples = sample_df[
        sample_df["source"].str.contains(group_a, case=False)
    ]["sample_id"].tolist()
    b_samples = sample_df[
        sample_df["source"].str.contains(group_b, case=False)
    ]["sample_id"].tolist()

    # 3) 差次的発現
    deg = geo_differential_expression(expr, a_samples, b_samples)
    deg.to_csv(output_dir / "deg_results.csv", index=False)

    print(f"GEO pipeline: {output_dir}")
    return {"samples": sample_df, "expression": expr, "deg": deg}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `geo` | GEO | データセット検索・サンプル情報・発現データ |

## パイプライン統合

```
ebi-databases → geo-expression → gene-expression-transcriptomics
  (ENA/EBI Search) (GEO データ)    (DESeq2/GTEx)
       │                 │                ↓
  literature-search ────┘         pathway-enrichment
  (PubMed/OpenAlex)   │          (KEGG/Reactome/GO)
                       ↓
                  multi-omics
                  (統合解析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/samples.csv` | サンプルメタデータ | → gene-expression-transcriptomics |
| `results/deg_results.csv` | 差次的発現結果 | → pathway-enrichment |
