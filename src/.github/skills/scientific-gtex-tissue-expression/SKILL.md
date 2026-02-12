---
name: scientific-gtex-tissue-expression
description: |
  GTEx 組織発現スキル。GTEx Portal REST API v2 による
  組織特異的遺伝子発現パターン解析・eQTL ルックアップ・
  多組織比較。直接 API (ToolUniverse 非連携)。
tu_tools: []
---

# Scientific GTEx Tissue Expression

GTEx (Genotype-Tissue Expression) Portal REST API v2 を活用した
組織特異的遺伝子発現解析・eQTL 検索・多組織比較パイプライン
を提供する。

## When to Use

- 遺伝子の組織特異的発現パターンを調べるとき
- 特定組織における eQTL (発現量的形質遺伝子座) を検索するとき
- 複数組織間で遺伝子発現レベルを比較するとき
- TPM (Transcripts Per Million) 発現データを取得するとき
- バリアントが遺伝子発現に与える影響を評価するとき
- 組織間の遺伝子共発現パターンを分析するとき

---

## Quick Start

## 1. 組織特異的遺伝子発現取得

```python
import requests
import pandas as pd

GTEX_BASE = "https://gtexportal.org/api/v2"


def gtex_gene_expression(gene_id, tissue=None):
    """
    GTEx — 組織別遺伝子発現 (中央値 TPM) 取得。

    Parameters:
        gene_id: str — 遺伝子シンボル or Ensembl ID
            (例: "BRCA1", "ENSG00000012048")
        tissue: str — 組織 ID (None で全組織)
            (例: "Breast_Mammary_Tissue")
    """
    url = f"{GTEX_BASE}/expression/medianGeneExpression"
    params = {
        "gencodeId": gene_id,
        "datasetId": "gtex_v8",
    }
    if tissue:
        params["tissueSiteDetailId"] = tissue

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("data", []):
        results.append({
            "gene_symbol": item.get("geneSymbol", ""),
            "gencode_id": item.get("gencodeId", ""),
            "tissue": item.get("tissueSiteDetailId", ""),
            "tissue_name": item.get("tissueSiteDetail", ""),
            "median_tpm": item.get("median", 0),
            "sample_count": item.get("numSamples", 0),
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("median_tpm", ascending=False)

    print(f"GTEx expression: {gene_id} → "
          f"{len(df)} tissues")
    return df


def gtex_top_tissues(gene_id, top_n=10):
    """
    GTEx — 発現量上位組織。

    Parameters:
        gene_id: str — 遺伝子シンボル or Ensembl ID
        top_n: int — 上位組織数
    """
    df = gtex_gene_expression(gene_id)
    top = df.head(top_n) if not df.empty else df
    print(f"GTEx top {top_n} tissues for {gene_id}:")
    for _, row in top.iterrows():
        print(f"  {row['tissue_name']}: "
              f"{row['median_tpm']:.2f} TPM "
              f"(n={row['sample_count']})")
    return top
```

## 2. eQTL 検索

```python
def gtex_eqtl_lookup(gene_id, tissue, variant_id=None):
    """
    GTEx — eQTL ルックアップ。

    Parameters:
        gene_id: str — 遺伝子シンボル or Ensembl ID
        tissue: str — 組織 ID
            (例: "Liver", "Whole_Blood")
        variant_id: str — バリアント ID (任意)
            (例: "rs12345")
    """
    url = f"{GTEX_BASE}/association/singleTissueEqtl"
    params = {
        "gencodeId": gene_id,
        "tissueSiteDetailId": tissue,
        "datasetId": "gtex_v8",
    }
    if variant_id:
        params["variantId"] = variant_id

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("data", []):
        results.append({
            "gene_symbol": item.get("geneSymbol", ""),
            "variant_id": item.get("variantId", ""),
            "tissue": tissue,
            "pvalue": item.get("pValue"),
            "nes": item.get("nes"),  # normalized effect size
            "maf": item.get("maf"),
            "ref": item.get("ref", ""),
            "alt": item.get("alt", ""),
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("pvalue")

    print(f"GTEx eQTL: {gene_id} in {tissue} → "
          f"{len(df)} associations")
    return df
```

## 3. 多組織比較

```python
def gtex_multi_gene_comparison(gene_ids, tissues=None):
    """
    GTEx — 複数遺伝子・複数組織の発現比較。

    Parameters:
        gene_ids: list[str] — 遺伝子リスト
        tissues: list[str] — 組織リスト (None で全組織)
    """
    all_data = []
    for gid in gene_ids:
        try:
            df = gtex_gene_expression(gid)
            if tissues:
                df = df[df["tissue"].isin(tissues)]
            all_data.append(df)
        except Exception as e:
            print(f"  Warning: {gid} — {e}")
            continue

    if not all_data:
        return pd.DataFrame()

    combined = pd.concat(all_data, ignore_index=True)

    # ピボットテーブル: 行=組織, 列=遺伝子, 値=TPM
    if not combined.empty:
        pivot = combined.pivot_table(
            index="tissue_name",
            columns="gene_symbol",
            values="median_tpm",
            aggfunc="first",
        )
        print(f"GTEx comparison: {len(gene_ids)} genes × "
              f"{len(pivot)} tissues")
        return pivot

    return combined
```

## 4. GTEx 統合パイプライン

```python
def gtex_pipeline(gene_ids, tissues=None,
                    output_dir="results"):
    """
    GTEx 統合パイプライン。

    Parameters:
        gene_ids: list[str] — 遺伝子リスト
        tissues: list[str] — 組織リスト (None で全組織)
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 全遺伝子の組織発現
    all_expr = []
    for gid in gene_ids:
        try:
            df = gtex_gene_expression(gid)
            df.to_csv(output_dir / f"expression_{gid}.csv",
                      index=False)
            all_expr.append(df)
        except Exception:
            continue

    # 2) 多組織比較マトリクス
    pivot = gtex_multi_gene_comparison(gene_ids, tissues)
    if isinstance(pivot, pd.DataFrame) and not pivot.empty:
        pivot.to_csv(output_dir / "expression_matrix.csv")

    # 3) eQTL 検索 (上位組織)
    eqtl_results = []
    for gid in gene_ids:
        if all_expr:
            top = all_expr[-1].head(3)
            for _, row in top.iterrows():
                try:
                    eqtl = gtex_eqtl_lookup(gid,
                                              row["tissue"])
                    eqtl_results.append(eqtl)
                except Exception:
                    continue
    if eqtl_results:
        eqtl_combined = pd.concat(eqtl_results,
                                    ignore_index=True)
        eqtl_combined.to_csv(output_dir / "eqtl_results.csv",
                              index=False)

    print(f"GTEx pipeline: {output_dir}")
    return {"expression": all_expr, "matrix": pivot}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| (direct) | GTEx Portal API v2 | 直接 REST API — TU 非連携 |

## パイプライン統合

```
gene-expression-transcriptomics → gtex-tissue-expression → variant-interpretation
  (DESeq2/edgeR 差分発現)     (組織別 TPM + eQTL)       (臨床変異評価)
          │                          │                        ↓
  arrayexpress-expression ──────────┘              gwas-catalog
  (ArrayExpress データ)          │               (GWAS 関連解析)
                                 ↓
                        disease-research
                        (疾患関連遺伝子)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/expression_*.csv` | 遺伝子別組織発現 | → disease-research |
| `results/expression_matrix.csv` | 多遺伝子比較 | → pathway-enrichment |
| `results/eqtl_results.csv` | eQTL 関連 | → variant-interpretation |
