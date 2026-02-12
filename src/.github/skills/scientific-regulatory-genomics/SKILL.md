---
name: scientific-regulatory-genomics
description: |
  レギュラトリーゲノミクススキル。RegulomeDB バリアント制御機能スコア、
  ReMap 転写因子結合マッピング、4D Nucleome (4DN) 三次元ゲノム構造
  解析の統合パイプライン。
---

# Scientific Regulatory Genomics

RegulomeDB / ReMap / 4D Nucleome を統合した
レギュラトリーゲノミクス (制御領域バリアント解析) パイプラインを提供する。

## When to Use

- 非コード領域バリアントの制御機能を評価するとき
- RegulomeDB で SNP の調節的影響をスコアリングするとき
- ReMap で転写因子結合部位のマッピングを確認するとき
- 4DN データから三次元ゲノム構造 (TAD/ループ) を解析するとき
- GWAS ヒットの制御メカニズムを解明するとき

---

## Quick Start

## 1. RegulomeDB バリアント制御スコア

```python
import requests
import pandas as pd

REGULOMEDB_API = "https://regulomedb.org/regulome-search"


def score_regulome_variants(variants):
    """
    RegulomeDB — 非コード領域バリアントの制御機能スコアリング。

    Parameters:
        variants: list — バリアントリスト (rsID or chr:pos 形式)
            e.g., ["rs12345", "chr1:109274570"]

    ToolUniverse:
        RegulomeDB_score_variant(variant=variant)
    """
    results = []
    for variant in variants:
        params = {"regions": variant, "genome": "GRCh38", "format": "json"}
        resp = requests.get(REGULOMEDB_API, params=params)
        if resp.status_code != 200:
            results.append({"variant": variant, "score": None, "error": True})
            continue

        data = resp.json()
        for hit in data.get("@graph", []):
            results.append({
                "variant": variant,
                "regulome_score": hit.get("regulome_score", {}).get("ranking", ""),
                "probability": hit.get("regulome_score", {}).get("probability", ""),
                "chrom": hit.get("chrom", ""),
                "start": hit.get("start", ""),
                "end": hit.get("end", ""),
                "dnase": hit.get("dnase", ""),
                "proteins_binding": hit.get("proteins_binding", []),
                "motifs": hit.get("motifs", []),
                "eqtls": hit.get("eqtls", []),
                "chromatin_state": hit.get("chromatin_state", {}),
            })

    df = pd.DataFrame(results)
    if not df.empty and "regulome_score" in df.columns:
        high_func = (df["regulome_score"].astype(str).str.match(r"^[12]")).sum()
        print(f"RegulomeDB: {len(variants)} variants scored, "
              f"{high_func} with high regulatory function (score 1-2)")
    return df
```

## 2. ReMap 転写因子結合マッピング

```python
REMAP_API = "https://remap.univ-amu.fr/api/v1"


def search_remap_binding(chrom, start, end, genome="hg38"):
    """
    ReMap — ゲノム領域の転写因子/コレギュレーター結合マッピング。

    Parameters:
        chrom: str — 染色体 (e.g., "chr1")
        start: int — 開始座標
        end: int — 終了座標
        genome: str — ゲノムアセンブリ ("hg38", "hg19", "mm10")

    ToolUniverse:
        ReMap_search_peaks(chrom=chrom, start=start, end=end)
        ReMap_get_tf_targets(tf_name=tf_name)
    """
    params = {
        "chrom": chrom,
        "start": start,
        "end": end,
        "genome": genome,
    }
    resp = requests.get(f"{REMAP_API}/peaks/search", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for peak in data.get("peaks", []):
        results.append({
            "tf_name": peak.get("tf_name", ""),
            "biotype": peak.get("biotype", ""),
            "cell_type": peak.get("cell_type", ""),
            "experiment": peak.get("experiment_accession", ""),
            "peak_start": peak.get("start", ""),
            "peak_end": peak.get("end", ""),
            "score": peak.get("score", ""),
        })

    df = pd.DataFrame(results)
    unique_tfs = df["tf_name"].nunique() if not df.empty else 0
    print(f"ReMap {chrom}:{start}-{end}: {len(df)} peaks, {unique_tfs} unique TFs")
    return df


def get_remap_tf_targets(tf_name, genome="hg38"):
    """
    ReMap — 特定転写因子の全結合部位取得。

    Parameters:
        tf_name: str — 転写因子名 (e.g., "TP53", "CTCF", "STAT3")
    """
    params = {"tf": tf_name, "genome": genome}
    resp = requests.get(f"{REMAP_API}/peaks/tf", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for peak in data.get("peaks", [])[:1000]:  # Limit for large TFs
        results.append({
            "chrom": peak.get("chrom", ""),
            "start": peak.get("start", ""),
            "end": peak.get("end", ""),
            "cell_type": peak.get("cell_type", ""),
            "score": peak.get("score", ""),
        })

    df = pd.DataFrame(results)
    print(f"ReMap TF '{tf_name}': {len(df)} binding sites")
    return df
```

## 3. 4D Nucleome (4DN) 三次元ゲノム構造

```python
FOURDN_API = "https://data.4dnucleome.org"


def search_4dn_experiments(query, experiment_type=None):
    """
    4D Nucleome ポータル — 三次元ゲノム実験データ検索。

    Parameters:
        query: str — 検索クエリ (細胞株名、タンパク質名等)
        experiment_type: str — 実験タイプ ("in situ Hi-C", "SPRITE", "GAM")

    ToolUniverse:
        FourDN_search_experiments(query=query)
    """
    params = {
        "searchTerm": query,
        "type": "ExperimentSetReplicate",
        "format": "json",
    }
    if experiment_type:
        params["experiment_type.display_title"] = experiment_type

    resp = requests.get(f"{FOURDN_API}/search/", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("@graph", []):
        results.append({
            "accession": item.get("accession", ""),
            "title": item.get("display_title", ""),
            "experiment_type": item.get("experiment_type", {}).get("display_title", ""),
            "biosource": item.get("biosource_summary", ""),
            "lab": item.get("lab", {}).get("display_title", ""),
            "status": item.get("status", ""),
        })

    df = pd.DataFrame(results)
    print(f"4DN search '{query}': {len(df)} experiment sets")
    return df
```

## 4. 制御バリアント統合解析パイプライン

```python
def regulatory_variant_pipeline(variants, genome="hg38"):
    """
    制御領域バリアント統合解析。

    Parameters:
        variants: list — バリアントリスト (rsID or chr:pos)
    """
    print("=" * 60)
    print("Regulatory Variant Analysis Pipeline")
    print("=" * 60)

    # Step 1: RegulomeDB scoring
    print("\n[1/3] RegulomeDB scoring...")
    regulome_df = score_regulome_variants(variants)

    # Step 2: ReMap TF binding for high-scoring variants
    print("\n[2/3] ReMap TF binding analysis...")
    remap_results = {}
    for _, row in regulome_df.iterrows():
        if row.get("chrom") and row.get("start"):
            chrom = row["chrom"]
            start = int(row["start"]) - 500
            end = int(row["end"]) + 500
            try:
                remap_df = search_remap_binding(chrom, start, end, genome)
                remap_results[row["variant"]] = remap_df
            except Exception as e:
                print(f"  ReMap error for {row['variant']}: {e}")

    # Step 3: Summary
    print("\n[3/3] Summary")
    summary = {
        "total_variants": len(variants),
        "regulome_scored": len(regulome_df),
        "high_regulatory": (
            regulome_df["regulome_score"].astype(str).str.match(r"^[12]")
        ).sum() if "regulome_score" in regulome_df.columns else 0,
        "remap_annotated": len(remap_results),
    }
    print(f"  Total: {summary['total_variants']}, "
          f"High regulatory: {summary['high_regulatory']}, "
          f"ReMap annotated: {summary['remap_annotated']}")

    return {"regulome": regulome_df, "remap": remap_results, "summary": summary}
```

---

## 利用可能ツール

| ToolUniverse カテゴリ | 主なツール |
|---|---|
| `regulomedb` | `RegulomeDB_score_variant` |
| `remap` | `ReMap_search_peaks`, `ReMap_get_tf_targets` |
| `fourdn_portal` | `FourDN_search_experiments` |

## パイプライン出力

| 出力ファイル | 説明 | 連携先スキル |
|---|---|---|
| `results/regulome_scores.csv` | バリアント制御スコア | → variant-interpretation, variant-effect-prediction |
| `results/remap_binding.csv` | TF 結合マッピング | → epigenomics-chromatin, disease-research |
| `results/4dn_contacts.json` | 3D ゲノム構造データ | → single-cell-genomics, epigenomics-chromatin |

## パイプライン統合

```
variant-interpretation ──→ regulatory-genomics ──→ epigenomics-chromatin
  (ACMG/AMP)               (RegulomeDB/ReMap/4DN)  (ChIP-seq/ATAC)
                                    │
                                    ├──→ disease-research (GWAS enhancer)
                                    ├──→ gene-expression (eQTL/制御)
                                    └──→ noncoding-rna (ncRNA 制御)
```
