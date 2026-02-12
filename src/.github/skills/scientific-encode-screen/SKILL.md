---
name: scientific-encode-screen
description: |
  ENCODE / ChIP-Atlas エピゲノムアトラススキル。ENCODE REST API
  実験・ファイル・バイオサンプル検索、SCREEN cis 制御エレメント、
  ChIP-Atlas エンリッチメント解析、エピゲノムアノテーション統合。
  ToolUniverse 連携: encode, chipatlas。
tu_tools:
  - key: encode
    name: ENCODE
    description: ENCODE プロジェクトの実験・バイオサンプル・ファイルデータ検索
  - key: chipatlas
    name: ChIP-Atlas
    description: ChIP-seq/ATAC-seq エンリッチメント解析・ピークブラウザ
---

# Scientific ENCODE / SCREEN / ChIP-Atlas

ENCODE REST API / SCREEN / ChIP-Atlas を活用したエピゲノム
アトラス統合解析パイプラインを提供する。

## When to Use

- ENCODE 実験データ (ChIP-seq/ATAC-seq/DNase) を検索するとき
- SCREEN cCRE (candidate cis-Regulatory Elements) を照会するとき
- ChIP-Atlas で転写因子結合・ヒストン修飾エンリッチメントを解析するとき
- バリアントの制御領域アノテーションを強化するとき
- エピゲノムデータを下流のクロマチン解析に統合するとき
- 組織/細胞型特異的エピゲノムプロファイルを比較するとき

---

## Quick Start

## 1. ENCODE 実験検索・ファイル取得

```python
import requests
import pandas as pd
import json

ENCODE_BASE = "https://www.encodeproject.org"
HEADERS = {"Accept": "application/json"}


def encode_search_experiments(assay_title=None, biosample=None,
                               target=None, organism="Homo sapiens",
                               limit=50):
    """
    ENCODE — 実験メタデータ検索。

    Parameters:
        assay_title: str — アッセイ種別 (例: "ChIP-seq", "ATAC-seq")
        biosample: str — バイオサンプル (例: "K562")
        target: str — ターゲット (例: "CTCF")
        organism: str — 生物種
        limit: int — 最大結果数
    """
    url = f"{ENCODE_BASE}/search/"
    params = {
        "type": "Experiment",
        "status": "released",
        "replicates.library.biosample.donor.organism.scientific_name": organism,
        "limit": limit,
        "format": "json",
    }
    if assay_title:
        params["assay_title"] = assay_title
    if biosample:
        params["biosample_ontology.term_name"] = biosample
    if target:
        params["target.label"] = target

    resp = requests.get(url, params=params, headers=HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for exp in data.get("@graph", []):
        results.append({
            "accession": exp.get("accession", ""),
            "assay": exp.get("assay_title", ""),
            "biosample": exp.get("biosample_summary", ""),
            "target": exp.get("target", {}).get("label", ""),
            "status": exp.get("status", ""),
            "lab": exp.get("lab", {}).get("title", ""),
            "date_released": exp.get("date_released", ""),
            "files_count": len(exp.get("files", [])),
        })

    df = pd.DataFrame(results)
    print(f"ENCODE: {len(df)} experiments found")
    return df


def encode_get_files(experiment_accession, file_type="bigWig",
                      output_type="signal p-value", assembly="GRCh38"):
    """
    ENCODE — 実験ファイル URL 取得。

    Parameters:
        experiment_accession: str — 実験アクセッション
        file_type: str — ファイルタイプ
        output_type: str — 出力タイプ
        assembly: str — アセンブリ
    """
    url = f"{ENCODE_BASE}/search/"
    params = {
        "type": "File",
        "dataset": f"/experiments/{experiment_accession}/",
        "file_format": file_type,
        "output_type": output_type,
        "assembly": assembly,
        "status": "released",
        "format": "json",
    }
    resp = requests.get(url, params=params, headers=HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    files = []
    for f in data.get("@graph", []):
        files.append({
            "accession": f.get("accession", ""),
            "file_format": f.get("file_format", ""),
            "output_type": f.get("output_type", ""),
            "assembly": f.get("assembly", ""),
            "href": ENCODE_BASE + f.get("href", ""),
            "file_size": f.get("file_size", 0),
            "biological_replicate": f.get("biological_replicates", []),
        })

    df = pd.DataFrame(files)
    print(f"ENCODE files ({experiment_accession}): {len(df)} files")
    return df
```

## 2. SCREEN cCRE 検索

```python
SCREEN_BASE = "https://api.wenglab.org/screen/v2"


def screen_ccre_search(gene=None, region=None, biosample=None,
                        assembly="GRCh38"):
    """
    SCREEN — candidate cis-Regulatory Element 検索。

    Parameters:
        gene: str — 遺伝子名 (例: "TP53")
        region: str — ゲノム領域 (例: "chr17:7668421-7687490")
        biosample: str — バイオサンプル名
        assembly: str — アセンブリ
    """
    url = f"{SCREEN_BASE}/search"
    query = {"assembly": assembly}

    if gene:
        query["gene"] = gene
    if region:
        chrom, pos = region.split(":")
        start, end = pos.split("-")
        query["coordinates"] = {
            "chromosome": chrom,
            "start": int(start),
            "end": int(end),
        }

    resp = requests.post(url, json=query, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    ccres = []
    for cre in data.get("data", []):
        ccres.append({
            "accession": cre.get("accession", ""),
            "chromosome": cre.get("chrom", ""),
            "start": cre.get("start", 0),
            "end": cre.get("end", 0),
            "ccre_class": cre.get("group", ""),
            "dnase_zscore": cre.get("dnase_zscore", None),
            "h3k4me3_zscore": cre.get("h3k4me3_zscore", None),
            "h3k27ac_zscore": cre.get("h3k27ac_zscore", None),
            "ctcf_zscore": cre.get("ctcf_zscore", None),
        })

    df = pd.DataFrame(ccres)
    print(f"SCREEN cCREs: {len(df)} elements")
    return df
```

## 3. ChIP-Atlas エンリッチメント解析

```python
CHIPATLAS_BASE = "https://chip-atlas.org/api"


def chipatlas_enrichment(gene_list, cell_type=None,
                          antigen_class="TFs and others",
                          genome="hg38", threshold=5):
    """
    ChIP-Atlas — 遺伝子リストのエンリッチメント解析。

    Parameters:
        gene_list: list[str] — 遺伝子リスト
        cell_type: str — 細胞型 (None = 全細胞型)
        antigen_class: str — 抗原クラス
        genome: str — ゲノムアセンブリ
        threshold: int — 距離閾値 (kb)
    """
    url = f"{CHIPATLAS_BASE}/enrichment"
    payload = {
        "genome": genome,
        "geneList": gene_list,
        "antigenClass": antigen_class,
        "distanceThreshold": threshold * 1000,
    }
    if cell_type:
        payload["cellType"] = cell_type

    resp = requests.post(url, json=payload, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for hit in data.get("results", []):
        results.append({
            "antigen": hit.get("antigen", ""),
            "cell_type": hit.get("cellType", ""),
            "experiment_id": hit.get("experimentId", ""),
            "p_value": hit.get("pValue", 1.0),
            "log_p": hit.get("logPValue", 0),
            "overlap_genes": hit.get("overlapGenes", 0),
            "total_peaks": hit.get("totalPeaks", 0),
        })

    df = pd.DataFrame(results)
    df = df.sort_values("p_value")
    print(f"ChIP-Atlas enrichment: {len(df)} TF/antigen hits")
    return df
```

## 4. ENCODE + SCREEN + ChIP-Atlas 統合パイプライン

```python
def encode_epigenome_pipeline(gene_name, biosample="K562",
                               output_dir="results"):
    """
    ENCODE/SCREEN/ChIP-Atlas エピゲノム統合パイプライン。

    Parameters:
        gene_name: str — 遺伝子名
        biosample: str — バイオサンプル
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) SCREEN cCRE
    ccres = screen_ccre_search(gene=gene_name)
    ccres.to_csv(output_dir / "screen_ccres.csv", index=False)

    # 2) ENCODE 実験
    experiments = encode_search_experiments(
        assay_title="ChIP-seq",
        biosample=biosample,
        target="H3K27ac",
    )
    experiments.to_csv(output_dir / "encode_experiments.csv", index=False)

    # 3) ChIP-Atlas エンリッチメント
    enrichment = chipatlas_enrichment(
        gene_list=[gene_name],
        cell_type=biosample,
    )
    enrichment.to_csv(output_dir / "chipatlas_enrichment.csv", index=False)

    print(f"ENCODE epigenome pipeline: {output_dir}")
    return {
        "ccres": ccres,
        "experiments": experiments,
        "enrichment": enrichment,
    }
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `encode` | ENCODE | 実験・バイオサンプル・ファイル検索 |
| `chipatlas` | ChIP-Atlas | 転写因子エンリッチメント解析 |

## パイプライン統合

```
regulatory-genomics → encode-screen → epigenomics-chromatin
  (RegulomeDB/ReMap)   (ENCODE/SCREEN)  (ChIP/ATAC bulk)
       │                     │                ↓
  variant-interpretation ───┘         scatac-signac
  (ACMG バリアント)     │              (scATAC-seq)
                        ↓
                  gene-regulatory-network
                  (GRN 推定)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/screen_ccres.csv` | cCRE アノテーション | → variant-interpretation |
| `results/encode_experiments.csv` | ENCODE 実験メタデータ | → epigenomics-chromatin |
| `results/chipatlas_enrichment.csv` | TF エンリッチメント | → gene-regulatory-network |
