---
name: scientific-arrayexpress-expression
description: |
  ArrayExpress 発現アーカイブスキル。BioStudies/ArrayExpress
  REST API によるマイクロアレイ・RNA-seq 発現実験検索・メタ
  データ取得・データ再解析。ToolUniverse 連携: arrayexpress。
tu_tools:
  - key: arrayexpress
    name: ArrayExpress
    description: ArrayExpress 発現実験検索・メタデータ・ファイル取得
---

# Scientific ArrayExpress Expression

EBI ArrayExpress / BioStudies REST API を活用した発現データ
アーカイブ検索・再解析パイプラインを提供する。

## When to Use

- ArrayExpress/BioStudies の発現実験を検索するとき
- マイクロアレイ/RNA-seq 発現データのメタデータを取得するとき
- SDRF サンプル情報テーブルを解析するとき
- E-MTAB/E-GEOD アクセッションからデータ再解析するとき
- 発現データアーカイブを横断検索するとき
- GEO と ArrayExpress の両方でデータを探すとき

---

## Quick Start

## 1. BioStudies 発現実験検索

```python
import requests
import pandas as pd

BIOSTUDIES_BASE = "https://www.ebi.ac.uk/biostudies/api/v1"
AE_BASE = "https://www.ebi.ac.uk/arrayexpress/json/v3"


def arrayexpress_search_experiments(query, organism=None,
                                       experiment_type=None,
                                       limit=50):
    """
    ArrayExpress — 発現実験検索 (BioStudies API)。

    Parameters:
        query: str — 検索クエリ (例: "breast cancer RNA-seq")
        organism: str — 生物種 (例: "Homo sapiens")
        experiment_type: str — 実験タイプ (例: "RNA-seq of coding RNA")
        limit: int — 最大結果数
    """
    url = f"{BIOSTUDIES_BASE}/search"
    params = {
        "query": query,
        "type": "study",
        "pageSize": limit,
    }
    if organism:
        params["organism"] = organism
    if experiment_type:
        params["experimenttype"] = experiment_type

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    hits = data.get("hits", [])
    results = []
    for h in hits:
        attrs = {a.get("name", ""): a.get("value", "")
                 for a in h.get("attributes", [])}
        results.append({
            "accession": h.get("accession", ""),
            "title": attrs.get("Title", h.get("title", "")),
            "organism": attrs.get("Organism", ""),
            "experiment_type": attrs.get("Experiment type", ""),
            "release_date": h.get("releaseDate", ""),
            "files_count": h.get("filesCount", 0),
            "links_count": h.get("linksCount", 0),
        })

    df = pd.DataFrame(results)
    print(f"ArrayExpress search: {len(df)} experiments "
          f"(query={query})")
    return df
```

## 2. 実験メタデータ・SDRF 取得

```python
def arrayexpress_get_experiment(accession):
    """
    ArrayExpress — 実験メタデータ & SDRF 取得。

    Parameters:
        accession: str — アクセッション (例: "E-MTAB-12345")
    """
    url = f"{BIOSTUDIES_BASE}/studies/{accession}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    # メタデータ
    attrs = {a.get("name", ""): a.get("value", "")
             for a in data.get("attributes", [])}
    metadata = {
        "accession": accession,
        "title": attrs.get("Title", ""),
        "description": attrs.get("Description", "")[:500],
        "organism": attrs.get("Organism", ""),
        "experiment_type": attrs.get("Experiment type", ""),
        "release_date": data.get("releaseDate", ""),
    }

    # ファイル一覧
    files = []
    for section in data.get("section", {}).get("files", []):
        if isinstance(section, list):
            for f in section:
                files.append({
                    "filename": f.get("path", ""),
                    "type": f.get("type", ""),
                    "size": f.get("size", 0),
                })
        elif isinstance(section, dict):
            files.append({
                "filename": section.get("path", ""),
                "type": section.get("type", ""),
                "size": section.get("size", 0),
            })

    files_df = pd.DataFrame(files)

    # SDRF 取得試行
    sdrf_url = (f"https://www.ebi.ac.uk/biostudies/files/"
                f"{accession}/{accession}.sdrf.txt")
    sdrf_df = pd.DataFrame()
    try:
        sdrf_resp = requests.get(sdrf_url, timeout=30)
        if sdrf_resp.status_code == 200:
            from io import StringIO
            sdrf_df = pd.read_csv(StringIO(sdrf_resp.text), sep="\t")
    except Exception:
        pass

    print(f"ArrayExpress {accession}: {len(files_df)} files, "
          f"{len(sdrf_df)} SDRF rows")
    return metadata, files_df, sdrf_df
```

## 3. 発現データダウンロード・処理

```python
def arrayexpress_download_matrix(accession, output_dir="results"):
    """
    ArrayExpress — 発現マトリクスダウンロード。

    Parameters:
        accession: str — アクセッション
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata, files_df, sdrf_df = arrayexpress_get_experiment(accession)

    # 処理済み発現ファイル検索
    expr_files = files_df[
        files_df["filename"].str.contains(
            r"processed|normalized|expression|counts",
            case=False, na=False)
    ]

    downloaded = []
    for _, frow in expr_files.iterrows():
        fname = frow["filename"]
        url = (f"https://www.ebi.ac.uk/biostudies/files/"
               f"{accession}/{fname}")
        try:
            resp = requests.get(url, timeout=120)
            if resp.status_code == 200:
                fpath = output_dir / fname.split("/")[-1]
                fpath.write_bytes(resp.content)
                downloaded.append(str(fpath))
        except Exception:
            continue

    # SDRF 保存
    if not sdrf_df.empty:
        sdrf_df.to_csv(output_dir / "sdrf.csv", index=False)

    print(f"ArrayExpress download: {len(downloaded)} files → "
          f"{output_dir}")
    return {
        "metadata": metadata,
        "files": downloaded,
        "sdrf": sdrf_df,
    }
```

## 4. ArrayExpress 統合パイプライン

```python
def arrayexpress_pipeline(query, organism="Homo sapiens",
                            output_dir="results"):
    """
    ArrayExpress 統合パイプライン。

    Parameters:
        query: str — 検索クエリ
        organism: str — 生物種
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 実験検索
    experiments = arrayexpress_search_experiments(
        query, organism=organism)
    experiments.to_csv(output_dir / "experiments.csv", index=False)

    # 2) トップ実験の詳細
    if not experiments.empty:
        top_acc = experiments.iloc[0]["accession"]
        metadata, files, sdrf = arrayexpress_get_experiment(top_acc)
        files.to_csv(output_dir / "experiment_files.csv", index=False)
        if not sdrf.empty:
            sdrf.to_csv(output_dir / "sdrf.csv", index=False)

    print(f"ArrayExpress pipeline: {output_dir}")
    return {"experiments": experiments}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `arrayexpress` | ArrayExpress | 発現実験検索・メタデータ・ファイル取得 |

## パイプライン統合

```
ebi-databases → arrayexpress-expression → gene-expression-transcriptomics
  (EBI Search)   (ArrayExpress/BioStudies)  (DESeq2/GSEA)
       │                   │                       ↓
  geo-expression ─────────┘               pathway-enrichment
  (GEO データ)        │                   (KEGG/Reactome)
                       ↓
                  multi-omics
                  (統合解析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/experiments.csv` | 実験一覧 | → geo-expression |
| `results/sdrf.csv` | サンプル情報 | → gene-expression-transcriptomics |
| `results/experiment_files.csv` | ファイルリスト | → data-preprocessing |
