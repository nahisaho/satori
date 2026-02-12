---
name: scientific-noncoding-rna
description: |
  非コード RNA (ncRNA) 解析スキル。Rfam RNA ファミリー検索、
  RNAcentral 統合 ncRNA データベース、共分散モデル、構造マッピング、
  系統樹解析パイプライン。
---

# Scientific Noncoding RNA

Rfam および RNAcentral を活用した ncRNA ファミリー検索、
配列アノテーション、構造予測パイプラインを提供する。

## When to Use

- RNA ファミリー (miRNA, lncRNA, rRNA, tRNA 等) を分類するとき
- Rfam 共分散モデルで RNA 配列を検索するとき
- RNAcentral で ncRNA のクロスリファレンスを取得するとき
- RNA 二次構造・構造マッピング情報を取得するとき
- RNA ファミリーの系統樹情報を調べるとき

---

## Quick Start

## 1. Rfam ファミリー検索

```python
import requests
import pandas as pd

RFAM_API = "https://rfam.org/family"


def get_rfam_family(rfam_acc):
    """
    Rfam RNA ファミリーの詳細情報を取得。

    Parameters:
        rfam_acc: str — Rfam accession (e.g., "RF00001") or ID

    ToolUniverse:
        Rfam_get_family(rfam_acc=rfam_acc)
        Rfam_id_to_accession(rfam_id=rfam_id)
    """
    url = f"https://rfam.org/family/{rfam_acc}?content-type=application/json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    info = data.get("rfam", {}).get("acc", {})
    desc = data.get("rfam", {}).get("description", "")

    print(f"Rfam {rfam_acc}: {data.get('rfam', {}).get('id', '?')}")
    return data
```

## 2. Rfam 配列検索 (Infernal cmscan)

```python
import time


def rfam_sequence_search(sequence, email=None):
    """
    Rfam に RNA 配列を投入し Infernal cmscan で
    マッチする RNA ファミリーを同定。

    Parameters:
        sequence: str — RNA sequence

    ToolUniverse:
        Rfam_search_sequence(sequence=sequence)
    """
    url = "https://rfam.org/search/sequence"

    payload = {
        "seq": sequence,
        "output": "json",
    }
    resp = requests.post(url, data=payload)
    resp.raise_for_status()

    # Async job → poll
    job_url = resp.json().get("resultURL", "")
    if not job_url:
        return resp.json()

    for _ in range(30):
        time.sleep(10)
        result = requests.get(job_url)
        if result.status_code == 200:
            data = result.json()
            if data.get("status", "") == "DONE":
                hits = data.get("hits", {}).get("hit", [])
                print(f"Rfam cmscan: {len(hits)} family hits")
                return hits

    print("Rfam cmscan: timeout")
    return []
```

## 3. Rfam 構造マッピング

```python
def get_rfam_structure_mapping(rfam_acc):
    """
    Rfam ファミリーの PDB 構造マッピング情報を取得。

    ToolUniverse:
        Rfam_get_structure_mapping(rfam_acc=rfam_acc)
        Rfam_get_covariance_model(rfam_acc=rfam_acc)
        Rfam_get_tree_data(rfam_acc=rfam_acc)
        Rfam_get_sequence_regions(rfam_acc=rfam_acc)
    """
    # Structure mapping
    url_struct = (
        f"https://rfam.org/family/{rfam_acc}/structures"
        "?content-type=application/json"
    )
    resp_s = requests.get(url_struct)
    structures = resp_s.json() if resp_s.status_code == 200 else []

    # Sequence regions
    url_regions = (
        f"https://rfam.org/family/{rfam_acc}/regions"
        "?content-type=application/json"
    )
    resp_r = requests.get(url_regions)
    regions = resp_r.json() if resp_r.status_code == 200 else []

    print(f"Rfam {rfam_acc}: {len(structures)} PDB structures, "
          f"{len(regions) if isinstance(regions, list) else '?'} regions")
    return structures, regions
```

## 4. RNAcentral ncRNA 検索

```python
RNACENTRAL_API = "https://rnacentral.org/api/v1"


def rnacentral_search(query, page_size=10):
    """
    RNAcentral で ncRNA を検索。

    Parameters:
        query: str — search term (gene name, accession, keyword)

    ToolUniverse:
        RNAcentral_search(query=query)
    """
    url = f"{RNACENTRAL_API}/rna/"
    params = {"query": query, "page_size": page_size}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    results = data.get("results", [])
    entries = []
    for r in results:
        entries.append({
            "rnacentral_id": r.get("rnacentral_id", ""),
            "description": r.get("description", ""),
            "rna_type": r.get("rna_type", ""),
            "length": r.get("length", 0),
            "num_xrefs": r.get("xref_count", 0),
        })

    df = pd.DataFrame(entries)
    print(f"RNAcentral '{query}': {data.get('count', 0)} total, "
          f"{len(df)} returned")
    return df


def rnacentral_get_by_accession(accession):
    """
    RNAcentral アクセッションから ncRNA 詳細情報を取得。

    ToolUniverse:
        RNAcentral_get_by_accession(accession=accession)
    """
    url = f"{RNACENTRAL_API}/rna/{accession}/"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    print(f"RNAcentral {accession}: {data.get('description', '')}")
    return data
```

## 5. ncRNA 統合解析パイプライン

```python
def ncRNA_integrated_search(sequence, rfam_acc=None):
    """
    配列ベースの ncRNA 統合解析。

    ToolUniverse (横断):
        Rfam_search_sequence(sequence) → Rfam_get_family(rfam_acc)
        RNAcentral_search(query)
    """
    pipeline = {"sequence_length": len(sequence)}

    # Step 1: Rfam family identification
    rfam_hits = rfam_sequence_search(sequence)
    pipeline["rfam_hits"] = len(rfam_hits) if isinstance(rfam_hits, list) else 0

    # Step 2: If Rfam family found, get details
    if rfam_hits and isinstance(rfam_hits, list) and len(rfam_hits) > 0:
        top_hit = rfam_hits[0]
        top_acc = top_hit.get("acc", rfam_acc or "")
        if top_acc:
            family = get_rfam_family(top_acc)
            pipeline["rfam_family"] = top_acc

    # Step 3: RNAcentral search
    rna_df = rnacentral_search(sequence[:30])  # truncate for search
    pipeline["rnacentral_hits"] = len(rna_df)

    print(f"ncRNA pipeline: Rfam={pipeline.get('rfam_family', 'none')}, "
          f"RNAcentral={pipeline['rnacentral_hits']} hits")
    return pipeline
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/rfam_family.json` | JSON |
| `results/rfam_cmscan_hits.json` | JSON |
| `results/rfam_structures.json` | JSON |
| `results/rnacentral_search.csv` | CSV |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| Rfam | `Rfam_get_family` | ファミリー情報 |
| Rfam | `Rfam_search_sequence` | 配列→ファミリー同定 |
| Rfam | `Rfam_get_covariance_model` | 共分散モデル |
| Rfam | `Rfam_get_structure_mapping` | PDB マッピング |
| Rfam | `Rfam_get_tree_data` | 系統樹 |
| Rfam | `Rfam_get_sequence_regions` | 配列領域 |
| Rfam | `Rfam_id_to_accession` | ID→アクセッション変換 |
| RNAcentral | `RNAcentral_search` | ncRNA 検索 |
| RNAcentral | `RNAcentral_get_by_accession` | 詳細取得 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-gene-expression-transcriptomics` | 転写産物解析 |
| `scientific-genome-sequence-tools` | 配列取得 |
| `scientific-structural-proteomics` | RNA 構造 |
| `scientific-biothings-idmapping` | ID マッピング |

### 依存パッケージ

`requests`, `pandas`
