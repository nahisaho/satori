---
name: scientific-protein-domain-family
description: |
  タンパク質ドメイン・ファミリー解析スキル。InterPro アノテーション検索、
  InterProScan によるシーケンスベースドメイン予測、Pfam/SMART/CDD
  ドメイン分類、ドメインアーキテクチャ可視化、ファミリー系統樹構築。
---

# Scientific Protein Domain & Family Analysis

InterPro / InterProScan を中心としたタンパク質ドメイン解析
およびファミリー分類パイプラインを提供する。

## When to Use

- 未知タンパク質のドメイン構成を同定するとき
- InterPro/Pfam アノテーションでファミリー分類するとき
- ドメインアーキテクチャを比較・可視化するとき
- InterProScan でバッチシーケンス解析するとき
- 進化的なドメイン保存性を評価するとき

---

## Quick Start

## 1. InterPro ドメイン検索

```python
import requests
import pandas as pd
import json


INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"


def search_interpro_domains(query, db_filter=None):
    """
    InterPro REST API でドメイン検索。

    Parameters:
        query: str — ドメイン名またはキーワード (e.g., "kinase", "SH3")
        db_filter: str — "pfam", "smart", "cdd", "prosite" など

    ToolUniverse:
        InterPro_search_domains(query=query)
    """
    url = f"{INTERPRO_API}/entry/interpro"
    params = {"search": query, "page_size": 25}

    if db_filter:
        url = f"{INTERPRO_API}/entry/{db_filter}"

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for entry in data.get("results", []):
        meta = entry.get("metadata", {})
        results.append({
            "accession": meta.get("accession", ""),
            "name": meta.get("name", ""),
            "type": meta.get("type", ""),
            "source_database": meta.get("source_database", "interpro"),
            "description": (meta.get("description", [{}])[0].get("text", "")[:200]
                          if meta.get("description") else ""),
            "member_databases": meta.get("member_databases", {}),
        })

    df = pd.DataFrame(results)
    print(f"InterPro search '{query}': {len(df)} entries found")
    return df
```

## 2. タンパク質のドメイン組成取得

```python
def get_protein_domains(uniprot_id):
    """
    UniProt タンパク質のドメインアノテーションを取得。

    Parameters:
        uniprot_id: str — UniProt accession (e.g., "P04637")

    ToolUniverse:
        InterPro_get_protein_domains(uniprot_id=uniprot_id)
    """
    url = f"{INTERPRO_API}/protein/uniprot/{uniprot_id}"
    params = {"page_size": 50}

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    domains = []
    for entry in data.get("results", []):
        meta = entry.get("metadata", {})
        # 各ドメインの位置情報
        for protein_info in entry.get("proteins", []):
            for loc in protein_info.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    domains.append({
                        "accession": meta.get("accession", ""),
                        "name": meta.get("name", ""),
                        "type": meta.get("type", ""),
                        "start": frag.get("start", 0),
                        "end": frag.get("end", 0),
                        "source": meta.get("source_database", "interpro"),
                    })

    df = pd.DataFrame(domains)
    if not df.empty:
        df = df.sort_values("start").reset_index(drop=True)

    print(f"Protein {uniprot_id}: {len(df)} domain annotations")
    return df
```

## 3. InterProScan シーケンスベース予測

```python
import time


IPRSCAN_API = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"


def submit_interproscan(sequence, email="user@example.com",
                          applications=None):
    """
    InterProScan REST API でシーケンスベースドメイン予測。

    Parameters:
        sequence: str — アミノ酸配列 (FASTA or raw)
        email: str — 結果通知メール
        applications: list — ["Pfam", "SMART", "CDD", "ProSitePatterns"]

    ToolUniverse:
        InterProScan_scan_sequence(sequence=sequence)
        InterProScan_get_job_status(job_id=job_id)
        InterProScan_get_job_results(job_id=job_id)
    """
    if applications is None:
        applications = ["Pfam", "SMART", "CDD", "ProSitePatterns",
                         "SUPERFAMILY", "Gene3D"]

    # 1. Submit job
    payload = {
        "email": email,
        "sequence": sequence,
        "appl": ",".join(applications),
        "goterms": "true",
        "pathways": "true",
    }
    resp = requests.post(f"{IPRSCAN_API}/run", data=payload)
    resp.raise_for_status()
    job_id = resp.text.strip()
    print(f"InterProScan job submitted: {job_id}")

    # 2. Poll for results
    max_wait = 600  # 10 min
    elapsed = 0
    while elapsed < max_wait:
        status_resp = requests.get(f"{IPRSCAN_API}/status/{job_id}")
        status = status_resp.text.strip()
        if status == "FINISHED":
            break
        elif status in ("ERROR", "FAILURE", "NOT_FOUND"):
            raise RuntimeError(f"InterProScan job {job_id}: {status}")
        time.sleep(30)
        elapsed += 30
        print(f"  Waiting... ({elapsed}s, status={status})")

    # 3. Get results
    result_resp = requests.get(f"{IPRSCAN_API}/result/{job_id}/json")
    result_resp.raise_for_status()
    results = result_resp.json()

    # Parse matches
    matches = []
    for res in results.get("results", [results]):
        for match in res.get("matches", []):
            sig = match.get("signature", {})
            for loc in match.get("locations", []):
                matches.append({
                    "accession": sig.get("accession", ""),
                    "name": sig.get("name", ""),
                    "description": sig.get("description", ""),
                    "database": sig.get("signatureLibraryRelease", {}).get("library", ""),
                    "start": loc.get("start", 0),
                    "end": loc.get("end", 0),
                    "score": loc.get("score", None),
                    "evalue": loc.get("evalue", None),
                    "interpro_accession": (
                        sig.get("entry", {}).get("accession", "")
                        if sig.get("entry") else ""
                    ),
                    "interpro_name": (
                        sig.get("entry", {}).get("name", "")
                        if sig.get("entry") else ""
                    ),
                })

    df = pd.DataFrame(matches)
    print(f"InterProScan: {len(df)} domain matches from {len(applications)} DBs")
    return df, job_id
```

## 4. ドメインアーキテクチャ可視化

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def visualize_domain_architecture(domains_df, protein_length=None,
                                    output="figures/domain_architecture.png"):
    """
    ドメインアーキテクチャ図を生成。

    Parameters:
        domains_df: DataFrame — columns: [name, start, end, source]
        protein_length: int — 全長 (None なら max(end) 使用)
        output: str — 出力ファイルパス
    """
    import os
    os.makedirs(os.path.dirname(output), exist_ok=True)

    if protein_length is None:
        protein_length = int(domains_df["end"].max()) + 10

    # 色パレット
    colors = plt.cm.Set3.colors
    unique_sources = domains_df["source"].unique() if "source" in domains_df.columns else [""]
    source_color = {s: colors[i % len(colors)] for i, s in enumerate(unique_sources)}

    fig, ax = plt.subplots(figsize=(14, 3))

    # バックボーン
    ax.plot([0, protein_length], [0.5, 0.5], "k-", linewidth=2)

    # ドメインボックス
    for _, row in domains_df.iterrows():
        start = row["start"]
        end = row["end"]
        source = row.get("source", "")
        color = source_color.get(source, "steelblue")

        rect = mpatches.FancyBboxPatch(
            (start, 0.2), end - start, 0.6,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="black", linewidth=1
        )
        ax.add_patch(rect)

        # ラベル
        mid = (start + end) / 2
        label = row.get("name", row.get("accession", ""))
        if len(label) > 12:
            label = label[:10] + ".."
        ax.text(mid, 0.5, label, ha="center", va="center",
                fontsize=7, fontweight="bold")

    ax.set_xlim(-10, protein_length + 10)
    ax.set_ylim(-0.5, 1.5)
    ax.set_xlabel("Residue position")
    ax.set_title("Domain Architecture")
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Domain architecture: {output}")
    return output
```

## 5. ドメインファミリー比較

```python
def compare_domain_architectures(protein_list):
    """
    複数タンパク質のドメインアーキテクチャを比較。

    Parameters:
        protein_list: list of str — UniProt accessions
    """
    all_domains = {}
    for uniprot_id in protein_list:
        try:
            df = get_protein_domains(uniprot_id)
            all_domains[uniprot_id] = df
        except Exception as e:
            print(f"  Warning: {uniprot_id} failed — {e}")
            all_domains[uniprot_id] = pd.DataFrame()

    # 共通ドメイン分析
    domain_sets = {}
    for prot, df in all_domains.items():
        if not df.empty:
            domain_sets[prot] = set(df["accession"].unique())
        else:
            domain_sets[prot] = set()

    # 全タンパク質に共通するドメイン
    if domain_sets:
        common = set.intersection(*domain_sets.values()) if domain_sets else set()
        all_unique = set.union(*domain_sets.values()) if domain_sets else set()
    else:
        common = set()
        all_unique = set()

    summary = {
        "proteins_analyzed": len(protein_list),
        "proteins_with_domains": sum(
            1 for df in all_domains.values() if not df.empty
        ),
        "common_domains": sorted(common),
        "total_unique_domains": len(all_unique),
        "per_protein": {
            prot: {"count": len(s), "domains": sorted(s)}
            for prot, s in domain_sets.items()
        },
    }

    print(f"Domain comparison: {len(protein_list)} proteins, "
          f"{len(common)} common domains, "
          f"{len(all_unique)} unique domains total")
    return summary, all_domains
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/interpro_search.csv` | CSV |
| `results/protein_domains.csv` | CSV |
| `results/interproscan_results.csv` | CSV |
| `results/domain_comparison.json` | JSON |
| `figures/domain_architecture.png` | PNG |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| InterPro | `InterPro_search_domains` | ドメイン検索 |
| InterPro | `InterPro_get_protein_domains` | タンパク質ドメイン取得 |
| InterPro | `InterPro_get_domain_details` | ドメイン詳細 |
| InterProScan | `InterProScan_scan_sequence` | 配列ベース予測 |
| InterProScan | `InterProScan_get_job_status` | ジョブステータス |
| InterProScan | `InterProScan_get_job_results` | 結果取得 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-protein-structure-analysis` | 3D 構造解析 |
| `scientific-protein-interaction-network` | PPI ネットワーク |
| `scientific-sequence-alignment` | 多重配列アラインメント |
| `scientific-phylogenetics` | 系統樹構築 |
| `scientific-gene-expression-transcriptomics` | 発現相関 |

### 依存パッケージ

`requests`, `pandas`, `matplotlib`, `json` (stdlib), `time` (stdlib)
