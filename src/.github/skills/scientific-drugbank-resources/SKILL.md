---
name: scientific-drugbank-resources
description: |
  DrugBank リソーススキル。DrugBank API を用いた薬剤記述・
  薬理情報・標的タンパク質・薬物相互作用検索。
  ToolUniverse 連携: drugbank。
tu_tools:
  - key: drugbank
    name: DrugBank
    description: 薬剤データベース API
---

# Scientific DrugBank Resources

DrugBank API を活用した薬剤記述・薬理情報 (MOA)・標的タンパク質
検索・薬物相互作用 (DDI) パイプラインを提供する。

## When to Use

- 薬剤の基本情報 (名前・分類・構造) を検索するとき
- 薬理メカニズム (MOA) を調べるとき
- 標的タンパク質から薬剤を逆引き検索するとき
- 薬物相互作用 (DDI) を確認するとき
- 薬剤の ADMET プロパティを取得するとき

---

## Quick Start

## 1. 薬剤検索・基本情報

```python
import requests
import pandas as pd

DRUGBANK_API = "https://api.drugbank.com/v1"


def drugbank_search(query, limit=25, api_key=None):
    """
    DrugBank — 薬剤テキスト検索。

    Parameters:
        query: str — 検索クエリ (例: "imatinib")
        limit: int — 最大結果数
        api_key: str — DrugBank API キー
    """
    headers = {}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    url = f"{DRUGBANK_API}/drugs"
    params = {"q": query, "per_page": limit}
    resp = requests.get(url, params=params,
                        headers=headers, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for d in data:
        rows.append({
            "drugbank_id": d.get("drugbank_id", ""),
            "name": d.get("name", ""),
            "cas_number": d.get("cas_number", ""),
            "drug_type": d.get("type", ""),
            "state": d.get("state", ""),
            "description": (d.get("description", "")
                            [:200]),
        })

    df = pd.DataFrame(rows)
    print(f"DrugBank search: '{query}' → {len(df)} drugs")
    return df


def drugbank_drug_detail(drugbank_id, api_key=None):
    """
    DrugBank — 薬剤詳細取得。

    Parameters:
        drugbank_id: str — DrugBank ID (例: "DB01254")
        api_key: str — DrugBank API キー
    """
    headers = {}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    url = f"{DRUGBANK_API}/drugs/{drugbank_id}"
    resp = requests.get(url, headers=headers, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    result = {
        "drugbank_id": data.get("drugbank_id", ""),
        "name": data.get("name", ""),
        "description": data.get("description", ""),
        "indication": data.get("indication", ""),
        "pharmacodynamics": data.get(
            "pharmacodynamics", ""),
        "mechanism_of_action": data.get(
            "mechanism_of_action", ""),
        "absorption": data.get("absorption", ""),
        "half_life": data.get("half_life", ""),
        "protein_binding": data.get(
            "protein_binding", ""),
        "molecular_weight": data.get(
            "average_molecular_weight", ""),
    }
    return result
```

## 2. 標的タンパク質検索

```python
def drugbank_targets(drugbank_id, api_key=None):
    """
    DrugBank — 薬剤の標的タンパク質取得。

    Parameters:
        drugbank_id: str — DrugBank ID
        api_key: str — DrugBank API キー
    """
    headers = {}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    url = f"{DRUGBANK_API}/drugs/{drugbank_id}/targets"
    resp = requests.get(url, headers=headers, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for t in data:
        polypeptide = t.get("polypeptide", {}) or {}
        rows.append({
            "drugbank_id": drugbank_id,
            "target_name": t.get("name", ""),
            "organism": t.get("organism", ""),
            "known_action": t.get("known_action", ""),
            "gene_name": polypeptide.get(
                "gene_name", ""),
            "uniprot_id": polypeptide.get(
                "external_identifiers", {}).get(
                "UniProtKB", ""),
        })

    df = pd.DataFrame(rows)
    print(f"DrugBank targets: {drugbank_id} "
          f"→ {len(df)} targets")
    return df
```

## 3. 薬物相互作用 (DDI)

```python
def drugbank_interactions(drugbank_id, api_key=None):
    """
    DrugBank — 薬物相互作用取得。

    Parameters:
        drugbank_id: str — DrugBank ID
        api_key: str — DrugBank API キー
    """
    headers = {}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    url = (f"{DRUGBANK_API}/drugs/"
           f"{drugbank_id}/drug_interactions")
    resp = requests.get(url, headers=headers, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for inter in data:
        rows.append({
            "drug_a": drugbank_id,
            "drug_b_id": inter.get(
                "drugbank_id", ""),
            "drug_b_name": inter.get("name", ""),
            "description": inter.get(
                "description", "")[:300],
        })

    df = pd.DataFrame(rows)
    print(f"DrugBank DDI: {drugbank_id} "
          f"→ {len(df)} interactions")
    return df
```

## 4. DrugBank 統合パイプライン

```python
def drugbank_pipeline(drug_name, api_key=None,
                         output_dir="results"):
    """
    DrugBank 統合パイプライン。

    Parameters:
        drug_name: str — 薬剤名 (例: "imatinib")
        api_key: str — DrugBank API キー
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 検索
    results = drugbank_search(drug_name,
                              api_key=api_key)
    results.to_csv(output_dir / "drugbank_search.csv",
                   index=False)

    if results.empty:
        print(f"DrugBank: '{drug_name}' not found")
        return {"search": results}

    db_id = results.iloc[0]["drugbank_id"]

    # 2) 詳細
    detail = drugbank_drug_detail(db_id,
                                  api_key=api_key)
    pd.DataFrame([detail]).to_csv(
        output_dir / "drugbank_detail.csv",
        index=False)

    # 3) 標的
    targets = drugbank_targets(db_id,
                               api_key=api_key)
    targets.to_csv(output_dir / "drugbank_targets.csv",
                   index=False)

    # 4) DDI
    ddi = drugbank_interactions(db_id,
                                api_key=api_key)
    ddi.to_csv(output_dir / "drugbank_ddi.csv",
               index=False)

    print(f"DrugBank pipeline: {drug_name} → {output_dir}")
    return {"detail": detail, "targets": targets,
            "ddi": ddi}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `drugbank` | DrugBank | 薬剤データベース API |

## パイプライン統合

```
drug-target-profiling → drugbank-resources → admet-pharmacokinetics
  (標的プロファイリング)   (DrugBank API)     (ADMET 予測)
          │                     │                  ↓
opentargets-genetics ──────────┘       compound-screening
  (OT 薬剤エビデンス)                  (ZINC 化合物検索)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/drugbank_search.csv` | 薬剤検索結果 | → drug-target-profiling |
| `results/drugbank_detail.csv` | 薬剤詳細 | → admet-pharmacokinetics |
| `results/drugbank_targets.csv` | 標的タンパク質 | → protein-interaction-network |
| `results/drugbank_ddi.csv` | 薬物相互作用 | → pharmacogenomics |
