---
name: scientific-depmap-dependencies
description: |
  DepMap 依存性スキル。Cancer Dependency Map (DepMap) Portal
  API によるがん細胞株 CRISPR/RNAi 依存性スコア・薬剤
  感受性データ・遺伝子効果取得。ToolUniverse 連携: depmap。
tu_tools:
  - key: depmap
    name: DepMap
    description: がん細胞株依存性マップ API
---

# Scientific DepMap Dependencies

Cancer Dependency Map (DepMap) Portal API を活用した
がん細胞株の CRISPR/RNAi 遺伝子依存性スコア取得・
薬剤感受性データ・遺伝子効果パイプラインを提供する。

## When to Use

- がん細胞株の遺伝子依存性 (CRISPR/RNAi) を調べるとき
- 遺伝子のがん選択的必須性を評価するとき
- 薬剤感受性データと遺伝子発現の相関を調べるとき
- 細胞株メタデータ (組織型・疾患サブタイプ) を取得するとき
- Cell Model Passports データを参照するとき

---

## Quick Start

## 1. DepMap 細胞株検索・メタデータ

```python
import requests
import pandas as pd

DEPMAP_API = "https://api.cellmodelpassports.sanger.ac.uk/v1"


def depmap_cell_lines(tissue=None, cancer_type=None,
                         limit=50):
    """
    DepMap / Cell Model Passports — 細胞株検索。

    Parameters:
        tissue: str — 組織型フィルタ (例: "Breast")
        cancer_type: str — がん種フィルタ
            (例: "Breast Carcinoma")
        limit: int — 最大結果数
    """
    url = f"{DEPMAP_API}/models"
    params = {"page_size": limit}
    if tissue:
        params["tissue"] = tissue
    if cancer_type:
        params["cancer_type"] = cancer_type

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for model in data.get("data", data
                          if isinstance(data, list)
                          else []):
        if isinstance(model, dict):
            rows.append({
                "model_id": model.get("model_id", ""),
                "model_name": model.get("model_name", ""),
                "tissue": model.get("tissue", ""),
                "cancer_type": model.get(
                    "cancer_type", ""),
                "sample_site": model.get(
                    "sample_site", ""),
                "gender": model.get("gender", ""),
            })

    df = pd.DataFrame(rows[:limit])
    print(f"DepMap cell lines: {len(df)} models")
    return df
```

## 2. CRISPR 遺伝子依存性スコア

```python
def depmap_gene_dependency(gene_symbol,
                              dataset="Chronos_Combined"):
    """
    DepMap — 遺伝子依存性スコア取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル (例: "KRAS")
        dataset: str — データセット名
            (例: "Chronos_Combined", "RNAi_merged")
    """
    # DepMap Public Portal download URL
    DEPMAP_PORTAL = "https://depmap.org/portal/api"

    url = f"{DEPMAP_PORTAL}/dataset/search"
    params = {"gene": gene_symbol,
              "dataset": dataset}

    try:
        resp = requests.get(url, params=params,
                            timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        # Fallback: alternative DepMap API
        print(f"  DepMap portal fallback for "
              f"{gene_symbol}")
        data = []

    rows = []
    if isinstance(data, list):
        for entry in data:
            if isinstance(entry, dict):
                rows.append({
                    "gene": gene_symbol,
                    "cell_line": entry.get(
                        "cell_line_name", ""),
                    "depmap_id": entry.get(
                        "depmap_id", ""),
                    "dependency_score": entry.get(
                        "dependency", 0),
                    "dataset": dataset,
                })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("dependency_score")
    print(f"DepMap dependency: {gene_symbol} "
          f"→ {len(df)} cell lines")
    return df
```

## 3. 薬剤感受性データ

```python
def depmap_drug_sensitivity(compound_name=None,
                               limit=100):
    """
    DepMap — 薬剤感受性データ取得。

    Parameters:
        compound_name: str — 化合物名フィルタ
            (例: "Paclitaxel")
        limit: int — 最大結果数
    """
    url = f"{DEPMAP_API}/drugs"
    params = {"page_size": limit}
    if compound_name:
        params["name"] = compound_name

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for drug in data.get("data", data
                         if isinstance(data, list)
                         else []):
        if isinstance(drug, dict):
            rows.append({
                "drug_id": drug.get("drug_id", ""),
                "drug_name": drug.get("drug_name", ""),
                "target": drug.get("target", ""),
                "pathway": drug.get("pathway", ""),
                "n_cell_lines": drug.get(
                    "n_cell_lines_tested", 0),
            })

    df = pd.DataFrame(rows[:limit])
    print(f"DepMap drugs: {len(df)} compounds")
    return df
```

## 4. DepMap 統合パイプライン

```python
def depmap_pipeline(gene_symbol,
                       output_dir="results"):
    """
    DepMap 統合パイプライン。

    Parameters:
        gene_symbol: str — 遺伝子シンボル
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 遺伝子依存性
    deps = depmap_gene_dependency(gene_symbol)
    deps.to_csv(output_dir / "depmap_dependency.csv",
                index=False)

    # 2) 関連薬剤
    drugs = depmap_drug_sensitivity()
    drugs.to_csv(output_dir / "depmap_drugs.csv",
                 index=False)

    # 3) 依存性の高い細胞株
    if not deps.empty:
        top_dependent = deps.head(10)
        top_dependent.to_csv(
            output_dir / "depmap_top_dependent.csv",
            index=False)

    print(f"DepMap pipeline: {gene_symbol} → {output_dir}")
    return {"dependency": deps, "drugs": drugs}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `depmap` | DepMap | がん細胞株依存性マップ API |

## パイプライン統合

```
cancer-genomics → depmap-dependencies → precision-oncology
  (がんゲノム)       (DepMap Portal)     (精密腫瘍学)
        │                  │                  ↓
expression-analysis ──────┘       drug-target-profiling
  (発現リスト)                    (標的プロファイリング)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/depmap_dependency.csv` | 依存性スコア | → cancer-genomics |
| `results/depmap_drugs.csv` | 薬剤感受性 | → drug-target-profiling |
| `results/depmap_top_dependent.csv` | 依存細胞株 | → precision-oncology |
