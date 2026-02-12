---
name: scientific-pharos-targets
description: |
  Pharos/TCRD ターゲットプロファイリングスキル。Illuminating
  the Druggable Genome (IDG) Pharos GraphQL API による
  ターゲット開発レベル (TDL) 分類・疾患関連・リガンド検索・
  既存薬候補スクリーニングパイプライン。
  ToolUniverse 連携: pharos。
tu_tools:
  - key: pharos
    name: Pharos
    description: IDG Pharos/TCRD ターゲットナレッジベース
---

# Scientific Pharos Targets

IDG (Illuminating the Druggable Genome) Pharos GraphQL API
を活用したターゲット TDL 分類・疾患-ターゲット関連・
リガンドアクティビティ検索・既存薬候補パイプラインを提供する。

## When to Use

- タンパク質標的の開発レベル (Tclin/Tchem/Tbio/Tdark) を確認するとき
- IDG ダークターゲットの知識ギャップを調査するとき
- ターゲットに関連する疾患・リガンド・PPI を検索するとき
- Tchem/Tbio ターゲットから創薬候補を探索するとき
- ターゲット間のドラッガビリティ比較を行うとき

---

## Quick Start

## 1. ターゲット TDL 分類検索

```python
import requests
import pandas as pd

PHAROS_API = "https://pharos-api.ncats.io/graphql"


def pharos_target(gene_symbol):
    """
    Pharos — ターゲット詳細情報取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル (例: "ACE2")
    """
    query = """
    query targetDetails($sym: String!) {
      target(q: {sym: $sym}) {
        name
        sym
        tdl
        fam
        novelty
        description
        uniprotId: accession
        diseaseCount
        ligandCount
        ppiCount
      }
    }
    """
    resp = requests.post(
        PHAROS_API,
        json={"query": query,
              "variables": {"sym": gene_symbol}},
        timeout=30)
    resp.raise_for_status()
    data = resp.json().get("data", {}).get("target", {})

    if not data:
        print(f"Pharos: {gene_symbol} not found")
        return {}

    print(f"Pharos: {gene_symbol} → TDL={data['tdl']}, "
          f"fam={data.get('fam', 'N/A')}, "
          f"diseases={data.get('diseaseCount', 0)}, "
          f"ligands={data.get('ligandCount', 0)}")
    return data


def pharos_target_list(gene_symbols):
    """
    Pharos — 複数ターゲット TDL バッチ取得。

    Parameters:
        gene_symbols: list[str] — 遺伝子シンボルリスト
    """
    results = []
    for sym in gene_symbols:
        data = pharos_target(sym)
        if data:
            results.append(data)

    df = pd.DataFrame(results)
    if not df.empty:
        tdl_counts = df["tdl"].value_counts()
        print(f"TDL distribution: "
              f"{tdl_counts.to_dict()}")
    return df
```

## 2. 疾患-ターゲット関連

```python
def pharos_target_diseases(gene_symbol, top_n=20):
    """
    Pharos — ターゲットの関連疾患取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル
        top_n: int — 上位件数
    """
    query = """
    query targetDiseases($sym: String!, $top: Int!) {
      target(q: {sym: $sym}) {
        sym
        diseases(top: $top) {
          name
          associationCount
          datasource_count
        }
      }
    }
    """
    resp = requests.post(
        PHAROS_API,
        json={"query": query,
              "variables": {"sym": gene_symbol,
                            "top": top_n}},
        timeout=30)
    resp.raise_for_status()
    target = resp.json().get("data", {}).get(
        "target", {})
    diseases = target.get("diseases", [])

    df = pd.DataFrame(diseases)
    print(f"Pharos diseases: {gene_symbol} → "
          f"{len(df)} associations")
    return df
```

## 3. リガンド・アクティビティ検索

```python
def pharos_target_ligands(gene_symbol, top_n=20):
    """
    Pharos — ターゲットのリガンド/活性化合物取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル
        top_n: int — 上位件数
    """
    query = """
    query targetLigands($sym: String!, $top: Int!) {
      target(q: {sym: $sym}) {
        sym
        ligands(top: $top) {
          name
          smiles
          isdrug
          actcnt
          activities {
            type
            value
            moa
          }
        }
      }
    }
    """
    resp = requests.post(
        PHAROS_API,
        json={"query": query,
              "variables": {"sym": gene_symbol,
                            "top": top_n}},
        timeout=30)
    resp.raise_for_status()
    target = resp.json().get("data", {}).get(
        "target", {})
    ligands = target.get("ligands", [])

    rows = []
    for lig in ligands:
        rows.append({
            "name": lig.get("name", ""),
            "smiles": lig.get("smiles", ""),
            "is_drug": lig.get("isdrug", False),
            "activity_count": lig.get("actcnt", 0),
        })

    df = pd.DataFrame(rows)
    n_drugs = df["is_drug"].sum() if not df.empty else 0
    print(f"Pharos ligands: {gene_symbol} → "
          f"{len(df)} compounds ({n_drugs} drugs)")
    return df
```

## 4. Pharos 統合パイプライン

```python
def pharos_pipeline(gene_symbols,
                      output_dir="results"):
    """
    Pharos 統合パイプライン。

    Parameters:
        gene_symbols: list[str] — 遺伝子シンボルリスト
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) TDL 分類
    targets = pharos_target_list(gene_symbols)
    targets.to_csv(output_dir / "pharos_targets.csv",
                   index=False)

    # 2) 疾患関連 (上位ターゲット)
    all_diseases = []
    for sym in gene_symbols[:10]:
        diseases = pharos_target_diseases(sym)
        diseases["target"] = sym
        all_diseases.append(diseases)
    if all_diseases:
        disease_df = pd.concat(all_diseases,
                               ignore_index=True)
        disease_df.to_csv(
            output_dir / "pharos_diseases.csv",
            index=False)

    # 3) リガンド
    all_ligands = []
    for sym in gene_symbols[:10]:
        ligands = pharos_target_ligands(sym)
        ligands["target"] = sym
        all_ligands.append(ligands)
    if all_ligands:
        ligand_df = pd.concat(all_ligands,
                              ignore_index=True)
        ligand_df.to_csv(
            output_dir / "pharos_ligands.csv",
            index=False)

    print(f"Pharos pipeline → {output_dir}")
    return {"targets": targets}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `pharos` | Pharos | IDG Pharos/TCRD ターゲットナレッジベース |

## パイプライン統合

```
drug-target-profiling → pharos-targets → drug-repurposing
  (DGIdb/標的同定)      (TDL/IDG)        (リポジショニング)
         │                   │                  ↓
  pharmacology-targets ─────┘      compound-screening
    (BindingDB/GtoPdb)             (ZINC ライブラリ)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/pharos_targets.csv` | TDL 分類結果 | → drug-target-profiling |
| `results/pharos_diseases.csv` | 疾患関連 | → disease-research |
| `results/pharos_ligands.csv` | リガンド/薬剤 | → drug-repurposing |
