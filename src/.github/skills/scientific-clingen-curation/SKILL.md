---
name: scientific-clingen-curation
description: |
  ClinGen 臨床ゲノム資源キュレーションスキル。ClinGen API に
  よる遺伝子-疾患バリディティ、臨床アクショナビリティ、
  投与量感受性、バリアントレベルエビデンス評価パイプライン。
  ToolUniverse 連携: clingen。
tu_tools:
  - key: clingen
    name: ClinGen
    description: ClinGen 臨床ゲノムリソース キュレーションデータ
---

# Scientific ClinGen Curation

ClinGen (Clinical Genome Resource) API を活用した
遺伝子-疾患バリディティ分類・臨床アクショナビリティ
スコアリング・投与量感受性評価・バリアントキュレーション
パイプラインを提供する。

## When to Use

- 遺伝子-疾患関連のエビデンスレベルを評価するとき
- 臨床アクショナビリティ (介入可能性) を判定するとき
- ハプロ不全/トリプロ感受性を評価するとき
- ClinGen キュレーション済みバリアント分類を取得するとき
- ACMG ガイドラインに基づくバリアント解釈を行うとき

---

## Quick Start

## 1. 遺伝子-疾患バリディティ

```python
import requests
import pandas as pd

CLINGEN_BASE = "https://search.clinicalgenome.org/kb"


def clingen_gene_validity(gene_symbol):
    """
    ClinGen — 遺伝子-疾患バリディティ分類取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル (例: "BRCA1")
    """
    url = (f"{CLINGEN_BASE}/gene-validity/"
           f"?search={gene_symbol}&format=json")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = data if isinstance(data, list) else \
        data.get("results", [])

    rows = []
    for item in results:
        rows.append({
            "gene": item.get("gene", {}).get(
                "symbol", gene_symbol),
            "disease": item.get("disease", {}).get(
                "label", ""),
            "classification": item.get(
                "classification", ""),
            "moi": item.get("moi", ""),
            "sop": item.get("sopVersion", ""),
        })

    df = pd.DataFrame(rows)
    print(f"ClinGen validity: {gene_symbol} → "
          f"{len(df)} gene-disease pairs")
    return df


def clingen_gene_validity_batch(gene_symbols):
    """
    ClinGen — 複数遺伝子バリディティバッチ取得。

    Parameters:
        gene_symbols: list[str] — 遺伝子シンボルリスト
    """
    all_results = []
    for sym in gene_symbols:
        df = clingen_gene_validity(sym)
        if not df.empty:
            all_results.append(df)
    if all_results:
        combined = pd.concat(all_results,
                             ignore_index=True)
        cls_dist = combined["classification"].value_counts()
        print(f"Validity distribution: "
              f"{cls_dist.to_dict()}")
        return combined
    return pd.DataFrame()
```

## 2. 投与量感受性

```python
def clingen_dosage_sensitivity(gene_symbol):
    """
    ClinGen — 投与量感受性 (haplo/triplo) 評価取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル
    """
    url = (f"{CLINGEN_BASE}/gene-dosage/"
           f"?search={gene_symbol}&format=json")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = data if isinstance(data, list) else \
        data.get("results", [])

    rows = []
    for item in results:
        rows.append({
            "gene": item.get("gene", {}).get(
                "symbol", gene_symbol),
            "haplo_score": item.get(
                "haploinsufficiency", {}).get(
                    "score", ""),
            "haplo_label": item.get(
                "haploinsufficiency", {}).get(
                    "label", ""),
            "triplo_score": item.get(
                "triplosensitivity", {}).get(
                    "score", ""),
            "triplo_label": item.get(
                "triplosensitivity", {}).get(
                    "label", ""),
        })

    df = pd.DataFrame(rows)
    print(f"ClinGen dosage: {gene_symbol} → "
          f"{len(df)} entries")
    return df
```

## 3. 臨床アクショナビリティ

```python
def clingen_actionability(gene_symbol):
    """
    ClinGen — 臨床アクショナビリティスコア取得。

    Parameters:
        gene_symbol: str — 遺伝子シンボル
    """
    url = (f"{CLINGEN_BASE}/actionability/"
           f"?search={gene_symbol}&format=json")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = data if isinstance(data, list) else \
        data.get("results", [])

    rows = []
    for item in results:
        rows.append({
            "gene": item.get("gene", {}).get(
                "symbol", gene_symbol),
            "disease": item.get("disease", {}).get(
                "label", ""),
            "classification": item.get(
                "classification", ""),
            "date": item.get("date", ""),
        })

    df = pd.DataFrame(rows)
    print(f"ClinGen actionability: {gene_symbol} → "
          f"{len(df)} entries")
    return df
```

## 4. ClinGen 統合パイプライン

```python
def clingen_pipeline(gene_symbols,
                       output_dir="results"):
    """
    ClinGen 統合キュレーションパイプライン。

    Parameters:
        gene_symbols: list[str] — 遺伝子シンボルリスト
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) Gene-disease validity
    validity_df = clingen_gene_validity_batch(
        gene_symbols)
    if not validity_df.empty:
        validity_df.to_csv(
            output_dir / "clingen_validity.csv",
            index=False)

    # 2) Dosage sensitivity
    dosage_results = []
    for sym in gene_symbols:
        dos = clingen_dosage_sensitivity(sym)
        if not dos.empty:
            dosage_results.append(dos)
    if dosage_results:
        dosage_df = pd.concat(dosage_results,
                              ignore_index=True)
        dosage_df.to_csv(
            output_dir / "clingen_dosage.csv",
            index=False)

    # 3) Actionability
    action_results = []
    for sym in gene_symbols:
        act = clingen_actionability(sym)
        if not act.empty:
            action_results.append(act)
    if action_results:
        action_df = pd.concat(action_results,
                              ignore_index=True)
        action_df.to_csv(
            output_dir / "clingen_actionability.csv",
            index=False)

    print(f"ClinGen pipeline → {output_dir}")
    return {"validity": validity_df}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `clingen` | ClinGen | ClinGen 臨床ゲノムリソース キュレーションデータ |

## パイプライン統合

```
variant-interpretation → clingen-curation → clinical-decision-support
  (ClinVar/ACMG)        (GDV/DOS/ACT)        (臨床判断支援)
         │                    │                      ↓
  variant-effect-prediction ─┘          pharmacogenomics
    (SpliceAI/CADD)                       (PGx 処方)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/clingen_validity.csv` | 遺伝子-疾患バリディティ | → genetic-counseling |
| `results/clingen_dosage.csv` | 投与量感受性 | → cnv-analysis |
| `results/clingen_actionability.csv` | 臨床介入可能性 | → precision-medicine |
