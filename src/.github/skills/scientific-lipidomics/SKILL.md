---
name: scientific-lipidomics
description: |
  リピドミクス解析スキル。LipidMAPS / SwissLipids / LION
  脂質データベース統合検索・脂質サブクラス分類・
  脂質 MS/MS スペクトル同定・脂質パスウェイエンリッチメント・
  脂質プロファイリングパイプライン。
  TU 外スキル (直接 REST API + Python ライブラリ)。
---

# Scientific Lipidomics

LipidMAPS / SwissLipids / LION 脂質データベースを統合した
脂質構造検索・サブクラス分類・MS/MS 同定・
脂質パスウェイエンリッチメント解析パイプラインを提供する。

## When to Use

- LC-MS/MS リピドミクスデータの脂質同定を行うとき
- LipidMAPS で脂質構造・サブクラスを検索するとき
- 脂質プロファイルの差次解析 (fold change/p-value) を行うとき
- LION エンリッチメントで脂質機能解析を行うとき
- 脂質パスウェイ (スフィンゴ脂質/リン脂質代謝) を可視化するとき

---

## Quick Start

## 1. LipidMAPS 構造検索

```python
import requests
import pandas as pd

LIPIDMAPS_API = "https://www.lipidmaps.org/rest"


def lipidmaps_search(name=None, formula=None,
                       mass=None, tolerance=0.01):
    """
    LipidMAPS — 脂質構造検索。

    Parameters:
        name: str | None — 脂質名 (部分一致)
        formula: str | None — 分子式
        mass: float | None — 精密質量
        tolerance: float — 質量誤差 (Da)
    """
    if mass is not None:
        url = (f"{LIPIDMAPS_API}/compound/lm_id/"
               f"mass/{mass}/{tolerance}")
    elif name:
        url = (f"{LIPIDMAPS_API}/compound/lm_id/"
               f"name/{name}")
    elif formula:
        url = (f"{LIPIDMAPS_API}/compound/lm_id/"
               f"formula/{formula}")
    else:
        print("LipidMAPS: provide name, formula, "
              "or mass")
        return pd.DataFrame()

    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    if isinstance(data, dict):
        data = [data]

    rows = []
    for item in data:
        rows.append({
            "lm_id": item.get("lm_id", ""),
            "name": item.get("name", ""),
            "sys_name": item.get(
                "systematic_name", ""),
            "formula": item.get("formula", ""),
            "mass": item.get("mass", 0),
            "main_class": item.get(
                "main_class", ""),
            "sub_class": item.get("sub_class", ""),
        })

    df = pd.DataFrame(rows)
    print(f"LipidMAPS: {len(df)} lipids found")
    return df


def lipidmaps_classification(lm_id):
    """
    LipidMAPS — 脂質分類階層取得。

    Parameters:
        lm_id: str — LipidMAPS ID (例: "LMFA01010001")
    """
    url = (f"{LIPIDMAPS_API}/compound/"
           f"lm_id/{lm_id}/all")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    classification = {
        "lm_id": data.get("lm_id", ""),
        "category": data.get("core", ""),
        "main_class": data.get("main_class", ""),
        "sub_class": data.get("sub_class", ""),
        "class_level4": data.get(
            "class_level4", ""),
        "smiles": data.get("smiles", ""),
        "inchi_key": data.get("inchi_key", ""),
    }

    print(f"LipidMAPS: {lm_id} → "
          f"{classification['main_class']} / "
          f"{classification['sub_class']}")
    return classification
```

## 2. 脂質差次解析

```python
import numpy as np
from scipy import stats


def lipid_differential_analysis(data, groups,
                                  fdr_threshold=0.05):
    """
    脂質差次解析 (Fold Change + t-test)。

    Parameters:
        data: pd.DataFrame — 脂質濃度行列
            (行=サンプル, 列=脂質)
        groups: list[int] — グループラベル (0 or 1)
        fdr_threshold: float — FDR 閾値
    """
    from statsmodels.stats.multitest import (
        multipletests)

    groups = np.array(groups)
    g0 = data[groups == 0]
    g1 = data[groups == 1]

    results = []
    for lipid in data.columns:
        mean0 = g0[lipid].mean()
        mean1 = g1[lipid].mean()
        fc = (mean1 / mean0 if mean0 > 0
              else np.inf)
        log2fc = np.log2(fc) if fc > 0 else 0
        _, pval = stats.ttest_ind(
            g0[lipid], g1[lipid])
        results.append({
            "lipid": lipid,
            "mean_ctrl": round(mean0, 4),
            "mean_case": round(mean1, 4),
            "fold_change": round(fc, 4),
            "log2FC": round(log2fc, 4),
            "pvalue": pval,
        })

    df = pd.DataFrame(results)
    _, fdr, _, _ = multipletests(
        df["pvalue"], method="fdr_bh")
    df["fdr"] = fdr
    df["significant"] = df["fdr"] < fdr_threshold

    n_sig = df["significant"].sum()
    n_up = ((df["significant"]) &
            (df["log2FC"] > 0)).sum()
    n_down = ((df["significant"]) &
              (df["log2FC"] < 0)).sum()
    print(f"Lipid DA: {n_sig} significant "
          f"({n_up} up, {n_down} down)")
    return df.sort_values("pvalue")
```

## 3. 脂質サブクラスエンリッチメント

```python
def lipid_subclass_enrichment(
        sig_lipids, all_lipids, class_map):
    """
    脂質サブクラスエンリッチメント (Fisher exact)。

    Parameters:
        sig_lipids: list[str] — 有意脂質リスト
        all_lipids: list[str] — 全脂質リスト
        class_map: dict — {lipid: subclass} マッピング
    """
    from scipy.stats import fisher_exact

    sig_set = set(sig_lipids)
    all_set = set(all_lipids)

    # サブクラス別集計
    subclasses = set(class_map.values())
    results = []
    for sc in subclasses:
        sc_all = {l for l, c in class_map.items()
                  if c == sc and l in all_set}
        sc_sig = sc_all & sig_set
        a = len(sc_sig)
        b = len(sig_set) - a
        c = len(sc_all) - a
        d = len(all_set) - a - b - c
        if a == 0:
            continue
        _, pval = fisher_exact(
            [[a, b], [c, d]],
            alternative="greater")
        results.append({
            "subclass": sc,
            "sig_in_class": a,
            "total_in_class": len(sc_all),
            "pvalue": pval,
            "ratio": round(a / len(sc_all), 3),
        })

    df = pd.DataFrame(results).sort_values("pvalue")
    print(f"Subclass enrichment: "
          f"{(df['pvalue'] < 0.05).sum()} "
          f"significant subclasses")
    return df
```

## 4. リピドミクス統合パイプライン

```python
def lipidomics_pipeline(data, groups,
                          output_dir="results"):
    """
    リピドミクス統合パイプライン。

    Parameters:
        data: pd.DataFrame — 脂質濃度行列
        groups: list[int] — グループラベル
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 差次解析
    da = lipid_differential_analysis(data, groups)
    da.to_csv(output_dir / "lipid_da.csv",
              index=False)

    # 2) LipidMAPS アノテーション
    annotations = []
    for lipid in data.columns[:30]:
        result = lipidmaps_search(name=lipid)
        if not result.empty:
            row = result.iloc[0].to_dict()
            row["query"] = lipid
            annotations.append(row)
    if annotations:
        ann_df = pd.DataFrame(annotations)
        ann_df.to_csv(
            output_dir / "lipid_annotations.csv",
            index=False)

    print(f"Lipidomics pipeline → {output_dir}")
    return {"da": da}
```

---

## パイプライン統合

```
metabolomics → lipidomics → pathway-enrichment
  (LC-MS 全代謝物)  (脂質特化)   (脂質代謝パスウェイ)
       │                │               ↓
  metabolomics-network ─┘     multi-omics
    (代謝物相関)               (オミクス統合)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/lipid_da.csv` | 差次脂質 | → biomarker-discovery |
| `results/lipid_annotations.csv` | LipidMAPS 注釈 | → pathway-enrichment |
