---
name: scientific-glycomics
description: |
  糖鎖構造解析スキル。GlyConnect / GlyGen / GlyCosmos
  糖鎖データベース統合検索・糖鎖構造描画・糖タンパク質
  グリコシル化部位予測・レクチンバインディング・
  糖鎖マスフラグメンテーション解析パイプライン。
  TU 外スキル (直接 REST API + Python ライブラリ)。
tu_tools:
  - key: glygen
    name: GlyGen
    description: 糖鎖構造・機能データベース検索
---

# Scientific Glycomics

GlyConnect / GlyGen / GlyCosmos 糖鎖データベースを統合した
糖鎖構造解析・糖タンパク質グリコサイト予測・レクチン特異性・
糖鎖 MS フラグメンテーション解析パイプラインを提供する。

## When to Use

- 糖鎖構造を GlyTouCan ID から検索・描画するとき
- タンパク質のグリコシル化部位を予測・マッピングするとき
- GlyGen/GlyConnect で糖鎖-タンパク質関連を検索するとき
- 糖鎖マススペクトルのフラグメンテーション解析を行うとき
- レクチン-糖鎖結合特異性を調査するとき

---

## Quick Start

## 1. GlyGen 糖鎖検索

```python
import requests
import pandas as pd

GLYGEN_API = "https://api.glygen.org"


def glygen_glycan_search(glycan_type=None,
                           mass_range=None):
    """
    GlyGen — 糖鎖検索。

    Parameters:
        glycan_type: str | None — 糖鎖タイプ
            ("N-linked", "O-linked", "GAG" 等)
        mass_range: tuple | None — (min_mass, max_mass)
    """
    query = {}
    if glycan_type:
        query["glycan_type"] = glycan_type
    if mass_range:
        query["mass"] = {
            "min": mass_range[0],
            "max": mass_range[1]}

    url = f"{GLYGEN_API}/glycan/search"
    resp = requests.post(url, json=query, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = data.get("results", [])
    rows = []
    for r in results:
        rows.append({
            "glytoucan_ac": r.get("glytoucan_ac", ""),
            "mass": r.get("mass", 0),
            "glycan_type": r.get("glycan_type", ""),
            "composition": r.get(
                "composition", ""),
        })

    df = pd.DataFrame(rows)
    print(f"GlyGen search: {len(df)} glycans found")
    return df


def glygen_glycan_detail(glytoucan_ac):
    """
    GlyGen — 糖鎖詳細情報取得。

    Parameters:
        glytoucan_ac: str — GlyTouCan アクセッション
    """
    url = f"{GLYGEN_API}/glycan/detail/{glytoucan_ac}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    info = {
        "glytoucan_ac": data.get("glytoucan_ac", ""),
        "mass": data.get("mass", 0),
        "glycan_type": data.get("glycan_type", ""),
        "iupac": data.get("iupac", ""),
        "glycoct": data.get("glycoct", ""),
        "species": [s.get("name", "")
                    for s in data.get("species", [])],
        "proteins": len(data.get("glycoprotein", [])),
    }

    print(f"GlyGen: {glytoucan_ac} → "
          f"type={info['glycan_type']}, "
          f"mass={info['mass']:.1f}, "
          f"proteins={info['proteins']}")
    return info
```

## 2. 糖タンパク質グリコサイト検索

```python
def glygen_protein_glycosylation(uniprot_ac):
    """
    GlyGen — タンパク質グリコシル化部位取得。

    Parameters:
        uniprot_ac: str — UniProt アクセッション
    """
    url = f"{GLYGEN_API}/protein/detail/{uniprot_ac}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    sites = data.get("glycosylation", [])
    rows = []
    for site in sites:
        rows.append({
            "position": site.get("position", 0),
            "type": site.get("type", ""),
            "glytoucan_ac": site.get(
                "glytoucan_ac", ""),
            "residue": site.get("residue", ""),
            "evidence": site.get("evidence", ""),
        })

    df = pd.DataFrame(rows)
    print(f"GlyGen glycosites: {uniprot_ac} → "
          f"{len(df)} sites")
    return df
```

## 3. 糖鎖 MS フラグメンテーション

```python
import numpy as np


def glycan_fragmentation(composition,
                           ion_type="[M+Na]+"):
    """
    糖鎖 MS フラグメンテーション予測。

    Parameters:
        composition: dict — 糖鎖組成
            例: {"Hex": 5, "HexNAc": 4, "Fuc": 1,
                 "NeuAc": 2}
        ion_type: str — イオン種
    """
    monosaccharide_mass = {
        "Hex": 162.0528,
        "HexNAc": 203.0794,
        "Fuc": 146.0579,
        "NeuAc": 291.0954,
        "NeuGc": 307.0903,
        "Pent": 132.0423,
    }

    adducts = {
        "[M+Na]+": 22.9892,
        "[M+H]+": 1.0073,
        "[M+K]+": 38.9632,
        "[M-H]-": -1.0073,
    }

    total_mass = 18.0106  # water
    for sugar, count in composition.items():
        if sugar in monosaccharide_mass:
            total_mass += (monosaccharide_mass[sugar]
                          * count)

    adduct = adducts.get(ion_type, 22.9892)
    precursor_mz = total_mass + adduct

    # Y-type fragments (reducing end)
    fragments = []
    for sugar, count in composition.items():
        if sugar not in monosaccharide_mass:
            continue
        for i in range(1, count + 1):
            loss = monosaccharide_mass[sugar] * i
            frag_mz = precursor_mz - loss
            fragments.append({
                "type": f"Y (loss {i}x{sugar})",
                "mz": round(frag_mz, 4),
                "loss": round(loss, 4),
            })

    df = pd.DataFrame(fragments).sort_values(
        "mz", ascending=False)
    print(f"Glycan fragmentation: "
          f"precursor={precursor_mz:.4f}, "
          f"{len(df)} fragments")
    return df
```

## 4. 糖鎖解析統合パイプライン

```python
def glycomics_pipeline(uniprot_ids,
                         output_dir="results"):
    """
    糖鎖解析統合パイプライン。

    Parameters:
        uniprot_ids: list[str] — UniProt ID リスト
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) グリコサイトマッピング
    all_sites = []
    for uid in uniprot_ids:
        sites = glygen_protein_glycosylation(uid)
        sites["protein"] = uid
        all_sites.append(sites)
    if all_sites:
        sites_df = pd.concat(all_sites,
                             ignore_index=True)
        sites_df.to_csv(
            output_dir / "glycosites.csv",
            index=False)

    # 2) 糖鎖詳細取得
    unique_glycans = set()
    for df in all_sites:
        if not df.empty:
            unique_glycans.update(
                df["glytoucan_ac"].dropna().unique())

    glycan_details = []
    for gac in list(unique_glycans)[:50]:
        if gac:
            detail = glygen_glycan_detail(gac)
            if detail:
                glycan_details.append(detail)
    if glycan_details:
        gdf = pd.DataFrame(glycan_details)
        gdf.to_csv(
            output_dir / "glycan_details.csv",
            index=False)

    print(f"Glycomics pipeline → {output_dir}")
    return {"sites": sites_df if all_sites else
            pd.DataFrame()}
```

---

## パイプライン統合

```
proteomics-mass-spectrometry → glycomics → pathway-enrichment
  (LC-MS/MS PTM 同定)          (糖鎖構造)   (糖鎖パスウェイ)
         │                        │               ↓
  protein-structure-analysis ────┘    immunoinformatics
    (糖鎖結合サイト構造)               (抗体グリコシル化)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/glycosites.csv` | グリコシル化部位 | → protein-structure-analysis |
| `results/glycan_details.csv` | 糖鎖詳細 | → pathway-enrichment |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `glygen` | GlyGen | 糖鎖構造・機能データベース検索 |
