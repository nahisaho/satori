---
name: scientific-rrna-taxonomy
description: |
  rRNA リファレンス・分類学スキル。SILVA SSU/LSU rRNA データベース・
  Greengenes2 系統分類・MGnify メタゲノム解析・QIIME2 分類器・
  scikit-bio 配列解析・系統分類パイプライン。
---

# Scientific rRNA Taxonomy

SILVA / Greengenes2 / MGnify を活用した rRNA リファレンスおよび
分類学的アノテーションパイプラインを提供する。16S/18S/ITS
アンプリコン配列の分類学的帰属と系統解析。

## When to Use

- 16S rRNA アンプリコン配列の分類学的帰属を行うとき
- SILVA/Greengenes2 リファレンスで分類器を訓練するとき
- MGnify からメタゲノム解析結果を取得するとき
- 18S/ITS 真核生物分類を行うとき
- ASV/OTU の分類学的コンセンサスを判定するとき
- QIIME2 カスタム分類器パイプラインを構築するとき

---

## Quick Start

## 1. SILVA rRNA リファレンス取得

```python
import requests
import pandas as pd
from pathlib import Path
from io import StringIO

SILVA_BASE = "https://www.arb-silva.de/api"


def download_silva_reference(version="138.1", subunit="SSU",
                              output_dir="references"):
    """
    SILVA rRNA リファレンス配列 & 分類取得。

    Parameters:
        version: str — SILVA バージョン
        subunit: str — "SSU" (16S/18S) or "LSU" (23S/28S)
        output_dir: str — 出力ディレクトリ
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # SILVA FTP から NR99 配列取得
    base_url = f"https://www.arb-silva.de/fileadmin/silva_databases/release_{version}/Exports"
    fasta_url = f"{base_url}/SILVA_{version}_{subunit}Ref_NR99_tax_silva.fasta.gz"
    tax_url = f"{base_url}/taxonomy/tax_slv_{subunit.lower()}_{version}.txt.gz"

    import urllib.request
    import gzip

    # 配列ダウンロード
    fasta_path = output_dir / f"silva_{version}_{subunit}_NR99.fasta.gz"
    if not fasta_path.exists():
        urllib.request.urlretrieve(fasta_url, str(fasta_path))
        print(f"Downloaded: {fasta_path}")

    # 分類辞書ダウンロード
    tax_path = output_dir / f"silva_{version}_{subunit}_taxonomy.txt.gz"
    if not tax_path.exists():
        urllib.request.urlretrieve(tax_url, str(tax_path))
        print(f"Downloaded: {tax_path}")

    # 配列数カウント
    n_seqs = 0
    with gzip.open(str(fasta_path), "rt") as f:
        for line in f:
            if line.startswith(">"):
                n_seqs += 1

    print(f"SILVA {version} {subunit}: {n_seqs} reference sequences")
    return {"fasta": str(fasta_path), "taxonomy": str(tax_path), "n_seqs": n_seqs}
```

## 2. Greengenes2 分類学取得

```python
def download_greengenes2(version="2024.09", output_dir="references"):
    """
    Greengenes2 分類・系統樹・配列取得。

    Parameters:
        version: str — GG2 バージョン
        output_dir: str — 出力ディレクトリ
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gg2_base = f"https://ftp.microbio.me/greengenes_release/{version}"
    files = {
        "taxonomy": f"{gg2_base}/taxonomy.tsv.gz",
        "backbone": f"{gg2_base}/gg2-backbone.nwk.gz",
        "seqs": f"{gg2_base}/gg2-seqs.fna.gz",
    }

    import urllib.request
    paths = {}
    for name, url in files.items():
        out_path = output_dir / f"gg2_{version}_{name}{Path(url).suffix}{Path(url).suffixes[-1] if len(Path(url).suffixes) > 1 else ''}"
        out_path = output_dir / Path(url).name
        if not out_path.exists():
            try:
                urllib.request.urlretrieve(url, str(out_path))
                print(f"Downloaded: {out_path}")
            except Exception as e:
                print(f"Warning: {name} download failed: {e}")
        paths[name] = str(out_path)

    print(f"Greengenes2 {version}: {len(paths)} files downloaded")
    return paths
```

## 3. MGnify メタゲノム解析結果取得

```python
MGNIFY_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"


def mgnify_study_search(query, biome=None, limit=25):
    """
    MGnify — メタゲノム研究検索。

    Parameters:
        query: str — 検索クエリ
        biome: str — バイオーム (例: "root:Environmental:Aquatic")
        limit: int — 最大取得数

    TU: mgnify
    """
    params = {"search": query, "page_size": limit}
    if biome:
        params["lineage"] = biome

    resp = requests.get(f"{MGNIFY_BASE}/studies", params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for study in data.get("data", []):
        attrs = study.get("attributes", {})
        results.append({
            "study_id": study["id"],
            "name": attrs.get("study-name", ""),
            "abstract": attrs.get("study-abstract", "")[:200],
            "biome": attrs.get("biome-name", ""),
            "samples_count": attrs.get("samples-count", 0),
        })

    df = pd.DataFrame(results)
    print(f"MGnify: '{query}' → {len(df)} studies")
    return df


def mgnify_taxonomy(analysis_id):
    """
    MGnify — 分類学的アノテーション結果取得。

    Parameters:
        analysis_id: str — MGnify 解析 ID

    TU: mgnify
    """
    url = f"{MGNIFY_BASE}/analyses/{analysis_id}/taxonomy/ssu"
    resp = requests.get(url, params={"page_size": 100}, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    taxa = []
    for entry in data.get("data", []):
        attrs = entry.get("attributes", {})
        taxa.append({
            "lineage": attrs.get("lineage", ""),
            "count": attrs.get("count", 0),
            "rank": attrs.get("hierarchy", {}).get("rank", ""),
        })

    df = pd.DataFrame(taxa)
    df = df.sort_values("count", ascending=False)
    print(f"MGnify taxonomy ({analysis_id}): {len(df)} taxa")
    return df
```

## 4. QIIME2 分類器パイプライン

```python
def qiime2_classify_sklearn(sequences_path, reference_seqs, reference_tax,
                             classifier_output="classifier.qza"):
    """
    QIIME2 scikit-learn 分類器訓練 & 分類。

    Parameters:
        sequences_path: str — 入力配列 (FASTA or QZA)
        reference_seqs: str — リファレンス配列パス
        reference_tax: str — リファレンス分類パス
        classifier_output: str — 分類器出力パス
    """
    import subprocess

    # 1) リファレンスインポート
    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Sequence]",
        "--input-path", reference_seqs,
        "--output-path", "ref-seqs.qza",
    ], check=True)

    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Taxonomy]",
        "--input-format", "HeaderlessTSVTaxonomyFormat",
        "--input-path", reference_tax,
        "--output-path", "ref-taxonomy.qza",
    ], check=True)

    # 2) 分類器訓練
    subprocess.run([
        "qiime", "feature-classifier", "fit-classifier-naive-bayes",
        "--i-reference-reads", "ref-seqs.qza",
        "--i-reference-taxonomy", "ref-taxonomy.qza",
        "--o-classifier", classifier_output,
    ], check=True)

    # 3) 分類実行
    subprocess.run([
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-classifier", classifier_output,
        "--i-reads", sequences_path,
        "--o-classification", "taxonomy.qza",
    ], check=True)

    print(f"QIIME2 classification complete: {classifier_output}")
    return "taxonomy.qza"
```

## 5. 分類学的コンセンサス解析

```python
import numpy as np


def taxonomy_consensus(classifications, confidence_threshold=0.8):
    """
    複数分類器のコンセンサス分類。

    Parameters:
        classifications: dict — {method: DataFrame(feature_id, taxon, confidence)}
        confidence_threshold: float — 信頼度閾値
    """
    all_features = set()
    for method_df in classifications.values():
        all_features.update(method_df["feature_id"].tolist())

    consensus = []
    for feat_id in all_features:
        taxa = {}
        for method, df in classifications.items():
            row = df[df["feature_id"] == feat_id]
            if len(row) > 0:
                taxa[method] = {
                    "taxon": row.iloc[0]["taxon"],
                    "confidence": row.iloc[0].get("confidence", 1.0),
                }

        # ランクごとのコンセンサス
        if taxa:
            lineages = [t["taxon"] for t in taxa.values()]
            confidences = [t["confidence"] for t in taxa.values()]

            # 分割してランク比較
            split_lineages = [l.split(";") for l in lineages]
            max_depth = max(len(sl) for sl in split_lineages)
            consensus_lineage = []

            for rank_idx in range(max_depth):
                rank_taxa = [sl[rank_idx] for sl in split_lineages
                             if rank_idx < len(sl)]
                most_common = max(set(rank_taxa), key=rank_taxa.count)
                agreement = rank_taxa.count(most_common) / len(rank_taxa)

                if agreement >= confidence_threshold:
                    consensus_lineage.append(most_common)
                else:
                    break

            consensus.append({
                "feature_id": feat_id,
                "consensus_taxon": ";".join(consensus_lineage),
                "depth": len(consensus_lineage),
                "methods_agree": len(taxa),
                "mean_confidence": np.mean(confidences),
            })

    df = pd.DataFrame(consensus)
    print(f"Consensus: {len(df)} features, "
          f"mean depth={df['depth'].mean():.1f}")
    return df
```

## 6. rRNA 分類統合パイプライン

```python
def rrna_taxonomy_pipeline(input_fasta, output_dir="results",
                            silva_version="138.1", use_greengenes=True):
    """
    SILVA + Greengenes2 統合 rRNA 分類パイプライン。

    Parameters:
        input_fasta: str — 入力 16S rRNA 配列
        output_dir: str — 出力ディレクトリ
        silva_version: str — SILVA バージョン
        use_greengenes: bool — GG2 も併用
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) リファレンスダウンロード
    silva_ref = download_silva_reference(
        version=silva_version, output_dir=str(output_dir / "refs")
    )

    refs = {"silva": silva_ref}
    if use_greengenes:
        gg2_ref = download_greengenes2(
            output_dir=str(output_dir / "refs")
        )
        refs["greengenes2"] = gg2_ref

    # 2) QIIME2 分類 (SILVA)
    silva_taxonomy = qiime2_classify_sklearn(
        input_fasta,
        silva_ref["fasta"],
        silva_ref["taxonomy"],
        classifier_output=str(output_dir / "silva_classifier.qza"),
    )

    # 3) MGnify 比較参照
    # (解析済みデータが MGnify にあれば取得)

    print(f"Pipeline complete: {len(refs)} references used")
    return refs
```

---

## パイプライン統合

```
microbiome-metagenomics → rrna-taxonomy → phylogenetics
  (DADA2 ASV パイプライン)  (SILVA/GG2 分類)  (ETE3 系統樹)
        │                       │                   ↓
environmental-ecology ─────────┘           population-genetics
  (α/β 多様性)              │              (Fst/ADMIXTURE)
                             ↓
                    pathway-enrichment
                    (微生物機能濃縮)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/taxonomy.csv` | 分類学的帰属結果 | → microbiome-metagenomics |
| `results/consensus.csv` | コンセンサス分類 | → phylogenetics |
| `results/refs/` | SILVA/GG2 リファレンス | — |
| `results/mgnify_taxonomy.csv` | MGnify 分類結果 | → environmental-ecology |

### 利用可能ツール (ToolUniverse SMCP)

| Config Key | ツール数 | 主要機能 |
|-----------|---------|---------|
| `mgnify` | 3+ | メタゲノム研究検索・分類学的プロファイル・機能アノテーション |
