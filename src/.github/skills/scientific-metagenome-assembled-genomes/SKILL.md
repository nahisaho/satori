---
name: scientific-metagenome-assembled-genomes
description: |
  メタゲノムアセンブルゲノム (MAG) 解析スキル。
  MetaBAT2 / CONCOCT / MaxBin2 ビニング・CheckM2 品質評価・
  GTDB-Tk 分類学的分類・dRep 脱重複・Prokka アノテーション・
  MAG アセンブリ品質レポートパイプライン。
  TU 外スキル (CLI ラッパー + Python ライブラリ)。
tu_tools:
  - key: mgnify
    name: MGnify
    description: メタゲノムアセンブル・MAG データ検索
---

# Scientific Metagenome-Assembled Genomes

メタゲノムリードから個別ゲノム (MAG) を再構築する
ビニング・品質評価・分類・アノテーションの
統合パイプラインを提供する。

## When to Use

- メタゲノムショットガンデータから MAG を再構築するとき
- コンティグビニング (MetaBAT2/CONCOCT/MaxBin2) を実行するとき
- CheckM/CheckM2 でゲノム完全性・コンタミネーションを評価するとき
- GTDB-Tk で MAG の分類学的位置づけを行うとき
- dRep で冗長な MAG を脱重複するとき
- Prokka/Bakta で MAG のアノテーションを行うとき

---

## Quick Start

## 1. MetaBAT2 ビニング

```python
import subprocess
import pandas as pd
from pathlib import Path


def run_metabat2(assembly_fasta, bam_file,
                   output_dir="metabat2_bins",
                   min_contig=2500):
    """
    MetaBAT2 — メタゲノムコンティグビニング。

    Parameters:
        assembly_fasta: str — アセンブリ FASTA
        bam_file: str — ソート済み BAM
        output_dir: str — 出力ディレクトリ
        min_contig: int — 最小コンティグ長
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 深度テーブル生成
    depth_file = out / "depth.txt"
    subprocess.run([
        "jgi_summarize_bam_contig_depths",
        "--outputDepth", str(depth_file),
        bam_file
    ], check=True)

    # MetaBAT2 実行
    subprocess.run([
        "metabat2",
        "-i", assembly_fasta,
        "-a", str(depth_file),
        "-o", str(out / "bin"),
        "-m", str(min_contig),
        "--seed", "42",
    ], check=True)

    bins = list(out.glob("bin.*.fa"))
    print(f"MetaBAT2: {len(bins)} bins generated")
    return bins
```

## 2. CheckM2 品質評価

```python
def run_checkm2(bin_dir, output_dir="checkm2_out",
                  threads=8):
    """
    CheckM2 — MAG 品質評価
    (完全性 / コンタミネーション / N50)。

    Parameters:
        bin_dir: str — ビンディレクトリ
        output_dir: str — 出力ディレクトリ
        threads: int — スレッド数
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    subprocess.run([
        "checkm2", "predict",
        "--input", bin_dir,
        "--output-directory", str(out),
        "--threads", str(threads),
        "-x", "fa",
    ], check=True)

    report = out / "quality_report.tsv"
    df = pd.read_csv(report, sep="\t")

    # MIMAG 基準による分類
    df["quality"] = df.apply(
        lambda r: (
            "high" if r["Completeness"] >= 90
            and r["Contamination"] < 5
            else "medium"
            if r["Completeness"] >= 50
            and r["Contamination"] < 10
            else "low"), axis=1)

    n_hq = (df["quality"] == "high").sum()
    n_mq = (df["quality"] == "medium").sum()
    n_lq = (df["quality"] == "low").sum()
    print(f"CheckM2: {n_hq} HQ, {n_mq} MQ, "
          f"{n_lq} LQ MAGs")
    return df


def filter_quality_mags(checkm_df,
                          min_completeness=50,
                          max_contamination=10):
    """
    品質基準によるMAGフィルタリング。

    Parameters:
        checkm_df: pd.DataFrame — CheckM2 結果
        min_completeness: float — 最小完全性 (%)
        max_contamination: float — 最大汚染 (%)
    """
    filtered = checkm_df[
        (checkm_df["Completeness"]
         >= min_completeness)
        & (checkm_df["Contamination"]
           <= max_contamination)
    ].copy()

    print(f"Filter: {len(filtered)}/"
          f"{len(checkm_df)} MAGs passed "
          f"(≥{min_completeness}% comp, "
          f"≤{max_contamination}% contam)")
    return filtered
```

## 3. GTDB-Tk 分類

```python
def run_gtdbtk(bin_dir, output_dir="gtdbtk_out",
                 threads=8):
    """
    GTDB-Tk — ゲノム分類学分類
    (GTDB taxonomy)。

    Parameters:
        bin_dir: str — フィルタ済みビンディレクトリ
        output_dir: str — 出力ディレクトリ
        threads: int — スレッド数
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    subprocess.run([
        "gtdbtk", "classify_wf",
        "--genome_dir", bin_dir,
        "--out_dir", str(out),
        "--cpus", str(threads),
        "-x", "fa",
    ], check=True)

    # 細菌/古細菌分類結果を統合
    results = []
    for domain in ["bac120", "ar53"]:
        tsv = (out / f"gtdbtk.{domain}."
               "summary.tsv")
        if tsv.exists():
            df = pd.read_csv(tsv, sep="\t")
            df["domain_marker"] = domain
            results.append(df)

    if results:
        combined = pd.concat(results,
                             ignore_index=True)
        print(f"GTDB-Tk: {len(combined)} MAGs "
              f"classified")
        return combined

    print("GTDB-Tk: no classification results")
    return pd.DataFrame()
```

## 4. dRep 脱重複

```python
def run_drep(bin_dir, output_dir="drep_out",
               ani_threshold=0.95):
    """
    dRep — MAG 脱重複 (ANI ベース)。

    Parameters:
        bin_dir: str — ビンディレクトリ
        output_dir: str — 出力ディレクトリ
        ani_threshold: float — ANI 閾値
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    subprocess.run([
        "dRep", "dereplicate",
        str(out),
        "-g", f"{bin_dir}/*.fa",
        "-sa", str(ani_threshold),
        "--ignoreGenomeQuality",
    ], check=True)

    derep = list(
        (out / "dereplicated_genomes").glob("*.fa"))
    print(f"dRep: {len(derep)} dereplicated MAGs "
          f"(ANI ≥ {ani_threshold})")
    return derep
```

## 5. MAG パイプライン統合

```python
def mag_pipeline(assembly_fasta, bam_file,
                   output_dir="mag_results",
                   threads=8):
    """
    MAG 統合パイプライン。

    Parameters:
        assembly_fasta: str — メタゲノムアセンブリ
        bam_file: str — ソート済み BAM
        output_dir: str — 出力ルート
        threads: int — スレッド数
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1) ビニング
    bins = run_metabat2(
        assembly_fasta, bam_file,
        str(out / "bins"))

    # 2) 品質評価
    checkm = run_checkm2(
        str(out / "bins"),
        str(out / "checkm2"),
        threads)

    # 3) フィルタリング (MIMAG medium+)
    quality = filter_quality_mags(checkm)

    # 4) GTDB-Tk 分類
    taxonomy = run_gtdbtk(
        str(out / "bins"),
        str(out / "gtdbtk"),
        threads)

    # 5) 脱重複
    derep = run_drep(
        str(out / "bins"),
        str(out / "drep"))

    print(f"MAG pipeline: {len(bins)} bins → "
          f"{len(quality)} QC passed → "
          f"{len(derep)} dereplicated")

    # === パイプライン連携用の構造化出力 ===
    results = Path(output_dir) / "results"
    results.mkdir(parents=True, exist_ok=True)

    # 1) MAG品質サマリーCSV (→ phylogenetics, environmental-ecology)
    quality.to_csv(results / "mag_quality_summary.csv", index=False)
    print(f"  ✔ MAG quality summary: {results / 'mag_quality_summary.csv'}")

    # 2) 分類学サマリーCSV (→ phylogenetics)
    if not taxonomy.empty:
        taxonomy.to_csv(results / "mag_taxonomy.csv", index=False)
        print(f"  ✔ MAG taxonomy: {results / 'mag_taxonomy.csv'}")

    # 3) 代表MAGをFASTAに統合 (→ phylogenetics, annotation)
    representative_fasta = results / "representative_mags.fasta"
    with open(representative_fasta, "w") as f:
        for mag_path in derep:
            mag_name = Path(mag_path).stem
            with open(mag_path) as mag_f:
                for line in mag_f:
                    if line.startswith(">"):
                        f.write(f">{mag_name}_{line[1:]}")
                    else:
                        f.write(line)
    print(f"  ✔ Representative MAGs FASTA: {representative_fasta}")

    # 4) パイプラインサマリーJSON
    import json
    pipeline_summary = {
        "total_bins": len(bins),
        "quality_passed": len(quality),
        "high_quality": int((quality["quality"] == "high").sum()) if "quality" in quality.columns else 0,
        "medium_quality": int((quality["quality"] == "medium").sum()) if "quality" in quality.columns else 0,
        "dereplicated": len(derep),
        "classified": len(taxonomy) if not taxonomy.empty else 0,
    }
    with open(results / "mag_pipeline_summary.json", "w") as f:
        json.dump(pipeline_summary, f, indent=2)
    print(f"  ✔ Pipeline summary: {results / 'mag_pipeline_summary.json'}")

    return {
        "bins": bins,
        "checkm": checkm,
        "quality": quality,
        "taxonomy": taxonomy,
        "dereplicated": derep,
    }
```

---

## パイプライン統合

```
microbiome-metagenomics → metagenome-assembled-genomes → environmental-ecology
  (メタゲノム組成解析)      (MAG 再構築)                   (生態系統合)
        │                        │                            ↓
  long-read-sequencing ─────────┘                    phylogenomics
    (ロングリードアセンブリ)                            (系統解析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `*_bins/bin.*.fa` | ビンゲノム | → dRep, GTDB-Tk |
| `checkm2_out/quality_report.tsv` | 品質レポート | → フィルタリング |
| `gtdbtk_out/*.summary.tsv` | 分類結果 | → phylogenomics |
| `drep_out/dereplicated_genomes/` | 脱重複 MAG | → environmental-ecology |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `mgnify` | MGnify | メタゲノムアセンブル・MAG データ検索 |
