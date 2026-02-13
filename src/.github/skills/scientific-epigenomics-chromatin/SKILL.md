---
name: scientific-epigenomics-chromatin
description: |
  エピゲノミクス・クロマチン生物学解析スキル。ChIP-seq ピーク呼び出し (MACS2/MACS3)、
  ATAC-seq ヌクレオソームフリー領域検出、DNA メチル化パターン解析 (WGBS/RRBS)、
  ヒストン修飾クロマチン状態モデリング (ChromHMM)、Hi-C 接触マップ・TAD 検出、
  転写因子結合サイト予測 (モチーフ濃縮)、差次結合解析 (DiffBind) を統合した
  計算エピゲノミクスパイプライン。ChIP-Atlas 43 万+実験との連携対応。
  ToolUniverse 連携: chipatlas。
tu_tools:
  - key: chipatlas
    name: ChIP-Atlas
    description: ChIP-Atlas エピゲノミクスエンリッチメント解析 (43万+実験)
---

# Scientific Epigenomics & Chromatin Biology

ChIP-seq・ATAC-seq・バイサルファイトシーケンシング・Hi-C データを対象に、
ピーク呼び出し→差次結合解析→クロマチン状態注釈→3D ゲノム構造解析の
統合エピゲノミクスパイプラインを提供する。

## When to Use

- ChIP-seq データからヒストン修飾・転写因子結合部位を同定するとき
- ATAC-seq でクロマチンアクセシビリティを評価するとき
- DNA メチル化（WGBS/RRBS）パターンを解析するとき
- Hi-C データから TAD/ループ・3D ゲノム構造を推定するとき
- 複数エピゲノムマークを統合してクロマチン状態を分類するとき

---

## Quick Start

## 1. ChIP-seq ピーク呼び出し (MACS2/MACS3)

```python
import subprocess
import pandas as pd
import numpy as np


def chipseq_peak_calling(treatment_bam, control_bam, genome_size="hs",
                         outdir="results/chipseq", name="sample",
                         peak_type="narrow", qvalue=0.05):
    """
    MACS2/MACS3 による ChIP-seq ピーク呼び出し。

    Parameters:
        treatment_bam: 処理群 BAM ファイル
        control_bam: コントロール BAM ファイル (Input/IgG)
        genome_size: 有効ゲノムサイズ (hs/mm/ce/dm or int)
        peak_type: "narrow" (TF) or "broad" (ヒストン修飾 H3K27me3 等)
        qvalue: FDR 閾値
    """
    import os
    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "macs3", "callpeak",
        "-t", treatment_bam,
        "-c", control_bam,
        "-g", str(genome_size),
        "--outdir", outdir,
        "-n", name,
        "-q", str(qvalue),
        "--keep-dup", "auto",
        "--call-summits",
    ]

    if peak_type == "broad":
        cmd.extend(["--broad", "--broad-cutoff", str(qvalue)])

    print(f"Running MACS3 peak calling ({peak_type} mode)...")
    subprocess.run(cmd, check=True)

    # ピークファイル読み込み
    suffix = "broadPeak" if peak_type == "broad" else "narrowPeak"
    peak_file = f"{outdir}/{name}_peaks.{suffix}"
    cols = ["chr", "start", "end", "name", "score", "strand",
            "signalValue", "pValue", "qValue"]
    if peak_type == "narrow":
        cols.append("summit")

    peaks = pd.read_csv(peak_file, sep="\t", header=None, names=cols)
    peaks["width"] = peaks["end"] - peaks["start"]

    print(f"  Called {len(peaks):,} {peak_type} peaks (q < {qvalue})")
    print(f"  Median peak width: {peaks['width'].median():.0f} bp")
    print(f"  Mean signal value: {peaks['signalValue'].mean():.2f}")

    return peaks


def chipseq_qc_metrics(peaks, frip_bam=None, total_reads=None):
    """
    ChIP-seq QC 指標の算出。

    Returns:
        dict: peak 数、中央値幅、FRiP (Fraction of Reads in Peaks)
    """
    metrics = {
        "n_peaks": len(peaks),
        "median_width_bp": float(peaks["width"].median()),
        "mean_signal": float(peaks["signalValue"].mean()),
        "mean_log10_qvalue": float(peaks["qValue"].mean()),
    }

    # ENCODE 品質基準との比較
    if metrics["n_peaks"] < 500:
        metrics["quality_flag"] = "LOW — < 500 peaks"
    elif metrics["n_peaks"] < 10000:
        metrics["quality_flag"] = "MODERATE"
    else:
        metrics["quality_flag"] = "HIGH"

    return metrics
```

## 2. ATAC-seq アクセシビリティ解析

```python
import numpy as np
import pandas as pd


def atacseq_nucleosome_free_regions(fragments_file, output_dir="results/atacseq"):
    """
    ATAC-seq フラグメントサイズ分布に基づくヌクレオソーム占有解析。

    フラグメントサイズによる分類:
    - < 150 bp: Nucleosome-Free Region (NFR)
    - 150-300 bp: Mono-nucleosome
    - 300-500 bp: Di-nucleosome
    - > 500 bp: Tri-nucleosome+
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    # フラグメントサイズ分布
    fragments = pd.read_csv(fragments_file, sep="\t",
                            names=["chr", "start", "end", "barcode", "count"])
    fragments["length"] = fragments["end"] - fragments["start"]

    # サイズカテゴリ分類
    bins = [0, 150, 300, 500, 10000]
    labels = ["NFR (<150)", "Mono-nuc (150-300)",
              "Di-nuc (300-500)", "Tri-nuc+ (>500)"]
    fragments["category"] = pd.cut(fragments["length"], bins=bins, labels=labels)

    size_dist = fragments["category"].value_counts(normalize=True)
    nfr_ratio = size_dist.get("NFR (<150)", 0)

    print(f"  Fragment size distribution:")
    for cat, pct in size_dist.items():
        print(f"    {cat}: {pct:.1%}")
    print(f"  NFR ratio: {nfr_ratio:.1%} (ENCODE target: >40%)")

    return fragments, size_dist


def atacseq_tss_enrichment(peaks, gene_gtf, window=2000):
    """
    TSS (Transcription Start Site) 周辺のシグナル濃縮スコア算出。
    TSS Enrichment Score > 7: 高品質 ATAC-seq データ (ENCODE 基準)。
    """
    from pybedtools import BedTool

    peaks_bt = BedTool.from_dataframe(
        peaks[["chr", "start", "end", "name", "score"]]
    )
    # TSS ± window bp のウィンドウ
    # (GTF parsing は省略 — 実運用時は pyranges/GTFparse 使用)

    print(f"  TSS enrichment window: ±{window} bp")
    print(f"  Total peaks overlapping TSS regions: calculated post-intersection")

    return peaks_bt
```

## 3. DNA メチル化パターン解析

```python
import numpy as np
import pandas as pd


def bisulfite_methylation_analysis(methylation_file, min_coverage=10,
                                    output_prefix="results/methylation"):
    """
    WGBS/RRBS バイサルファイトシーケンシングデータのメチル化解析。

    入力: Bismark methylation extractor 出力 (CpG context)
    処理:
    1. カバレッジフィルタリング
    2. メチル化レベル算出 (β 値)
    3. CpG アイランド/ショア/シェルフ注釈
    4. 差次メチル化領域 (DMR) 検出
    """
    import os
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    # Bismark 出力の読み込み
    df = pd.read_csv(methylation_file, sep="\t",
                     names=["chr", "pos", "strand", "count_m", "count_u"])
    df["coverage"] = df["count_m"] + df["count_u"]
    df["beta"] = df["count_m"] / df["coverage"]

    # カバレッジフィルタ
    n_before = len(df)
    df = df[df["coverage"] >= min_coverage].copy()
    print(f"  Coverage filter (≥{min_coverage}x): {n_before:,} → {len(df):,} CpGs")

    # グローバルメチル化統計
    mean_beta = df["beta"].mean()
    median_beta = df["beta"].median()
    print(f"  Global methylation: mean β = {mean_beta:.3f}, median β = {median_beta:.3f}")

    # メチル化状態分類
    df["status"] = pd.cut(df["beta"],
                          bins=[0, 0.2, 0.8, 1.0],
                          labels=["hypo", "intermediate", "hyper"])
    status_counts = df["status"].value_counts(normalize=True)
    print(f"  Hypomethylated (β<0.2): {status_counts.get('hypo', 0):.1%}")
    print(f"  Intermediate (0.2≤β≤0.8): {status_counts.get('intermediate', 0):.1%}")
    print(f"  Hypermethylated (β>0.8): {status_counts.get('hyper', 0):.1%}")

    return df


def detect_dmrs(group1_betas, group2_betas, positions, min_cpgs=5,
                delta_beta_cutoff=0.2, pvalue_cutoff=0.05):
    """
    差次メチル化領域 (DMR) 検出。

    Parameters:
        group1_betas, group2_betas: n_cpgs × n_samples メチル化マトリクス
        min_cpgs: DMR 内の最小 CpG 数
        delta_beta_cutoff: Δβ 閾値
    """
    from scipy.stats import mannwhitneyu

    results = []
    mean_g1 = group1_betas.mean(axis=1)
    mean_g2 = group2_betas.mean(axis=1)
    delta_beta = mean_g2 - mean_g1

    for i in range(len(positions)):
        stat, pval = mannwhitneyu(
            group1_betas[i, :], group2_betas[i, :], alternative="two-sided"
        )
        results.append({
            "chr": positions[i]["chr"],
            "pos": positions[i]["pos"],
            "delta_beta": float(delta_beta[i]),
            "pvalue": pval,
            "mean_group1": float(mean_g1[i]),
            "mean_group2": float(mean_g2[i]),
        })

    df = pd.DataFrame(results)

    # 多重検定補正
    from statsmodels.stats.multitest import multipletests
    df["padj"] = multipletests(df["pvalue"], method="fdr_bh")[1]

    # DMR フィルタ
    sig = df[(df["padj"] < pvalue_cutoff) &
             (df["delta_beta"].abs() >= delta_beta_cutoff)]
    print(f"  Significant DMCs (Δβ≥{delta_beta_cutoff}, FDR<{pvalue_cutoff}): {len(sig):,}")

    return df, sig
```

## 4. クロマチン状態モデリング (ChromHMM)

```python
import subprocess
import pandas as pd
import numpy as np


def chromhmm_learn_model(binarized_dir, output_dir, n_states=15,
                         assembly="hg38"):
    """
    ChromHMM によるクロマチン状態モデリング。

    複数のヒストン修飾マーク (H3K4me1/me3, H3K27ac, H3K27me3,
    H3K36me3, H3K9me3 等) を入力として、ゲノムをクロマチン状態に分類。

    Roadmap Epigenomics 15-state モデル:
    1-TssA, 2-TssAFlnk, 3-TxFlnk, 4-Tx, 5-TxWk,
    6-EnhG, 7-Enh, 8-ZNF/Rpts, 9-Het, 10-TssBiv,
    11-BivFlnk, 12-EnhBiv, 13-ReprPC, 14-ReprPCWk, 15-Quies
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        "java", "-mx8G", "-jar", "ChromHMM.jar", "LearnModel",
        "-b", "200",
        binarized_dir, output_dir, str(n_states), assembly
    ]
    print(f"Running ChromHMM LearnModel with {n_states} states...")
    subprocess.run(cmd, check=True)

    # 遷移確率マトリクス読み込み
    trans_file = f"{output_dir}/transitions_{n_states}.txt"
    if os.path.exists(trans_file):
        trans = pd.read_csv(trans_file, sep="\t", index_col=0)
        print(f"  Transition matrix: {trans.shape}")

    # エミッション確率読み込み
    emit_file = f"{output_dir}/emissions_{n_states}.txt"
    if os.path.exists(emit_file):
        emit = pd.read_csv(emit_file, sep="\t", index_col=0)
        print(f"  Emission matrix: {emit.shape}")

    return {"n_states": n_states, "output_dir": output_dir}


def annotate_chromatin_states(segments_bed, state_labels=None):
    """
    ChromHMM セグメンテーション結果のゲノム注釈。

    Parameters:
        segments_bed: ChromHMM 出力の *_segments.bed
        state_labels: 状態番号→機能ラベルのマッピング辞書
    """
    default_labels = {
        "E1": "Active TSS", "E2": "Flanking Active TSS",
        "E3": "Transcription at gene 5'/3'", "E4": "Strong Transcription",
        "E5": "Weak Transcription", "E6": "Genic Enhancers",
        "E7": "Enhancers", "E8": "ZNF genes & Repeats",
        "E9": "Heterochromatin", "E10": "Bivalent/Poised TSS",
        "E11": "Flanking Bivalent TSS/Enh", "E12": "Bivalent Enhancer",
        "E13": "Repressed PolyComb", "E14": "Weak Repressed PolyComb",
        "E15": "Quiescent/Low",
    }
    labels = state_labels or default_labels

    segments = pd.read_csv(segments_bed, sep="\t",
                           names=["chr", "start", "end", "state"])
    segments["width"] = segments["end"] - segments["start"]
    segments["label"] = segments["state"].map(labels)

    # ゲノムカバレッジ統計
    total_bp = segments["width"].sum()
    state_coverage = segments.groupby("label")["width"].sum() / total_bp
    print("  Chromatin state genome coverage:")
    for label, pct in state_coverage.sort_values(ascending=False).items():
        print(f"    {label}: {pct:.1%}")

    return segments, state_coverage
```

## 5. Hi-C 3D ゲノム構造解析

```python
import numpy as np
import pandas as pd


def hic_contact_matrix_analysis(cool_file, resolution=10000,
                                 chromosome="chr1"):
    """
    Hi-C 接触マップ解析 (.cool/.mcool 形式)。

    1. ICE 正規化
    2. A/B コンパートメント同定 (PCA)
    3. TAD 呼び出し (Insulation Score)
    """
    import cooler

    # クールファイル読み込み
    clr = cooler.Cooler(f"{cool_file}::resolutions/{resolution}")
    matrix = clr.matrix(balance=True).fetch(chromosome)

    print(f"  Contact matrix shape: {matrix.shape}")
    print(f"  Resolution: {resolution:,} bp")
    print(f"  Non-zero entries: {np.count_nonzero(~np.isnan(matrix)):,}")

    return matrix


def call_tads_insulation_score(matrix, resolution=10000, window_size=500000):
    """
    Insulation Score 法による TAD (Topologically Associating Domain) 呼び出し。

    Parameters:
        window_size: Insulation window サイズ (bp)
    """
    window_bins = window_size // resolution

    n = matrix.shape[0]
    insulation = np.zeros(n)

    for i in range(window_bins, n - window_bins):
        submat = matrix[i - window_bins:i, i:i + window_bins]
        insulation[i] = np.nanmean(submat)

    # log2 正規化
    mean_val = np.nanmean(insulation[insulation > 0])
    log_insulation = np.log2(insulation / mean_val + 1e-10)

    # TAD 境界 = Insulation Score の極小値
    from scipy.signal import argrelextrema
    minima = argrelextrema(log_insulation, np.less, order=5)[0]

    tad_boundaries = minima * resolution
    n_tads = len(tad_boundaries) - 1

    print(f"  Found {len(tad_boundaries)} TAD boundaries")
    print(f"  Estimated {n_tads} TADs")
    print(f"  Mean TAD size: {np.diff(tad_boundaries).mean() / 1e6:.2f} Mb")

    return log_insulation, tad_boundaries


def ab_compartment_analysis(matrix, resolution=100000):
    """
    Hi-C データからの A/B コンパートメント同定。

    A コンパートメント: euchromatin, 活性, 遺伝子リッチ
    B コンパートメント: heterochromatin, 不活性, 遺伝子プア
    """
    from sklearn.decomposition import PCA

    # O/E (Observed/Expected) マトリクス
    matrix_clean = np.nan_to_num(matrix, nan=0.0)
    expected = np.zeros_like(matrix_clean)
    for d in range(matrix_clean.shape[0]):
        diag_vals = np.diag(matrix_clean, d)
        mean_val = np.mean(diag_vals) if len(diag_vals) > 0 else 0
        np.fill_diagonal(expected[d:, :], mean_val)
        np.fill_diagonal(expected[:, d:], mean_val)

    oe_matrix = matrix_clean / (expected + 1e-10)

    # 相関マトリクス → PCA
    corr_matrix = np.corrcoef(oe_matrix)
    corr_matrix = np.nan_to_num(corr_matrix)

    pca = PCA(n_components=2)
    components = pca.fit_transform(corr_matrix)
    pc1 = components[:, 0]

    # A/B 分類 (PC1 正 = A, 負 = B)
    compartment = np.where(pc1 > 0, "A", "B")
    a_frac = np.mean(compartment == "A")

    print(f"  A compartment: {a_frac:.1%}")
    print(f"  B compartment: {1 - a_frac:.1%}")
    print(f"  PC1 variance explained: {pca.explained_variance_ratio_[0]:.1%}")

    return pc1, compartment
```

## 6. 転写因子モチーフ濃縮解析

```python
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact


def motif_enrichment_analysis(peak_sequences, background_sequences,
                               jaspar_db="JASPAR2024_CORE_vertebrates",
                               pvalue_cutoff=0.01):
    """
    ピーク領域における転写因子結合モチーフの濃縮解析。

    Parameters:
        peak_sequences: FASTA ファイル (ピーク中心 ±250 bp)
        background_sequences: ランダムゲノム領域 FASTA
        jaspar_db: JASPAR データベースバージョン
    """
    from Bio import motifs

    # JASPAR PWM スキャン (概念的コード)
    results = []

    # Homer / MEME-ChIP 呼び出し (実運用)
    print(f"  Scanning {jaspar_db} motifs against peak sequences...")
    print("  (Using FIMO from MEME Suite for motif scanning)")

    # Fisher exact test による濃縮
    # peak_hits / peak_total vs bg_hits / bg_total
    # for each motif in JASPAR database

    return results


def differential_binding_analysis(sample_sheet, peaks_dir,
                                   contrast=("Treatment", "Control"),
                                   fdr_cutoff=0.05, fold_change_cutoff=2):
    """
    DiffBind による差次結合解析。

    Parameters:
        sample_sheet: DiffBind サンプルシート CSV
        contrast: (treatment, control) 比較群
        fdr_cutoff: FDR 閾値
        fold_change_cutoff: log2FC 閾値
    """
    # R/rpy2 経由で DiffBind を呼び出し
    import subprocess
    r_script = f"""
    library(DiffBind)
    samples <- read.csv("{sample_sheet}")
    dba <- dba(sampleSheet=samples)
    dba <- dba.count(dba)
    dba <- dba.contrast(dba, categories=DBA_CONDITION)
    dba <- dba.analyze(dba)
    db_sites <- dba.report(dba, th={fdr_cutoff}, fold={np.log2(fold_change_cutoff)})
    write.csv(as.data.frame(db_sites), "results/diffbind_results.csv")
    """

    print(f"  Running DiffBind: {contrast[0]} vs {contrast[1]}")
    print(f"  FDR cutoff: {fdr_cutoff}, log2FC cutoff: ±{np.log2(fold_change_cutoff):.1f}")

    return r_script
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/chipseq/{name}_peaks.narrowPeak` | BED/narrowPeak |
| `results/chipseq/{name}_peaks.broadPeak` | BED/broadPeak |
| `results/atacseq/fragment_size_dist.csv` | CSV |
| `results/methylation/dmr_results.csv` | CSV |
| `results/chromhmm/emissions_{n}.txt` | TSV |
| `results/hic/tad_boundaries.bed` | BED |
| `results/hic/compartments.csv` | CSV |
| `results/diffbind_results.csv` | CSV |
| `figures/chromatin_state_heatmap.png` | PNG |
| `figures/hic_contact_map.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| ChIP-Atlas | `ChIPAtlas_enrichment_analysis` | TF/ヒストン修飾エンリッチメント解析 |
| ChIP-Atlas | `ChIPAtlas_get_experiments` | 実験メタデータ取得 (43 万+実験) |
| ChIP-Atlas | `ChIPAtlas_get_peak_data` | ピークコールデータ取得 |
| ChIP-Atlas | `ChIPAtlas_search_datasets` | データセット検索 (抗原/細胞種) |
| 4DN | `FourDN_search_data` | Hi-C/ChIA-PET 3D ゲノムデータ検索 |
| JASPAR | `jaspar_search_matrices` | 転写因子結合モチーフ (PWM) 検索 |
| JASPAR | `jaspar_get_matrix` | PWM (Position Weight Matrix) 取得 |
| JASPAR | `jaspar_list_collections` | JASPAR コレクション一覧 |
| SCREEN | `SCREEN_get_regulatory_elements` | cCRE (候補シス調節エレメント) 取得 |
| ENCODE | `ENCODE_search_experiments` | ENCODE ChIP-seq/ATAC-seq 実験検索 |
| ENCODE | `ENCODE_get_experiment` | ENCODE 実験詳細取得 |
| ENCODE | `ENCODE_list_files` | ENCODE ファイル一覧 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-single-cell-genomics` | scATAC-seq 連携 |
| `scientific-sequence-analysis` | ゲノム配列操作 |
| `scientific-bioinformatics` | BAM/VCF 処理 |
| `scientific-population-genetics` | eQTL・調節バリアント |
| `scientific-gene-expression-transcriptomics` | 発現-エピゲノム統合 |

### 依存パッケージ

`macs3`, `cooler`, `pybedtools`, `deeptools`, `scikit-learn`, `scipy`, `pandas`, `numpy`, `biopython`
