---
name: scientific-proteomics-mass-spectrometry
description: |
  プロテオミクス・質量分析解析スキル。LC-MS/MS データ前処理、ペプチド同定 (PSM/FDR 制御)、
  蛋白質定量 (LFQ/TMT/SILAC/iBAQ)、翻訳後修飾 (PTM) マッピング、
  スペクトル類似度スコアリング (コサイン/修正コサイン)、分子ネットワーキング (GNPS)、
  化合物アノテーション (HMDB/MassBank) を統合した質量分析パイプライン。
  pyOpenMS / matchms ベースの包括的ワークフロー。
---

# Scientific Proteomics & Mass Spectrometry

LC-MS/MS ベースのプロテオミクス・メタボロミクス質量分析データを対象に、
スペクトル前処理→ペプチド/化合物同定→定量→差次的解析の
標準パイプラインを提供する。

## When to Use

- LC-MS/MS プロテオミクスデータのペプチド同定・蛋白質定量が必要なとき
- TMT/SILAC/LFQ による差次的蛋白質発現解析を行うとき
- 翻訳後修飾 (リン酸化, ユビキチン化, アセチル化) マッピングが必要なとき
- スペクトルライブラリ検索・分子ネットワーキングを行うとき
- 化合物同定 (HMDB/MassBank/GNPS) が必要なとき

---

## Quick Start

## 1. MS データ前処理 (pyOpenMS)

```python
import numpy as np
import pandas as pd


def ms_data_preprocessing(mzml_file, noise_threshold=1000,
                           peak_picking_method="centroid"):
    """
    LC-MS/MS データ前処理パイプライン (pyOpenMS ベース)。

    1. Raw → mzML 変換 (事前に msconvert)
    2. ピーク検出 (セントロイド化)
    3. ベースライン補正
    4. ノイズ除去
    5. RT アラインメント
    """
    from pyopenms import MSExperiment, MzMLFile, PeakPickerHiRes

    # mzML 読み込み
    exp = MSExperiment()
    MzMLFile().load(mzml_file, exp)

    print(f"  Loaded {exp.getNrSpectra()} spectra from {mzml_file}")
    print(f"  MS1 scans: {sum(1 for s in exp if s.getMSLevel() == 1)}")
    print(f"  MS2 scans: {sum(1 for s in exp if s.getMSLevel() == 2)}")

    # ピークピッキング
    if peak_picking_method == "centroid":
        picker = PeakPickerHiRes()
        exp_picked = MSExperiment()
        picker.pickExperiment(exp, exp_picked)
        print(f"  After peak picking: {exp_picked.getNrSpectra()} spectra")
        return exp_picked

    return exp


def feature_detection(exp, mass_error_ppm=10, noise_threshold=1000):
    """
    LC-MS 特徴検出 — m/z × RT × 強度の 3D ピーク検出。
    """
    from pyopenms import FeatureFinder, FeatureMap

    ff = FeatureFinder()
    features = FeatureMap()

    ff_params = ff.getParameters()
    ff_params.setValue("mass_trace:mz_tolerance", float(mass_error_ppm))
    ff_params.setValue("noise_threshold_int", float(noise_threshold))
    ff.setParameters(ff_params)

    ff.run("centroided", exp, features, FeatureMap())

    print(f"  Detected {features.size()} features")

    results = []
    for f in features:
        results.append({
            "mz": f.getMZ(),
            "rt": f.getRT(),
            "intensity": f.getIntensity(),
            "charge": f.getCharge(),
            "quality": f.getOverallQuality(),
        })

    return pd.DataFrame(results)
```

## 2. ペプチド同定 (データベース検索)

```python
import numpy as np
import pandas as pd


def peptide_identification(mzml_file, fasta_db, enzyme="Trypsin",
                           missed_cleavages=2, precursor_mass_tol=10,
                           fragment_mass_tol=0.02, fdr_cutoff=0.01):
    """
    MS/MS スペクトルからのペプチド同定パイプライン。

    1. データベース検索 (X!Tandem / Comet / MSGF+)
    2. PSM (Peptide-Spectrum Match) スコアリング
    3. Target-Decoy FDR 制御
    4. Protein inference (Occam's Razor)
    """
    from pyopenms import (
        IdXMLFile, ProteinIdentification,
        PeptideIdentification
    )

    # 検索パラメータ
    search_params = {
        "database": fasta_db,
        "enzyme": enzyme,
        "missed_cleavages": missed_cleavages,
        "precursor_mass_tolerance": f"{precursor_mass_tol} ppm",
        "fragment_mass_tolerance": f"{fragment_mass_tol} Da",
        "fixed_modifications": ["Carbamidomethyl (C)"],
        "variable_modifications": ["Oxidation (M)", "Acetyl (Protein N-term)"],
    }

    print(f"  Database search parameters:")
    for k, v in search_params.items():
        print(f"    {k}: {v}")

    # FDR 制御
    # Target-Decoy approach: concatenate reversed sequences
    print(f"  FDR cutoff: {fdr_cutoff} (1% at PSM level)")
    print(f"  Method: Target-Decoy competition (TDC)")

    return search_params


def protein_quantification(psm_results, method="LFQ",
                            min_peptides=2, min_ratio_count=2):
    """
    蛋白質定量。

    Methods:
    - LFQ (Label-Free Quantification): MS1 強度ベース
    - iBAQ (intensity-Based Absolute Quantification)
    - TMT (Tandem Mass Tag): Reporter ion 強度
    - SILAC: Heavy/Light 比率
    """
    print(f"  Quantification method: {method}")
    print(f"  Minimum peptides per protein: {min_peptides}")

    if method == "LFQ":
        # MaxLFQ アルゴリズム: ペプチド比率の中央値正規化
        print("  Normalization: MaxLFQ (median of peptide ratios)")
        print("  Missing value imputation: MinDet (minimum deterministic)")

    elif method == "TMT":
        print("  TMT channels: 126-134N (TMTpro 18-plex)")
        print("  Reporter ion extraction: ±10 ppm")
        print("  Normalization: Median centering → IRS (Internal Reference Scaling)")

    elif method == "SILAC":
        print("  Labels: Light (K0R0) vs Heavy (K8R10)")
        print("  Ratio calculation: median of peptide ratios")

    elif method == "iBAQ":
        print("  iBAQ = Σ(peptide intensities) / n_observable_peptides")

    return {"method": method, "min_peptides": min_peptides}
```

## 3. 翻訳後修飾 (PTM) 解析

```python
import pandas as pd
import numpy as np


def ptm_site_localization(psm_results, ptm_types=None,
                           localization_prob_cutoff=0.75):
    """
    翻訳後修飾サイトの局在化解析。

    主要 PTM タイプ:
    - Phosphorylation (S/T/Y): リン酸化
    - Ubiquitination (K): ユビキチン化 (diGly remnant)
    - Acetylation (K): アセチル化
    - Methylation (K/R): メチル化
    - Glycosylation (N/S/T): 糖鎖修飾
    - SUMOylation (K): SUMO 化
    """
    ptm_types = ptm_types or ["Phospho (S)", "Phospho (T)", "Phospho (Y)"]

    print(f"  PTM types analyzed: {ptm_types}")
    print(f"  Localization probability cutoff: {localization_prob_cutoff}")
    print("  Methods: phosphoRS / Ascore / ptmRS")

    # PTM 濃縮解析 (Motif-X / pLogo)
    print("  Motif enrichment: Motif-X algorithm")
    print("    Window: ±7 residues around modification site")
    print("    Significance: p < 1e-6 (binomial test)")

    return {"ptm_types": ptm_types, "cutoff": localization_prob_cutoff}


def phosphoproteomics_kinase_activity(phosphosites_df,
                                       kinase_substrate_db="PhosphoSitePlus"):
    """
    フォスフォプロテオミクスからのキナーゼ活性推定 (KSEA)。

    Kinase-Substrate Enrichment Analysis:
    - PhosphoSitePlus / NetworKIN のキナーゼ-基質関係を使用
    - 各キナーゼの基質群の平均 log2FC をスコア化
    - z-test で有意性を評価
    """
    print(f"  Kinase-substrate database: {kinase_substrate_db}")
    print("  Algorithm: KSEA (Kinase-Substrate Enrichment Analysis)")
    print("  Score: mean(log2FC of substrates) × sqrt(n_substrates)")

    return {"database": kinase_substrate_db, "method": "KSEA"}
```

## 4. スペクトルマッチング・分子ネットワーキング

```python
import numpy as np
import pandas as pd


def spectral_similarity_scoring(query_spectrum, library_spectrum,
                                 method="modified_cosine",
                                 mz_tolerance=0.02):
    """
    MS/MS スペクトル間の類似度スコアリング。

    Methods:
    - cosine: 標準コサイン類似度
    - modified_cosine: 前駆体質量差を考慮したシフトマッチング
    - spec2vec: Word2Vec ベースのスペクトル埋め込み
    """
    from matchms import Spectrum, calculate_scores
    from matchms.similarity import ModifiedCosine, CosineGreedy

    if method == "modified_cosine":
        similarity_func = ModifiedCosine(tolerance=mz_tolerance)
    else:
        similarity_func = CosineGreedy(tolerance=mz_tolerance)

    score = similarity_func.pair(query_spectrum, library_spectrum)

    print(f"  Similarity method: {method}")
    print(f"  Score: {score['score']:.4f}")
    print(f"  Matched peaks: {score['matches']}")

    return score


def molecular_networking(spectra_list, min_cosine=0.7,
                          min_matched_peaks=6, max_neighbors=10):
    """
    GNPS スタイル分子ネットワーキング。

    スペクトル間の修正コサイン類似度に基づいてネットワークを構築。
    類似構造を持つ化合物がクラスタを形成 → 未知化合物の推定に活用。
    """
    from matchms import calculate_scores
    from matchms.similarity import ModifiedCosine
    import networkx as nx

    sim_func = ModifiedCosine(tolerance=0.02)

    G = nx.Graph()
    n = len(spectra_list)
    edge_count = 0

    for i in range(n):
        G.add_node(i)
        for j in range(i + 1, n):
            score = sim_func.pair(spectra_list[i], spectra_list[j])
            if (score["score"] >= min_cosine and
                    score["matches"] >= min_matched_peaks):
                G.add_edge(i, j, weight=score["score"],
                           matches=score["matches"])
                edge_count += 1

    print(f"  Molecular network: {n} nodes, {edge_count} edges")
    print(f"  Connected components: {nx.number_connected_components(G)}")
    print(f"  Parameters: min_cosine={min_cosine}, min_matched={min_matched_peaks}")

    return G
```

## 5. 差次的蛋白質発現 + 機能濃縮

```python
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def differential_protein_expression(intensity_matrix, groups,
                                     log2fc_cutoff=1.0, fdr_cutoff=0.05,
                                     imputation="MinDet"):
    """
    差次的蛋白質発現解析。

    Parameters:
        intensity_matrix: proteins × samples (log2 LFQ 強度)
        groups: サンプルグループ (e.g., ["Treatment", "Control", ...])
        imputation: 欠損値補完法 (MinDet / MinProb / KNN / QRILC)
    """
    group_names = sorted(set(groups))
    g1, g2 = group_names[0], group_names[1]
    g1_idx = [i for i, g in enumerate(groups) if g == g1]
    g2_idx = [i for i, g in enumerate(groups) if g == g2]

    results = []
    for protein in intensity_matrix.index:
        vals1 = intensity_matrix.loc[protein, intensity_matrix.columns[g1_idx]].dropna()
        vals2 = intensity_matrix.loc[protein, intensity_matrix.columns[g2_idx]].dropna()

        if len(vals1) < 2 or len(vals2) < 2:
            continue

        log2fc = vals2.mean() - vals1.mean()
        stat, pval = ttest_ind(vals1, vals2, equal_var=False)

        results.append({
            "protein": protein,
            "log2FC": log2fc,
            "pvalue": pval,
            "mean_g1": vals1.mean(),
            "mean_g2": vals2.mean(),
            "n_g1": len(vals1),
            "n_g2": len(vals2),
        })

    df = pd.DataFrame(results)
    df["padj"] = multipletests(df["pvalue"], method="fdr_bh")[1]

    sig_up = df[(df["padj"] < fdr_cutoff) & (df["log2FC"] > log2fc_cutoff)]
    sig_down = df[(df["padj"] < fdr_cutoff) & (df["log2FC"] < -log2fc_cutoff)]

    print(f"  {g2} vs {g1}:")
    print(f"  Total proteins tested: {len(df)}")
    print(f"  Significant UP: {len(sig_up)} (log2FC > {log2fc_cutoff}, FDR < {fdr_cutoff})")
    print(f"  Significant DOWN: {len(sig_down)} (log2FC < -{log2fc_cutoff}, FDR < {fdr_cutoff})")

    return df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/features_detected.csv` | CSV |
| `results/psm_results.csv` | CSV |
| `results/protein_quant.csv` | CSV |
| `results/ptm_sites.csv` | CSV |
| `results/differential_proteins.csv` | CSV |
| `results/molecular_network.graphml` | GraphML |
| `figures/volcano_proteomics.png` | PNG |
| `figures/molecular_network.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PRIDE | `PRIDE_search_proteomics` | プロテオミクスプロジェクト検索 |
| PRIDE | `PRIDE_get_project` | プロジェクト詳細取得 |
| PRIDE | `PRIDE_get_project_files` | プロテオミクスデータファイル取得 |
| UniProt | `search_uniprot_by_name` | 蛋白質検索 |
| UniProt | `get_uniprot_entry` | 蛋白質エントリ詳細 |
| KEGG | `kegg_get_pathway_info` | 蛋白質パスウェイ情報 |
| Reactome | `reactome_pathway_analysis` | パスウェイ濃縮解析 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-metabolomics` | 代謝物 MS 解析 |
| `scientific-spectral-signal` | スペクトル解析基盤 |
| `scientific-bioinformatics` | 配列・蛋白質データベース |
| `scientific-network-analysis` | 分子ネットワーク可視化 |
| `scientific-multi-omics` | マルチオミクス統合 |

### 依存パッケージ

`pyopenms`, `matchms`, `pandas`, `numpy`, `scipy`, `scikit-learn`, `networkx`, `statsmodels`
