# SATORI（悟り）— Agent Skills for Science

**SATORI** は、科学データ解析のための **GitHub Copilot Agent Skills** コレクションです。

[![npm version](https://img.shields.io/npm/v/@nahisaho/satori)](https://www.npmjs.com/package/@nahisaho/satori)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENCE)

## Overview

このディレクトリには、Exp-01〜13 で蓄積した科学データ解析技法を Agent Skills として体系化した **76 個**のスキルを格納しています。Copilot がプロンプトの文脈に応じて適切なスキルを自動ロードし、各実験で確立した解析パターンを再利用します。42 のスキルは [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で 1,200 以上の外部科学データベースツールとも連携可能です。

### パイプラインフロー

```
hypothesis-pipeline → pipeline-scaffold → academic-writing → critical-review
  (仮説定義)         (解析実行)         (草稿作成)         (レビュー・修正)
                                                                       ↓
  paper-quality ← revision-tracker ← peer-review-response ← [査読結果受領]
  (品質評価)    (改訂追跡)       (査読対応)         (ジャーナル)
```

**ドメイン特化パイプライン（J-O）**

```
research-methodology → grant-writing → hypothesis-pipeline    ← [O→A 研究計画]
                                              ↓
drug-target-profiling → admet-pharmacokinetics ─→ drug-repurposing
  (標的同定)           (ADMET/PK 評価)        (リポジショニング)
        ↓                                            ↓
protein-structure-analysis → protein-design → lab-automation
  (構造解析)               (de novo 設計)   (実験自動化)
                                                     ↓
variant-interpretation → clinical-decision-support → presentation-design
  (バリアント解釈)       (臨床意思決定)            (学会発表)
```

各ステップで生成されるファイルが次のステップに自動的に引き継がれます：

**先端計算・医用パイプライン（P-S）**

```
pharmacovigilance ← admet-pharmacokinetics       ← [P 安全性監視]
  (市販後安全性)     (前臨床 ADMET)                   ↓
precision-oncology → clinical-decision-support → medical-imaging
  (腫瘍プロファイル)  (臨床意思決定)            (画像診断)
        ↓                                           ↓
disease-research → variant-interpretation      → deep-learning
  (疾患-遺伝子)    (バリアント解釈)            (DL フレームワーク)
                                                     ↓
quantum-computing → bayesian-statistics → graph-neural-networks
  (量子計算)       (ベイズ推論)          (GNN 分子予測)
                                               ↓
                  explainable-ai ← deep-learning ← [R 先端計算]
                  (XAI 説明可能性)  (DL パイプライン)
```

**次世代オミクス・疫学パイプライン（T-Z）**

```
single-cell-genomics → spatial-transcriptomics     ← [T シングルセル・空間]
  (scRNA-seq QC)       (Visium/MERFISH)
        ↓                     ↓
immunoinformatics → infectious-disease             ← [U 免疫・感染症]
  (エピトープ予測)    (AMR・系統解析)
        ↓                     ↓
microbiome-metagenomics → environmental-ecology    ← [V マイクロバイオーム・環境]
  (16S/メタゲノム)         (SDM・生物多様性)
        ↓                     ↓
systems-biology           population-genetics       ← [W+Y モデル・集団]
  (SBML/FBA/GRN)           (Fst/ADMIXTURE)
        ↓                     ↓
epidemiology-public-health → text-mining-nlp       ← [X+Z 疫学・NLP]
  (RR/OR/空間クラスタ)       (NER/KG/BERTopic)
```

| フェーズ | 生成ファイル | 参照先 |
|---|---|---|
| 仮説立案 | `docs/hypothesis.{md,json}`, `docs/workflow_design.{md,json}` | → scaffold, writing |
| 解析実行 | `results/analysis_summary.json`, `figures/*.png` | → writing |
| 草稿作成 | `manuscript/manuscript.md` | → critical-review, citation-checker |
| レビュー | `manuscript/review_report.{md,json}`, `manuscript/manuscript_revised.md` | → latex-formatter |
| 引用検証 | `manuscript/citation_report.json` | → latex-formatter |
| SI 生成 | `manuscript/supplementary.md`, `manuscript/si_crossref_report.json` | → latex-formatter |
| 査読対応 | `manuscript/response_to_reviewers.md`, `manuscript/response_mapping.json` | → revision-tracker |
| 改訂追跡 | `manuscript/manuscript_tracked.md`, `manuscript/revision_summary.json` | → paper-quality |
| 品質評価 | `manuscript/quality_report.json` | → latex-formatter |
| LaTeX 変換 | `manuscript/manuscript.tex`, `manuscript/references.bib` | — |
| 標的プロファイリング | `results/target_profile_report.md`, `results/target_profile.json` | → admet-pk, protein-structure, drug-repurposing |
| ADMET/PK 評価 | `results/admet_profile.json`, `results/pk_model.json` | → drug-repurposing, clinical-decision |
| ドラッグリポジショニング | `results/repurposing_candidates.json`, `results/network_proximity.json` | → clinical-decision |
| 構造解析 | `results/structure_analysis.json`, `results/binding_sites.json` | → protein-design, cheminformatics |
| タンパク質設計 | `results/design_candidates.json`, `results/esm_scores.json` | → lab-automation, admet-pk |
| バリアント解釈 | `results/variant_classification.json`, `results/pgx_report.json` | → clinical-decision |
| 臨床意思決定 | `results/clinical_recommendation.json`, `results/trial_matches.json` | → presentation, writing |
| 実験室自動化 | `protocols/protocol.py`, `results/qc_report.json` | → data-preprocessing |
| プレゼンテーション | `presentation/slides.md`, `presentation/poster.tex` | — |
| 研究方法論 | `docs/methodology_design.md`, `docs/study_design.json` | → grant-writing, doe |
| グラント | `grants/specific_aims.md`, `grants/research_strategy.md` | → hypothesis-pipeline |
| ファーマコビジランス | `results/pv_signal_report.{md,json}`, `figures/pv_temporal_trend.png` | → clinical-decision |
| 精密腫瘍学 | `results/mtb_report.{md,json}`, `results/variant_actionability.json` | → clinical-decision, writing |
| 疾患研究 | `results/disease_research_report.{md,json}`, `results/gwas_significant_loci.json` | → variant-interpretation |
| 量子計算 | `results/quantum_result.json`, `figures/quantum_convergence.png` | → bayesian, cheminformatics |
| GNN | `results/gnn_predictions.json`, `figures/gnn_training_curve.png` | → drug-target, admet |
| ベイズ統計 | `results/bayesian_summary.json`, `figures/bayesian_trace.png` | → doe, meta-analysis |
| 説明可能 AI | `results/xai_report.json`, `figures/shap_summary.png` | → clinical-decision |
| 深層学習 | `results/dl_training_log.json`, `models/model.onnx` | → GNN, medical-imaging |
| 医用イメージング | `results/imaging_report.{md,json}`, `results/radiomics_features.json` | → precision-oncology |
| scRNA-seq 解析 | `results/sc_markers.json`, `figures/umap_clusters.png`, `results/rna_velocity.json` | → spatial, systems-biology |
| 空間トランスクリプトミクス | `results/spatial_domains.json`, `figures/spatial_svg_map.png` | → single-cell, systems-biology |
| 免疫情報学 | `results/epitope_candidates.json`, `results/tcr_diversity.json` | → infectious-disease, drug-target |
| 感染症ゲノミクス | `results/amr_report.json`, `results/mlst_profile.json`, `results/sir_simulation.json` | → epidemiology, microbiome |
| マイクロバイオーム | `results/asv_table.json`, `results/diversity_metrics.json`, `results/da_results.json` | → environmental-ecology |
| 環境生態学 | `results/sdm_predictions.json`, `results/biodiversity_indices.json` | → microbiome, text-mining |
| システム生物学 | `results/sbml_timecourse.json`, `results/fba_fluxes.json`, `results/grn_edges.json` | → multi-omics, network-analysis |
| 疫学・公衆衛生 | `results/epi_risk_measures.json`, `results/spatial_clusters.json`, `results/dag_analysis.json` | → survival-clinical, causal-inference |
| 集団遺伝学 | `results/pop_structure.json`, `results/fst_matrix.json`, `results/selection_scan.json` | → disease-research, variant-interpretation |
| テキストマイニング | `results/ner_entities.json`, `results/knowledge_graph.json`, `results/topic_model.json` | → deep-research, meta-analysis |
| 神経電気生理学 | `results/spike_sorting.json`, `results/eeg_erp.json`, `results/connectivity.json` | → biosignal, deep-learning |
| プロテオミクス | `results/protein_quant.csv`, `results/ptm_sites.json`, `results/molecular_network.json` | → multi-omics, network-analysis |
| トランスクリプトミクス | `results/deseq2_results.csv`, `results/gsea/`, `figures/volcano_rnaseq.png` | → bioinformatics, multi-omics |
| 計算材料科学 | `results/structure.cif`, `figures/phase_diagram.png`, `figures/band_structure.png` | → quantum-computing |
| 臨床試験解析 | `results/clinical_trials.csv`, `results/competitive_landscape.json` | → survival-clinical, meta-analysis |
| ラボデータ管理 | `results/benchling_sequences.json`, `results/dnanexus_workflow_output.json` | → bioinformatics, lab-automation |
| 科学図式 | `figures/consort_flow.svg`, `figures/nn_architecture.svg`, `figures/pathway.md` | → presentation, writing |
| 規制科学 | `results/fda_orange_book.json`, `results/510k_clearances.csv`, `results/patent_search.csv` | → pharmacovigilance, clinical-trials |
| 薬理ゲノミクス | `results/pgx_report.json`, `results/cpic_recommendations.json` | → variant-interpretation, clinical-decision |
| エピゲノミクス | `results/peak_calls.bed`, `results/dmr_results.csv`, `results/chromatin_states.bed` | → single-cell, multi-omics |

### ToolUniverse MCP ツール連携

42 のスキル（HIGH 13 + MEDIUM 9 + 新規 20）は、[ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP サーバー経由で 1,200 以上の外部科学ツールを利用可能です。各 SKILL.md 内の `### 利用可能ツール` セクションに対応ツールが記載されています。

```
SATORI Skill (方法論・判断)        ToolUniverse SMCP (データ取得・計算)
┌──────────────────────┐        ┌─────────────────────────────┐
│ pharmacovigilance    │───MCP──│ FAERS, FDA Labels, DailyMed │
│ precision-oncology   │───MCP──│ OncoKB, CIViC, COSMIC, GDC  │
│ disease-research     │───MCP──│ OpenTargets, HPO, Monarch    │
│ drug-target-profiling│───MCP──│ UniProt, ChEMBL, DGIdb       │
│ variant-interpretation│──MCP──│ ClinVar, gnomAD, ClinGen     │
│ admet-pharmacokinetics│──MCP──│ ADMET-AI, PubChem, ChEMBL    │
│ ... (42 skills total) │       │ ... (1,200+ tools)           │
└──────────────────────┘        └─────────────────────────────┘
```

スキルは **26 の中区分**に分類されています。

| 中区分 | スキル数 | 概要 |
|---|:---:|---|
| A. 基盤・ワークフロー | 13 | パイプライン構築・前処理・データ生成・図表・執筆・仮説立案・批判的レビュー・SI 生成・LaTeX 変換・引用検証・査読対応・改訂追跡・論文品質 |
| B. 統計・探索的解析 | 3 | EDA・仮説検定・次元削減 |
| C. 機械学習・モデリング | 3 | 回帰・分類・特徴量重要度 |
| D. 実験計画・プロセス最適化 | 2 | DOE・応答曲面法・ベイズ最適化 |
| E. 信号・スペクトル・時系列 | 4 | スペクトル解析・生体信号・時系列分解・神経電気生理学 |
| F. 生命科学・オミクス | 7 | バイオインフォ・メタボロ・ゲノム配列・マルチオミクス・ネットワーク・プロテオミクス・トランスクリプトミクス |
| G. 化学・材料・イメージング | 4 | ケモインフォ・材料特性評価・画像形態解析・計算材料科学 |
| H. 臨床・疫学・メタ科学 | 4 | 生存解析・因果推論・メタアナリシス・臨床試験解析 |
| I. Deep Research | 1 | 科学文献深層リサーチ・エビデンス階層評価・ソース追跡・交差検証 |
| J. 創薬・ファーマコロジー | 3 | 標的プロファイリング・ADMET/PK・ドラッグリポジショニング |
| K. 構造生物学・タンパク質工学 | 2 | PDB/AlphaFold 構造解析・de novo タンパク質設計 |
| L. 精密医療・臨床意思決定 | 2 | 変異解釈 (ACMG/AMP)・エビデンスベース臨床意思決定 |
| M. 実験室自動化・データ管理 | 2 | 液体ハンドリング・プロトコル管理・ELN/LIMS 連携・ラボデータ管理 |
| N. 科学プレゼンテーション・図式 | 2 | 科学スライド・ポスター・ワークフロー図・科学図式 |
| O. 研究計画・グラント・規制 | 3 | 助成金申請書・研究方法論・倫理審査・規制科学 |
| P. ファーマコビジランス・薬理ゲノミクス | 2 | FAERS 不均衡分析・MedDRA 階層・安全性シグナル検出・PGx 代謝型 |
| Q. 腫瘍学・疾患研究 | 2 | 精密腫瘍学 (CIViC/OncoKB)・疾患-遺伝子関連 (GWAS/Orphanet) |
| R. 量子・先端計算 | 5 | 量子計算・GNN・ベイズ統計・説明可能 AI・深層学習 |
| S. 医用イメージング | 1 | DICOM/NIfTI・WSI 病理画像・Radiomics・MONAI |
| T. シングルセル・空間・エピゲノミクス | 3 | scRNA-seq・Visium・MERFISH・CELLxGENE・RNA velocity・エピゲノミクス |
| U. 免疫・感染症 | 2 | 免疫情報学・MHC 結合予測・病原体ゲノミクス・AMR・IEDB |
| V. マイクロバイオーム・環境 | 2 | 16S/メタゲノム・α/β 多様性・SDM・OBIS・GBIF |
| W. システム生物学 | 1 | SBML シミュレーション・FBA・GRN 推定・BioModels |
| X. 疫学・公衆衛生 | 1 | リスク指標 (RR/OR)・年齢標準化・空間疫学・WHO・CDC |
| Y. 集団遺伝学 | 1 | HWE・PCA/ADMIXTURE・Fst・選択スキャン・gnomAD・GWAS |
| Z. 科学テキストマイニング | 1 | NER・関係抽出・知識グラフ・BERTopic・PubTator |

---

## Skills 一覧

### A. 基盤・ワークフロー（13 種）

全 Exp に共通する横断的な基盤スキル。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 1 | [scientific-pipeline-scaffold](scientific-pipeline-scaffold/SKILL.md) | パイプライン雛形・ディレクトリ構造・StepLogger・JSON 要約 | 全 Exp |
| 2 | [scientific-data-preprocessing](scientific-data-preprocessing/SKILL.md) | 欠損値補完・エンコーディング・スケーリング・外れ値処理 | 全 Exp |
| 3 | [scientific-data-simulation](scientific-data-simulation/SKILL.md) | 物理/化学/生物ベースの合成データ生成 | 06-09, 12, 13 |
| 4 | [scientific-publication-figures](scientific-publication-figures/SKILL.md) | 論文品質図表・rcParams・カラーパレット・マルチパネル | 10, 11-13 |
| 5 | [scientific-academic-writing](scientific-academic-writing/SKILL.md) | 学術論文執筆・ジャーナル別テンプレート・Cover Letter・査読対応 | 汎用 |
| 6 | [scientific-hypothesis-pipeline](scientific-hypothesis-pipeline/SKILL.md) | プロンプトから仮説立案・PICO/PECO 構造化・解析パイプライン自動生成 | 汎用 |
| 7 | [scientific-critical-review](scientific-critical-review/SKILL.md) | 草稿の批判的レビュー・考察深化・論理検証・修正案生成 | 汎用 |
| 8 | [scientific-supplementary-generator](scientific-supplementary-generator/SKILL.md) | Supplementary Information 自動生成・SI 図表整理・本文-SI 相互参照検証 | 汎用 |
| 9 | [scientific-latex-formatter](scientific-latex-formatter/SKILL.md) | Markdown→LaTeX 変換・ジャーナルテンプレート適用・BibTeX 生成 | 汎用 |
| 10 | [scientific-citation-checker](scientific-citation-checker/SKILL.md) | 引用文献の自動検索・網羅性チェック・整合性検証・重複検出 | 汎用 |
| 11 | [scientific-peer-review-response](scientific-peer-review-response/SKILL.md) | 査読コメント構造化・ポイントバイポイント回答・リバッタルレター生成 | 汎用 |
| 12 | [scientific-revision-tracker](scientific-revision-tracker/SKILL.md) | 改訂履歴追跡・差分管理・変更マークアップ・トレーサビリティ検証 | 汎用 |
| 13 | [scientific-paper-quality](scientific-paper-quality/SKILL.md) | 可読性スコア・構造バランス・語彙品質・ジャーナル適合性・再現性チェック | 汎用 |

### B. 統計・探索的解析（3 種）

データの理解・検定・次元削減を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 14 | [scientific-eda-correlation](scientific-eda-correlation/SKILL.md) | 探索的データ解析・相関ヒートマップ・分布可視化 | 02, 12, 13 |
| 15 | [scientific-statistical-testing](scientific-statistical-testing/SKILL.md) | 仮説検定・多重比較・エンリッチメント・ベイズ推論 | 03, 04, 06, 07 |
| 16 | [scientific-pca-tsne](scientific-pca-tsne/SKILL.md) | PCA / t-SNE / UMAP 次元削減・クラスタリング | 02, 03, 07, 11, 13 |

### C. 機械学習・モデリング（3 種）

教師あり学習と特徴量解釈を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 17 | [scientific-ml-regression](scientific-ml-regression/SKILL.md) | マルチターゲット回帰・モデル比較・レーダーチャート | 05, 12, 13 |
| 18 | [scientific-ml-classification](scientific-ml-classification/SKILL.md) | 分類 ML・ROC・PR 曲線・混同行列・PDP・Volcano | 03, 05 |
| 19 | [scientific-feature-importance](scientific-feature-importance/SKILL.md) | Tree-based & Permutation 特徴量重要度・PDP | 05, 12, 13 |

### D. 実験計画・プロセス最適化（2 種）

実験設計と応答曲面最適化を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 20 | [scientific-doe](scientific-doe/SKILL.md) | 田口直交表・CCD/Box-Behnken・ANOVA 因子効果・ベイズ最適化 | 汎用 |
| 21 | [scientific-process-optimization](scientific-process-optimization/SKILL.md) | 応答曲面法 (ML-RSM)・パレート最適化・プロセスウィンドウ | 12, 13 |

### E. 信号・スペクトル・時系列（4 種）

波形・周波数領域・神経電気生理学の解析を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 22 | [scientific-spectral-signal](scientific-spectral-signal/SKILL.md) | スペクトル前処理・フィルタリング・ピーク検出 | 11 |
| 23 | [scientific-biosignal-processing](scientific-biosignal-processing/SKILL.md) | ECG R波/HRV・EEG バンドパワー/ERP・EMG バースト・Poincaré | 08 |
| 24 | [scientific-time-series](scientific-time-series/SKILL.md) | STL 分解・SARIMA 予測・変化点検出・FFT 周期解析・Granger 因果 | 汎用 |
| 67 | [scientific-neuroscience-electrophysiology](scientific-neuroscience-electrophysiology/SKILL.md) | SpikeInterface/Kilosort4 スパイクソート・MNE EEG/ERP・NeuroKit2 HRV/EDA・脳機能結合 | 汎用 |

### F. 生命科学・オミクス（7 種）

バイオ・オミクス・ネットワーク解析を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 25 | [scientific-bioinformatics](scientific-bioinformatics/SKILL.md) | scRNA-seq・PPI ネットワーク・バルク RNA-seq | 01, 04 |
| 26 | [scientific-metabolomics](scientific-metabolomics/SKILL.md) | PLS-DA/VIP スコア・Pareto スケーリング・パスウェイ濃縮 | 07 |
| 27 | [scientific-sequence-analysis](scientific-sequence-analysis/SKILL.md) | RSCU/CAI コドン解析・アラインメント・系統樹・ORF/CpG 島 | 09 |
| 28 | [scientific-multi-omics](scientific-multi-omics/SKILL.md) | CCA 正準相関・SNF ネットワーク融合・パスウェイ統合・マルチオミクスクラスタ | 汎用 |
| 29 | [scientific-network-analysis](scientific-network-analysis/SKILL.md) | ネットワーク構築・中心性・コミュニティ・PSP パス図 | 04, 07, 13 |
| 68 | [scientific-proteomics-mass-spectrometry](scientific-proteomics-mass-spectrometry/SKILL.md) | pyOpenMS LC-MS/MS・ペプチド ID・タンパク質定量・PTM・GNPS 分子ネットワーク | 汎用 |
| 69 | [scientific-gene-expression-transcriptomics](scientific-gene-expression-transcriptomics/SKILL.md) | GEO データ取得・PyDESeq2 差次発現・GTEx 組織発現/eQTL・GSEA | 汎用 |

### G. 化学・材料・イメージング（4 種）

化学構造・材料特性評価・画像形態解析・計算材料科学を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 30 | [scientific-cheminformatics](scientific-cheminformatics/SKILL.md) | RDKit 分子記述子・Tanimoto・構造アラート・Lipinski | 02, 05 |
| 31 | [scientific-materials-characterization](scientific-materials-characterization/SKILL.md) | Thornton-Anders SZM・XRD Scherrer・Tauc プロット | 11, 12, 13 |
| 32 | [scientific-image-analysis](scientific-image-analysis/SKILL.md) | Otsu/Watershed セグメンテーション・粒径分布・GLCM テクスチャ・蛍光合成 | 汎用 |
| 70 | [scientific-computational-materials](scientific-computational-materials/SKILL.md) | pymatgen 結晶構造・Materials Project・相図・バンド構造/DOS・VASP/QE I/O | 汎用 |

### H. 臨床・疫学・メタ科学（4 種）

臨床試験・因果推論・メタアナリシス・臨床試験解析を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 33 | [scientific-survival-clinical](scientific-survival-clinical/SKILL.md) | Kaplan-Meier・Cox PH・検出力分析・安全性解析 | 03, 06 |
| 34 | [scientific-causal-inference](scientific-causal-inference/SKILL.md) | PSM 傾向スコア・IPW・DID・RDD・DAG 共変量選択・Rosenbaum 感度分析 | 汎用 |
| 35 | [scientific-meta-analysis](scientific-meta-analysis/SKILL.md) | 固定/ランダム効果モデル・Forest/Funnel プロット・Egger 検定・サブグループ | 汎用 |
| 71 | [scientific-clinical-trials-analytics](scientific-clinical-trials-analytics/SKILL.md) | ClinicalTrials.gov API v2 検索・競合ランドスケープ・AE/アウトカム抽出 | 汎用 |

### I. Deep Research（1 種）

科学文献の反復的深層リサーチを担うスキル。SHIKIGAMI の WebResearcher パラダイムを科学研究に適応。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 36 | [scientific-deep-research](scientific-deep-research/SKILL.md) | SHIKIGAMI 準拠 Think→Search→Evaluate→Synthesize 反復サイクル・学術 DB 検索・エビデンス階層評価・ソース追跡・交差検証・ハルシネーション防止 | 汎用 |

### J. 創薬・ファーマコロジー（3 種）

ドラッグディスカバリーの標的評価・薬物動態・リポジショニングを担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 37 | [scientific-drug-target-profiling](scientific-drug-target-profiling/SKILL.md) | 9-path 標的プロファイリング・TDL 分類・ドラッガビリティ評価・競合ランドスケープ | 汎用 |
| 38 | [scientific-admet-pharmacokinetics](scientific-admet-pharmacokinetics/SKILL.md) | 5 段階 ADMET パイプライン・Lipinski/Veber ルール・CYP 予測・PK モデリング | 汎用 |
| 39 | [scientific-drug-repurposing](scientific-drug-repurposing/SKILL.md) | 7 戦略ドラッグリポジショニング・ネットワーク近接解析・多基準候補スコアリング | 汎用 |

### K. 構造生物学・タンパク質工学（2 種）

タンパク質構造解析と de novo 設計を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 40 | [scientific-protein-structure-analysis](scientific-protein-structure-analysis/SKILL.md) | PDB/AlphaFold 構造検索・品質評価 (pLDDT/R-factor)・結合サイト検出 | 汎用 |
| 41 | [scientific-protein-design](scientific-protein-design/SKILL.md) | ESM-2 変異スキャン・RFdiffusion/ProteinMPNN de novo 設計・バインダー/酵素設計 | 汎用 |

### L. 精密医療・臨床意思決定（2 種）

バリアント解釈とエビデンスベース臨床判断を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 42 | [scientific-variant-interpretation](scientific-variant-interpretation/SKILL.md) | ACMG/AMP 28 基準・薬理ゲノミクス (CPIC)・OncoKB 体細胞変異レベル | 汎用 |
| 43 | [scientific-clinical-decision-support](scientific-clinical-decision-support/SKILL.md) | GRADE エビデンス枠組・精密腫瘍学ワークフロー・臨床試験マッチング | 汎用 |

### M. 実験室自動化・データ管理（2 種）

ラボ実験の自動化とデータ管理を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 44 | [scientific-lab-automation](scientific-lab-automation/SKILL.md) | PyLabRobot/Opentrons プロトコル・SOP テンプレート・ELN/LIMS 連携・QC 検証 | 汎用 |
| 72 | [scientific-lab-data-management](scientific-lab-data-management/SKILL.md) | Benchling ELN/DNA 設計・DNAnexus PaaS・OMERO バイオイメージング・Protocols.io | 汎用 |

### N. 科学プレゼンテーション・図式（2 種）

学会発表用スライド・ポスター・科学図式のデザインを担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 45 | [scientific-presentation-design](scientific-presentation-design/SKILL.md) | 15 スライド構成テンプレート・tikzposter・matplotlib ワークフロー図・アクセシビリティ | 汎用 |
| 73 | [scientific-scientific-schematics](scientific-scientific-schematics/SKILL.md) | CONSORT フロー図・NN アーキテクチャ図・パスウェイ図・TikZ/SVG | 汎用 |

### O. 研究計画・グラント・規制（3 種）

助成金申請書と研究方法論・規制科学の設計を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 46 | [scientific-grant-writing](scientific-grant-writing/SKILL.md) | NIH Specific Aims テンプレート・JSPS 科研費・予算計画・Budget Justification | 汎用 |
| 47 | [scientific-research-methodology](scientific-research-methodology/SKILL.md) | SCAMPER/TRIZ ブレインストーミング・研究デザインマトリクス・FINER 基準・IRB 倫理チェック | 汎用 |
| 74 | [scientific-regulatory-science](scientific-regulatory-science/SKILL.md) | FDA Orange Book/医療機器 510(k)/ISO 13485 QMS/CAPA/USPTO 特許検索 | 汎用 |

### P. ファーマコビジランス・薬理ゲノミクス（2 種）

市販後医薬品安全性監視と薬理ゲノミクスのためのシグナル検出・定量評価を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 48 | [scientific-pharmacovigilance](scientific-pharmacovigilance/SKILL.md) | FAERS 不均衡分析 (PRR/ROR/IC/EBGM)・MedDRA 階層・時系列トレンド・Naranjo 因果評価 | 汎用 |
| 75 | [scientific-pharmacogenomics](scientific-pharmacogenomics/SKILL.md) | PharmGKB/CPIC ガイドライン・Star アレル・代謝型・FDA PGx バイオマーカー | 汎用 |

### Q. 腫瘍学・疾患研究（2 種）

精密腫瘍学と疾患-遺伝子関連研究を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 49 | [scientific-precision-oncology](scientific-precision-oncology/SKILL.md) | CIViC/OncoKB/cBioPortal 統合・TMB/MSI 判定・AMP Tiering・MTB レポート | 汎用 |
| 50 | [scientific-disease-research](scientific-disease-research/SKILL.md) | GWAS Catalog・DisGeNET GDA・Orphanet/OMIM/HPO 表現型マッチング・PRS 算出 | 汎用 |

### R. 量子・先端計算（5 種）

量子計算・GNN・ベイズ統計・XAI・深層学習など次世代計算手法を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 51 | [scientific-quantum-computing](scientific-quantum-computing/SKILL.md) | Qiskit/Cirq/PennyLane VQE・QAOA・量子 ML・QuTiP 量子ダイナミクス | 汎用 |
| 52 | [scientific-graph-neural-networks](scientific-graph-neural-networks/SKILL.md) | PyG GCN/GAT/GIN・TorchDrug 分子特性予測・知識グラフ推論・Scaffold Split | 汎用 |
| 53 | [scientific-bayesian-statistics](scientific-bayesian-statistics/SKILL.md) | PyMC/Stan 階層ベイズ・MCMC 診断・PPC・WAIC/LOO-CV モデル比較・ベイズ最適化 | 汎用 |
| 54 | [scientific-explainable-ai](scientific-explainable-ai/SKILL.md) | SHAP/LIME/Captum 特徴量寄与・反実仮想説明・公平性監査・DeepSHAP | 汎用 |
| 55 | [scientific-deep-learning](scientific-deep-learning/SKILL.md) | Lightning/timm/Transformers・CNN/ViT/BERT Fine-tune・Optuna HPO・ONNX エクスポート | 汎用 |

### S. 医用イメージング（1 種）

DICOM・WSI 等の医用画像の解析・セグメンテーションを担うスキル。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 56 | [scientific-medical-imaging](scientific-medical-imaging/SKILL.md) | DICOM/NIfTI 処理・MONAI U-Net/SwinUNETR・WSI パッチ抽出・Radiomics・3D 可視化 | 汎用 |

### T. シングルセル・空間・エピゲノミクス（3 種）

scRNA-seq・空間トランスクリプトミクス・エピゲノミクスの解析パイプラインを担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 57 | [scientific-single-cell-genomics](scientific-single-cell-genomics/SKILL.md) | scRNA-seq QC・Scanpy Leiden クラスタリング・DEG・RNA velocity・CellChat 細胞間通信 | 汎用 |
| 58 | [scientific-spatial-transcriptomics](scientific-spatial-transcriptomics/SKILL.md) | Visium/MERFISH 前処理・Squidpy SVG 検出・空間ドメイン・cell2location デコンボリューション | 汎用 |
| 76 | [scientific-epigenomics-chromatin](scientific-epigenomics-chromatin/SKILL.md) | ChIP-seq MACS2/3・ATAC-seq・WGBS DMR・ChromHMM・Hi-C TAD・モチーフ濃縮 | 汎用 |

### U. 免疫・感染症（2 種）

免疫情報学・病原体ゲノミクスの解析パイプラインを担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 59 | [scientific-immunoinformatics](scientific-immunoinformatics/SKILL.md) | MHC-I/II 結合予測・B 細胞エピトープ・TCR/BCR レパトア多様性・抗体 CDR 解析・ワクチン候補ランキング | 汎用 |
| 60 | [scientific-infectious-disease](scientific-infectious-disease/SKILL.md) | 病原体 WGS QC・AMR 遺伝子検出・MLST 型別・系統解析 (IQ-TREE)・SIR/SEIR 数理モデル | 汎用 |

### V. マイクロバイオーム・環境（2 種）

マイクロバイオーム解析と環境・生態系モデリングを担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 61 | [scientific-microbiome-metagenomics](scientific-microbiome-metagenomics/SKILL.md) | DADA2 ASV パイプライン・MetaPhlAn/Kraken2・α/β 多様性・ANCOM-BC 差次的存在量・HUMAnN 機能プロファイリング | 汎用 |
| 62 | [scientific-environmental-ecology](scientific-environmental-ecology/SKILL.md) | SDM (MaxEnt/RF/GBM)・生物多様性指数・群集序列化 (NMDS/CCA)・保全優先順位ランキング | 汎用 |

### W. システム生物学（1 種）

SBML 動的シミュレーション・代謝フラックス・遺伝子制御ネットワーク推定を担うスキル。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 63 | [scientific-systems-biology](scientific-systems-biology/SKILL.md) | SBML/RoadRunner シミュレーション・FBA/pFBA (cobrapy)・GRN 推定 (GENIE3)・Sobol 感度解析 | 汎用 |

### X. 疫学・公衆衛生（1 種）

疫学的リスク指標算出・標準化・空間疫学・DAG 交絡分析を担うスキル。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 64 | [scientific-epidemiology-public-health](scientific-epidemiology-public-health/SKILL.md) | RR/OR/RD/NNT/AF リスク指標・直接/間接年齢標準化・LISA/Getis-Ord 空間クラスタリング・DAG バックドア基準 | 汎用 |

### Y. 集団遺伝学（1 種）

集団構造推定・分化指標・自然選択検出を担うスキル。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 65 | [scientific-population-genetics](scientific-population-genetics/SKILL.md) | PLINK2 QC・HWE 検定・PCA/ADMIXTURE・Weir-Cockerham Fst・iHS/Tajima's D 選択スキャン | 汎用 |

### Z. 科学テキストマイニング（1 種）

科学文献からの情報抽出・知識グラフ構築・トピックモデリングを担うスキル。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 66 | [scientific-text-mining-nlp](scientific-text-mining-nlp/SKILL.md) | BioBERT/SciSpaCy NER・関係抽出・知識グラフ構築 (Louvain)・BERTopic トピックモデリング・引用ネットワーク分析 | 汎用 |

---

## インストール

```bash
# npx でワンコマンドインストール
npx @nahisaho/satori init

# またはグローバルインストール
npm install -g @nahisaho/satori
satori init
```

`.github/skills/` がカレントディレクトリにコピーされ、Copilot Agent Mode で即座に利用できます。

---

## 使い方

### GitHub Copilot Agent Mode / Copilot CLI での利用

Skills は `.github/skills/` に配置されているため、Copilot が自動的に検出します。
プロンプトの文脈に応じて、関連する Skill の `SKILL.md` がエージェントに注入されます。

```
# 例：相関解析を依頼すると scientific-eda-correlation が自動ロードされる
> ZnO 薄膜のプロセスパラメータと膜物性の相関ヒートマップを作成して

# 例：分類モデルを依頼すると scientific-ml-classification が自動ロードされる
> がん遺伝子発現データで Random Forest と SVM の ROC 比較をして

# 例：実験計画法を依頼すると scientific-doe が自動ロードされる
> 3因子の田口L9直交表を作成して主効果プロットを描画して

# 例：メタアナリシスを依頼すると scientific-meta-analysis が自動ロードされる
> 5本のRCT論文からランダム効果モデルでForestプロットを作成して
```

### ディレクトリ構造

```
.github/skills/
├── README.md
│
│── [A] 基盤・ワークフロー
│   ├── scientific-pipeline-scaffold/
│   ├── scientific-data-preprocessing/
│   ├── scientific-data-simulation/
│   ├── scientific-publication-figures/
│   ├── scientific-academic-writing/
│   │   └── assets/   ← ジャーナル別テンプレート 7 種
│   ├── scientific-hypothesis-pipeline/
│   ├── scientific-critical-review/
│   ├── scientific-supplementary-generator/
│   ├── scientific-latex-formatter/
│   ├── scientific-citation-checker/
│   ├── scientific-peer-review-response/
│   ├── scientific-revision-tracker/
│   └── scientific-paper-quality/
│
│── [B] 統計・探索的解析
│   ├── scientific-eda-correlation/
│   ├── scientific-statistical-testing/
│   └── scientific-pca-tsne/
│
│── [C] 機械学習・モデリング
│   ├── scientific-ml-regression/
│   ├── scientific-ml-classification/
│   └── scientific-feature-importance/
│
│── [D] 実験計画・プロセス最適化
│   ├── scientific-doe/
│   └── scientific-process-optimization/
│
│── [E] 信号・スペクトル・時系列
│   ├── scientific-spectral-signal/
│   ├── scientific-biosignal-processing/
│   ├── scientific-time-series/
│   └── scientific-neuroscience-electrophysiology/
│
│── [F] 生命科学・オミクス
│   ├── scientific-bioinformatics/
│   ├── scientific-metabolomics/
│   ├── scientific-sequence-analysis/
│   ├── scientific-multi-omics/
│   ├── scientific-network-analysis/
│   ├── scientific-proteomics-mass-spectrometry/
│   └── scientific-gene-expression-transcriptomics/
│
│── [G] 化学・材料・イメージング
│   ├── scientific-cheminformatics/
│   ├── scientific-materials-characterization/
│   ├── scientific-image-analysis/
│   └── scientific-computational-materials/
│
├── [H] 臨床・疫学・メタ科学
│   ├── scientific-survival-clinical/
│   ├── scientific-causal-inference/
│   ├── scientific-meta-analysis/
│   └── scientific-clinical-trials-analytics/
│
├── [I] Deep Research
│   └── scientific-deep-research/
│
├── [J] 創薬・ファーマコロジー
│   ├── scientific-drug-target-profiling/
│   ├── scientific-admet-pharmacokinetics/
│   └── scientific-drug-repurposing/
│
├── [K] 構造生物学・タンパク質工学
│   ├── scientific-protein-structure-analysis/
│   └── scientific-protein-design/
│
├── [L] 精密医療・臨床意思決定
│   ├── scientific-variant-interpretation/
│   └── scientific-clinical-decision-support/
│
├── [M] 実験室自動化・データ管理
│   ├── scientific-lab-automation/
│   └── scientific-lab-data-management/
│
├── [N] 科学プレゼンテーション・図式
│   ├── scientific-presentation-design/
│   └── scientific-scientific-schematics/
│
└── [O] 研究計画・グラント・規制
    ├── scientific-grant-writing/
    ├── scientific-research-methodology/
    └── scientific-regulatory-science/
│
├── [P] ファーマコビジランス・薬理ゲノミクス
│   ├── scientific-pharmacovigilance/
│   └── scientific-pharmacogenomics/
│
├── [Q] 腫瘍学・疾患研究
│   ├── scientific-precision-oncology/
│   └── scientific-disease-research/
│
├── [R] 量子・先端計算
│   ├── scientific-quantum-computing/
│   ├── scientific-graph-neural-networks/
│   ├── scientific-bayesian-statistics/
│   ├── scientific-explainable-ai/
│   └── scientific-deep-learning/
│
└── [S] 医用イメージング
    └── scientific-medical-imaging/
│
│── [T] シングルセル・空間・エピゲノミクス
│   ├── scientific-single-cell-genomics/
│   ├── scientific-spatial-transcriptomics/
│   └── scientific-epigenomics-chromatin/
│
│── [U] 免疫・感染症
│   ├── scientific-immunoinformatics/
│   └── scientific-infectious-disease/
│
│── [V] マイクロバイオーム・環境
│   ├── scientific-microbiome-metagenomics/
│   └── scientific-environmental-ecology/
│
│── [W] システム生物学
│   └── scientific-systems-biology/
│
│── [X] 疫学・公衆衛生
│   └── scientific-epidemiology-public-health/
│
│── [Y] 集団遺伝学
│   └── scientific-population-genetics/
│
└── [Z] 科学テキストマイニング
    └── scientific-text-mining-nlp/
```

> 注: 実際のファイルシステム上ではすべてのスキルディレクトリは `.github/skills/` 直下にフラットに配置されています。上記の中区分グルーピングは論理的な分類です。

---

## 参考

- [SATORI 使い方ガイド](../../docs/qiita-satori-guide.md)
- [GitHub Copilot Agent Skills ドキュメント](https://docs.github.com/en/copilot/concepts/agents/about-agent-skills)
- [Agent Skills オープン標準](https://github.com/agentskills/agentskills)
