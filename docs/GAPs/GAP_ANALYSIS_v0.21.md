# GAP ANALYSIS — SATORI v0.21.0 → v0.22.0+

## Current State (v0.21.0)

| Metric | Value |
|--------|-------|
| Total Skills | 160 |
| TU-linked Skills | 124 |
| Non-TU Skills | 36 |
| Categories | 26 (A-Z) |
| TU Coverage | 77.5% |

### Category Distribution (v0.21.0)

| Cat | Name | Skills | TU-linked |
|-----|------|--------|-----------|
| A | 基盤・ワークフロー | 17 | 4 |
| B | 統計・探索的解析 | 4 | 0 |
| C | 機械学習・モデリング | 3 | 0 |
| D | 実験計画・プロセス最適化 | 2 | 0 |
| E | 信号・スペクトル・時系列 | 4 | 0 |
| F | 生命科学・オミクス | 26 | 21 |
| G | 化学・材料・イメージング | 9 | 4 |
| H | 臨床・疫学・メタ科学 | 6 | 1 |
| I | Deep Research・文献検索 | 4 | 3 |
| J | 創薬・ファーマコロジー | 9 | 8 |
| K | 構造生物学・タンパク質工学 | 7 | 6 |
| L | 精密医療・臨床意思決定 | 6 | 5 |
| M | 実験室自動化・データ管理 | 2 | 0 |
| N | 科学プレゼンテーション・図式 | 2 | 0 |
| O | 研究計画・グラント・規制 | 3 | 1 |
| P | ファーマコビジランス・薬理ゲノミクス | 3 | 3 |
| Q | 腫瘍学・疾患研究 | 10 | 9 |
| R | 量子・先端計算 | 7 | 0 |
| S | 医用イメージング | 1 | 0 |
| T | シングルセル・空間・エピゲノミクス | 13 | 5 |
| U | 免疫・感染症 | 2 | 1 |
| V | マイクロバイオーム・環境 | 8 | 4 |
| W | システム生物学 | 4 | 2 |
| X | 疫学・公衆衛生 | 3 | 1 |
| Y | 集団遺伝学 | 2 | 1 |
| Z | 科学テキストマイニング | 3 | 1 |

---

## Phase 13 (v0.21.0) Summary — Completed

### Track A: 既存 8 スキルに TU key 追加 (+7 newly TU-linked)
1. human-cell-atlas → cellxgene_census (既存 hca_tools に追加、カウント変更なし)
2. drug-repurposing → pharos (+1)
3. pharmacogenomics → fda_pharmacogenomic_biomarkers (+1)
4. human-protein-atlas → hpa (+1)
5. variant-effect-prediction → spliceai, cadd (+1)
6. gtex-tissue-expression → gtex_v2 (+1)
7. biothings-idmapping → biothings (+1)
8. protein-structure-analysis → proteinsplus (+1)

### Track B: 新規 6 スキル (#155-#160、うち 3 件 TU 連携)
- #155 scientific-cellxgene-census (T) → cellxgene_census
- #156 scientific-pharos-targets (J) → pharos
- #157 scientific-clingen-curation (L) → clingen
- #158 scientific-clinical-nlp (Z) → no TU
- #159 scientific-hgnc-nomenclature (F) → no TU
- #160 scientific-metabolomics-network (F) → no TU

---

## Phase 14 Roadmap (v0.22.0 候補)

### Priority A: 既存スキルへの TU key 追加候補

| Skill | Candidate TU Key | Tools | Notes |
|-------|------------------|-------|-------|
| scientific-epigenomics-chromatin | chipatlas | 4 tools | ChIP-Atlas エンリッチメント (encode-screen にも連携) |
| scientific-metabolomics | metabolomics_workbench | 6 tools | MWB REST API メタボロームデータ |
| scientific-systems-biology | biomodels | existing | BioModels SBML (metabolic-modeling と共有) |
| scientific-infectious-disease | card | TBD | CARD AMR データベース (要検証) |
| scientific-environmental-ecology | gbif | TBD | GBIF 生物多様性 (marine-ecology と共有) |

### Priority B: 新規スキル候補

| # | Candidate Skill | Category | TU Key | Description |
|---|----------------|----------|--------|-------------|
| 161 | scientific-proteomics-atlas | F | hpa, proteinsplus | プロテオミクスアトラス統合 (HPA+ProteinsPlus+構造) |
| 162 | scientific-glycomics | F | None | 糖鎖構造解析・GlyConnect/GlyGen |
| 163 | scientific-lipidomics | F | None | リピドミクス・LipidMAPS/SwissLipids |
| 164 | scientific-metagenome-assembled-genomes | V | mgnify | MAG 解析・ビニング・品質評価 |
| 165 | scientific-crispr-design | M | None | CRISPR gRNA 設計・オフターゲット予測 |
| 166 | scientific-clinical-pharmacology | P | None | 臨床薬理学モデリング・PopPK/PBPK |

### Priority C: パイプライン強化

| Pipeline | Current | Proposed Enhancement |
|----------|---------|---------------------|
| Monarch → GDC → precision-oncology | 3-step | + pharos-targets → drug-repurposing (5-step) |
| STITCH → STRING → drug-target | 3-step | + clingen-curation → clinical-decision (5-step) |
| human-cell-atlas → cellxgene-census → single-cell | 3-step | 完成 (v0.21.0) |
| biothings-idmapping → hgnc-nomenclature → genome-sequence | 3-step | 完成 (v0.21.0) |
| metabolomics → metabolomics-network → pathway-enrichment | 3-step | 完成 (v0.21.0) |
| text-mining-nlp → clinical-nlp → clinical-reporting | 3-step | 完成 (v0.21.0) |

---

## Remaining Non-TU Skills (36)

以下は ToolUniverse に対応する TU key が存在しない/未確認のスキル:

### 計算・ML 系 (B-E, R) — 20 スキル
TU 連携不要（Python ライブラリベース）: statistical-testing, eda-correlation, pca-tsne, ml-regression, ml-classification, feature-importance, doe, process-optimization, spectral-signal, biosignal-processing, time-series, neuroscience-electrophysiology, quantum-computing, graph-neural-networks, bayesian-statistics, explainable-ai, deep-learning, healthcare-ai, reinforcement-learning, symbolic-mathematics

### 実験室・プレゼン・研究計画 (M-N) — 4 スキル
TU 連携不要: lab-automation, lab-data-management, presentation-design, scientific-schematics

### ライティング・基盤 (A) — 9 スキル
TU 連携不要: pipeline-scaffold, data-preprocessing, data-simulation, publication-figures, academic-writing, hypothesis-pipeline, critical-review, supplementary-generator, latex-formatter, citation-checker, peer-review-response, revision-tracker, paper-quality (一部は crossref-metadata 等で TU 連携済み)

### 新規 Non-TU (v0.21.0) — 3 スキル
clinical-nlp, hgnc-nomenclature, metabolomics-network

---

## Milestone Projection

| Version | Skills | TU-linked | Coverage | Phase |
|---------|--------|-----------|----------|-------|
| v0.15.0 | 120 | 65 | 54.2% | 7 |
| v0.16.0 | 124 | 69 | 55.6% | 8 |
| v0.17.0 | 129 | 74 | 57.4% | 9 |
| v0.18.0 | 135 | 80 | 59.3% | 10 |
| v0.19.0 | 148 | 99 | 66.9% | 11 |
| v0.20.0 | 154 | 114 | 74.0% | 12 |
| **v0.21.0** | **160** | **124** | **77.5%** | **13** |
| v0.22.0 (proj) | ~166 | ~130 | ~78.3% | 14 |
| v0.25.0 (goal) | ~180 | ~145 | ~80.6% | — |
