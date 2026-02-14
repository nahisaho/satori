# SATORI スキルギャップ分析

## 比較対象リポジトリ
- **mims-harvard/ToolUniverse** — 1000+ 科学ツール (8カテゴリ, 84 DB, 281 API)
- **K-Dense-AI/claude-scientific-skills** — 140+ 科学スキル (~17カテゴリ)

## 分析基準
- SATORI 既存47スキル (15カテゴリ A–O) を基準に、上記2リポから**未カバー領域**を抽出
- 既に十分カバー済みの機能は**除外** (例: UniProt連携, PubMed検索, Reactome/KEGGパスウェイ等)

---

## Tier 1: 高優先度ギャップ (新スキル推奨)

### 1. ファーマコビジランス / 医薬品安全性
| 項目 | 内容 |
|------|------|
| **ソース** | ToolUniverse (`tooluniverse-pharmacovigilance`) |
| **重要度** | **H** |
| **種別** | 🆕 スタンドアロン |
| **既存との重複** | `admet-pharmacokinetics` は前臨床 ADMET。市販後安全性は未カバー |
| **カバー内容** | FAERS有害事象報告、不均衡分析 (ROR/PRR/IC)、MedDRA階層集約、時系列トレンド、薬物比較、FDA添付文書警告、PharmGKB薬理ゲノミクス統合、ICD-10マッピング |
| **推奨スキル名** | `scientific-pharmacovigilance` |
| **カテゴリ** | J (創薬) に追加、or 新カテゴリ P |

### 2. 精密腫瘍学 (Precision Oncology)
| 項目 | 内容 |
|------|------|
| **ソース** | ToolUniverse (`tooluniverse-precision-oncology`) |
| **重要度** | **H** |
| **種別** | 🆕 スタンドアロン |
| **既存との重複** | `variant-interpretation` は一般的。がん特化ではない |
| **カバー内容** | CIViC (Clinical Interpretation of Variants in Cancer)、OncoKB治療アクション可能性、COSMIC体細胞変異、GDC/TCGA腫瘍データ、cBioPortalクロススタディ解析、DepMap標的バリデーション、HPA発現バリデーション |
| **推奨スキル名** | `scientific-precision-oncology` |
| **カテゴリ** | L (精密医療) に追加 |

### 3. 疾患研究パイプライン
| 項目 | 内容 |
|------|------|
| **ソース** | ToolUniverse (`tooluniverse-disease-research`) |
| **重要度** | **H** |
| **種別** | 🆕 スタンドアロン |
| **既存との重複** | 部分的 (`variant-interpretation`, `survival-clinical`, `network-analysis`) |
| **カバー内容** | OpenTargets疾患-標的マッピング、GWAS Catalog (13ツール) 統合、ClinVar病原性変異、GTEx組織発現、HPA蛋白発現、GEO公開データセット、疫学統計、治療オプション統合レポート生成 |
| **推奨スキル名** | `scientific-disease-research` |
| **カテゴリ** | H (臨床・疫学) に追加 |

### 4. 量子コンピューティング
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (4スキル: `qiskit`, `cirq`, `pennylane`, `qutip`) |
| **重要度** | **H** |
| **種別** | 🆕 新カテゴリ |
| **既存との重複** | **なし** — SATORI に量子計算の機能は皆無 |
| **カバー内容** | IBM Qiskit (量子回路構築、Estimator/Sampler、VQE)、Google Cirq (ノイズモデリング、パラメータスイープ)、PennyLane (量子ML、自動微分、ハイブリッド量子-古典モデル)、QuTiP (開放量子系、密度行列、マスター方程式)。量子化学 (分子ハミルトニアン、基底状態エネルギー計算) を含む |
| **推奨スキル名** | `scientific-quantum-computing` (統合) or 個別4スキル |
| **カテゴリ** | 🆕 新カテゴリ P: 量子計算・シミュレーション |

### 5. グラフニューラルネットワーク (GNN)
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`torch-geometric`) |
| **重要度** | **H** |
| **種別** | 🆕 スタンドアロン |
| **既存との重複** | `network-analysis` はグラフ理論。GNN は未カバー |
| **カバー内容** | GCN/GAT/GraphSAGE/GIN 層、分子グラフ畳み込み (SchNet, DimeNet)、ヘテロジニアスグラフ、グラフプーリング、GNN Explainability、PyTorch Geometric データセット (ZINC, MoleculeNet 等) |
| **推奨スキル名** | `scientific-graph-neural-networks` |
| **カテゴリ** | C (機械学習) に追加 |

### 6. 医用画像・デジタル病理学
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (3スキル: `pydicom`, `histolab`, `pathml`) |
| **重要度** | **H** |
| **種別** | 🆕 スタンドアロン |
| **既存との重複** | `image-analysis` は一般画像処理。DICOM/WSI/病理特化ではない |
| **カバー内容** | DICOM読み書き・メタデータ処理 (pydicom)、WSI組織タイル抽出 (histolab)、デジタル病理パイプライン — 細胞セグメンテーション、組織分類、細胞グラフ構築、空間解析、PyG統合 (PathML) |
| **推奨スキル名** | `scientific-medical-imaging` |
| **カテゴリ** | 🆕 新カテゴリ Q: 医用画像、or G (イメージング) 拡張 |

### 7. ベイズ統計モデリング
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`pymc`) |
| **重要度** | **H** |
| **種別** | 🆕 スタンドアロン or `statistical-testing` 大幅拡張 |
| **既存との重複** | `statistical-testing` は頻度主義。ベイズ推論は未カバー |
| **カバー内容** | 確率的プログラミング、MCMC (NUTS, Metropolis)、変分推論、階層モデル、ベイジアンネットワーク、事後分布推定、収束診断、モデル比較 (WAIC, LOO) |
| **推奨スキル名** | `scientific-bayesian-modeling` |
| **カテゴリ** | B (統計) に追加 |

### 8. 説明可能AI (XAI)
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`shap`) |
| **重要度** | **H** |
| **種別** | 🔄 `feature-importance` の大幅拡張 or 🆕 スタンドアロン |
| **既存との重複** | `feature-importance` は基本的な特徴量重要度。SHAP方法論全体は未カバー |
| **カバー内容** | TreeExplainer (XGBoost/LightGBM/RF)、DeepExplainer (TF/PyTorch)、KernelExplainer (モデル非依存)、Waterfall/Beeswarm/Bar/Scatter/Force/Heatmapプロット、公平性分析、モデル比較、production XAI |
| **推奨スキル名** | `scientific-explainable-ai` |
| **カテゴリ** | C (機械学習) に追加 |

### 9. 深層学習フレームワーク
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`pytorch-lightning`, `transformers`) |
| **重要度** | **H** |
| **種別** | 🆕 スタンドアロン |
| **既存との重複** | `ml-regression`/`ml-classification` は scikit-learn レベル。DL は未カバー |
| **カバー内容** | PyTorch Lightning (訓練ループ抽象化、分散訓練、混合精度)、HuggingFace Transformers (NLP/CV/Audio、ファインチューニング、tokenization)、Transfer Learning |
| **推奨スキル名** | `scientific-deep-learning` |
| **カテゴリ** | C (機械学習) に追加 |

---

## Tier 2: 中優先度ギャップ (拡張 or 条件付き新スキル)

### 10. プロテオミクス・質量分析
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`matchms`, `pyopenms`) |
| **重要度** | **M** |
| **種別** | 🆕 スタンドアロン |
| **既存との重複** | `spectral-signal` は一般スペクトル、`metabolomics` は代謝物。MS特化は未カバー |
| **カバー内容** | マススペクトルマッチング、スペクトルライブラリ検索、ペプチド/蛋白同定、定量プロテオミクス、pyOpenMS LC-MS/MSワークフロー |
| **推奨スキル名** | `scientific-proteomics-ms` |
| **カテゴリ** | F (生命科学) に追加 |

### 11. 代謝モデリング (制約ベース再構築)
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`cobrapy`) |
| **重要度** | **M** |
| **種別** | 🔄 `metabolomics` or `network-analysis` の拡張 |
| **既存との重複** | `metabolomics` は代謝物プロファイリング。FBA/制約ベースモデリングは別領域 |
| **カバー内容** | フラックスバランス解析 (FBA)、フラックス変動解析 (FVA)、最小培地予測、遺伝子ノックアウトシミュレーション、代謝ネットワーク可視化 |
| **推奨スキル名** | `scientific-metabolic-modeling` or `metabolomics` 拡張 |

### 12. 神経科学ツール
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`neurokit2`, `neuropixels-analysis`) |
| **重要度** | **M** |
| **種別** | 🔄 `biosignal-processing` の拡張 |
| **既存との重複** | `biosignal-processing` は一般生体信号。神経科学特化ツールは不足 |
| **カバー内容** | NeuroKit2 (ECG/EDA/EMG/EEG解析、HRV)、Neuropixels (大規模神経記録、スパイクソーティング、アセンブリ検出)、脳結合性解析 |
| **推奨スキル名** | `biosignal-processing` 内に統合 or `scientific-neuroscience` |

### 13. 多目的最適化
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`pymoo`) |
| **重要度** | **M** |
| **種別** | 🔄 `process-optimization` or `doe` の拡張 |
| **既存との重複** | `process-optimization` は単目的が主。パレート最適化は未カバー |
| **カバー内容** | NSGA-II/NSGA-III、パレートフロント、多目的進化アルゴリズム、制約付き最適化、分解ベース手法 (MOEA/D) |
| **推奨スキル名** | `process-optimization` 拡張 |

### 14. 離散事象シミュレーション
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`simpy`) |
| **重要度** | **M** |
| **種別** | 🆕 スタンドアロン or `process-optimization` 拡張 |
| **既存との重複** | プロセスシミュレーション特化は未カバー |
| **カバー内容** | 待ち行列モデリング、リソーススケジューリング、臨床試験シミュレーション、製造プロセスシミュレーション、確率的イベント駆動モデル |
| **推奨スキル名** | `scientific-process-simulation` |

### 15. 大規模データ処理
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`dask`, `polars`, `vaex`) |
| **重要度** | **M** |
| **種別** | 🔄 `data-preprocessing` の拡張 |
| **既存との重複** | `data-preprocessing` は pandas レベルか。分散/Out-of-core は未確認 |
| **カバー内容** | Dask 分散計算、Polars 高速DataFrame、Vaex Out-of-core処理、並列ETLパイプライン |
| **推奨スキル名** | `data-preprocessing` 内に統合 |

### 16. 科学スキマティクス生成
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`scientific-schematics`) |
| **重要度** | **M** |
| **種別** | 🔄 `publication-figures` or `presentation-design` の拡張 |
| **既存との重複** | `publication-figures` はmatplotlib/seabornプロット。概念図・フローチャートは別 |
| **カバー内容** | CONSORT フローチャート、神経ネットワークアーキテクチャ図、生物学的パスウェイ図、システムアーキテクチャ図、TikZ/SVG生成、AI生成モード (Nano Banana Pro + Gemini Review) |
| **推奨スキル名** | `scientific-schematics` or `publication-figures` 拡張 |

### 17. ELN/LIMS 統合
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`benchling-integration`, `labarchive-integration`, `protocolsio-integration`) |
| **重要度** | **M** |
| **種別** | 🔄 `lab-automation` の拡張 |
| **既存との重複** | `lab-automation` はロボティクス操作。ELN/LIMSは別領域 |
| **カバー内容** | Benchling API (DNA配列、レジストリ)、LabArchives (ノートブック管理)、Protocols.io (プロトコル共有・バージョン管理) |
| **推奨スキル名** | `lab-automation` 拡張 or `scientific-eln-integration` |

### 18. クラウド科学計算基盤
| 項目 | 内容 |
|------|------|
| **ソース** | K-Dense (`modal`, `dnanexus`, `latchbio-integration`, `omero`) |
| **重要度** | **M** |
| **種別** | 🔄 `pipeline-scaffold` の拡張 |
| **既存との重複** | `pipeline-scaffold` はローカルパイプライン。クラウド基盤は未カバー |
| **カバー内容** | Modal (サーバーレスGPU計算)、DNAnexus (ゲノミクスプラットフォーム)、LatchBio (バイオインフォパイプライン — bulk RNA-seq, DESeq2, AlphaFold等のVerified Workflows)、OMERO (バイオイメージングデータ管理) |
| **推奨スキル名** | `scientific-cloud-infrastructure` |

---

## Tier 3: 低優先度ギャップ (ニッチ or 強い重複)

| # | 機能 | ソース | 重要度 | 種別 | 既存との重複 | 備考 |
|---|------|--------|--------|------|------------|------|
| 19 | **天文学** (Astropy) | K-Dense | L | 新規 | なし | 天文ユーザー対象時のみ |
| 20 | **CFD** (FluidSim) | K-Dense | L | 新規 | なし | 工学ユーザー対象時のみ |
| 21 | **地理空間解析** (GeoPandas) | K-Dense | L | 新規 | なし | 環境科学ユーザー対象時のみ |
| 22 | **規制・品質管理** (ISO 13485) | K-Dense | L | 新規 | なし | 医療機器業界向けニッチ |
| 23 | **フローサイトメトリー** (FlowIO) | K-Dense | L | 拡張 | `biosignal-processing` 部分的 | 免疫学特化 |
| 24 | **強化学習** (SB3, PufferLib) | K-Dense | L | 新規 | なし | 科学応用は限定的 |
| 25 | **生存分析ライブラリ** (scikit-survival) | K-Dense | L | 拡張 | `survival-clinical` **強い重複** | 既存で十分な可能性 |
| 26 | **時系列ML** (aeon) | K-Dense | L | 拡張 | `time-series` **強い重複** | 既存で十分な可能性 |
| 27 | **ポスターデザイン** (LaTeX/PPTX) | K-Dense | L | 拡張 | `presentation-design` 部分的 | 低付加価値 |
| 28 | **UMAP次元削減** (umap-learn) | K-Dense | L | 拡張 | `pca-tsne` **強い重複** | tSNE → UMAP拡張のみ |
| 29 | **シンボリック数学** (SymPy) | K-Dense | L | 拡張 | `statistical-testing` 部分的 | 数式処理ニーズ限定 |
| 30 | **GRN推論** (Arboreto) | K-Dense | L | 拡張 | `network-analysis` 部分的 | GENIE3/GRNBoost2特化 |

---

## サマリー: 推奨アクション優先順位

### Phase 1: 即時追加推奨 (9 スキル)

| 順位 | 推奨スキル名 | カテゴリ | 根拠 |
|------|-------------|---------|------|
| 1 | `scientific-pharmacovigilance` | J→拡張 | ToolUniverse に完成度の高いスキルが存在。FAERS/FDA統合は創薬サイクル完結に必須 |
| 2 | `scientific-precision-oncology` | L→拡張 | がん研究は最大の研究領域。CIViC/OncoKB/cBioPortal統合は高価値 |
| 3 | `scientific-disease-research` | H→拡張 | GWAS Catalog + OpenTargets + ClinVar の疾患横断パイプラインは汎用性が高い |
| 4 | `scientific-quantum-computing` | 🆕 P | 量子計算は科学計算の次世代基盤。VQE/量子化学は材料・創薬に直結 |
| 5 | `scientific-graph-neural-networks` | C→拡張 | 分子性質予測、蛋白相互作用、知識グラフ推論に不可欠 |
| 6 | `scientific-medical-imaging` | G→拡張 | DICOM/WSI/病理AIは臨床研究の急成長分野 |
| 7 | `scientific-bayesian-modeling` | B→拡張 | ベイズ推論は現代統計の主流手法。頻度主義のみでは不十分 |
| 8 | `scientific-explainable-ai` | C→拡張 | 規制要件としてのXAIニーズ急増。SHAP は事実上の標準 |
| 9 | `scientific-deep-learning` | C→拡張 | Transformers/Lightning は現代MLの基盤。scikit-learnのみでは不足 |

### Phase 2: 計画的拡張推奨 (9 スキル/拡張)

| 順位 | 推奨アクション | 対象既存スキル | 根拠 |
|------|---------------|---------------|------|
| 10 | `scientific-proteomics-ms` 新規 | — | オミクス解析の重要な欠落ピース |
| 11 | `scientific-metabolic-modeling` or 拡張 | `metabolomics` | COBRApy/FBA は系統生物学の必須ツール |
| 12 | 神経科学ツール統合 | `biosignal-processing` | NeuroKit2/Neuropixels 追加 |
| 13 | 多目的最適化追加 | `process-optimization` | PyMOO パレート最適化 |
| 14 | 離散事象シミュレーション | `process-optimization` or 新規 | SimPy で臨床試験/製造プロセスモデリング |
| 15 | 大規模データ処理拡張 | `data-preprocessing` | Dask/Polars/Vaex 統合 |
| 16 | スキマティクス生成 | `publication-figures` | CONSORT/パスウェイ/アーキテクチャ図 |
| 17 | ELN/LIMS統合 | `lab-automation` | Benchling/LabArchives/Protocols.io |
| 18 | クラウド基盤統合 | `pipeline-scaffold` | Modal/DNAnexus/LatchBio |

---

## カバレッジマトリクス

```
                          SATORI    ToolUniverse    K-Dense
─────────────────────────────────────────────────────────────
前臨床 ADMET              ✅         ✅              ✅
市販後安全性(PV)          ❌         ✅✅            ❌
精密腫瘍学                ❌         ✅✅            ❌
疾患研究パイプライン      ❌         ✅✅            ❌
量子コンピューティング    ❌         ❌              ✅✅
GNN                       ❌         ❌              ✅✅
医用画像/病理             △          ❌              ✅✅
ベイズ統計                ❌         ❌              ✅✅
XAI (SHAP等)              △          ❌              ✅✅
深層学習フレームワーク    ❌         ❌              ✅✅
プロテオミクス/MS         ❌         ❌              ✅
代謝モデリング(FBA)       ❌         ❌              ✅
神経科学特化              △          ❌              ✅
多目的最適化              ❌         ❌              ✅
バイオインフォマティクス  ✅         ✅              ✅
分子動力学/構造           ✅         ✅              ✅
ネットワーク解析          ✅         ✅              ✅
臨床試験/メタ分析         ✅         ✅              ✅
学術文書/論文             ✅         ❌              ✅
実験計画(DOE)             ✅         ❌              ❌
因果推論                  ✅         ❌              ❌
─────────────────────────────────────────────────────────────
✅✅ = 非常に充実  ✅ = カバー  △ = 部分的  ❌ = 未カバー
```

---

## 結論

SATORI は学術文書・統計・実験計画・因果推論で独自の強みを持つが、以下3つの大きなギャップ領域が存在:

1. **臨床/トランスレーショナル領域** — ファーマコビジランス、精密腫瘍学、疾患研究パイプライン (ToolUniverse が先行)
2. **次世代計算手法** — 量子計算、GNN、深層学習、ベイズ推論、XAI (K-Dense が先行)
3. **専門イメージング/オミクス** — 医用画像・病理AI、プロテオミクス/MS (K-Dense が先行)

Phase 1 の 9 スキルを追加することで、既存47 → 56 スキルとなり、2リポとの主要ギャップを解消できる。
