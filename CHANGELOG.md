# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.27.1] - 2026-02-15

### Changed
- **ドキュメント統一**: v0.26.0 → v0.27.0 対応
  - `docs/SATORI_REVERSE_INDEX.md`: パイプライン 26 → 50、TU キー統計更新
  - `docs/SATORI_PIPELINE_EXAMPLES.md`: タイトル・バージョン更新
  - `docs/qiita/SATORI_PIPELINE_EXAMPLES_QIITA.md`: v0.26.0 → v0.27.0、50 パイプライン記載
  - `docs/qiita/SATORI_REVERSE_INDEX_QIITA.md`: v0.26.0 → v0.27.0、50 パイプライン記載

### Documentation
- v0.27.0 の新機能（skill search/info、50 パイプライン）をドキュメント全体で反映
- Qiita 記事をバージョンアップ（日本語テクニカルライティング）

## [0.27.0] - 2026-02-15

### Added
- **`satori skill search` コマンド**: キーワードでスキルを検索
  - 対話的検索: 名前・説明・TU キーでマッチング
  -複数キーワード対応（スペース区切り）
  - トップ 10 件を表示し、詳細は `satori skill info` で確認可能
  - 使用例: `satori skill search 創薬`, `satori skill search machine learning`

- **`satori skill info` コマンド**: スキルの詳細情報を表示
  - スキル説明（説明フィールドから抽出）
  - When to Use セクション
  - ToolUniverse 連携ツール（tu_tools セクション）
  - 関連パイプライン（最大 5 件）
  - ファイル位置
  - 使用例: `satori skill info deep-learning`, `satori skill info drug-target-profiling`

- **パイプラインカタログの大幅拡張**: 26 → 50 パイプライン
  - 既存ドメインパイプライン: 26 
  - クロスドメインパイプライン: 15 (A-O)
    - ゲノム創薬統合（A）
    - AI 駆動臨床意思決定（B）
    - 研究自動化（C）
    - マルチオミクス疾患解明（D）
    - 個別化薬物療法（E）
    - バイオインフォマティクス完全（F）
    - がん精密医療 End-to-End（G）
    - マルチオミクス縦断統合（H）
    - 環境メタボ・マイクロバイオーム One Health（I）
    - AI 駆動マテリアルズインフォマティクス（J）
    - 研究ライフサイクル完全自動化（K）
    - AI 駆動エビデンス合成（L）
    - がんマルチレイヤーゲノム創薬（M）
    - 臨床→規制→出版バリューチェーン（N）
    - シングルセルプロテオーム統合（O）
  - インダストリーパイプライン: 5 (Ind-1～Ind-5)
    - 製薬企業レギュラトリー
    - 農業バイオテクノロジー
    - 臨床検査室ワークフロー
    - 食品安全・毒性評価
    - 法医・公衆衛生
  - メソドロジーパイプライン: 4 (M-α～M-δ)
    - ベイズ推論ワークフロー
    - 因果推論パイプライン
    - 時系列予測パイプライン
    - テキストマイニング・NLP

- **GitHub コミュニティステンプレート**
  - `.github/ISSUE_TEMPLATE/bug_report.md`: バグ報告テンプレート
  - `.github/ISSUE_TEMPLATE/feature_request.md`: 機能リクエストテンプレート
  - `.github/PULL_REQUEST_TEMPLATE.md`: プルリクエストテンプレート

- **CODE_OF_CONDUCT.md**: Contributor Covenant ベースの行動規範
  - コミュニティの価値観（敬意・包括性・透明性・協力）
  - 許容されない行為の明示
  - 問題報告と調査プロセス
  - 処分の段階

### Changed
- **`satori pipeline list` コマンド**: 50 パイプラインを分類表示
  - ドメインパイプライン (26) セクション
  - クロスドメインパイプライン (15) セクション
  - インダストリーパイプライン (5) セクション
  - メソドロジーパイプライン (4) セクション
  - 各セクション内でパイプラインを整理

- **パイプライン推奨システム**: 拡張キーワード＆スコアリング
  - クロスドメイン・インダストリー・メソドロジーパイプラインも対象
  - `satori pipeline suggest` でインタラクティブ推奨が改善

- **ヘルプテキスト更新**: `skill search` / `skill info` コマンド追加

- **統計レポート更新**
  - パイプライン数: 26 → 50
  - ユニーク TU キー: 325 (v0.26.0 より増加)

### Fixed
- テスト: `satori stats` でパイプライン数が 50 と表示されることを確認
- テスト: `satori init --dry-run` のファイル数マッチング正規表現改善

### Testing
- 12 新規 unit テスト追加（`skill search` 4 件、`skill info` 6 件、エッジケース 2 件）
- 統合テスト更新: パイプライン カウント検証を 26 → 50 に更新
- 全テスト: 1,382 テスト (1,383 → 1,382 テスト実行、36 unit + 11 integration + 1,333 validation + 2 tu-coverage)

## [0.26.0] - 2026-02-15

### Added
- **テストフレームワーク導入**: Vitest テストフレームワークを導入（devDependencies）
  - `vitest.config.js`: テスト設定ファイル（unit / integration / validation の 3 カテゴリ）
  - `npm test`: 全テスト実行、`npm run test:unit` / `npm run test:validation` で個別実行
- **CLI ユニットテスト** (`tests/unit/cli.test.js`): 24 テストケース
  - `--version` / `-v`: バージョン出力検証
  - `help` / `--help` / `-h`: ヘルプ表示検証
  - `init --dry-run`: ドライラン動作検証
  - `init`: 実際のインストール・`--force` 上書き検証
  - `pipeline list`: 26 パイプライン一覧表示検証
  - `pipeline (不正サブコマンド)`: エラーハンドリング検証
  - `validate`: 190 スキル検証 / `--verbose` 出力検証（3 テスト）
  - `stats`: 統計表示・スキル数・パイプライン数・TU カバレッジ・バージョン検証（5 テスト）
  - 不明コマンド・引数なし: エッジケース検証
- **SKILL.md フォーマット検証テスト** (`tests/validation/skill-format.test.js`): 190 × 7 = 1,330 テストケース
  - 全 190 スキルの YAML Frontmatter（name + description）整合性検証
  - `## When to Use` / `## Quick Start` セクション存在検証
  - コードブロック（python/markdown/json）存在検証
  - Frontmatter name ↔ ディレクトリ名一致検証
  - TU 参照統計（80+ スキルが TU 連携）
- **TU カバレッジ分析テスト** (`tests/validation/tu-coverage.test.js`)
  - 全 190 スキルの ToolUniverse 連携状況を自動検出・レポート生成
  - カテゴリ別 TU 未連携スキル一覧出力
  - カバレッジベースライン: 54.2%（自動検出ベース）
- **GitHub Actions CI** (`.github/workflows/ci.yml`)
  - Node.js 18/20/22 マトリクスでの全テスト実行
  - SKILL.md フォーマット検証ジョブ
  - CLI スモークテストジョブ（version / help / init --dry-run / pipeline list / validate / stats）
  - Biome lint チェックジョブ
- **`engines` フィールド追加**: `"node": ">=18"` を package.json に明記
- **`satori validate` コマンド**: 全 190 SKILL.md の自動検証（Frontmatter / セクション / コードブロック）
  - `--verbose` オプション対応
  - パス/フェイル数の集計表示
  - 不合格スキルの問題点一覧出力
- **`satori stats` コマンド**: プロジェクト統計ダッシュボード
  - スキル総数 / パイプライン数 / TU 連携率 / ユニーク TU キー数 / コードブロック数
- **Biome v2.3 リンター/フォーマッター導入** (`biome.json`)
  - `@biomejs/biome` を devDependencies に追加
  - `npm run lint` / `npm run lint:fix` / `npm run format` スクリプト追加
  - recommended ルールセット + `useConst` / `noUnusedVariables` / `noUnusedImports` 有効化
  - シングルクォート / トレイリングカンマ / セミコロン / 120 桁幅

### Fixed
- `scientific-critical-review/SKILL.md`: 1 行目の不正な ` ```skill ` ラッパーを除去（Frontmatter 解析不能の原因）
- `scientific-hypothesis-pipeline/SKILL.md`: 同上の不正ラッパーを除去

### Changed
- `package.json`: `scripts.test` を `"node bin/satori.js init --dry-run"` → `"vitest run"` に変更
- `package.json`: 10 の npm scripts 追加（test:watch / test:unit / test:validation / test:coverage / validate:skills / validate:tu / lint / lint:fix / format）
- `.gitignore`: `coverage/` / OS 一時ファイル / IDE 設定ファイルの除外ルール追加
- **README 更新**: CI バッジ追加、CLI コマンド一覧テーブル、開発セクション（テスト・リント・CI 構成）追加
- **CONTRIBUTING.md 新設**: コントリビューションガイド（SKILL.md フォーマット要件・テスト構成・コードスタイル・コミット規約・CI 説明）
- **統合テスト** (`tests/integration/cli-integration.test.js`): 11 テストケース
  - `validate` ↔ `stats` 出力整合性検証（スキル数一致）
  - `stats` 出力とファイルシステム実数との照合
  - `pipeline list` 26 パイプライン完全性・スキル連鎖セパレータ検証
  - `validate --verbose` 全スキル ✔ マーク検証
  - `init --dry-run` ファイル数とファイルシステム実数の照合
  - `help` 全コマンド記載検証
- **GAP ANALYSIS v0.26** (`docs/GAPs/GAP_ANALYSIS_v0.26.md`): v0.25→v0.26 差分分析・7 ギャップ中 6 解消記録・TU ベースライン確立・v0.27 ロードマップ
- **TU カバレッジ 100% 達成**: 全 87 未連携スキルに ToolUniverse 連携セクション + frontmatter `tu_tools` を一括追加
  - カバレッジ: 54.2% (103/190) → **100.0% (190/190)**
  - 追加スクリプト: `scripts/add-tu-keys.js`（ドメイン別 TU key マッピング定義）
  - 使用 TU key 例: `openml` (ML), `biotools` (統計/基盤), `papers_with_code` (先端計算), `cellxgene` (シングルセル), `tcia` (医用画像) 等
  - TU カバレッジテストのアサーションを `>=50%` → `100%` に引き上げ
- **pre-commit hooks** (`husky` + `lint-staged`): コミット時に Biome 自動チェック
  - `husky` v9 + `lint-staged` v16 を devDependencies に追加
  - `.husky/pre-commit`: `npx lint-staged` 実行
  - `"prepare": "husky"` スクリプト自動追加
  - `lint-staged` 設定: `bin/`, `tests/`, `vitest.config.js` に `biome check --write` 適用
- **npm publish ワークフロー** (`.github/workflows/publish.yml`)
  - GitHub Release 作成時に自動で `npm publish --provenance --access public`
  - テスト → validate → lint → publish の安全なパイプライン
  - npm provenance (SLSA) 対応

---

## [0.25.5] - 2026-02-14

### Added
- **`satori pipeline suggest`** コマンド: キーワード入力でインタラクティブにパイプラインを推薦
- **`satori pipeline list`** コマンド: 全 26 ドメインパイプラインの一覧表示
- **パイプラインドキュメント** (`docs/SATORI_PIPELINE_EXAMPLES.md`): 50 パイプライン（26 ドメイン + 15 クロスドメイン + 5 産業特化 + 4 方法論特化）の完全ガイド
- **逆引き辞典** (`docs/SATORI_REVERSE_INDEX.md`): スキル名・TU キー・成果物パスからの逆引き検索
- **エラーハンドリング指針** (Section 5): チェックポイント戦略・フォールバック戦略・ロギング標準
- **共有スキルカスタマイズ指針** (Section 6): `publication-figures`/`statistical-testing`/`data-preprocessing` のパイプライン別パラメータ推奨

### Changed
- **`data-preprocessing`**: チェックポイント永続化を追加 — `preprocessed_data.csv` と `preprocessing_summary.json` を `results/` に自動保存
- **`metagenome-assembled-genomes`**: 構造化出力を追加 — `mag_quality_summary.csv`, `mag_taxonomy.csv`, `representative_mags.fasta`, `mag_pipeline_summary.json`

### Fixed
- Pipeline #21 ステップ詳細: `clinical-standards` の TU ツール欄を `—` → `loinc`, `icd` に修正
- Pipeline #1, #8, #18 の「期待される成果物」に新規出力ファイルを追記

### Docs
- `docs/GAPs/` にギャップ分析ドキュメントを整理・移動

---

## [0.25.0] - 2026-02-13

### Added
- **Phase 17: 8 新スキル — ML/AI 先端手法・統計シミュレーション・適応実験・放射線 AI 拡張** (182→190 スキル、TU 連携 131 維持)

#### R. 量子・先端計算（9→11 種、+2 新スキル）
- **scientific-federated-learning** スキル (#183): Flower FL パイプライン・FedAvg/FedProx 集約戦略・Opacus DP-SGD 差分プライバシー・Dirichlet 非 IID 分割
- **scientific-neural-architecture-search** スキル (#184): Optuna NAS 構造探索・Pareto 多目的最適化 (精度 vs サイズ)・探索空間定義・枝刈り

#### C. 機械学習・モデリング（9→11 種、+2 新スキル）
- **scientific-semi-supervised-learning** スキル (#185): Self-Training パイプライン・Label Propagation/Spreading・Pseudo-Labeling 品質評価・ラベル効率化
- **scientific-multi-task-learning** スキル (#186): Hard Parameter Sharing MTL・GradNorm 動的タスクバランシング・マルチ出力回帰/分類

#### B. 統計・探索的解析（9→11 種、+2 新スキル）
- **scientific-statistical-simulation** スキル (#187): Monte Carlo シミュレーション・Bootstrap BCa 信頼区間・Permutation Test・統計的検出力分析
- **scientific-streaming-analytics** スキル (#188): River オンライン学習 (HT/ARF)・ストリーミング異常検知 (Z-score/IQR/EWMA)・概念ドリフト検出 (ADWIN/DDM)

#### S. 医用イメージング（1→2 種、+1 新スキル）
- **scientific-radiology-ai** スキル (#189): MONAI CADe/CADx パイプライン・CT/MRI 分類・Grad-CAM 説明可能性・構造化放射線レポート・AI-RADS

#### D. 実験計画・プロセス最適化（2→3 種、+1 新スキル）
- **scientific-adaptive-experiments** スキル (#190): Thompson Sampling/UCB バンディット・Wald SPRT 逐次検定・ベイズ適応用量探索 (CRM)

### Changed
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリーを更新
- カテゴリ数: 26 (A-Z) — 変更なし
- カテゴリ展開: R(9→11), C(9→11), B(9→11), S(1→2), D(2→3)
- パイプライン統合: deep-learning→federated-learning→model-monitoring, automl→neural-architecture-search→deep-learning, active-learning→semi-supervised-learning→ml-classification, deep-learning→multi-task-learning→feature-importance, doe→statistical-simulation→statistical-testing, anomaly-detection→streaming-analytics→model-monitoring, medical-imaging→radiology-ai→clinical-report, doe→adaptive-experiments→statistical-testing

---

## [0.24.0] - 2025-07-27

### Added
- **Phase 16: 8 新スキル — ML/AI インフラ・データ解析・可視化スキル拡張** (174→182 スキル、TU 連携 131 維持)

#### C. 機械学習・モデリング（6→9 種、+3 新スキル）
- **scientific-anomaly-detection** スキル (#175): Isolation Forest/LOF/OCSVM アンサンブル異常検知・Autoencoder 異常検出・SPC 管理図 (Individuals-MR/CUSUM)
- **scientific-causal-ml** スキル (#176): DoWhy 因果推論フレームワーク・EconML Double ML/Causal Forest・S/T/X-Learner メタラーナー・CATE 推定
- **scientific-model-monitoring** スキル (#177): データドリフト検出 (KS/PSI/Wasserstein)・性能劣化検出 (スライディングウィンドウ)・A/B テストモデル比較 (Bootstrap CI)

#### E. 信号・スペクトル・時系列（4→5 種、+1 新スキル）
- **scientific-time-series-forecasting** スキル (#178): Prophet/NeuralProphet ML 予測・時系列特徴量エンジニアリング (ラグ/ローリング)・バックテストフレームワーク

#### B. 統計・探索的解析（6→9 種、+3 新スキル）
- **scientific-data-profiling** スキル (#179): ydata-profiling 自動 EDA・データ品質スコア (5 次元)・Great Expectations バリデーション
- **scientific-geospatial-analysis** スキル (#180): GeoPandas 地理空間処理・Moran's I/LISA 空間自己相関・Kriging 補間・Folium インタラクティブ地図
- **scientific-network-visualization** スキル (#181): NetworkX グラフ構築・Louvain/Leiden コミュニティ検出・中心性分析 (5 指標)・PyVis インタラクティブネットワーク

#### N. 科学プレゼンテーション・図式（3→4 種、+1 新スキル）
- **scientific-reproducible-reporting** スキル (#182): Quarto 科学文書テンプレート・Jupyter Book 多章構成・Papermill パラメトリック実行・nbconvert 一括変換

### Changed
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリーを更新
- カテゴリ数: 26 (A-Z) — 変更なし
- カテゴリ展開: B(6→9), C(6→9), E(4→5), N(3→4)
- パイプライン統合: eda→data-profiling→eda-correlation, eda→anomaly-detection→ml-classification, causal-inference→causal-ml→feature-importance, time-series→time-series-forecasting→model-monitoring, environmental-geodata→geospatial-analysis→advanced-visualization, eda-correlation→network-visualization→advanced-visualization, interactive-dashboard→reproducible-reporting→presentation-design, ensemble-methods→model-monitoring→anomaly-detection

---

## [0.23.0] - 2025-07-27

### Added
- **Phase 15: 8 新スキル — ML/AI・データ解析・可視化スキル拡張** (166→174 スキル、TU 連携 131 維持)

#### C. 機械学習・モデリング（3→6 種、+3 新スキル）
- **scientific-active-learning** スキル (#167): 不確実性サンプリング・QBC・バッチ AL・能動学習ループ・停止基準・戦略比較
- **scientific-automl** スキル (#168): Optuna HPO・マルチモデル AutoML・自動特徴量エンジニアリング・AutoML レポート
- **scientific-ensemble-methods** スキル (#169): XGBoost/LightGBM/CatBoost 比較・Stacking OOF・Voting アンサンブル・多様性評価

#### R. 量子・先端計算（7→9 種、+2 新スキル）
- **scientific-transfer-learning** スキル (#170): Vision/NLP ファインチューニング・Few-shot (Prototypical Net)・知識蒸留・ドメイン適応
- **scientific-uncertainty-quantification** スキル (#171): Conformal Prediction・MC Dropout・深層アンサンブル不確実性・Calibration 評価・ECE

#### B. 統計・探索的解析（4→6 種、+2 新スキル）
- **scientific-missing-data-analysis** スキル (#172): 欠損パターン診断 (MCAR/MAR/MNAR)・Little's MCAR テスト・MICE 多重代入・KNN/MissForest 補完・Rubin's Rules
- **scientific-advanced-visualization** スキル (#173): Plotly 3D・Altair 宣言的可視化・Parallel Coordinates・Radar・出版品質図 (Nature/Science style)・アニメーション

#### N. 科学プレゼンテーション・図式（2→3 種、+1 新スキル）
- **scientific-interactive-dashboard** スキル (#174): Streamlit/Dash/Panel/Voilà 科学データダッシュボード・パラメータ探索 UI・ウィジェット連動

### Changed
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリーを更新
- カテゴリ数: 26 (A-Z) — 変更なし
- カテゴリ展開: B(4→6), C(3→6), N(2→3), R(7→9)
- パイプライン統合: eda→missing-data→ml, eda→automl→ensemble→uncertainty-quantification→explainable-ai, advanced-visualization→interactive-dashboard→presentation-design

---

## [0.22.0] - 2025-07-27

### Added
- **Phase 14: 6 新スキル + 8 既存 TU key 追加** — ToolUniverse & K-Dense ギャップ分析に基づく v0.22.0 科学スキル拡張 (160→166 スキル、TU 連携 124→131)

#### Track A: 既存スキル TU key 追加 (8 件、うち 6 件が新規 TU 連携化)
- **scientific-epigenomics-chromatin** に `tu_tools: chipatlas` (ChIP-Atlas エピゲノミクスエンリッチメント解析) を追加 ★新規 TU 連携
- **scientific-metabolomics** に `tu_tools: metabolomics_workbench` (Metabolomics Workbench REST API メタボロームデータ・RefMet) を追加 (既存 hmdb に併記)
- **scientific-systems-biology** に `tu_tools: bigg_models, complex_portal, wikipathways` (ゲノムスケール代謝モデル、EBI タンパク質複合体、コミュニティパスウェイ) を追加 ★新規 TU 連携
- **scientific-immunoinformatics** に `tu_tools: imgt, sabdab, therasabdab` (国際免疫遺伝学情報システム、構造抗体データベース、治療用抗体構造データベース) を追加 (既存 iedb に併記)
- **scientific-public-health-data** に `tu_tools: nhanes, medlineplus, odphp` (全米健康栄養調査、NLM 健康情報 API、Healthy People ガイドライン) を追加 ★新規 TU 連携
- **scientific-epidemiology-public-health** に `tu_tools: who_gho` (WHO Global Health Observatory 健康統計 API) を追加 ★新規 TU 連携
- **scientific-model-organism-db** に `tu_tools: impc, mpd` (国際マウス表現型解析コンソーシアム、Mouse Phenome Database) を追加 ★新規 TU 連携
- **scientific-environmental-ecology** に `tu_tools: gbif` (地球規模生物多様性情報ファシリティ) を追加 ★新規 TU 連携

#### Track B: 新規スキル (6 件、うち 1 件 TU 連携)
- **F. 生命科学・オミクス（26→28 種）**
  - **scientific-glycomics** スキル (#161): GlyGen/GlyConnect/GlyCosmos 糖鎖データベース統合・糖タンパク質部位検索・MS 断片化予測
  - **scientific-lipidomics** スキル (#162): LipidMAPS/SwissLipids 脂質構造検索・サブクラス分類・脂質差次解析・脂質エンリッチメント
- **V. マイクロバイオーム・環境（8→9 種）**
  - **scientific-metagenome-assembled-genomes** スキル (#163): MetaBAT2/CONCOCT ビニング・CheckM2 品質評価・GTDB-Tk 分類・MAG パイプライン
- **M. 実験室自動化・データ管理（2→3 種）**
  - **scientific-crispr-design** スキル (#164): CRISPR gRNA 設計・Cas9/Cas12a PAM 検索・オフターゲットスコアリング・sgRNA ライブラリ構築
- **P. ファーマコビジランス・薬理ゲノミクス（3→4 種）**
  - **scientific-clinical-pharmacology** スキル (#165): PopPK NLME・PBPK シミュレーション・TDM 投与量最適化・Emax PD モデリング
- **H. 臨床・疫学・メタ科学（6→7 種）**
  - **scientific-clinical-standards** スキル (#166): LOINC/ICD-10/ICD-11 臨床標準コード検索・FHIR R4 マッピング・用語相互運用・loinc, icd = TU ツール連携

### Changed
- ToolUniverse SMCP 連携スキル: 124→131 (+1 新規 TU 連携: clinical-standards + 6 既存新規 TU 連携化: epigenomics-chromatin, systems-biology, public-health-data, epidemiology-public-health, model-organism-db, environmental-ecology + 2 既存 TU key 追加: metabolomics, immunoinformatics)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新
- カテゴリ数: 26 (A-Z) — 変更なし
- カテゴリ展開: F(26→28), H(6→7), M(2→3), P(3→4), V(8→9)

## [0.21.0] - 2025-07-26

### Added
- **Phase 13: 6 新スキル + 8 既存 TU key 追加** — ToolUniverse & K-Dense ギャップ分析に基づく v0.21.0 科学スキル拡張 (154→160 スキル、TU 連携 114→124)

#### Track A: 既存スキル TU key 追加 (8 件)
- **scientific-human-cell-atlas** に `tu_tools: cellxgene_census` (CELLxGENE Census 大規模シングルセルアトラス API) を追加 (既存 hca_tools に併記)
- **scientific-drug-repurposing** に `tu_tools: pharos` (IDG Pharos/TCRD ターゲットナレッジベース) を追加
- **scientific-pharmacogenomics** に `tu_tools: fda_pharmacogenomic_biomarkers` (FDA 薬理ゲノミクスバイオマーカーテーブル) を追加
- **scientific-human-protein-atlas** に `tu_tools: hpa` (組織/細胞タンパク質発現・RNA 発現・がん予後) を追加
- **scientific-variant-effect-prediction** に `tu_tools: spliceai, cadd` (スプライシング効果予測、統合アノテーション枯渇スコア) を追加
- **scientific-gtex-tissue-expression** に `tu_tools: gtex_v2` (GTEx Portal REST API v2 組織特異的発現・eQTL) を追加 — 空配列から更新
- **scientific-biothings-idmapping** に `tu_tools: biothings` (MyGene/MyVariant/MyChem 統合アノテーション API) を追加
- **scientific-protein-structure-analysis** に `tu_tools: proteinsplus` (タンパク質結合部位検出・構造解析ツール群) を追加

#### Track B: 新規スキル (6 件、うち 3 件 TU 連携)
- **T. シングルセル・空間・エピゲノミクス（12→13 種）**
  - **scientific-cellxgene-census** スキル (#155): CELLxGENE Census API 大規模シングルセルアトラスデータアクセス・細胞型分布・遺伝子発現マトリクス・cellxgene_census = TU ツール連携
- **J. 創薬・ファーマコロジー（8→9 種）**
  - **scientific-pharos-targets** スキル (#156): Pharos/TCRD IDG GraphQL API ターゲット TDL 分類・疾患関連・リガンドアクティビティ検索・pharos = TU ツール連携
- **L. 精密医療・臨床意思決定（5→6 種）**
  - **scientific-clingen-curation** スキル (#157): ClinGen API 遺伝子-疾患バリディティ・投与量感受性・臨床アクショナビリティ・clingen = TU ツール連携
- **Z. 科学テキストマイニング（2→3 種）**
  - **scientific-clinical-nlp** スキル (#158): MedSpaCy/scispaCy 臨床テキスト NER・否定文検出 (NegEx)・セクション分類・ICD-10/SNOMED-CT リンキング
- **F. 生命科学・オミクス（24→26 種）**
  - **scientific-hgnc-nomenclature** スキル (#159): HGNC REST API 遺伝子命名法・公式シンボル検索・エイリアス/旧シンボル解決・遺伝子ファミリー
  - **scientific-metabolomics-network** スキル (#160): 代謝物相関ネットワーク構築 (GGM/WGCNA)・KEGG パスウェイグラフ・ハブ代謝物同定・パスウェイエンリッチメント

### Changed
- ToolUniverse SMCP 連携スキル: 114→124 (+3 新規 TU 連携: cellxgene_census, pharos, clingen + 7 既存 TU key 追加: pharos, fda_pharmacogenomic_biomarkers, hpa, spliceai/cadd, gtex_v2, biothings, proteinsplus)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新
- カテゴリ数: 26 (A-Z) — 変更なし
- カテゴリ展開: F(24→26), J(8→9), L(5→6), T(12→13), Z(2→3)

## [0.20.0] - 2025-07-25

### Added
- **Phase 12: 6 新スキル + 12 既存 TU key 追加** — ToolUniverse & K-Dense ギャップ分析に基づく v0.20.0 科学スキル拡張 (148→154 スキル、TU 連携 99→114)

#### Track A: 既存スキル TU key 追加 (12 件)
- **scientific-cancer-genomics** に `tu_tools: cosmic, cbioportal` (がん体細胞変異カタログ、がんゲノミクスポータル) を追加
- **scientific-precision-oncology** に `tu_tools: oncokb` (精密腫瘍学アノテーション) を追加
- **scientific-immunoinformatics** に `tu_tools: iedb` (免疫エピトープデータベース) を追加
- **scientific-microbiome-metagenomics** に `tu_tools: mgnify` (EBI メタゲノミクス解析プラットフォーム) を追加
- **scientific-cell-line-resources** に `tu_tools: cellosaurus` (細胞株データベース ExPASy) を追加
- **scientific-string-network-api** に `tu_tools: ppi` (タンパク質・化学物質相互作用ネットワーク) を追加
- **scientific-metabolic-modeling** に `tu_tools: biomodels` (SBML モデルリポジトリ EBI) を追加
- **scientific-chembl-assay-mining** に `tu_tools: chembl` (創薬生理活性データベース EBI) を追加
- **scientific-pharmacology-targets** に `tu_tools: bindingdb, gtopdb, brenda` (結合親和性・GtoPdb・酵素動態) を追加
- **scientific-drug-target-profiling** に `tu_tools: dgidb` (薬物-遺伝子相互作用データベース) を追加
- **scientific-admet-pharmacokinetics** に `tu_tools: pubchem` (化合物・物質・生理活性アッセイ) を追加
- **scientific-metabolomics** に `tu_tools: hmdb` (ヒトメタボロームデータベース) を追加

#### Track B: 新規スキル (6 件、うち 3 件 TU 連携)
- **G. 化学・材料・イメージング（8→9 種）**
  - **scientific-stitch-chemical-network** スキル (#154): STITCH REST API 化学物質-タンパク質相互作用検索・信頼度スコアリング・ネットワーク薬理学・ポリファーマコロジー解析・stitch = TU ツール連携
- **H. 臨床・疫学・メタ科学（5→6 種）**
  - **scientific-biobank-cohort** スキル (#151): UK Biobank/BBJ/All of Us 大規模コホートデータ解析・フェノタイプ辞書検索・GWAS サマリー統計処理・PheWAS パイプライン
- **Q. 腫瘍学・疾患研究（8→10 種）**
  - **scientific-monarch-ontology** スキル (#149): Monarch Initiative REST API 疾患-遺伝子-表現型オントロジー・HPO フェノタイピング・エンティティ検索・monarch = TU ツール連携
  - **scientific-gdc-portal** スキル (#150): NCI Genomic Data Commons REST API プロジェクト横断検索・ケースメタデータ・体細胞変異 (SSM)・gdc = TU ツール連携
- **T. シングルセル・空間・エピゲノミクス（11→12 種）**
  - **scientific-spatial-multiomics** スキル (#152): MERFISH/CODEX 空間マルチオミクス統合・マルチモーダルアライメント・空間共検出解析・近傍グラフ空間コミュニティ検出
- **W. システム生物学（3→4 種）**
  - **scientific-metabolic-flux** スキル (#153): 13C/15N 安定同位体代謝フラックス解析・EMU モデリング・MID フィッティング・フラックス推定パイプライン

### Changed
- ToolUniverse SMCP 連携スキル: 99→114 (+3 新規 TU 連携: monarch, gdc, stitch + 12 既存 TU key 追加: cosmic, cbioportal, oncokb, iedb, mgnify, cellosaurus, ppi, biomodels, chembl, bindingdb/gtopdb/brenda, dgidb, pubchem, hmdb)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新
- カテゴリ数: 26 (A-Z) — 変更なし
- カテゴリ展開: G(8→9), H(5→6), Q(8→10), T(11→12), W(3→4)

## [0.19.0] - 2025-07-24

### Added
- **Phase 11: 8 新スキル + 6 既存 TU key 追加** — ToolUniverse & K-Dense ギャップ分析に基づく v0.19.0 科学スキル拡張 (140→148 スキル、TU 連携 85→99)

#### Track A: 既存スキル TU key 追加 (6 件)
- **scientific-rare-disease-genetics** に `tu_tools: orphanet` (希少疾患分類・遺伝子関連・有病率データベース) を追加
- **scientific-disease-research** に `tu_tools: disgenet` (疾患-遺伝子関連スコア GDA データベース) を追加
- **scientific-protein-interaction-network** に `tu_tools: intact` (分子相互作用データベース EBI) を追加
- **scientific-metabolomics-databases** に `tu_tools: metacyc` (代謝パスウェイ・反応・化合物データベース) を追加
- **scientific-compound-screening** に `tu_tools: zinc` (購入可能化合物データベース) を追加
- **scientific-variant-interpretation** に `tu_tools: clinvar` (臨床的バリアント解釈データベース) を追加

#### Track B: 新規 TU 連携スキル (8 件)
- **F. 生命科学・オミクス（22→24 種）**
  - **scientific-uniprot-proteome** スキル (#141): UniProt REST API プロテオーム検索・Swiss-Prot フィルタ・エントリ詳細・非同期 ID マッピング・ドメイン/特徴抽出・uniprot = TU ツール連携
  - **scientific-reactome-pathways** スキル (#144): Reactome Content Service REST API パスウェイ検索・UniProt→パスウェイマッピング・参加者取得・reactome = TU ツール連携
- **J. 創薬・ファーマコロジー（7→8 種）**
  - **scientific-drugbank-resources** スキル (#146): DrugBank API 薬剤テキスト検索・薬理 MOA・標的タンパク質逆引き・薬物相互作用 (DDI)・drugbank = TU ツール連携
- **K. 構造生物学・タンパク質工学（6→7 種）**
  - **scientific-rcsb-pdb-search** スキル (#142): RCSB PDB Search API + Data API テキスト/解像度/実験手法検索・エントリメタデータ・リガンド情報・rcsb_pdb, rcsb_search = TU ツール連携
- **L. 精密医療・臨床意思決定（3→5 種）**
  - **scientific-civic-evidence** スキル (#147): CIViC REST API がんバリアント臨床解釈・エビデンスアイテム・遺伝子サマリー・アサーション・civic = TU ツール連携
  - **scientific-gnomad-variants** スキル (#148): gnomAD GraphQL API 集団アレル頻度・遺伝子制約 (pLI/LOEUF)・リージョンクエリ・gnomad = TU ツール連携
- **Q. 腫瘍学・疾患研究（6→8 種）**
  - **scientific-opentargets-genetics** スキル (#143): Open Targets Platform GraphQL API 標的-疾患アソシエーション・薬剤エビデンス・L2G 遺伝的関連・opentarget = TU ツール連携
  - **scientific-depmap-dependencies** スキル (#145): DepMap Portal API / Sanger Cell Model Passports CRISPR/RNAi 遺伝子依存性・薬剤感受性・depmap = TU ツール連携

### Changed
- ToolUniverse SMCP 連携スキル: 85→99 (+8 新規: uniprot, rcsb_pdb/rcsb_search, opentarget, reactome, depmap, drugbank, civic, gnomad + 6 既存追加: orphanet, disgenet, intact, metacyc, zinc, clinvar)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新
- カテゴリ数: 26 (A-Z) — 変更なし

## [0.18.0] - 2026-02-23

### Added
- **8 新スキル (7 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense ギャップ分析に基づく Phase 10 科学スキル拡張 (132→140 スキル、TU 連携 79→85)
- **A. 基盤・ワークフロー（16→17 種）**
  - **scientific-crossref-metadata** スキル (#139): CrossRef REST API DOI 解決・論文メタデータ・引用数・ジャーナル情報・助成金レジストリ検索・crossref = TU ツール連携
- **F. 生命科学・オミクス（20→22 種）**
  - **scientific-arrayexpress-expression** スキル (#135): BioStudies/ArrayExpress REST API 発現実験検索・SDRF サンプルメタデータ・発現マトリクスダウンロード・データ再解析統合パイプライン・arrayexpress = TU ツール連携
  - **scientific-gtex-tissue-expression** スキル (#137): GTEx Portal REST API v2 組織特異的遺伝子発現 (中央値 TPM)・eQTL ルックアップ・多組織比較マトリクス・TU 外スキル (直接 API)
- **I. Deep Research・文献検索（3→4 種）**
  - **scientific-semantic-scholar** スキル (#136): Semantic Scholar Academic Graph API 論文検索・著者プロファイル (h-index)・引用/被引用グラフ・TLDR 自動要約・年次引用傾向分析・semantic_scholar = TU ツール連携
- **K. 構造生物学・タンパク質工学（5→6 種）**
  - **scientific-alphafold-structures** スキル (#134): AlphaFold Protein Structure Database REST API 構造予測取得・pLDDT 残基レベル信頼度プロファイリング・PAE マトリクスドメイン境界推定・バッチパイプライン・alphafold = TU ツール連携
- **P. ファーマコビジランス・薬理ゲノミクス（2→3 種）**
  - **scientific-pharmgkb-pgx** スキル (#138): PharmGKB REST API 臨床アノテーション・薬物遺伝子関連・投与量ガイドライン (CPIC/DPWG)・スターアレル表現型対応・pharmgkb = TU ツール連携
- **Q. 腫瘍学・疾患研究（5→6 種）**
  - **scientific-icgc-cancer-data** スキル (#140): ICGC DCC API 国際がんゲノムプロジェクト検索・ドナー/検体・体細胞変異 (SSM) 検索・がん種統計・遺伝子変異頻度・TU 外スキル (直接 API)
- **Y. 集団遺伝学（1→2 種）**
  - **scientific-gwas-catalog** スキル (#133): NHGRI-EBI GWAS Catalog REST API 関連解析 (P 値フィルタリング)・研究メタデータ検索・バリアント PheWAS (多形質解析)・統合パイプライン・gwas = TU ツール連携

### Changed
- ToolUniverse SMCP 連携スキル: 79→85 (+6 スキル: gwas-catalog [gwas]、alphafold-structures [alphafold]、arrayexpress-expression [arrayexpress]、semantic-scholar [semantic_scholar]、pharmgkb-pgx [pharmgkb]、crossref-metadata [crossref])
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新

## [0.17.0] - 2026-02-22

### Added
- **8 新スキル (4 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense ギャップ分析に基づく Phase 9 科学スキル拡張 (124→132 スキル、TU 連携 74→79)
- **F. 生命科学・オミクス（18→20 種）**
  - **scientific-geo-expression** スキル (#127): GEO REST API (E-utilities) データセット検索・GSE 発現マトリクス取得 (GEOparse)・差次的発現解析 (t-test/Wilcoxon/BH-FDR)・発現プロファイリング統合パイプライン・geo (3) = 3 TU ツール連携
  - **scientific-parasite-genomics** スキル (#132): VEuPathDB (PlasmoDB/VectorBase/ToxoDB/TriTrypDB) 遺伝子検索・機能アノテーション (GO/InterPro)・薬剤標的候補スクリーニング (キナーゼ/プロテアーゼ/トランスポーター)・TU 外スキル (直接 API)
- **T. シングルセル・空間・エピゲノミクス（8→11 種）**
  - **scientific-encode-screen** スキル (#125): ENCODE REST API 実験/ファイル検索・SCREEN cCRE 候補シス制御エレメント検索・ChIP-Atlas エンリッチメント解析・制御ゲノミクス統合パイプライン・encode (5)/chipatlas (1) = 6 TU ツール連携
  - **scientific-human-cell-atlas** スキル (#126): HCA Data Portal プロジェクト/ファイル検索・CELLxGENE Census 大規模アトラスクエリ・細胞タイプ組成分析・統合アトラスパイプライン・hca_tools (1) = 1 TU ツール連携
  - **scientific-squidpy-advanced** スキル (#131): Squidpy 空間自己相関 (Moran's I/Geary's C)・空間共起解析・近傍エンリッチメント・空間ニッチ同定 (KMeans)・中心性スコア・K-Dense: squidpy-advanced
- **V. マイクロバイオーム・環境（6→8 種）**
  - **scientific-environmental-geodata** スキル (#128): SoilGrids REST API (ISRIC) 土壌特性取得・WorldClim バイオクリマティック変数 (BIO1-19)・SDM 環境変数スタック・TU 外スキル (直接 API)
  - **scientific-paleobiology** スキル (#129): Paleobiology Database (PBDB) REST API 化石産出記録/分類群検索・地質年代多様性曲線・古地理サマリ・paleobiology (1) = 1 TU ツール連携
- **W. システム生物学（2→3 種）**
  - **scientific-metabolic-atlas** スキル (#130): Metabolic Atlas/Human-GEM REST API 代謝反応/代謝産物検索・代謝ネットワーク (NetworkX DiGraph) 構築・ハブ代謝産物同定・K-Dense: metabolic-atlas

### Changed
- ToolUniverse SMCP 連携スキル: 74→79 (+5 スキル: geo-expression [geo]、encode-screen [encode, chipatlas]、human-cell-atlas [hca_tools]、paleobiology [paleobiology])
- K-Dense 新規カバー: squidpy-advanced, metabolic-atlas (+2)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新

## [0.16.0] - 2026-02-21

### Added
- **8 新スキル (5 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense ギャップ分析に基づく Phase 8 科学スキル拡張 (116→124 スキル、TU 連携 70→74)
- **A. 基盤・ワークフロー（15→16 種）**
  - **scientific-data-submission** スキル (#119): GenBank/SRA/GEO/BioProject/BioSample データ投稿・Aspera アップロード・GEO SOFT フォーマット・FAIR 原則準拠チェックリスト (F/A/I/R スコアリング)・TU 外スキル (直接 API)
- **J. 創薬・ファーマコロジー（6→7 種）**
  - **scientific-nci60-screening** スキル (#120): CellMiner NCI-60 薬剤活性 (GI50/TGI/LC50)・組織別応答パターン・薬剤-分子マーカー相関 (Pearson/Bonferroni)・DepMap CRISPR/RNAi 依存性統合・TU 外スキル (直接 API)
- **T. シングルセル・空間・エピゲノミクス（6→8 種）**
  - **scientific-scatac-signac** スキル (#123): SnapATAC2 scATAC-seq 前処理・ピークコーリング (MACS2/Tile)・chromVAR モチーフエンリッチメント・Gene Activity スコア・RNA+ATAC WNN マルチモーダル統合 (muon)・K-Dense: signac
  - **scientific-gpu-singlecell** スキル (#124): rapids-singlecell GPU 前処理・cuML GPU PCA/kNN・cuGraph Leiden/Louvain クラスタリング・GPU UMAP/t-SNE・CPU vs GPU ベンチマーク・大規模 (>1M cells) パイプライン・K-Dense: rapids-singlecell
- **V. マイクロバイオーム・環境（3→6 種）**
  - **scientific-rrna-taxonomy** スキル (#118): SILVA SSU/LSU rRNA リファレンスダウンロード (v138.1)・Greengenes2 分類体系・MGnify メタゲノム検索/分類・QIIME2 sklearn 分類器・マルチ分類器コンセンサス・mgnify (1) = 1 TU ツール連携
  - **scientific-plant-biology** スキル (#121): Plant Reactome 代謝パスウェイ検索・TAIR Arabidopsis 遺伝子/発現情報・Ensembl Plants 種間オーソログ・Phytozome 比較ゲノミクス・TU 外スキル (直接 API)
  - **scientific-marine-ecology** スキル (#122): OBIS 海洋生物出現記録/チェックリスト・WoRMS 分類学検索/階層・GBIF 種検索/出現記録・FishBase 魚類データ・obis (1)/worms (1)/gbif (1) = 3 TU ツール連携
- **X. 疫学・公衆衛生（2→3 種）**
  - **scientific-toxicology-env** スキル (#117): CTD 化学物質-遺伝子/疾患相互作用・CompTox Dashboard Tox21/ToxCast アッセイ・T3DB 毒素検索・EPA IRIS リスク評価・毒性パスウェイ統合解析・TU 外スキル (直接 API)

### Changed
- ToolUniverse SMCP 連携スキル: 70→74 (+4 スキル: rrna-taxonomy [mgnify]、marine-ecology [obis, worms, gbif])
- K-Dense 新規カバー: signac, rapids-singlecell (+2)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新

## [0.15.0] - 2026-02-20

### Added
- **10 新スキル (3 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense ギャップ分析に基づく Phase 7 科学スキル拡張 (106→116 スキル、TU 連携 66→70)
- **F. 生命科学・オミクス（14→18 種）**
  - **scientific-ensembl-genomics** スキル (#108): Ensembl REST API 遺伝子検索・配列取得・VEP バリアント効果予測 (SIFT/PolyPhen/CADD)・クロスリファレンス ID 変換・ホモロジー検索 (オーソログ/パラログ)・制御領域・遺伝子ツリー・ensembl (16) = 16 TU ツール連携
  - **scientific-string-network-api** スキル (#109): STRING v12 PPI ネットワーク・BioGRID 実験的 PPI・STITCH 化学-タンパク質相互作用・NetworkX トポロジー解析・コミュニティ検出・機能濃縮・ppi (2)/stitch (3) = 5 TU ツール連携
  - **scientific-expression-comparison** スキル (#110): EBI Expression Atlas ベースライン発現・差次的発現検索・実験メタデータ・組織横断比較マトリクス・ヒートマップ・expression_atlas (4) = 4 TU ツール連携
  - **scientific-model-organism-db** スキル (#111): FlyBase/WormBase/ZFIN/RGD/MGI 直接 REST API・種間オーソログ検索・表現型データ・TU 外スキル (直接 API)
- **G. 化学・材料・イメージング（4→8 種）**
  - **scientific-chembl-assay-mining** スキル (#107): ChEMBL REST API ターゲット検索・バイオアクティビティデータ・アッセイ検索・SAR 解析・選択性プロファイリング・類似性検索・ATC 分類・構造アラート・ChEMBL (19) = 19 TU ツール連携
  - **scientific-md-simulation** スキル (#112): MDAnalysis トラジェクトリ解析・RMSD/RMSF/回転半径/水素結合/SASA・OpenFF 力場パラメータ化・K-Dense: mdanalysis, openff
  - **scientific-advanced-imaging** スキル (#114): Cellpose 深層学習セルセグメンテーション・CellProfiler 形態学的プロファイリング・Cell Painting 5ch 解析・napari 3D 可視化・K-Dense: cellprofiler, napari, cellpose
  - **scientific-deep-chemistry** スキル (#115): DeepChem GCN/MPNN/AttentiveFP 分子特性予測・MoleculeNet ベンチマーク・ChemBERTa 分子表現学習・K-Dense: deepchem
- **T. シングルセル・空間・エピゲノミクス（4→6 種）**
  - **scientific-perturbation-analysis** スキル (#113): pertpy 摂動解析・Augur 摂動応答性スコアリング・scGen 摂動予測・scIB 統合ベンチマーク・摂動シグネチャ抽出・K-Dense: pertpy, scib
  - **scientific-scvi-integration** スキル (#116): scVI 変分オートエンコーダ統合・scANVI 半教師有りアノテーション・totalVI CITE-seq 結合解析・SOLO ダブレット検出・確率的差次的発現・K-Dense: scvi-tools

### Changed
- ToolUniverse SMCP 連携スキル: 66→70 (+4 スキル、~44 新規 TU ツール統合)
- K-Dense 新規カバー: mdanalysis, openff, cellprofiler, napari, cellpose, deepchem, pertpy, scib, scvi-tools (+9)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新

## [0.14.0] - 2026-02-19

### Added
- **10 新スキル (9 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense ギャップ分析に基づく Phase 6 ドメイン拡張 (96→106 スキル、TU 連携 59→66)
- **パイプライン統合アーキテクチャ**を README.md に新設 — データアクセス・パイプライン統合図でスキル間連携フローを明示
- **B. 統計・探索的解析（3→4 種）**
  - **scientific-symbolic-mathematics** スキル (#105): SymPy ODE ソルバー・記号微積分・線形代数・PK 区画モデル・LaTeX エクスポート・K-Dense: sympy
- **F. 生命科学・オミクス（12→14 種）**
  - **scientific-ontology-enrichment** スキル (#99): EFO 形質オントロジー・OLS API オントロジー横断検索・Enrichr 遺伝子セット解析・UMLS メタシソーラスクロスウォーク・efo (2)/ols (2)/enrichr (1)/umls (1) = 6 TU ツール連携
  - **scientific-ebi-databases** スキル (#100): EBI Search・ENA Browser・BioStudies・dbfetch・MetaboLights・ebi_search (2)/ena_browser (1)/biostudies (1)/dbfetch (1)/metabolights (1) = 6 TU ツール連携・K-Dense: ena-database
- **I. Deep Research・文献検索（2→3 種）**
  - **scientific-preprint-archive** スキル (#97): bioRxiv/medRxiv/arXiv/PMC/DOAJ/Unpaywall/HAL/CORE/Zenodo/OpenAIRE/OSF/Fatcat/DBLP 横断プレプリント検索・biorxiv (2)/medrxiv (1)/arxiv (1)/pmc (2)/doaj (1)/unpaywall (1)/hal (1)/core (1)/zenodo (1)/openaire (1)/osf_preprints (1)/fatcat (1)/dblp (1) = 16 TU ツール連携
- **Q. 腫瘍学・疾患研究（4→5 種）**
  - **scientific-cell-line-resources** スキル (#101): Cellosaurus 細胞株同定・STR プロファイル検証・HeLa 汚染チェック・cellosaurus (3) = 3 TU ツール連携
- **R. 量子・先端計算（6→7 種）**
  - **scientific-reinforcement-learning** スキル (#104): Stable-Baselines3 RL エージェント訓練・Gymnasium カスタム環境・PufferLib マルチエージェント・分子設計 RL・実験最適化ベイズ RL・K-Dense: stable-baselines3, pufferlib
- **T. シングルセル・空間・エピゲノミクス（3→4 種）**
  - **scientific-regulatory-genomics** スキル (#102): RegulomeDB バリアントスコアリング・ReMap 転写因子 ChIP-seq アトラス・4DN Portal 三次元ゲノム構造・regulomedb (1)/remap (1)/fourdn_portal (1) = 3 TU ツール連携
- **V. マイクロバイオーム・環境（2→3 種）**
  - **scientific-phylogenetics** スキル (#103): ETE3 系統樹構築・scikit-bio Faith's PD/UniFrac・分岐年代推定・祖先配列復元・K-Dense: etetoolkit, scikit-bio
- **X. 疫学・公衆衛生（1→2 種）**
  - **scientific-public-health-data** スキル (#98): NHANES XPT データアクセス・MedlinePlus 健康情報・RxNorm 医薬品正規化・ODPHP ガイドライン・Health Disparities API・nhanes (2)/medlineplus (2)/rxnorm (2)/health_disparities (2)/odphp (2)/guidelines_tools (2) = 12 TU ツール連携
- **Z. 科学テキストマイニング（1→2 種）**
  - **scientific-biomedical-pubtator** スキル (#106): PubTator3 バイオメディカル NER (Gene/Disease/Chemical/Species/Variant)・エンティティ関係抽出・知識グラフ構築・pubtator (2) = 2 TU ツール連携

### Changed
- ToolUniverse SMCP 連携スキル: 59→66 (+7 スキル、~48 新規 TU ツール統合、~33 新規 TU config keys)
- K-Dense 新規カバー: ena-database, etetoolkit, stable-baselines3, pufferlib, sympy (+5、推定 ~123/140)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図・パイプライン統合図を更新

## [0.13.0] - 2026-02-18

### Added
- **10 新スキル (7 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense-AI/claude-scientific-skills ギャップ分析に基づく Phase 5 ドメイン拡張 (86→96 スキル、TU 連携 50→59)
- **A. 基盤・ワークフロー（14→15 種）**
  - **scientific-biothings-idmapping** スキル (#91): MyGene.info 遺伝子アノテーション・MyVariant.info 変異アノテーション・MyChem.info 化合物アノテーション・クロス DB ID マッピング (Entrez↔Ensembl↔UniProt↔RefSeq)・バッチ統合アノテーション・biothings (7) = 7 TU ツール連携
- **F. 生命科学・オミクス（9→12 種）**
  - **scientific-human-protein-atlas** スキル (#88): HPA 組織/細胞タンパク質発現・RNA 発現プロファイル (HPA/GTEx/FANTOM5)・がん予後バイオマーカー・細胞内局在・タンパク質相互作用・hpa (14) = 14 TU ツール連携
  - **scientific-genome-sequence-tools** スキル (#90): dbSNP rsID 変異アレル頻度・BLAST 相同性検索 (blastn/blastp)・NCBI Nucleotide 配列フェッチ・GDC がんゲノミクス体細胞変異/CNV/発現・dbsnp (3)/blast (2)/ncbi_nucleotide (3)/gdc (7) = 15 TU ツール連携
  - **scientific-noncoding-rna** スキル (#92): Rfam RNA ファミリー検索/Infernal cmscan・RNAcentral ncRNA 統合 DB・共分散モデル取得・構造マッピング・系統樹・rfam (7)/rnacentral (2) = 9 TU ツール連携
- **J. 創薬・ファーマコロジー（4→6 種）**
  - **scientific-pharmacology-targets** スキル (#89): BindingDB 結合親和性・GPCRdb GPCR プロファイリング・GtoPdb 薬理学データ・BRENDA 酵素阻害剤・Pharos/TCRD 未解明ターゲット TDL 分類・bindingdb (4)/gpcrdb (5)/gtopdb (8)/brenda (1)/pharos (4) = 22 TU ツール連携
  - **scientific-compound-screening** スキル (#94): ZINC 購入可能化合物ライブラリ検索・SMILES 類似性検索・Lipinski フィルタ・VS 前処理パイプライン・zinc (4) = 4 TU ツール連携
- **K. 構造生物学・タンパク質工学（4→5 種）**
  - **scientific-structural-proteomics** スキル (#93): EMDB クライオ EM 構造・PDBe 構造品質/実験/二次構造・Proteins API ドメイン/変異/ゲノムマッピング・Complex Portal 複合体・DeepGO 機能予測・EVE 変異効果スコア・emdb (7)/pdbe_api (10)/proteins_api (10)/complex_portal (2)/deepgo (1)/eve (2) = 32 TU ツール連携
- **Q. 腫瘍学・疾患研究（3→4 種）**
  - **scientific-rare-disease-genetics** スキル (#87): OMIM 遺伝子-疾患マッピング・Orphanet 希少疾患分類・DisGeNET GDA スコア・IMPC マウス表現型・統合希少疾患遺伝子プロファイリング・omim (4)/orphanet (5)/disgenet (5)/impc (4) = 18 TU ツール連携
- **R. 量子・先端計算（5→6 種）**
  - **scientific-healthcare-ai** スキル (#96): PyHealth 臨床 ML パイプライン (MIMIC-III)・Transformer/RETAIN/GRU モデル・FlowIO フローサイトメトリー FCS 解析・医療コードマッピング (ICD-10/SNOMED/ATC)・K-Dense 参照 (pyhealth/flowio)
- **W. システム生物学（1→2 種）**
  - **scientific-metabolic-modeling** スキル (#95): BiGG Models ゲノムスケール代謝モデル検索/反応/代謝物・BioModels SBML リポジトリ検索/取得・統合代謝モデル探索・bigg_models (7)/biomodels_tools (5) = 12 TU ツール連携

### Changed
- ToolUniverse SMCP 連携スキル: 50→59 (+9 スキル、133+ 新規 TU ツール統合)
- README.md: 全カウント・カテゴリ表・スキル一覧・ディレクトリツリー・TU 連携図を更新

## [0.12.0] - 2026-02-17

### Added
- **10 新スキル (8 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense-AI/claude-scientific-skills ギャップ分析に基づく Phase 4 ドメイン拡張 (76→86 スキル)
- **A. 基盤・ワークフロー（13→14 種）**
  - **scientific-systematic-review** スキル (#84): PRISMA 2020 系統的レビュー・マルチ DB 検索戦略・タイトル/抄録スクリーニング・RoB 2/ROBINS-I/NOS バイアスリスク評価・PRISMA フロー図自動生成
- **F. 生命科学・オミクス（7→9 種）**
  - **scientific-pathway-enrichment** スキル (#77): ORA (超幾何検定)・GSEA (running sum)・KEGG/Reactome/GO/WikiPathways 統合パスウェイ濃縮解析・統合ヒートマップ・kegg (5)/reactome (20)/go (5)/wikipathways (2)/pathway_commons (2) = 34 TU ツール連携
  - **scientific-metabolomics-databases** スキル (#82): HMDB XML 検索・MetaCyc パスウェイ・Metabolomics Workbench REST API・m/z ベース多 DB 同定・RefMet 標準化・hmdb (3)/metacyc (4)/metabolomics_workbench (6) = 13 TU ツール連携
- **H. 臨床・疫学・メタ科学（4→5 種）**
  - **scientific-clinical-reporting** スキル (#85): SOAP ノート生成・バイオマーカープロファイルレポート・CPIC 準拠ファーマコゲノミクスレポート・HL7 FHIR DiagnosticReport・HTML/LaTeX/FHIR JSON 出力
- **I. Deep Research・文献検索（1→2 種）**
  - **scientific-literature-search** スキル (#78): PubMed E-utilities・Semantic Scholar API・OpenAlex・EuropePMC・CrossRef メタデータ・引用ネットワーク構築・pubmed (6)/europepmc (7)/semantic_scholar (2)/openalex (8)/crossref (6) = 29 TU ツール連携
- **J. 創薬・ファーマコロジー（3→4 種）**
  - **scientific-molecular-docking** スキル (#83): AutoDock Vina (Python API + CLI)・DiffDock AI 推論・受容体/リガンド前処理・バーチャルスクリーニングパイプライン
- **K. 構造生物学・タンパク質工学（2→4 種）**
  - **scientific-protein-interaction-network** スキル (#79): STRING v12 API・IntAct REST・STITCH 化学-タンパク質ネットワーク・PPI トポロジー解析 (中心性/コミュニティ検出)・ネットワーク可視化・intact (8)/ppi (2)/stitch (3)/humanbase (1) = 14 TU ツール連携
  - **scientific-protein-domain-family** スキル (#86): InterPro REST API ドメイン検索・InterProScan シーケンスベース予測・ドメインアーキテクチャ可視化・ファミリー比較・interpro (3)/interproscan (3) = 6 TU ツール連携
- **L. 精密医療・臨床意思決定（2→3 種）**
  - **scientific-variant-effect-prediction** スキル (#80): AlphaMissense 病原性スコア (>0.564 pathogenic)・CADD PHRED スコア (≥20 deleterious)・SpliceAI デルタスコア 4 チャネル・コンセンサス病原性判定・alphamissense (3)/cadd (3)/spliceai (3) = 9 TU ツール連携
- **Q. 腫瘍学・疾患研究（2→3 種）**
  - **scientific-cancer-genomics** スキル (#81): COSMIC CGC 変異検索・cBioPortal REST API がん研究照会・DepMap 遺伝子依存性 (Chronos スコア)・変異シグネチャー解析 (NMF)・cosmic (2)/cbioportal (6)/depmap (5) = 13 TU ツール連携

### Changed
- **README.md**: スキル数を 76→86 に更新、カテゴリ数 26 (変更なし, 既存カテゴリ拡張)、ToolUniverse 連携スキル数を 42→50 に更新、各カテゴリのスキル数・説明・ディレクトリ構造・ファイル I/O テーブルを更新
- **docs/GAP_ANALYSIS_v0.12.md**: v0.12.0 ギャップ分析ドキュメント新規追加

## [0.11.0] - 2026-02-16

### Added
- **10 新スキル (既存 9 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense-AI/claude-scientific-skills ギャップ分析に基づく Phase 3 ドメイン拡張 (66→76 スキル)
- **E. 信号・スペクトル・時系列（3→4 種）**
  - **scientific-neuroscience-electrophysiology** スキル (#67): SpikeInterface/Kilosort4 スパイクソート・Allen Institute 品質基準・MNE-Python EEG (前処理/ERP/マイクロステート)・NeuroKit2 ECG HRV/EDA 解析・脳機能結合 (wPLI/グラフ理論)
- **F. 生命科学・オミクス（5→7 種）**
  - **scientific-proteomics-mass-spectrometry** スキル (#68): pyOpenMS LC-MS/MS 前処理・ペプチド ID (FDR)・タンパク質定量 (LFQ/TMT/SILAC/iBAQ)・PTM マッピング・matchms スペクトル類似度・GNPS 分子ネットワーク・PRIDE/UniProt/KEGG/Reactome ツール連携
  - **scientific-gene-expression-transcriptomics** スキル (#69): GEO データセット取得 (GEOparse)・PyDESeq2 差次発現解析・GTEx 組織発現/eQTL 照会・GSEA/ORA 遺伝子セット濃縮解析・geo/GTEx/ExpressionAtlas/ArrayExpress (24 ツール) 連携
- **G. 化学・材料・イメージング（3→4 種）**
  - **scientific-computational-materials** スキル (#70): pymatgen 結晶構造操作・対称性解析・Materials Project API 照会・相図凸包解析・電子バンド構造/DOS 可視化・VASP/QE 入出力
- **H. 臨床・疫学・メタ科学（3→4 種）**
  - **scientific-clinical-trials-analytics** スキル (#71): ClinicalTrials.gov API v2 多基準検索・試験詳細取得・競合ランドスケープ解析・AE/アウトカム抽出・バルクエクスポート・ClinicalTrials (7)/FDA (2) ツール連携
- **M. 実験室自動化・データ管理（1→2 種）**
  - **scientific-lab-data-management** スキル (#72): Benchling ELN/DNA 設計/レジストリ・DNAnexus ゲノミクス PaaS・OMERO バイオイメージング管理・Protocols.io プロトコル共有
- **N. 科学プレゼンテーション・図式（1→2 種）**
  - **scientific-scientific-schematics** スキル (#73): CONSORT フロー図・NN アーキテクチャ図・分子パスウェイ図 (Mermaid)・TikZ 出版品質ベクター図
- **O. 研究計画・グラント・規制（2→3 種）**
  - **scientific-regulatory-science** スキル (#74): FDA Orange Book 承認履歴/特許/排他性・510(k) 医療機器クリアランス・ISO 13485 設計管理/CAPA・USPTO 特許検索・FDA_OrangeBook (6)/FAERS ツール連携
- **P. ファーマコビジランス・薬理ゲノミクス（1→2 種）**
  - **scientific-pharmacogenomics** スキル (#75): PharmGKB 遺伝子-薬物相互作用・CPIC ガイドライン・Star アレルアノテーション・代謝型判定 (PM/IM/NM/RM/UM)・FDA PGx バイオマーカー・PharmGKB (7)/FDA (3)/OpenTargets (1) ツール連携
- **T. シングルセル・空間・エピゲノミクス（2→3 種）**
  - **scientific-epigenomics-chromatin** スキル (#76): ChIP-seq MACS2/MACS3 ピークコール・ATAC-seq NFR 検出・WGBS/RRBS DMR 検出・ChromHMM 15 状態モデル・Hi-C TAD/A-B コンパートメント・DiffBind 差次的結合・ChIPAtlas (4)/JASPAR (3)/SCREEN (1)/ENCODE (3) ツール連携

### Changed
- **README.md**: スキル数を 66→76 に更新、カテゴリ数 26 (変更なし, 既存カテゴリ拡張)、ToolUniverse 連携スキル数を 32→42 に更新、各カテゴリのスキル数・説明・ディレクトリ構造・ファイル I/O テーブルを更新
- **docs/GAP_ANALYSIS_v0.11.md**: v0.11.0 ギャップ分析ドキュメント新規追加

## [0.10.0] - 2026-02-15

### Added
- **7 新カテゴリ (T-Z)・10 新スキル** を追加 — ToolUniverse ギャップ分析に基づくドメイン拡張
- **T. シングルセル・空間オミクス（2 種）**
  - **scientific-single-cell-genomics** スキル (#57): scRNA-seq QC (Scanpy)・Leiden クラスタリング・DEG (Wilcoxon)・RNA velocity (scVelo)・CellChat 細胞間通信・CELLxGENE/HCA ツール連携
  - **scientific-spatial-transcriptomics** スキル (#58): Visium/MERFISH 前処理 (Squidpy)・SVG 検出 (Moran's I)・空間ドメイン検出・cell2location デコンボリューション・空間 L-R 解析
- **U. 免疫・感染症（2 種）**
  - **scientific-immunoinformatics** スキル (#59): MHC-I/II 結合予測 (MHCflurry)・B 細胞エピトープ・TCR/BCR レパトア多様性・抗体 CDR 解析 (ANARCI)・ワクチン候補ランキング・IEDB/IMGT/SAbDab ツール連携
  - **scientific-infectious-disease** スキル (#60): 病原体 WGS QC・AMR 遺伝子検出 (ResFinder/RGI)・MLST/cgMLST 型別・系統解析 (IQ-TREE)・SIR/SEIR 数理モデル・CDC/EUHealthInfo ツール連携
- **V. マイクロバイオーム・環境（2 種）**
  - **scientific-microbiome-metagenomics** スキル (#61): DADA2 ASV パイプライン・MetaPhlAn/Kraken2 分類プロファイリング・α/β 多様性 (Shannon/Bray-Curtis/PERMANOVA)・ANCOM-BC 差次的存在量・HUMAnN 機能プロファイリング・MGnify/KEGG ツール連携
  - **scientific-environmental-ecology** スキル (#62): SDM (MaxEnt/RF/GBM)・生物多様性指数 (α/β/γ)・群集序列化 (NMDS/CCA/RDA)・保全優先順位ランキング・OBIS/GBIF ツール連携
- **W. システム生物学（1 種）**
  - **scientific-systems-biology** スキル (#63): SBML/RoadRunner 動的シミュレーション・定常状態/ヤコビアン解析・FBA/pFBA (cobrapy)・GRN 推定 (GENIE3)・Sobol 感度解析 (SALib)・BioModels/Reactome/BiGG ツール連携
- **X. 疫学・公衆衛生（1 種）**
  - **scientific-epidemiology-public-health** スキル (#64): RR/OR/RD/NNT/AF リスク指標・直接/間接年齢標準化 (SMR)・LISA/Getis-Ord 空間クラスタリング・DAG バックドア基準交絡分析・WHO/CDC/EUHealthInfo ツール連携
- **Y. 集団遺伝学（1 種）**
  - **scientific-population-genetics** スキル (#65): PLINK2 遺伝子型 QC・HWE 検定・PCA/ADMIXTURE 集団構造推定・Weir-Cockerham Fst・iHS/Tajima's D 選択スキャン・gnomAD/GWAS Catalog ツール連携
- **Z. 科学テキストマイニング（1 種）**
  - **scientific-text-mining-nlp** スキル (#66): BioBERT/SciSpaCy NER・関係抽出 (PPI/DDI/GDA)・知識グラフ構築 (Louvain)・BERTopic トピックモデリング・引用ネットワーク分析 (PageRank/HITS)・PubTator/PubMed/EuropePMC ツール連携

### Changed
- **README.md**: スキル数を 56→66 に更新、カテゴリ数を 19→26 に拡張 (T-Z 追加)、ToolUniverse 連携スキル数を 22→32 に更新、次世代オミクス・疫学パイプラインフロー図を追加、ディレクトリ構造に 7 新カテゴリを追加

## [0.9.0] - 2026-02-14

### Added
- **ToolUniverse MCP ツール連携**: 22 スキル（HIGH 13 + MEDIUM 9）の SKILL.md に `### 利用可能ツール` セクションを追加 — [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で 1,200 以上の外部科学データベースツール（ClinVar, gnomAD, OncoKB, CIViC, FAERS, UniProt, ChEMBL, PubMed 等）への参照を記載
  - **HIGH（13 スキル）**: variant-interpretation, precision-oncology, disease-research, drug-target-profiling, drug-repurposing, pharmacovigilance, clinical-decision-support, protein-structure-analysis, admet-pharmacokinetics, protein-design, bioinformatics, multi-omics, metabolomics
  - **MEDIUM（9 スキル）**: deep-research, cheminformatics, sequence-analysis, citation-checker, meta-analysis, network-analysis, graph-neural-networks, survival-clinical, grant-writing

### Changed
- **README.md**: Overview に ToolUniverse MCP 連携の説明を追加、ファイル入出力セクションに「ToolUniverse MCP ツール連携」サブセクションを追加（アーキテクチャ図付き）

## [0.8.0] - 2026-02-13

### Added
- **4 新カテゴリ (P-S)・9 新スキル** を追加 — ToolUniverse (mims-harvard) と claude-scientific-skills (K-Dense-AI) のギャップ分析 Tier-1 に基づく MECE 設計
- **P. ファーマコビジランス（1 種）**
  - **scientific-pharmacovigilance** スキル (#48): FAERS 不均衡分析 (PRR/ROR/IC/EBGM)・MedDRA 階層構造 (LLT→PT→HLT→HLGT→SOC)・時系列トレンド・人口統計層別化・Naranjo 因果評価スケール・安全性シグナルレポート生成
- **Q. 腫瘍学・疾患研究（2 種）**
  - **scientific-precision-oncology** スキル (#49): CIViC GraphQL API・OncoKB Annotation API・cBioPortal TCGA 変異頻度・TMB/MSI 判定・AMP/ASCO/CAP Tiering (I-IV)・分子腫瘍ボード (MTB) レポート生成
  - **scientific-disease-research** スキル (#50): GWAS Catalog (EBI) 検索・DisGeNET GDA スコア・Orphanet 希少疾患カタログ・HPO 表現型マッチング・Polygenic Risk Score (PRS) 算出
- **R. 量子・先端計算（5 種）**
  - **scientific-quantum-computing** スキル (#51): Qiskit VQE (H₂ 基底エネルギー)・PennyLane 量子 ML (VQC)・Cirq ノイズモデリング・QAOA MaxCut・QuTiP Jaynes-Cummings ダイナミクス
  - **scientific-graph-neural-networks** スキル (#52): PyG GCN/GAT/GIN 分子特性予測・TorchDrug MoleculeNet ベンチマーク・知識グラフ推論 (TransE/RotatE)・Scaffold Split・GNNExplainer
  - **scientific-bayesian-statistics** スキル (#53): PyMC/Stan 階層ベイズ回帰・NUTS MCMC サンプリング・事後予測チェック (PPC)・WAIC/LOO-CV モデル比較・Gaussian Process ベイズ最適化
  - **scientific-explainable-ai** スキル (#54): TreeSHAP/KernelSHAP/DeepSHAP 特徴量寄与分解・LIME ローカル説明・Captum Integrated Gradients・反実仮想説明・公平性監査
  - **scientific-deep-learning** スキル (#55): PyTorch Lightning 構造化パイプライン・timm 事前学習 CNN/ViT・HuggingFace Transformer Fine-tune・Optuna HPO (Hyperband Pruning)・ONNX/TorchScript エクスポート
- **S. 医用イメージング（1 種）**
  - **scientific-medical-imaging** スキル (#56): DICOM/NIfTI 読み込み・HU ウィンドウ処理・Isotropic リサンプリング・MONAI U-Net/SwinUNETR 3D セグメンテーション・WSI パッチ抽出 (openslide)・PyRadiomics 特徴量抽出

### Changed
- **README.md**: スキル数を 47→56 に更新、カテゴリ数を 15→19 に拡張 (P-S 追加)、先端計算・医用パイプラインフロー図を追加、ディレクトリ構造に 4 新カテゴリを追加

## [0.7.1] - 2026-02-13

### Fixed
- J-O カテゴリ 11 スキル全てに `## References` セクション（Output Files + 参照スキル テーブル）を追加 — パイプライン連携の欠落を修正

## [0.7.0] - 2026-02-13

### Added
- **6 新カテゴリ (J-O)・11 新スキル** を追加 — ToolUniverse (mims-harvard) と claude-scientific-skills (K-Dense-AI) のギャップ分析に基づく MECE 設計
- **J. 創薬・ファーマコロジー（3 種）**
  - **scientific-drug-target-profiling** スキル (#37): 9-path 標的プロファイリング (Identity→Safety)・TDL 分類 (Tclin/Tchem/Tbio/Tdark)・ドラッガビリティマトリクス (SM/Ab/PROTAC/ASO/遺伝子治療)・競合ランドスケープ
  - **scientific-admet-pharmacokinetics** スキル (#38): 5 段階 ADMET パイプライン (吸収→分布→代謝→排泄→毒性)・Lipinski/Veber/Ghose/QED ドラッグライクネス・PAINS/Brenk 構造アラート・CYP 阻害/基質予測・1-コンパートメント PK モデル
  - **scientific-drug-repurposing** スキル (#39): 7 戦略リポジショニング (Target/Compound/Disease/Mechanism/Network/Phenotype/Structure)・ネットワーク近接解析 (Guney et al.)・6 因子重み付き多基準スコアリング
- **K. 構造生物学・タンパク質工学（2 種）**
  - **scientific-protein-structure-analysis** スキル (#40): PDB Search API v2・AlphaFold DB 取得・品質評価 (Resolution/R-factor/pLDDT)・結合サイト検出 (fpocket/PDBe)
  - **scientific-protein-design** スキル (#41): ESM-2 変異スキャン (ΔLL)・de novo 設計パイプライン (RFdiffusion→ProteinMPNN→ESMFold)・バインダー/スキャフォールド/酵素設計・検証基準 (pLDDT>70, pTM>0.5)
- **L. 精密医療・臨床意思決定（2 種）**
  - **scientific-variant-interpretation** スキル (#42): ACMG/AMP 28 基準分類アルゴリズム・薬理ゲノミクス (CYP2D6/2C19/2C9/DPYD/TPMT/HLA-B, CPIC ガイドライン)・体細胞変異 OncoKB レベル・in silico 予測 (CADD/REVEL/AlphaMissense/SpliceAI)
  - **scientific-clinical-decision-support** スキル (#43): GRADE エビデンス枠組 (1a-5)・推奨グレード (A-D)・精密腫瘍学ワークフロー (OncoKB)・ClinicalTrials.gov API マッチング・リスク-ベネフィット分析
- **M. 実験室自動化（1 種）**
  - **scientific-lab-automation** スキル (#44): PyLabRobot ユニバーサル API・Opentrons OT-2 プロトコル・SOP テンプレート・Protocols.io API・ELN/LIMS 連携・QC 検証ワークフロー
- **N. 科学プレゼンテーション（1 種）**
  - **scientific-presentation-design** スキル (#45): 15 スライド口頭発表テンプレート (時間配分付)・tikzposter テンプレート・matplotlib ワークフロー図・3m ルール・アクセシビリティ (コントラスト比 ≥4.5)
- **O. 研究計画・グラント（2 種）**
  - **scientific-grant-writing** スキル (#46): NIH Specific Aims 1 ページテンプレート・Research Strategy (Significance/Innovation/Approach 12 頁)・JSPS 科研費フォーマット・予算計画・Budget Justification
  - **scientific-research-methodology** スキル (#47): SCAMPER/Cross-Domain/TRIZ ブレインストーミング・仮説型分類・研究デザインマトリクス・FINER 基準・バイアス分類 (6 種)・IRB/倫理チェックリスト

### Changed
- **README.md**: スキル数を 36→47 に更新、カテゴリ数を 9→15 に拡張 (J-O 追加)、ディレクトリ構造に 6 新カテゴリを追加

## [0.6.0] - 2026-02-12

### Added
- **scientific-deep-research** スキル (#36): SHIKIGAMI の WebResearcher パラダイムを科学研究に適応した深層リサーチスキル
  - Think→Search→Evaluate→Synthesize 反復サイクル（最大 15 ラウンド）
  - 学術データベース検索統合（PubMed, Google Scholar, arXiv, Semantic Scholar, CiNii, J-STAGE 等）
  - エビデンス階層評価（Level 1a〜5 + プレプリント）
  - ソース信頼性スコアリング（IF, h-index, サンプルサイズ, 統計手法, 再現性, 出版年）
  - ハルシネーション防止マーキング（✅ 検証済 / 📎 単一ソース / ⚠️ 未査読 / ❓ AI推定 / 🔄 古いデータ / ⚡ 矛盾）
  - 交差検証（数値乖離率判定・内容矛盾処理）
  - 品質ゲート（Phase 1→2, 2→3, 3→4, 4 完了）
  - PICO/PECO/SPIDER 構造化・検索戦略設計
  - PRISMA 2020 フローテンプレート（Systematic Review 用）
  - 批判的評価チェックリスト（Cochrane RoB / NOS / AMSTAR-2 準拠）
  - 日英並列検索必須・学術ドメイン優先
  - 分野別検索戦略テンプレート（生命科学・材料科学・CS/AI）
  - 出力ファイル: research_report.md, evidence_table.json, search_log.md, source_registry.json, prisma_flow.md

### Changed
- **README.md**: スキル数を 35→36 に更新、カテゴリ I (Deep Research) を追加、ディレクトリ構造に [I] Deep Research を追加

## [0.5.2] - 2026-02-12

### Fixed
- **scientific-academic-writing** スキル: SATORI バージョン参照ルールを追加 — 論文原稿内でバージョン番号をハードコードせず `package.json` から動的取得するよう明記
- **scientific-academic-writing** スキル: AI 使用開示 (AI Usage Disclosure) セクションを追加 — ジャーナル別記載場所・テンプレート・Authors セクション禁止ルールを整備

## [0.5.0] - 2026-02-12

### Added
- **scientific-peer-review-response** スキル (#11): 査読コメント構造化 (Major/Minor/Editorial 自動分類)・ポイントバイポイント回答レター生成・改訂版カバーレター・コメント-改訂マッピング
- **scientific-revision-tracker** スキル (#12): 改訂履歴追跡・バージョンスナップショット・差分検出・変更マークアップ (Markdown/LaTeX/プレーンテキスト)・トレーサビリティ検証
- **scientific-paper-quality** スキル (#13): 可読性スコア (Flesch-Kincaid/Gunning Fog)・IMRAD バランス分析・語彙品質 (冗長表現/過剰主張検出)・ジャーナル適合性チェック・再現可能性チェック・総合品質スコアカード

### Changed
- **README.md**: スキル数を 32→35 に更新、カテゴリ A を 10→13 種に拡張、パイプラインフローに査読対応・改訂追跡・品質評価フェーズを追加
- スキル番号を再付番 (B-H カテゴリ: #11-32 → #14-35)

## [0.4.0] - 2025-06-13

### Added
- **scientific-supplementary-generator** スキル (#8): Supplementary Information 自動生成・SI 図表整理・本文-SI 相互参照検証・ジャーナル別 SI フォーマット対応 (Nature/Science/ACS/IEEE/Elsevier)
- **scientific-latex-formatter** スキル (#9): Markdown→LaTeX 変換・ジャーナルテンプレート適用 (Nature/Science/ACS/IEEE/Elsevier/RevTeX)・数式/図表/引用の LaTeX 構文変換・BibTeX 生成
- **scientific-citation-checker** スキル (#10): 引用文献の自動検索・網羅性チェック・整合性検証 (孤立参考文献/未解決引用/重複検出)・引用カバレッジ分析

### Changed
- **README.md**: スキル数を 29→32 に更新、カテゴリ A を 7→10 種に拡張、パイプラインフローテーブルに引用検証・SI 生成・LaTeX 変換フェーズを追加
- スキル番号を再付番 (B-H カテゴリ: #8-29 → #11-32)

## [0.3.0] - 2026-02-11

### Added
- **scientific-critical-review** スキル (#29): 草稿の批判的レビュー・5パスレビューワークフロー・7層考察深化フレームワーク・修正案生成
- **scientific-hypothesis-pipeline**: 仮説定義の永続化 (`docs/hypothesis.{md,json}`)、ワークフロー設計の永続化 (`docs/workflow_design.{md,json}`)、`load_hypothesis()` / `load_workflow_design()` 関数
- **scientific-critical-review**: レビュー結果の永続化 (`manuscript/review_report.{md,json}`, `manuscript/manuscript_revised.md`, `manuscript/revision_diff.md`)、`run_review_pipeline()` 統合関数
- **scientific-pipeline-scaffold**: Step 0 (仮説・ワークフロー保存) / Step 8 (論文執筆・レビュー) を追加、参照スキル表を追加
- **scientific-academic-writing**: 参照スキル表・パイプラインフロー図を追加、`docs/hypothesis.md` / `results/analysis_summary.json` からの自動参照を明記
- 4スキル間の双方向相互参照を完備（hypothesis-pipeline ↔ critical-review ↔ academic-writing ↔ pipeline-scaffold）

### Fixed
- **scientific-critical-review**: 全 save 関数のパス方式を `Path()` (CWD相対) から `BASE_DIR` ベースに統一
- **README.md**: パイプラインフロー図追加、インストール方法追加、ディレクトリ構造に hypothesis-pipeline / critical-review を追加

## [0.2.0] - 2025-06-12

### Added
- **scientific-academic-writing**: 図の埋め込みワークフロー（セクション 8）を追加
- 全 6 ジャーナルテンプレートに `![Figure N](figures/...)` 構文を追加
  - IMRaD 標準、Nature、Science、ACS、IEEE、Elsevier

## [0.1.0] - 2025-06-12

### Added
- 初回リリース: 27 スキルを含む Agent Skills コレクション
- CLI ツール (`npx @nahisaho/satori init`) による `.github/skills/` のインストール
- `--force` / `--dry-run` オプション対応
- 8 カテゴリ (A-H) のスキル分類体系
- 7 種のジャーナルテンプレート (IMRaD, Nature, Science, ACS, IEEE, Elsevier, Qiita)

[0.8.0]: https://github.com/nahisaho/satori/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/nahisaho/satori/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/nahisaho/satori/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/nahisaho/satori/compare/v0.5.2...v0.6.0
[0.5.2]: https://github.com/nahisaho/satori/compare/v0.5.0...v0.5.2
[0.5.0]: https://github.com/nahisaho/satori/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/nahisaho/satori/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/nahisaho/satori/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/nahisaho/satori/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/nahisaho/satori/releases/tag/v0.1.0
