# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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

[0.6.0]: https://github.com/nahisaho/satori/compare/v0.5.2...v0.6.0
[0.5.2]: https://github.com/nahisaho/satori/compare/v0.5.0...v0.5.2
[0.5.0]: https://github.com/nahisaho/satori/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/nahisaho/satori/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/nahisaho/satori/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/nahisaho/satori/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/nahisaho/satori/releases/tag/v0.1.0
