# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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

[0.4.0]: https://github.com/nahisaho/satori/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/nahisaho/satori/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/nahisaho/satori/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/nahisaho/satori/releases/tag/v0.1.0
