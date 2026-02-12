---
name: scientific-deep-research
description: |
  科学文献の深層リサーチスキル。SHIKIGAMI の WebResearcher パラダイム
  （Think→Report→Action 反復サイクル）を科学研究に特化させた実装。
  学術データベース検索、エビデンス階層評価、ソース追跡、交差検証、
  ハルシネーション防止を統合した反復的深層調査を提供。
  「文献調査して」「先行研究を調べて」「systematic review して」で発火。
---

# Scientific Deep Research

科学文献の反復的深層リサーチを支援するスキル。SHIKIGAMI の WebResearcher
パラダイム（Think→Report→Action サイクル）を科学研究コンテキストに適応し、
学術データベース検索・エビデンス階層評価・ソース追跡・交差検証を統合する。

## When to Use

- 科学研究テーマの先行研究調査を行うとき
- Systematic Review / Scoping Review の文献検索を行うとき
- 特定の研究トピックの最新動向を網羅的に調査するとき
- 技術比較・手法比較のための情報収集が必要なとき
- Grant proposal / Research proposal の背景調査をするとき
- 論文 Introduction のための文献レビューを行うとき

## Quick Start

### 1. Deep Research ワークフロー概要

```
Phase 1: Research Question Definition（研究課題の定義）
    ↓
Phase 2: Deep Research（反復的深層リサーチ）
    ┌──────────────────────────────────────────┐
    │  Think → Search → Evaluate → Synthesize  │ ← 反復サイクル
    │         (最大 15 ラウンド)                │
    └──────────────────────────────────────────┘
    ↓
Phase 3: Evidence Synthesis（エビデンス統合）
    ↓
Phase 4: Research Report（リサーチレポート生成）
```

### 2. WebResearcher パラダイム（科学研究適応版）

SHIKIGAMIの Think→Report→Action サイクルを科学研究に最適化:

```
┌─────────────────────────────────────────────────────┐
│  Round N                                            │
│  ┌─────────────────────────────────────────────┐   │
│  │ Workspace (前ラウンドの Synthesis + 新文献)   │   │
│  └─────────────────────────────────────────────┘   │
│           ↓                                         │
│  ┌───────────┐  ┌───────────┐  ┌───────────────┐  │
│  │  Think    │→│  Search   │→│  Evaluate      │  │
│  │ (知識     │  │ (学術DB   │  │ (エビデンス   │  │
│  │  ギャップ │  │  検索)    │  │  評価)        │  │
│  │  分析)    │  │           │  │               │  │
│  └───────────┘  └───────────┘  └───────┬───────┘  │
│                                         ↓          │
│                                ┌───────────────┐   │
│                                │  Synthesize   │   │
│                                │ (進化する     │   │
│                                │  レビュー)    │   │
│                                └───────────────┘   │
└─────────────────────────────────────────────────────┘
```

---

## Phase 1: Research Question Definition

### PICO/PECO 構造化

リサーチクエスチョンを構造化し、検索戦略の基盤を構築する。

```markdown
## Research Question の構造化

### PICO (介入研究)
- **P** (Population): 対象集団・研究対象
- **I** (Intervention): 介入・手法・技術
- **C** (Comparison): 比較対象・既存手法
- **O** (Outcome): 評価指標・アウトカム

### PECO (観察研究)
- **P** (Population): 対象集団
- **E** (Exposure): 曝露・要因
- **C** (Comparison): 比較群
- **O** (Outcome): アウトカム

### SPIDER (質的研究)
- **S** (Sample): サンプル
- **PI** (Phenomenon of Interest): 関心現象
- **D** (Design): 研究デザイン
- **E** (Evaluation): 評価方法
- **R** (Research type): 研究種別
```

### 検索戦略の設計

```markdown
## 検索戦略テンプレート

### 1. キーワード設計
| 概念 | 日本語キーワード | 英語キーワード | MeSH / 統制語彙 |
|------|-----------------|---------------|-----------------|
| P    | [対象]          | [population]  | [MeSH term]     |
| I/E  | [介入/曝露]     | [intervention]| [MeSH term]     |
| C    | [比較]          | [comparison]  | [MeSH term]     |
| O    | [アウトカム]    | [outcome]     | [MeSH term]     |

### 2. ブール演算子の組み合わせ
(P terms) AND (I terms) AND (O terms)

### 3. 包含・除外基準
| 基準 | 包含 | 除外 |
|------|------|------|
| 出版年 | [YYYY]-present | < [YYYY] |
| 言語 | English, Japanese | Others |
| 研究デザイン | [RCT, Cohort, etc.] | [Case report, etc.] |
| 出版種別 | Original article, Review | Editorial, Letter |
```

---

## Phase 2: Deep Research（反復的深層リサーチ）

### Think フェーズ（知識ギャップ分析）

各ラウンドの冒頭で現在の知識状態を評価し、次のアクションを決定する。

```markdown
## Think フェーズのチェックリスト

### 知識状態の評価
- [ ] 現在収集済みの文献数とカバレッジ
- [ ] 未解決の研究課題・サブクエスチョン
- [ ] エビデンスが不足している領域
- [ ] 矛盾する知見の有無

### 次アクションの決定
- SEARCH: 知識ギャップが大きい → 追加検索
- VISIT: 有望な論文を発見 → 詳細取得
- VERIFY: 重要な知見 → 別ソースで交差検証
- COMPLETE: 十分な情報収集 → 統合フェーズへ
```

**重要**: Think フェーズの内部推論は次ラウンドに引き継がない
（コンテキスト汚染防止、SHIKIGAMI 準拠）。

### Search フェーズ（学術データベース検索）

#### 検索対象データベース

| データベース | 対象分野 | 検索方法 | 優先度 |
|-------------|---------|----------|--------|
| **PubMed / MEDLINE** | 生命科学・医学 | MeSH + Free text | P0 (生命科学) |
| **Google Scholar** | 全分野 | Free text | P0 (全般) |
| **arXiv** | 物理・CS・数学・生物 | Free text | P0 (物理・CS) |
| **Semantic Scholar** | 全分野 (AI 解析付き) | API / Free text | P1 |
| **Web of Science** | 全分野 (IF 付き) | Topic search | P1 |
| **Scopus** | 全分野 | Title-Abs-Key | P1 |
| **CrossRef** | DOI メタデータ | DOI lookup | P2 |
| **CiNii** | 日本語文献 | Free text | P1 (日本語) |
| **J-STAGE** | 日本語学術誌 | Free text | P1 (日本語) |
| **ChemRxiv** | 化学 | Free text | P2 (化学) |
| **bioRxiv / medRxiv** | プレプリント | Free text | P2 (生命科学) |

#### 検索ルール（必須）

```markdown
## 検索の実行規則

### 日英並列検索（必須）
- 全てのクエリは日本語・英語の両方で実行する
- 例: ["ZnO薄膜 スパッタリング 特性", "ZnO thin film sputtering properties"]

### 学術ドメイン優先
- 以下のドメインを優先的に評価:
  .ac.jp, .edu, .gov, .go.jp, .org (学術機関)
  pubmed.ncbi.nlm.nih.gov
  scholar.google.com
  arxiv.org
  doi.org

### 検索式の構造化
- ブール演算子 (AND, OR, NOT) を活用
- ワイルドカード (*) で語幹展開
- フレーズ検索 ("exact phrase") で精度向上
- フィルタ: 出版年、研究デザイン、言語

### 失敗時のフォールバック
- 404/5xx エラー → Wayback Machine → Archive.today
- 検索結果 0 件 → クエリ簡略化 → 同義語展開 → 上位概念へ
```

### Evaluate フェーズ（エビデンス評価）

収集した各文献・情報源のエビデンスレベルと信頼性を評価する。

#### エビデンス階層（Evidence Hierarchy）

```markdown
## エビデンスレベル分類

| Level | 研究デザイン | 信頼度 | 記号 |
|-------|-------------|--------|------|
| 1a | Systematic Review of RCTs | ★★★★★ | [SR] |
| 1b | Individual RCT | ★★★★☆ | [RCT] |
| 2a | Systematic Review of Cohort | ★★★★☆ | [SR-C] |
| 2b | Individual Cohort Study | ★★★☆☆ | [Cohort] |
| 3a | Systematic Review of Case-Control | ★★★☆☆ | [SR-CC] |
| 3b | Individual Case-Control Study | ★★☆☆☆ | [CC] |
| 4 | Case Series / Cross-sectional | ★★☆☆☆ | [CS] |
| 5 | Expert Opinion / Narrative Review | ★☆☆☆☆ | [EO] |
| — | Preprint (not peer-reviewed) | ☆☆☆☆☆ | [PP] |
```

#### ソース信頼性スコアリング

```markdown
## 信頼性評価マトリクス

| 評価項目 | スコア範囲 | 重み |
|---------|-----------|------|
| **ジャーナル IF / ランク** | 0-25 | ×2 |
| **著者の h-index / 実績** | 0-25 | ×1.5 |
| **サンプルサイズ n** | 0-25 | ×1.5 |
| **統計手法の適切性** | 0-25 | ×1 |
| **再現性・データ公開** | 0-25 | ×1 |
| **出版年の新しさ** | 0-25 | ×0.5 |

### 総合信頼度スコア
- 80-100: ★★★★★ 高信頼
- 60-79:  ★★★★☆ 信頼
- 40-59:  ★★★☆☆ 中程度
- 20-39:  ★★☆☆☆ 低信頼
- 0-19:   ★☆☆☆☆ 要注意

### 信頼度不足時のアクション
- スコア 60 未満 → 交差検証必須（別ソースで確認）
- プレプリント → 「⚠️ 未査読」マーキング必須
- 5 年以上前 → 「⚠️ 要更新確認」マーキング推奨
```

#### 批判的評価チェックリスト

```markdown
## 文献の批判的評価 (Critical Appraisal)

### RCT 評価 (Cochrane Risk of Bias 準拠)
- [ ] ランダム化の方法は適切か
- [ ] 割り付けの隠蔽化はされているか
- [ ] ブラインドは適切か（参加者・評価者）
- [ ] 欠測データの処理は適切か
- [ ] アウトカムの選択的報告はないか
- [ ] その他のバイアスリスクはないか

### 観察研究評価 (Newcastle-Ottawa Scale 準拠)
- [ ] 対象群の選択は適切か
- [ ] 群間の比較可能性はあるか
- [ ] アウトカムの測定は適切か
- [ ] フォローアップは十分か

### 系統的レビュー評価 (AMSTAR-2 準拠)
- [ ] プロトコルは事前登録されているか
- [ ] 包括的な文献検索が行われているか
- [ ] バイアスリスク評価が行われているか
- [ ] メタアナリシスの手法は適切か
```

### Synthesize フェーズ（進化するレビュー）

各ラウンドで収集した情報を統合し、進化するリサーチレポートを更新する。

```markdown
## Evolving Research Report（進化するリサーチレポート）

### Round N の更新内容
| 要素 | 内容 |
|------|------|
| **Rationale（根拠）** | 研究課題に直接関連する具体的知見の特定 |
| **Evidence（証拠）** | 最も関連性の高い情報の抽出（原文保持優先） |
| **Summary（要約）** | 現時点での知見の統合・整理 |
| **Confidence（確信度）** | 0-100% — 十分な情報が集まったかの指標 |
| **Gaps（ギャップ）** | 残存する知識ギャップの一覧 |

### 完了判定基準
- [ ] 主要な研究課題が回答されている
- [ ] 最低 5 件の査読済み文献を収集
- [ ] 日英両方で検索が実行されている
- [ ] 全てのデータに出典 URL/DOI が付与されている
- [ ] 信頼度スコア 60% 以上のソースが 3 件以上
- [ ] 矛盾する知見がある場合、交差検証が完了している
```

---

## Phase 3: Evidence Synthesis（エビデンス統合）

### PRISMA フロー（Systematic Review 用）

```markdown
## PRISMA 2020 フローテンプレート

### Identification（同定）
- データベース検索: n = ___
  - PubMed: n = ___
  - Google Scholar: n = ___
  - arXiv: n = ___
  - その他: n = ___
- その他の情報源: n = ___
  - 引用文献追跡: n = ___
  - 被引用文献追跡: n = ___
  - 手動検索: n = ___

### Screening（スクリーニング）
- 重複除去後: n = ___
- タイトル/抄録スクリーニング: n = ___
  - 除外: n = ___ (理由: ___)
- 全文スクリーニング: n = ___
  - 除外: n = ___ (理由: ___)

### Included（採択）
- 質的統合に含めた研究: n = ___
- 量的統合（メタアナリシス）に含めた研究: n = ___
```

### 統合分析テンプレート

```markdown
## エビデンス統合テーブル

| # | 著者 (年) | 研究デザイン | n | 主要結果 | 信頼度 | DOI |
|---|----------|-------------|---|---------|--------|-----|
| 1 | Author (2024) | [RCT] | 100 | [key finding] | ★★★★☆ | 10.xxxx |
| 2 | Author (2023) | [Cohort] | 500 | [key finding] | ★★★☆☆ | 10.xxxx |
| ... | | | | | | |

## 知見の統合
### 一致する知見
- [複数の研究で支持されている知見を記述]

### 矛盾する知見
- [文献 A] は [知見 X] を報告しているが、[文献 B] は [知見 Y] を報告
  - 考えられる理由: [方法論の差異 / 対象集団の違い / ...]

### エビデンスギャップ
- [十分な研究がない領域を特定]
```

---

## Phase 4: Research Report（リサーチレポート生成）

### レポートテンプレート

```markdown
# Deep Research Report: [研究テーマ]

**生成日**: YYYY-MM-DD
**リサーチラウンド数**: N 回
**収集文献数**: M 件
**主要データベース**: [PubMed, Google Scholar, ...]

---

## Executive Summary

[3-5 文で主要な知見を要約]

## 1. Research Question

### 1.1 構造化された研究課題
[PICO/PECO 形式で記述]

### 1.2 検索戦略
[実行した検索式・フィルタの記録]

## 2. Findings

### 2.1 [サブトピック 1]
[知見の記述 + 引用]

### 2.2 [サブトピック 2]
[知見の記述 + 引用]

### 2.3 [サブトピック 3]
[知見の記述 + 引用]

## 3. Evidence Summary

### 3.1 エビデンステーブル
[統合分析テーブル]

### 3.2 エビデンスの強さ
[GRADE アプローチによる評価]

## 4. Knowledge Gaps & Future Directions

[特定されたギャップと今後の研究方向]

## 5. References

[全引用文献リスト — DOI/URL 必須]
```

---

## ソース追跡・ハルシネーション防止

### ソース追跡ルール（必須）

```markdown
## 全データにソース情報を付与する

### 必須項目
| 項目 | 説明 | 例 |
|------|------|---|
| **著者** | 筆頭著者 et al. | Smith et al. |
| **年** | 出版年 | 2024 |
| **ジャーナル** | 掲載誌名 | Nature |
| **DOI** | Digital Object Identifier | 10.1038/s41586-... |
| **URL** | アクセス用 URL | https://doi.org/10.1038/... |
| **アクセス日** | 情報取得日 | 2026-02-12 |
| **エビデンスレベル** | 研究デザイン | [RCT] |
| **信頼度スコア** | 評価結果 | ★★★★☆ (78/100) |
```

### ハルシネーション防止マーキング

```markdown
## 情報の確度マーキング（必須）

| マーカー | 意味 | 使用条件 |
|---------|------|---------|
| ✅ | 検証済み（2+ ソースで確認） | 複数の独立したソースで確認 |
| 📎 | 単一ソース | 1 つのソースからの情報 |
| ⚠️ | 未査読・要確認 | プレプリント、灰色文献 |
| ❓ | AI による推定・補完 | ソースなし、推測を含む |
| 🔄 | 古いデータ（5 年以上） | 最新情報で要更新確認 |
| ⚡ (?) | 情報源間の矛盾 | 複数ソースで矛盾する情報 |

### 使用例
- ✅ ZnO 薄膜のバンドギャップは 3.37 eV である [Smith 2024, Tanaka 2023]
- 📎 新規触媒の変換効率は 95% と報告されている [Lee 2024]
- ⚠️ 量子コンピュータによる最適化で 10x 高速化（プレプリント）[arXiv:2024.xxxxx]
- ❓ この方法は他の金属酸化物にも適用可能と考えられる
- ⚡ (?) 市場規模は $50B [Report A] vs $35B [Report B]
```

---

## 交差検証（Cross-Validation）

```markdown
## 交差検証の実行ルール

### 必須条件
- 重要な数値データ → 最低 2 ソースで確認
- 主要な主張 → 独立した研究で確認
- 矛盾発見時 → 第 3 のソースで検証

### 数値データの矛盾処理
| 乖離率 | 重大度 | アクション |
|--------|--------|----------|
| < 5% | 低 | 範囲として記述 (e.g., 3.35-3.37 eV) |
| 5-20% | 中 | ⚡ (?) マーク、測定条件の差異を注記 |
| > 20% | 高 | ⚡ (?) マーク、追加調査、条件差異の分析 |

### 内容・見解の矛盾処理
1. 両方の見解を併記する
2. 各見解のエビデンスレベルを明示する
3. 矛盾の考えられる要因を分析する
   - 研究デザインの違い
   - 対象集団の違い
   - 測定方法の違い
   - 出版バイアスの可能性
```

---

## 品質ゲート（Quality Gates）

```markdown
## Deep Research 完了チェックリスト

### Phase 1→2 ゲート（研究課題定義→リサーチ開始）
- [ ] PICO/PECO が構造化されている
- [ ] 検索キーワードが日英両方で定義されている
- [ ] 包含・除外基準が明確である

### Phase 2→3 ゲート（リサーチ→統合）
- [ ] 最低 5 件の査読済み文献を収集
- [ ] 全てのデータに DOI/URL が付与されている
- [ ] 日英両方で検索が実行されている
- [ ] 信頼度スコア 60% 以上のソースが 3 件以上
- [ ] 主要な知見に対する交差検証が完了

### Phase 3→4 ゲート（統合→レポート）
- [ ] エビデンステーブルが完成している
- [ ] 矛盾する知見が整理されている
- [ ] 知識ギャップが特定されている

### Phase 4 完了ゲート（レポート完成）
- [ ] 全引用に DOI/URL が記載されている
- [ ] ハルシネーション防止マーカーが付与されている
- [ ] PRISMA フロー（該当する場合）が作成されている
- [ ] Executive Summary が本文と整合している
```

---

## 検索戦略テンプレート（分野別）

### 生命科学・医学

```markdown
## PubMed 検索戦略

### MeSH + Free text の組み合わせ
("Drug Delivery Systems"[MeSH] OR "drug delivery"[tiab])
AND
("Nanoparticles"[MeSH] OR nanoparticle*[tiab])
AND
("Neoplasms"[MeSH] OR cancer[tiab] OR tumor[tiab])

### フィルタ
- Publication date: 2020/01/01 - present
- Article type: Clinical Trial, Randomized Controlled Trial, Review
- Language: English, Japanese
```

### 材料科学・化学

```markdown
## Google Scholar / Scopus 検索戦略

### キーワード組み合わせ
("ZnO thin film" OR "zinc oxide film")
AND
(sputtering OR "physical vapor deposition" OR PVD)
AND
(optical OR electrical OR "band gap")

### 年代フィルタ: 2020-present
### 引用数ソート: 被引用数降順で上位 50 件を確認
```

### コンピュータサイエンス・AI

```markdown
## arXiv / Semantic Scholar 検索戦略

### キーワード
("large language model" OR LLM)
AND
("retrieval augmented generation" OR RAG)
AND
(evaluation OR benchmark)

### カテゴリフィルタ: cs.CL, cs.AI, cs.IR
### ソート: 最新順 + 被引用数
```

---

## 他スキルとの連携

### パイプライン上の位置づけ

```
deep-research → hypothesis-pipeline → pipeline-scaffold → academic-writing
  (文献調査)    (仮説定義)          (解析実行)         (草稿作成)
```

### 連携スキル

| スキル | 連携内容 |
|--------|---------|
| `scientific-hypothesis-pipeline` | Deep Research の結果から仮説を定式化 |
| `scientific-academic-writing` | 収集した文献を Introduction / Discussion に反映 |
| `scientific-citation-checker` | 引用の網羅性・整合性を検証 |
| `scientific-critical-review` | リサーチレポートの批判的レビュー |
| `scientific-meta-analysis` | 収集した文献データでメタアナリシスを実行 |
| `scientific-statistical-testing` | 数値データの統計的検証 |

---

## Output Files

| ファイル | 形式 | 内容 |
|---------|------|------|
| `research/research_report.md` | Markdown | Deep Research レポート |
| `research/evidence_table.json` | JSON | エビデンス統合テーブル |
| `research/search_log.md` | Markdown | 検索実行ログ（検索式・結果数・日時） |
| `research/source_registry.json` | JSON | ソース追跡レジストリ |
| `research/prisma_flow.md` | Markdown | PRISMA フロー（該当時） |

### research_report.md のフォーマット

```markdown
---
title: "[研究テーマ]"
date: "YYYY-MM-DD"
rounds: N
sources_count: M
databases: ["PubMed", "Google Scholar", ...]
confidence_score: XX%
---

# Deep Research Report: [研究テーマ]
...
```

### evidence_table.json のスキーマ

```json
{
  "research_question": "...",
  "search_strategy": {
    "databases": ["PubMed", "Google Scholar"],
    "date_range": "2020-2026",
    "keywords": {
      "ja": ["キーワード1", "キーワード2"],
      "en": ["keyword1", "keyword2"]
    }
  },
  "sources": [
    {
      "id": "S001",
      "authors": "Smith et al.",
      "year": 2024,
      "title": "...",
      "journal": "Nature",
      "doi": "10.1038/...",
      "url": "https://doi.org/...",
      "evidence_level": "1b",
      "study_design": "RCT",
      "sample_size": 100,
      "key_findings": ["..."],
      "credibility_score": 85,
      "accessed_date": "2026-02-12"
    }
  ],
  "synthesis": {
    "consistent_findings": ["..."],
    "conflicting_findings": ["..."],
    "evidence_gaps": ["..."]
  }
}
```

### source_registry.json のスキーマ

```json
{
  "registry_version": "1.0",
  "created_at": "2026-02-12T00:00:00+09:00",
  "sources": [
    {
      "id": "S001",
      "url": "https://doi.org/...",
      "doi": "10.1038/...",
      "accessed_at": "2026-02-12T10:30:00+09:00",
      "retrieval_method": "search",
      "query_used": "ZnO thin film sputtering",
      "database": "Google Scholar",
      "status": "active",
      "credibility_score": 85,
      "evidence_level": "2b",
      "cross_validated": true,
      "cross_validation_sources": ["S003", "S007"]
    }
  ]
}
```

## References

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PubMed | `PubMed_search_articles` | 文献検索 |
| PubMed | `PubMed_get_cited_by` | 被引用構造分析 |
| EuropePMC | `EuropePMC_search_articles` | ヨーロッパ文献検索 |
| Crossref | `Crossref_search_works` | 出版情報検索 |
| OpenAlex | `OpenAlex_Guidelines_Search` | オープンアクセス文献 |
| ArXiv | `ArXiv_search_papers` | プレプリント検索 |
