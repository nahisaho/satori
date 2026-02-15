---
name: scientific-academic-writing
description: |
  科学技術・学術論文の執筆スキル。IMRaD 標準、Nature/Science 系、ACS 系、IEEE 系、
  Elsevier 系のジャーナル形式に対応した論文構成・セクション設計・文章パターンを提供。
  「論文を書いて」「Abstract を作成して」「Methods セクションを書いて」で発火。
  assets/ に主要ジャーナル形式の Markdown テンプレートを同梱。
tu_tools:
  - key: crossref
    name: Crossref
    description: 論文メタデータ検索・引用情報取得
---

# Scientific Academic Writing

科学技術・学術論文の執筆を支援するスキル。ジャーナル形式に応じた構成テンプレート、
セクション別の文章パターン、引用・図表参照の規約を提供する。

## SATORI バージョン参照ルール

**重要**: 論文原稿内で SATORI のバージョンに言及する場合（Authors, Acknowledgements,
AI Disclosure 等）、バージョン番号をハードコードしてはならない。
必ず `package.json` の `version` フィールドを読み取り、実際の値を使用すること。

```
✅ 正しい例: GitHub Copilot Agent (SATORI v0.5.1)  ← package.json から取得
❌ 誤った例: GitHub Copilot Agent (SATORI v0.2.0)  ← 古いバージョンのハードコード
```

取得手順:
1. ワークスペースの `package.json` を読み取る
2. `version` フィールドの値を取得する
3. 原稿内で `SATORI v{version}` の形式で記載する

## When to Use

- 学術論文の草稿を作成するとき
- Abstract / Introduction / Methods / Results / Discussion を執筆するとき
- 特定ジャーナル形式に合わせた論文構成が必要なとき
- Cover Letter / Response to Reviewers を作成するとき
- 既存原稿のセクション構成を見直すとき

## Quick Start

### 1. テンプレート選択ガイド

```
ジャーナル形式の選択フロー:

Q1: 投稿先は？
├── Nature / Nature Communications / Nature Materials → nature_article.md
├── Science / Science Advances                       → science_research_article.md
├── ACS Nano / JACS / Chem. Mater.                  → acs_article.md
├── IEEE Trans. / IEEE Access                        → ieee_transactions.md
├── Elsevier 系 (Acta Mater., etc.)                 → elsevier_article.md
├── Qiita 技術記事                                   → qiita_technical_article.md
└── その他 / 不明                                    → imrad_standard.md (最も汎用的)
```

### 2. 論文構成の基本原則

```markdown
## 論文構成の CARS モデル (Create A Research Space)

### Introduction の 3 ステップ:
1. **Establishing a territory** — 研究分野の重要性・先行研究の概観
2. **Establishing a niche** — 先行研究のギャップ・未解決の問題
3. **Occupying the niche** — 本研究の目的・アプローチ・主要な発見

### Discussion の構成:
1. 主要な発見の要約（Results の再述ではなく解釈）
2. 先行研究との比較・位置づけ
3. メカニズムの考察・理論的意味
4. 研究の限界 (Limitations)
5. 将来の展望 (Future perspectives)
6. 結論 (Conclusion) — 一部ジャーナルでは独立セクション
```

### 3. セクション別ライティングパターン

#### Abstract（構造化抄録）

```markdown
## Abstract テンプレート（250 words 以内が一般的）

**Background/Context**: [研究分野] において [課題/問題] は依然として [重要な課題] である。

**Objective/Purpose**: 本研究では、[手法/アプローチ] を用いて [研究目的] を達成することを
目的とした。

**Methods**: [材料/データ] に対して [実験手法/解析手法] を適用した。
[主要なパラメータ/条件] は [値] とした。

**Results**: [主要な定量的結果を 2-3 文で記述]。
[統計的有意性: p < 0.05, 効果量, 信頼区間などを含む]。

**Conclusions**: これらの結果は [解釈/意義] を示唆しており、
[応用/今後の展望] に有用である。

**Keywords**: keyword1, keyword2, keyword3, keyword4, keyword5
```

#### Introduction

```markdown
## Introduction テンプレート（800-1500 words）

### 第 1 段落: 研究分野の重要性
[研究分野] は [産業/科学的意義] において重要な役割を果たしている [ref1, ref2]。
特に [特定のトピック] は [理由] から注目を集めている。

### 第 2 段落: 先行研究の概観
これまでに [研究グループ A] は [手法/発見] を報告しており [ref3]、
[研究グループ B] は [別の手法/発見] を示した [ref4]。
さらに [研究グループ C] による [成果] は [意義] を明らかにした [ref5, ref6]。

### 第 3 段落: ギャップの特定
しかしながら、[未解決の課題 1] や [未解決の課題 2] は依然として明らかにされていない。
特に [具体的なギャップ] については [理由] から十分な検討がなされていない。

### 第 4 段落: 本研究の目的
本研究では、[アプローチ/手法] を用いて [目的] を明らかにすることを目指した。
具体的には、(1) [サブ目的 1]、(2) [サブ目的 2]、(3) [サブ目的 3] を検討した。
```

#### Methods / Experimental

```markdown
## Methods テンプレート

### 2.1 Materials / Datasets
[材料名] (純度 XX%, [メーカー名], [国]) を使用した。
[データセット名] は [出典/生成方法] から取得した (n = XXX)。

### 2.2 Experimental Procedure / Data Processing
[装置名] ([型番], [メーカー名]) を用いて [条件] で [操作] を行った。
[パラメータ 1] は [範囲/値]、[パラメータ 2] は [範囲/値] とした。

### 2.3 Characterization / Analysis
[分析手法] ([装置型番]) を用いて [測定対象] を評価した。
測定条件は [条件の詳細] とした。

### 2.4 Statistical Analysis
統計解析には [ソフトウェア名] (ver. X.X) を使用した。
群間比較には [検定名] を用い、有意水準は p < 0.05 とした。
多重比較には [補正法] を適用した。
```

#### Results

```markdown
## Results テンプレート

### 3.1 [実験/解析 1 の結果]
**図表の導入**: Figure 1a に [測定対象] の [可視化内容] を示す。
**定量的記述**: [パラメータ] は [値 ± SD] (n = XX) であった。
**比較・傾向**: [条件 A] と比較して [条件 B] では [XX]% の [増加/減少] が
観察された (p = X.XXX, [検定名])。
**図表の参照**: この傾向は Figure 1b の [内容] からも確認される。

### 注意: Results での禁止事項
- ❌ データの解釈・推測（→ Discussion で記述）
- ❌ 先行研究との比較（→ Discussion で記述）
- ❌ 方法の説明（→ Methods で記述）
- ✅ 客観的な事実と数値のみ記述
```

#### Discussion

```markdown
## Discussion テンプレート

### 第 1 段落: 主要な発見の要約
本研究の主要な発見は以下の 3 点である:
(1) [発見 1 — 最も重要な結果]、
(2) [発見 2]、(3) [発見 3]。

### 第 2 段落: 先行研究との比較
[発見 1] は [先行研究の結果] と一致しており [refX]、
[メカニズム/理論] を支持する結果である。
一方、[発見 2] については [先行研究 Y] の報告 ([refY]) と
異なる結果が得られた。この差異は [考えられる理由] に
起因すると考えられる。

### 第 3 段落: メカニズムの考察
[観察された現象] のメカニズムとして、
[仮説/モデル] が考えられる (Figure X)。
[根拠 1] および [根拠 2] がこの解釈を支持している。

### 第 4 段落: 限界と展望
本研究にはいくつかの限界がある。
第一に [限界 1]、第二に [限界 2] である。
今後 [将来の研究方向] を検討することで、
[期待される成果] が得られると考える。

### 第 5 段落: 結論
以上の結果から、[主要な結論] が示された。
本研究は [分野] における [貢献/意義] を提供するものである。
```

### 4. 図表の参照規約

```markdown
## 図表参照のジャーナル別書式

| ジャーナル | 図の参照 | 表の参照 | 補足図 |
|---|---|---|---|
| Nature 系 | Fig. 1a | Table 1 | Extended Data Fig. 1 |
| Science 系 | Fig. 1A | Table 1 | fig. S1 |
| ACS 系 | Figure 1a | Table 1 | Figure S1 |
| IEEE 系 | Fig. 1(a) | TABLE I | — |
| Elsevier 系 | Fig. 1(a) | Table 1 | Fig. S1 |
| IMRaD 標準 | Figure 1 | Table 1 | Supplementary Figure S1 |

## 図表キャプションの書き方
**Figure キャプション**: 一文目は図全体の説明（太字）。
二文目以降でパネルごとの説明 (a) ..., (b) ..., (c) ...

**Table キャプション**: 表の上に配置。内容を一文で要約。
略語は表の下部に脚注として記載。
```

### 5. 引用・参考文献の書式

```markdown
## 引用スタイル一覧

### Nature 系（番号順）
本文: "... has been reported¹."
文献: 1. Author, A. B., Author, C. D. & Author, E. F.
       Title of article. *Journal* **vol**, pages (year).

### Science 系（番号順）
本文: "... has been reported (1)."
文献: 1. A. B. Author, C. D. Author, E. F. Author,
       Title of article. *Journal* **vol**, pages (year).

### ACS 系（番号順、上付き）
本文: "... has been reported.¹"
文献: (1) Author, A. B.; Author, C. D.; Author, E. F.
       Title of Article. *Journal* **year**, *vol*, pages.

### IEEE 系（角括弧番号）
本文: "... has been reported [1]."
文献: [1] A. B. Author, C. D. Author, and E. F. Author,
       "Title of article," *Journal*, vol. X, no. Y, pp. XX-YY, year.

### Elsevier 系（著者名-年）
本文: "... has been reported (Author et al., 2024)."
文献: Author, A.B., Author, C.D., Author, E.F., 2024.
       Title of article. Journal Vol, pages.
```

### 6. Cover Letter テンプレート

```markdown
Dear Editor,

We are pleased to submit our manuscript entitled "[タイトル]" for
consideration for publication in [ジャーナル名].

[1-2 文で研究背景と動機を記述]

In this study, we [研究の主要なアプローチを記述].
Our key findings include:
(1) [主要な発見 1],
(2) [主要な発見 2], and
(3) [主要な発見 3].

These results [研究の意義・インパクトを記述], which we believe will be of
significant interest to the readership of [ジャーナル名].

This manuscript has not been published or submitted elsewhere.
All authors have approved the manuscript and agree with its submission
to [ジャーナル名]. We declare no competing interests.

We suggest the following reviewers:
1. Prof. [Name], [Affiliation] ([email])
2. Prof. [Name], [Affiliation] ([email])
3. Prof. [Name], [Affiliation] ([email])

Thank you for your consideration.

Sincerely,
[Corresponding Author Name]
[Affiliation]
[Email]
```

### 7. Response to Reviewers テンプレート

```markdown
# Response to Reviewers

We thank the reviewers for their constructive comments, which have
significantly improved our manuscript. Below we address each comment
point by point. Reviewer comments are shown in **bold**, and our
responses follow each comment. Changes in the revised manuscript
are highlighted in blue.

---

## Reviewer 1

**Comment 1**: [レビュアーのコメントをそのまま引用]

**Response**: We thank the reviewer for this insightful comment.
[回答の内容]
We have revised the manuscript accordingly (page X, lines YY-ZZ).

> **Revised text**: "[修正箇所の文章を引用]"

---

**Comment 2**: [コメント]

**Response**: [回答]

---

## Reviewer 2

**Comment 1**: [コメント]

**Response**: [回答]
```

### 8. 図の埋め込みワークフロー

論文作成時、解析で生成した図を原稿内に埋め込む。`figures/` ディレクトリに
保存された画像ファイルを自動的に検出し、Markdown の画像構文で本文に挿入する。

```markdown
## 図の埋め込みルール

### 基本手順
1. `figures/` ディレクトリを走査し、生成済みの画像ファイル一覧を取得する
2. 各図を Results セクションの該当箇所に `![Figure N](figures/filename.png)` で埋め込む
3. Figure Captions / Figure Legends セクションにも同じ画像を配置する
4. Supplementary の図は `![Figure SN](figures/filename.png)` で埋め込む

### 画像埋め込み構文

#### 本文中（Results セクション）
- 結果の記述と共に図を直接埋め込む
- キャプションは図の直下に記述する

![Figure 1](figures/fig1_overview.png)
**Figure 1.** [図の説明]

#### マルチパネル図
- composite figure は 1 つの画像として埋め込む
- パネルの説明はキャプション内で (a), (b), (c) で記述する

![Figure 2](figures/fig2_composite.png)
**Figure 2.** [図全体の説明。]
(a) [パネル a]。(b) [パネル b]。(c) [パネル c]。

### ジャーナル別の図参照 + 埋め込み

| ジャーナル | 埋め込み例 |
|---|---|
| Nature 系 | `![Fig. 1](figures/fig1.png)` + `**Fig. 1 \| [タイトル].**` |
| Science 系 | `![Fig. 1](figures/fig1.png)` + `**Fig. 1. [タイトル].**` |
| ACS 系 | `![Figure 1](figures/figure1.png)` + `**Figure 1.** [説明]` |
| IEEE 系 | `![Fig. 1](figures/fig1.png)` + `**Fig. 1.** [説明]` |
| Elsevier 系 | `![Fig. 1](figures/fig1.png)` + `**Fig. 1.** [説明]` |
| Qiita | `![代替テキスト](./figures/Fig1_description.png)` |

### figures/ ディレクトリの走査

論文執筆時、以下の手順で図を収集する:
1. ワークスペース内の `figures/` ディレクトリを検索
2. `.png`, `.svg`, `.pdf` ファイルを一覧取得
3. ファイル名から図番号・内容を推定（例: `violin_hardness.png` → 硬さの Violin Plot）
4. Results の記述順に合わせて Figure 番号を割り当て
5. 本文中の該当箇所に画像を埋め込み、キャプションを生成
```

### 9. AI 使用開示 (AI Usage Disclosure)

多くのジャーナルが生成 AI の使用開示を義務化している。
AI はオーサーシップの要件を満たさないため、**著者として記載しない**。

```markdown
## AI 使用開示の記載場所（ジャーナル別）

| ジャーナル | 記載場所 | ポリシー |
|---|---|---|
| Nature 系 | Methods セクション末尾 | 必須 (2023年〜) |
| Science 系 | Acknowledgements | 必須 |
| ACS 系 | Methods or Acknowledgements | 推奨 |
| IEEE 系 | Acknowledgements | 推奨 |
| Elsevier 系 | 専用 AI Disclosure セクション | 必須 (2024年〜) |
```

#### AI 使用開示テンプレート

**注意**: バージョン番号は package.json から動的に取得すること。

```markdown
## Methods セクション末尾（Nature 系）

**Use of AI tools**: This study used GitHub Copilot Agent with SATORI
skills (v{package.json の version}) for [具体的な用途: data analysis /
figure generation / manuscript drafting / literature review].
All AI-generated content was reviewed and verified by the authors,
who take full responsibility for the content of this publication.

## Acknowledgements（Science / ACS / IEEE 系）

The authors acknowledge the use of GitHub Copilot Agent with SATORI
skills (v{package.json の version}) for [具体的な用途].
All outputs were critically reviewed and validated by the authors.

## AI Disclosure セクション（Elsevier 系）

During the preparation of this work, the authors used GitHub Copilot
Agent with SATORI skills (v{package.json の version}) for [具体的な用途].
After using this tool, the authors reviewed and edited the content as
needed and take full responsibility for the content of the publication.
```

#### Authors セクションのルール

```markdown
## AI をオーサーに含める場合の注意

❌ 禁止: AI ツールをオーサーとして記載してはならない
   (ICMJE ガイドライン: オーサーシップには説明責任が必要)

❌ 誤り:
   Author A¹, GitHub Copilot Agent (SATORI v0.2.0)¹, Author B²

✅ 正しい:
   Author A¹*, Author B²
   （AI 使用は Methods / Acknowledgements / AI Disclosure に記載）
```

### 10. Supplementary Information 構成

```markdown
## Supplementary Information テンプレート

# Supplementary Information for:
# [論文タイトル]

[著者名]

## Supplementary Figures

**Figure S1.** [キャプション]

**Figure S2.** [キャプション]

## Supplementary Tables

**Table S1.** [キャプション]

## Supplementary Methods

### S1. [追加実験手法の詳細]

### S2. [追加解析手法の詳細]

## Supplementary References

[SI 内でのみ引用した文献]
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `crossref` | Crossref | 論文メタデータ検索・引用情報取得 |

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `manuscript/manuscript.md` | Markdown 原稿 |
| `manuscript/figures/` | 図表ディレクトリ |
| `manuscript/supplementary.md` | 補足情報 |
| `manuscript/cover_letter.md` | カバーレター |
| `manuscript/response_to_reviewers.md` | レビュー回答 |

### テンプレートファイル (assets/)

| ファイル | 対応ジャーナル |
|---|---|
| `assets/imrad_standard.md` | IMRaD 標準形式 |
| `assets/nature_article.md` | Nature / Nature Communications |
| `assets/science_research_article.md` | Science / Science Advances |
| `assets/acs_article.md` | ACS Nano / JACS / Chem. Mater. |
| `assets/ieee_transactions.md` | IEEE Transactions |
| `assets/elsevier_article.md` | Elsevier 系ジャーナル |
| `assets/qiita_technical_article.md` | Qiita 技術記事（AI for Science シリーズ） |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-hypothesis-pipeline` | `docs/hypothesis.md` を読み込み、Introduction の研究仮説・Methods の解析計画に自動反映 |
| `scientific-critical-review` | 草稿完成後のセルフレビュー。投稿前に必ず実行を推奨 |
| `scientific-publication-figures` | `figures/` の図表を原稿に埋め込み |
| `scientific-pipeline-scaffold` | `results/analysis_summary.json` から数値・統計結果を参照 |
| `scientific-statistical-testing` | 統計的主張の正確な記述に使用 |

### パイプライン上の位置づけ

```
hypothesis-pipeline → pipeline-scaffold → academic-writing → critical-review
  (仮説定義)         (解析実行)         (草稿作成)         (レビュー・修正)
```

論文執筆時には以下のファイルを参照する:
- `docs/hypothesis.md` — 仮説定義（Introduction / Methods に反映）
- `docs/workflow_design.md` — 解析ワークフロー（Methods に反映）
- `results/analysis_summary.json` — 解析結果サマリー（Results に反映）
- `figures/*.png` — 図表（本文中に `![Figure N](figures/...)` で埋め込み）
