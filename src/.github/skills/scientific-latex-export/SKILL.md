---
name: scientific-latex-export
description: |
  実験結果を論文形式（LaTeX / IMRaD）にエクスポートするスキル。
  Introduction・Materials & Methods・Results・Discussion の構造で
  出版準備用の原稿を自動生成する。
  「論文にして」「LaTeX出力」「出版準備」で発火。
---

# Scientific LaTeX Export

実験結果を学術論文形式（IMRaD: Introduction, Materials & Methods,
Results, Discussion）に構造化し、LaTeX ソースを生成するスキル。

## When to Use

- 実験結果を論文形式にまとめるとき
- 学会発表用のアブストラクトを作成するとき
- 共同研究者と出版用原稿を共有するとき
- arXiv / bioRxiv へのプレプリント投稿を準備するとき

## Quick Start

### Step 1: 実験データの収集
- 実験のメッセージ履歴を読解
- 生成された成果物（CSV, 画像, レポート）を収集
- 使用したスキル・ツールの記録

### Step 2: IMRaD 構造への変換
1. **Title & Abstract** — 研究の要約（250語以内）
2. **Introduction** — 背景・目的・仮説
3. **Materials & Methods** — 実験条件・手法・解析手法
4. **Results** — データ・図表の記述
5. **Discussion** — 結果の解釈・先行研究との比較・限界
6. **References** — 引用文献（BibTeX形式）

### Step 3: 図表の整形
- matplotlib / R の図を publication-quality に調整
- 表を LaTeX tabular 形式に変換
- キャプションの作成

### Step 4: LaTeX ソース生成

## Output Format

```markdown
paper.tex       # LaTeX メインファイル
paper.bib       # BibTeX 参考文献
figures/        # 図ファイル
```

### LaTeX テンプレート構造:
```latex
\documentclass[12pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{graphicx, amsmath, booktabs, hyperref}
\usepackage[japanese]{babel}  % 日本語対応

\title{...}
\author{...}
\date{\today}

\begin{document}
\maketitle
\begin{abstract} ... \end{abstract}

\section{Introduction}
\section{Materials and Methods}
\section{Results}
\section{Discussion}
\section*{Acknowledgments}

\bibliographystyle{unsrt}
\bibliography{paper}
\end{document}
```

## Examples

- 「この実験結果を論文にして」
- 「IMRaD形式でLaTeX出力して」
- 「bioRxiv投稿用のプレプリントを準備して」

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `crossref` | Crossref | 参考文献メタデータ・DOI 解決 |
