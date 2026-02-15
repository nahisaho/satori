---
name: scientific-latex-formatter
description: |
  Markdown 原稿を LaTeX 形式に変換し、ジャーナル指定のテンプレート（.cls/.sty）に
  適合するフォーマッティングを行うスキル。数式・図表・引用・相互参照の LaTeX 構文変換、
  ジャーナル別スタイル適用、コンパイル可能な .tex ファイル生成を担う。
  「LaTeX に変換して」「投稿用 TeX を作って」「ジャーナルフォーマットにして」で発火。
tu_tools:
  - key: crossref
    name: Crossref
    description: 参考文献メタデータ・DOI 解決
---

# Scientific LaTeX Formatter

Markdown 原稿を LaTeX に変換し、ジャーナル指定テンプレートに適合する
投稿用 .tex ファイルを生成するスキル。

## When to Use

- Markdown で執筆した原稿を LaTeX 形式に変換するとき
- ジャーナル指定の LaTeX テンプレート（revtex4, elsarticle 等）に合わせるとき
- 数式・図表・参照の LaTeX 構文を正確に生成したいとき
- 投稿用の .tex + .bib + figures/ 一式を準備するとき
- Supplementary Information も LaTeX で整形するとき

## Quick Start

## 1. 変換ワークフロー

```
Markdown 原稿 (manuscript/manuscript.md)
  ├─ Phase 1: Markdown → LaTeX 基本変換
  │   ├─ セクション見出し (#, ##, ###) → \section, \subsection
  │   ├─ 太字・斜体 → \textbf, \textit
  │   ├─ リスト → \begin{itemize/enumerate}
  │   └─ コードブロック → \begin{lstlisting} or verbatim
  ├─ Phase 2: 学術要素の変換
  │   ├─ 図の埋め込み → \begin{figure}...\includegraphics
  │   ├─ 表 → \begin{table}...\begin{tabular}
  │   ├─ 数式 ($...$, $$...$$) → \( \), \[ \], equation 環境
  │   ├─ 引用 [1], [2-5] → \cite{key1}, \cite{key1,key2,...}
  │   └─ 相互参照 (Figure 1, Table 2) → \ref{fig:1}, \ref{tab:2}
  ├─ Phase 3: ジャーナルテンプレート適用
  │   ├─ documentclass の設定
  │   ├─ タイトル・著者・所属の構造化
  │   ├─ Abstract 環境の適用
  │   └─ 参考文献スタイル (bibliographystyle) の設定
  └─ Phase 4: ファイル出力
      ├─ manuscript/manuscript.tex
      ├─ manuscript/references.bib
      └─ manuscript/figures/ (コピー)
```

## 2. ジャーナル別 LaTeX テンプレート

```python
JOURNAL_TEMPLATES = {
    "nature": {
        "documentclass": r"\documentclass{nature}",
        "packages": [
            r"\usepackage{graphicx}",
            r"\usepackage{amsmath}",
            r"\usepackage{natbib}",
        ],
        "bib_style": "naturemag",
        "figure_width": r"0.8\textwidth",
        "abstract_env": "abstract",
        "title_cmd": r"\title{%s}",
        "author_cmd": r"\author{%s}",
        "affiliation_cmd": r"\affiliation{%s}",
    },
    "science": {
        "documentclass": r"\documentclass[12pt]{article}",
        "packages": [
            r"\usepackage{graphicx}",
            r"\usepackage{amsmath}",
            r"\usepackage[numbers,sort&compress]{natbib}",
            r"\usepackage{setspace}",
            r"\doublespacing",
        ],
        "bib_style": "Science",
        "figure_width": r"0.9\textwidth",
        "abstract_env": "abstract",
        "title_cmd": r"\title{%s}",
        "author_cmd": r"\author{%s}",
        "affiliation_cmd": None,
    },
    "acs": {
        "documentclass": r"\documentclass[journal=jacsat,manuscript=article]{achemso}",
        "packages": [
            r"\usepackage{graphicx}",
            r"\usepackage{amsmath}",
        ],
        "bib_style": "achemso",
        "figure_width": r"3.25in",
        "abstract_env": "abstract",
        "title_cmd": r"\title{%s}",
        "author_cmd": r"\author{%s}",
        "affiliation_cmd": r"\affiliation{%s}",
    },
    "ieee": {
        "documentclass": r"\documentclass[conference]{IEEEtran}",
        "packages": [
            r"\usepackage{graphicx}",
            r"\usepackage{amsmath}",
            r"\usepackage{cite}",
        ],
        "bib_style": "IEEEtran",
        "figure_width": r"\columnwidth",
        "abstract_env": "abstract",
        "title_cmd": r"\title{%s}",
        "author_cmd": r"\author{%s}",
        "affiliation_cmd": None,
    },
    "elsevier": {
        "documentclass": r"\documentclass[preprint,12pt]{elsarticle}",
        "packages": [
            r"\usepackage{graphicx}",
            r"\usepackage{amsmath}",
            r"\usepackage{lineno}",
            r"\modulolinenumbers[5]",
        ],
        "bib_style": "elsarticle-num",
        "figure_width": r"\textwidth",
        "abstract_env": "abstract",
        "title_cmd": r"\title{%s}",
        "author_cmd": r"\author{%s}",
        "affiliation_cmd": r"\address{%s}",
    },
    "revtex": {
        "documentclass": r"\documentclass[aps,prl,twocolumn,superscriptaddress]{revtex4-2}",
        "packages": [
            r"\usepackage{graphicx}",
            r"\usepackage{amsmath}",
        ],
        "bib_style": "apsrev4-2",
        "figure_width": r"\columnwidth",
        "abstract_env": "abstract",
        "title_cmd": r"\title{%s}",
        "author_cmd": r"\author{%s}",
        "affiliation_cmd": r"\affiliation{%s}",
    },
}
```

## 3. Markdown → LaTeX 変換エンジン

```python
import re
from pathlib import Path


def md_to_latex(md_text, journal_format="elsevier"):
    """
    Markdown テキストを LaTeX に変換する。

    Args:
        md_text: str — Markdown 原稿テキスト
        journal_format: str — ジャーナル形式

    Returns:
        str: LaTeX テキスト
    """
    tex = md_text

    # === セクション見出し ===
    tex = re.sub(r'^### (.+)$', r'\\subsubsection{\1}', tex, flags=re.MULTILINE)
    tex = re.sub(r'^## (.+)$', r'\\subsection{\1}', tex, flags=re.MULTILINE)
    tex = re.sub(r'^# (.+)$', r'\\section{\1}', tex, flags=re.MULTILINE)

    # === インライン書式 ===
    tex = re.sub(r'\*\*\*(.+?)\*\*\*', r'\\textbf{\\textit{\1}}', tex)
    tex = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\1}', tex)
    tex = re.sub(r'\*(.+?)\*', r'\\textit{\1}', tex)

    # === 数式 ===
    # ディスプレイ数式 $$...$$ → \[ \] (equation 環境にラベルが必要な場合は手動)
    tex = re.sub(r'\$\$(.+?)\$\$', r'\\[\1\\]', tex, flags=re.DOTALL)
    # インライン数式は LaTeX でもそのまま使える ($...$)

    # === 図の埋め込み ===
    tex = _convert_figures(tex, journal_format)

    # === 表の変換 ===
    tex = _convert_tables(tex)

    # === 引用の変換 ===
    # [1] → \cite{ref1}, [1, 2] → \cite{ref1,ref2}, [1-3] → \cite{ref1,ref2,ref3}
    tex = re.sub(r'\[(\d+(?:[-,\s]\d+)*)\]', _convert_citation, tex)

    # === 箇条書きの変換 ===
    tex = _convert_lists(tex)

    # === 特殊文字エスケープ ===
    # & と % はすでに LaTeX 構文で使っている場合があるので慎重に
    tex = tex.replace('&amp;', r'\&')
    tex = tex.replace('%', r'\%')

    return tex


def _convert_figures(tex, journal_format):
    """Markdown の図埋め込みを LaTeX figure 環境に変換する。"""
    tmpl = JOURNAL_TEMPLATES.get(journal_format, JOURNAL_TEMPLATES["elsevier"])
    fig_width = tmpl["figure_width"]

    def fig_replacer(m):
        alt_text = m.group(1)
        img_path = m.group(2)
        # ファイル名から拡張子を除去 (LaTeX は拡張子なしでも可)
        img_base = Path(img_path).stem
        # Figure 番号を alt_text から抽出
        fig_num = re.search(r'(\d+)', alt_text)
        label = f"fig:{fig_num.group(1)}" if fig_num else f"fig:{img_base}"

        return (
            f"\\begin{{figure}}[htbp]\n"
            f"  \\centering\n"
            f"  \\includegraphics[width={fig_width}]{{{img_path}}}\n"
            f"  \\caption{{{alt_text}}}\n"
            f"  \\label{{{label}}}\n"
            f"\\end{{figure}}"
        )

    tex = re.sub(r'!\[([^\]]*)\]\(([^)]+)\)', fig_replacer, tex)
    return tex


def _convert_tables(tex):
    """Markdown テーブルを LaTeX tabular 環境に変換する。"""
    lines = tex.split('\n')
    result = []
    in_table = False
    table_lines = []

    for line in lines:
        if re.match(r'^\|.*\|$', line.strip()):
            if not in_table:
                in_table = True
                table_lines = []
            table_lines.append(line.strip())
        else:
            if in_table:
                result.append(_table_to_latex(table_lines))
                in_table = False
                table_lines = []
            result.append(line)

    if in_table:
        result.append(_table_to_latex(table_lines))

    return '\n'.join(result)


def _table_to_latex(table_lines):
    """Markdown テーブル行のリストを LaTeX tabular に変換する。"""
    # ヘッダー行
    headers = [c.strip() for c in table_lines[0].strip('|').split('|')]
    n_cols = len(headers)
    col_spec = '|' + '|'.join(['c'] * n_cols) + '|'

    latex = [
        r"\begin{table}[htbp]",
        r"  \centering",
        f"  \\begin{{tabular}}{{{col_spec}}}",
        r"    \hline",
        "    " + " & ".join(f"\\textbf{{{h}}}" for h in headers) + r" \\",
        r"    \hline",
    ]

    # データ行（セパレータ行 |---|---| をスキップ）
    for line in table_lines[2:]:
        cells = [c.strip() for c in line.strip('|').split('|')]
        latex.append("    " + " & ".join(cells) + r" \\")

    latex.extend([
        r"    \hline",
        r"  \end{tabular}",
        r"  \caption{[テーブルのキャプション]}",
        r"  \label{tab:label}",
        r"\end{table}",
    ])

    return '\n'.join(latex)


def _convert_citation(m):
    """引用番号を \\cite{} 形式に変換する。"""
    raw = m.group(1)
    # "1-3" → "ref1,ref2,ref3", "1, 2" → "ref1,ref2"
    nums = []
    for part in re.split(r'[,\s]+', raw):
        if '-' in part:
            start, end = part.split('-')
            nums.extend(range(int(start), int(end) + 1))
        else:
            nums.append(int(part))
    keys = ','.join(f"ref{n}" for n in sorted(set(nums)))
    return f"\\cite{{{keys}}}"


def _convert_lists(tex):
    """Markdown リストを LaTeX 環境に変換する。"""
    lines = tex.split('\n')
    result = []
    in_list = False
    list_type = None

    for line in lines:
        ul_match = re.match(r'^(\s*)[-*]\s+(.+)$', line)
        ol_match = re.match(r'^(\s*)\d+\.\s+(.+)$', line)

        if ul_match and not in_list:
            in_list = True
            list_type = "itemize"
            result.append(r"\begin{itemize}")
            result.append(f"  \\item {ul_match.group(2)}")
        elif ul_match and in_list and list_type == "itemize":
            result.append(f"  \\item {ul_match.group(2)}")
        elif ol_match and not in_list:
            in_list = True
            list_type = "enumerate"
            result.append(r"\begin{enumerate}")
            result.append(f"  \\item {ol_match.group(2)}")
        elif ol_match and in_list and list_type == "enumerate":
            result.append(f"  \\item {ol_match.group(2)}")
        else:
            if in_list:
                result.append(f"\\end{{{list_type}}}")
                in_list = False
                list_type = None
            result.append(line)

    if in_list:
        result.append(f"\\end{{{list_type}}}")

    return '\n'.join(result)
```

## 4. 完全な .tex ファイル生成

```python
def generate_tex_file(manuscript_path, journal_format="elsevier",
                       title="", authors="", affiliations="",
                       abstract_text="", keywords=None, filepath=None):
    """
    Markdown 原稿から完全な .tex ファイルを生成する。

    Args:
        manuscript_path: Path — 入力 Markdown ファイル
        journal_format: str — ジャーナル形式
        title, authors, affiliations: str — メタデータ
        abstract_text: str — アブストラクト
        keywords: list[str] — キーワード
        filepath: Path — 出力先（デフォルト: manuscript/manuscript.tex）
    """
    from pathlib import Path

    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "manuscript.tex"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    tmpl = JOURNAL_TEMPLATES.get(journal_format, JOURNAL_TEMPLATES["elsevier"])

    # Markdown 本文を読み込み・変換
    with open(manuscript_path, "r", encoding="utf-8") as f:
        md_text = f.read()

    body_tex = md_to_latex(md_text, journal_format)

    # プリアンブル構築
    preamble = [tmpl["documentclass"]]
    preamble.extend(tmpl["packages"])
    preamble.append("")

    # タイトル・著者
    doc_begin = [r"\begin{document}", ""]
    if title:
        doc_begin.append(tmpl["title_cmd"] % title)
    if authors:
        doc_begin.append(tmpl["author_cmd"] % authors)
    if affiliations and tmpl.get("affiliation_cmd"):
        doc_begin.append(tmpl["affiliation_cmd"] % affiliations)
    doc_begin.append(r"\maketitle")

    # アブストラクト
    if abstract_text:
        doc_begin.extend([
            "",
            f"\\begin{{{tmpl['abstract_env']}}}",
            abstract_text,
            f"\\end{{{tmpl['abstract_env']}}}",
        ])

    # キーワード
    if keywords:
        doc_begin.extend([
            "",
            r"\begin{keyword}",
            " \\sep ".join(keywords),
            r"\end{keyword}",
        ])

    # ドキュメント末尾
    doc_end = [
        "",
        f"\\bibliographystyle{{{tmpl['bib_style']}}}",
        r"\bibliography{references}",
        "",
        r"\end{document}",
    ]

    # 結合
    full_tex = '\n'.join(preamble) + '\n' + '\n'.join(doc_begin) + '\n\n'
    full_tex += body_tex + '\n'
    full_tex += '\n'.join(doc_end) + '\n'

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(full_tex)

    print(f"  → LaTeX ファイルを保存: {filepath}")
    return filepath
```

## 5. BibTeX ファイル生成

```python
def generate_bib_file(references, filepath=None):
    """
    引用リストから BibTeX ファイルを生成する。

    Args:
        references: list[dict] — 各参照のメタデータ
            [{"key": "ref1", "type": "article", "author": "...",
              "title": "...", "journal": "...", "year": "...", ...}]
        filepath: Path — 出力先（デフォルト: manuscript/references.bib）
    """
    from pathlib import Path

    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "references.bib"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    content = ""
    for ref in references:
        ref_type = ref.get("type", "article")
        key = ref.get("key", "unknown")
        fields = []

        for field_name in ["author", "title", "journal", "year", "volume",
                            "number", "pages", "doi", "publisher", "booktitle",
                            "url", "note"]:
            if field_name in ref and ref[field_name]:
                fields.append(f"  {field_name} = {{{ref[field_name]}}}")

        content += f"@{ref_type}{{{key},\n"
        content += ",\n".join(fields)
        content += "\n}\n\n"

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"  → BibTeX ファイルを保存: {filepath} ({len(references)} 件)")
    return filepath
```

## 6. LaTeX 変換パイプライン統合

```python
def run_latex_pipeline(manuscript_path, journal_format="elsevier",
                        title="", authors="", affiliations="",
                        abstract_text="", keywords=None, references=None):
    """
    Markdown → LaTeX 変換パイプラインを実行する。

    出力ファイル:
        manuscript/manuscript.tex  — メイン LaTeX ファイル
        manuscript/references.bib  — BibTeX 参照データベース
    """
    print("=" * 60)
    print("LaTeX Formatter Pipeline")
    print("=" * 60)

    # Phase 1-3: Markdown → LaTeX 変換 + テンプレート適用
    print("\n[Phase 1-3] Markdown → LaTeX 変換...")
    tex_path = generate_tex_file(
        manuscript_path, journal_format=journal_format,
        title=title, authors=authors, affiliations=affiliations,
        abstract_text=abstract_text, keywords=keywords,
    )

    # Phase 4: BibTeX 生成
    if references:
        print("\n[Phase 4] BibTeX ファイル生成...")
        generate_bib_file(references)

    print("\n" + "=" * 60)
    print("LaTeX 変換完了！")
    print(f"  コンパイル: pdflatex manuscript.tex && bibtex manuscript && pdflatex manuscript.tex × 2")
    print("=" * 60)

    return tex_path
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `crossref` | Crossref | 参考文献メタデータ・DOI 解決 |

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `manuscript/manuscript.tex` | LaTeX メインファイル | 変換完了時 |
| `manuscript/references.bib` | BibTeX 参照データベース | 変換完了時 |

### 対応ジャーナルテンプレート

| テンプレート | documentclass | bib_style |
|---|---|---|
| Nature 系 | `nature` | `naturemag` |
| Science 系 | `article` (12pt) | `Science` |
| ACS 系 | `achemso` | `achemso` |
| IEEE 系 | `IEEEtran` | `IEEEtran` |
| Elsevier 系 | `elsarticle` | `elsarticle-num` |
| APS/PRL 系 | `revtex4-2` | `apsrev4-2` |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-academic-writing` | 入力: `manuscript/manuscript.md` (Markdown 原稿) |
| `scientific-supplementary-generator` | SI の LaTeX 変換にも対応 |
| `scientific-citation-checker` | 引用キーの整合性を検証後に BibTeX 生成 |
| `scientific-publication-figures` | `figures/` の図パスを `\includegraphics` に変換 |

```
