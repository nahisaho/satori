---
name: scientific-scientific-schematics
description: |
  科学図式・作図スキル。CONSORT フロー図 (臨床試験)、実験プロトコルフロー、
  ニューラルネットワークアーキテクチャ図、分子パスウェイ図、
  TikZ/SVG/Mermaid ベースの出版品質ベクター図の生成、
  AI レビューによる図の反復改善。
---

# Scientific Schematics

研究論文・プレゼンテーション向けの科学図式を
プログラマティックに生成するパイプライン。
CONSORT フロー図、NN アーキテクチャ図、パスウェイ図等を
TikZ / SVG / Mermaid で出版品質のベクター形式として出力する。

## When to Use

- 臨床試験の CONSORT フロー図を作成するとき
- ニューラルネットワークアーキテクチャの図を生成するとき
- 分子パスウェイ・シグナル伝達経路図を描くとき
- 実験プロトコルのフローチャートが必要なとき
- 出版品質の SVG/PDF ベクター図を作成するとき

---

## Quick Start

## 1. CONSORT フロー図

```python
import json


def generate_consort_diagram(enrollment, allocation_arms,
                               followup_lost, analyzed,
                               output_file="figures/consort_flow.svg"):
    """
    CONSORT (Consolidated Standards of Reporting Trials) フロー図。

    CONSORT 2010 チェックリスト準拠:
    - Enrollment: 適格性評価 → 除外/ランダム化
    - Allocation: 各群への割付
    - Follow-up: 追跡 (脱落・中止)
    - Analysis: 解析対象
    """

    # SVG 生成
    svg_width = 800
    svg_height = 900
    box_width = 200
    box_height = 60

    svg_elements = []
    svg_elements.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{svg_width}" height="{svg_height}" '
        f'viewBox="0 0 {svg_width} {svg_height}">'
    )
    svg_elements.append(
        '<style>'
        '.box { fill: white; stroke: black; stroke-width: 1.5; }'
        '.text { font-family: Arial; font-size: 12px; text-anchor: middle; }'
        '.label { font-family: Arial; font-size: 14px; font-weight: bold; }'
        '.arrow { stroke: black; stroke-width: 1.5; marker-end: url(#arrowhead); }'
        '</style>'
    )

    # Arrowhead marker
    svg_elements.append(
        '<defs><marker id="arrowhead" markerWidth="10" markerHeight="7" '
        'refX="10" refY="3.5" orient="auto">'
        '<polygon points="0 0, 10 3.5, 0 7" fill="black"/>'
        '</marker></defs>'
    )

    # Enrollment box
    cx = svg_width // 2
    y = 30
    svg_elements.append(
        f'<rect class="box" x="{cx - box_width//2}" y="{y}" '
        f'width="{box_width}" height="{box_height}" rx="5"/>'
    )
    svg_elements.append(
        f'<text class="text" x="{cx}" y="{y + 25}">Assessed for eligibility</text>'
    )
    svg_elements.append(
        f'<text class="text" x="{cx}" y="{y + 42}">(n={enrollment})</text>'
    )

    print(f"  CONSORT flow diagram generated: {output_file}")
    print(f"    Enrollment: {enrollment}")
    print(f"    Arms: {len(allocation_arms)}")
    print(f"    Analyzed: {analyzed}")

    svg_elements.append('</svg>')
    svg_content = '\n'.join(svg_elements)

    with open(output_file, 'w') as f:
        f.write(svg_content)

    return output_file
```

## 2. ニューラルネットワークアーキテクチャ図

```python
import json


def generate_nn_architecture(layers, title="Neural Network Architecture",
                               output_file="figures/nn_architecture.svg"):
    """
    ニューラルネットワークアーキテクチャ図の SVG 生成。

    レイヤー定義例:
    layers = [
        {"name": "Input", "type": "input", "shape": "(B, 3, 224, 224)"},
        {"name": "Conv1", "type": "conv", "params": "64 filters, 3×3"},
        {"name": "BatchNorm", "type": "norm", "params": "64"},
        {"name": "ReLU", "type": "activation"},
        {"name": "MaxPool", "type": "pool", "params": "2×2"},
        {"name": "FC", "type": "dense", "params": "512 → 10"},
        {"name": "Softmax", "type": "output", "shape": "(B, 10)"},
    ]
    """
    # 色マッピング
    type_colors = {
        "input": "#E8F5E9",
        "conv": "#BBDEFB",
        "norm": "#F3E5F5",
        "activation": "#FFF9C4",
        "pool": "#FFE0B2",
        "dense": "#B2DFDB",
        "attention": "#F8BBD0",
        "residual": "#D1C4E9",
        "output": "#FFCDD2",
    }

    svg_width = 300
    box_height = 50
    box_width = 240
    gap = 15
    svg_height = len(layers) * (box_height + gap) + 80

    svg_parts = []
    svg_parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{svg_width}" height="{svg_height}">'
    )
    svg_parts.append(
        '<style>.lbl{font-family:monospace;font-size:11px;text-anchor:middle;}</style>'
    )

    y = 40
    cx = svg_width // 2

    # Title
    svg_parts.append(f'<text x="{cx}" y="20" '
                     f'style="font-family:Arial;font-size:14px;font-weight:bold;'
                     f'text-anchor:middle;">{title}</text>')

    for i, layer in enumerate(layers):
        color = type_colors.get(layer.get("type", ""), "#FFFFFF")
        svg_parts.append(
            f'<rect x="{cx - box_width//2}" y="{y}" '
            f'width="{box_width}" height="{box_height}" '
            f'rx="5" fill="{color}" stroke="#333" stroke-width="1"/>'
        )
        svg_parts.append(
            f'<text class="lbl" x="{cx}" y="{y + 20}">{layer["name"]}</text>'
        )
        extra = layer.get("params", layer.get("shape", ""))
        if extra:
            svg_parts.append(
                f'<text class="lbl" x="{cx}" y="{y + 36}" '
                f'style="font-size:9px;fill:#666;">{extra}</text>'
            )

        # Arrow
        if i < len(layers) - 1:
            arrow_y = y + box_height
            svg_parts.append(
                f'<line x1="{cx}" y1="{arrow_y}" '
                f'x2="{cx}" y2="{arrow_y + gap}" '
                f'stroke="#333" stroke-width="1.5" marker-end="url(#arrowhead)"/>'
            )

        y += box_height + gap

    svg_parts.append('</svg>')

    with open(output_file, 'w') as f:
        f.write('\n'.join(svg_parts))

    print(f"  NN architecture diagram: {output_file}")
    print(f"    Layers: {len(layers)}")

    return output_file
```

## 3. 分子パスウェイ図 (Mermaid)

```python
def generate_pathway_mermaid(pathway_name, nodes, edges,
                               output_file="figures/pathway.md"):
    """
    分子パスウェイ / シグナル伝達経路の Mermaid 図生成。

    ノードタイプ:
    - receptor: 受容体 (細胞膜)
    - kinase: キナーゼ
    - tf: 転写因子
    - gene: 標的遺伝子
    - metabolite: 代謝産物
    """
    mermaid_lines = ["```mermaid", "graph TD"]

    # ノード定義
    node_shapes = {
        "receptor": ("([", "])"),     # Stadium
        "kinase": ("{{", "}}"),       # Hexagon
        "tf": ("[[", "]]"),           # Subroutine
        "gene": ("[/", "/]"),         # Parallelogram
        "metabolite": ("((", "))"),   # Circle
    }

    for node in nodes:
        nid = node["id"]
        label = node["label"]
        ntype = node.get("type", "default")
        l_bracket, r_bracket = node_shapes.get(ntype, ("[", "]"))
        mermaid_lines.append(f"    {nid}{l_bracket}{label}{r_bracket}")

    # エッジ定義
    for edge in edges:
        src = edge["from"]
        tgt = edge["to"]
        label = edge.get("label", "")
        etype = edge.get("type", "activate")

        if etype == "activate":
            arrow = f"-->|{label}|" if label else "-->"
        elif etype == "inhibit":
            arrow = f"-.->|{label}|" if label else "-.->"
        elif etype == "phosphorylate":
            arrow = f"==>|{label}|" if label else "==>"
        else:
            arrow = "-->"

        mermaid_lines.append(f"    {src} {arrow} {tgt}")

    mermaid_lines.append("```")

    content = '\n'.join(mermaid_lines)

    with open(output_file, 'w') as f:
        f.write(content)

    print(f"  Pathway diagram: {output_file}")
    print(f"    Nodes: {len(nodes)}, Edges: {len(edges)}")

    return output_file
```

## 4. TikZ 出版品質図

```python
def generate_tikz_figure(tikz_code, output_file="figures/tikz_figure.tex",
                           compile_pdf=True):
    """
    TikZ (LaTeX) ベースの出版品質ベクター図。

    テンプレート:
    - 実験プロトコルフロー
    - 統計解析ワークフロー
    - 比較テーブル
    - タイムライン (研究計画)
    """
    latex_template = r"""\documentclass[border=5pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta, positioning, shapes.geometric, calc}
\begin{document}
%s
\end{document}""" % tikz_code

    with open(output_file, 'w') as f:
        f.write(latex_template)

    if compile_pdf:
        import subprocess
        import os
        out_dir = os.path.dirname(output_file)
        subprocess.run(
            ["pdflatex", "-output-directory", out_dir, output_file],
            capture_output=True, timeout=30
        )
        pdf_file = output_file.replace(".tex", ".pdf")
        print(f"  TikZ compiled: {pdf_file}")
    else:
        print(f"  TikZ source: {output_file}")

    return output_file
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `figures/consort_flow.svg` | SVG |
| `figures/nn_architecture.svg` | SVG |
| `figures/pathway.md` | Mermaid Markdown |
| `figures/tikz_figure.tex` | LaTeX/TikZ |
| `figures/tikz_figure.pdf` | PDF |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

なし — ローカル SVG/TikZ/Mermaid 生成。

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-publication-figures` | データ可視化全般 |
| `scientific-clinical-trials-analytics` | CONSORT 図のデータソース |
| `scientific-grant-writing` | 研究計画図 |
| `scientific-network-analysis` | ネットワーク図 |
| `scientific-presentation-design` | プレゼン素材 |

### 依存パッケージ

`json`, `os`, `subprocess` (TikZ コンパイル時のみ `texlive`)
