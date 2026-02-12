---
name: scientific-presentation-design
description: |
  科学プレゼンテーション・ポスター・模式図設計スキル。学会発表スライド、
  LaTeX/PPTX ポスター、科学模式図（ワークフロー図・メカニズム図）、
  ビジュアルアブストラクトの作成を支援。claude-scientific-skills の
  Scientific Communication カテゴリ（slides, posters, schematics）を統合。
  「学会スライドを作成して」「ポスターのレイアウトを設計して」で発火。
---

# Scientific Presentation Design

科学プレゼンテーションの設計・制作スキル。学会発表スライド、
ポスター（LaTeX/PPTX）、科学模式図（pathway/mechanism 図）、
ビジュアルアブストラクトの作成を支援する。

## When to Use

- 学会発表用のスライドデッキを作成するとき
- 学会ポスターのレイアウトを設計するとき
- 論文の Graphical Abstract を作成するとき
- メカニズム図・ワークフロー図を描くとき
- 研究成果のビジュアルサマリーを作成するとき

## Quick Start

### プレゼンテーション設計パイプライン

```
Input: 研究内容 / 対象学会 / 発表形式
    ↓
Step 1: Structure Design
  - ストーリーライン構築
  - セクション構成
  - 時間配分（oral の場合）
    ↓
Step 2: Visual Design
  - カラーパレット選択
  - フォント設定
  - レイアウトグリッド
    ↓
Step 3: Content Creation
  - テキストブロック
  - 図表配置
  - 模式図作成
    ↓
Step 4: Polish
  - 整合性チェック
  - アクセシビリティ確認
  - エクスポート
    ↓
Output: Slides / Poster / Schematic
```

---

## Phase 1: 学会スライド設計

### 標準構成テンプレート

```markdown
## Oral Presentation Structure (15 min)

### Slide 1: Title (0:00 - 0:30)
- タイトル / 著者 / 所属
- 1 sentence hook

### Slides 2-3: Background (0:30 - 2:30)
- 研究の文脈・先行研究
- Knowledge gap の明示

### Slide 4: Research Question (2:30 - 3:30)
- 明確な RQ / 仮説
- 本研究の目的

### Slides 5-7: Methods (3:30 - 6:00)
- 実験デザイン
- 解析手法
- ワークフロー図

### Slides 8-11: Results (6:00 - 11:00)
- 主要結果（1 slide = 1 message）
- 図表は大きく、テキストは最小限
- 統計値を明示

### Slide 12: Discussion (11:00 - 13:00)
- 結果の解釈
- 先行研究との比較
- 限界

### Slide 13: Conclusions (13:00 - 14:00)
- Take-home message（3 点以内）
- Future directions

### Slide 14: Acknowledgments (14:00 - 14:30)
- 共同研究者 / Funding

### Slide 15: Q&A / Backup
```

### スライドデザインルール

```python
SLIDE_DESIGN_RULES = {
    "fonts": {
        "title": {"family": "Arial/Helvetica", "size": 28, "weight": "bold"},
        "body": {"family": "Arial/Helvetica", "size": 20, "weight": "normal"},
        "caption": {"family": "Arial/Helvetica", "size": 16, "weight": "normal"},
    },
    "colors": {
        "scientific_blue": ["#1B365D", "#2E5090", "#4A7FB5", "#8CB4D6", "#C5DAE9"],
        "nature_green": ["#1B4332", "#2D6A4F", "#40916C", "#74C69D", "#B7E4C7"],
        "warm_red": ["#6A040F", "#9D0208", "#D00000", "#DC2F02", "#E85D04"],
    },
    "layout": {
        "margin": "0.5 inch all sides",
        "max_bullets": 5,
        "max_words_per_slide": 40,
        "figure_ratio": 0.6,  # 図の面積が 60% 以上
    },
    "accessibility": {
        "min_contrast_ratio": 4.5,
        "colorblind_safe_palette": True,
        "alt_text_for_figures": True,
    },
}
```

---

## Phase 2: ポスター設計

### LaTeX ポスターテンプレート

```latex
\\documentclass[25pt, a0paper, portrait]{tikzposter}
\\usepackage[utf8]{inputenc}

\\title{\\parbox{\\linewidth}{\\centering
  Research Title Goes Here: \\\\
  A Comprehensive Study
}}
\\author{Author A$^1$, Author B$^{1,2}$}
\\institute{$^1$Institution, $^2$Institution}

\\usetheme{Default}
\\usecolorstyle{Britain}

\\begin{document}
\\maketitle

\\begin{columns}
\\column{0.33}
\\block{Background}{
  % 背景・先行研究
}
\\block{Methods}{
  % 手法（ワークフロー図推奨）
}

\\column{0.34}
\\block{Results}{
  % 主要結果・図表
}

\\column{0.33}
\\block{Discussion}{
  % 考察
}
\\block{Conclusions}{
  % 結論（3 bullet points）
}
\\block{References}{
  % 引用文献（最小限）
}
\\end{columns}

\\end{document}
```

### ポスターレイアウト原則

```markdown
## Poster Design Principles

### Layout
- 3-column or 2-column layout
- 読み順: 上→下、左→右
- セクション間の明確な区切り
- ヘッダー（タイトル + 著者）は上部 15%

### Typography
- Title: 72-96 pt
- Section headers: 48-60 pt
- Body text: 28-36 pt
- 3 メートル離れて読めること

### Figures
- 高解像度 (300+ DPI)
- 大きく（ポスター面積の 40-50%）
- 自己説明的なキャプション
- カラーバリアフリー

### Content
- Introduction: 150-200 words
- Methods: 図で示す（テキスト最小限）
- Results: 図中心
- Conclusions: 3-5 bullet points
- References: 5-10 citations max
```

---

## Phase 3: 科学模式図

### Mermaid ベースの模式図

```markdown
## 模式図タイプ

### 1. ワークフロー図
ステップの流れ・パイプラインを視覚化

### 2. メカニズム図
分子・シグナル経路のメカニズムを図示

### 3. 実験デザイン図
実験群の配置・処理条件を視覚化

### 4. Graphical Abstract
論文全体を 1 枚に要約
```

```python
# matplotlib ベースの模式図生成
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def create_workflow_schematic(steps, title="Workflow"):
    """
    ワークフロー模式図の自動生成。
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, len(steps) * 1.5 + 2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(steps) * 2 + 1)
    ax.axis('off')

    colors = plt.cm.Blues(np.linspace(0.3, 0.8, len(steps)))

    for i, (step_name, description) in enumerate(steps):
        y = len(steps) * 2 - i * 2
        # Box
        rect = mpatches.FancyBboxPatch(
            (1, y - 0.5), 8, 1.2,
            boxstyle="round,pad=0.1",
            facecolor=colors[i], edgecolor='#333333', linewidth=1.5
        )
        ax.add_patch(rect)
        ax.text(5, y + 0.1, step_name, ha='center', va='center',
                fontsize=14, fontweight='bold', color='white')
        ax.text(5, y - 0.25, description, ha='center', va='center',
                fontsize=10, color='white', style='italic')
        # Arrow
        if i < len(steps) - 1:
            ax.annotate('', xy=(5, y - 0.7), xytext=(5, y - 0.5),
                        arrowprops=dict(arrowstyle='->', color='#333333', lw=2))

    ax.set_title(title, fontsize=18, fontweight='bold', pad=20)
    plt.tight_layout()
    return fig
```

---

## Report Template

```markdown
# Presentation Design: [Title]

**Format**: [Oral / Poster / Graphical Abstract]
**Venue**: [Conference / Journal]
**Date**: [date]

## 1. Story Arc
## 2. Design Specifications
| Element | Value |
|---------|-------|
| Format | |
| Dimensions | |
| Color palette | |
| Font family | |

## 3. Section Outline
## 4. Figures Required
## 5. Generated Files
- [ ] slides.pptx / slides.tex
- [ ] poster.pdf / poster.tex
- [ ] figures/*.svg
```

---

## Completeness Checklist

- [ ] ストーリーライン: 明確な論理の流れ
- [ ] デザイン: カラーパレット選択、フォント統一
- [ ] 図表: 高解像度、カラーバリアフリー
- [ ] テキスト: 最小限、1 slide = 1 message
- [ ] アクセシビリティ: コントラスト比 ≥ 4.5
- [ ] エクスポート: PDF + ソースファイル

## Best Practices

1. **1 Slide = 1 Message**: 情報を詰め込みすぎない
2. **図 > テキスト**: 視覚的に伝える
3. **カラーバリアフリー**: viridis/cividis 系パレットを推奨
4. **3 メートルルール**: ポスターは遠距離から読めること
5. **10-20-30 ルール**: 10 slides, 20 min, 30pt font (Kawasaki)
