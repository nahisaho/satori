---
name: scientific-grant-writing
description: |
  研究グラント（助成金申請書）執筆スキル。NIH R01/R21、NSF、JSPS 科研費、
  ERC 等のフォーマットに対応。Specific Aims、Research Strategy、予算計画、
  Biosketch の構造化作成を支援。
  「グラント申請書を書いて」「科研費を作成して」「Specific Aims を書いて」で発火。
---

# Scientific Grant Writing

研究助成金申請書の執筆支援スキル。主要ファンディング機関
（NIH、NSF、JSPS、ERC）のフォーマットに準拠した申請書の
構造化・ドラフト作成を支援する。

## When to Use

- 研究費申請書のドラフトを作成するとき
- Specific Aims ページを構成するとき
- Research Strategy を構造化するとき
- 予算計画（Budget Justification）を作成するとき
- ファンディング機関のフォーマットに適合させるとき

## Quick Start

### グラント執筆パイプライン

```
Step 1: Funder Selection & Format
  - ファンディング機関の選択
  - 申請カテゴリ（NIH R01, R21, JSPS 基盤A, etc.）
  - ページ/文字制限の確認
    ↓
Step 2: Specific Aims
  - Opening paragraph (知識のギャップ)
  - Aims (2-3 aims, 具体的・測定可能)
  - Impact statement
    ↓
Step 3: Research Strategy
  - Significance
  - Innovation
  - Approach (per aim)
    ↓
Step 4: Supporting Documents
  - Budget & Justification
  - Biosketch / 研究業績リスト
  - Facilities & Equipment
  - Letters of Support
    ↓
Output: Complete Grant Application
```

---

## Phase 1: Specific Aims ページ

### NIH Specific Aims テンプレート

```markdown
## Specific Aims (1 page)

### Opening Paragraph
[問題の重要性 → 現在の知識 → 知識のギャップ → なぜそのギャップが重要か]

**句型テンプレート:**
"[Disease/Problem] affects [X million] people annually, resulting in [consequence].
Current approaches [limitation]. Our preliminary data show [finding],
suggesting that [hypothesis]. However, [gap]. This gap is critical because [reason]."

### Long-term Goal & Objective
"The long-term goal of this project is to [broad vision].
The objective of this application is to [specific, achievable goal].
Our central hypothesis is that [testable hypothesis],
based on [preliminary data / rationale]."

### Specific Aims

**Aim 1: [Action verb] + [measurable outcome]**
[2-3 sentences describing the aim, approach, and expected outcome]
Working hypothesis: [specific to this aim]

**Aim 2: [Action verb] + [measurable outcome]**
[2-3 sentences describing the aim, approach, and expected outcome]
Working hypothesis: [specific to this aim]

**Aim 3 (optional): [Action verb] + [measurable outcome]**
[2-3 sentences]

### Impact Statement
"Upon completion, this project will [expected outcome].
This contribution is significant because [reason].
The proposed research is innovative because [novel aspect]."
```

---

## Phase 2: Research Strategy

### Significance / Innovation / Approach

```markdown
## Research Strategy (12 pages for NIH R01)

### A. Significance (2-3 pages)
#### A.1 Importance of the Problem
#### A.2 Scientific Premise
  - Rigor of prior research (strengths & weaknesses)
  - How this proposal addresses gaps
#### A.3 Impact if Successful

### B. Innovation (1-2 pages)
#### B.1 Conceptual Innovation
  - New paradigm / model / theory
#### B.2 Technical Innovation
  - New methods / tools / approaches
#### B.3 Application Innovation
  - New applications / fields

### C. Approach (7-8 pages)
#### C.1 Overview & Rationale
#### C.2 Preliminary Data
  - Figures with quality data
  - Demonstrates feasibility

#### C.3 Aim 1: [Title]
  - Rationale
  - Experimental Design
    - Variables (independent, dependent, controlled)
    - Sample size justification (power analysis)
    - Controls (positive, negative)
  - Methods
  - Expected Outcomes
  - Potential Problems & Alternative Strategies
  - Timeline

#### C.4 Aim 2: [Title]
  [同様の構成]

#### C.5 Rigor & Reproducibility
  - Biological variables (sex, age, strain)
  - Statistical approach
  - Data management plan
  - Authentication of key resources

#### C.6 Timeline
| Activity | Year 1 | Year 2 | Year 3 | Year 4 | Year 5 |
|----------|--------|--------|--------|--------|--------|
| Aim 1    | ████   | ██     |        |        |        |
| Aim 2    |        | ██     | ████   | ██     |        |
| Aim 3    |        |        |        | ██     | ████   |
```

---

## Phase 3: JSPS 科研費テンプレート

### 研究目的・研究計画・研究業績

```markdown
## JSPS 科研費 様式

### 1. 研究目的、研究方法など（概要）
[400字以内: 研究の全体像を簡潔に記述]

### 2. 本研究の学術的背景
- 学術的「問い」は何か
- 国内外の研究動向と本研究の位置づけ
- これまでの研究成果

### 3. 研究の目的と手法
[研究目的 → 具体的な方法 → 期待される成果]

### 4. 研究の特色・独創性
- なぜこの方法を選択するか
- 先行研究との差別化
- 技術的新規性

### 5. 本研究で何をどこまで明らかにしようとするのか
[具体的達成目標]

### 6. 研究計画
| 年度 | 研究内容 | 予算（千円）|
|------|----------|-----------|

### 7. 研究業績
[ORCID / Google Scholar から自動取得可能]
```

---

## Phase 4: Budget Planning

### 予算計画テンプレート

```python
def create_budget_plan(total_budget, duration_years, personnel_pct=0.60):
    """
    研究費の予算計画テンプレート生成。
    NIH モジュラー予算に対応。
    """
    budget = {
        "total": total_budget,
        "duration": duration_years,
        "annual": total_budget / duration_years,
        "categories": {
            "personnel": total_budget * personnel_pct,
            "equipment": total_budget * 0.10,
            "supplies": total_budget * 0.15,
            "travel": total_budget * 0.05,
            "other": total_budget * 0.05,
            "indirect": total_budget * 0.05,
        }
    }
    return budget

BUDGET_JUSTIFICATION_TEMPLATE = """
## Budget Justification

### A. Senior/Key Personnel
- PI ([Name]): X calendar months effort @ [salary]
  Justification: [role in project]

### B. Other Personnel
- Postdoctoral researcher: 12 calendar months @ [salary]
  Justification: [specific tasks]

### C. Equipment (>$5,000)
- [Equipment name]: $[cost]
  Justification: [why needed, no alternative]

### D. Travel
- Domestic conference: $[cost] x [#]
  Justification: [dissemination]
- International collaboration: $[cost]

### E. Supplies
- Reagents: $[cost]
- Computing: $[cost]

### F. Other Direct Costs
"""
```

---

## Report Template

```markdown
# Grant Application: [Title]

**Funder**: [NIH / NSF / JSPS / ERC]
**Category**: [R01 / 基盤B / etc.]
**PI**: [Name]
**Date**: [date]

## Document Checklist
- [ ] Specific Aims (1 page)
- [ ] Research Strategy (12 pages)
- [ ] Budget & Justification
- [ ] Biosketch
- [ ] Facilities & Equipment
- [ ] Data Management Plan
- [ ] Letters of Support
```

---

## Completeness Checklist

- [ ] Specific Aims: ギャップ → 仮説 → Aims → インパクト
- [ ] Research Strategy: Significance + Innovation + Approach
- [ ] Preliminary Data: ≥2 figures showing feasibility
- [ ] Power Analysis: サンプルサイズの根拠
- [ ] Alternative Strategy: 各 Aim に Plan B を記載
- [ ] Budget: カテゴリ別の合理的配分
- [ ] Timeline: 全 Aim の Gantt チャート

## Best Practices

1. **Reviewer の視点で書く**: 専門外の審査員にも伝わる明快さ
2. **Specific Aims は独立させる**: 1 Aim の失敗が全体に波及しない設計
3. **Preliminary Data を重視**: Feasibility の証拠が採択率を大きく左右
4. **ページ制限の 95% を使う**: 短すぎるのも減点対象
5. **Broader Impacts を忘れない**: 社会的インパクトの記述（NSF）

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `grants/specific_aims.md` | Specific Aims ページ（Markdown） | Aims 策定完了時 |
| `grants/research_strategy.md` | Research Strategy（Markdown） | 全セクション完了時 |
| `grants/budget.json` | 予算計画（JSON） | 予算見積完了時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-research-methodology` | ← 研究デザイン・方法論の提供 |
| `scientific-hypothesis-pipeline` | ← 仮説定義・課題設定 |
| `scientific-deep-research` | ← 文献レビュー・先行研究調査 |
| `scientific-academic-writing` | ← アカデミックライティングスキル |
