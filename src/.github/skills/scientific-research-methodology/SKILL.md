---
name: scientific-research-methodology
description: |
  研究方法論・研究デザインスキル。体系的な研究計画策定、ブレインストーミング
  フレームワーク、批判的思考法、研究倫理・IRB、先行研究評価、
  クロスドメイン着想法を含む研究者のメタスキル群。
  「研究計画を立てて」「ブレインストーミングして」「研究デザインを設計して」で発火。
---

# Scientific Research Methodology

研究方法論のメタスキル。体系的な研究デザイン、仮説生成フレームワーク、
批判的思考、研究倫理、先行研究の質的評価、クロスドメインアプローチを
統合的に提供する。

## When to Use

- 新規研究プロジェクトの計画を体系的に立てるとき
- ブレインストーミングでアイデアを構造化するとき
- 研究デザイン（観察/実験/準実験）を選択するとき
- 先行研究の質を系統的に評価するとき
- 倫理審査（IRB/倫理委員会）の準備をするとき
- 研究の妥当性と信頼性を確保する方法を検討するとき

## Quick Start

### 研究計画策定パイプライン

```
Phase 1: Ideation（着想）
  - ブレインストーミング
  - 文献ギャップ分析
  - クロスドメイン着想
    ↓
Phase 2: Research Question（研究課題の定義）
  - PICO/FINER 基準
  - 仮説の形式化
  - スコープの設定
    ↓
Phase 3: Study Design（研究デザイン）
  - デザイン選択
  - サンプルサイズ計算
  - バイアスコントロール
    ↓
Phase 4: Methodology（方法論）
  - データ収集計画
  - 解析計画（SAP）
  - 倫理審査
    ↓
Phase 5: Validity & Rigor（妥当性・厳密性）
  - 内的妥当性
  - 外的妥当性
  - 再現性計画
```

---

## Phase 1: ブレインストーミングフレームワーク

### 構造化ブレインストーミング

```python
BRAINSTORMING_FRAMEWORKS = {
    "SCAMPER": {
        "S": "Substitute（代替）: 何を別のもので置き換えられるか？",
        "C": "Combine（結合）: 何と何を組み合わせられるか？",
        "A": "Adapt（適応）: 他分野のアイデアを適用できるか？",
        "M": "Modify/Magnify（修正/拡大）: 何を変更・拡大できるか？",
        "P": "Put to other uses（転用）: 別の用途はないか？",
        "E": "Eliminate（除去）: 何を取り除けるか？",
        "R": "Reverse/Rearrange（逆転/再配列）: 順序・構造を変えられるか？",
    },
    "Cross_Domain": {
        "description": "他分野のアナロジーから着想するフレームワーク",
        "steps": [
            "1. 自分野の問題を抽象化する",
            "2. 類似の抽象構造を持つ他分野を探す",
            "3. 他分野の解決策を自分野に翻訳する",
            "4. 翻訳した解を具体化・評価する",
        ],
        "examples": [
            "ニューラルネットワーク（CS）→ 遺伝子調節ネットワーク（生物学）",
            "品質管理（製造）→ 臨床試験のモニタリング（医学）",
            "進化アルゴリズム（CS）→ 分子最適化（創薬）",
        ]
    },
    "TRIZ": {
        "description": "矛盾解決の体系的アプローチ",
        "principles": [
            "Segmentation（分割）",
            "Extraction（抽出）",
            "Local quality（局所的性質）",
            "Asymmetry（非対称）",
            "Merging（統合）",
        ]
    },
}
```

### 仮説タイプ分類

```markdown
## Hypothesis Types

| Type | Format | Example |
|------|--------|---------|
| **Directional** | X increases Y | 運動は BMI を低下させる |
| **Non-directional** | X affects Y | 運動は BMI に影響する |
| **Null (H₀)** | X does not affect Y | 運動は BMI に影響しない |
| **Mechanistic** | X causes Y via Z | 運動は代謝率向上を介して BMI を低下させる |
| **Interaction** | X modifies Z's effect on Y | 年齢が運動の BMI 効果を修飾する |
```

---

## Phase 2: 研究デザイン選択

### デザイン分類マトリクス

```python
STUDY_DESIGNS = {
    "experimental": {
        "RCT": {
            "description": "ランダム化比較試験",
            "strength": "因果推論が最も強い",
            "when": "介入効果の評価",
            "bias_control": ["randomization", "blinding", "allocation concealment"],
            "evidence_level": "1b",
        },
        "crossover": {
            "description": "クロスオーバー試験",
            "strength": "被験者内比較（個人差を排除）",
            "when": "慢性疾患、washout が可能",
        },
        "factorial": {
            "description": "要因配置試験",
            "strength": "複数介入の交互作用",
            "when": "2+ 介入を同時に評価",
        },
    },
    "observational": {
        "cohort": {
            "description": "コホート研究",
            "strength": "発生率・リスク比算出可能",
            "when": "曝露→アウトカムの前向き追跡",
            "evidence_level": "2b",
        },
        "case_control": {
            "description": "症例対照研究",
            "strength": "稀な疾患に適する",
            "when": "稀な疾患の危険因子探索",
            "evidence_level": "3b",
        },
        "cross_sectional": {
            "description": "横断研究",
            "strength": "有病率・関連の把握",
            "when": "スナップショット評価",
            "evidence_level": "4",
        },
    },
    "qualitative": {
        "grounded_theory": "データからの理論構築",
        "phenomenology": "体験の本質の記述",
        "ethnography": "文化・社会的文脈の理解",
        "case_study": "特定事例の深い記述",
    },
}

def recommend_design(research_question):
    """
    研究課題のタイプからデザインを推奨。
    """
    rq_type = research_question.get("type")
    recommendations = {
        "causation": ["RCT", "quasi_experimental"],
        "association": ["cohort", "case_control"],
        "prevalence": ["cross_sectional"],
        "prediction": ["cohort", "validation_study"],
        "exploration": ["qualitative", "mixed_methods"],
    }
    return recommendations.get(rq_type, ["cross_sectional"])
```

---

## Phase 3: FINER 基準による評価

### 研究課題の質的評価

```markdown
## FINER Criteria

| Criterion | Question | Score (1-5) |
|-----------|----------|-------------|
| **F**easible | 限られた時間・予算・人員で実行可能か？ | |
| **I**nteresting | 研究者コミュニティにとって興味深いか？ | |
| **N**ovel | 新しい知見を生み出すか？ | |
| **E**thical | 倫理的に問題がないか？ | |
| **R**elevant | 科学的・社会的に意義があるか？ | |
| **Total** | | /25 |

### 判定基準
- 20-25: Excellent — すぐに計画開始
- 15-19: Good — 一部修正で実行可能
- 10-14: Fair — 大幅なリフレーミングが必要
- <10: Poor — 根本的な再検討が必要
```

---

## Phase 4: バイアスコントロール

### バイアス分類と対策

```python
BIAS_TYPES = {
    "selection_bias": {
        "description": "対象者選択の偏り",
        "countermeasures": ["randomization", "stratification", "matching"],
    },
    "information_bias": {
        "description": "情報収集の偏り",
        "countermeasures": ["blinding", "standardized_instruments", "validated_scales"],
    },
    "confounding": {
        "description": "交絡因子による偏り",
        "countermeasures": ["randomization", "stratification", "regression_adjustment",
                           "propensity_score", "instrumental_variables"],
    },
    "performance_bias": {
        "description": "介入実施の偏り",
        "countermeasures": ["double_blinding", "placebo_control"],
    },
    "attrition_bias": {
        "description": "追跡脱落の偏り",
        "countermeasures": ["ITT_analysis", "sensitivity_analysis", "IPCW"],
    },
    "reporting_bias": {
        "description": "報告の偏り（p-hacking、HARKing）",
        "countermeasures": ["pre_registration", "registered_report", "SAP"],
    },
}
```

---

## Phase 5: 研究倫理

### IRB / 倫理審査チェックリスト

```markdown
## Ethics Review Checklist

### 基本原則（ヘルシンキ宣言 / ベルモントレポート）
- [ ] **Respect for persons**: インフォームドコンセント
- [ ] **Beneficence**: リスクの最小化、利益の最大化
- [ ] **Justice**: 公正な対象者選択

### 書類チェック
- [ ] 研究計画書（プロトコル）
- [ ] 同意書・説明文書
- [ ] 個人情報保護計画
- [ ] 利益相反申告
- [ ] データ管理計画
- [ ] 有害事象報告計画

### 特殊な配慮が必要な場合
- [ ] ヒト組織・検体の使用
- [ ] 脆弱な集団（子ども、認知症、囚人）
- [ ] 遺伝情報の取り扱い
- [ ] 国際共同研究
```

---

## Report Template

```markdown
# Research Methodology Plan: [Title]

**PI**: [Name] | **Date**: [date]

## 1. Research Question
## 2. FINER Evaluation
## 3. Study Design
## 4. Methodology
### 4.1 Participants / Samples
### 4.2 Interventions / Exposures
### 4.3 Outcomes
### 4.4 Data Collection
### 4.5 Statistical Analysis Plan
## 5. Bias Control Strategy
## 6. Ethics Considerations
## 7. Timeline
## 8. Feasibility Assessment
```

---

## Completeness Checklist

- [ ] 研究課題: PICO + FINER 評価
- [ ] デザイン: 研究タイプの明示的選択と根拠
- [ ] サンプルサイズ: 正式な検出力分析
- [ ] バイアス: 少なくとも 3 タイプの対策を明記
- [ ] 統計解析計画 (SAP): Pre-registration に対応
- [ ] 倫理: IRB/倫理審査の準備

## Best Practices

1. **Pre-registration**: 主要解析を事前登録 (OSF, ClinicalTrials.gov)
2. **SAP を先に書く**: データ収集前に統計解析計画を完成
3. **Pilot study**: 本研究前にパイロットで feasibility 確認
4. **バイアスを列挙**: 排除できないバイアスは limitation に明記
5. **再現性を設計**: プロトコル・データ・コードの公開計画を含める
