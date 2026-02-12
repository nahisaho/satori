---
name: scientific-variant-interpretation
description: |
  遺伝子バリアント臨床解釈スキル。ClinVar / gnomAD / COSMIC / ACMG ガイドラインに
  基づく病原性評価、薬理ゲノミクス（PharmGKB/ClinPGx）、バリアント-表現型相関の
  エビデンスグレーディング。ToolUniverse の Variant Interpretation パラダイムを統合。
  「バリアントの病原性を評価して」「pharmacogenomics 解析して」で発火。
tu_tools:
  - key: clinvar
    name: ClinVar
    description: 臨床的バリアント解釈データベース
---

# Scientific Variant Interpretation

遺伝子バリアントの臨床的解釈スキル。ACMG/AMP ガイドラインに準拠した
病原性分類、薬理ゲノミクスによる薬物応答予測、体細胞変異の
ドライバー/パッセンジャー判定を統合的に実行する。

## When to Use

- SNV/Indel の病原性を ACMG 基準で分類するとき
- 薬理ゲノミクス (PGx) による薬物応答予測を行うとき
- がん体細胞変異のドライバー判定を行うとき
- 希少疾患の原因バリアント同定を支援するとき
- バリアント-表現型相関のエビデンスを評価するとき

## Quick Start

### バリアント解釈パイプライン

```
Input: Variant (HGVS notation / rsID / genomic coordinates)
    ↓
Step 1: Annotation
  - ClinVar 臨床意義
  - gnomAD 集団頻度
  - InterVar ACMG 自動分類
    ↓
Step 2: Population Frequency
  - gnomAD allele frequency
  - 民族・集団別頻度
  - 稀少性判定 (AF < 0.01 or 0.001)
    ↓
Step 3: Functional Impact
  - CADD / REVEL / AlphaMissense スコア
  - 保存度 (PhyloP, GERP++)
  - スプライス予測 (SpliceAI)
    ↓
Step 4: Clinical Evidence
  - ClinVar submissions
  - COSMIC (体細胞変異)
  - Literature evidence  (PubMed)
    ↓
Step 5: ACMG Classification
  - 28 基準の系統的適用
  - 5 段階分類 (P, LP, VUS, LB, B)
    ↓
Output: Variant Interpretation Report
```

---

## Phase 1: ACMG/AMP 分類

### 28 基準の系統的評価

```python
ACMG_CRITERIA = {
    # Pathogenic - Very Strong
    "PVS1": "Null variant (nonsense, frameshift, canonical splice) in LOF-intolerant gene",
    # Pathogenic - Strong
    "PS1": "Same amino acid change as established pathogenic variant",
    "PS2": "De novo (confirmed parentage) in patient with disease",
    "PS3": "Functional study shows damaging effect",
    "PS4": "Prevalence in affected >> controls (OR > 5)",
    # Pathogenic - Moderate
    "PM1": "In mutational hot spot or critical functional domain",
    "PM2": "Absent from controls (or extremely rare in gnomAD)",
    "PM3": "Detected in trans with pathogenic variant (recessive)",
    "PM4": "Protein length change (in-frame del/ins in non-repeat region)",
    "PM5": "Novel missense at same position as established pathogenic",
    "PM6": "Assumed de novo (parentage not confirmed)",
    # Pathogenic - Supporting
    "PP1": "Co-segregation with disease in multiple family members",
    "PP2": "Missense in gene with low rate of benign missense",
    "PP3": "Computational evidence supports deleterious effect",
    "PP4": "Patient phenotype highly specific for gene",
    "PP5": "Reputable source reports as pathogenic",
    # Benign - Stand-alone
    "BA1": "Allele frequency > 5% in any population",
    # Benign - Strong
    "BS1": "Allele frequency greater than expected for disorder",
    "BS2": "Observed in healthy adult (for early-onset/penetrant disorder)",
    "BS3": "Functional study shows no damaging effect",
    "BS4": "Lack of segregation in affected family members",
    # Benign - Supporting
    "BP1": "Missense in gene where only truncating cause disease",
    "BP2": "Observed in trans with pathogenic variant (dominant)",
    "BP3": "In-frame del/ins in repetitive region",
    "BP4": "Computational evidence suggests no impact",
    "BP5": "Variant found in case with alternate molecular basis",
    "BP6": "Reputable source reports as benign",
    "BP7": "Synonymous with no splicing impact predicted",
}

def classify_acmg(criteria_met):
    """
    適用された ACMG 基準から最終分類を導出。
    Richards et al., Genetics in Medicine 2015
    """
    pathogenic_criteria = [c for c in criteria_met if c.startswith(("PVS", "PS", "PM", "PP"))]
    benign_criteria = [c for c in criteria_met if c.startswith(("BA", "BS", "BP"))]

    pvs = [c for c in pathogenic_criteria if c.startswith("PVS")]
    ps = [c for c in pathogenic_criteria if c.startswith("PS")]
    pm = [c for c in pathogenic_criteria if c.startswith("PM")]
    pp = [c for c in pathogenic_criteria if c.startswith("PP")]

    ba = [c for c in benign_criteria if c.startswith("BA")]
    bs = [c for c in benign_criteria if c.startswith("BS")]
    bp = [c for c in benign_criteria if c.startswith("BP")]

    # Pathogenic Rules
    if (len(pvs) >= 1 and (len(ps) >= 1 or len(pm) >= 2 or
        (len(pm) == 1 and len(pp) == 1) or len(pp) >= 2)):
        return "Pathogenic"
    if len(ps) >= 2:
        return "Pathogenic"

    # Likely Pathogenic
    if len(pvs) >= 1 and len(pm) == 1:
        return "Likely Pathogenic"
    if len(ps) >= 1 and (len(pm) >= 1 or len(pm) >= 2):
        return "Likely Pathogenic"

    # Benign
    if len(ba) >= 1:
        return "Benign"
    if len(bs) >= 2:
        return "Benign"

    # Likely Benign
    if len(bs) >= 1 and len(bp) >= 1:
        return "Likely Benign"
    if len(bp) >= 2:
        return "Likely Benign"

    return "Variant of Uncertain Significance (VUS)"
```

---

## Phase 2: 薬理ゲノミクス (PGx)

### PharmGKB / CPIC ガイドライン

```python
PGX_GENES = {
    "CYP2D6": {
        "drugs": ["codeine", "tramadol", "tamoxifen", "atomoxetine"],
        "phenotypes": ["Poor Metabolizer", "Intermediate", "Normal", "Ultrarapid"],
    },
    "CYP2C19": {
        "drugs": ["clopidogrel", "voriconazole", "escitalopram", "omeprazole"],
        "phenotypes": ["Poor Metabolizer", "Intermediate", "Normal", "Rapid", "Ultrarapid"],
    },
    "CYP2C9": {
        "drugs": ["warfarin", "phenytoin", "celecoxib"],
        "phenotypes": ["Poor Metabolizer", "Intermediate", "Normal"],
    },
    "DPYD": {
        "drugs": ["fluorouracil", "capecitabine"],
        "phenotypes": ["Poor Metabolizer", "Intermediate", "Normal"],
    },
    "TPMT": {
        "drugs": ["azathioprine", "mercaptopurine", "thioguanine"],
        "phenotypes": ["Poor Metabolizer", "Intermediate", "Normal"],
    },
    "HLA-B": {
        "variants": {"*57:01": "abacavir hypersensitivity", "*58:01": "allopurinol SJS/TEN"},
    },
}

def pgx_recommendation(gene, phenotype, drug):
    """
    CPIC ガイドラインに基づく薬物用量推奨。
    """
    recommendations = {
        "avoid": "代替薬の使用を推奨",
        "reduce_dose": "標準用量の 25-50% に減量",
        "standard_dose": "標準用量で開始",
        "increase_dose": "標準用量で効果不十分の可能性、増量を検討",
    }
    # 実際の推奨は CPIC テーブルから取得
    return recommendations
```

---

## Phase 3: 体細胞変異解釈

### OncoKB / COSMIC エビデンスレベル

```markdown
## Somatic Variant Interpretation

### Oncogenicity Classification
| Level | Description |
|-------|-------------|
| Oncogenic | Functionally validated driver |
| Likely Oncogenic | Strong computational/indirect evidence |
| VUS | Insufficient evidence |
| Likely Benign | Evidence against oncogenicity |
| Benign | Confirmed passenger |

### Therapeutic Actionability (OncoKB Levels)
| Level | Description |
|-------|-------------|
| 1 | FDA-approved, same tumor type |
| 2 | Standard care, different tumor type |
| 3A | Clinical evidence in same tumor type |
| 3B | Clinical evidence in different tumor type |
| 4 | Preclinical evidence |
| R1 | Resistance to approved therapy |
| R2 | Resistance to investigational therapy |
```

---

## Report Template

```markdown
# Variant Interpretation Report

**Variant**: [HGVS notation]
**Gene**: [gene symbol]
**Date**: [date]

## 1. Variant Summary
| Feature | Value |
|---------|-------|
| Genomic location | |
| Transcript | |
| Protein change | |
| Variant type | |

## 2. Population Frequency
| Database | Frequency | Population |
|----------|-----------|------------|

## 3. In Silico Predictions
| Tool | Score | Prediction |
|------|-------|------------|
| CADD | | |
| REVEL | | |
| AlphaMissense | | |
| SpliceAI | | |

## 4. Clinical Evidence
### 4.1 ClinVar
### 4.2 Literature
### 4.3 Functional Studies

## 5. ACMG Classification
| Criterion | Applied | Evidence |
|-----------|---------|----------|
**Final Classification**: [P/LP/VUS/LB/B]

## 6. Pharmacogenomic Implications
（該当する場合）

## 7. Treatment Implications
（がん体細胞変異の場合）

## 8. Recommendations
```

---

## Completeness Checklist

- [ ] バリアント注釈: HGVS、rsID、ゲノム座標の正規化
- [ ] 集団頻度: gnomAD 全集団 + 民族別
- [ ] In silico 予測: CADD + REVEL + 少なくとも 1 追加ツール
- [ ] ClinVar: 全サブミッションの確認
- [ ] ACMG 分類: 28 基準の系統的評価
- [ ] PGx: 該当遺伝子の場合 CPIC ガイドライン参照

## Best Practices

1. **ACMG は系統的に**: 全 28 基準を 1 つずつ評価し、根拠を明記
2. **集団頻度はサブ集団で確認**: 特定民族でのみ高頻度な場合がある
3. **機能研究を重視**: 計算予測よりも実験的エビデンスが優先
4. **ClinVar は星評価を確認**: 3-4 星のエントリが最も信頼性が高い
5. **PGx は CPIC レベルを確認**: Level A/B のみ臨床実装

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/variant_report.md` | バリアント解釈レポート（Markdown） | 全解析完了時 |
| `results/variant_classification.json` | ACMG/AMP 分類データ（JSON） | 分類完了時 |
| `results/pgx_report.json` | 薬理ゲノミクスレポート（JSON） | PGx 評価完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| ClinVar | `clinvar_search_variants` | バリアントの病原性分類検索 |
| gnomAD | `gnomad_get_gene_constraints` | 遺伝子制約メトリクス（pLI / LOEUF） |
| ClinGen | `ClinGen_get_gene_validity` | 遺伝子-疾患の妥当性評価 |
| AlphaMissense | `AlphaMissense_get_variant_score` | ミスセンス病原性予測スコア |
| PharmGKB | `PharmGKB_search_variants` | 薬理ゲノミクスバリアント検索 |
| CADD | `CADD_get_variant_score` | バリアント有害性スコア |
| MyVariant | `MyVariant_get_variant_annotation` | 統合バリアントアノテーション |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-bioinformatics` | ← ゲノムデータ・バリアントコール |
| `scientific-sequence-analysis` | ← 配列コンテキスト・保存度情報 |
| `scientific-data-preprocessing` | ← バリアントデータの前処理・正規化 |
| `scientific-clinical-decision-support` | → バリアント解釈結果の臨床意思決定 |
| `scientific-academic-writing` | → 研究成果の論文化 |
| `scientific-pharmacogenomics` | ← Star アレル・代謝型・薬理ゲノミクス |
