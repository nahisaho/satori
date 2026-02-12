---
name: scientific-drug-target-profiling
description: |
  創薬ターゲットプロファイリングスキル。ToolUniverse / Open Targets / ChEMBL / UniProt
  を活用したドラッグターゲットインテリジェンス。ドラッガビリティ評価、安全性プロファイリング、
  ターゲット-疾患アソシエーション、競合パイプライン分析を統合的に実行。
  「ターゲット評価して」「druggability 分析して」「標的タンパク質を調べて」で発火。
---

# Scientific Drug Target Profiling

創薬ターゲットの包括的プロファイリングスキル。ToolUniverse（mims-harvard）の
Target Intelligence Gatherer パラダイムに準拠し、9 つの並行リサーチパスで
ターゲットを多角的に評価する。

## When to Use

- 創薬ターゲット候補のドラッガビリティを評価するとき
- ターゲットタンパク質の安全性プロファイルを確認するとき
- ターゲット-疾患アソシエーションの強度を定量化するとき
- 既知リガンド・ケミカルプローブを網羅的に収集するとき
- 競合パイプライン（臨床段階の化合物）を調査するとき

## Quick Start

### 1. ターゲットプロファイリング 9 パス戦略

```
PATH 1: Identity Resolution
  Gene Symbol → UniProt → Ensembl → ChEMBL Target ID

PATH 2: Basic Protein Information
  UniProt Entry → Function, Localization, Domains

PATH 3: Structural Biology
  PDB Structures → AlphaFold → Binding Sites

PATH 4: Function & Pathways
  GO Terms → Reactome → KEGG → Pathway Context

PATH 5: Expression Profile
  GTEx → HPA → Tissue Specificity → Single-cell

PATH 6: Genetic Variation & Disease
  ClinVar → gnomAD → GWAS → Constraint Scores

PATH 7: Drug Interactions & Druggability
  ChEMBL Activities → DGIdb → Known Drugs → Probes

PATH 8: Literature & Research Landscape
  PubMed → OpenAlex → Publication Trends

PATH 9: Safety & Toxicology
  Essential Genes → Phenotypes → Off-target Risk
```

---

## Phase 1: Identity Resolution

### ID マッピングチェーン

```python
def resolve_target_ids(query):
    """
    ターゲット名/Gene Symbol から主要 ID を解決する。
    """
    ids = {
        "query": query,
        "uniprot_accession": None,     # e.g., P04637
        "ensembl_id": None,            # e.g., ENSG00000141510
        "entrez_id": None,             # e.g., 7157
        "chembl_target_id": None,      # e.g., CHEMBL3927
        "hgnc_symbol": None,           # e.g., TP53
        "open_targets_id": None,       # = Ensembl ID
    }

    # Step 1: UniProt 検索（Gene Name → Accession）
    # Step 2: Ensembl ID 取得（UniProt xref）
    # Step 3: ChEMBL Target ID（UniProt → ChEMBL mapping）
    # Step 4: Cross-validation（各 DB 間で一致確認）

    return ids
```

> **名前衝突の検出**: 同名の遺伝子が複数生物種に存在する場合、
> UniProt taxonomy filter で human (9606) を優先する。

---

## Phase 2: Druggability Assessment

### ドラッガビリティ 3 軸評価

```
┌─────────────────────────────────────────┐
│         Druggability Matrix             │
├─────────────┬───────────┬───────────────┤
│  Modality   │  Metric   │   Threshold   │
├─────────────┼───────────┼───────────────┤
│ Small Mol   │ Pocket?   │ ≥1 druggable  │
│ Antibody    │ Surface?  │ extracellular │
│ PROTAC      │ E3 dist   │ ≤30 Å         │
│ ASO/siRNA   │ mRNA expr │ detectable    │
│ Gene Therapy│ LOF/GOF   │ disease link  │
└─────────────┴───────────┴───────────────┘
```

### Target Development Level (TDL) 分類

```python
def classify_tdl(target_data):
    """
    Pharos TDL 分類に準拠したターゲット分類。
    Tclin: 承認薬あり
    Tchem: 高活性化合物あり（ChEMBL）
    Tbio: 生物学的機能が判明
    Tdark: 情報不足
    """
    if target_data.get("approved_drugs"):
        return "Tclin"
    elif target_data.get("chembl_activities_count", 0) > 0:
        potent = [a for a in target_data["activities"]
                  if a.get("pchembl_value", 0) >= 6.0]
        if potent:
            return "Tchem"
    if target_data.get("go_annotations") or target_data.get("publications", 0) > 5:
        return "Tbio"
    return "Tdark"
```

---

## Phase 3: Safety Profiling

### 安全性評価チェックリスト

```markdown
## Safety Assessment

### Genetic Constraint
- [ ] pLI score: ___  (>0.9 = highly constrained, LOF intolerant)
- [ ] LOEUF: ___  (<0.35 = constrained)
- [ ] Missense Z-score: ___  (>3.09 = missense-constrained)

### Essential Gene Analysis
- [ ] DepMap dependency score: ___  (<-0.5 = broadly essential)
- [ ] Mouse knockout phenotype: ___
- [ ] Lethal phenotype: ___  (YES/NO)

### Expression Breadth
- [ ] Tissue specificity index (tau): ___
- [ ] Ubiquitously expressed: ___  (risk for on-target toxicity)
- [ ] Brain/Heart/Liver expression: ___  (safety-critical organs)

### Off-target Risk
- [ ] Paralog count: ___
- [ ] Closest paralog similarity: ___  %
- [ ] Shared binding site features: ___
```

---

## Phase 4: Disease Association

### エビデンスグレーディング

```python
EVIDENCE_TIERS = {
    "T1": "Genetic + Clinical (GWAS + ClinVar pathogenic)",
    "T2": "Strong biological (functional studies + animal models)",
    "T3": "Associative (expression correlation + network guilt-by-association)",
    "T4": "Computational prediction only",
}

def grade_disease_association(target_id, disease_id, evidence_sources):
    """
    ターゲット-疾患アソシエーションのエビデンス評価。
    Open Targets overall_association_score + 追加エビデンスで T1-T4 判定。
    """
    score = evidence_sources.get("open_targets_score", 0)
    has_gwas = evidence_sources.get("gwas_significance", False)
    has_clinvar = evidence_sources.get("clinvar_pathogenic", False)
    has_functional = evidence_sources.get("functional_study", False)

    if has_gwas and has_clinvar:
        return "T1", score
    elif has_functional or score > 0.7:
        return "T2", score
    elif score > 0.3:
        return "T3", score
    else:
        return "T4", score
```

---

## Phase 5: Competitive Landscape

### パイプラインマッピング

```markdown
## Competitive Intelligence

### Known Drugs (Approved)
| Drug | Mechanism | Indication | Approval Year |
|------|-----------|------------|---------------|

### Clinical Pipeline
| Compound | Phase | Sponsor | Indication | NCT ID |
|----------|-------|---------|------------|--------|

### Chemical Probes
| Probe | Potency | Selectivity | Source |
|-------|---------|-------------|--------|

### Patent Landscape
| Patent Family | Assignee | Filing Date | Key Claims |
|---------------|----------|-------------|------------|
```

---

## Report Template

### ターゲットインテリジェンスレポート

```markdown
# Target Intelligence Report: [TARGET NAME]

**Generated**: [Date] | **Analyst**: SATORI Drug Target Profiling

## 1. Executive Summary
[2-3 sentences: target name, key disease links (with evidence tier), druggability verdict]

## 2. Target Identifiers
| Database | ID | Verified |
|----------|----|----------|
| UniProt  |    | ✓/✗     |
| Ensembl  |    | ✓/✗     |
| ChEMBL   |    | ✓/✗     |

## 3. Protein Biology
### 3.1 Function & Localization
### 3.2 Domain Architecture
### 3.3 Pathway Context

## 4. Structural Biology
### 4.1 Experimental Structures (PDB)
### 4.2 AlphaFold Prediction
### 4.3 Binding Sites & Pockets

## 5. Expression Profile
### 5.1 Tissue Expression (GTEx/HPA)
### 5.2 Disease-specific Expression

## 6. Disease Associations
[Table with evidence tiers T1-T4]

## 7. Druggability Assessment
### 7.1 TDL Classification
### 7.2 Modality Assessment
### 7.3 Tractability Score

## 8. Known Ligands & Drugs
### 8.1 Approved Drugs
### 8.2 Clinical Candidates
### 8.3 Chemical Probes & Tool Compounds

## 9. Safety Profile
### 9.1 Genetic Constraint
### 9.2 Essential Gene Status
### 9.3 Off-target Risk

## 10. Competitive Landscape

## 11. Recommendations
### 11.1 Go/No-Go Assessment
### 11.2 Suggested Modality
### 11.3 Key Experiments Needed

## 12. Data Sources & Methodology
```

---

## Completeness Checklist

### 必須項目

- [ ] ID Resolution: UniProt, Ensembl, ChEMBL, HGNC の 4ID 確認
- [ ] Druggability: TDL 分類 + 少なくとも 2 モダリティの評価
- [ ] Safety: pLI + LOEUF + DepMap の 3 指標
- [ ] Disease: 上位 5 疾患の T1-T4 グレード
- [ ] Literature: 主要レビュー論文 ≥3 件引用
- [ ] Competitive: 承認薬と Ph3 候補を網羅

## Best Practices

1. **ID は必ず Cross-validate**: UniProt と Ensembl の双方向マッピングで確認
2. **名前衝突を検出**: Gene Symbol の同名異義語を taxonomy filter で排除
3. **Evidence Tier を明記**: すべての疾患アソシエーションに T1-T4 を付与
4. **安全性を最初に確認**: Essential gene の場合は早期に Go/No-Go 判断
5. **構造情報を優先**: PDB 実験構造 > AlphaFold 予測 > Homology model

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/target_profile_report.md` | ターゲットプロファイルレポート（Markdown） | 全解析完了時 |
| `results/target_profile.json` | 構造化プロファイルデータ（JSON） | 全解析完了時 |
| `results/druggability_matrix.json` | ドラッガビリティマトリクス（JSON） | Druggability 評価完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| UniProt | `UniProt_get_entry_by_accession` | タンパク質エントリ取得 |
| UniProt | `UniProt_get_function_by_accession` | タンパク質機能情報 |
| ChEMBL | `ChEMBL_get_target` | ターゲット情報取得 |
| ChEMBL | `ChEMBL_get_target_activities` | ターゲット活性データ |
| OpenTargets | `OpenTargets_get_associated_targets_by_disease_efoId` | 疾患-ターゲット関連 |
| DGIdb | `DGIdb_get_gene_druggability` | ドラッガビリティ評価 |
| DGIdb | `DGIdb_get_drug_gene_interactions` | 薬物-遺伝子相互作用 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-hypothesis-pipeline` | ← 仮説定義からターゲット同定への入力 |
| `scientific-deep-research` | ← 文献深層調査で標的エビデンス収集 |
| `scientific-bioinformatics` | ← ゲノム・プロテオームデータ提供 |
| `scientific-network-analysis` | ← PPI ネットワーク・パスウェイ情報 |
| `scientific-admet-pharmacokinetics` | → ターゲットに対する化合物の ADMET 評価 |
| `scientific-protein-structure-analysis` | → ターゲットタンパク質の構造解析 |
| `scientific-drug-repurposing` | → ターゲットベースのリポジショニング |
| `scientific-academic-writing` | → 研究成果の論文化 |
