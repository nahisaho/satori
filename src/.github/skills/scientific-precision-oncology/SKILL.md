---
name: scientific-precision-oncology
description: |
  精密腫瘍学スキル。CIViC・OncoKB・cBioPortal・COSMIC・GDC/TCGA を統合し、
  腫瘍ゲノムプロファイリング・分子標的選定・バイオマーカー評価・治療推奨を支援。
  「がんゲノム解析して」「腫瘍プロファイリングして」「OncoKB で検索して」で発火。
---

# Scientific Precision Oncology

精密腫瘍学（Precision Oncology）のための統合解析スキル。
腫瘍ゲノムデータベース（CIViC, OncoKB, cBioPortal, COSMIC, GDC/TCGA）を
横断的に活用し、分子特性に基づく治療戦略の立案を支援する。

## When to Use

- 腫瘍体細胞変異のアクショナビリティ評価
- CIViC / OncoKB エビデンスレベル検索
- cBioPortal / TCGA 変異頻度・共起解析
- バイオマーカー駆動の治療推奨
- 分子標的治療のエビデンス統合
- がん種横断的なドライバー変異分析

## Quick Start

### 精密腫瘍学パイプライン

```
Phase 1: Tumor Profiling
  - 体細胞変異・CNV・融合遺伝子の同定
  - TMB (Tumor Mutational Burden) 算出
  - MSI (Microsatellite Instability) 判定
    ↓
Phase 2: Variant Annotation
  - CIViC GraphQL API クエリ
  - OncoKB Annotation API
  - COSMIC / cBioPortal 変異頻度
    ↓
Phase 3: Actionability Assessment
  - OncoKB Evidence Level (1-4, R1-R2)
  - CIViC Evidence Rating (A-E)
  - AMP/ASCO/CAP Tiering (I-IV)
    ↓
Phase 4: Treatment Selection
  - 分子標的薬マッチング
  - 併用療法候補の同定
  - 耐性メカニズムの評価
    ↓
Phase 5: Clinical Trial Matching
  - ClinicalTrials.gov API 検索
  - 適格基準の自動マッチング
  - バイオマーカー駆動試験の優先順位付け
    ↓
Phase 6: Molecular Tumor Board Report
  - 統合ゲノムレポート生成
  - 治療推奨サマリー
  - エビデンステーブル
```

## Workflow

### 1. CIViC (Clinical Interpretation of Variants in Cancer)

```python
import requests
import pandas as pd

# === CIViC GraphQL API ===
CIVIC_URL = "https://civicdb.org/api/graphql"

def query_civic_variant(gene, variant_name):
    """CIViC で遺伝子バリアントのエビデンスを検索"""
    query = """
    query($gene: String!) {
      genes(name: $gene) {
        nodes {
          name
          description
          variants {
            nodes {
              name
              variantTypes { name }
              evidenceItems {
                nodes {
                  status
                  evidenceType
                  evidenceLevel
                  evidenceDirection
                  significance
                  disease { name }
                  therapies { name }
                  source { citation }
                }
              }
            }
          }
        }
      }
    }
    """
    resp = requests.post(CIVIC_URL, json={"query": query, "variables": {"gene": gene}})
    data = resp.json()["data"]["genes"]["nodes"]

    results = []
    for g in data:
        for v in g["variants"]["nodes"]:
            if variant_name.upper() in v["name"].upper():
                for ev in v["evidenceItems"]["nodes"]:
                    if ev["status"] == "accepted":
                        results.append({
                            "gene": g["name"],
                            "variant": v["name"],
                            "type": ev["evidenceType"],
                            "level": ev["evidenceLevel"],
                            "direction": ev["evidenceDirection"],
                            "significance": ev["significance"],
                            "disease": ev["disease"]["name"] if ev["disease"] else "",
                            "therapies": ", ".join(t["name"] for t in ev["therapies"]),
                            "citation": ev["source"]["citation"] if ev["source"] else "",
                        })
    return pd.DataFrame(results)

civic_results = query_civic_variant("BRAF", "V600E")
print(f"CIViC evidence items: {len(civic_results)}")
print(civic_results[["gene", "variant", "level", "significance", "therapies"]].head(10))
```

### 2. OncoKB Annotation

```python
ONCOKB_URL = "https://www.oncokb.org/api/v1"
ONCOKB_TOKEN = "YOUR_ONCOKB_TOKEN"  # oncokb.org で取得

def annotate_oncokb(gene, variant, tumor_type=None):
    """OncoKB でバリアントをアノテーション"""
    headers = {"Authorization": f"Bearer {ONCOKB_TOKEN}"}

    # Variant annotation
    params = {"hugoSymbol": gene, "alteration": variant}
    if tumor_type:
        params["tumorType"] = tumor_type

    resp = requests.get(f"{ONCOKB_URL}/annotate/mutations/byHGVSg",
                        headers=headers, params=params)
    data = resp.json()

    return {
        "gene": gene,
        "variant": variant,
        "oncogenic": data.get("oncogenic", ""),
        "mutation_effect": data.get("mutationEffect", {}).get("knownEffect", ""),
        "highest_sensitive_level": data.get("highestSensitiveLevel", ""),
        "highest_resistance_level": data.get("highestResistanceLevel", ""),
        "treatments": [
            {
                "drugs": ", ".join(d["drugName"] for d in t.get("drugs", [])),
                "level": t.get("level", ""),
                "indication": t.get("levelAssociatedCancerType", {}).get("name", ""),
            }
            for t in data.get("treatments", [])
        ],
    }

# OncoKB Evidence Levels
ONCOKB_LEVELS = {
    "LEVEL_1": "FDA-recognized biomarker (same indication)",
    "LEVEL_2": "Standard care biomarker (same indication)",
    "LEVEL_3A": "Compelling clinical evidence (same indication)",
    "LEVEL_3B": "Standard care or compelling evidence (different indication)",
    "LEVEL_4": "Compelling biological evidence",
    "LEVEL_R1": "Standard care resistance biomarker",
    "LEVEL_R2": "Compelling clinical resistance evidence",
}
```

### 3. cBioPortal 解析

```python
CBIOPORTAL_URL = "https://www.cbioportal.org/api"

def query_cbioportal_mutations(gene, study_ids=None):
    """cBioPortal で変異頻度を取得"""
    if study_ids is None:
        # TCGA PanCancer Atlas studies
        resp = requests.get(f"{CBIOPORTAL_URL}/studies",
                            params={"keyword": "tcga_pan_can_atlas"})
        study_ids = [s["studyId"] for s in resp.json()]

    all_mutations = []
    for study_id in study_ids:
        # Get molecular profile
        profiles = requests.get(
            f"{CBIOPORTAL_URL}/molecular-profiles",
            params={"studyId": study_id}
        ).json()

        mut_profiles = [p for p in profiles if p["molecularAlterationType"] == "MUTATION_EXTENDED"]
        if not mut_profiles:
            continue

        profile_id = mut_profiles[0]["molecularProfileId"]

        # Get mutations
        mutations = requests.get(
            f"{CBIOPORTAL_URL}/molecular-profiles/{profile_id}/mutations",
            params={"entrezGeneId": gene_to_entrez(gene)}
        ).json()

        for m in mutations:
            all_mutations.append({
                "study": study_id,
                "sample_id": m.get("sampleId", ""),
                "protein_change": m.get("proteinChange", ""),
                "mutation_type": m.get("mutationType", ""),
                "variant_type": m.get("variantType", ""),
            })

    df = pd.DataFrame(all_mutations)

    # 変異頻度集計
    if not df.empty:
        freq = df["protein_change"].value_counts().head(20)
        print(f"Top mutations for {gene}:")
        print(freq)

    return df

def gene_to_entrez(gene_symbol):
    """Hugo symbol → Entrez ID 変換"""
    mapping = {"BRAF": 673, "EGFR": 1956, "KRAS": 3845, "TP53": 7157, "PIK3CA": 5290}
    return mapping.get(gene_symbol, 0)
```

### 4. TMB / MSI 算出

```python
def calculate_tmb(mutations_df, exome_size_mb=38.0):
    """
    Tumor Mutational Burden (TMB) 算出
    TMB = nonsynonymous mutations / exome size (Mb)
    """
    # Nonsynonymous のみ
    nonsynonymous = mutations_df[
        mutations_df["mutation_type"].isin([
            "Missense_Mutation", "Nonsense_Mutation",
            "Frame_Shift_Del", "Frame_Shift_Ins",
            "In_Frame_Del", "In_Frame_Ins",
            "Splice_Site", "Translation_Start_Site",
        ])
    ]

    tmb = len(nonsynonymous) / exome_size_mb

    # TMB カテゴリ (FDA: ≥10 mut/Mb = TMB-High → Pembrolizumab)
    if tmb >= 20:
        category = "TMB-Very High"
    elif tmb >= 10:
        category = "TMB-High"
    elif tmb >= 5:
        category = "TMB-Intermediate"
    else:
        category = "TMB-Low"

    return {"tmb": round(tmb, 2), "category": category,
            "nonsynonymous_count": len(nonsynonymous)}


def assess_msi(microsatellite_loci_results):
    """
    MSI (Microsatellite Instability) 判定
    Bethesda Panel: BAT25, BAT26, D2S123, D5S346, D17S250
    """
    unstable_count = sum(1 for r in microsatellite_loci_results if r["status"] == "unstable")
    total = len(microsatellite_loci_results)

    if unstable_count >= 2:
        status = "MSI-H"  # High
    elif unstable_count == 1:
        status = "MSI-L"  # Low
    else:
        status = "MSS"    # Stable

    return {"status": status, "unstable_loci": unstable_count, "total_loci": total}
```

### 5. 分子腫瘍ボードレポート生成

```python
import json

def generate_mtb_report(patient_id, variants, civic_data, oncokb_data,
                        tmb_result, msi_result, output_dir="results"):
    """Molecular Tumor Board (MTB) レポート生成"""
    report = {
        "patient_id": patient_id,
        "report_date": pd.Timestamp.now().isoformat(),
        "genomic_profile": {
            "tmb": tmb_result,
            "msi": msi_result,
            "variants": variants,
        },
        "actionable_findings": [],
        "clinical_trials": [],
    }

    # Actionability 統合
    for v in variants:
        finding = {
            "gene": v["gene"],
            "variant": v["variant"],
            "oncokb_level": oncokb_data.get(f"{v['gene']}_{v['variant']}", {}).get("highest_sensitive_level", ""),
            "civic_evidence": [],
            "therapies": [],
        }

        # CIViC エビデンス統合
        civic_match = civic_data[
            (civic_data["gene"] == v["gene"]) &
            (civic_data["variant"].str.contains(v["variant"], case=False))
        ]
        for _, ev in civic_match.iterrows():
            finding["civic_evidence"].append({
                "level": ev["level"],
                "therapies": ev["therapies"],
                "significance": ev["significance"],
            })

        report["actionable_findings"].append(finding)

    # AMP/ASCO/CAP Tier 分類
    for finding in report["actionable_findings"]:
        level = finding["oncokb_level"]
        if level in ["LEVEL_1", "LEVEL_2"]:
            finding["amp_tier"] = "Tier I (Strong clinical significance)"
        elif level in ["LEVEL_3A", "LEVEL_3B"]:
            finding["amp_tier"] = "Tier II (Potential clinical significance)"
        elif level == "LEVEL_4":
            finding["amp_tier"] = "Tier III (Unknown significance)"
        else:
            finding["amp_tier"] = "Tier IV (Benign/likely benign)"

    # JSON 出力
    with open(f"{output_dir}/mtb_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    # Markdown レポート
    md = f"# Molecular Tumor Board Report\n\n"
    md += f"**Patient**: {patient_id} | **Date**: {report['report_date']}\n\n"
    md += f"## Genomic Profile\n\n"
    md += f"| Metric | Result |\n|---|---|\n"
    md += f"| TMB | {tmb_result['tmb']} mut/Mb ({tmb_result['category']}) |\n"
    md += f"| MSI | {msi_result['status']} |\n\n"
    md += f"## Actionable Findings\n\n"
    md += "| Gene | Variant | OncoKB Level | AMP Tier | Therapies |\n|---|---|---|---|---|\n"
    for f_ in report["actionable_findings"]:
        therapies = "; ".join(set(
            e["therapies"] for e in f_["civic_evidence"] if e["therapies"]
        ))
        md += f"| {f_['gene']} | {f_['variant']} | {f_['oncokb_level']} | {f_['amp_tier']} | {therapies} |\n"

    with open(f"{output_dir}/mtb_report.md", "w") as f_out:
        f_out.write(md)

    return report
```

---

## Best Practices

1. **多データベースクロスバリデーション**: CIViC + OncoKB + COSMIC の一致を重視
2. **がん種特異的解釈**: 同じ変異でもがん種により臨床的意義が異なる
3. **エビデンスレベルの階層**: Level 1 > 2 > 3A > 3B > 4 の優先順位に従う
4. **耐性変異を見逃さない**: 一次耐性・獲得耐性の両方を評価
5. **TMB/MSI を免疫療法判断の補助指標に**: TMB-High (≥10) → Pembrolizumab 適応
6. **VUS (Variant of Unknown Significance) の扱い**: 機能予測ツールを補助的に使用

## Completeness Checklist

- [ ] 体細胞変異リスト確定（SNV/Indel/CNV/Fusion）
- [ ] CIViC エビデンス取得
- [ ] OncoKB アノテーション完了
- [ ] cBioPortal 変異頻度確認
- [ ] TMB / MSI 判定
- [ ] AMP/ASCO/CAP Tier 分類
- [ ] 治療推奨サマリー（分子標的薬 + 免疫療法）
- [ ] 臨床試験マッチング
- [ ] MTB レポート（JSON + Markdown）生成

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/mtb_report.json` | 分子腫瘍ボードレポート（JSON） | プロファイリング完了時 |
| `results/mtb_report.md` | MTB レポート（Markdown） | レポート生成時 |
| `results/variant_actionability.json` | バリアント臨床的意義（JSON） | アノテーション完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| OncoKB | `OncoKB_annotate_variant` | 体細胞変異の臨床的アノテーション |
| OncoKB | `OncoKB_get_cancer_genes` | がん遺伝子リスト取得 |
| CIViC | `civic_search_evidence_items` | 臨床エビデンス検索 |
| CIViC | `civic_get_variant` | バリアント臨床解釈 |
| COSMIC | `COSMIC_get_mutations_by_gene` | 体細胞変異頻度データ |
| GDC | `GDC_get_mutation_frequency` | TCGA 変異頻度 |
| ClinicalTrials | `search_clinical_trials` | 腫瘍学臨床試験マッチング |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-variant-interpretation` | ← 生殖細胞系変異の ACMG 分類 |
| `scientific-clinical-decision-support` | → 治療推奨の臨床意思決定反映 |
| `scientific-bioinformatics` | ← RNA-seq 発現データ解析 |
| `scientific-network-analysis` | ← シグナル経路解析・ドライバー予測 |
| `scientific-drug-target-profiling` | ← 標的ドラッガビリティ評価 |
| `scientific-disease-research` | ← がん種の疫学・遺伝的背景 |
| `scientific-deep-research` | ← 腫瘍学最新文献リサーチ |
| `scientific-pharmacogenomics` | ← PGx 代謝型・投与量調整 |
