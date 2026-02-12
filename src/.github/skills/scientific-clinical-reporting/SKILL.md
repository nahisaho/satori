---
name: scientific-clinical-reporting
description: |
  臨床レポート自動生成スキル。検査結果サマリー (SOAP ノート)、バイオマーカー
  プロファイルレポート、薬理ゲノミクスレポート、臨床試験要約を構造化テンプレート
  (PDF/LaTeX/HTML) で出力。HL7 FHIR DiagnosticReport 形式にも対応。
---

# Scientific Clinical Reporting

臨床データから構造化レポートを自動生成するパイプラインを提供する。

## When to Use

- 検査結果を SOAP ノート形式でまとめるとき
- バイオマーカープロファイルレポートを作成するとき
- ファーマコゲノミクスレポート (CPIC ガイドライン準拠) が必要なとき
- 臨床試験の CSR (Clinical Study Report) サマリーを生成するとき
- HL7 FHIR DiagnosticReport 形式で出力するとき

---

## Quick Start

## 1. SOAP ノート生成

```python
import json
from datetime import datetime


def generate_soap_note(patient_data, findings, assessment, plan):
    """
    SOAP ノート形式の臨床レポートを生成。

    Parameters:
        patient_data: dict — {"id": "...", "age": 45, "sex": "M", ...}
        findings: dict — {"subjective": [...], "objective": [...]}
        assessment: list — 評価・診断リスト
        plan: list — 治療計画リスト
    """
    soap = {
        "report_type": "SOAP_Note",
        "generated_at": datetime.now().isoformat(),
        "patient": {
            "id": patient_data.get("id", "ANON"),
            "age": patient_data.get("age"),
            "sex": patient_data.get("sex"),
        },
        "S": {  # Subjective
            "chief_complaint": findings.get("chief_complaint", ""),
            "history": findings.get("subjective", []),
        },
        "O": {  # Objective
            "vitals": findings.get("vitals", {}),
            "lab_results": findings.get("lab_results", []),
            "imaging": findings.get("imaging", []),
            "physical_exam": findings.get("objective", []),
        },
        "A": {  # Assessment
            "diagnoses": assessment,
            "differential": findings.get("differential", []),
        },
        "P": {  # Plan
            "treatment": plan,
            "follow_up": findings.get("follow_up", ""),
            "referrals": findings.get("referrals", []),
        },
    }

    print(f"SOAP note: patient={soap['patient']['id']}, "
          f"diagnoses={len(assessment)}, plans={len(plan)}")
    return soap
```

## 2. バイオマーカープロファイルレポート

```python
import pandas as pd


def biomarker_profile_report(biomarkers_df, reference_ranges=None):
    """
    バイオマーカープロファイルレポート生成。

    Parameters:
        biomarkers_df: DataFrame — columns: [marker, value, unit, specimen]
        reference_ranges: dict — {"marker": {"low": x, "high": y, "unit": "..."}}
    """
    if reference_ranges is None:
        reference_ranges = {
            "CEA": {"low": 0, "high": 5.0, "unit": "ng/mL"},
            "AFP": {"low": 0, "high": 10.0, "unit": "ng/mL"},
            "CA19-9": {"low": 0, "high": 37.0, "unit": "U/mL"},
            "CA125": {"low": 0, "high": 35.0, "unit": "U/mL"},
            "PSA": {"low": 0, "high": 4.0, "unit": "ng/mL"},
            "HER2": {"low": 0, "high": 1, "unit": "IHC score"},
            "Ki-67": {"low": 0, "high": 14, "unit": "%"},
            "PD-L1 TPS": {"low": 0, "high": 1, "unit": "%"},
        }

    results = []
    for _, row in biomarkers_df.iterrows():
        marker = row["marker"]
        value = float(row["value"])
        ref = reference_ranges.get(marker, {})

        status = "normal"
        if ref:
            if value > ref.get("high", float("inf")):
                status = "HIGH"
            elif value < ref.get("low", float("-inf")):
                status = "LOW"

        results.append({
            "marker": marker,
            "value": value,
            "unit": row.get("unit", ref.get("unit", "")),
            "reference": f"{ref.get('low', '?')}-{ref.get('high', '?')}",
            "status": status,
        })

    report_df = pd.DataFrame(results)
    abnormal = report_df[report_df["status"] != "normal"]

    report = {
        "report_type": "Biomarker_Profile",
        "total_markers": len(report_df),
        "abnormal_count": len(abnormal),
        "results": report_df.to_dict("records"),
        "summary": (
            f"{len(abnormal)} of {len(report_df)} markers outside reference range"
            if len(abnormal) > 0
            else "All markers within reference range"
        ),
    }

    print(f"Biomarker profile: {len(report_df)} markers, "
          f"{len(abnormal)} abnormal")
    return report
```

## 3. ファーマコゲノミクスレポート

```python
def pharmacogenomics_report(genotypes, medications):
    """
    CPIC ガイドライン準拠のファーマコゲノミクスレポート。

    Parameters:
        genotypes: dict — {"CYP2D6": "*1/*4", "CYP2C19": "*1/*2", ...}
        medications: list — ["codeine", "clopidogrel", ...]
    """
    # CPIC phenotype マッピング (簡略)
    cpic_phenotypes = {
        "CYP2D6": {
            "*1/*1": "Normal Metabolizer",
            "*1/*4": "Intermediate Metabolizer",
            "*4/*4": "Poor Metabolizer",
            "*1/*2xN": "Ultrarapid Metabolizer",
        },
        "CYP2C19": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Intermediate Metabolizer",
            "*2/*2": "Poor Metabolizer",
            "*1/*17": "Rapid Metabolizer",
            "*17/*17": "Ultrarapid Metabolizer",
        },
    }

    # 推奨アクション (簡略)
    drug_gene_map = {
        "codeine": {"gene": "CYP2D6", "action": {
            "Poor Metabolizer": "AVOID — use alternative analgesic",
            "Ultrarapid Metabolizer": "AVOID — toxicity risk",
            "Intermediate Metabolizer": "Use with caution, consider alternative",
        }},
        "clopidogrel": {"gene": "CYP2C19", "action": {
            "Poor Metabolizer": "Use alternative antiplatelet (e.g., prasugrel)",
            "Intermediate Metabolizer": "Consider alternative antiplatelet",
        }},
    }

    recommendations = []
    for drug in medications:
        entry = drug_gene_map.get(drug, {})
        gene = entry.get("gene", "Unknown")
        genotype = genotypes.get(gene, "Unknown")
        phenotype_map = cpic_phenotypes.get(gene, {})
        phenotype = phenotype_map.get(genotype, "Indeterminate")

        action = entry.get("action", {}).get(phenotype, "Standard dosing")
        recommendations.append({
            "drug": drug,
            "gene": gene,
            "genotype": genotype,
            "phenotype": phenotype,
            "recommendation": action,
            "cpic_level": "A" if drug in drug_gene_map else "N/A",
        })

    report = {
        "report_type": "Pharmacogenomics",
        "genotypes_tested": genotypes,
        "medications_queried": medications,
        "recommendations": recommendations,
    }

    print(f"PGx report: {len(genotypes)} genes, {len(medications)} drugs, "
          f"{len(recommendations)} recommendations")
    return report
```

## 4. 構造化レポート出力 (LaTeX/HTML)

```python
def export_clinical_report(report, output_format="html",
                            output_path="reports/clinical_report"):
    """
    臨床レポートを LaTeX/HTML/FHIR JSON 形式で出力。

    Parameters:
        report: dict — SOAP, Biomarker, PGx レポート
        output_format: "html", "latex", "fhir_json"
        output_path: str — 出力先パス (拡張子なし)
    """
    import os
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    report_type = report.get("report_type", "Clinical")

    if output_format == "html":
        filepath = f"{output_path}.html"
        html_parts = [
            "<!DOCTYPE html><html><head>",
            f"<title>{report_type} Report</title>",
            "<style>body{font-family:Arial;margin:2em;}"
            "table{border-collapse:collapse;width:100%;}"
            "td,th{border:1px solid #ddd;padding:8px;}</style>",
            "</head><body>",
            f"<h1>{report_type} Report</h1>",
        ]

        if report_type == "SOAP_Note":
            for section in ["S", "O", "A", "P"]:
                html_parts.append(f"<h2>{section}</h2>")
                html_parts.append(f"<pre>{json.dumps(report.get(section, {}), indent=2, ensure_ascii=False)}</pre>")

        elif report_type == "Biomarker_Profile":
            html_parts.append("<table><tr><th>Marker</th><th>Value</th>"
                              "<th>Reference</th><th>Status</th></tr>")
            for r in report.get("results", []):
                status_color = "red" if r["status"] != "normal" else "green"
                html_parts.append(
                    f"<tr><td>{r['marker']}</td><td>{r['value']} {r['unit']}</td>"
                    f"<td>{r['reference']}</td>"
                    f"<td style='color:{status_color}'>{r['status']}</td></tr>"
                )
            html_parts.append("</table>")

        html_parts.append("</body></html>")
        with open(filepath, "w") as f:
            f.write("\n".join(html_parts))

    elif output_format == "fhir_json":
        filepath = f"{output_path}.fhir.json"
        fhir = {
            "resourceType": "DiagnosticReport",
            "status": "final",
            "category": [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                                       "code": "LAB"}]}],
            "code": {"text": report_type},
            "issued": report.get("generated_at", datetime.now().isoformat()),
            "result": [],
        }
        with open(filepath, "w") as f:
            json.dump(fhir, f, indent=2)

    elif output_format == "latex":
        filepath = f"{output_path}.tex"
        with open(filepath, "w") as f:
            f.write(f"\\documentclass{{article}}\n")
            f.write(f"\\title{{{report_type} Report}}\n")
            f.write("\\begin{document}\n\\maketitle\n")
            f.write(f"Report type: {report_type}\n")
            f.write("\\end{document}\n")

    print(f"Report exported: {filepath}")
    return filepath
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `reports/soap_note.json` | JSON |
| `reports/biomarker_profile.json` | JSON |
| `reports/pgx_report.json` | JSON |
| `reports/clinical_report.html` | HTML |
| `reports/clinical_report.tex` | LaTeX |
| `reports/clinical_report.fhir.json` | FHIR JSON |

### 利用可能ツール

> 本スキルは ToolUniverse ツールに直接依存しない。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| — | — | — |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-variant-interpretation` | バリアント解釈レポート |
| `scientific-variant-effect-prediction` | バリアント病原性スコア |
| `scientific-pharmacogenomics` | PGx ガイドライン |
| `scientific-precision-oncology` | 精密腫瘍学レポート |
| `scientific-disease-research` | 疾患情報統合 |

### 依存パッケージ

`pandas`, `json` (stdlib), `datetime` (stdlib)
