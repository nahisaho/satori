---
name: scientific-data-submission
description: |
  科学データ登録・アーカイブスキル。GenBank/SRA 配列登録・
  ENA 配列アーカイブ・GEO 発現データ登録・BioProject/BioSample
  メタデータ管理・FAIR 原則準拠データ共有。
---

# Scientific Data Submission

GenBank / SRA / ENA / GEO / BioProject を活用した科学データの
登録・アーカイブパイプラインを提供する。FAIR 原則に準拠した
配列データ・発現データ・メタデータの公開準備。

## When to Use

- 配列データを GenBank/DDBJ/ENA に登録するとき
- RNA-seq/WGS データを SRA にアーカイブするとき
- GEO にマイクロアレイ/RNA-seq 発現データを登録するとき
- BioProject/BioSample でメタデータを構造化するとき
- 論文投稿時にデータアクセッション番号が必要なとき
- FAIR 原則 (Findable, Accessible, Interoperable, Reusable) に準拠するとき

---

## Quick Start

## 1. BioProject/BioSample メタデータ作成

```python
import json
import pandas as pd
from pathlib import Path
from datetime import date


def create_bioproject_metadata(title, description, organism,
                                data_type="Genome Sequencing",
                                relevance="Medical"):
    """
    BioProject メタデータ XML/JSON 生成。

    Parameters:
        title: str — プロジェクトタイトル
        description: str — プロジェクト説明
        organism: str — 生物種名
        data_type: str — データ種別
        relevance: str — 関連分野
    """
    bioproject = {
        "Project": {
            "ProjectID": {"ArchiveID": {"accession": "PRJNA_PENDING"}},
            "Descriptor": {
                "Title": title,
                "Description": description,
                "Relevance": relevance,
            },
            "ProjectType": {
                "ProjectTypeSubmission": {
                    "Target": {
                        "Organism": {"OrganismName": organism},
                    },
                    "Method": {"MethodType": data_type},
                    "Objectives": {"Data": {"DataType": data_type}},
                }
            },
        }
    }

    print(f"BioProject metadata created:")
    print(f"  Title: {title}")
    print(f"  Organism: {organism}")
    print(f"  Data type: {data_type}")
    return bioproject


def create_biosample_table(samples, organism, package="Generic"):
    """
    BioSample TSV テンプレート生成。

    Parameters:
        samples: list[dict] — サンプル情報
        organism: str — 生物種名
        package: str — BioSample パッケージ
    """
    required_fields = [
        "sample_name", "organism", "collection_date",
        "geo_loc_name", "tissue", "description",
    ]

    rows = []
    for s in samples:
        row = {
            "sample_name": s.get("name", ""),
            "organism": organism,
            "collection_date": s.get("date", str(date.today())),
            "geo_loc_name": s.get("location", "not collected"),
            "tissue": s.get("tissue", "not applicable"),
            "description": s.get("description", ""),
        }
        row.update({k: v for k, v in s.items()
                     if k not in ["name", "date", "location"]})
        rows.append(row)

    df = pd.DataFrame(rows)
    print(f"BioSample table: {len(df)} samples, package='{package}'")
    return df
```

## 2. GenBank 配列登録準備

```python
def prepare_genbank_submission(sequences, annotations, output_dir="submission"):
    """
    GenBank 配列登録用 .sqn ファイル準備。

    Parameters:
        sequences: dict — {seq_id: sequence_string}
        annotations: dict — {seq_id: {gene, product, organism, ...}}
        output_dir: str — 出力ディレクトリ
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # FASTA 生成
    fasta_path = output_dir / "sequences.fsa"
    with open(fasta_path, "w") as f:
        for seq_id, seq in sequences.items():
            ann = annotations.get(seq_id, {})
            organism = ann.get("organism", "Unknown organism")
            f.write(f">{seq_id} [organism={organism}]\n")
            # 80文字折返し
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

    # Feature Table 生成
    tbl_path = output_dir / "sequences.tbl"
    with open(tbl_path, "w") as f:
        for seq_id, ann in annotations.items():
            f.write(f">Feature {seq_id}\n")
            if "gene" in ann:
                f.write(f"1\t{len(sequences[seq_id])}\tgene\n")
                f.write(f"\t\t\tgene\t{ann['gene']}\n")
            if "product" in ann:
                f.write(f"1\t{len(sequences[seq_id])}\tCDS\n")
                f.write(f"\t\t\tproduct\t{ann['product']}\n")

    # Template 生成
    template = {
        "source": {
            "organism": list(annotations.values())[0].get("organism", ""),
            "mol_type": "genomic DNA",
        },
        "submitter": {
            "name": "AutoSubmission",
        },
    }

    template_path = output_dir / "template.json"
    with open(template_path, "w") as f:
        json.dump(template, f, indent=2)

    print(f"GenBank submission prepared: {len(sequences)} sequences")
    print(f"  FASTA: {fasta_path}")
    print(f"  Feature Table: {tbl_path}")
    return {"fasta": str(fasta_path), "tbl": str(tbl_path)}
```

## 3. SRA メタデータ & アップロード

```python
def prepare_sra_metadata(samples, library_strategy="WGS",
                          library_source="GENOMIC",
                          platform="ILLUMINA",
                          instrument_model="Illumina NovaSeq 6000"):
    """
    SRA メタデータ TSV 生成。

    Parameters:
        samples: list[dict] — {biosample, title, file_r1, file_r2}
        library_strategy: str — WGS/RNA-Seq/AMPLICON/etc.
        library_source: str — GENOMIC/TRANSCRIPTOMIC/etc.
        platform: str — ILLUMINA/OXFORD_NANOPORE/etc.
    """
    rows = []
    for s in samples:
        rows.append({
            "biosample_accession": s.get("biosample", "SAMN_PENDING"),
            "library_ID": s.get("library_id", s.get("title", "")),
            "title": s.get("title", ""),
            "library_strategy": library_strategy,
            "library_source": library_source,
            "library_selection": s.get("selection", "RANDOM"),
            "library_layout": "paired" if s.get("file_r2") else "single",
            "platform": platform,
            "instrument_model": instrument_model,
            "filetype": "fastq",
            "filename": s.get("file_r1", ""),
            "filename2": s.get("file_r2", ""),
        })

    df = pd.DataFrame(rows)
    print(f"SRA metadata: {len(df)} runs, strategy={library_strategy}")
    return df


def sra_upload_ascp(files, destination, ascp_key=None):
    """
    Aspera (ascp) による SRA データ高速アップロード。

    Parameters:
        files: list — アップロードファイルリスト
        destination: str — SRA アップロード先
        ascp_key: str — Aspera SSH キーパス
    """
    import subprocess

    if ascp_key is None:
        ascp_key = Path.home() / ".aspera/connect/etc/asperaweb_id_dsa.openssh"

    for f in files:
        cmd = [
            "ascp", "-i", str(ascp_key),
            "-QT", "-l", "300m", "-k", "1",
            str(f), destination,
        ]
        print(f"Uploading: {f}")
        subprocess.run(cmd, check=True)

    print(f"SRA upload complete: {len(files)} files")
```

## 4. GEO 発現データ登録

```python
def prepare_geo_submission(expression_matrix, sample_metadata,
                            platform="GPL16791", output_dir="geo_submission"):
    """
    GEO SOFT 形式サブミッション準備。

    Parameters:
        expression_matrix: pd.DataFrame — 遺伝子 × サンプルマトリクス
        sample_metadata: pd.DataFrame — サンプルメタデータ
        platform: str — GEO プラットフォーム ID
        output_dir: str — 出力ディレクトリ
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # SOFT テンプレート
    soft_path = output_dir / "submission.soft"
    with open(soft_path, "w") as f:
        # Series section
        f.write("^SERIES\n")
        f.write("!Series_title = \n")
        f.write("!Series_summary = \n")
        f.write(f"!Series_platform_id = {platform}\n")

        # Sample sections
        for col in expression_matrix.columns:
            meta = sample_metadata[sample_metadata["sample_id"] == col]
            f.write(f"\n^SAMPLE = {col}\n")
            f.write(f"!Sample_title = {col}\n")
            if len(meta) > 0:
                for key, val in meta.iloc[0].items():
                    if key != "sample_id":
                        f.write(f"!Sample_characteristics_ch1 = {key}: {val}\n")

    # Matrix file
    matrix_path = output_dir / "expression_matrix.txt"
    expression_matrix.to_csv(matrix_path, sep="\t")

    # Raw data files list
    raw_files_path = output_dir / "raw_files.txt"
    with open(raw_files_path, "w") as f:
        for col in expression_matrix.columns:
            f.write(f"{col}.fastq.gz\n")

    print(f"GEO submission: {expression_matrix.shape[1]} samples, "
          f"{expression_matrix.shape[0]} genes")
    return {"soft": str(soft_path), "matrix": str(matrix_path)}
```

## 5. FAIR データ検証チェックリスト

```python
def fair_checklist(submission_package):
    """
    FAIR 原則準拠チェックリスト。

    Parameters:
        submission_package: dict — 登録パッケージ情報
    """
    checks = {
        "Findable": {
            "F1_persistent_id": bool(submission_package.get("accession")),
            "F2_metadata_rich": bool(submission_package.get("metadata")),
            "F3_id_in_metadata": True,
            "F4_searchable_registry": bool(submission_package.get("repository")),
        },
        "Accessible": {
            "A1_retrievable_protocol": bool(submission_package.get("access_url")),
            "A1_1_open_protocol": True,
            "A2_metadata_persists": True,
        },
        "Interoperable": {
            "I1_formal_language": bool(submission_package.get("format")),
            "I2_fair_vocabularies": bool(submission_package.get("ontology_terms")),
            "I3_qualified_references": bool(submission_package.get("references")),
        },
        "Reusable": {
            "R1_usage_license": bool(submission_package.get("license")),
            "R1_1_community_standards": bool(submission_package.get("standard")),
            "R1_2_provenance": bool(submission_package.get("methods")),
        },
    }

    total = 0
    passed = 0
    for principle, items in checks.items():
        for check, status in items.items():
            total += 1
            if status:
                passed += 1

    score = passed / total * 100 if total > 0 else 0
    print(f"FAIR checklist: {passed}/{total} ({score:.0f}%)")
    for principle, items in checks.items():
        n_pass = sum(items.values())
        print(f"  {principle}: {n_pass}/{len(items)}")

    return checks
```

---

## パイプライン統合

```
bioinformatics → data-submission → literature-search
  (解析完了)       (データ登録)       (論文投稿時)
       │                │                  ↓
lab-data-management ───┘           academic-writing
  (Benchling/OMERO)    │            (論文執筆)
                       ↓
                 ebi-databases
                 (ENA/BioStudies 連携)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `submission/sequences.fsa` | GenBank 登録用 FASTA | → bioinformatics |
| `submission/sra_metadata.tsv` | SRA メタデータ | → ebi-databases |
| `geo_submission/submission.soft` | GEO SOFT テンプレート | → gene-expression |
| `submission/fair_report.json` | FAIR チェックリスト結果 | → academic-writing |
