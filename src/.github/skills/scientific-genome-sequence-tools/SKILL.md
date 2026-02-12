---
name: scientific-genome-sequence-tools
description: |
  ゲノム配列解析総合スキル。Ensembl ゲノムブラウザ、dbSNP 変異データ、
  BLAST 相同性検索、NCBI Nucleotide 配列取得、GDC がんゲノミクスデータの
  統合パイプライン。
---

# Scientific Genome Sequence Tools

公的ゲノムデータベース (Ensembl, dbSNP, BLAST, NCBI, GDC) を横断した
配列検索・変異アノテーション・がんゲノミクスパイプラインを提供する。

## When to Use

- ゲノム配列・エクソン構造を Ensembl から取得するとき
- rsID から変異のアレル頻度を調べるとき
- BLAST で塩基/アミノ酸配列の相同性検索を行うとき
- NCBI Nucleotide から配列をフェッチするとき
- GDC がんゲノミクスデータ (体細胞変異, CNV, 発現) を取得するとき

---

## Quick Start

## 1. dbSNP 変異情報取得

```python
import requests
import pandas as pd


def get_dbsnp_variant(rsid):
    """
    dbSNP から rsID ベースの変異情報 (アレル頻度含む) を取得。

    Parameters:
        rsid: str — e.g. "rs7412"

    ToolUniverse:
        dbsnp_get_variant_by_rsid(rsid=rsid)
        dbsnp_get_frequencies(rsid=rsid)
        dbsnp_search_by_gene(gene_symbol=gene_symbol)
    """
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid.lstrip('rs')}"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    # Extract primary info
    info = {
        "rsid": f"rs{data.get('refsnp_id', '')}",
        "create_date": data.get("create_date", ""),
        "update_date": data.get("update_date", ""),
    }

    # Allele frequencies
    alleles = data.get("primary_snapshot_data", {}).get(
        "allele_annotations", []
    )
    freq_data = []
    for allele in alleles:
        for freq_entry in allele.get("frequency", []):
            freq_data.append({
                "study": freq_entry.get("study_name", ""),
                "allele": freq_entry.get("allele", ""),
                "count": freq_entry.get("allele_count", 0),
                "total": freq_entry.get("total_count", 0),
            })

    df_freq = pd.DataFrame(freq_data)
    print(f"dbSNP {info['rsid']}: {len(df_freq)} frequency entries")
    return info, df_freq
```

## 2. BLAST 相同性検索

```python
import time


def blast_search(sequence, program="blastn", database="nt", max_hits=10):
    """
    NCBI BLAST REST API で相同性検索。

    Parameters:
        sequence: str — query sequence (nucleotide or protein)
        program: str — "blastn", "blastp", "blastx", "tblastn"
        database: str — "nt", "nr", "refseq_rna", etc.

    ToolUniverse:
        BLAST_nucleotide_search(sequence=sequence, database=database)
        BLAST_protein_search(sequence=sequence, database=database)
    """
    put_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": sequence,
        "FORMAT_TYPE": "JSON2",
        "HITLIST_SIZE": max_hits,
    }
    resp = requests.post(put_url, data=params)
    resp.raise_for_status()

    # Extract RID
    import re
    rid_match = re.search(r"RID = (\S+)", resp.text)
    if not rid_match:
        raise ValueError("BLAST RID not found")
    rid = rid_match.group(1)
    print(f"BLAST submitted: RID={rid}")

    # Poll for results
    for _ in range(60):
        time.sleep(10)
        check = requests.get(put_url, params={
            "CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"
        })
        if "Status=WAITING" not in check.text:
            break

    return check.json() if check.headers.get(
        "Content-Type", ""
    ).startswith("application/json") else check.text
```

## 3. NCBI Nucleotide 配列フェッチ

```python
def fetch_ncbi_sequence(accession, rettype="fasta"):
    """
    NCBI Nucleotide (E-utilities) から配列を取得。

    Parameters:
        accession: str — NCBI accession (e.g., "NM_000546.6")
        rettype: str — "fasta", "gb", "gbwithparts"

    ToolUniverse:
        NCBI_search_nucleotide(query=query)
        NCBI_fetch_accessions(accessions=accessions)
        NCBI_get_sequence(accession=accession)
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": accession,
        "rettype": rettype,
        "retmode": "text",
    }
    resp = requests.get(url, params=params)
    resp.raise_for_status()

    print(f"NCBI Nucleotide '{accession}': {len(resp.text)} chars ({rettype})")
    return resp.text
```

## 4. GDC がんゲノミクスデータ

```python
def get_gdc_mutations(gene_symbol, project_id=None):
    """
    NCI GDC (Genomic Data Commons) から体細胞変異データを取得。

    Parameters:
        gene_symbol: str — e.g. "TP53"
        project_id: str | None — e.g. "TCGA-BRCA"

    ToolUniverse:
        GDC_get_ssm_by_gene(gene_symbol=gene_symbol)
        GDC_get_mutation_frequency(project_id=project_id)
        GDC_get_gene_expression(gene_id=gene_id, project_id=project_id)
        GDC_get_cnv_data(gene_id=gene_id)
        GDC_list_projects()
        GDC_search_cases(filters=filters)
        GDC_list_files(filters=filters)
    """
    url = "https://api.gdc.cancer.gov/ssms"
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {
                "field": "consequence.transcript.gene.symbol",
                "value": [gene_symbol],
            }},
        ],
    }
    if project_id:
        filters["content"].append({
            "op": "in",
            "content": {
                "field": "cases.project.project_id",
                "value": [project_id],
            },
        })

    import json
    params = {
        "filters": json.dumps(filters),
        "fields": ("ssm_id,consequence.transcript.gene.symbol,"
                    "consequence.transcript.aa_change,"
                    "consequence.transcript.consequence_type,"
                    "genomic_dna_change"),
        "size": 100,
        "format": "json",
    }
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    hits = resp.json().get("data", {}).get("hits", [])

    results = []
    for hit in hits:
        for csq in hit.get("consequence", []):
            tx = csq.get("transcript", {})
            results.append({
                "ssm_id": hit.get("ssm_id", ""),
                "gene": tx.get("gene", {}).get("symbol", ""),
                "aa_change": tx.get("aa_change", ""),
                "consequence_type": tx.get("consequence_type", ""),
                "genomic_dna_change": hit.get("genomic_dna_change", ""),
            })

    df = pd.DataFrame(results)
    print(f"GDC SSMs '{gene_symbol}'"
          f"{f' ({project_id})' if project_id else ''}: {len(df)} mutations")
    return df
```

## 5. 統合ゲノム変異パイプライン

```python
def integrated_variant_pipeline(rsid, gene_symbol=None):
    """
    dbSNP + GDC を統合したゲノム変異解析パイプライン。

    ToolUniverse (横断):
        dbsnp_get_variant_by_rsid(rsid) → GDC_get_ssm_by_gene(gene_symbol)
    """
    pipeline_result = {"rsid": rsid}

    # Step 1: dbSNP
    info, freq_df = get_dbsnp_variant(rsid)
    pipeline_result["dbsnp"] = info

    # Step 2: GDC somatic mutations (if gene provided)
    if gene_symbol:
        gdc_df = get_gdc_mutations(gene_symbol)
        pipeline_result["gdc_mutation_count"] = len(gdc_df)
        pipeline_result["gdc_top_consequences"] = (
            gdc_df["consequence_type"].value_counts().head(5).to_dict()
            if not gdc_df.empty else {}
        )

    print(f"Integrated variant: {rsid}"
          f" | GDC={pipeline_result.get('gdc_mutation_count', 'N/A')}")
    return pipeline_result
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/dbsnp_variant.json` | JSON |
| `results/dbsnp_frequencies.csv` | CSV |
| `results/blast_results.json` | JSON |
| `results/ncbi_sequence.fasta` | FASTA |
| `results/gdc_mutations.csv` | CSV |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| dbSNP | `dbsnp_get_variant_by_rsid` | rsID 変異情報 |
| dbSNP | `dbsnp_get_frequencies` | アレル頻度 |
| dbSNP | `dbsnp_search_by_gene` | 遺伝子→変異 |
| BLAST | `BLAST_nucleotide_search` | 核酸相同性検索 |
| BLAST | `BLAST_protein_search` | タンパク質相同性検索 |
| NCBI | `NCBI_search_nucleotide` | 配列検索 |
| NCBI | `NCBI_fetch_accessions` | アクセッション取得 |
| NCBI | `NCBI_get_sequence` | 配列フェッチ |
| GDC | `GDC_get_ssm_by_gene` | 体細胞変異 |
| GDC | `GDC_get_mutation_frequency` | 変異頻度 |
| GDC | `GDC_get_gene_expression` | 発現データ |
| GDC | `GDC_get_cnv_data` | CNV データ |
| GDC | `GDC_list_projects` | プロジェクト一覧 |
| GDC | `GDC_search_cases` | 症例検索 |
| GDC | `GDC_list_files` | ファイル一覧 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-variant-interpretation` | 変異アノテーション |
| `scientific-population-genetics` | 集団遺伝学 |
| `scientific-cancer-genomics` | がんゲノミクス |
| `scientific-rare-disease-genetics` | 希少疾患遺伝学 |
| `scientific-biothings-idmapping` | ID マッピング |

### 依存パッケージ

`requests`, `pandas`
