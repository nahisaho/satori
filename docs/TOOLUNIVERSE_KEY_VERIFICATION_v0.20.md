# ToolUniverse Key 検証レポート — SATORI v0.20.0 計画用

> **生成日**: 2025-07  
> **リポジトリ**: mims-harvard/ToolUniverse (main branch)  
> **参照ファイル**: `src/tooluniverse/default_config.py`, `src/tooluniverse/tools/__init__.py`, `docs/tools/*.rst`  
> **検証方法**: GitHub リポジトリ直接検索 + raw ファイル取得

---

## 1. 検証サマリー

| 分類 | 件数 |
|------|------|
| 候補キー総数 | 19 |
| **EXISTS** (完全一致) | 12 |
| **EXISTS (キー名相違)** | 5 |
| **NOT FOUND** | 2 |

---

## 2. Priority A — 検証結果 (6 件)

| # | 候補キー | 実際の TU Key | ステータス | ツール数 | 主要ツール |
|---|---------|--------------|-----------|---------|-----------|
| 1 | `cosmic` | `"cosmic"` | **EXISTS** | 2 | `COSMIC_get_mutations_by_gene`, `COSMIC_search_mutations` |
| 2 | `cbioportal` | `"cbioportal"` | **EXISTS** | ~12 | `cBioPortal_get_cancer_studies`, `cBioPortal_get_mutations`, etc. |
| 3 | `oncokb` | `"oncokb"` | **EXISTS** | ~5 | `OncoKB_annotate_mutations`, `OncoKB_get_gene_info`, etc. |
| 4 | `iedb` | `"iedb_tools"` | **EXISTS (キー名相違)** | ~8 | `iedb_search_epitopes`, `iedb_search_bcell`, `iedb_get_epitope_antigens`, `iedb_get_epitope_mhc`, `iedb_search_mhc`, `iedb_search_antigens`, `iedb_search_references`, `iedb_get_epitope_references` |
| 5 | `card` | — | **NOT FOUND** | 0 | CARD (Comprehensive Antibiotic Resistance Database) はTUに未実装 |
| 6 | `mgnify` | `"mgnify"` | **EXISTS** | 2 | `MGnify_search_studies`, `MGnify_list_analyses` |

### Priority A 注意事項

- **`iedb`**: 正しいキー名は `"iedb_tools"` (not `"iedb"`)。スキル定義で `tu_keys` に指定する際は `iedb_tools` を使用すること。
- **`card`**: `default_config.py` に一切記載なし。CARD API ツールは ToolUniverse リポジトリ全体で実装が確認できず。v0.20.0 では **スキップ推奨**。将来 TU 側で実装された場合に追加を検討。

---

## 3. Priority B — 検証結果 (13 件)

| # | 候補キー | 実際の TU Key | ステータス | ツール数 | 主要ツール |
|---|---------|--------------|-----------|---------|-----------|
| 7 | `hgnc` | — | **NOT FOUND** | 0 | HGNC (HUGO Gene Nomenclature) はTUに未実装 |
| 8 | `hmdb` | `"hmdb"` | **EXISTS** | 3 | `HMDB_search`, `HMDB_get_metabolite`, `HMDB_get_diseases` |
| 9 | `chembl` | `"ChEMBL"` | **EXISTS (キー名相違)** | ~20+ | `ChEMBL_search_molecule`, `ChEMBL_get_activity`, `ChEMBL_get_assay`, `ChEMBL_get_target`, `ChEMBL_search_tissue`, etc. |
| 10 | `pubchem` | `"pubchem"` | **EXISTS** | ~12+ | `PubChem_get_compound`, `PubChem_search_by_name`, `PubChem_get_assay`, etc. |
| 11 | `ncbi` / `entrez` | `"ncbi_nucleotide"` | **EXISTS (キー名相違)** | ~5+ | `NCBI_search_nucleotide`, `NCBI_get_sequence`, `NCBI_fetch_accessions`, etc. |
| 12 | `bindingdb` | `"bindingdb"` | **EXISTS** | 4 | `BindingDB_get_ligands_by_pdb`, `BindingDB_get_ligands_by_uniprot`, `BindingDB_get_ligands_by_uniprots`, `BindingDB_get_targets_by_compound` |
| 13 | `brenda` | `"brenda"` | **EXISTS** | ~2+ | `BRENDA_get_inhibitors`, etc. |
| 14 | `dgidb` | `"dgidb"` | **EXISTS** | 4 | `DGIdb_get_drug_gene_interactions`, `DGIdb_get_drug_info`, `DGIdb_get_gene_druggability`, `DGIdb_get_gene_info` |
| 15 | `string_db` / `string` | `"ppi"` | **EXISTS (キー名相違)** | 7+ | `STRING_get_protein_interactions`, `STRING_get_interaction_partners`, `STRING_get_network`, `STRING_map_identifiers`, `STRING_ppi_enrichment`, `STRING_functional_enrichment`, `BioGRID_get_interactions` |
| 16 | `biomodels` | `"biomodels_tools"` | **EXISTS (キー名相違)** | 4 | `BioModels_search_parameters`, `BioModels_get_model`, `BioModels_list_files`, `BioModels_download_model` |
| 17 | `monarch` | `"monarch"` | **EXISTS** | ~3 | `Monarch_search_gene`, `Monarch_get_gene_diseases`, `Monarch_get_gene_phenotypes` |
| 18 | `cellosaurus` | `"cellosaurus"` | **EXISTS** | ~3 | `Cellosaurus_search`, `Cellosaurus_get_cell_line`, `Cellosaurus_query_converter` |
| 19 | `gtopdb` | `"gtopdb"` | **EXISTS** | 8 | `GtoPdb_get_target`, `GtoPdb_get_targets`, `GtoPdb_get_ligand`, `GtoPdb_get_disease`, `GtoPdb_get_target_interactions`, `GtoPdb_list_diseases`, `GtoPdb_list_ligands`, `GtoPdb_search_interactions` |

### Priority B 注意事項

- **`ChEMBL`**: 大文字始まり。`tu_keys` には `ChEMBL` と正確に指定必要。
- **`ppi`**: STRING + BioGRID の統合キー。さらに `"biogrid"` が別キーとして独立存在 (`biogrid_tools.json`)。STRING 単体を使いたい場合も `ppi` で可。
- **`ncbi_nucleotide`**: NCBI ヌクレオチド配列データベース専用。PubMed は `"pubmed"` 、ClinVar は `"clinvar"` 、dbSNP は `"dbsnp"` でそれぞれ別キー。汎用 NCBI Entrez キーは存在しない。
- **`biomodels_tools`**: `"biomodels"` ではなく `"biomodels_tools"` が正しい。
- **`hgnc`**: TU 全体で未実装。遺伝子命名は UniProt/Ensembl/NCBI の各ツール内で間接的にサポート。

---

## 4. キー名相違マッピング一覧

スキル定義の `tu_keys` で正しいキー名を使用するための参照表:

| 想定キー名 | 正しい TU Key | 備考 |
|-----------|--------------|------|
| `iedb` | `iedb_tools` | アンダースコア付き |
| `chembl` | `ChEMBL` | 大文字始まり |
| `ncbi` / `entrez` | `ncbi_nucleotide` | ヌクレオチド専用 |
| `string_db` / `string` | `ppi` | STRING + BioGRID 統合 |
| `biomodels` | `biomodels_tools` | `_tools` サフィックス |

---

## 5. 未使用 TU キー — v0.20.0 向け新規発見

以下は `default_config.py` に存在するが、SATORI v0.19.0 のスキルマッピングでまだ十分活用されていない可能性のあるキー:

### Tier 1: 高インパクト (科学的重要度 ★★★)

| TU Key | データベース名 | ツール数 | 推奨カテゴリ | 根拠 |
|--------|--------------|---------|------------|------|
| `"stitch"` | STITCH Chemical-Protein Interactions | ~3+ | J (創薬) | STRING の化合物版。Drug-target discovery に有用 |
| `"pharos"` | Pharos/TCRD (NIH IDG) | ~3+ | J (創薬) | NIH Illuminating the Druggable Genome。未研究タンパク質の創薬ターゲット探索 |
| `"biothings"` | BioThings (MyGene, MyVariant, MyChem) | ~5+ | F (生命科学) | 統合的遺伝子/バリアント/化合物 API。HGNC の代替としても機能 |
| `"cellxgene_census"` | CELLxGENE Census | 7 | F (生命科学) | 単一細胞 RNA-seq アトラス。CELLxGENE_download_h5ad ~ CELLxGENE_get_presence_matrix |
| `"expression_atlas"` | Expression Atlas (EBI GXA) | ~3+ | F (生命科学) | 遺伝子発現アトラス。baseline + differential expression |
| `"fda_pharmacogenomic_biomarkers"` | FDA PGx Biomarkers | ~2+ | P (薬理ゲノミクス) | FDA 承認の薬理ゲノミクスバイオマーカー |
| `"gdc"` | GDC (Genomic Data Commons) | ~7 | Q (腫瘍学) | NCI ゲノミクスデータ。GDC_get_gene_expression, GDC_get_mutation_frequency 等 |

### Tier 2: 中インパクト (科学的重要度 ★★☆)

| TU Key | データベース名 | ツール数 | 推奨カテゴリ | 根拠 |
|--------|--------------|---------|------------|------|
| `"chipatlas"` | ChIP-Atlas | ~4 | F (生命科学) | エピゲノミクス/クロマチン解析 |
| `"fourdn"` | 4DN Data Portal | ~4 | F (生命科学) | 3D ゲノム構造 |
| `"metabolomics_workbench"` | Metabolomics Workbench | ~3+ | F (生命科学) | メタボロミクスデータリポジトリ |
| `"bigg_models"` | BiGG Models | ~3+ | F (生命科学) | ゲノムスケール代謝モデル |
| `"rfam"` | Rfam | ~2+ | F (生命科学) | RNA ファミリーデータベース |
| `"enamine"` | Enamine Make-on-Demand | ~3+ | G (化学) | 化合物ベンダー/仮想スクリーニング |
| `"emolecules"` | eMolecules | ~2+ | G (化学) | 化合物ベンダーアグリゲーター |
| `"nvidia_nim"` | NVIDIA NIM Healthcare | ~3+ | K (構造生物) | 構造予測、分子ドッキング、ゲノミクス AI |
| `"therasabdab"` | Thera-SAbDab | 3 | K (構造生物) | 治療用抗体構造データベース |
| `"impc"` | IMPC | ~4 | Q (疾患研究) | マウス表現型コンソーシアム。遺伝子機能解析 |

### Tier 3: ニッチだが有用 (★☆☆)

| TU Key | データベース名 | ツール数 | 備考 |
|--------|--------------|---------|------|
| `"sabdab"` | SAbDab | ~3+ | 構造抗体データベース |
| `"imgt"` | IMGT | ~3 | 免疫遺伝学 |
| `"gpcrdb"` | GPCRdb | ~5 | GPCR 受容体データベース |
| `"swissdock"` | SwissDock | 3 | 分子ドッキング |
| `"proteinsplus"` | ProteinsPlus | ~3+ | タンパク質-リガンドドッキング |
| `"ols"` | OLS (Ontology Lookup) | ~2+ | オントロジー検索 |
| `"hca_tools"` | Human Cell Atlas | 2 | ヒト細胞アトラス |

---

## 6. v0.20.0 推奨優先順位

### Phase 12 実装推奨 — Track A: 既存スキルへの TU Key 追加

| スキル (既存) | 追加 TU Key | ツール数増 | 根拠 |
|-------------|------------|----------|------|
| `scientific-immunology-databases` | `iedb_tools` | +8 | エピトープ検索・MHC結合予測。免疫学研究の標準 |
| `scientific-protein-interaction-network` | `ppi` | +7 | STRING/BioGRID による PPI ネットワーク強化 |
| `scientific-cancer-genomics` | `cosmic`, `cbioportal`, `oncokb` | +19 | 腫瘍ゲノミクス3大データベース統合 |
| `scientific-enzyme-kinetics` | `brenda` | +2 | 酵素速度論パラメータ |
| `scientific-metabolomics-databases` | `hmdb` | +3 | ヒトメタボロームデータベース |
| `scientific-pharmacology-targets` | `gtopdb`, `dgidb` | +12 | 薬理ターゲット + Drug-Gene 相互作用 |
| `scientific-compound-screening` | `bindingdb` | +4 | タンパク質-リガンド結合アフィニティ |

### Phase 12 実装推奨 — Track B: 新規スキル

| # | 提案スキル名 | Cat | TU Keys | ツール数 | 根拠 |
|---|------------|-----|---------|---------|------|
| 149 | `scientific-metagenomics-microbiome` | F | `mgnify` | ~2 | メタゲノミクス/マイクロバイオーム研究 |
| 150 | `scientific-ncbi-nucleotide` | F | `ncbi_nucleotide` | ~5 | ヌクレオチド配列検索・取得 |
| 151 | `scientific-cell-line-databases` | F | `cellosaurus` | ~3 | 細胞株データベース検索 |
| 152 | `scientific-gene-disease-monarch` | L | `monarch` | ~3 | 遺伝子-疾患-表現型の統合リソース |
| 153 | `scientific-systems-biology-models` | F | `biomodels_tools` | 4 | 生体系計算モデル |
| 154 | `scientific-single-cell-atlas` | F | `cellxgene_census` | 7 | 単一細胞 RNA-seq アトラス |
| 155 | `scientific-pharmacogenomics-biomarkers` | P | `fda_pharmacogenomic_biomarkers`, `pharmgkb` | ~5+ | 薬理ゲノミクスバイオマーカー |
| 156 | `scientific-druggable-genome` | J | `pharos`, `stitch` | ~6+ | NIH IDG + 化合物-タンパク質相互作用 |

### 実装しないキー（v0.20.0）

| キー | 理由 |
|------|------|
| `card` | TU に未実装。将来追加される可能性あり |
| `hgnc` | TU に未実装。`biothings` (MyGene) が代替手段 |
| `ChEMBL`, `pubchem` | 既に v0.11〜v0.18 で活用済み |

---

## 7. ツール総数概算

| 区分 | ツール数 |
|------|---------|
| Priority A (EXISTS) | ~29 tools (5 keys) |
| Priority B (EXISTS) | ~72+ tools (11 keys) |
| 新規発見 Tier 1 | ~30+ tools (7 keys) |
| 新規発見 Tier 2 | ~32+ tools (10 keys) |
| **合計利用可能** | **~163+ tools** |

---

*END OF REPORT*
