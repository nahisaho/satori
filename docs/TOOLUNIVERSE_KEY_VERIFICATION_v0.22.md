# ToolUniverse Key Verification Report — SATORI v0.22.0 Phase 14

**Date**: 2025-07-17  
**Repository**: [mims-harvard/ToolUniverse](https://github.com/mims-harvard/ToolUniverse)  
**Total TU Tools**: 1,265+ (per `__init__.py` header)  
**Source Files Analyzed**: `default_config.py` (lines 0–451), `_lazy_registry_static.py`, `tools/__init__.py` (`__all__`), individual `*_tool.py` files  
**Verification Method**: GitHub リポジトリ直接検索 × 9 回 + ローカル SKILL.md frontmatter 読取り × 16 ファイル

---

## 1. PRIORITY A — 検証結果 (4 件)

| # | 候補キー | ステータス | ツール数 | Config 行 |
|---|---------|-----------|---------|-----------|
| 1 | `chipatlas` | ✅ **VERIFIED** | 4 | `default_config.py:268` |
| 2 | `metabolomics_workbench` | ✅ **VERIFIED** | 6 | `default_config.py:253` |
| 3 | `card` | ❌ **NOT FOUND** | 0 | — |
| 4 | `gbif` | ✅ **VERIFIED** | 2 | `default_config.py:224` |

### ✅ `chipatlas` — ChIP-Atlas エピゲノミクス

- **Config**: `"chipatlas": os.path.join(current_dir, "data", "chipatlas_tools.json")`
- **Implementation**: `chipatlas_tool.py` (class `ChIPAtlasTool`)
- **Registry**: `_lazy_registry_static.py` → `"ChIPAtlasTool": "chipatlas_tool"`
- **Tools (4)**:

| Tool Function | Description |
|---|---|
| `ChIPAtlas_enrichment_analysis` | Enrichment analysis for epigenomic datasets |
| `ChIPAtlas_get_experiments` | Get experiment metadata (433K+ experiments) |
| `ChIPAtlas_get_peak_data` | Get peak data for specific experiments |
| `ChIPAtlas_search_datasets` | Search ChIP-Atlas datasets by condition |

- **SATORI 割当先候補**: `scientific-epigenomics-chromatin` (現在 tu_tools なし)

### ✅ `metabolomics_workbench` — Metabolomics Workbench

- **Config**: `"metabolomics_workbench": os.path.join(current_dir, "data", "metabolomics_workbench_tools.json")`
- **Implementation**: `metabolomics_workbench_tool.py` (class `MetabolomicsWorkbenchTool`)
- **Registry**: `"MetabolomicsWorkbenchTool": "metabolomics_workbench_tool"`
- **Tools (6)**:

| Tool Function | Description |
|---|---|
| `MetabolomicsWorkbench_get_study` | Get study metadata & results |
| `MetabolomicsWorkbench_get_refmet_info` | RefMet standardized nomenclature |
| `MetabolomicsWorkbench_search_by_exact_mass` | Exact mass search |
| `MetabolomicsWorkbench_search_by_mz` | m/z search for metabolite identification |
| `MetabolomicsWorkbench_search_compound_by_name` | Compound name search |
| `MetabolomicsWorkbench_get_compound_by_pubchem_cid` | PubChem CID lookup |

- **SATORI 割当先候補**: `scientific-metabolomics` (現在 `hmdb` のみ) または `scientific-metabolomics-databases` (現在 `metacyc` のみ)

### ❌ `card` — NOT FOUND

- **検索**: 2回の GitHub 検索 (`"CARD antimicrobial resistance database AMR tools"` + `"CARD antibiotic resistance"`) で未検出
- **v0.20 レポートとの一貫性**: v0.20 でも NOT FOUND (Priority A #5)
- **`default_config.py` 全キーリストに不在**
- **結論**: CARD (Comprehensive Antibiotic Resistance Database) は **ToolUniverse に未実装**。v0.22.0 でもスキップ。
- **代替**: AMR 関連解析は `ncbi_nucleotide` (GenBank 耐性遺伝子配列) + ローカル ABRicate/AMRFinderPlus で対応

### ✅ `gbif` — Global Biodiversity Information Facility

- **Config**: `"gbif": os.path.join(current_dir, "data", "gbif_tools.json")`
- **Implementation**: `gbif_tool.py` (class `GBIFTool`)
- **Tools (2)**:

| Tool Function | Description |
|---|---|
| `GBIF_search_species` | Species taxonomy search |
| `GBIF_search_occurrences` | Biodiversity occurrence records |

- **SATORI 現状**: **既に `scientific-marine-ecology` に割当済** (`obis`, `worms`, `gbif`)
- **`scientific-environmental-ecology` への重複割当**: 可能（海洋+陸域の両方で利用）

---

## 2. DOMAIN-SPECIFIC 新発見 — 4 ドメイン探索結果

### 2-A. 臨床・疫学 (Clinical / Epidemiology)

| TU Key | Tools | Implementation | SATORI 現状 | 推奨先スキル |
|--------|-------|----------------|-------------|-------------|
| **`nhanes`** | 2 | `NHANESTool` class, `nhanes_tool.py` | ❌ 未割当 | `scientific-public-health-data` |
| **`who_gho`** | 2 | `WHOGHOQueryTool` class, `who_gho_tool.py` | ❌ 未割当 | `scientific-epidemiology-public-health` |
| **`medlineplus`** | 5 | `medlineplus_tool.py` | ❌ 未割当 | `scientific-public-health-data` |
| **`loinc`** | 4 | `LOINCTool` class, `loinc_tool.py` | ❌ 未割当 | `scientific-clinical-decision-support` |
| **`odphp`** | 3+ | ODPHP health objectives | ❌ 未割当 | `scientific-public-health-data` |
| **`icd`** | 3+ | ICD-10/11 disease codes | ❌ 未割当 | `scientific-clinical-decision-support` |
| **`umls`** | varies | UMLS medical terminologies | ❌ 未割当 | `scientific-clinical-nlp` |
| **`euhealth`** | varies | EU health surveillance | ❌ 未割当 | `scientific-epidemiology-public-health` |
| `clingen` | 7 | `ClinGenTool` class | ✅ 割当済 (`scientific-clingen-curation`) | — |
| `rxnorm` | varies | Drug name normalization | description にあるが tu_tools 未設定 | `scientific-public-health-data` |

**特筆**: `scientific-public-health-data` は description で NHANES / MedlinePlus / RxNorm / ODPHP を言及しているが **tu_tools が未設定**。Phase 14 で最優先追加。

### 2-B. 免疫学 (Immunology)

| TU Key | Tools | Implementation | SATORI 現状 | 推奨先スキル |
|--------|-------|----------------|-------------|-------------|
| **`imgt`** | 3 | `IMGTTool` class, `imgt_tool.py` | ❌ 未割当 | `scientific-immunoinformatics` |
| **`sabdab`** | 3 | `SAbDabTool` class, `sabdab_tool.py` | ❌ 未割当 | `scientific-immunoinformatics` |
| **`therasabdab`** | 3 | `TheraSAbDabTool` class, `therasabdab_tool.py` (config L328) | ❌ 未割当 | `scientific-immunoinformatics` |
| `iedb` | 8 | `iedb_tools` (key名注意) | ✅ 割当済 (`scientific-immunoinformatics`) | — |

**特筆**: `scientific-immunoinformatics` は description で IEDB/IMGT/SAbDab すべてを言及しているが、tu_tools に `iedb` しか設定されていない。`imgt`, `sabdab`, `therasabdab` を追加すべき。

### 2-C. 環境・生態学 (Environmental / Ecology)

| TU Key | Tools | Implementation | SATORI 現状 | 推奨先スキル |
|--------|-------|----------------|-------------|-------------|
| **`impc`** | 4 | `IMPCTool` class, `impc_tool.py` (config L336) | ❌ 未割当 | `scientific-model-organism-db` |
| **`mpd`** | 1 | `MPDRESTTool` class, `mpd_tools.json` | ❌ 未割当 | `scientific-model-organism-db` |
| `gbif` | 2 | `GBIFTool` | ✅ 割当済 (`scientific-marine-ecology`) | 重複割当可 |
| `obis` | 2 | `OBISTool` | ✅ 割当済 (`scientific-marine-ecology`) | — |
| `worms` | 1 | `WoRMSRESTTool` | ✅ 割当済 (`scientific-marine-ecology`) | — |
| `paleobiology` | 1 | `PaleobiologyRESTTool` | ✅ 割当済 (`scientific-paleobiology`) | — |

**特筆**: `scientific-environmental-ecology` は description で OBIS/GBIF を言及しているが tu_tools 未設定。`gbif` と `obis` の重複割当、または `impc` + `mpd` の新規割当を検討。

### 2-D. システム生物学 (Systems Biology)

| TU Key | Tools | Implementation | SATORI 現状 | 推奨先スキル |
|--------|-------|----------------|-------------|-------------|
| **`bigg_models`** | 7 | BiGG genome-scale metabolic models | ❌ 未割当 | `scientific-systems-biology` or `scientific-metabolic-modeling` |
| **`complex_portal`** | 2 | `ComplexPortalTool` class (config L338) | ❌ 未割当 | `scientific-systems-biology` |
| **`wikipathways`** | 2 | WikiPathways pathway search | ❌ 未割当 | `scientific-systems-biology` or `scientific-pathway-enrichment` |
| `metacyc` | 4 | `MetaCycTool` | ✅ 割当済 (`scientific-metabolomics-databases`) | — |
| `biomodels` | 4 | `BioModelsTool` | ✅ 割当済 (`scientific-metabolic-modeling`) | — |

**特筆**: `scientific-systems-biology` は description で BioModels / Reactome / KEGG / **BiGG** を言及しているが tu_tools 未設定。`bigg_models`, `complex_portal`, `wikipathways` を割当すべき。

---

## 3. 既存 SKILL.md Frontmatter ステータス (8 対象ファイル)

| # | スキル | tu_tools | `---` 行 | 推奨追加キー |
|---|--------|----------|---------|-------------|
| 1 | `scientific-epigenomics-chromatin` | ❌ なし | L9 | `chipatlas`, `fourdn`, `regulomedb`, `jaspar` |
| 2 | `scientific-metabolomics` | ✅ `hmdb` | L12 | `metabolomics_workbench` |
| 3 | `scientific-systems-biology` | ❌ なし | L8 | `bigg_models`, `complex_portal`, `wikipathways` |
| 4 | `scientific-infectious-disease` | ❌ なし | L8 | —（CARD 未実装、候補なし） |
| 5 | `scientific-environmental-ecology` | ❌ なし | L7 | `gbif`(重複割当), `impc`, `mpd` |
| 6 | `scientific-population-genetics` | ❌ なし | L7 | —（`gwas`, `gnomad`, `dbsnp` は別スキルで割当済/候補） |
| 7 | `scientific-toxicology-env` | ❌ なし | L7 | —（CTD/ToxCast 未実装） |
| 8 | `scientific-microbiome-metagenomics` | ✅ `mgnify` | L13 | —（完了） |

### 追加候補マトリクス

```
scientific-epigenomics-chromatin     → chipatlas (4), fourdn (3), regulomedb (1), jaspar (1)
scientific-metabolomics              → metabolomics_workbench (6)  ← PRIORITY
scientific-systems-biology           → bigg_models (7), complex_portal (2), wikipathways (2)
scientific-infectious-disease        → (no TU candidate available)
scientific-environmental-ecology     → gbif (2, 重複可), impc (4), mpd (1)
scientific-population-genetics       → (gwas/gnomad/dbsnp は別スキルに割当推奨)
scientific-toxicology-env            → (no TU candidate available)
scientific-microbiome-metagenomics   → (already complete)
```

---

## 4. 追加発見 — Description-tu_tools 不整合スキル

以下のスキルは description で ToolUniverse データベースを言及しているが **tu_tools が未設定** (Phase 14 での整合化を推奨):

| スキル | Description 記載 DB | 推奨 tu_tools |
|--------|---------------------|--------------|
| **`scientific-public-health-data`** | NHANES, MedlinePlus, RxNorm, ODPHP, Health Disparities | `nhanes`, `medlineplus`, `rxnorm`, `odphp`, `health_disparities` |
| **`scientific-epidemiology-public-health`** | WHO/CDC/EU 公衆衛生 | `who_gho`, `euhealth` |
| **`scientific-immunoinformatics`** | IEDB/IMGT/SAbDab | `imgt`, `sabdab`, `therasabdab` (iedb 割当済) |
| **`scientific-systems-biology`** | BioModels/Reactome/KEGG/BiGG | `bigg_models`, `complex_portal`, `wikipathways` |
| **`scientific-environmental-ecology`** | OBIS/GBIF | `gbif`, `obis` (重複割当) |
| **`scientific-metabolic-modeling`** | BiGG Models | `bigg_models` (biomodels 割当済) |

---

## 5. Phase 14 推奨アクション TOP 10

| Rank | TU Key | ツール数 | 割当先スキル | アクション |
|------|--------|---------|-------------|-----------|
| 1 | `chipatlas` | 4 | `scientific-epigenomics-chromatin` | **新規 tu_tools 追加** |
| 2 | `metabolomics_workbench` | 6 | `scientific-metabolomics` | **既存 tu_tools に追記** |
| 3 | `nhanes` | 2 | `scientific-public-health-data` | **新規 tu_tools 追加** |
| 4 | `medlineplus` | 5 | `scientific-public-health-data` | 同上 (バッチ追加) |
| 5 | `imgt` | 3 | `scientific-immunoinformatics` | **既存 tu_tools に追記** |
| 6 | `sabdab` | 3 | `scientific-immunoinformatics` | 同上 (バッチ追加) |
| 7 | `bigg_models` | 7 | `scientific-systems-biology` | **新規 tu_tools 追加** |
| 8 | `complex_portal` | 2 | `scientific-systems-biology` | 同上 (バッチ追加) |
| 9 | `who_gho` | 2 | `scientific-epidemiology-public-health` | **新規 tu_tools 追加** |
| 10 | `impc` | 4 | `scientific-model-organism-db` or `scientific-environmental-ecology` | **新規 tu_tools 追加** |

### NOT FOUND 最終リスト

| 候補 | ステータス | 代替手段 |
|------|-----------|---------|
| **CARD** (AMR) | ❌ TU 未実装 (v0.20, v0.22 連続不在) | ローカル ABRicate / AMRFinderPlus + `ncbi_nucleotide` |
| **CTD** (Comparative Toxicogenomics) | ❌ TU 未実装 | `disgenet` + `stitch` で代替 |
| **ToxCast / Tox21** | ❌ TU 未実装 | EPA CompTox API (カスタム実装要) |
| **HGNC** | ❌ TU 未実装 | `biothings` (MyGene) で代替 |
| **IUCN Red List** | ❌ TU 未実装 | IUCN API (カスタム実装要) |

---

## 6. default_config.py 全キーインベントリ (v0.22 時点)

`default_config.py` lines 0–347 から確認されたカテゴリキー (全 ~135 キー):

```
special_tools, tool_finder, opentarget, fda_drug_label, monarch, clinical_trials,
fda_drug_adverse_event, fda_drug_adverse_event_detail, ChEMBL, EuropePMC,
semantic_scholar, pubtator, EFO, Enrichr, HumanBase, OpenAlex, literature_search,
arxiv, crossref, simbad, dblp, pubmed, ncbi_nucleotide, ncbi_sra, doaj, unpaywall,
biorxiv, medrxiv, hal, core, pmc, zenodo, openaire, osf_preprints, fatcat,
wikidata_sparql, wikipedia, dbpedia, agents, smolagents, tool_discovery_agents,
web_search_tools, package_discovery_tools, pypi_package_inspector_tools,
drug_discovery_agents, dataset, health_disparities, hpa, reactome, pubchem,
medlineplus, rxnorm, loinc, uniprot, cellosaurus,
software_bioinformatics, software_genomics, software_single_cell,
software_structural_biology, software_cheminformatics, software_machine_learning,
software_visualization, software_scientific_computing, software_physics_astronomy,
software_earth_sciences, software_image_processing, software_neuroscience,
visualization_protein_3d, visualization_molecule_2d, visualization_molecule_3d,
interpro, ebi_search, intact, metabolights, proteins_api, arrayexpress, biostudies,
dbfetch, pdbe_api, ena_browser, blast, cbioportal, regulomedb, jaspar, remap,
screen, pride, emdb, sasbdb, gtopdb, mpd, worms, paleobiology,
go, compose, python_executor, idmap, disease_target_score,
mcp_auto_loader_uspto_downloader, uspto, xml, mcp_auto_loader_boltz, url,
file_download, rcsb_pdb, rcsb_search, tool_composition, embedding, gwas, admetai,
alphafold, output_summarization, odphp, who_gho, umls, icd, euhealth, markitdown,
guidelines, kegg, ensembl, clinvar, geo, dbsnp, gnomad, gbif, obis, wikipathways,
rnacentral, encode, gtex, biomodels_tools, biothings,
fda_pharmacogenomic_biomarkers, metabolomics_workbench, pharmgkb, dgidb, stitch,
civic, cellxgene_census, chipatlas, fourdn, gtex_v2, rfam, bigg_models, ppi,
biogrid, nvidia_nim, cosmic, oncokb, omim, orphanet, disgenet, bindingdb, gpcrdb,
brenda, hmdb, metacyc, zinc, enamine, emolecules, sabdab, imgt, depmap,
interproscan, eve, therasabdab, deepgo, clingen, spliceai, impc, complex_portal,
expression_atlas, proteinsplus, swissdock, pharos, alphamissense, cadd,
nhanes, ols, enrichr, string, faers, admet_ai
```

**v0.21 → v0.22 で新たに確認されたキー**: `nhanes`, `ols`, `enrichr`, `string`, `faers`, `admet_ai` (末尾に追加されている)

---

*Generated by SATORI ToolUniverse Key Verification Workflow*
