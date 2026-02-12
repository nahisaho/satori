# ToolUniverse Key Verification Report — SATORI v0.21.0 Phase 13

**Date**: 2025-07-16  
**Repository**: [mims-harvard/ToolUniverse](https://github.com/mims-harvard/ToolUniverse)  
**Total TU Tools**: 1,265 (per `__init__.py` header)  
**Source Files Analyzed**: `default_config.py`, `_lazy_registry_static.py`, `tools/__init__.py` (`__all__`), individual `*_tool.py` files

---

## 1. PRIORITY A — VERIFIED KEYS

### ✅ `cellxgene_census`
- **Config**: `default_config.py:266` → `cellxgene_census_tools.json`
- **Implementation**: `cellxgene_census_tool.py` (class `CELLxGENECensusTool`)
- **Tools (7)**:
  | Tool Function | Description |
  |---|---|
  | `CELLxGENE_get_census_versions` | List available census data versions |
  | `CELLxGENE_get_cell_metadata` | Query cell metadata by organism/tissue/cell type |
  | `CELLxGENE_get_gene_metadata` | Gene-level metadata from census |
  | `CELLxGENE_get_expression_data` | Gene expression matrix queries |
  | `CELLxGENE_get_presence_matrix` | Gene presence across datasets |
  | `CELLxGENE_get_embeddings` | Cell embeddings (scVI, etc.) |
  | `CELLxGENE_download_h5ad` | Download H5AD files |
- **Bonus**: `get_cellxgene_census_info` (PackageTool for library info)

### ✅ `pharos`
- **Config**: `default_config.py:316` → `pharos_tools.json`
- **Implementation**: `pharos_tool.py` (class `PharosTool`)
- **Tools (4)**:
  | Tool Function | Description |
  |---|---|
  | `Pharos_get_target` | Get target by gene/UniProt with TDL classification |
  | `Pharos_search_targets` | Search targets with TDL filter (Tdark/Tbio/Tchem/Tclin) |
  | `Pharos_get_tdl_summary` | TDL level descriptions & statistics |
  | `Pharos_get_disease_targets` | Disease-associated targets with TDL |

### ✅ `fda_pharmacogenomic_biomarkers`
- **Config**: `default_config.py:249` → `fda_pharmacogenomic_biomarkers_tools.json`
- **Implementation**: `fda_pharmacogenomic_biomarkers_tool.py` (class `FDAPharmacogenomicBiomarkersTool`)
- **Tools (1 dedicated + 2 related)**:
  | Tool Function | Description |
  |---|---|
  | `fda_pharmacogenomic_biomarkers` | Scrape FDA PGx biomarker table, filter by drug/biomarker |
  | `FDA_get_pharmacogenomics_info_by_drug_name` | FDA label PGx info (in `fda_drug_label` category) |
  | `FDA_get_drug_name_by_pharmacogenomics` | Reverse lookup by PGx biomarker |

---

## 2. PRIORITY B — SEARCH RESULTS

| Candidate | Status | Notes |
|---|---|---|
| **PROSITE** | ❌ NOT FOUND | No dedicated `prosite` key. InterPro covers domain/motif queries. |
| **HGNC** | ❌ NOT FOUND | No dedicated `hgnc` key. Gene nomenclature via MyGene (`biothings`) and Ensembl. |
| **MetaboAnalyst** | ❌ NOT FOUND | No dedicated key. Metabolomics via `metabolomics_workbench` and `hmdb`. |
| **Clinical NLP** | ❌ NOT FOUND | No dedicated key. `MedicalTermNormalizer` agent exists in `agents` category. |

---

## 3. NEW DISCOVERIES — Unused TU Keys with 3+ Tools

These keys exist in ToolUniverse `default_config.py` but are **NOT** in SATORI's current `already-used` list. Ordered by scientific relevance to SATORI's mission.

### Tier 1: HIGH PRIORITY (direct scientific value, 4+ tools)

| # | TU Key | Tools | Description |
|---|---|---|---|
| 1 | **`gtex_v2`** | 12 | GTEx Portal V2 — tissue-specific gene expression, eQTLs (V11 Jan 2026) |
| 2 | **`clingen`** | 8 | ClinGen — gene-disease validity, dosage sensitivity, actionability, variant classification |
| 3 | **`hpa`** | 13+ | Human Protein Atlas — expression, subcellular location, cancer prognostics, biological processes |
| 4 | **`gwas`** | 10+ | GWAS Catalog — associations, SNPs, studies, traits, variants |
| 5 | **`biothings`** | 7 | BioThings (MyGene/MyVariant/MyChem) — gene annotation, variant annotation, chemical annotation |
| 6 | **`metabolomics_workbench`** | 6 | Metabolomics Workbench — compound search, study data, RefMet, exact mass search |
| 7 | **`proteinsplus`** | 5 | ProteinsPlus — binding site prediction (DoGSiteScorer), interaction diagrams (PLIP), structure quality |
| 8 | **`gpcrdb`** | 5 | GPCRdb — GPCR protein info, structures, mutations, ligands |
| 9 | **`medlineplus`** | 5 | MedlinePlus — genetics conditions, gene info, topic search, connect lookup |
| 10 | **`chipatlas`** | 4 | ChIP-Atlas — epigenomic landscapes, enrichment analysis, peak data (433K+ experiments) |
| 11 | **`impc`** | 4 | IMPC — mouse knockout phenotypes, gene summary, phenotype hits |
| 12 | **`ncbi_sra`** | 4 | NCBI SRA — sequence read archive search, run info, download URLs, BioSample links |
| 13 | **`loinc`** | 4 | LOINC — lab test codes, forms, answer lists |
| 14 | **`biostudies`** | 4 | EBI BioStudies — study search, files, collections |
| 15 | **`dbfetch`** | 4 | EBI DbFetch — entry fetch, batch fetch, database/format listing |

### Tier 2: MEDIUM PRIORITY (solid value, 3 tools)

| # | TU Key | Tools | Description |
|---|---|---|---|
| 16 | **`spliceai`** | 3 | SpliceAI — deep learning splice prediction, Pangolin, max delta scores |
| 17 | **`swissdock`** | 3 | SwissDock — molecular docking (AutoDock Vina), job management |
| 18 | **`therasabdab`** | 3 | Thera-SAbDab — therapeutic antibody sequences, target search |
| 19 | **`cadd`** | 3 | CADD — combined annotation dependent depletion (variant pathogenicity scoring) |
| 20 | **`interproscan`** | 3 | InterProScan — protein domain/family prediction via sequence scanning |
| 21 | **`fourdn`** | 3 | 4D Nucleome — 3D genome organization, experiment/file metadata |
| 22 | **`imgt`** | 3 | IMGT — immunogenetics (gene info, sequences) |
| 23 | **`dbsnp`** | 3 | dbSNP — variant frequencies, rsID lookup, gene-based search |
| 24 | **`pride`** | 3 | PRIDE — proteomics data, project info, files |
| 25 | **`geo`** | 3 | GEO — gene expression omnibus (dataset info, sample info, search) |
| 26 | **`icd`** | 3+ | ICD-10/ICD-11 — disease code search, hierarchy, entity lookup |

### Tier 3: SPECIALIZED (2 tools or niche domain)

| # | TU Key | Tools | Description |
|---|---|---|---|
| 27 | **`complex_portal`** | 2 | EBI Complex Portal — curated protein complexes (CORUM) |
| 28 | **`nvidia_nim`** | varies | NVIDIA NIM — AlphaFold2, Boltz docking, genomics |
| 29 | **`omim`** | varies | OMIM — Mendelian inheritance in man |
| 30 | **`eve`** | 1+ | EVE — evolutionary variant effect predictions |
| 31 | **`alphamissense`** | 1+ | AlphaMissense — DeepMind pathogenicity predictions |
| 32 | **`deepgo`** | 1 | DeepGO — protein function prediction from sequence |
| 33 | **`disease_target_score`** | many | Disease-target score aggregation (cancer_biomarkers, cancer_gene_census, chembl, europepmc, eva, expression_atlas, genomics_england) |
| 34 | **`sabdab`** | varies | SAbDab — structural antibody database |
| 35 | **`enamine`** | varies | Enamine — make-on-demand compound vendor |
| 36 | **`emolecules`** | varies | eMolecules — vendor aggregator |
| 37 | **`bigg_models`** | varies | BiGG — genome-scale metabolic models |
| 38 | **`who_gho`** | 2 | WHO GHO — global health observatory data |
| 39 | **`euhealth`** | varies | EU Health Info — surveillance, vaccination, obesity, etc. |

### Additional Non-Database Keys (Tools/Utilities)

| TU Key | Description |
|---|---|
| `proteins_api` | EBI Proteins API (UniProt REST alternative) |
| `pdbe_api` | PDBe API (PDB Europe) |
| `blast` | NCBI BLAST sequence alignment |
| `jaspar` | JASPAR transcription factor binding matrices |
| `rxnorm` | RxNorm drug name normalization |
| `umls` | UMLS medical terminologies |
| `pubmed` | PubMed article search & retrieval |
| `gtex` | GTEx v1 (older, prefer `gtex_v2`) |
| `screen` | Screen tools |
| `sasbdb` | SASBDB small-angle scattering |
| `emdb` | EMDB electron microscopy data bank |
| `health_disparities` | SVI, county rankings |
| `wikipathways` | WikiPathways pathway search/get |
| `odphp` | ODPHP health objectives |
| `fatcat` | Fatcat Scholar (open access articles) |
| `gbif` | GBIF biodiversity occurrences |
| `obis` | OBIS ocean biodiversity |
| `worms` | WoRMS marine species taxonomy |
| `paleobiology` | Paleobiology fossils |
| `mpd` | Mouse Phenome Database |

---

## 4. SUMMARY & RECOMMENDATIONS

### Priority A Result: All 3 VERIFIED ✅
All three priority candidates exist as registered TU keys with dedicated tool implementations.

### Top 10 Recommendations for Phase 13

Based on tool count, scientific breadth, and complementarity to SATORI's existing skills:

| Rank | TU Key | Tools | Rationale |
|---|---|---|---|
| 1 | `gtex_v2` | 12 | Tissue expression is fundamental; V2 API with eQTL support |
| 2 | `clingen` | 8 | Gene-disease validity + actionability fills clinical genetics gap |
| 3 | `hpa` | 13+ | Protein Atlas is critical for expression/localization/cancer |
| 4 | `biothings` | 7 | MyGene/MyVariant/MyChem are universal bio-annotation tools |
| 5 | `gwas` | 10+ | Genome-wide association data essential for genetics skills |
| 6 | `metabolomics_workbench` | 6 | Fills metabolomics gap alongside HMDB |
| 7 | `chipatlas` | 4 | Epigenomics data (ChIP-seq/ATAC-seq) is unique capability |
| 8 | `proteinsplus` | 5 | Binding site analysis complements structure tools |
| 9 | `spliceai` | 3 | Splice variant pathogenicity is high clinical value |
| 10 | `cadd` | 3 | Variant pathogenicity scoring is foundational for clinical genomics |

### NOT FOUND in ToolUniverse
The following requested tools have **no** dedicated TU key:
- **PROSITE** — Use `interpro` for motif/domain queries
- **HGNC** — Use `biothings` (MyGene) for gene nomenclature
- **MetaboAnalyst** — Use `metabolomics_workbench` + `hmdb`
- **Clinical NLP** — No dedicated tools; agents category has `MedicalTermNormalizer`
- **Pfam/CATH/SCOP** — Covered by `interpro` and `interproscan`
- **CTD** — Not found (use `disgenet` + `stitch` for chemical-disease/gene)
- **ToxCast** — Not found
- **LINCS** — Not found
- **BioAssay** — PubChem BioAssay tools exist within `pubchem` (already used)

---

## 5. FULL default_config.py KEY INVENTORY

Complete list of all TU category keys from `default_config.py` (lines 0-345):

```
special_tools, tool_finder, opentarget, fda_drug_label, monarch, clinical_trials,
fda_drug_adverse_event, fda_drug_adverse_event_detail, ChEMBL, EuropePMC,
semantic_scholar, pubtator, EFO, Enrichr, HumanBase, OpenAlex, literature_search,
arxiv, crossref, simbad, dblp, pubmed, ncbi_nucleotide, ncbi_sra, doaj, unpaywall,
biorxiv, medrxiv, hal, core, pmc, zenodo, openaire, osf_preprints, fatcat,
wikidata_sparql, wikipedia, dbpedia, agents, smolagents, tool_discovery_agents,
web_search_tools, package_discovery_tools, pypi_package_inspector_tools,
drug_discovery_agents, dataset, health_disparities, hpa, reactome, pubchem,
medlineplus, rxnorm, loinc, uniprot, cellosaurus, software_bioinformatics,
software_scientific_computing, software_physics_astronomy, software_earth_sciences,
software_image_processing, software_neuroscience, visualization_molecule_2d,
visualization_molecule_3d, interpro, ebi_search, intact, metabolights, proteins_api,
arrayexpress, biostudies, dbfetch, pdbe_api, ena_browser, blast, cbioportal,
regulomedb, jaspar, remap, screen, pride, emdb, sasbdb, gtopdb, mpd, worms,
paleobiology, go, compose, python_executor, idmap, disease_target_score,
mcp_auto_loader_uspto_downloader, uspto, xml, mcp_auto_loader_boltz, url,
file_download, rcsb_pdb, rcsb_search, tool_composition, embedding, gwas, admetai,
alphafold, output_summarization, odphp, who_gho, umls, icd, euhealth, markitdown,
guidelines, kegg, ensembl, clinvar, geo, dbsnp, gnomad, gbif, obis, wikipathways,
rnacentral, encode, gtex, biomodels_tools, biothings, fda_pharmacogenomic_biomarkers,
metabolomics_workbench, pharmgkb, dgidb, stitch, civic, cellxgene_census, chipatlas,
fourdn, gtex_v2, rfam, bigg_models, ppi, biogrid, nvidia_nim, cosmic, oncokb, omim,
orphanet, disgenet, bindingdb, gpcrdb, brenda, hmdb, metacyc, zinc, enamine,
emolecules, sabdab, imgt, depmap, interproscan, eve, therasabdab, deepgo, clingen,
spliceai, impc, complex_portal, expression_atlas, proteinsplus, swissdock, pharos,
alphamissense, cadd, nhanes, ols, enrichr, string, faers, admet_ai
```

**Total registered category keys**: ~130+
