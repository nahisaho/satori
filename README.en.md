# SATORI — Agent Skills for Science

**SATORI** is a collection of **GitHub Copilot Agent Skills** for scientific data analysis.

[![npm version](https://img.shields.io/npm/v/@nahisaho/satori)](https://www.npmjs.com/package/@nahisaho/satori)
[![CI](https://github.com/nahisaho/satori/actions/workflows/ci.yml/badge.svg)](https://github.com/nahisaho/satori/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENCE)

## Overview

This directory contains **190 skills** that systematize scientific data analysis techniques accumulated through Exp-01–13 as Agent Skills. Copilot automatically loads the appropriate skill based on the prompt context, reusing established analysis patterns from each experiment. All 190 skills can integrate with over 1,200 external scientific database tools via [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP.

### Pipeline Flow

```
hypothesis-pipeline → pipeline-scaffold → academic-writing → critical-review
  (Hypothesis Def.)    (Analysis Exec.)    (Draft Writing)    (Review & Revision)
                                                                       ↓
  paper-quality ← revision-tracker ← peer-review-response ← [Peer Review Received]
  (Quality Eval.)  (Revision Tracking) (Review Response)      (Journal)
```

**Domain-Specific Pipelines (J–O)**

```
research-methodology → grant-writing → hypothesis-pipeline    ← [O→A Research Plan]
        │                    ↑                   ↓
regulatory-science ───┘    scientific-schematics
  (FDA/ISO/Patents)       (Research Diagrams)
                                              ↓
drug-target-profiling → admet-pharmacokinetics ─→ drug-repurposing
  (Target ID)            (ADMET/PK Evaluation)    (Repositioning)
        ↓                        │                    ↓
protein-structure-analysis → protein-design → lab-automation
  (Structure Analysis)      (de novo Design)  (Experiment Automation)
        │                                             │
protein-interaction-network    molecular-docking    lab-data-management
  (PPI: STRING/IntAct)        (Vina/DiffDock)    (Benchling/DNAnexus/OMERO)
protein-domain-family                                    ↓
  (InterPro/Pfam)                              ↓
variant-interpretation → variant-effect-prediction → clinical-decision-support → presentation-design
  (Variant Interpret.)  (AlphaMissense/CADD/SpliceAI) (Clinical Decision)       (Conference Present.)
        ↑                                             ↓
pharmacogenomics                             clinical-reporting
  (PGx Metabolizer)                          (SOAP/FHIR Reports)
```

Files generated at each step are automatically passed to the next step:

**Advanced Computing & Medical Pipelines (P–S)**

```
pharmacovigilance ← admet-pharmacokinetics       ← [P Safety Surveillance]
  (Post-market Safety)  (Preclinical ADMET)              ↓
        │                                   regulatory-science
        │                                   (FDA/ISO/Patents)
        ↓                                          ↓
precision-oncology → clinical-decision-support → medical-imaging
  (Tumor Profiling)   (Clinical Decision)        (Imaging Diagnosis)
        ↓              ↑                          ↓
cancer-genomics → disease-research → variant-interpretation      → deep-learning
  (COSMIC/DepMap)  (Disease-Gene)    (Variant Interpret.)         (DL Frameworks)
                        ↑                          ↓
              pharmacogenomics              neuroscience-electrophysiology
              (PGx Metabolizer)            (Spike Sorting/EEG/HRV)
                                                     ↓
quantum-computing → bayesian-statistics → graph-neural-networks
  (Quantum Comp.)    (Bayesian Inference)  (GNN Molecular Pred.)
        │                                     ↓
computational-materials  explainable-ai ← deep-learning ← [R Advanced Computing]
  (pymatgen/VASP)        (XAI Explainability) (DL Pipeline)
```

**Next-Gen Omics & Epidemiology Pipelines (T–Z)**

```
single-cell-genomics → spatial-transcriptomics     ← [T Single-Cell & Spatial]
  (scRNA-seq QC)       (Visium/MERFISH)
        ↓                     ↓
epigenomics-chromatin → gene-expression-transcriptomics
  (ChIP-seq/ATAC/WGBS)   (GEO/GTEx/DESeq2)
        │                     │
        └─────────────┬─────┘
                          ↓
proteomics-mass-spectrometry → multi-omics      ← [F Omics Integration]
  (LC-MS/MS/PTM/GNPS)         (Integrated Analysis)
        │                          ↓
        │              pathway-enrichment ← metabolomics-databases
        │              (KEGG/Reactome/GO)  (HMDB/MetaCyc/MWB)
                                    ↓
immunoinformatics → infectious-disease             ← [U Immune & Infectious]
  (Epitope Prediction)  (AMR & Phylogenetics)
        ↓                     ↓
microbiome-metagenomics → environmental-ecology    ← [V Microbiome & Environment]
  (16S/Metagenomics)       (SDM & Biodiversity)
        ↓                     ↓
systems-biology           population-genetics       ← [W+Y Modeling & Population]
  (SBML/FBA/GRN)           (Fst/ADMIXTURE)
        ↓                     ↓
epidemiology-public-health → text-mining-nlp        ← [X+Z Epidemiology & NLP]
  (RR/OR/Spatial Clusters)    (NER/KG/BERTopic)
        ↑                          ↓
clinical-trials-analytics    literature-search → systematic-review
  (ClinicalTrials.gov)       (PubMed/OpenAlex)   (PRISMA 2020)
```

**Data Access & Pipeline Integration (v0.14.0 New)**

```
preprint-archive ───→ literature-search ───→ systematic-review
  (bioRxiv/arXiv/CORE)  (PubMed/OpenAlex)    (PRISMA 2020)
        │                       ↓
        │              biomedical-pubtator ───→ text-mining-nlp
        │              (PubTator NER)           (KG Construction)
        │                       ↓
        └──────────────→ deep-research
                        (Evidence Synthesis)

ontology-enrichment ──→ disease-research ───→ variant-interpretation
  (EFO/OLS/UMLS/Enrichr) (GWAS/DisGeNET)      (ACMG/AMP)
        │                       ↓
        │              regulatory-genomics ──→ epigenomics-chromatin
        │              (RegulomeDB/ReMap/4DN)  (ChIP-seq/ATAC)
        ↓
public-health-data ──→ epidemiology-public-health → clinical-decision-support
  (NHANES/RxNorm/CDC)   (RR/OR/Spatial Clusters)    (GRADE Evidence)

ebi-databases ────→ genome-sequence-tools → bioinformatics
  (ENA/EBI Search)    (BLAST/NCBI)           (scRNA/PPI)
        │
        └──→ metabolomics-databases ──→ metabolic-modeling
             (MetaboLights/HMDB)        (BiGG/BioModels)

cell-line-resources ──→ cancer-genomics ──→ precision-oncology
  (Cellosaurus STR)     (COSMIC/DepMap)     (MTB Report)

phylogenetics ────→ microbiome-metagenomics → environmental-ecology
  (ETE3/scikit-bio)   (UniFrac/16S)            (SDM/Biodiversity)

reinforcement-learning → doe → lab-automation
  (SB3/PufferLib)       (DOE) (Robot Control)

symbolic-mathematics ──→ systems-biology ──→ admet-pharmacokinetics
  (SymPy Analytical)     (SBML/FBA)         (PK Modeling)
```

| Phase | Generated Files | Referenced By |
|---|---|---|
| Hypothesis Formulation | `docs/hypothesis.{md,json}`, `docs/workflow_design.{md,json}` | → scaffold, writing |
| Analysis Execution | `results/analysis_summary.json`, `figures/*.png` | → writing |
| Draft Writing | `manuscript/manuscript.md` | → critical-review, citation-checker |
| Review | `manuscript/review_report.{md,json}`, `manuscript/manuscript_revised.md` | → latex-formatter |
| Citation Verification | `manuscript/citation_report.json` | → latex-formatter |
| SI Generation | `manuscript/supplementary.md`, `manuscript/si_crossref_report.json` | → latex-formatter |
| Peer Review Response | `manuscript/response_to_reviewers.md`, `manuscript/response_mapping.json` | → revision-tracker |
| Revision Tracking | `manuscript/manuscript_tracked.md`, `manuscript/revision_summary.json` | → paper-quality |
| Quality Evaluation | `manuscript/quality_report.json` | → latex-formatter |
| LaTeX Conversion | `manuscript/manuscript.tex`, `manuscript/references.bib` | — |
| Target Profiling | `results/target_profile_report.md`, `results/target_profile.json` | → admet-pk, protein-structure, drug-repurposing |
| ADMET/PK Evaluation | `results/admet_profile.json`, `results/pk_model.json` | → drug-repurposing, clinical-decision |
| Drug Repositioning | `results/repurposing_candidates.json`, `results/network_proximity.json` | → clinical-decision |
| Structure Analysis | `results/structure_analysis.json`, `results/binding_sites.json` | → protein-design, cheminformatics |
| Protein Design | `results/design_candidates.json`, `results/esm_scores.json` | → lab-automation, admet-pk |
| Variant Interpretation | `results/variant_classification.json`, `results/pgx_report.json` | → clinical-decision |
| Clinical Decision | `results/clinical_recommendation.json`, `results/trial_matches.json` | → presentation, writing |
| Lab Automation | `protocols/protocol.py`, `results/qc_report.json` | → data-preprocessing |
| Presentation | `presentation/slides.md`, `presentation/poster.tex` | — |
| Research Methodology | `docs/methodology_design.md`, `docs/study_design.json` | → grant-writing, doe |
| Grant Writing | `grants/specific_aims.md`, `grants/research_strategy.md` | → hypothesis-pipeline |
| Pharmacovigilance | `results/pv_signal_report.{md,json}`, `figures/pv_temporal_trend.png` | → clinical-decision |
| Precision Oncology | `results/mtb_report.{md,json}`, `results/variant_actionability.json` | → clinical-decision, writing |
| Disease Research | `results/disease_research_report.{md,json}`, `results/gwas_significant_loci.json` | → variant-interpretation |
| Quantum Computing | `results/quantum_result.json`, `figures/quantum_convergence.png` | → bayesian, cheminformatics |
| GNN | `results/gnn_predictions.json`, `figures/gnn_training_curve.png` | → drug-target, admet |
| Bayesian Statistics | `results/bayesian_summary.json`, `figures/bayesian_trace.png` | → doe, meta-analysis |
| Explainable AI | `results/xai_report.json`, `figures/shap_summary.png` | → clinical-decision |
| Deep Learning | `results/dl_training_log.json`, `models/model.onnx` | → GNN, medical-imaging |
| Medical Imaging | `results/imaging_report.{md,json}`, `results/radiomics_features.json` | → precision-oncology |
| scRNA-seq Analysis | `results/sc_markers.json`, `figures/umap_clusters.png`, `results/rna_velocity.json` | → spatial, systems-biology |
| Spatial Transcriptomics | `results/spatial_domains.json`, `figures/spatial_svg_map.png` | → single-cell, systems-biology |
| Immunoinformatics | `results/epitope_candidates.json`, `results/tcr_diversity.json` | → infectious-disease, drug-target |
| Infectious Disease Genomics | `results/amr_report.json`, `results/mlst_profile.json`, `results/sir_simulation.json` | → epidemiology, microbiome |
| Microbiome | `results/asv_table.json`, `results/diversity_metrics.json`, `results/da_results.json` | → environmental-ecology |
| Environmental Ecology | `results/sdm_predictions.json`, `results/biodiversity_indices.json` | → microbiome, text-mining |
| Systems Biology | `results/sbml_timecourse.json`, `results/fba_fluxes.json`, `results/grn_edges.json` | → multi-omics, network-analysis |
| Epidemiology & Public Health | `results/epi_risk_measures.json`, `results/spatial_clusters.json`, `results/dag_analysis.json` | → survival-clinical, causal-inference |
| Population Genetics | `results/pop_structure.json`, `results/fst_matrix.json`, `results/selection_scan.json` | → disease-research, variant-interpretation |
| Text Mining | `results/ner_entities.json`, `results/knowledge_graph.json`, `results/topic_model.json` | → deep-research, meta-analysis |
| Neuroelectrophysiology | `results/spike_sorting.json`, `results/eeg_erp.json`, `results/connectivity.json` | → biosignal, deep-learning |
| Proteomics | `results/protein_quant.csv`, `results/ptm_sites.json`, `results/molecular_network.json` | → multi-omics, network-analysis |
| Transcriptomics | `results/deseq2_results.csv`, `results/gsea/`, `figures/volcano_rnaseq.png` | → bioinformatics, multi-omics |
| Computational Materials | `results/structure.cif`, `figures/phase_diagram.png`, `figures/band_structure.png` | → quantum-computing |
| Clinical Trials Analytics | `results/clinical_trials.csv`, `results/competitive_landscape.json` | → survival-clinical, meta-analysis |
| Lab Data Management | `results/benchling_sequences.json`, `results/dnanexus_workflow_output.json` | → bioinformatics, lab-automation |
| Scientific Schematics | `figures/consort_flow.svg`, `figures/nn_architecture.svg`, `figures/pathway.md` | → presentation, writing |
| Regulatory Science | `results/fda_orange_book.json`, `results/510k_clearances.csv`, `results/patent_search.csv` | → pharmacovigilance, clinical-trials |
| Pharmacogenomics | `results/pgx_report.json`, `results/cpic_recommendations.json` | → variant-interpretation, clinical-decision |
| Epigenomics | `results/peak_calls.bed`, `results/dmr_results.csv`, `results/chromatin_states.bed` | → single-cell, multi-omics |
| Pathway Enrichment | `results/ora_results.csv`, `results/gsea_results.csv`, `figures/enrichment_heatmap.png` | → gene-expression, metabolomics, multi-omics |
| Literature Search | `results/literature_search.csv`, `results/citation_network.json` | → deep-research, systematic-review |
| PPI Network | `results/ppi_network.json`, `figures/ppi_network.png` | → drug-target, network-analysis |
| Variant Effect Prediction | `results/variant_predictions.csv`, `results/consensus_pathogenicity.json` | → variant-interpretation, cancer-genomics |
| Cancer Genomics | `results/cosmic_mutations.csv`, `results/mutation_signatures.json` | → precision-oncology, variant-interpretation |
| Metabolite DB | `results/hmdb_search.csv`, `results/mz_identification.csv` | → metabolomics, pathway-enrichment |
| Molecular Docking | `results/docking_results.csv`, `results/vina_poses.pdbqt` | → drug-target, protein-structure |
| Systematic Review | `results/screening_records.csv`, `figures/prisma_flow.mmd` | → literature-search, meta-analysis |
| Clinical Reporting | `reports/clinical_report.html`, `reports/clinical_report.fhir.json` | → variant-interpretation, clinical-decision |
| Domain/Family | `results/interproscan_results.csv`, `figures/domain_architecture.png` | → protein-structure, protein-interaction |

### ToolUniverse MCP Tool Integration

131 skills (HIGH 13 + MEDIUM 9 + Phase 3: 20 + Phase 4: 8 + Phase 5: 9 + Phase 6: 7 + Phase 7: 4 + Phase 8: 4 + Phase 9: 5 + Phase 10: 6 + Phase 11: 8 new + 6 existing + Phase 12: 3 new + 12 existing key additions + Phase 13: 3 new + 7 existing key additions + Phase 14: 1 new + 6 existing key additions) can access over 1,200 external scientific tools via the [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP server. Corresponding tools are listed in the `### Available Tools` section within each SKILL.md.

```
SATORI Skill (Methodology/Judgment)    ToolUniverse SMCP (Data Retrieval/Computation)
┌──────────────────────┐        ┌─────────────────────────────┐
│ pharmacovigilance    │───MCP──│ FAERS, FDA Labels, DailyMed │
│ precision-oncology   │───MCP──│ OncoKB, CIViC, COSMIC, GDC  │
│ disease-research     │───MCP──│ OpenTargets, HPO, Monarch    │
│ drug-target-profiling│───MCP──│ UniProt, ChEMBL, DGIdb       │
│ variant-interpretation│──MCP──│ ClinVar, gnomAD, ClinGen     │
│ admet-pharmacokinetics│──MCP──│ ADMET-AI, PubChem, ChEMBL    │
│ pathway-enrichment   │───MCP──│ KEGG, Reactome, GO, WP       │
│ literature-search    │───MCP──│ PubMed, EuropePMC, OpenAlex  │
│ protein-interaction  │───MCP──│ STRING, IntAct, STITCH       │
│ variant-effect-pred  │───MCP──│ AlphaMissense, CADD, SpliceAI│
│ cancer-genomics      │───MCP──│ COSMIC, cBioPortal, DepMap   │
│ metabolomics-databases│──MCP──│ HMDB, MetaCyc, MWB           │
│ protein-domain-family│───MCP──│ InterPro, InterProScan       │
│ systematic-review    │───MCP──│ PubMed (shared)              │
│ rare-disease-genetics│───MCP──│ OMIM, Orphanet, DisGeNET     │
│ human-protein-atlas  │───MCP──│ HPA tissue/RNA/cancer        │
│ pharmacology-targets │───MCP──│ BindingDB, GPCRdb, GtoPdb    │
│ genome-sequence-tools│───MCP──│ dbSNP, BLAST, NCBI, GDC      │
│ biothings-idmapping  │───MCP──│ MyGene, MyVariant, MyChem    │
│ noncoding-rna        │───MCP──│ Rfam, RNAcentral             │
│ structural-proteomics│───MCP──│ EMDB, PDBe, Proteins API     │
│ compound-screening   │───MCP──│ ZINC                         │
│ metabolic-modeling   │───MCP──│ BiGG Models, BioModels       │
│ preprint-archive     │───MCP──│ bioRxiv, arXiv, CORE, Zenodo │
│ public-health-data   │───MCP──│ NHANES, MedlinePlus, RxNorm  │
│ ontology-enrichment  │───MCP──│ EFO, OLS, Enrichr, UMLS      │
│ ebi-databases        │───MCP──│ EBI Search, ENA, MetaboLights│
│ cell-line-resources  │───MCP──│ Cellosaurus                  │
│ regulatory-genomics  │───MCP──│ RegulomeDB, ReMap, 4DN       │
│ biomedical-pubtator  │───MCP──│ PubTator NER                 │
│ chembl-assay-mining  │───MCP──│ ChEMBL Assay/Activity/Target  │
│ ensembl-genomics     │───MCP──│ Ensembl REST, VEP             │
│ string-network-api   │───MCP──│ STRING, BioGRID, STITCH       │
│ expression-comparison│───MCP──│ Expression Atlas              │
│ rrna-taxonomy        │───MCP──│ MGnify metagenomics          │
│ marine-ecology       │───MCP──│ OBIS, WoRMS, GBIF            │
│ encode-screen        │───MCP──│ ENCODE, SCREEN, ChIP-Atlas    │
│ human-cell-atlas     │───MCP──│ HCA Data Portal, CELLxGENE    │
│ geo-expression       │───MCP──│ GEO (NCBI E-utilities)        │
│ paleobiology         │───MCP──│ Paleobiology Database (PBDB)  │
│ gwas-catalog           │───MCP──│ GWAS Catalog (EBI)           │
│ alphafold-structures   │───MCP──│ AlphaFold DB API             │
│ arrayexpress-expression│───MCP──│ ArrayExpress, BioStudies     │
│ semantic-scholar       │───MCP──│ Semantic Scholar Graph API   │
│ pharmgkb-pgx          │───MCP──│ PharmGKB, CPIC Guidelines    │
│ crossref-metadata     │───MCP──│ CrossRef DOI/Metadata        │
│ uniprot-proteome     │───MCP──│ UniProt REST API              │
│ rcsb-pdb-search      │───MCP──│ RCSB PDB Search/Data API      │
│ opentargets-genetics  │───MCP──│ Open Targets GraphQL           │
│ reactome-pathways    │───MCP──│ Reactome Content Service       │
│ depmap-dependencies  │───MCP──│ DepMap Portal, Cell Model      │
│ drugbank-resources   │───MCP──│ DrugBank API                   │
│ civic-evidence       │───MCP──│ CIViC REST API                 │
│ gnomad-variants      │───MCP──│ gnomAD GraphQL API             │
│ monarch-ontology   │───MCP──│ Monarch Initiative API        │
│ gdc-portal         │───MCP──│ NCI GDC REST API              │
│ stitch-chemical-net│───MCP──│ STITCH Chemical-Protein       │
│ drug-repurposing   │───MCP──│ Pharos IDG Targets            │
│ pharmacogenomics   │───MCP──│ FDA PGx Biomarkers            │
│ gtex-tissue-expr   │───MCP──│ GTEx v2 REST API              │
│ protein-structure  │───MCP──│ ProteinsPlus Binding Sites    │
│ cellxgene-census   │───MCP──│ CELLxGENE Census API          │
│ pharos-targets     │───MCP──│ Pharos GraphQL API            │
│ clingen-curation   │───MCP──│ ClinGen Validity/Dosage       │
│ ... (131 skills total)│      │ ... (1,200+ tools)           │
└──────────────────────┘        └─────────────────────────────┘
```

Skills are classified into **26 categories**.

| Category | # Skills | Overview |
|---|:---:|---|
| A. Foundation & Workflow | 17 | Pipeline construction, preprocessing, data generation, figures, writing, hypothesis formulation, critical review, SI generation, LaTeX conversion, citation verification, peer review response, revision tracking, paper quality, systematic review, BioThings ID mapping, data submission, CrossRef metadata |
| B. Statistics & Exploratory Analysis | 11 | EDA, hypothesis testing, dimensionality reduction, symbolic math, missing data analysis, advanced visualization, data profiling, geospatial analysis, network visualization, statistical simulation, streaming analytics |
| C. Machine Learning & Modeling | 11 | Regression, classification, feature importance, active learning, AutoML, ensemble learning, anomaly detection, causal ML, model monitoring, semi-supervised learning, multi-task learning |
| D. Experimental Design & Process Optimization | 3 | DOE, response surface methodology, Bayesian optimization, adaptive experimental design |
| E. Signal, Spectral & Time Series | 5 | Spectral analysis, biosignal, time series decomposition, neuroelectrophysiology, ML time series forecasting |
| F. Life Sciences & Omics | 28 | Bioinformatics, metabolomics, genome sequence, multi-omics, network, proteomics, transcriptomics, pathway enrichment, metabolite DB, HPA, genome sequence tools, noncoding RNA, ontology, EBI DB group, Ensembl genomics, STRING/BioGRID PPI, expression comparison, model organism DB, GEO expression profiles, parasite genomics, ArrayExpress expression archive, GTEx tissue expression, UniProt proteome, Reactome pathways, HGNC nomenclature, metabolic network, glycomics, lipidomics |
| G. Chemistry, Materials & Imaging | 9 | Cheminformatics, materials characterization, image morphology, computational materials, ChEMBL assay mining, MD simulation, advanced imaging, deep chemistry, STITCH chemical-protein network |
| H. Clinical, Epidemiology & Meta-Science | 7 | Survival analysis, causal inference, meta-analysis, clinical trials analytics, clinical reporting, biobank large-scale cohort, clinical standard terminology |
| I. Deep Research & Literature Search | 4 | Scientific literature deep research, evidence hierarchy evaluation, multi-DB literature search, citation network, preprint cross-search, Semantic Scholar academic graph |
| J. Drug Discovery & Pharmacology | 9 | Target profiling, ADMET/PK, drug repositioning, molecular docking, pharmacological targets, compound screening, NCI-60 screening, DrugBank resources, Pharos targets |
| K. Structural Biology & Protein Engineering | 7 | PDB/AlphaFold structure analysis, de novo protein design, PPI network, domain/family, structural proteomics, AlphaFold DB structure prediction, RCSB PDB structure search |
| L. Precision Medicine & Clinical Decision | 6 | Variant interpretation (ACMG/AMP), evidence-based clinical decision, variant effect prediction, CIViC clinical evidence, gnomAD variants, ClinGen curation |
| M. Lab Automation & Data Management | 3 | Liquid handling, protocol management, ELN/LIMS integration, lab data management, CRISPR gRNA design |
| N. Scientific Presentation & Schematics | 4 | Scientific slides, posters, workflow diagrams, scientific schematics, interactive dashboards, reproducible reports |
| O. Research Planning, Grants & Regulation | 3 | Grant applications, research methodology, ethics review, regulatory science |
| P. Pharmacovigilance & Pharmacogenomics | 4 | FAERS disproportionality analysis, MedDRA hierarchy, safety signal detection, PGx metabolizer, PharmGKB clinical annotations, clinical pharmacology PopPK/PBPK |
| Q. Oncology & Disease Research | 10 | Precision oncology (CIViC/OncoKB), disease-gene association (GWAS/Orphanet), cancer genomics (COSMIC/DepMap), rare disease genetics, cell line resources, ICGC cancer genome data, Open Targets genetics, DepMap dependencies, Monarch ontology, GDC portal |
| R. Quantum & Advanced Computing | 11 | Quantum computing, GNN, Bayesian statistics, explainable AI, deep learning, healthcare AI, reinforcement learning, transfer learning, uncertainty quantification, federated learning, NAS |
| S. Medical Imaging | 2 | DICOM/NIfTI, WSI pathology images, radiomics, MONAI, radiology AI |
| T. Single-Cell, Spatial & Epigenomics | 13 | scRNA-seq, Visium, MERFISH, CELLxGENE, RNA velocity, epigenomics, regulatory genomics, perturbation analysis, scVI integration, scATAC-seq/Signac, GPU single-cell, ENCODE/SCREEN, Human Cell Atlas, advanced Squidpy spatial analysis, spatial multi-omics, CELLxGENE Census |
| U. Immune & Infectious Disease | 2 | Immunoinformatics, MHC binding prediction, pathogen genomics, AMR, IEDB |
| V. Microbiome & Environment | 9 | 16S/metagenomics, α/β diversity, SDM, OBIS, GBIF, phylogenetics, rRNA taxonomy, plant biology, marine ecology, environmental geodata, paleobiology, MAG reconstruction |
| W. Systems Biology | 4 | SBML simulation, FBA, GRN inference, BioModels, metabolic modeling, Metabolic Atlas, metabolic flux analysis |
| X. Epidemiology & Public Health | 3 | Risk measures (RR/OR), age standardization, spatial epidemiology, WHO, CDC, public health data, environmental toxicology |
| Y. Population Genetics | 2 | HWE, PCA/ADMIXTURE, Fst, selection scan, gnomAD, GWAS, GWAS Catalog |
| Z. Scientific Text Mining | 3 | NER, relation extraction, knowledge graph, BERTopic, PubTator bio-annotation, clinical NLP |

---

## Skills List

### A. Foundation & Workflow (17 Skills)

Cross-cutting foundational skills common to all experiments.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 1 | [scientific-pipeline-scaffold](scientific-pipeline-scaffold/SKILL.md) | Pipeline scaffolding, directory structure, StepLogger, JSON summary | All Exp |
| 2 | [scientific-data-preprocessing](scientific-data-preprocessing/SKILL.md) | Missing value imputation, encoding, scaling, outlier handling | All Exp |
| 3 | [scientific-data-simulation](scientific-data-simulation/SKILL.md) | Physics/chemistry/biology-based synthetic data generation | 06-09, 12, 13 |
| 4 | [scientific-publication-figures](scientific-publication-figures/SKILL.md) | Publication-quality figures, rcParams, color palettes, multi-panel | 10, 11-13 |
| 5 | [scientific-academic-writing](scientific-academic-writing/SKILL.md) | Scientific paper writing, journal-specific templates, cover letters, peer review response | General |
| 6 | [scientific-hypothesis-pipeline](scientific-hypothesis-pipeline/SKILL.md) | Hypothesis formulation from prompts, PICO/PECO structuring, auto-generation of analysis pipeline | General |
| 7 | [scientific-critical-review](scientific-critical-review/SKILL.md) | Critical review of drafts, discussion deepening, logic verification, revision suggestions | General |
| 8 | [scientific-supplementary-generator](scientific-supplementary-generator/SKILL.md) | Auto-generation of Supplementary Information, SI figure/table organization, main text-SI cross-reference verification | General |
| 9 | [scientific-latex-formatter](scientific-latex-formatter/SKILL.md) | Markdown→LaTeX conversion, journal template application, BibTeX generation | General |
| 10 | [scientific-citation-checker](scientific-citation-checker/SKILL.md) | Automated citation search, coverage check, consistency verification, duplicate detection | General |
| 11 | [scientific-peer-review-response](scientific-peer-review-response/SKILL.md) | Structured peer review comments, point-by-point responses, rebuttal letter generation | General |
| 12 | [scientific-revision-tracker](scientific-revision-tracker/SKILL.md) | Revision history tracking, diff management, change markup, traceability verification | General |
| 13 | [scientific-paper-quality](scientific-paper-quality/SKILL.md) | Readability score, structural balance, vocabulary quality, journal fit, reproducibility check | General |
| 84 | [scientific-systematic-review](scientific-systematic-review/SKILL.md) | PRISMA 2020 systematic review, multi-DB search strategy, screening, risk of bias assessment | General |
| 91 | [scientific-biothings-idmapping](scientific-biothings-idmapping/SKILL.md) | BioThings API (MyGene/MyVariant/MyChem) cross-database ID mapping & annotation | General |
| 119 | [scientific-data-submission](scientific-data-submission/SKILL.md) | GenBank/SRA/GEO/BioProject/BioSample data submission, FAIR principles compliance | General |
| 139 | [scientific-crossref-metadata](scientific-crossref-metadata/SKILL.md) | CrossRef REST API DOI resolution, paper metadata, citation counts, journal information | General |

### B. Statistics & Exploratory Analysis (11 Skills)

Skills for data understanding, testing, dimensionality reduction, symbolic math, missing data analysis, advanced visualization, data profiling, geospatial analysis, network visualization, statistical simulation, and streaming analytics.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 14 | [scientific-eda-correlation](scientific-eda-correlation/SKILL.md) | Exploratory data analysis, correlation heatmaps, distribution visualization | 02, 12, 13 |
| 15 | [scientific-statistical-testing](scientific-statistical-testing/SKILL.md) | Hypothesis testing, multiple comparisons, enrichment, Bayesian inference | 03, 04, 06, 07 |
| 16 | [scientific-pca-tsne](scientific-pca-tsne/SKILL.md) | PCA / t-SNE / UMAP dimensionality reduction, clustering | 02, 03, 07, 11, 13 |
| 105 | [scientific-symbolic-mathematics](scientific-symbolic-mathematics/SKILL.md) | SymPy analytical calculus, ODE solving, linear algebra, symbolic scientific modeling | General |
| 172 | [scientific-missing-data-analysis](scientific-missing-data-analysis/SKILL.md) | Missing pattern diagnosis (MCAR/MAR/MNAR), Little's MCAR test, MICE multiple imputation, KNN/MissForest imputation | General |
| 173 | [scientific-advanced-visualization](scientific-advanced-visualization/SKILL.md) | Plotly 3D, Altair declarative visualization, parallel coordinates, publication-quality figures, animation | General |
| 179 | [scientific-data-profiling](scientific-data-profiling/SKILL.md) | ydata-profiling automated EDA, data quality score (5 dimensions), Great Expectations validation | General |
| 180 | [scientific-geospatial-analysis](scientific-geospatial-analysis/SKILL.md) | GeoPandas geospatial processing, Moran's I/LISA spatial autocorrelation, Kriging interpolation, Folium maps | General |
| 181 | [scientific-network-visualization](scientific-network-visualization/SKILL.md) | NetworkX graph construction, Louvain/Leiden community detection, centrality analysis, PyVis interactive | General |
| 187 | [scientific-statistical-simulation](scientific-statistical-simulation/SKILL.md) | Monte Carlo simulation, Bootstrap inference (BCa), Permutation Test, power analysis | General |
| 188 | [scientific-streaming-analytics](scientific-streaming-analytics/SKILL.md) | River online learning, streaming anomaly detection, concept drift detection (ADWIN/DDM) | General |

### C. Machine Learning & Modeling (11 Skills)

Skills for supervised learning, feature interpretation, active learning, AutoML, ensemble learning, anomaly detection, causal ML, model monitoring, semi-supervised learning, and multi-task learning.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 17 | [scientific-ml-regression](scientific-ml-regression/SKILL.md) | Multi-target regression, model comparison, radar charts | 05, 12, 13 |
| 18 | [scientific-ml-classification](scientific-ml-classification/SKILL.md) | Classification ML, ROC, PR curves, confusion matrix, PDP, Volcano | 03, 05 |
| 19 | [scientific-feature-importance](scientific-feature-importance/SKILL.md) | Tree-based & permutation feature importance, PDP | 05, 12, 13 |
| 167 | [scientific-active-learning](scientific-active-learning/SKILL.md) | Uncertainty sampling, QBC, batch AL, active learning loop, stopping criteria | General |
| 168 | [scientific-automl](scientific-automl/SKILL.md) | Optuna HPO, multi-model AutoML, automated feature engineering, AutoML report | General |
| 169 | [scientific-ensemble-methods](scientific-ensemble-methods/SKILL.md) | XGBoost/LightGBM/CatBoost comparison, stacking OOF, voting, ensemble diversity | General |
| 175 | [scientific-anomaly-detection](scientific-anomaly-detection/SKILL.md) | Isolation Forest/LOF/OCSVM ensemble anomaly detection, autoencoder anomaly, SPC control charts | General |
| 176 | [scientific-causal-ml](scientific-causal-ml/SKILL.md) | DoWhy causal inference, EconML Double ML/Causal Forest, S/T/X-Learner meta-learners | General |
| 177 | [scientific-model-monitoring](scientific-model-monitoring/SKILL.md) | Data drift detection (KS/PSI/Wasserstein), performance degradation detection, A/B test statistics | General |
| 185 | [scientific-semi-supervised-learning](scientific-semi-supervised-learning/SKILL.md) | Self-Training, Label Propagation, Pseudo-Labeling quality evaluation, label efficiency | General |
| 186 | [scientific-multi-task-learning](scientific-multi-task-learning/SKILL.md) | Hard/Soft parameter sharing MTL, GradNorm dynamic task balancing, PCGrad | General |

### D. Experimental Design & Process Optimization (3 Skills)

Skills for experimental design, response surface optimization, and adaptive experimental design.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 20 | [scientific-doe](scientific-doe/SKILL.md) | Taguchi orthogonal arrays, CCD/Box-Behnken, ANOVA factor effects, Bayesian optimization | General |
| 21 | [scientific-process-optimization](scientific-process-optimization/SKILL.md) | Response surface methodology (ML-RSM), Pareto optimization, process windows | 12, 13 |
| 190 | [scientific-adaptive-experiments](scientific-adaptive-experiments/SKILL.md) | Thompson Sampling/UCB bandits, SPRT sequential testing, Bayesian adaptive dose finding | General |

### E. Signal, Spectral & Time Series (5 Skills)

Skills for waveform, frequency domain, and neuroelectrophysiology analysis.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 22 | [scientific-spectral-signal](scientific-spectral-signal/SKILL.md) | Spectral preprocessing, filtering, peak detection | 11 |
| 23 | [scientific-biosignal-processing](scientific-biosignal-processing/SKILL.md) | ECG R-wave/HRV, EEG band power/ERP, EMG burst, Poincaré | 08 |
| 24 | [scientific-time-series](scientific-time-series/SKILL.md) | STL decomposition, SARIMA forecasting, change point detection, FFT periodicity analysis, Granger causality | General |
| 67 | [scientific-neuroscience-electrophysiology](scientific-neuroscience-electrophysiology/SKILL.md) | SpikeInterface/Kilosort4 spike sorting, MNE EEG/ERP, NeuroKit2 HRV/EDA, brain functional connectivity | General |
| 178 | [scientific-time-series-forecasting](scientific-time-series-forecasting/SKILL.md) | Prophet/NeuralProphet ML forecasting, time series feature engineering, backtesting framework | General |

### F. Life Sciences & Omics (28 Skills)

Skills for bioinformatics, omics, network analysis, ontology, EBI databases, genomics, PPI, expression comparison, model organism DB, GEO expression profiles, parasite genomics, ArrayExpress expression archive, GTEx tissue expression, UniProt proteome, Reactome pathways, HGNC nomenclature, metabolic network, glycomics, and lipidomics.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 25 | [scientific-bioinformatics](scientific-bioinformatics/SKILL.md) | scRNA-seq, PPI network, bulk RNA-seq | 01, 04 |
| 26 | [scientific-metabolomics](scientific-metabolomics/SKILL.md) | PLS-DA/VIP scores, Pareto scaling, pathway enrichment | 07 |
| 27 | [scientific-sequence-analysis](scientific-sequence-analysis/SKILL.md) | RSCU/CAI codon analysis, alignment, phylogenetic trees, ORF/CpG islands | 09 |
| 28 | [scientific-multi-omics](scientific-multi-omics/SKILL.md) | CCA canonical correlation, SNF network fusion, pathway integration, multi-omics clustering | General |
| 29 | [scientific-network-analysis](scientific-network-analysis/SKILL.md) | Network construction, centrality, community, PSP path diagrams | 04, 07, 13 |
| 68 | [scientific-proteomics-mass-spectrometry](scientific-proteomics-mass-spectrometry/SKILL.md) | pyOpenMS LC-MS/MS, peptide ID, protein quantification, PTM, GNPS molecular networking | General |
| 69 | [scientific-gene-expression-transcriptomics](scientific-gene-expression-transcriptomics/SKILL.md) | GEO data retrieval, PyDESeq2 differential expression, GTEx tissue expression/eQTL, GSEA | General |
| 77 | [scientific-pathway-enrichment](scientific-pathway-enrichment/SKILL.md) | ORA/GSEA pathway enrichment analysis, KEGG/Reactome/GO/WikiPathways integration | General |
| 82 | [scientific-metabolomics-databases](scientific-metabolomics-databases/SKILL.md) | HMDB/MetaCyc/Metabolomics Workbench metabolite DB search, m/z identification | General |
| 88 | [scientific-human-protein-atlas](scientific-human-protein-atlas/SKILL.md) | HPA tissue/cell protein expression, RNA expression, cancer prognosis, subcellular localization | General |
| 90 | [scientific-genome-sequence-tools](scientific-genome-sequence-tools/SKILL.md) | Ensembl/dbSNP/BLAST/NCBI Nucleotide/GDC genome sequence analysis | General |
| 92 | [scientific-noncoding-rna](scientific-noncoding-rna/SKILL.md) | Rfam RNA families, RNAcentral ncRNA, covariance models, structure mapping | General |
| 99 | [scientific-ontology-enrichment](scientific-ontology-enrichment/SKILL.md) | EFO/OLS/Enrichr/UMLS ontology search, gene set enrichment, term mapping | General |
| 100 | [scientific-ebi-databases](scientific-ebi-databases/SKILL.md) | EBI Search/ENA/BioStudies/dbfetch/MetaboLights integrated data access | General |
| 108 | [scientific-ensembl-genomics](scientific-ensembl-genomics/SKILL.md) | Ensembl REST API genome analysis, VEP variant effect prediction, homology search, regulatory regions | General |
| 109 | [scientific-string-network-api](scientific-string-network-api/SKILL.md) | STRING v12/BioGRID/STITCH PPI network, chemical-protein interactions, topology analysis | General |
| 110 | [scientific-expression-comparison](scientific-expression-comparison/SKILL.md) | EBI Expression Atlas expression comparison, baseline/differential expression, cross-tissue heatmaps | General |
| 111 | [scientific-model-organism-db](scientific-model-organism-db/SKILL.md) | FlyBase/WormBase/ZFIN/RGD/MGI model organism databases, cross-species ortholog search | General |
| 127 | [scientific-geo-expression](scientific-geo-expression/SKILL.md) | GEO REST API expression profiles, matrix retrieval, differential expression analysis | General |
| 132 | [scientific-parasite-genomics](scientific-parasite-genomics/SKILL.md) | PlasmoDB/VectorBase/ToxoDB parasite genomics, drug target identification | General |
| 135 | [scientific-arrayexpress-expression](scientific-arrayexpress-expression/SKILL.md) | ArrayExpress/BioStudies REST API expression experiment search, SDRF metadata, data reanalysis | General |
| 137 | [scientific-gtex-tissue-expression](scientific-gtex-tissue-expression/SKILL.md) | GTEx Portal REST API v2 tissue-specific expression, eQTL, multi-tissue comparison | General |
| 141 | [scientific-uniprot-proteome](scientific-uniprot-proteome/SKILL.md) | UniProt REST API proteome search, ID mapping, domain/feature extraction | General |
| 144 | [scientific-reactome-pathways](scientific-reactome-pathways/SKILL.md) | Reactome Content Service pathway search, UniProt mapping, participant retrieval | General |
| 159 | [scientific-hgnc-nomenclature](scientific-hgnc-nomenclature/SKILL.md) | HGNC REST API gene nomenclature, official symbol search, alias resolution, gene families | General |
| 160 | [scientific-metabolomics-network](scientific-metabolomics-network/SKILL.md) | Metabolite correlation network construction, KEGG pathway graph, hub metabolites, enrichment | General |
| 161 | [scientific-glycomics](scientific-glycomics/SKILL.md) | GlyGen/GlyConnect glycan database integration, glycoprotein site search, MS fragmentation prediction | General |
| 162 | [scientific-lipidomics](scientific-lipidomics/SKILL.md) | LipidMAPS/SwissLipids lipid structure search, subclass classification, differential lipid analysis, lipid enrichment | General |

### G. Chemistry, Materials & Imaging (9 Skills)

Skills for chemical structure, materials characterization, image morphology analysis, computational materials science, ChEMBL assay mining, MD simulation, advanced imaging, deep chemistry, and STITCH chemical-protein networks.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 30 | [scientific-cheminformatics](scientific-cheminformatics/SKILL.md) | RDKit molecular descriptors, Tanimoto, structural alerts, Lipinski | 02, 05 |
| 31 | [scientific-materials-characterization](scientific-materials-characterization/SKILL.md) | Thornton-Anders SZM, XRD Scherrer, Tauc plot | 11, 12, 13 |
| 32 | [scientific-image-analysis](scientific-image-analysis/SKILL.md) | Otsu/Watershed segmentation, particle size distribution, GLCM texture, fluorescence synthesis | General |
| 70 | [scientific-computational-materials](scientific-computational-materials/SKILL.md) | pymatgen crystal structure, Materials Project, phase diagrams, band structure/DOS, VASP/QE I/O | General |
| 107 | [scientific-chembl-assay-mining](scientific-chembl-assay-mining/SKILL.md) | ChEMBL REST API assay mining, bioactivity, SAR analysis, selectivity profiling | General |
| 112 | [scientific-md-simulation](scientific-md-simulation/SKILL.md) | MDAnalysis/OpenFF molecular dynamics simulation, RMSD/RMSF/SASA/hydrogen bond analysis | General |
| 114 | [scientific-advanced-imaging](scientific-advanced-imaging/SKILL.md) | Cellpose segmentation, CellProfiler morphological profiling, napari 3D visualization | General |
| 115 | [scientific-deep-chemistry](scientific-deep-chemistry/SKILL.md) | DeepChem GCN/MPNN/AttentiveFP molecular property prediction, MoleculeNet, ChemBERTa | General |
| 154 | [scientific-stitch-chemical-network](scientific-stitch-chemical-network/SKILL.md) | STITCH chemical-protein interaction network, network pharmacology, polypharmacology | General |

### H. Clinical, Epidemiology & Meta-Science (7 Skills)

Skills for clinical trials, causal inference, meta-analysis, clinical trials analytics, biobank large-scale cohorts, and clinical standard terminology.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 33 | [scientific-survival-clinical](scientific-survival-clinical/SKILL.md) | Kaplan-Meier, Cox PH, power analysis, safety analysis | 03, 06 |
| 34 | [scientific-causal-inference](scientific-causal-inference/SKILL.md) | PSM propensity score, IPW, DID, RDD, DAG covariate selection, Rosenbaum sensitivity analysis | General |
| 35 | [scientific-meta-analysis](scientific-meta-analysis/SKILL.md) | Fixed/random effects models, forest/funnel plots, Egger test, subgroup analysis | General |
| 71 | [scientific-clinical-trials-analytics](scientific-clinical-trials-analytics/SKILL.md) | ClinicalTrials.gov API v2 search, competitive landscape, AE/outcome extraction | General |
| 85 | [scientific-clinical-reporting](scientific-clinical-reporting/SKILL.md) | SOAP notes, biomarker reports, pharmacogenomics, FHIR JSON | General |
| 151 | [scientific-biobank-cohort](scientific-biobank-cohort/SKILL.md) | UK Biobank/BBJ/All of Us large-scale cohort, GWAS summary statistics, PheWAS | General |
| 166 | [scientific-clinical-standards](scientific-clinical-standards/SKILL.md) | LOINC/ICD-10/ICD-11 clinical standard code search, FHIR R4 mapping, terminology interoperability | General |

### I. Deep Research & Literature Search (4 Skills)

Skills for iterative deep research of scientific literature, multi-DB literature search, preprint cross-search, and Semantic Scholar academic graph.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 36 | [scientific-deep-research](scientific-deep-research/SKILL.md) | SHIKIGAMI-compliant Think→Search→Evaluate→Synthesize iterative cycle, academic DB search, evidence hierarchy evaluation, source tracking, cross-validation, hallucination prevention | General |
| 78 | [scientific-literature-search](scientific-literature-search/SKILL.md) | PubMed/Semantic Scholar/OpenAlex/EuropePMC/CrossRef multi-DB search, citation networks | General |
| 97 | [scientific-preprint-archive](scientific-preprint-archive/SKILL.md) | bioRxiv/medRxiv/arXiv/PMC/CORE/Zenodo/OpenAIRE/Unpaywall preprint & OA cross-search | General |
| 136 | [scientific-semantic-scholar](scientific-semantic-scholar/SKILL.md) | Semantic Scholar Academic Graph API paper search, citation graph, author profiles, TLDR | General |

### J. Drug Discovery & Pharmacology (9 Skills)

Skills for drug discovery target evaluation, pharmacokinetics, repositioning, pharmacological targets, compound screening, NCI-60 screening, DrugBank resources, and Pharos targets.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 37 | [scientific-drug-target-profiling](scientific-drug-target-profiling/SKILL.md) | 9-path target profiling, TDL classification, druggability assessment, competitive landscape | General |
| 38 | [scientific-admet-pharmacokinetics](scientific-admet-pharmacokinetics/SKILL.md) | 5-stage ADMET pipeline, Lipinski/Veber rules, CYP prediction, PK modeling | General |
| 39 | [scientific-drug-repurposing](scientific-drug-repurposing/SKILL.md) | 7-strategy drug repositioning, network proximity analysis, multi-criteria candidate scoring | General |
| 83 | [scientific-molecular-docking](scientific-molecular-docking/SKILL.md) | AutoDock Vina/DiffDock molecular docking, virtual screening | General |
| 89 | [scientific-pharmacology-targets](scientific-pharmacology-targets/SKILL.md) | BindingDB/GPCRdb/GtoPdb/BRENDA/Pharos pharmacological target profiling | General |
| 94 | [scientific-compound-screening](scientific-compound-screening/SKILL.md) | ZINC compound library search, virtual screening preprocessing | General |
| 120 | [scientific-nci60-screening](scientific-nci60-screening/SKILL.md) | NCI-60/CellMiner/DepMap cancer cell line drug response screening | General |
| 146 | [scientific-drugbank-resources](scientific-drugbank-resources/SKILL.md) | DrugBank API drug information, pharmacological MOA, target proteins, drug interactions | General |
| 156 | [scientific-pharos-targets](scientific-pharos-targets/SKILL.md) | Pharos/TCRD IDG target TDL classification, disease associations, ligand search | General |

### K. Structural Biology & Protein Engineering (7 Skills)

Skills for protein structure analysis, design, PPI networks, domain analysis, structural proteomics, AlphaFold DB structure prediction, and RCSB PDB structure search.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 40 | [scientific-protein-structure-analysis](scientific-protein-structure-analysis/SKILL.md) | PDB/AlphaFold structure search, quality evaluation (pLDDT/R-factor), binding site detection | General |
| 41 | [scientific-protein-design](scientific-protein-design/SKILL.md) | ESM-2 mutation scan, RFdiffusion/ProteinMPNN de novo design, binder/enzyme design | General |
| 79 | [scientific-protein-interaction-network](scientific-protein-interaction-network/SKILL.md) | STRING/IntAct/STITCH PPI network, topology analysis, community detection | General |
| 86 | [scientific-protein-domain-family](scientific-protein-domain-family/SKILL.md) | InterPro/InterProScan domain prediction, family classification, architecture visualization | General |
| 93 | [scientific-structural-proteomics](scientific-structural-proteomics/SKILL.md) | EMDB/PDBe/Proteins API/Complex Portal/DeepGO/EVE structural proteomics | General |
| 134 | [scientific-alphafold-structures](scientific-alphafold-structures/SKILL.md) | AlphaFold DB REST API structure prediction retrieval, pLDDT confidence, PAE analysis | General |
| 142 | [scientific-rcsb-pdb-search](scientific-rcsb-pdb-search/SKILL.md) | RCSB PDB Search/Data API structure search, metadata, ligand information | General |

### L. Precision Medicine & Clinical Decision (6 Skills)

Skills for variant interpretation, evidence-based clinical judgment, CIViC clinical evidence, gnomAD variants, and ClinGen curation.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 42 | [scientific-variant-interpretation](scientific-variant-interpretation/SKILL.md) | ACMG/AMP 28 criteria, pharmacogenomics (CPIC), OncoKB somatic mutation levels | General |
| 43 | [scientific-clinical-decision-support](scientific-clinical-decision-support/SKILL.md) | GRADE evidence framework, precision oncology workflow, clinical trial matching | General |
| 80 | [scientific-variant-effect-prediction](scientific-variant-effect-prediction/SKILL.md) | AlphaMissense/CADD/SpliceAI variant effect prediction, consensus pathogenicity assessment | General |
| 147 | [scientific-civic-evidence](scientific-civic-evidence/SKILL.md) | CIViC REST API cancer variant clinical interpretation, evidence, assertions | General |
| 148 | [scientific-gnomad-variants](scientific-gnomad-variants/SKILL.md) | gnomAD GraphQL population allele frequencies, gene constraint (pLI/LOEUF), region queries | General |
| 157 | [scientific-clingen-curation](scientific-clingen-curation/SKILL.md) | ClinGen gene-disease validity, dosage sensitivity, clinical actionability | General |

### M. Lab Automation & Data Management (3 Skills)

Skills for lab experiment automation, data management, and CRISPR gRNA design.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 44 | [scientific-lab-automation](scientific-lab-automation/SKILL.md) | PyLabRobot/Opentrons protocol, SOP templates, ELN/LIMS integration, QC verification | General |
| 72 | [scientific-lab-data-management](scientific-lab-data-management/SKILL.md) | Benchling ELN/DNA design, DNAnexus PaaS, OMERO bioimaging, Protocols.io | General |
| 164 | [scientific-crispr-design](scientific-crispr-design/SKILL.md) | CRISPR gRNA design, Cas9/Cas12a PAM search, off-target scoring, sgRNA library construction | General |

### N. Scientific Presentation & Schematics (4 Skills)

Skills for designing conference slides, posters, scientific schematics, interactive dashboards, and reproducible reports.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 45 | [scientific-presentation-design](scientific-presentation-design/SKILL.md) | 15-slide structure template, tikzposter, matplotlib workflow diagrams, accessibility | General |
| 73 | [scientific-scientific-schematics](scientific-scientific-schematics/SKILL.md) | CONSORT flow diagrams, NN architecture diagrams, pathway diagrams, TikZ/SVG | General |
| 174 | [scientific-interactive-dashboard](scientific-interactive-dashboard/SKILL.md) | Streamlit/Dash/Panel scientific data dashboards, parameter exploration UI, real-time analysis | General |
| 182 | [scientific-reproducible-reporting](scientific-reproducible-reporting/SKILL.md) | Quarto scientific documents, Jupyter Book multi-chapter composition, Papermill parametric execution, nbconvert auto-conversion | General |

### O. Research Planning, Grants & Regulation (3 Skills)

Skills for grant applications, research methodology, and regulatory science design.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 46 | [scientific-grant-writing](scientific-grant-writing/SKILL.md) | NIH Specific Aims template, JSPS KAKENHI, budget planning, Budget Justification | General |
| 47 | [scientific-research-methodology](scientific-research-methodology/SKILL.md) | SCAMPER/TRIZ brainstorming, research design matrix, FINER criteria, IRB ethics check | General |
| 74 | [scientific-regulatory-science](scientific-regulatory-science/SKILL.md) | FDA Orange Book, medical device 510(k), ISO 13485 QMS, CAPA, USPTO patent search | General |

### P. Pharmacovigilance & Pharmacogenomics (4 Skills)

Skills for post-market drug safety surveillance, pharmacogenomics, PharmGKB clinical annotations, and clinical pharmacology modeling signal detection and quantitative evaluation.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 48 | [scientific-pharmacovigilance](scientific-pharmacovigilance/SKILL.md) | FAERS disproportionality analysis (PRR/ROR/IC/EBGM), MedDRA hierarchy, temporal trends, Naranjo causality assessment | General |
| 75 | [scientific-pharmacogenomics](scientific-pharmacogenomics/SKILL.md) | PharmGKB/CPIC guidelines, star alleles, metabolizer types, FDA PGx biomarkers | General |
| 138 | [scientific-pharmgkb-pgx](scientific-pharmgkb-pgx/SKILL.md) | PharmGKB REST API clinical annotations, drug-gene associations, dosing guidelines | General |
| 165 | [scientific-clinical-pharmacology](scientific-clinical-pharmacology/SKILL.md) | PopPK NLME, PBPK simulation, TDM dose optimization, Emax PD modeling | General |

### Q. Oncology & Disease Research (10 Skills)

Skills for precision oncology, disease-gene association research, cancer genomics, rare disease genetics, cell line resources, ICGC cancer genome data, Open Targets genetics, DepMap dependencies, Monarch ontology, and GDC portal.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 49 | [scientific-precision-oncology](scientific-precision-oncology/SKILL.md) | CIViC/OncoKB/cBioPortal integration, TMB/MSI assessment, AMP Tiering, MTB report | General |
| 50 | [scientific-disease-research](scientific-disease-research/SKILL.md) | GWAS Catalog, DisGeNET GDA, Orphanet/OMIM/HPO phenotype matching, PRS calculation | General |
| 81 | [scientific-cancer-genomics](scientific-cancer-genomics/SKILL.md) | COSMIC/cBioPortal/DepMap cancer genomics, mutational signature analysis | General |
| 87 | [scientific-rare-disease-genetics](scientific-rare-disease-genetics/SKILL.md) | OMIM/Orphanet/DisGeNET/IMPC rare disease genetics, integrated analysis | General |
| 101 | [scientific-cell-line-resources](scientific-cell-line-resources/SKILL.md) | Cellosaurus cell line search, STR profile verification, contamination detection | General |
| 140 | [scientific-icgc-cancer-data](scientific-icgc-cancer-data/SKILL.md) | ICGC DCC API international cancer genome data, somatic mutations, cancer type statistics | General |
| 143 | [scientific-opentargets-genetics](scientific-opentargets-genetics/SKILL.md) | Open Targets Platform GraphQL target-disease associations, drug evidence, L2G | General |
| 145 | [scientific-depmap-dependencies](scientific-depmap-dependencies/SKILL.md) | DepMap Portal CRISPR/RNAi gene dependencies, drug sensitivity | General |
| 149 | [scientific-monarch-ontology](scientific-monarch-ontology/SKILL.md) | Monarch Initiative disease-gene-phenotype ontology, HPO, entity search | General |
| 150 | [scientific-gdc-portal](scientific-gdc-portal/SKILL.md) | NCI Genomic Data Commons REST API, project/case/SSM search | General |

### R. Quantum & Advanced Computing (11 Skills)

Skills for next-generation computational methods including quantum computing, GNN, Bayesian statistics, XAI, deep learning, healthcare AI, reinforcement learning, transfer learning, uncertainty quantification, federated learning, and NAS.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 51 | [scientific-quantum-computing](scientific-quantum-computing/SKILL.md) | Qiskit/Cirq/PennyLane VQE, QAOA, quantum ML, QuTiP quantum dynamics | General |
| 52 | [scientific-graph-neural-networks](scientific-graph-neural-networks/SKILL.md) | PyG GCN/GAT/GIN, TorchDrug molecular property prediction, knowledge graph reasoning, scaffold split | General |
| 53 | [scientific-bayesian-statistics](scientific-bayesian-statistics/SKILL.md) | PyMC/Stan hierarchical Bayes, MCMC diagnostics, PPC, WAIC/LOO-CV model comparison, Bayesian optimization | General |
| 54 | [scientific-explainable-ai](scientific-explainable-ai/SKILL.md) | SHAP/LIME/Captum feature attribution, counterfactual explanations, fairness audit, DeepSHAP | General |
| 55 | [scientific-deep-learning](scientific-deep-learning/SKILL.md) | Lightning/timm/Transformers, CNN/ViT/BERT fine-tune, Optuna HPO, ONNX export | General |
| 96 | [scientific-healthcare-ai](scientific-healthcare-ai/SKILL.md) | PyHealth clinical ML pipeline, flow cytometry, EHR processing | General |
| 104 | [scientific-reinforcement-learning](scientific-reinforcement-learning/SKILL.md) | Stable-Baselines3/PufferLib RL agent training, molecular design/experiment optimization RL | General |
| 170 | [scientific-transfer-learning](scientific-transfer-learning/SKILL.md) | Vision/NLP fine-tuning, few-shot, knowledge distillation, domain adaptation | General |
| 171 | [scientific-uncertainty-quantification](scientific-uncertainty-quantification/SKILL.md) | Conformal Prediction, MC Dropout, deep ensemble, calibration, ECE | General |
| 183 | [scientific-federated-learning](scientific-federated-learning/SKILL.md) | Flower FL pipeline, FedAvg/FedProx aggregation, DP-SGD differential privacy, non-IID partitioning | General |
| 184 | [scientific-neural-architecture-search](scientific-neural-architecture-search/SKILL.md) | Optuna NAS, Pareto multi-objective search (accuracy vs size), search space definition, pruning | General |

### S. Medical Imaging (2 Skills)

Skills for medical image analysis, segmentation, and radiology AI for DICOM, WSI, and other medical images.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 56 | [scientific-medical-imaging](scientific-medical-imaging/SKILL.md) | DICOM/NIfTI processing, MONAI U-Net/SwinUNETR, WSI patch extraction, radiomics, 3D visualization | General |
| 189 | [scientific-radiology-ai](scientific-radiology-ai/SKILL.md) | MONAI CADe/CADx pipeline, CT/MRI classification, Grad-CAM explainability, structured reports | General |

### T. Single-Cell, Spatial & Epigenomics (13 Skills)

Skills for scRNA-seq, spatial transcriptomics, epigenomics, regulatory genomics, perturbation analysis, scVI integration, scATAC-seq, GPU single-cell, ENCODE/SCREEN, Human Cell Atlas, advanced Squidpy spatial analysis, spatial multi-omics, and CELLxGENE Census analysis pipelines.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 57 | [scientific-single-cell-genomics](scientific-single-cell-genomics/SKILL.md) | scRNA-seq QC, Scanpy Leiden clustering, DEG, RNA velocity, CellChat cell-cell communication | General |
| 58 | [scientific-spatial-transcriptomics](scientific-spatial-transcriptomics/SKILL.md) | Visium/MERFISH preprocessing, Squidpy SVG detection, spatial domains, cell2location deconvolution | General |
| 76 | [scientific-epigenomics-chromatin](scientific-epigenomics-chromatin/SKILL.md) | ChIP-seq MACS2/3, ATAC-seq, WGBS DMR, ChromHMM, Hi-C TAD, motif enrichment | General |
| 102 | [scientific-regulatory-genomics](scientific-regulatory-genomics/SKILL.md) | RegulomeDB/ReMap/4DN regulatory region variants, 3D genome structure analysis | General |
| 113 | [scientific-perturbation-analysis](scientific-perturbation-analysis/SKILL.md) | pertpy/Augur/scIB perturbation analysis, CRISPR screens, scGen perturbation prediction, integration benchmarks | General |
| 116 | [scientific-scvi-integration](scientific-scvi-integration/SKILL.md) | scVI/scANVI/totalVI/SOLO single-cell integration, semi-supervised annotation, CITE-seq | General |
| 123 | [scientific-scatac-signac](scientific-scatac-signac/SKILL.md) | Signac/SnapATAC2 scATAC-seq, motif analysis, Gene Activity, RNA+ATAC multimodal integration | General |
| 124 | [scientific-gpu-singlecell](scientific-gpu-singlecell/SKILL.md) | rapids-singlecell/cuML/cuGraph GPU-accelerated single-cell analysis | General |
| 125 | [scientific-encode-screen](scientific-encode-screen/SKILL.md) | ENCODE REST API experiment/file search, SCREEN cCRE, ChIP-Atlas enrichment | General |
| 126 | [scientific-human-cell-atlas](scientific-human-cell-atlas/SKILL.md) | HCA Data Portal projects/files, CELLxGENE Census large-scale atlas | General |
| 131 | [scientific-squidpy-advanced](scientific-squidpy-advanced/SKILL.md) | Squidpy spatial autocorrelation, co-occurrence analysis, neighborhood enrichment, niche identification | General |
| 152 | [scientific-spatial-multiomics](scientific-spatial-multiomics/SKILL.md) | MERFISH/CODEX spatial multi-omics integration, co-detection analysis, spatial community detection | General |
| 155 | [scientific-cellxgene-census](scientific-cellxgene-census/SKILL.md) | CELLxGENE Census API large-scale single-cell atlas, cell type distribution, gene expression | General |

### U. Immune & Infectious Disease (2 Skills)

Skills for immunoinformatics and pathogen genomics analysis pipelines.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 59 | [scientific-immunoinformatics](scientific-immunoinformatics/SKILL.md) | MHC-I/II binding prediction, B-cell epitopes, TCR/BCR repertoire diversity, antibody CDR analysis, vaccine candidate ranking | General |
| 60 | [scientific-infectious-disease](scientific-infectious-disease/SKILL.md) | Pathogen WGS QC, AMR gene detection, MLST typing, phylogenetic analysis (IQ-TREE), SIR/SEIR mathematical models | General |

### V. Microbiome & Environment (9 Skills)

Skills for microbiome analysis, environmental/ecosystem modeling, phylogenetics, rRNA taxonomy, plant biology, marine ecology, environmental geodata, paleobiology, and MAG reconstruction.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 61 | [scientific-microbiome-metagenomics](scientific-microbiome-metagenomics/SKILL.md) | DADA2 ASV pipeline, MetaPhlAn/Kraken2, α/β diversity, ANCOM-BC differential abundance, HUMAnN functional profiling | General |
| 62 | [scientific-environmental-ecology](scientific-environmental-ecology/SKILL.md) | SDM (MaxEnt/RF/GBM), biodiversity indices, community ordination (NMDS/CCA), conservation priority ranking | General |
| 103 | [scientific-phylogenetics](scientific-phylogenetics/SKILL.md) | ETE3/scikit-bio phylogenetic tree construction, Faith's PD, UniFrac, divergence time estimation, ancestral sequence reconstruction | General |
| 118 | [scientific-rrna-taxonomy](scientific-rrna-taxonomy/SKILL.md) | SILVA/Greengenes2/MGnify rRNA reference, taxonomy, consensus classification | General |
| 121 | [scientific-plant-biology](scientific-plant-biology/SKILL.md) | Plant Reactome/TAIR/Ensembl Plants plant metabolic pathways, cross-species comparison | General |
| 122 | [scientific-marine-ecology](scientific-marine-ecology/SKILL.md) | OBIS/WoRMS/GBIF/FishBase marine biodiversity, distribution analysis | General |
| 128 | [scientific-environmental-geodata](scientific-environmental-geodata/SKILL.md) | SoilGrids/WorldClim environmental geospatial data, species distribution model environmental variables | General |
| 129 | [scientific-paleobiology](scientific-paleobiology/SKILL.md) | PBDB fossil occurrence records, taxon search, geological age diversity curves | General |
| 163 | [scientific-metagenome-assembled-genomes](scientific-metagenome-assembled-genomes/SKILL.md) | MetaBAT2/CONCOCT binning, CheckM2 quality assessment, GTDB-Tk classification, MAG pipeline | General |

### W. Systems Biology (4 Skills)

Skills for SBML dynamic simulation, metabolic flux analysis, and gene regulatory network inference.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 63 | [scientific-systems-biology](scientific-systems-biology/SKILL.md) | SBML/RoadRunner simulation, FBA/pFBA (cobrapy), GRN inference (GENIE3), Sobol sensitivity analysis | General |
| 95 | [scientific-metabolic-modeling](scientific-metabolic-modeling/SKILL.md) | BiGG Models/BioModels genome-scale metabolic models, reactions, metabolite search | General |
| 130 | [scientific-metabolic-atlas](scientific-metabolic-atlas/SKILL.md) | Metabolic Atlas/Human-GEM metabolic reactions, metabolite search, network analysis | General |
| 153 | [scientific-metabolic-flux](scientific-metabolic-flux/SKILL.md) | ¹³C/¹⁵N stable isotope metabolic flux analysis, EMU modeling, MID fitting | General |

### X. Epidemiology & Public Health (3 Skills)

Skills for epidemiological risk measure calculation, standardization, spatial epidemiology, DAG confounding analysis, public health data access, and environmental toxicology.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 64 | [scientific-epidemiology-public-health](scientific-epidemiology-public-health/SKILL.md) | RR/OR/RD/NNT/AF risk measures, direct/indirect age standardization, LISA/Getis-Ord spatial clustering, DAG backdoor criteria | General |
| 98 | [scientific-public-health-data](scientific-public-health-data/SKILL.md) | NHANES/MedlinePlus/RxNorm/ODPHP public health data access, health disparities API | General |
| 117 | [scientific-toxicology-env](scientific-toxicology-env/SKILL.md) | CTD/Tox21/ToxCast/T3DB/EPA IRIS environmental toxicology, chemical health impact analysis | General |

### Y. Population Genetics (2 Skills)

Skills for population structure estimation, differentiation metrics, and natural selection detection, plus GWAS Catalog.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 65 | [scientific-population-genetics](scientific-population-genetics/SKILL.md) | PLINK2 QC, HWE testing, PCA/ADMIXTURE, Weir-Cockerham Fst, iHS/Tajima's D selection scan | General |
| 133 | [scientific-gwas-catalog](scientific-gwas-catalog/SKILL.md) | NHGRI-EBI GWAS Catalog REST API association analysis, study search, PheWAS | General |

### Z. Scientific Text Mining (3 Skills)

Skills for information extraction from scientific literature, knowledge graph construction, topic modeling, biomedical NER, and clinical NLP.

| # | Skill | Description | Ref. Exp |
|---|---|---|---|
| 66 | [scientific-text-mining-nlp](scientific-text-mining-nlp/SKILL.md) | BioBERT/SciSpaCy NER, relation extraction, knowledge graph construction (Louvain), BERTopic topic modeling, citation network analysis | General |
| 106 | [scientific-biomedical-pubtator](scientific-biomedical-pubtator/SKILL.md) | PubTator3 biomedical NER, entity relation extraction, knowledge graph construction | General |
| 158 | [scientific-clinical-nlp](scientific-clinical-nlp/SKILL.md) | MedSpaCy/scispaCy clinical text NER, negation detection, section classification, UMLS linking | General |

---

## Installation

```bash
# One-command install with npx
npx @nahisaho/satori init

# Or global install
npm install -g @nahisaho/satori
satori init
```

`.github/skills/` is copied to the current directory and becomes immediately available in Copilot Agent Mode.

### CLI Commands

| Command | Description |
|---|---|
| `satori init [--force] [--dry-run]` | Install `.github/skills/` to the current directory |
| `satori pipeline suggest` | Interactively recommend pipelines by keyword |
| `satori pipeline list` | List all 26 domain pipelines |
| `satori validate [--verbose]` | Auto-validate all 190 SKILL.md files |
| `satori stats` | Display skill count, TU coverage, and other statistics |
| `satori help` | Display help |
| `satori --version` | Display version |

---

## Usage

### Using with GitHub Copilot Agent Mode / Copilot CLI

Skills are placed in `.github/skills/`, so Copilot automatically detects them.
Based on prompt context, the relevant Skill's `SKILL.md` is injected into the agent.

```
# Example: Requesting correlation analysis auto-loads scientific-eda-correlation
> Create a correlation heatmap of ZnO thin film process parameters and film properties

# Example: Requesting a classification model auto-loads scientific-ml-classification
> Compare Random Forest and SVM ROC for cancer gene expression data

# Example: Requesting DOE auto-loads scientific-doe
> Create a Taguchi L9 orthogonal array for 3 factors and draw main effects plots

# Example: Requesting meta-analysis auto-loads scientific-meta-analysis
> Create a forest plot with a random effects model from 5 RCT papers
```

### Directory Structure

```
.github/skills/
├── README.md
│
│── [A] Foundation & Workflow
│   ├── scientific-pipeline-scaffold/
│   ├── scientific-data-preprocessing/
│   ├── scientific-data-simulation/
│   ├── scientific-publication-figures/
│   ├── scientific-academic-writing/
│   │   └── assets/   ← 7 journal-specific templates
│   ├── scientific-hypothesis-pipeline/
│   ├── scientific-critical-review/
│   ├── scientific-supplementary-generator/
│   ├── scientific-latex-formatter/
│   ├── scientific-citation-checker/
│   ├── scientific-peer-review-response/
│   ├── scientific-revision-tracker/
│   ├── scientific-paper-quality/
│   ├── scientific-systematic-review/
│   ├── scientific-biothings-idmapping/
│   ├── scientific-data-submission/
│   └── scientific-crossref-metadata/
│
│── [B] Statistics & Exploratory Analysis
│   ├── scientific-eda-correlation/
│   ├── scientific-statistical-testing/
│   ├── scientific-pca-tsne/
│   ├── scientific-symbolic-mathematics/
│   ├── scientific-missing-data-analysis/
│   ├── scientific-advanced-visualization/
│   ├── scientific-data-profiling/
│   ├── scientific-geospatial-analysis/
│   ├── scientific-network-visualization/
│   ├── scientific-statistical-simulation/
│   └── scientific-streaming-analytics/
│
│── [C] Machine Learning & Modeling
│   ├── scientific-ml-regression/
│   ├── scientific-ml-classification/
│   ├── scientific-feature-importance/
│   ├── scientific-active-learning/
│   ├── scientific-automl/
│   ├── scientific-ensemble-methods/
│   ├── scientific-anomaly-detection/
│   ├── scientific-causal-ml/
│   ├── scientific-model-monitoring/
│   ├── scientific-semi-supervised-learning/
│   └── scientific-multi-task-learning/
│
│── [D] Experimental Design & Process Optimization
│   ├── scientific-doe/
│   ├── scientific-process-optimization/
│   └── scientific-adaptive-experiments/
│
│── [E] Signal, Spectral & Time Series
│   ├── scientific-spectral-signal/
│   ├── scientific-biosignal-processing/
│   ├── scientific-time-series/
│   ├── scientific-neuroscience-electrophysiology/
│   └── scientific-time-series-forecasting/
│
│── [F] Life Sciences & Omics
│   ├── scientific-bioinformatics/
│   ├── scientific-metabolomics/
│   ├── scientific-sequence-analysis/
│   ├── scientific-multi-omics/
│   ├── scientific-network-analysis/
│   ├── scientific-proteomics-mass-spectrometry/
│   ├── scientific-gene-expression-transcriptomics/
│   ├── scientific-pathway-enrichment/
│   ├── scientific-metabolomics-databases/
│   ├── scientific-human-protein-atlas/
│   ├── scientific-genome-sequence-tools/
│   ├── scientific-noncoding-rna/
│   ├── scientific-ontology-enrichment/
│   ├── scientific-ebi-databases/
│   ├── scientific-ensembl-genomics/
│   ├── scientific-string-network-api/
│   ├── scientific-expression-comparison/
│   ├── scientific-model-organism-db/
│   ├── scientific-geo-expression/
│   ├── scientific-parasite-genomics/
│   ├── scientific-arrayexpress-expression/
│   ├── scientific-gtex-tissue-expression/
│   ├── scientific-uniprot-proteome/
│   ├── scientific-reactome-pathways/
│   ├── scientific-hgnc-nomenclature/
│   ├── scientific-metabolomics-network/
│   ├── scientific-glycomics/
│   └── scientific-lipidomics/
│
│── [G] Chemistry, Materials & Imaging
│   ├── scientific-cheminformatics/
│   ├── scientific-materials-characterization/
│   ├── scientific-image-analysis/
│   ├── scientific-computational-materials/
│   ├── scientific-chembl-assay-mining/
│   ├── scientific-md-simulation/
│   ├── scientific-advanced-imaging/
│   ├── scientific-deep-chemistry/
│   └── scientific-stitch-chemical-network/
│
├── [H] Clinical, Epidemiology & Meta-Science
│   ├── scientific-survival-clinical/
│   ├── scientific-causal-inference/
│   ├── scientific-meta-analysis/
│   ├── scientific-clinical-trials-analytics/
│   ├── scientific-clinical-reporting/
│   ├── scientific-biobank-cohort/
│   └── scientific-clinical-standards/
│
├── [I] Deep Research & Literature Search
│   ├── scientific-deep-research/
│   ├── scientific-literature-search/
│   ├── scientific-preprint-archive/
│   └── scientific-semantic-scholar/
│
├── [J] Drug Discovery & Pharmacology
│   ├── scientific-drug-target-profiling/
│   ├── scientific-admet-pharmacokinetics/
│   ├── scientific-drug-repurposing/
│   ├── scientific-molecular-docking/
│   ├── scientific-pharmacology-targets/
│   ├── scientific-compound-screening/
│   ├── scientific-nci60-screening/
│   ├── scientific-drugbank-resources/
│   └── scientific-pharos-targets/
│
├── [K] Structural Biology & Protein Engineering
│   ├── scientific-protein-structure-analysis/
│   ├── scientific-protein-design/
│   ├── scientific-protein-interaction-network/
│   ├── scientific-protein-domain-family/
│   ├── scientific-structural-proteomics/
│   ├── scientific-alphafold-structures/
│   └── scientific-rcsb-pdb-search/
│
├── [L] Precision Medicine & Clinical Decision
│   ├── scientific-variant-interpretation/
│   ├── scientific-clinical-decision-support/
│   ├── scientific-variant-effect-prediction/
│   ├── scientific-civic-evidence/
│   ├── scientific-gnomad-variants/
│   └── scientific-clingen-curation/
│
├── [M] Lab Automation & Data Management
│   ├── scientific-lab-automation/
│   ├── scientific-lab-data-management/
│   └── scientific-crispr-design/
│
├── [N] Scientific Presentation & Schematics
│   ├── scientific-presentation-design/
│   ├── scientific-scientific-schematics/
│   ├── scientific-interactive-dashboard/
│   └── scientific-reproducible-reporting/
│
└── [O] Research Planning, Grants & Regulation
    ├── scientific-grant-writing/
    ├── scientific-research-methodology/
    └── scientific-regulatory-science/
│
├── [P] Pharmacovigilance & Pharmacogenomics
│   ├── scientific-pharmacovigilance/
│   ├── scientific-pharmacogenomics/
│   ├── scientific-pharmgkb-pgx/
│   └── scientific-clinical-pharmacology/
│
├── [Q] Oncology & Disease Research
│   ├── scientific-precision-oncology/
│   ├── scientific-disease-research/
│   ├── scientific-cancer-genomics/
│   ├── scientific-rare-disease-genetics/
│   ├── scientific-cell-line-resources/
│   ├── scientific-icgc-cancer-data/
│   ├── scientific-opentargets-genetics/
│   ├── scientific-depmap-dependencies/
│   ├── scientific-monarch-ontology/
│   └── scientific-gdc-portal/
│
├── [R] Quantum & Advanced Computing
│   ├── scientific-quantum-computing/
│   ├── scientific-graph-neural-networks/
│   ├── scientific-bayesian-statistics/
│   ├── scientific-explainable-ai/
│   ├── scientific-deep-learning/
│   ├── scientific-healthcare-ai/
│   ├── scientific-reinforcement-learning/
│   ├── scientific-transfer-learning/
│   ├── scientific-uncertainty-quantification/
│   ├── scientific-federated-learning/
│   └── scientific-neural-architecture-search/
│
├── [S] Medical Imaging
│   ├── scientific-medical-imaging/
│   └── scientific-radiology-ai/
│
│── [T] Single-Cell, Spatial & Epigenomics
│   ├── scientific-single-cell-genomics/
│   ├── scientific-spatial-transcriptomics/
│   ├── scientific-epigenomics-chromatin/
│   ├── scientific-regulatory-genomics/
│   ├── scientific-perturbation-analysis/
│   ├── scientific-scvi-integration/
│   ├── scientific-scatac-signac/
│   ├── scientific-gpu-singlecell/
│   ├── scientific-encode-screen/
│   ├── scientific-human-cell-atlas/
│   ├── scientific-squidpy-advanced/
│   ├── scientific-spatial-multiomics/
│   └── scientific-cellxgene-census/
│
│── [U] Immune & Infectious Disease
│   ├── scientific-immunoinformatics/
│   └── scientific-infectious-disease/
│
│── [V] Microbiome & Environment
│   ├── scientific-microbiome-metagenomics/
│   ├── scientific-environmental-ecology/
│   ├── scientific-phylogenetics/
│   ├── scientific-rrna-taxonomy/
│   ├── scientific-plant-biology/
│   ├── scientific-marine-ecology/
│   ├── scientific-environmental-geodata/
│   ├── scientific-paleobiology/
│   └── scientific-metagenome-assembled-genomes/
│
│── [W] Systems Biology
│   ├── scientific-systems-biology/
│   ├── scientific-metabolic-modeling/
│   ├── scientific-metabolic-atlas/
│   └── scientific-metabolic-flux/
│
│── [X] Epidemiology & Public Health
│   ├── scientific-epidemiology-public-health/
│   ├── scientific-public-health-data/
│   └── scientific-toxicology-env/
│
│── [Y] Population Genetics
│   ├── scientific-population-genetics/
│   └── scientific-gwas-catalog/
│
└── [Z] Scientific Text Mining
    ├── scientific-text-mining-nlp/
    ├── scientific-biomedical-pubtator/
    └── scientific-clinical-nlp/
```

> Note: In the actual file system, all skill directories are placed flat under `.github/skills/`. The category groupings above are logical classifications.

---

## Development

### Prerequisites

- Node.js >= 18

### Setup

```bash
git clone https://github.com/nahisaho/satori.git
cd satori
npm install
```

### Testing

```bash
npm test                  # Run all tests (1,359 tests)
npm run test:unit         # CLI unit tests only (24 tests)
npm run test:validation   # SKILL.md validation tests only (1,335 tests)
npm run test:coverage     # Run with coverage
npm run test:watch        # Watch mode
```

### Linting & Formatting

```bash
npm run lint              # Biome lint check
npm run lint:fix          # Auto-fix
npm run format            # Code formatting
```

### CI

GitHub Actions runs the following 4 jobs on push/PR to the `main` branch:

- **test**: All tests across Node.js 18/20/22 matrix
- **validate-skills**: All SKILL.md format validation
- **cli-smoke**: CLI command smoke tests
- **lint**: Code quality check with Biome

---

## References

- [GitHub Copilot Agent Skills Documentation](https://docs.github.com/en/copilot/concepts/agents/about-agent-skills)
- [Agent Skills Open Standard](https://github.com/agentskills/agentskills)
