#!/usr/bin/env node

const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const COMMAND = process.argv[2];
const SUBCOMMAND = process.argv[3];
const FLAGS = process.argv.slice(3);

const PACKAGE_ROOT = path.resolve(__dirname, '..');
const SOURCE_DIR = path.join(PACKAGE_ROOT, 'src', '.github');
const CUSTOM_PIPELINES_PATH = path.join(process.env.HOME || os.homedir(), '.satori', 'custom-pipelines.json');

function copyDirSync(src, dest) {
  fs.mkdirSync(dest, { recursive: true });
  for (const entry of fs.readdirSync(src, { withFileTypes: true })) {
    const srcPath = path.join(src, entry.name);
    const destPath = path.join(dest, entry.name);
    if (entry.isDirectory()) {
      copyDirSync(srcPath, destPath);
    } else {
      fs.copyFileSync(srcPath, destPath);
    }
  }
}

function countFiles(dir) {
  let count = 0;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    if (entry.isDirectory()) {
      count += countFiles(path.join(dir, entry.name));
    } else {
      count++;
    }
  }
  return count;
}

function init() {
  const dryRun = FLAGS.includes('--dry-run');
  const force = FLAGS.includes('--force');
  const targetDir = path.join(process.cwd(), '.github');

  if (!fs.existsSync(SOURCE_DIR)) {
    console.error('Error: source directory not found:', SOURCE_DIR);
    process.exit(1);
  }

  const fileCount = countFiles(SOURCE_DIR);

  if (fs.existsSync(targetDir) && !force) {
    console.error(`Error: ${targetDir} already exists.`);
    console.error('Use --force to overwrite.');
    process.exit(1);
  }

  if (dryRun) {
    console.log('[dry-run] Would copy:');
    console.log(`  ${SOURCE_DIR}`);
    console.log(`  -> ${targetDir}`);
    console.log(`  (${fileCount} files)`);
    return;
  }

  copyDirSync(SOURCE_DIR, targetDir);
  console.log(`✔ Installed .github/ (${fileCount} files) into ${targetDir}`);
}

function showHelp() {
  console.log(`
SATORI — Agent Skills for Science

Usage:
  satori init [--force] [--dry-run]   Install .github/ skills into current directory
  satori skill search <query>         Search skills by keyword
  satori skill info <name>            Show detailed skill information
  satori skill recommend <name>       Get similar/related skills
  satori pipeline suggest             Interactive pipeline recommendation
  satori pipeline list                List all available pipelines
  satori pipeline custom <action>     Manage custom pipelines
  satori docs generate [--preview]    Generate docs/ and docs/qiita/ files
  satori validate [--verbose]         Validate all SKILL.md files
  satori stats                        Show skill/TU coverage statistics
  satori help                         Show this help message
  satori --version, -v                Show version number

Options:
  --force     Overwrite existing .github/ directory
  --dry-run   Preview what would be installed without making changes
  --verbose   Show detailed validation output
  --preview   Show docs generation summary without writing files

Custom Pipelines:
  satori pipeline custom list         List custom pipelines
  satori pipeline custom add <path>   Add custom pipeline from file
  satori pipeline custom remove <id>  Remove custom pipeline
`);
}

// ── Pipeline Suggest ──

const PIPELINES = [
  {
    id: 1,
    name: '仮説検証→論文化',
    domain: 'general',
    keywords: ['仮説', '統計', '論文', 'hypothesis'],
    skills:
      'hypothesis-engine → data-preprocessing → statistical-testing → ml-classification → publication-figures → academic-writing → critical-review',
  },
  {
    id: 2,
    name: 'バリアント→臨床',
    domain: 'genomics',
    keywords: ['バリアント', 'variant', 'VCF', 'WGS', 'WES'],
    skills:
      'variant-interpretation → pharmacogenomics → precision-oncology → clinical-decision-support → clinical-reporting',
  },
  {
    id: 3,
    name: 'トランスクリプトーム',
    domain: 'genomics',
    keywords: ['RNA-seq', 'トランスクリプトーム', 'DEG', '発現'],
    skills: 'rnaseq-analysis → pathway-enrichment → network-analysis → publication-figures',
  },
  {
    id: 4,
    name: 'エピジェネティクス',
    domain: 'genomics',
    keywords: ['エピゲノム', 'ChIP-seq', 'ATAC-seq', 'メチル化'],
    skills: 'epigenomics-chromatin → regulatory-genomics → noncoding-rna → gene-regulation',
  },
  {
    id: 5,
    name: 'AlphaFold 構造解析',
    domain: 'structural',
    keywords: ['AlphaFold', 'タンパク質構造', '3D', 'protein structure'],
    skills: 'alphafold-structures → protein-structure-analysis → molecular-docking',
  },
  {
    id: 6,
    name: 'エビデンス合成',
    domain: 'literature',
    keywords: ['メタアナリシス', 'systematic review', '文献', 'エビデンス'],
    skills:
      'deep-research → literature-search → meta-analysis → evidence-synthesis → academic-writing → critical-review',
  },
  {
    id: 7,
    name: '創薬パイプライン',
    domain: 'pharma',
    keywords: ['創薬', 'drug discovery', 'ADMET', 'ドッキング'],
    skills:
      'drug-target-profiling → compound-screening → molecular-docking → admet-pharmacokinetics → drug-repurposing',
  },
  {
    id: 8,
    name: 'ML/XAI パイプライン',
    domain: 'ml',
    keywords: ['機械学習', 'ML', 'SHAP', 'XAI', '予測モデル'],
    skills:
      'data-preprocessing → ml-classification → ml-regression → explainable-ai → fairness-bias → publication-figures',
  },
  {
    id: 9,
    name: '環境・生態学',
    domain: 'ecology',
    keywords: ['生態', '生物多様性', 'SDM', '環境', 'ecology'],
    skills: 'environmental-ecology → biodiversity-conservation → species-distribution → time-series-forecasting',
  },
  {
    id: 10,
    name: '計算材料科学',
    domain: 'materials',
    keywords: ['材料', 'materials', 'DFT', '物性'],
    skills: 'computational-materials → cheminformatics → molecular-dynamics → ml-regression',
  },
  {
    id: 11,
    name: '医薬品安全性',
    domain: 'pharma',
    keywords: ['有害事象', 'ファーマコビジランス', '安全性', 'adverse'],
    skills: 'pharmacovigilance → pharmacogenomics → regulatory-science',
  },
  {
    id: 12,
    name: '希少疾患',
    domain: 'clinical',
    keywords: ['希少疾患', 'rare disease', 'Orphanet'],
    skills: 'rare-disease-genetics → gene-panel-design → variant-interpretation → clinical-reporting',
  },
  {
    id: 13,
    name: 'がんゲノミクス',
    domain: 'oncology',
    keywords: ['がん', 'cancer', 'TMB', '体細胞変異'],
    skills: 'cancer-genomics → precision-oncology → biomarker-discovery → clinical-reporting',
  },
  {
    id: 14,
    name: 'GWAS・集団遺伝学',
    domain: 'genomics',
    keywords: ['GWAS', '集団遺伝学', 'population genetics', 'biobank'],
    skills: 'biobank-cohort → population-genetics → statistical-testing → publication-figures',
  },
  {
    id: 15,
    name: 'シングルセル',
    domain: 'genomics',
    keywords: ['シングルセル', 'single-cell', 'scRNA-seq', '空間トランスクリプトーム'],
    skills: 'cellxgene-census → scvi-integration → spatial-transcriptomics → gene-regulation',
  },
  {
    id: 16,
    name: 'プロテオミクス',
    domain: 'omics',
    keywords: ['プロテオミクス', 'proteomics', '質量分析'],
    skills: 'proteomics → protein-structure-analysis → network-analysis',
  },
  {
    id: 17,
    name: 'メタボロミクス',
    domain: 'omics',
    keywords: ['メタボロミクス', 'metabolomics', '代謝物', '脂質'],
    skills: 'metabolomics → lipidomics → systems-biology → network-analysis',
  },
  {
    id: 18,
    name: 'マイクロバイオーム',
    domain: 'ecology',
    keywords: ['マイクロバイオーム', 'metagenome', '16S', '腸内細菌'],
    skills: 'microbiome-metagenomics → metagenome-assembled-genomes → phylogenetics → environmental-ecology',
  },
  {
    id: 19,
    name: 'パスウェイ・KG',
    domain: 'systems',
    keywords: ['パスウェイ', 'ナレッジグラフ', 'knowledge graph', 'pathway'],
    skills: 'gene-id-mapping → pathway-enrichment → ontology-integration → network-analysis → knowledge-graph',
  },
  {
    id: 20,
    name: '農業・食品',
    domain: 'agriculture',
    keywords: ['農業', '食品', 'agriculture', 'food safety'],
    skills: 'agricultural-science → food-science-nutrition → environmental-ecology',
  },
  {
    id: 21,
    name: '臨床情報学',
    domain: 'clinical',
    keywords: ['臨床', 'EHR', 'FHIR', 'OMOP', '電子カルテ'],
    skills: 'clinical-standards → clinical-nlp → clinical-reporting → healthcare-ai → survival-clinical',
  },
  {
    id: 22,
    name: 'ロボティクス・IoT',
    domain: 'engineering',
    keywords: ['ロボティクス', 'IoT', 'ロボット', 'robotics'],
    skills: 'robotics-automation → lab-automation → lab-data-management → interactive-dashboard',
  },
  {
    id: 23,
    name: '実験計画・統計',
    domain: 'general',
    keywords: ['実験計画', 'DOE', '検出力', 'サンプルサイズ'],
    skills:
      'experimental-design → experiment-fork → statistical-testing → experiment-template → reproducibility-assessment → publication-figures',
  },
  {
    id: 24,
    name: '科学的可視化',
    domain: 'general',
    keywords: ['可視化', 'visualization', 'ダッシュボード', 'dashboard'],
    skills: 'publication-figures → interactive-dashboard',
  },
  {
    id: 25,
    name: '学術出版',
    domain: 'literature',
    keywords: ['論文投稿', 'journal', 'グラント', 'grant'],
    skills: 'academic-writing → peer-review → latex-export → critical-review → citation-network',
  },
  {
    id: 26,
    name: '科学教育',
    domain: 'education',
    keywords: ['教育', 'education', 'カリキュラム'],
    skills: 'science-education → reproducibility-assessment',
  },
  // ── クロスドメインパイプライン ──
  {
    id: 'A',
    name: 'ゲノム創薬統合',
    domain: 'cross-domain',
    keywords: ['ゲノム創薬', 'GWAS', '創薬ターゲット', 'drug target', 'biobank'],
    skills:
      'biobank-cohort → population-genetics → drug-target-profiling → compound-screening → molecular-docking → admet-pharmacokinetics',
  },
  {
    id: 'B',
    name: 'AI 駆動臨床意思決定',
    domain: 'cross-domain',
    keywords: ['臨床AI', '予後予測', 'SHAP', '患者', 'clinical AI'],
    skills: 'clinical-decision-support → healthcare-ai → explainable-ai → pharmacovigilance → regulatory-science',
  },
  {
    id: 'C',
    name: '研究自動化',
    domain: 'cross-domain',
    keywords: ['研究自動化', '論文化', '仮説', 'research automation'],
    skills:
      'deep-research → hypothesis-pipeline → experiment-template → pipeline-scaffold → data-preprocessing → statistical-testing → publication-figures → academic-writing → systematic-review',
  },
  {
    id: 'D',
    name: 'マルチオミクス疾患解明',
    domain: 'cross-domain',
    keywords: ['マルチオミクス', '疾患', 'scRNA-seq', 'GRN', 'multi-omics'],
    skills:
      'single-cell-genomics → spatial-transcriptomics → disease-research → systems-biology → multi-omics → network-analysis',
  },
  {
    id: 'E',
    name: '個別化薬物療法',
    domain: 'cross-domain',
    keywords: ['個別化医療', 'PGx', 'Star アレル', '投与量最適化', 'pharmacogenomics'],
    skills:
      'variant-interpretation → pharmacogenomics → drug-target-profiling → admet-pharmacokinetics → clinical-decision-support → pharmacovigilance',
  },
  {
    id: 'F',
    name: 'バイオインフォマティクス完全',
    domain: 'cross-domain',
    keywords: ['バイオインフォマティクス', 'FASTQ', '配列解析', '統合パイプライン'],
    skills:
      'bioinformatics → single-cell-genomics → biobank-cohort → multi-omics → population-genetics → systems-biology → hypothesis-pipeline → academic-writing',
  },
  {
    id: 'G',
    name: 'がん精密医療 End-to-End',
    domain: 'cross-domain',
    keywords: ['がん精密医療', 'GDC', 'DepMap', '精密腫瘍学', 'TCGA'],
    skills:
      'gdc-portal → cancer-genomics → depmap-dependencies → civic-evidence → pharos-targets → compound-screening → precision-oncology → clinical-decision-support → healthcare-ai → survival-clinical',
  },
  {
    id: 'H',
    name: 'マルチオミクス縦断統合',
    domain: 'cross-domain',
    keywords: ['縦断統合', 'エピゲノム', 'プロテオーム', 'パスウェイ', 'VEP'],
    skills:
      'genome-sequence-tools → bioinformatics → variant-effect-prediction → epigenomics-chromatin → regulatory-genomics → cellxgene-census → scvi-integration → uniprot-proteome → alphafold-structures → protein-interaction-network → pathway-enrichment → reactome-pathways → network-visualization',
  },
  {
    id: 'I',
    name: '環境メタボ・マイクロバイオーム One Health',
    domain: 'cross-domain',
    keywords: ['One Health', '環境メタボ', '土壌', '微生物群集', 'SDM'],
    skills:
      'environmental-ecology → environmental-geodata → geospatial-analysis → microbiome-metagenomics → metagenome-assembled-genomes → phylogenetics → metabolomics-databases → metabolomics-network → metabolic-modeling → toxicology-env → publication-figures',
  },
  {
    id: 'J',
    name: 'AI 駆動マテリアルズインフォマティクス',
    domain: 'cross-domain',
    keywords: ['マテリアルズインフォマティクス', 'GNN', '能動学習', 'Materials Project', '材料探索'],
    skills:
      'computational-materials → cheminformatics → automl → graph-neural-networks → uncertainty-quantification → active-learning → doe → bayesian-statistics → adaptive-experiments → materials-characterization → advanced-visualization',
  },
  {
    id: 'K',
    name: '研究ライフサイクル完全自動化',
    domain: 'cross-domain',
    keywords: ['研究ライフサイクル', 'ラボ自動化', 'LIMS', 'ダッシュボード', 'グラント'],
    skills:
      'lab-automation → lab-data-management → streaming-analytics → model-monitoring → data-profiling → advanced-visualization → interactive-dashboard → scientific-schematics → reproducible-reporting → audit-report → paper-quality → peer-review → latex-formatter → latex-export → peer-review-response → grant-writing → preprint-archive',
  },
  {
    id: 'L',
    name: 'AI 駆動エビデンス合成',
    domain: 'cross-domain',
    keywords: ['エビデンス合成AI', 'DL文献', 'AutoML', 'スクリーニング'],
    skills:
      'deep-research → literature-search → text-mining-nlp → deep-learning → transfer-learning → automl → meta-analysis → explainable-ai → systematic-review → academic-writing',
  },
  {
    id: 'M',
    name: 'がんマルチレイヤーゲノム創薬',
    domain: 'cross-domain',
    keywords: ['がんゲノム創薬', 'ICGC', 'ChEMBL', 'エピゲノム'],
    skills:
      'gdc-portal → cancer-genomics → icgc-cancer-data → ensembl-genomics → variant-effect-prediction → epigenomics-chromatin → gwas-catalog → pharos-targets → chembl-assay-mining → compound-screening',
  },
  {
    id: 'N',
    name: '臨床→規制→出版バリューチェーン',
    domain: 'cross-domain',
    keywords: ['バリューチェーン', 'EHR', '規制報告', '学術出版', 'HL7'],
    skills:
      'clinical-standards → clinical-nlp → clinical-reporting → healthcare-ai → pharmacovigilance → regulatory-science → audit-report → reproducible-reporting → paper-quality → latex-formatter → peer-review-response',
  },
  {
    id: 'O',
    name: 'シングルセルプロテオーム統合',
    domain: 'cross-domain',
    keywords: ['シングルセルプロテオーム', '質量分析', '代謝モデル', 'MOFA+'],
    skills:
      'single-cell-genomics → spatial-transcriptomics → proteomics-mass-spectrometry → structural-proteomics → alphafold-structures → metabolomics-databases → metabolic-modeling → systems-biology → multi-omics',
  },
  // ── インダストリーパイプライン ──
  {
    id: 'Ind-1',
    name: '製薬企業レギュラトリー',
    domain: 'industry',
    keywords: ['製薬', 'CTD', 'レギュラトリー', '規制申請', 'regulatory'],
    skills:
      'drug-target-profiling → molecular-docking → admet-pharmacokinetics → clinical-trials-analytics → pharmacovigilance → regulatory-science → reproducible-reporting → paper-quality',
  },
  {
    id: 'Ind-2',
    name: '農業バイオテクノロジー',
    domain: 'industry',
    keywords: ['農業バイオ', '土壌微生物', 'CRISPR', '圃場', 'ゲノム編集'],
    skills:
      'environmental-ecology → microbiome-metagenomics → geospatial-analysis → plant-biology → crispr-design → gene-expression-transcriptomics → doe → publication-figures',
  },
  {
    id: 'Ind-3',
    name: '臨床検査室ワークフロー',
    domain: 'industry',
    keywords: ['臨床検査', 'NGS', 'ACMG', 'PGx', '臨床レポート'],
    skills:
      'genome-sequence-tools → variant-interpretation → pharmacogenomics → clinical-decision-support → clinical-standards → clinical-nlp → clinical-reporting',
  },
  {
    id: 'Ind-4',
    name: '食品安全・毒性評価',
    domain: 'industry',
    keywords: ['食品安全', '毒性', '残留農薬', 'フードセーフティ', 'food safety'],
    skills:
      'microbiome-metagenomics → rrna-taxonomy → metabolomics-databases → metabolomics-network → toxicology-env → data-profiling → regulatory-science → publication-figures',
  },
  {
    id: 'Ind-5',
    name: '法医・公衆衛生',
    domain: 'industry',
    keywords: ['法医学', '公衆衛生', 'アウトブレイク', 'サーベイランス', 'forensic'],
    skills:
      'variant-interpretation → population-genetics → infectious-disease → phylogenetics → immunoinformatics → epidemiology-public-health → public-health-data → biobank-cohort',
  },
  // ── メソドロジーパイプライン ──
  {
    id: 'M-α',
    name: 'ベイズ推論ワークフロー',
    domain: 'methodology',
    keywords: ['ベイズ', 'MCMC', '事後分布', 'Bayesian', '事前分布'],
    skills:
      'data-preprocessing → bayesian-statistics → statistical-simulation → uncertainty-quantification → doe → adaptive-experiments',
  },
  {
    id: 'M-β',
    name: '因果推論パイプライン',
    domain: 'methodology',
    keywords: ['因果推論', 'DAG', '傾向スコア', 'CATE', 'causal'],
    skills:
      'data-preprocessing → missing-data-analysis → causal-inference → causal-ml → explainable-ai → statistical-testing → publication-figures',
  },
  {
    id: 'M-γ',
    name: '時系列予測パイプライン',
    domain: 'methodology',
    keywords: ['時系列', 'Prophet', 'ARIMA', 'LSTM', '異常検知', 'forecasting'],
    skills:
      'data-preprocessing → time-series → time-series-forecasting → anomaly-detection → streaming-analytics → model-monitoring',
  },
  {
    id: 'M-δ',
    name: 'テキストマイニング・NLP',
    domain: 'methodology',
    keywords: ['テキストマイニング', 'NLP', 'PubTator', '引用ネットワーク', 'NER'],
    skills:
      'deep-research → literature-search → text-mining-nlp → biomedical-pubtator → clinical-nlp → semantic-scholar → citation-checker',
  },
];

// ── Synonym Dictionary ──
const SYNONYM_DICT = {
  // 機械学習・AI
  ml: ['machine learning', '機械学習', 'ML'],
  ai: ['artificial intelligence', '人工知能', 'AI', 'AI'],
  dl: ['deep learning', '深層学習', 'DL'],
  'neural network': ['NN', 'ニューラルネットワーク'],

  // バイオインフォマティクス
  bioinfo: ['バイオインフォマティクス', 'bioinformatics'],
  genomics: ['ゲノミクス', 'ゲノム', 'genomics'],
  seq: ['シーケンシング', 'sequencing'],
  rna: ['RNA', 'RNA-seq', 'トランスクリプトーム'],
  protein: ['プロテイン', 'タンパク質'],

  // 創薬・化学
  'drug discovery': ['創薬', '創薬', 'drug-discovery'],
  admet: ['ADMET', '薬物動態'],
  docking: ['ドッキング', 'molecular docking'],
  cheminformatics: ['ケモインフォマティクス'],

  // データ分析
  'data analysis': ['データ解析', 'data analysis'],
  statistics: ['統計', '統計学'],
  visualization: ['可視化', 'ビジュアル'],
  pipeline: ['パイプライン'],

  // 医療・臨床
  clinical: ['臨床', 'クリニカル'],
  precision: ['精密医療', '精密'],
  oncology: ['腫瘍学', 'がん'],
};

function normalizeKeyword(keyword) {
  const lower = keyword.toLowerCase().trim();

  // 同義語チェック
  for (const [key, synonyms] of Object.entries(SYNONYM_DICT)) {
    for (const syn of synonyms) {
      if (lower.includes(syn.toLowerCase()) || syn.toLowerCase().includes(lower)) {
        return key;
      }
    }
  }

  return lower;
}

function pipelineSuggest() {
  const readline = require('node:readline');
  const rl = readline.createInterface({ input: process.stdin, output: process.stdout });

  const ask = (q) => new Promise((resolve) => rl.question(q, resolve));

  (async () => {
    console.log('\n🔬 SATORI Pipeline Suggest — インタラクティブパイプライン推薦\n');
    console.log('研究内容を入力すると、最適なパイプラインを提案します。');
    console.log('(Ctrl+C で終了)\n');

    const input = await ask('何を解析しますか？ キーワードや研究テーマを入力してください:\n> ');
    const query = input.toLowerCase();
    const normalizedQuery = normalizeKeyword(query);

    // Score each pipeline by keyword match (with synonym support)
    const scored = PIPELINES.map((p) => {
      let score = 0;
      for (const kw of p.keywords) {
        const normalizedKw = normalizeKeyword(kw);
        // 完全一致：2点、同義語マッチ：1.5点、部分一致：1点
        if (query.includes(kw.toLowerCase())) score += 2;
        else if (normalizedQuery === normalizedKw) score += 1.5;
        else if (normalizedKw.includes(normalizedQuery) || query.includes(normalizedKw)) score += 1;
      }
      // Partial match on name
      if (query.includes(p.name.toLowerCase()) || p.name.toLowerCase().includes(query)) score += 0.5;
      return { ...p, score };
    })
      .filter((p) => p.score > 0)
      .sort((a, b) => b.score - a.score);

    console.log('');
    if (scored.length === 0) {
      console.log('❌ 該当するパイプラインが見つかりませんでした。');
      console.log('');
      console.log('利用可能なキーワード例:');
      console.log('  遺伝子/バリアント, 創薬/ADMET, RNA-seq, がん, 機械学習/ML,');
      console.log('  メタボロミクス, マイクロバイオーム, 環境/生態, 材料, 臨床/EHR,');
      console.log('  文献/メタアナリシス, 可視化, 論文, AlphaFold, シングルセル');
      console.log('');
      console.log('全パイプライン一覧は `satori pipeline list` で確認できます。');
    } else {
      console.log(`✅ ${scored.length} 件のパイプラインが見つかりました:\n`);
      const top = scored.slice(0, 5);
      for (const p of top) {
        console.log(`  📋 Pipeline #${p.id}: ${p.name}`);
        console.log(`     スキル連鎖: ${p.skills}`);
        console.log('');
      }
      if (scored.length > 5) {
        console.log(`  ... 他 ${scored.length - 5} 件`);
      }
      console.log('詳細は docs/SATORI_PIPELINE_EXAMPLES.md を参照してください。');
      console.log('全パイプライン一覧は `satori pipeline list` で確認できます。');
    }

    rl.close();
  })();
}

function pipelineList() {
  const domain = PIPELINES.filter((p) => typeof p.id === 'number');
  const cross = PIPELINES.filter((p) => p.domain === 'cross-domain');
  const industry = PIPELINES.filter((p) => p.domain === 'industry');
  const methodology = PIPELINES.filter((p) => p.domain === 'methodology');

  console.log(`\n📋 SATORI パイプライン一覧 (全 ${PIPELINES.length} パイプライン)\n`);

  console.log('── ドメインパイプライン (26) ──\n');
  for (const p of domain) {
    console.log(`  #${String(p.id).padStart(2, ' ')}  ${p.name}`);
    console.log(`       ${p.skills}`);
    console.log('');
  }

  console.log('── クロスドメインパイプライン (15) ──\n');
  for (const p of cross) {
    console.log(`  #${p.id}   ${p.name}`);
    console.log(`       ${p.skills}`);
    console.log('');
  }

  console.log('── インダストリーパイプライン (5) ──\n');
  for (const p of industry) {
    console.log(`  #${p.id}  ${p.name}`);
    console.log(`       ${p.skills}`);
    console.log('');
  }

  console.log('── メソドロジーパイプライン (4) ──\n');
  for (const p of methodology) {
    console.log(`  #${p.id}  ${p.name}`);
    console.log(`       ${p.skills}`);
    console.log('');
  }

  console.log('詳細は docs/SATORI_PIPELINE_EXAMPLES.md を参照してください。');
}

function showVersion() {
  const pkg = require(path.join(PACKAGE_ROOT, 'package.json'));
  console.log(pkg.version);
}

// ── Validate ──

function parseFrontmatter(content) {
  const match = content.match(/^---\n([\s\S]*?)\n---/);
  if (!match) return null;
  const yaml = match[1];
  const name = yaml.match(/^name:\s*(.+)$/m)?.[1]?.trim();
  const hasDescription = /^description:/m.test(yaml);
  return { name, hasDescription };
}

function validate() {
  const verbose = FLAGS.includes('--verbose');
  const skillsDir = path.join(SOURCE_DIR, 'skills');

  if (!fs.existsSync(skillsDir)) {
    console.error('Error: skills directory not found:', skillsDir);
    process.exit(1);
  }

  const dirs = fs
    .readdirSync(skillsDir)
    .filter((d) => d.startsWith('scientific-'))
    .sort();
  let pass = 0;
  let fail = 0;
  const errors = [];

  for (const dir of dirs) {
    const filePath = path.join(skillsDir, dir, 'SKILL.md');
    const issues = [];

    if (!fs.existsSync(filePath)) {
      issues.push('SKILL.md が見つかりません');
    } else {
      const content = fs.readFileSync(filePath, 'utf-8');
      const fm = parseFrontmatter(content);

      if (!fm) issues.push('YAML Frontmatter がありません');
      else {
        if (!fm.name) issues.push('Frontmatter に name がありません');
        else if (fm.name !== dir) issues.push(`name 不一致: "${fm.name}" (期待値: "${dir}")`);
        if (!fm.hasDescription) issues.push('Frontmatter に description がありません');
      }

      if (!/^# .+$/m.test(content)) issues.push('H1 タイトルがありません');
      if (!/^## When to Use/m.test(content)) issues.push('## When to Use セクションがありません');
      if (!/^## Quick Start/m.test(content)) issues.push('## Quick Start セクションがありません');
      if (!/```(?:python|markdown|json)/.test(content)) issues.push('コードブロックがありません');
    }

    if (issues.length === 0) {
      pass++;
      if (verbose) console.log(`  ✔ ${dir}`);
    } else {
      fail++;
      errors.push({ dir, issues });
      if (verbose) {
        console.log(`  ✘ ${dir}`);
        for (const issue of issues) console.log(`      - ${issue}`);
      }
    }
  }

  console.log(`\n📋 SKILL.md 検証結果: ${pass} pass / ${fail} fail (全 ${dirs.length} スキル)`);

  if (errors.length > 0 && !verbose) {
    console.log('\n問題のあるスキル:');
    for (const e of errors) {
      console.log(`  ✘ ${e.dir}: ${e.issues.join(', ')}`);
    }
  }

  if (fail > 0) {
    console.log('\n詳細は --verbose オプションで確認してください。');
    process.exit(1);
  } else {
    console.log('\n✔ 全スキルの検証に成功しました。');
  }
}

// ── Stats ──

function collectStats() {
  const skillsDir = path.join(SOURCE_DIR, 'skills');

  if (!fs.existsSync(skillsDir)) {
    console.error('Error: skills directory not found:', skillsDir);
    process.exit(1);
  }

  const dirs = fs
    .readdirSync(skillsDir)
    .filter((d) => d.startsWith('scientific-'))
    .sort();
  const totalSkills = dirs.length;
  let tuLinked = 0;
  let totalCodeBlocks = 0;
  const tuPattern = /ToolUniverse|利用可能ツール|SMCP/i;
  const tuKeyPattern = /`([A-Z][a-zA-Z]*_[a-z]+_[a-z_]+)`/g;
  const allTuKeys = new Set();

  for (const dir of dirs) {
    const filePath = path.join(skillsDir, dir, 'SKILL.md');
    if (!fs.existsSync(filePath)) continue;
    const content = fs.readFileSync(filePath, 'utf-8');

    if (tuPattern.test(content)) tuLinked++;

    const codeBlocks = content.match(/```(?:python|markdown|json)/g);
    if (codeBlocks) totalCodeBlocks += codeBlocks.length;

    for (const m of content.matchAll(tuKeyPattern)) {
      allTuKeys.add(m[1]);
    }
  }

  const pipelineBreakdown = {
    domain: PIPELINES.filter((p) => typeof p.id === 'number').length,
    cross: PIPELINES.filter((p) => p.domain === 'cross-domain').length,
    industry: PIPELINES.filter((p) => p.domain === 'industry').length,
    methodology: PIPELINES.filter((p) => p.domain === 'methodology').length,
  };

  const pkg = require(path.join(PACKAGE_ROOT, 'package.json'));
  const date = new Date().toISOString().slice(0, 10);
  const coverage = ((tuLinked / totalSkills) * 100).toFixed(1);

  return {
    version: pkg.version,
    date,
    totalSkills,
    pipelinesCount: PIPELINES.length,
    pipelineBreakdown,
    tuLinked,
    tuKeysCount: allTuKeys.size,
    totalCodeBlocks,
    coverage,
  };
}

function stats() {
  const summary = collectStats();

  console.log(`
📊 SATORI v${summary.version} — 統計

  スキル総数:          ${summary.totalSkills}
  パイプライン数:      ${summary.pipelinesCount}
  TU 連携スキル:       ${summary.tuLinked} (${summary.coverage}%)
  TU 未連携:           ${summary.totalSkills - summary.tuLinked}
  ユニーク TU キー:    ${summary.tuKeysCount}
  コードブロック総数:  ${summary.totalCodeBlocks}
`);
}

// ── Docs Generate ──

function updateReverseIndexDoc(content, summary) {
  const pipelineLine = `| パイプライン数 | **${summary.pipelinesCount}** (ドメイン ${summary.pipelineBreakdown.domain} + クロスドメイン ${summary.pipelineBreakdown.cross} + インダストリー ${summary.pipelineBreakdown.industry} + メソドロジー ${summary.pipelineBreakdown.methodology}) |`;

  return content
    .replace(/\n\| SATORI バージョン \| \*\*v[^*]+\*\* \|/g, `\n| SATORI バージョン | **v${summary.version}** |`)
    .replace(/\n\| 生成日 \| .* \|/g, `\n| 生成日 | ${summary.date} |`)
    .replace(/\n\| スキル数 \| \*\*\d+\*\* \|/g, `\n| スキル数 | **${summary.totalSkills}** |`)
    .replace(/\n\| パイプライン数 \| .* \|/g, `\n${pipelineLine}`)
    .replace(
      /\n\| ToolUniverse 連携スキル数 \| \*\*\d+\*\* \|/g,
      `\n| ToolUniverse 連携スキル数 | **${summary.tuLinked}** |`,
    )
    .replace(
      /\n\| ToolUniverse キー数（ユニーク） \| \*\*\d+\*\* \|/g,
      `\n| ToolUniverse キー数（ユニーク） | **${summary.tuKeysCount}** |`,
    );
}

function updatePipelineExamplesDoc(content, summary) {
  const pipelineLine = `| 掲載パイプライン数 | ${summary.pipelineBreakdown.domain} ドメイン + ${summary.pipelineBreakdown.cross} クロスドメイン + ${summary.pipelineBreakdown.industry} インダストリー + ${summary.pipelineBreakdown.methodology} メソドロジー = **${summary.pipelinesCount}** |`;

  return content
    .replace(
      /^> \*\*SATORI v[^*]+\*\* — .*$/m,
      `> **SATORI v${summary.version}** — ${summary.totalSkills} スキル + ${summary.pipelinesCount} パイプラインの連携レシピ集`,
    )
    .replace(/\n\| 生成日 \| .* \|/g, `\n| 生成日 | ${summary.date} |`)
    .replace(/\n\| 対象バージョン \| .* \|/g, `\n| 対象バージョン | v${summary.version} |`)
    .replace(/\n\| 掲載パイプライン数 \| .* \|/g, `\n${pipelineLine}`)
    .replace(
      /\n\| スキル総数 \| \d+ \(.*\) \|/g,
      `\n| スキル総数 | ${summary.totalSkills} (\`src/.github/skills/scientific-*/SKILL.md\`) |`,
    )
    .replace(/\n\| ToolUniverse キー数 \| .* \|/g, `\n| ToolUniverse キー数 | ${summary.tuKeysCount} (ユニーク) |`);
}

function updateQiitaReverseIndexDoc(content, summary) {
  const title = `title: 【SATORI v${summary.version}】${summary.totalSkills}スキル×${summary.pipelinesCount}パイプライン逆引き辞書 完全索引`;
  const intro = `**[SATORI](https://github.com/nahisaho/satori)** は **GitHub Copilot** 上で動作する、${summary.totalSkills} の専門スキルと ${summary.pipelinesCount} の統合パイプライン（${summary.pipelineBreakdown.domain} ドメイン + ${summary.pipelineBreakdown.cross} クロスドメイン + ${summary.pipelineBreakdown.industry} インダストリー + ${summary.pipelineBreakdown.methodology} メソドロジー）により、仮説構築から論文出版まで、あらゆる科学研究ワークフローを自動化するフレームワークです。`;

  return content.replace(/^title: .+$/m, title).replace(/^\*\*\[SATORI\][\s\S]*?フレームワークです。$/m, intro);
}

function updateQiitaPipelineExamplesDoc(content, summary) {
  const title = `title: 【SATORI v${summary.version}】${summary.totalSkills}スキル×${summary.pipelinesCount}パイプラインで実現する科学研究自動化 完全ガイド`;
  const intro = `**[SATORI](https://github.com/nahisaho/satori)** は **GitHub Copilot** 上で動作する、${summary.totalSkills} の専門スキルを組み合わせて構築した ${summary.pipelinesCount} 個のパイプライン（${summary.pipelineBreakdown.domain} ドメイン + ${summary.pipelineBreakdown.cross} クロスドメイン + ${summary.pipelineBreakdown.industry} インダストリー + ${summary.pipelineBreakdown.methodology} メソドロジー）により、仮説構築から論文出版まで、あらゆる科学研究ワークフローを自動化するフレームワークです。`;

  return content.replace(/^title: .+$/m, title).replace(/^\*\*\[SATORI\][\s\S]*?フレームワークです。$/m, intro);
}

function applyDocUpdate(filePath, updater, summary, preview) {
  if (!fs.existsSync(filePath)) {
    console.error(`Error: ドキュメントが見つかりません: ${filePath}`);
    return false;
  }
  const content = fs.readFileSync(filePath, 'utf-8');
  const updated = updater(content, summary);
  if (updated !== content && !preview) {
    fs.writeFileSync(filePath, updated);
  }
  return updated !== content;
}

function docsGenerate() {
  const summary = collectStats();
  const docsFlags = process.argv.slice(4);
  const preview = docsFlags.includes('--preview');

  const targets = [
    { path: path.join(PACKAGE_ROOT, 'docs', 'SATORI_REVERSE_INDEX.md'), update: updateReverseIndexDoc },
    { path: path.join(PACKAGE_ROOT, 'docs', 'SATORI_PIPELINE_EXAMPLES.md'), update: updatePipelineExamplesDoc },
    {
      path: path.join(PACKAGE_ROOT, 'docs', 'qiita', 'SATORI_REVERSE_INDEX_QIITA.md'),
      update: updateQiitaReverseIndexDoc,
    },
    {
      path: path.join(PACKAGE_ROOT, 'docs', 'qiita', 'SATORI_PIPELINE_EXAMPLES_QIITA.md'),
      update: updateQiitaPipelineExamplesDoc,
    },
  ];

  let updatedCount = 0;
  for (const target of targets) {
    if (applyDocUpdate(target.path, target.update, summary, preview)) {
      updatedCount++;
    }
  }

  if (preview) {
    console.log(`✔ docs generate (preview) 完了: ${updatedCount} 件更新予定`);
  } else {
    console.log(`✔ docs generate 完了: ${updatedCount} 件更新`);
  }
}

// ── Skill Search / Info ──

function loadAllSkills() {
  const skillsDir = path.join(SOURCE_DIR, 'skills');
  if (!fs.existsSync(skillsDir)) {
    console.error('Error: skills directory not found:', skillsDir);
    process.exit(1);
  }
  const dirs = fs
    .readdirSync(skillsDir)
    .filter((d) => d.startsWith('scientific-'))
    .sort();
  const skills = [];
  for (const dir of dirs) {
    const filePath = path.join(skillsDir, dir, 'SKILL.md');
    if (!fs.existsSync(filePath)) continue;
    const content = fs.readFileSync(filePath, 'utf-8');
    const fm = parseFrontmatter(content);
    const descMatch = content.match(/^description:\s*\|?\s*\n([\s\S]*?)(?=\n\w|\n---)/m);
    const description = descMatch ? descMatch[1].replace(/^\s+/gm, '').trim() : fm?.hasDescription ? '' : '';
    const tuKeyPattern = /`([A-Z][a-zA-Z]*_[a-z]+_[a-z_]+)`/g;
    const tuKeys = [];
    for (const m of content.matchAll(tuKeyPattern)) {
      tuKeys.push(m[1]);
    }
    const h1Match = content.match(/^# (.+)$/m);
    const title = h1Match ? h1Match[1].trim() : dir;
    skills.push({ dir, name: fm?.name || dir, title, description, content, tuKeys });
  }
  return skills;
}

function skillSearch() {
  const query = process.argv.slice(4).join(' ').toLowerCase();
  if (!query) {
    console.error('Error: 検索クエリを指定してください。');
    console.log('Usage: satori skill search <query>');
    process.exit(1);
  }

  const skills = loadAllSkills();
  const scored = skills
    .map((s) => {
      let score = 0;
      // 名前の完全一致
      if (s.name.toLowerCase() === query) score += 10;
      // 名前に含まれる
      else if (s.name.toLowerCase().includes(query)) score += 5;
      // タイトルに含まれる
      if (s.title.toLowerCase().includes(query)) score += 3;
      // 説明に含まれる
      if (s.description.toLowerCase().includes(query)) score += 2;
      // TU キーに含まれる
      for (const k of s.tuKeys) {
        if (k.toLowerCase().includes(query)) score += 1;
      }
      return { ...s, score };
    })
    .filter((s) => s.score > 0)
    .sort((a, b) => b.score - a.score);

  console.log(`\n🔍 "${process.argv.slice(4).join(' ')}" の検索結果\n`);
  if (scored.length === 0) {
    console.log('❌ 該当するスキルが見つかりませんでした。');
    console.log('');
    console.log('ヒント: 英語名（例: deep-learning, cancer-genomics）や');
    console.log('       日本語キーワード（例: 創薬, 機械学習）で検索してみてください。');
  } else {
    const top = scored.slice(0, 10);
    for (const s of top) {
      const desc = s.description ? s.description.split('\n')[0].substring(0, 60) : '';
      console.log(`  📖 ${s.name}`);
      if (desc) console.log(`     ${desc}`);
      console.log('');
    }
    if (scored.length > 10) {
      console.log(`  ... 他 ${scored.length - 10} 件`);
    }
    console.log(`合計 ${scored.length} 件がヒットしました。`);
    console.log('詳細は `satori skill info <name>` で確認できます。');
  }
}

function skillInfo() {
  const name = process.argv[4];
  if (!name) {
    console.error('Error: スキル名を指定してください。');
    console.log('Usage: satori skill info <name>');
    console.log('スキル検索は `satori skill search <query>` を使ってください。');
    process.exit(1);
  }

  const skillsDir = path.join(SOURCE_DIR, 'skills');
  // scientific- プレフィックスを自動補完
  const dirName = name.startsWith('scientific-') ? name : `scientific-${name}`;
  const filePath = path.join(skillsDir, dirName, 'SKILL.md');

  if (!fs.existsSync(filePath)) {
    console.error(`Error: スキル "${name}" が見つかりません。`);
    console.log('');
    // 部分一致候補を提示
    if (fs.existsSync(skillsDir)) {
      const dirs = fs
        .readdirSync(skillsDir)
        .filter((d) => d.startsWith('scientific-') && d.includes(name))
        .slice(0, 5);
      if (dirs.length > 0) {
        console.log('もしかして:');
        for (const d of dirs) {
          console.log(`  - ${d.replace('scientific-', '')}`);
        }
      }
    }
    process.exit(1);
  }

  const content = fs.readFileSync(filePath, 'utf-8');
  const fm = parseFrontmatter(content);
  const h1Match = content.match(/^# (.+)$/m);
  const title = h1Match ? h1Match[1].trim() : dirName;

  // 説明抽出
  const descMatch = content.match(/^description:\s*\|?\s*\n([\s\S]*?)(?=\n\w|\n---)/m);
  const description = descMatch ? descMatch[1].replace(/^\s+/gm, '').trim() : '';

  // When to Use セクション抽出
  const whenMatch = content.match(/^## When to Use\s*\n([\s\S]*?)(?=\n## )/m);
  const whenToUse = whenMatch ? whenMatch[1].trim() : '';

  // TU ツール
  const tuKeyPattern = /`([A-Z][a-zA-Z]*_[a-z]+_[a-z_]+)`/g;
  const tuKeys = new Set();
  for (const m of content.matchAll(tuKeyPattern)) {
    tuKeys.add(m[1]);
  }

  // tu_tools from frontmatter
  const tuToolMatches = content.match(/^tu_tools:\s*\n([\s\S]*?)(?=\n---|\n[a-z])/m);
  const tuToolNames = [];
  if (tuToolMatches) {
    const toolLines = tuToolMatches[1].matchAll(/name:\s*(.+)/g);
    for (const m of toolLines) {
      tuToolNames.push(m[1].trim());
    }
  }

  // 関連パイプライン
  const shortName = dirName.replace('scientific-', '');
  const relatedPipelines = PIPELINES.filter((p) => p.skills.includes(shortName));

  console.log(`\n📖 ${title}`);
  console.log(`   名前: ${fm?.name || dirName}`);
  if (description) {
    console.log(`   説明: ${description}`);
  }
  console.log('');

  if (whenToUse) {
    console.log('── When to Use ──');
    console.log(whenToUse);
    console.log('');
  }

  if (tuKeys.size > 0 || tuToolNames.length > 0) {
    console.log('── ToolUniverse 連携 ──');
    if (tuToolNames.length > 0) {
      for (const t of tuToolNames) {
        console.log(`  🔧 ${t}`);
      }
    }
    if (tuKeys.size > 0) {
      console.log(`  TU キー: ${[...tuKeys].join(', ')}`);
    }
    console.log('');
  }

  if (relatedPipelines.length > 0) {
    console.log('── 関連パイプライン ──');
    for (const p of relatedPipelines.slice(0, 5)) {
      console.log(`  📋 #${p.id}: ${p.name}`);
    }
    if (relatedPipelines.length > 5) {
      console.log(`  ... 他 ${relatedPipelines.length - 5} 件`);
    }
    console.log('');
  }

  console.log(`ファイル: src/.github/skills/${dirName}/SKILL.md`);
}

function skillRecommend() {
  const name = process.argv[4];
  if (!name) {
    console.error('Error: スキル名を指定してください。');
    console.log('Usage: satori skill recommend <name>');
    process.exit(1);
  }

  const skillsDir = path.join(SOURCE_DIR, 'skills');
  const dirName = name.startsWith('scientific-') ? name : `scientific-${name}`;
  const filePath = path.join(skillsDir, dirName, 'SKILL.md');

  if (!fs.existsSync(filePath)) {
    console.error(`Error: スキル "${name}" が見つかりません。`);
    process.exit(1);
  }

  // スキル名から短縮形を取得
  const shortName = dirName.replace('scientific-', '');

  // 全パイプラインから、このスキルを使用するパイプラインを検出
  const usedInPipelines = PIPELINES.filter((p) => p.skills.includes(shortName));

  // これらのパイプラインで使用される他のスキルをカウント
  const skillCooccurrence = {};
  for (const p of usedInPipelines) {
    const skillNames = p.skills.split(' → ').map((s) => s.trim());
    for (const sk of skillNames) {
      if (sk !== shortName && !skillCooccurrence[sk]) {
        skillCooccurrence[sk] = 0;
      }
      if (sk !== shortName) {
        skillCooccurrence[sk]++;
      }
    }
  }

  // スコア降順にソート
  const recommended = Object.entries(skillCooccurrence)
    .map(([skill, count]) => ({ skill, count }))
    .sort((a, b) => b.count - a.count)
    .slice(0, 5);

  console.log(`\n🎯 "${name}" に関連するスキル\n`);

  if (usedInPipelines.length > 0) {
    console.log(`このスキルが使用されるパイプライン: ${usedInPipelines.length} 件\n`);
    console.log('関連スキル:');
    for (let i = 0; i < recommended.length; i++) {
      const { skill, count } = recommended[i];
      console.log(`  ${i + 1}. ${skill} (${count} パイプラインで併用)`);
    }
    console.log('');
    console.log('詳細は `satori skill info <related-skill>` で確認できます。');
  } else {
    console.log('❌ このスキルが使用されるパイプラインが見つかりませんでした。');
    console.log('すべてのパイプラインは `satori pipeline list` で確認できます。');
  }
}

// ── Custom Pipeline Management ──

function loadCustomPipelines() {
  if (!fs.existsSync(CUSTOM_PIPELINES_PATH)) {
    return [];
  }
  try {
    const content = fs.readFileSync(CUSTOM_PIPELINES_PATH, 'utf-8');
    const data = JSON.parse(content);
    return data.customPipelines || [];
  } catch (err) {
    console.error('Warning: Failed to load custom pipelines:', err.message);
    return [];
  }
}

function saveCustomPipelines(pipelines) {
  const dir = path.dirname(CUSTOM_PIPELINES_PATH);
  fs.mkdirSync(dir, { recursive: true });
  fs.writeFileSync(CUSTOM_PIPELINES_PATH, JSON.stringify({ customPipelines: pipelines }, null, 2));
}

function pipelineCustom() {
  const action = process.argv[4];
  const customPipelines = loadCustomPipelines();

  if (action === 'list') {
    listCustomPipelines(customPipelines);
  } else if (action === 'add') {
    addCustomPipeline(customPipelines);
  } else if (action === 'remove') {
    removeCustomPipeline(customPipelines);
  } else {
    console.error(`Unknown custom pipeline action: ${action || '(none)'}`);
    console.log('Usage: satori pipeline custom list | add <file> | remove <id>');
    process.exit(1);
  }
}

function listCustomPipelines(pipelines) {
  if (pipelines.length === 0) {
    console.log('\n📋 カスタムパイプライン: 0 件\n');
    console.log('新しいカスタムパイプラインを追加するには:');
    console.log('  satori pipeline custom add <file>');
    return;
  }

  console.log(`\n📋 カスタムパイプライン一覧 (${pipelines.length} 件)\n`);
  for (const p of pipelines) {
    console.log(`  🔧 [${p.id}] ${p.name}`);
    console.log(`     スキル連鎖: ${p.skills}`);
    console.log('');
  }
}

function addCustomPipeline(pipelines) {
  const filePath = process.argv[5];
  if (!filePath) {
    console.error('Error: ファイルパスを指定してください。');
    console.log('Usage: satori pipeline custom add <file>');
    process.exit(1);
  }

  if (!fs.existsSync(filePath)) {
    console.error(`Error: ファイルが見つかりません: ${filePath}`);
    process.exit(1);
  }

  try {
    const content = fs.readFileSync(filePath, 'utf-8');
    const pipelineData = JSON.parse(content);

    // バリデーション
    if (!pipelineData.id || !pipelineData.name || !pipelineData.skills) {
      console.error('Error: パイプラインは id, name, skills を含む必要があります。');
      process.exit(1);
    }

    // 重複チェック
    if (pipelines.some((p) => p.id === pipelineData.id)) {
      console.error(`Error: ID "${pipelineData.id}" は既に存在します。`);
      process.exit(1);
    }

    pipelines.push(pipelineData);
    saveCustomPipelines(pipelines);
    console.log(`✔ カスタムパイプライン "${pipelineData.name}" を追加しました。`);
  } catch (err) {
    console.error('Error: パイプラインファイルの解析に失敗しました:', err.message);
    process.exit(1);
  }
}

function removeCustomPipeline(pipelines) {
  const id = process.argv[5];
  if (!id) {
    console.error('Error: パイプライン ID を指定してください。');
    console.log('Usage: satori pipeline custom remove <id>');
    process.exit(1);
  }

  const index = pipelines.findIndex((p) => p.id === id);
  if (index === -1) {
    console.error(`Error: パイプライン "${id}" が見つかりません。`);
    process.exit(1);
  }

  const removed = pipelines.splice(index, 1)[0];
  saveCustomPipelines(pipelines);
  console.log(`✔ カスタムパイプライン "${removed.name}" を削除しました。`);
}

switch (COMMAND) {
  case 'init':
    init();
    break;
  case 'skill':
    if (SUBCOMMAND === 'search') {
      skillSearch();
    } else if (SUBCOMMAND === 'info') {
      skillInfo();
    } else if (SUBCOMMAND === 'recommend') {
      skillRecommend();
    } else {
      console.error(`Unknown skill subcommand: ${SUBCOMMAND || '(none)'}`);
      console.log('Usage: satori skill search <query> | satori skill info <name> | satori skill recommend <name>');
      process.exit(1);
    }
    break;
  case 'pipeline':
    if (SUBCOMMAND === 'suggest') {
      pipelineSuggest();
    } else if (SUBCOMMAND === 'list') {
      pipelineList();
    } else if (SUBCOMMAND === 'custom') {
      pipelineCustom();
    } else {
      console.error(`Unknown pipeline subcommand: ${SUBCOMMAND || '(none)'}`);
      console.log('Usage: satori pipeline suggest | satori pipeline list | satori pipeline custom list|add|remove');
      process.exit(1);
    }
    break;
  case 'docs':
    if (SUBCOMMAND === 'generate') {
      docsGenerate();
    } else {
      console.error(`Unknown docs subcommand: ${SUBCOMMAND || '(none)'}`);
      console.log('Usage: satori docs generate [--preview]');
      process.exit(1);
    }
    break;
  case 'validate':
    validate();
    break;
  case 'stats':
    stats();
    break;
  case 'help':
  case '--help':
  case '-h':
  case undefined:
    showHelp();
    break;
  case '--version':
  case '-v':
  case 'version':
    showVersion();
    break;
  default:
    console.error(`Unknown command: ${COMMAND}`);
    showHelp();
    process.exit(1);
}
