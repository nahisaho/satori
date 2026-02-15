#!/usr/bin/env node

const fs = require('node:fs');
const path = require('node:path');

const COMMAND = process.argv[2];
const SUBCOMMAND = process.argv[3];
const FLAGS = process.argv.slice(3);

const PACKAGE_ROOT = path.resolve(__dirname, '..');
const SOURCE_DIR = path.join(PACKAGE_ROOT, 'src', '.github');

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
  satori pipeline suggest             Interactive pipeline recommendation
  satori pipeline list                List all available pipelines
  satori validate [--verbose]         Validate all SKILL.md files
  satori stats                        Show skill/TU coverage statistics
  satori help                         Show this help message
  satori --version, -v                Show version number

Options:
  --force     Overwrite existing .github/ directory
  --dry-run   Preview what would be installed without making changes
  --verbose   Show detailed validation output
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
    skills: 'experimental-design → statistical-testing → reproducibility-assessment → publication-figures',
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
    skills: 'academic-writing → critical-review → citation-network',
  },
  {
    id: 26,
    name: '科学教育',
    domain: 'education',
    keywords: ['教育', 'education', 'カリキュラム'],
    skills: 'science-education → reproducibility-assessment',
  },
];

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

    // Score each pipeline by keyword match
    const scored = PIPELINES.map((p) => {
      let score = 0;
      for (const kw of p.keywords) {
        if (query.includes(kw.toLowerCase())) score += 2;
      }
      // Partial match on name
      if (query.includes(p.name.toLowerCase()) || p.name.toLowerCase().includes(query)) score += 1;
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
    }

    rl.close();
  })();
}

function pipelineList() {
  console.log('\n📋 SATORI パイプライン一覧 (26 ドメインパイプライン)\n');
  for (const p of PIPELINES) {
    console.log(`  #${String(p.id).padStart(2, ' ')}  ${p.name}`);
    console.log(`       ${p.skills}`);
    console.log('');
  }
  console.log('クロスドメイン (15), 産業特化 (5), 方法論特化 (4) パイプラインは');
  console.log('docs/SATORI_PIPELINE_EXAMPLES.md を参照してください。');
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

function stats() {
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

  const coverage = ((tuLinked / totalSkills) * 100).toFixed(1);
  const pkg = require(path.join(PACKAGE_ROOT, 'package.json'));

  console.log(`
📊 SATORI v${pkg.version} — 統計

  スキル総数:          ${totalSkills}
  パイプライン数:      ${PIPELINES.length}
  TU 連携スキル:       ${tuLinked} (${coverage}%)
  TU 未連携:           ${totalSkills - tuLinked}
  ユニーク TU キー:    ${allTuKeys.size}
  コードブロック総数:  ${totalCodeBlocks}
`);
}

switch (COMMAND) {
  case 'init':
    init();
    break;
  case 'pipeline':
    if (SUBCOMMAND === 'suggest') {
      pipelineSuggest();
    } else if (SUBCOMMAND === 'list') {
      pipelineList();
    } else {
      console.error(`Unknown pipeline subcommand: ${SUBCOMMAND || '(none)'}`);
      console.log('Usage: satori pipeline suggest | satori pipeline list');
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
