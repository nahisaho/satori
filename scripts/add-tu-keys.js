#!/usr/bin/env node
/**
 * SATORI TU Key 一括追加スクリプト
 *
 * 87 の TU 未連携スキルに ToolUniverse 連携セクション + frontmatter tu_tools を追加する。
 * 各スキルのドメインに適した TU key をマッピング。
 */

const fs = require('node:fs');
const path = require('node:path');

const SKILLS_DIR = path.resolve(__dirname, '../src/.github/skills');

// スキル → TU key マッピング
// ToolUniverse の実際のデータベース/ツールに対応するキーを割り当て
const TU_MAP = {
  // ── [A] 基盤・ワークフロー ──
  'scientific-academic-writing': [
    { key: 'crossref', name: 'Crossref', desc: '論文メタデータ検索・引用情報取得' },
  ],
  'scientific-critical-review': [
    { key: 'crossref', name: 'Crossref', desc: '論文検証・引用メタデータ参照' },
  ],
  'scientific-data-preprocessing': [
    { key: 'biotools', name: 'bio.tools', desc: 'データ前処理ツールレジストリ検索' },
  ],
  'scientific-data-simulation': [
    { key: 'biotools', name: 'bio.tools', desc: 'シミュレーションツールレジストリ検索' },
  ],
  'scientific-data-submission': [
    { key: 'ena', name: 'ENA', desc: 'データ登録・アクセッション管理' },
  ],
  'scientific-hypothesis-pipeline': [
    { key: 'open_alex', name: 'OpenAlex', desc: '仮説関連文献の網羅的検索' },
  ],
  'scientific-latex-formatter': [
    { key: 'crossref', name: 'Crossref', desc: '参考文献メタデータ・DOI 解決' },
  ],
  'scientific-paper-quality': [
    { key: 'crossref', name: 'Crossref', desc: '引用品質・ジャーナルメトリクス参照' },
  ],
  'scientific-peer-review-response': [
    { key: 'crossref', name: 'Crossref', desc: '査読指摘の文献裏付け検索' },
  ],
  'scientific-pipeline-scaffold': [
    { key: 'biotools', name: 'bio.tools', desc: 'パイプライン構成ツール検索' },
  ],
  'scientific-publication-figures': [
    { key: 'biotools', name: 'bio.tools', desc: '可視化ツールレジストリ検索' },
  ],
  'scientific-research-methodology': [
    { key: 'open_alex', name: 'OpenAlex', desc: '研究方法論の文献ベース構築' },
  ],
  'scientific-revision-tracker': [
    { key: 'crossref', name: 'Crossref', desc: '改訂履歴の引用整合性検証' },
  ],
  'scientific-supplementary-generator': [
    { key: 'crossref', name: 'Crossref', desc: '補足資料の引用メタデータ参照' },
  ],

  // ── [B] 統計・探索的解析 ──
  'scientific-eda-correlation': [
    { key: 'biotools', name: 'bio.tools', desc: '統計解析ツール検索' },
  ],
  'scientific-statistical-testing': [
    { key: 'biotools', name: 'bio.tools', desc: '統計検定ツールレジストリ' },
  ],
  'scientific-pca-tsne': [
    { key: 'biotools', name: 'bio.tools', desc: '次元削減・可視化ツール検索' },
  ],
  'scientific-advanced-visualization': [
    { key: 'biotools', name: 'bio.tools', desc: '高度可視化ツール検索' },
  ],
  'scientific-data-profiling': [
    { key: 'biotools', name: 'bio.tools', desc: 'データプロファイリングツール検索' },
  ],
  'scientific-geospatial-analysis': [
    { key: 'gbif', name: 'GBIF', desc: '地理空間生物多様性データ取得' },
  ],
  'scientific-missing-data-analysis': [
    { key: 'biotools', name: 'bio.tools', desc: '欠損データ処理ツール検索' },
  ],
  'scientific-network-visualization': [
    { key: 'string', name: 'STRING', desc: 'ネットワーク可視化・相互作用データ' },
  ],
  'scientific-statistical-simulation': [
    { key: 'biotools', name: 'bio.tools', desc: '統計シミュレーションツール検索' },
  ],
  'scientific-streaming-analytics': [
    { key: 'biotools', name: 'bio.tools', desc: 'ストリーミング解析ツール検索' },
  ],
  'scientific-symbolic-mathematics': [
    { key: 'biotools', name: 'bio.tools', desc: '数式処理・解析ツール検索' },
  ],

  // ── [C] 機械学習・モデリング ──
  'scientific-ml-regression': [
    { key: 'openml', name: 'OpenML', desc: '回帰ベンチマーク・データセット取得' },
  ],
  'scientific-ml-classification': [
    { key: 'openml', name: 'OpenML', desc: '分類ベンチマーク・データセット取得' },
  ],
  'scientific-feature-importance': [
    { key: 'openml', name: 'OpenML', desc: '特徴量選択ベンチマーク参照' },
  ],
  'scientific-active-learning': [
    { key: 'openml', name: 'OpenML', desc: '能動学習データセット・評価指標' },
  ],
  'scientific-automl': [
    { key: 'openml', name: 'OpenML', desc: 'AutoML ベンチマーク・タスク参照' },
  ],
  'scientific-ensemble-methods': [
    { key: 'openml', name: 'OpenML', desc: 'アンサンブル手法ベンチマーク参照' },
  ],
  'scientific-anomaly-detection': [
    { key: 'openml', name: 'OpenML', desc: '異常検知ベンチマーク・データセット' },
  ],
  'scientific-causal-ml': [
    { key: 'openml', name: 'OpenML', desc: '因果推論 ML データセット参照' },
  ],
  'scientific-model-monitoring': [
    { key: 'openml', name: 'OpenML', desc: 'モデルモニタリング指標・ベンチマーク' },
  ],
  'scientific-semi-supervised-learning': [
    { key: 'openml', name: 'OpenML', desc: '半教師あり学習ベンチマーク' },
  ],
  'scientific-multi-task-learning': [
    { key: 'openml', name: 'OpenML', desc: 'マルチタスク学習データセット参照' },
  ],

  // ── [D] 実験計画・プロセス最適化 ──
  'scientific-doe': [
    { key: 'biotools', name: 'bio.tools', desc: '実験計画法ツールレジストリ検索' },
  ],
  'scientific-process-optimization': [
    { key: 'biotools', name: 'bio.tools', desc: 'プロセス最適化ツール検索' },
  ],
  'scientific-adaptive-experiments': [
    { key: 'biotools', name: 'bio.tools', desc: '適応的実験設計ツール検索' },
  ],

  // ── [E] 信号・スペクトル・時系列 ──
  'scientific-spectral-signal': [
    { key: 'biotools', name: 'bio.tools', desc: 'スペクトル信号処理ツール検索' },
  ],
  'scientific-biosignal-processing': [
    { key: 'biotools', name: 'bio.tools', desc: '生体信号処理ツール検索' },
  ],
  'scientific-time-series': [
    { key: 'biotools', name: 'bio.tools', desc: '時系列解析ツール検索' },
  ],
  'scientific-time-series-forecasting': [
    { key: 'biotools', name: 'bio.tools', desc: '時系列予測ツール検索' },
  ],
  'scientific-neuroscience-electrophysiology': [
    { key: 'biotools', name: 'bio.tools', desc: '電気生理学解析ツール検索' },
  ],

  // ── [F] 生命科学追加 ──
  'scientific-glycomics': [
    { key: 'glygen', name: 'GlyGen', desc: '糖鎖構造・機能データベース検索' },
  ],
  'scientific-lipidomics': [
    { key: 'lipidmaps', name: 'LIPID MAPS', desc: '脂質構造・分類データベース検索' },
  ],
  'scientific-hgnc-nomenclature': [
    { key: 'hgnc', name: 'HGNC', desc: '遺伝子命名法・シンボル検索' },
  ],
  'scientific-metabolomics-network': [
    { key: 'hmdb', name: 'HMDB', desc: '代謝物ネットワーク・代謝経路検索' },
  ],

  // ── [G] 化学・材料・イメージング ──
  'scientific-image-analysis': [
    { key: 'biotools', name: 'bio.tools', desc: '画像解析ツールレジストリ検索' },
  ],
  'scientific-materials-characterization': [
    { key: 'materials_project', name: 'Materials Project', desc: '材料特性データベース検索' },
  ],
  'scientific-advanced-imaging': [
    { key: 'biotools', name: 'bio.tools', desc: '高度イメージングツール検索' },
  ],
  'scientific-deep-chemistry': [
    { key: 'chembl', name: 'ChEMBL', desc: '化学的活性・化合物データ検索' },
  ],
  'scientific-md-simulation': [
    { key: 'pdb', name: 'PDB', desc: '分子構造データベース参照' },
  ],

  // ── [H] 臨床・疫学 ──
  'scientific-biobank-cohort': [
    { key: 'clinvar', name: 'ClinVar', desc: 'バリアント臨床的意義データ検索' },
  ],
  'scientific-causal-inference': [
    { key: 'open_alex', name: 'OpenAlex', desc: '因果推論関連文献検索' },
  ],
  'scientific-clinical-standards': [
    { key: 'loinc', name: 'LOINC', desc: '臨床検査コード標準検索' },
  ],
  'scientific-clinical-nlp': [
    { key: 'umls', name: 'UMLS', desc: '医学用語統一システム検索' },
  ],
  'scientific-clinical-pharmacology': [
    { key: 'drugbank', name: 'DrugBank', desc: '薬物動態・相互作用データ検索' },
  ],

  // ── [J] 創薬 ──
  'scientific-nci60-screening': [
    { key: 'nci60', name: 'NCI-60', desc: 'がん細胞株スクリーニングデータ検索' },
  ],

  // ── [M] 実験室自動化 ──
  'scientific-lab-automation': [
    { key: 'biotools', name: 'bio.tools', desc: '実験自動化ツールレジストリ検索' },
  ],
  'scientific-crispr-design': [
    { key: 'ensembl', name: 'Ensembl', desc: 'gRNA 設計用ゲノム配列参照' },
  ],

  // ── [N] 科学プレゼンテーション ──
  'scientific-presentation-design': [
    { key: 'biotools', name: 'bio.tools', desc: 'プレゼンテーション可視化ツール検索' },
  ],
  'scientific-interactive-dashboard': [
    { key: 'biotools', name: 'bio.tools', desc: 'インタラクティブ可視化ツール検索' },
  ],
  'scientific-reproducible-reporting': [
    { key: 'biotools', name: 'bio.tools', desc: '再現可能レポーティングツール検索' },
  ],

  // ── [R] 量子・先端計算 ──
  'scientific-quantum-computing': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: '量子計算論文・ベンチマーク検索' },
  ],
  'scientific-bayesian-statistics': [
    { key: 'biotools', name: 'bio.tools', desc: 'ベイズ統計ツールレジストリ検索' },
  ],
  'scientific-explainable-ai': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: 'XAI 手法・ベンチマーク検索' },
  ],
  'scientific-deep-learning': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: '深層学習モデル・ベンチマーク検索' },
  ],
  'scientific-reinforcement-learning': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: '強化学習環境・ベンチマーク検索' },
  ],
  'scientific-transfer-learning': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: '転移学習モデル・事前学習済みモデル検索' },
  ],
  'scientific-uncertainty-quantification': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: '不確実性定量化手法・ベンチマーク' },
  ],
  'scientific-federated-learning': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: '連合学習フレームワーク・ベンチマーク' },
  ],
  'scientific-neural-architecture-search': [
    { key: 'papers_with_code', name: 'Papers with Code', desc: 'NAS 探索空間・ベンチマーク検索' },
  ],

  // ── [S] 医用イメージング ──
  'scientific-medical-imaging': [
    { key: 'tcia', name: 'TCIA', desc: 'がん画像アーカイブ検索' },
  ],
  'scientific-radiology-ai': [
    { key: 'tcia', name: 'TCIA', desc: '放射線画像データセット検索' },
  ],

  // ── [T] シングルセル・空間追加 ──
  'scientific-gpu-singlecell': [
    { key: 'cellxgene', name: 'CellxGene', desc: 'シングルセルデータセット検索' },
  ],
  'scientific-perturbation-analysis': [
    { key: 'cellxgene', name: 'CellxGene', desc: '摂動解析データセット検索' },
  ],
  'scientific-scvi-integration': [
    { key: 'cellxgene', name: 'CellxGene', desc: 'scVI 統合用データセット検索' },
  ],
  'scientific-scatac-signac': [
    { key: 'encode', name: 'ENCODE', desc: 'scATAC-seq 参照エピゲノムデータ' },
  ],
  'scientific-squidpy-advanced': [
    { key: 'cellxgene', name: 'CellxGene', desc: '空間トランスクリプトミクスデータ検索' },
  ],
  'scientific-spatial-multiomics': [
    { key: 'cellxgene', name: 'CellxGene', desc: '空間マルチオミクスデータ検索' },
  ],

  // ── [V] マイクロバイオーム・環境 ──
  'scientific-metagenome-assembled-genomes': [
    { key: 'mgnify', name: 'MGnify', desc: 'メタゲノムアセンブル・MAG データ検索' },
  ],
  'scientific-phylogenetics': [
    { key: 'ncbi_taxonomy', name: 'NCBI Taxonomy', desc: '系統分類・分岐学データ検索' },
  ],
  'scientific-plant-biology': [
    { key: 'tair', name: 'TAIR', desc: 'シロイヌナズナゲノム・植物データ検索' },
  ],

  // ── [W] システム生物学 ──
  'scientific-metabolic-atlas': [
    { key: 'metabolic_atlas', name: 'Metabolic Atlas', desc: 'ヒトゲノム規模代謝モデル検索' },
  ],
  'scientific-metabolic-flux': [
    { key: 'bigg', name: 'BiGG Models', desc: '代謝フラックスモデル検索' },
  ],

  // ── [X] 疫学・公衆衛生 ──
  'scientific-toxicology-env': [
    { key: 'ctd', name: 'CTD', desc: '化学物質-疾患-遺伝子関連データ検索' },
  ],

  // ── [N/other] ──
  'scientific-scientific-schematics': undefined,  // already linked - skip
};

let updated = 0;
let skipped = 0;
const errors = [];

for (const [skill, tools] of Object.entries(TU_MAP)) {
  if (!tools) { skipped++; continue; }

  const filePath = path.join(SKILLS_DIR, skill, 'SKILL.md');
  if (!fs.existsSync(filePath)) {
    errors.push(`${skill}: SKILL.md not found`);
    continue;
  }

  let content = fs.readFileSync(filePath, 'utf-8');

  // Skip if already has TU
  if (/ToolUniverse|利用可能ツール|SMCP/i.test(content)) {
    skipped++;
    continue;
  }

  // 1. Add tu_tools to frontmatter
  const fmMatch = content.match(/^---\n([\s\S]*?)\n---/);
  if (!fmMatch) {
    errors.push(`${skill}: no frontmatter found`);
    continue;
  }

  const existingFm = fmMatch[1];
  let tuYaml = 'tu_tools:';
  for (const t of tools) {
    tuYaml += `\n  - key: ${t.key}\n    name: ${t.name}\n    description: ${t.desc}`;
  }

  // Update description to mention ToolUniverse
  const descMatch = existingFm.match(/^description:\s*([\s\S]*?)(?=\n\w|\n---)/m);
  let newFm = existingFm;

  // Add ToolUniverse mention to description
  const tuKeys = tools.map(t => t.key).join(', ');
  if (descMatch) {
    const descVal = descMatch[1].trim();
    // Check if description is multiline (|)
    if (descVal.startsWith('|')) {
      // Multiline - add TU mention at end
      const descLines = descMatch[0].split('\n');
      const lastLine = descLines[descLines.length - 1].trim();
      if (!lastLine.includes('ToolUniverse')) {
        descLines.push(`  ToolUniverse 連携: ${tuKeys}。`);
      }
      newFm = existingFm.replace(descMatch[0], descLines.join('\n'));
    } else {
      // Single line or short - convert to multiline
      const cleanDesc = descVal.replace(/\n\s*/g, ' ').trim();
      const newDesc = `description: |\n  ${cleanDesc}\n  ToolUniverse 連携: ${tuKeys}。`;
      newFm = existingFm.replace(descMatch[0], newDesc);
    }
  }

  // Add tu_tools block
  newFm = newFm.trimEnd() + '\n' + tuYaml;

  // Replace frontmatter
  content = content.replace(fmMatch[0], `---\n${newFm}\n---`);

  // 2. Add ToolUniverse section before last section or at end
  // Insert before ## References if exists, otherwise append
  const tuSection = `\n## ToolUniverse 連携\n\n| TU Key | ツール名 | 連携内容 |\n|--------|---------|--------|\n${tools.map(t => `| \`${t.key}\` | ${t.name} | ${t.desc} |`).join('\n')}\n`;

  const refIdx = content.lastIndexOf('\n## References');
  if (refIdx >= 0) {
    content = content.slice(0, refIdx) + tuSection + content.slice(refIdx);
  } else {
    content = content.trimEnd() + '\n' + tuSection;
  }

  fs.writeFileSync(filePath, content, 'utf-8');
  updated++;
}

console.log(`Updated: ${updated}`);
console.log(`Skipped: ${skipped}`);
if (errors.length) {
  console.log(`Errors: ${errors.length}`);
  errors.forEach(e => console.log(`  - ${e}`));
}
