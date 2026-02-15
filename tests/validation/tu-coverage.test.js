/**
 * TU (ToolUniverse) カバレッジ分析テスト
 *
 * 全 190 スキルの TU 連携状況を解析し、
 * カバレッジ統計とギャップを出力する。
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';

const SKILLS_DIR = path.resolve(__dirname, '../../src/.github/skills');

function getSkillDirs() {
  return fs
    .readdirSync(SKILLS_DIR)
    .filter((d) => d.startsWith('scientific-'))
    .sort();
}

function readSkillMd(skillName) {
  const filePath = path.join(SKILLS_DIR, skillName, 'SKILL.md');
  if (!fs.existsSync(filePath)) return null;
  return fs.readFileSync(filePath, 'utf-8');
}

/** TU キーパターン検出: 大文字始まり + アンダースコア + 動詞 のパターン */
const TU_KEY_PATTERN = /`?([A-Z][a-zA-Z]*(?:_[a-zA-Z]+)+)\s*\(/g;

/** ToolUniverse セクション見出し */
const TU_SECTION_PATTERN = /ToolUniverse|利用可能ツール|SMCP/i;

/** TU キー抽出 */
function extractTuKeys(content) {
  const keys = new Set();

  // パターン 1: 関数呼び出し形式
  const calls = content.matchAll(TU_KEY_PATTERN);
  for (const m of calls) {
    const key = m[1];
    // 汎用的ではないもの（大文字始まり＋動詞パターン）を TU キーとみなす
    if (
      /^[A-Z]/.test(key) &&
      /_(?:get|search|fetch|list|create|analyze|predict|calculate|classify|extract|validate|convert|download|upload|query|find|submit|check|compare|run|generate|cluster|score|recommend|retrieve|process|compute|detect|identify|filter|map|annotate|enrich|measure|simulate|evaluate|assess|screen|profile|parse|transform|build|train|optimize|select|rank|assign|aggregate|normalize|integrate|visualize|align|assemble|call|index|merge|count|phase|impute|quantify|infer|embed)_/i.test(
        key,
      )
    ) {
      keys.add(key);
    }
  }

  // パターン 2: バッククォート内の TU キー（テーブル内）
  const backtickKeys = content.matchAll(/`([A-Z][a-zA-Z]*_[a-z]+_[a-z_]+)`/g);
  for (const m of backtickKeys) {
    keys.add(m[1]);
  }

  return Array.from(keys);
}

/** スキルのカテゴリ推定 (ディレクトリ名から) */
function inferCategory(skillName) {
  const name = skillName.replace('scientific-', '');
  const categories = {
    A: /literature|citation|crossref|preprint|semantic-scholar|deep-research/,
    B: /statist|eda-|bayesian|data-profiling|geospatial|network-vis|simulation|streaming/,
    C: /ml-|deep-learning|automl|ensemble|transfer|active-learn|explainable|feature-|anomaly|causal-ml|model-monitor|semi-supervised|multi-task|neural-arch|federated|reinforcement|graph-neural/,
    D: /doe|process-opt|adaptive-exp/,
    E: /signal|spectral|biosignal|time-series/,
    F: /chemi|compound|admet|molecular-dock|deep-chem|drug/,
    G: /alphafold|protein-structure|protein-design|structural-proteo|rcsb|protein-domain|protein-interaction/,
    H: /gene-express|rna|expression-comp|geo-|arrayexpress|gtex|encode|regulation|noncoding/,
    I: /variant|gnomad|clingen|civic|pharmaco|precision-onco|rare-disease|cancer/,
    J: /genome-seq|sequence-anal|crispr|ensembl|bioinformatics/,
    K: /proteomics|uniprot|mass-spec/,
    L: /metabol|lipid|glyco/,
    M: /microbiome|metagenom|rrna/,
    N: /publication-fig|present|schematics|latex|interactive-dash|advanced-vis/,
    O: /data-preprocess|missing-data|data-sim|data-sub|pipeline-scaffold|reproducib/,
    P: /pathway|ontology|gene-id|network-analy|reactome|knowledge/,
    Q: /image-analy|advanced-imag|medical-imag|radiology/,
    R: /quantum|symbolic|gpu-single|uncertainty|pca-tsne/,
    S: /clinical|healthcare|survival/,
    T: /environ|marine|paleobio|plant-bio|agricultural|food/,
    U: /epidemiol|public-health|infectious|immunoinfor|pharmacovig/,
    V: /material|computational-mat/,
    W: /robot|lab-auto|lab-data/,
    X: /text-min|academic-writ|critical-rev|peer-review|supplement|revision|paper-qual|systematic|meta-analy|hypothesis|grant|research-method/,
    Y: /education/,
    Z: /regulatory|toxicol|cell-line|model-organ|depmap|nci60|biobank|cohort|population|single-cell|cellxgene|scvi|spatial|squidpy|scatac|human-cell|human-protein|multi-omics|systems-bio|stitch|string-net|opentargets|pharos|pharmgkb|drugbank|ebi-|biothings|hgnc|monarch|disease|gdc|icgc|metabolic-atlas|metabolic-flux|metabolic-model|chembl|biomedical|clinical-trial|clinical-pharm|clinical-decision|neuroscience/,
  };
  for (const [cat, pattern] of Object.entries(categories)) {
    if (pattern.test(name)) return cat;
  }
  return '?';
}

// =============================================================
describe('TU カバレッジ分析', () => {
  const skillDirs = getSkillDirs();
  const linked = [];
  const unlinked = [];

  for (const skillName of skillDirs) {
    const content = readSkillMd(skillName);
    if (!content) {
      unlinked.push({ name: skillName, category: inferCategory(skillName), keys: [] });
      continue;
    }

    const hasTuSection = TU_SECTION_PATTERN.test(content);
    const tuKeys = extractTuKeys(content);

    if (hasTuSection || tuKeys.length > 0) {
      linked.push({ name: skillName, category: inferCategory(skillName), keys: tuKeys });
    } else {
      unlinked.push({ name: skillName, category: inferCategory(skillName), keys: [] });
    }
  }

  const totalTuKeys = new Set(linked.flatMap((s) => s.keys));
  const coverage = ((linked.length / skillDirs.length) * 100).toFixed(1);

  it(`TU カバレッジレポートを生成する (現在: ${coverage}%)`, () => {
    console.log(`\n${'='.repeat(60)}`);
    console.log('📊 SATORI TU カバレッジレポート');
    console.log('='.repeat(60));
    console.log(`\n  総スキル数:      ${skillDirs.length}`);
    console.log(`  TU 連携済み:     ${linked.length}`);
    console.log(`  TU 未連携:       ${unlinked.length}`);
    console.log(`  カバレッジ:      ${coverage}%`);
    console.log(`  ユニーク TU キー: ${totalTuKeys.size}`);

    console.log('\n--- TU 未連携スキル一覧 ---');
    const byCat = {};
    for (const s of unlinked) {
      if (!byCat[s.category]) byCat[s.category] = [];
      byCat[s.category].push(s.name.replace('scientific-', ''));
    }
    for (const [cat, skills] of Object.entries(byCat).sort()) {
      console.log(`\n  [${cat}] (${skills.length} 件)`);
      for (const s of skills) {
        console.log(`    - ${s}`);
      }
    }
    console.log(`\n${'='.repeat(60)}`);

    // 全スキルが TU 連携済み
    expect(linked.length).toBe(skillDirs.length);
  });

  it('TU カバレッジが 100% (v0.26.0: 54.2% → 100.0%)', () => {
    expect(Number(coverage)).toBe(100);
  });
});
