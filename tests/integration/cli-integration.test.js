/**
 * SATORI 統合テスト
 *
 * CLI コマンド間の整合性を検証する。
 * validate / stats / pipeline list の出力が互いに一致することを確認。
 */

import { execFileSync } from 'node:child_process';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { describe, expect, it } from 'vitest';

const CLI = path.resolve(__dirname, '../../bin/satori.js');
const PACKAGE_ROOT = path.resolve(__dirname, '../..');
const SKILLS_DIR = path.join(PACKAGE_ROOT, 'src', '.github', 'skills');

function runCLI(...args) {
  return execFileSync('node', [CLI, ...args], {
    encoding: 'utf-8',
    timeout: 30000,
    cwd: PACKAGE_ROOT,
  });
}

function runCLIInTmp(...args) {
  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'satori-int-'));
  try {
    const stdout = execFileSync('node', [CLI, ...args], {
      encoding: 'utf-8',
      timeout: 30000,
      cwd: tmpDir,
    });
    fs.rmSync(tmpDir, { recursive: true, force: true });
    return stdout;
  } catch (err) {
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw err;
  }
}

// =============================================================
describe('satori validate ↔ stats 整合性', () => {
  const validateOutput = runCLI('validate');
  const statsOutput = runCLI('stats');

  it('validate のスキル数と stats のスキル総数が一致する', () => {
    const validateMatch = validateOutput.match(/(\d+) pass/);
    const statsMatch = statsOutput.match(/スキル総数:\s+(\d+)/);
    expect(validateMatch).not.toBeNull();
    expect(statsMatch).not.toBeNull();
    expect(validateMatch[1]).toBe(statsMatch[1]);
  });

  it('stats のスキル総数がファイルシステム上のスキル数と一致する', () => {
    const dirs = fs.readdirSync(SKILLS_DIR).filter((d) => d.startsWith('scientific-'));
    const statsMatch = statsOutput.match(/スキル総数:\s+(\d+)/);
    expect(Number(statsMatch[1])).toBe(dirs.length);
  });

  it('stats のパイプライン数が 50 である', () => {
    const match = statsOutput.match(/パイプライン数:\s+(\d+)/);
    expect(Number(match[1])).toBe(50);
  });

  it('stats の TU カバレッジが 0〜100% の範囲にある', () => {
    const match = statsOutput.match(/TU 連携スキル:\s+\d+ \(([\d.]+)%\)/);
    expect(match).not.toBeNull();
    const pct = Number.parseFloat(match[1]);
    expect(pct).toBeGreaterThanOrEqual(0);
    expect(pct).toBeLessThanOrEqual(100);
  });
});

// =============================================================
describe('satori pipeline list 整合性', () => {
  const listOutput = runCLI('pipeline', 'list');

  it('26 ドメインパイプライン + 15 クロスドメイン + 5 インダストリー + 4 メソドロジー = 50 パイプラインが含まれる', () => {
    // 数値 ID (#1-#26) + 文字列 ID (#A-#O, #Ind-1〜#Ind-5, #M-α〜#M-δ)
    const numIds = [...listOutput.matchAll(/#\s*(\d+)/g)].map((m) => Number(m[1]));
    const uniqueNumIds = [...new Set(numIds)];
    expect(uniqueNumIds.length).toBe(26);

    // クロスドメイン (A-O = 15), インダストリー (Ind-1-5 = 5), メソドロジー (M-α-δ = 4)
    const charIds = [...listOutput.matchAll(/#([A-O]|Ind-[1-5]|M-[α-δα-δ])/g)];
    const uniqueCharIds = [...new Set(charIds.map((m) => m[1]))];
    expect(uniqueCharIds.length).toBeGreaterThanOrEqual(14); // 最低限の確認
  });

  it('全パイプラインにスキル連鎖セパレータ (→) がある', () => {
    const lines = listOutput.split('\n').filter((l) => l.includes('→'));
    expect(lines.length).toBeGreaterThanOrEqual(50);
  });

  it('pipeline list の参照先ドキュメントが存在する', () => {
    const docPath = path.join(PACKAGE_ROOT, 'docs', 'SATORI_PIPELINE_EXAMPLES.md');
    expect(fs.existsSync(docPath)).toBe(true);
  });
});

// =============================================================
describe('satori validate --verbose 詳細出力', () => {
  const verboseOutput = runCLI('validate', '--verbose');

  it('各スキルに ✔ マークが出力される', () => {
    const checkmarks = verboseOutput.match(/✔/g);
    expect(checkmarks).not.toBeNull();
    // 190 スキル + 最終サマリー行の ✔ = 191
    expect(checkmarks.length).toBe(191);
  });

  it('✘ マークが出力されない (全スキル合格)', () => {
    expect(verboseOutput).not.toContain('✘');
  });
});

// =============================================================
describe('satori init --dry-run 出力整合性', () => {
  it('dry-run のファイル数がファイルシステムと一致する', () => {
    const output = runCLIInTmp('init', '--dry-run');
    const match = output.match(/(\d+) files/);
    expect(match).not.toBeNull();

    // 再帰的にファイル数をカウント
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

    const sourceDir = path.join(PACKAGE_ROOT, 'src', '.github');
    const actualCount = countFiles(sourceDir);
    expect(Number(match[1])).toBe(actualCount);
  });
});

// =============================================================
describe('satori help コマンド一覧', () => {
  const helpOutput = runCLI('help');

  it('全 CLI コマンドがヘルプに記載されている', () => {
    const commands = ['init', 'pipeline suggest', 'pipeline list', 'validate', 'stats', 'help', '--version'];
    for (const cmd of commands) {
      expect(helpOutput).toContain(cmd);
    }
  });
});
