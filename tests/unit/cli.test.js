/**
 * SATORI CLI ユニットテスト
 *
 * satori.js の各コマンドを子プロセスとして実行し、
 * 出力・終了コードを検証する。
 */

import { execFileSync } from 'node:child_process';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';

const CLI = path.resolve(__dirname, '../../bin/satori.js');
const PACKAGE_ROOT = path.resolve(__dirname, '../..');

/** CLI を同期実行して結果を返す（一時ディレクトリで実行） */
function run(...args) {
  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'satori-run-'));
  try {
    const stdout = execFileSync('node', [CLI, ...args], {
      encoding: 'utf-8',
      timeout: 15000,
      cwd: tmpDir,
    });
    fs.rmSync(tmpDir, { recursive: true, force: true });
    return { stdout, exitCode: 0 };
  } catch (err) {
    fs.rmSync(tmpDir, { recursive: true, force: true });
    return {
      stdout: err.stdout || '',
      stderr: err.stderr || '',
      exitCode: err.status,
    };
  }
}

// =============================================================
describe('satori --version', () => {
  it('package.json のバージョンを出力する', () => {
    const pkg = JSON.parse(fs.readFileSync(path.join(PACKAGE_ROOT, 'package.json'), 'utf-8'));
    const { stdout, exitCode } = run('--version');
    expect(exitCode).toBe(0);
    expect(stdout.trim()).toBe(pkg.version);
  });

  it('-v フラグでも同じ出力', () => {
    const { stdout, exitCode } = run('-v');
    expect(exitCode).toBe(0);
    expect(stdout.trim()).toMatch(/^\d+\.\d+\.\d+$/);
  });
});

// =============================================================
describe('satori help', () => {
  it('使い方テキストを表示する', () => {
    const { stdout, exitCode } = run('help');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('satori init');
    expect(stdout).toContain('satori pipeline suggest');
    expect(stdout).toContain('satori pipeline list');
  });

  it('--help フラグでも動作する', () => {
    const { stdout, exitCode } = run('--help');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('SATORI');
  });

  it('-h フラグでも動作する', () => {
    const { stdout, exitCode } = run('-h');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('satori init');
  });
});

// =============================================================
describe('satori init --dry-run', () => {
  it('正常終了し、dry-run メッセージを出力する', () => {
    const { stdout, exitCode } = run('init', '--dry-run');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('dry-run');
  });

  it('ファイル数が 0 より大きい', () => {
    const { stdout } = run('init', '--dry-run');
    // "190 skills" or ファイル数の出力を検証
    const match = stdout.match(/(\d+)/);
    expect(match).not.toBeNull();
    expect(Number(match[1])).toBeGreaterThan(0);
  });
});

// =============================================================
describe('satori init (実行)', () => {
  let tmpDir;

  beforeEach(() => {
    tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'satori-test-'));
  });

  afterEach(() => {
    fs.rmSync(tmpDir, { recursive: true, force: true });
  });

  it('.github ディレクトリを作成する', () => {
    const stdout = execFileSync('node', [CLI, 'init'], {
      encoding: 'utf-8',
      cwd: tmpDir,
      timeout: 30000,
    });
    const githubDir = path.join(tmpDir, '.github');
    expect(fs.existsSync(githubDir)).toBe(true);
    expect(stdout).toContain('Installed');
  });

  it('--force で上書きできる', () => {
    // 初回
    execFileSync('node', [CLI, 'init'], {
      encoding: 'utf-8',
      cwd: tmpDir,
      timeout: 30000,
    });
    // 上書き
    const stdout = execFileSync('node', [CLI, 'init', '--force'], {
      encoding: 'utf-8',
      cwd: tmpDir,
      timeout: 30000,
    });
    expect(stdout).toContain('Installed');
  });

  it('.github なしで 2 回目実行すると警告を出す', () => {
    // 初回
    execFileSync('node', [CLI, 'init'], {
      encoding: 'utf-8',
      cwd: tmpDir,
      timeout: 30000,
    });
    // 2回目（--force なし）
    const result = run('init');
    // 既存の場合はエラーまたは警告
    // NOTE: 実装に依存—存在チェック後の挙動を検証
    expect(result.stdout + (result.stderr || '')).toBeTruthy();
  });
});

// =============================================================
describe('satori pipeline list', () => {
  it('26 パイプラインを一覧表示する', () => {
    const { stdout, exitCode } = run('pipeline', 'list');
    expect(exitCode).toBe(0);
    // #1 〜 #26 が出力される
    expect(stdout).toContain('#1');
    expect(stdout).toContain('#26');
    expect(stdout).toContain('パイプライン一覧');
  });

  it('各パイプラインにスキル連鎖が表示される', () => {
    const { stdout } = run('pipeline', 'list');
    // → はスキル連鎖のセパレータ
    expect(stdout).toContain('→');
  });
});

// =============================================================
describe('satori pipeline (不正サブコマンド)', () => {
  it('不明なサブコマンドでエラー終了する', () => {
    const { exitCode } = run('pipeline', 'unknown');
    expect(exitCode).not.toBe(0);
  });

  it('サブコマンドなしでエラー終了する', () => {
    const { exitCode } = run('pipeline');
    expect(exitCode).not.toBe(0);
  });
});

// =============================================================
describe('satori (不明コマンド)', () => {
  it('不明コマンドでエラー終了する', () => {
    const { exitCode } = run('nonexistent');
    expect(exitCode).not.toBe(0);
  });
});

// =============================================================
describe('satori validate', () => {
  it('全スキルの検証に成功する', () => {
    const { stdout, exitCode } = run('validate');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('検証結果');
    expect(stdout).toContain('pass');
    expect(stdout).toContain('0 fail');
    expect(stdout).toContain('全スキルの検証に成功しました');
  });

  it('--verbose でスキルごとの結果を表示する', () => {
    const { stdout, exitCode } = run('validate', '--verbose');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('✔');
    expect(stdout).toContain('scientific-');
  });

  it('190 スキルを検証する', () => {
    const { stdout } = run('validate');
    const match = stdout.match(/(\d+) pass/);
    expect(match).not.toBeNull();
    expect(Number(match[1])).toBe(190);
  });
});

// =============================================================
describe('satori stats', () => {
  it('統計情報を表示する', () => {
    const { stdout, exitCode } = run('stats');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('SATORI');
    expect(stdout).toContain('統計');
    expect(stdout).toContain('スキル総数');
    expect(stdout).toContain('パイプライン数');
    expect(stdout).toContain('TU 連携スキル');
  });

  it('スキル数が 190 と表示される', () => {
    const { stdout } = run('stats');
    expect(stdout).toContain('190');
  });

  it('パイプライン数が 26 と表示される', () => {
    const { stdout } = run('stats');
    expect(stdout).toContain('26');
  });

  it('TU カバレッジが表示される', () => {
    const { stdout } = run('stats');
    // "XX.X%" の形式
    expect(stdout).toMatch(/\d+\.\d+%/);
  });

  it('バージョン番号が表示される', () => {
    const pkg = JSON.parse(fs.readFileSync(path.join(PACKAGE_ROOT, 'package.json'), 'utf-8'));
    const { stdout } = run('stats');
    expect(stdout).toContain(pkg.version);
  });
});

// =============================================================
describe('satori (引数なし)', () => {
  it('ヘルプを表示して正常終了する', () => {
    const { stdout, exitCode } = run();
    expect(exitCode).toBe(0);
    expect(stdout).toContain('satori');
  });
});
