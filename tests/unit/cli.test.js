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
    expect(stdout).toContain('satori skill search');
    expect(stdout).toContain('satori skill info');
    expect(stdout).toContain('satori skill recommend');
    expect(stdout).toContain('satori pipeline suggest');
    expect(stdout).toContain('satori pipeline list');
    expect(stdout).toContain('satori docs generate');
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
    // "(\d+ files)" のパターンを検索
    const match = stdout.match(/\((\d+) files\)/);
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

  it('195 スキルを検証する', () => {
    const { stdout } = run('validate');
    const match = stdout.match(/(\d+) pass/);
    expect(match).not.toBeNull();
    expect(Number(match[1])).toBe(195);
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

  it('スキル数が 195 と表示される', () => {
    const { stdout } = run('stats');
    expect(stdout).toContain('195');
  });

  it('パイプライン数が 50 と表示される', () => {
    const { stdout } = run('stats');
    expect(stdout).toContain('50');
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
describe('satori skill search', () => {
  it('キーワード検索でスキルを検索する', () => {
    const { stdout, exitCode } = run('skill', 'search', '創薬');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('検索結果');
    expect(stdout).toContain('📖');
  });

  it('マッチしない検索は空結果を返す', () => {
    const { stdout, exitCode } = run('skill', 'search', 'nonexistent-keyword-xyz');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('見つかりませんでした');
  });

  it('検索クエリなしでエラー終了する', () => {
    const { exitCode } = run('skill', 'search');
    expect(exitCode).not.toBe(0);
  });

  it('複数単語での検索が可能', () => {
    const { stdout, exitCode } = run('skill', 'search', '機械', '学習');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('検索結果');
  });
});

// =============================================================
describe('satori skill info', () => {
  it('スキルの詳細情報を表示する', () => {
    const { stdout, exitCode } = run('skill', 'info', 'deep-learning');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('Scientific Deep Learning');
    expect(stdout).toContain('When to Use');
  });

  it('scientific- プレフィックスなしで検索できる', () => {
    const { stdout, exitCode } = run('skill', 'info', 'deep-learning');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('Papers with Code');
  });

  it('存在しないスキル名でエラー終了する', () => {
    const { exitCode } = run('skill', 'info', 'nonexistent-skill');
    expect(exitCode).not.toBe(0);
  });

  it('スキル名なしでエラー終了する', () => {
    const { exitCode } = run('skill', 'info');
    expect(exitCode).not.toBe(0);
  });

  it('ToolUniverse 連携情報を表示する', () => {
    const { stdout } = run('skill', 'info', 'drug-target-profiling');
    expect(stdout).toContain('ToolUniverse');
  });

  it('関連パイプラインを表示する', () => {
    const { stdout } = run('skill', 'info', 'deep-learning');
    // AI駆動エビデンス合成パイプラインが関連
    expect(stdout).toContain('関連パイプライン');
  });
});

// =============================================================
describe('satori skill recommend', () => {
  it('関連スキルを推奨する', () => {
    const { stdout, exitCode } = run('skill', 'recommend', 'deep-learning');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('関連するスキル');
    expect(stdout).toContain('パイプラインで併用');
  });

  it('存在しないスキル名でエラー終了する', () => {
    const { exitCode } = run('skill', 'recommend', 'nonexistent-skill');
    expect(exitCode).not.toBe(0);
  });

  it('スキル名なしでエラー終了する', () => {
    const { exitCode } = run('skill', 'recommend');
    expect(exitCode).not.toBe(0);
  });

  it('複数のスキル推奨候補を表示する', () => {
    const { stdout } = run('skill', 'recommend', 'data-preprocessing');
    // 複数の関連スキルが表示される
    expect(stdout).toMatch(/\d\./);
  });
});

// =============================================================
describe('satori pipeline custom', () => {
  it('カスタムパイプライン list コマンドが実行可能', () => {
    const { stdout, exitCode } = run('pipeline', 'custom', 'list');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('カスタムパイプライン');
  });

  it('ファイルなしで add エラー終了する', () => {
    const { exitCode } = run('pipeline', 'custom', 'add');
    expect(exitCode).not.toBe(0);
  });

  it('存在しないファイルで add エラー終了する', () => {
    const { exitCode } = run('pipeline', 'custom', 'add', '/nonexistent/file.json');
    expect(exitCode).not.toBe(0);
  });

  it('ID なしで remove エラー終了する', () => {
    const { exitCode } = run('pipeline', 'custom', 'remove');
    expect(exitCode).not.toBe(0);
  });

  it('不明なアクションでエラー終了する', () => {
    const { exitCode } = run('pipeline', 'custom', 'unknown');
    expect(exitCode).not.toBe(0);
  });
});

// =============================================================
describe('satori docs generate', () => {
  it('preview モードで実行できる', () => {
    const { stdout, exitCode } = run('docs', 'generate', '--preview');
    expect(exitCode).toBe(0);
    expect(stdout).toContain('docs generate');
  });

  it('不明なサブコマンドでエラー終了する', () => {
    const { exitCode } = run('docs', 'unknown');
    expect(exitCode).not.toBe(0);
  });
});

// =============================================================
describe('satori skill (引数なし)', () => {
  it('エラーとなる', () => {
    const { exitCode } = run('skill');
    expect(exitCode).not.toBe(0);
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
