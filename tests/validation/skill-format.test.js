/**
 * SKILL.md フォーマット検証テスト
 *
 * 全 190 スキルの SKILL.md が以下を満たすこと検証:
 *   1. YAML Frontmatter (name + description)
 *   2. # タイトル (H1)
 *   3. ## When to Use セクション
 *   4. ## Quick Start セクション
 *   5. Python コードブロックが 1 つ以上
 *   6. ファイル名が SKILL.md
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';

const SKILLS_DIR = path.resolve(__dirname, '../../src/.github/skills');
const EXPECTED_SKILL_COUNT = 190;

/** スキルディレクトリ一覧を取得 */
function getSkillDirs() {
  return fs
    .readdirSync(SKILLS_DIR)
    .filter((d) => d.startsWith('scientific-'))
    .sort();
}

/** SKILL.md を読み込み */
function readSkillMd(skillName) {
  const filePath = path.join(SKILLS_DIR, skillName, 'SKILL.md');
  if (!fs.existsSync(filePath)) return null;
  return fs.readFileSync(filePath, 'utf-8');
}

/** YAML Frontmatter をパース（簡易） */
function parseFrontmatter(content) {
  const match = content.match(/^---\n([\s\S]*?)\n---/);
  if (!match) return null;
  const yaml = match[1];
  const name = yaml.match(/^name:\s*(.+)$/m)?.[1]?.trim();
  const hasDescription = /^description:/m.test(yaml);
  return { name, hasDescription };
}

// =============================================================
describe('SKILL.md — スキル数の検証', () => {
  const skillDirs = getSkillDirs();

  it(`スキルディレクトリが ${EXPECTED_SKILL_COUNT} 個存在する`, () => {
    expect(skillDirs.length).toBe(EXPECTED_SKILL_COUNT);
  });

  it('全ディレクトリが scientific- プレフィックスを持つ', () => {
    for (const d of skillDirs) {
      expect(d).toMatch(/^scientific-/);
    }
  });
});

// =============================================================
describe('SKILL.md — 各スキルのフォーマット検証', () => {
  const skillDirs = getSkillDirs();

  for (const skillName of skillDirs) {
    describe(skillName, () => {
      const content = readSkillMd(skillName);

      it('SKILL.md が存在する', () => {
        expect(content).not.toBeNull();
      });

      if (content === null) return; // ファイルなしなら以降スキップ

      it('YAML Frontmatter に name と description がある', () => {
        const fm = parseFrontmatter(content);
        expect(fm).not.toBeNull();
        expect(fm.name).toBeTruthy();
        expect(fm.hasDescription).toBe(true);
      });

      it('Frontmatter の name がディレクトリ名と一致する', () => {
        const fm = parseFrontmatter(content);
        expect(fm?.name).toBe(skillName);
      });

      it('H1 タイトル (#) が存在する', () => {
        expect(content).toMatch(/^# .+$/m);
      });

      it('## When to Use セクションがある', () => {
        expect(content).toMatch(/^## When to Use/m);
      });

      it('## Quick Start セクションがある', () => {
        expect(content).toMatch(/^## Quick Start/m);
      });

      it('コードブロックが 1 つ以上含まれる (python/markdown/json)', () => {
        expect(content).toMatch(/```(?:python|markdown|json)/);
      });
    });
  }
});

// =============================================================
describe('SKILL.md — TU 連携の統計', () => {
  const skillDirs = getSkillDirs();

  it('TU 参照を持つスキルが 80 以上ある', () => {
    let tuCount = 0;
    for (const skillName of skillDirs) {
      const content = readSkillMd(skillName);
      if (!content) continue;
      // ToolUniverse 参照パターン（docstring 内 / テーブル内 / セクション見出し）
      const hasTU =
        /ToolUniverse/i.test(content) ||
        /`[A-Z][a-zA-Z]+_(?:get|search|fetch|list|create|analyze|predict|calculate|classify|extract|validate|convert|download|upload|query|find|submit|check|compare)_/i.test(
          content,
        );
      if (hasTU) tuCount++;
    }
    expect(tuCount).toBeGreaterThanOrEqual(80);
  });
});
