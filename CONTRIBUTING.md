# Contributing to SATORI

SATORI へのコントリビューションを歓迎します！以下のガイドラインに従ってください。

## 前提条件

- Node.js >= 18
- Git

## セットアップ

```bash
git clone https://github.com/nahisaho/satori.git
cd satori
npm install
```

## 開発フロー

1. `main` ブランチから feature ブランチを作成
2. 変更を実装
3. テストとリントを実行
4. Pull Request を作成

```bash
git checkout -b feature/your-feature-name

# 作業後
npm test              # 全テスト実行
npm run lint          # リントチェック
npm run validate:skills  # SKILL.md フォーマット検証
```

## SKILL.md フォーマット

新しいスキルを追加する場合、以下のフォーマットに準拠してください:

```markdown
---
name: scientific-your-skill-name
description: スキルの概要（1行）
---

# スキルのタイトル

## When to Use

- このスキルを使うべき場面の説明

## Quick Start

```python
# 最小実行例
```
```

### 必須要件

| 要件 | 詳細 |
|---|---|
| ディレクトリ名 | `scientific-` プレフィックス必須 |
| Frontmatter `name` | ディレクトリ名と完全一致 |
| Frontmatter `description` | 必須 |
| `## When to Use` セクション | 必須 |
| `## Quick Start` セクション | 必須 |
| コードブロック | `python`, `markdown`, `json` のいずれか 1 つ以上 |

フォーマット検証は以下で実行できます:

```bash
npm run validate:skills   # 全 190 スキルの自動検証
satori validate --verbose # CLI からの詳細検証
```

## テスト

```bash
npm test                  # 全テスト (1,359 ケース)
npm run test:unit         # CLI ユニットテスト (24 ケース)
npm run test:validation   # SKILL.md 検証テスト (1,335 ケース)
npm run test:coverage     # カバレッジレポート付き
```

### テスト構成

```
tests/
├── unit/
│   └── cli.test.js          # CLI コマンド検証
├── validation/
│   ├── skill-format.test.js  # SKILL.md フォーマット検証
│   └── tu-coverage.test.js   # TU カバレッジ分析
└── integration/              # パイプライン E2E テスト
```

新しい機能を追加する場合は、対応するテストも追加してください。

## コードスタイル

Biome を使用してリント・フォーマットを統一しています:

```bash
npm run lint          # チェックのみ
npm run lint:fix      # 自動修正
npm run format        # フォーマット
```

### ルール概要

- シングルクォート
- セミコロンあり
- トレイリングカンマ (`all`)
- インデント: 2 スペース
- 行幅: 120 文字
- `const` 推奨 (`let` は再代入がある場合のみ)

## コミットメッセージ

以下の形式を推奨:

```
<type>: <summary>

<body>
```

| type | 用途 |
|---|---|
| `feat` | 新機能 |
| `fix` | バグ修正 |
| `docs` | ドキュメント変更 |
| `test` | テスト追加・修正 |
| `chore` | ビルド・設定変更 |
| `refactor` | リファクタリング |

例: `feat: add satori validate CLI command`

## Pull Request

- タイトルに変更の概要を記載
- 関連 Issue がある場合はリンク
- 全テストがパスしていることを確認
- CHANGELOG.md に変更内容を追記

## CI

PR 作成時に以下の GitHub Actions ジョブが自動実行されます:

- **test**: Node.js 18/20/22 でのテスト
- **validate-skills**: SKILL.md フォーマット検証
- **cli-smoke**: CLI コマンドのスモークテスト
- **lint**: Biome リントチェック

全ジョブがパスしないとマージできません。

## ライセンス

MIT License — 詳細は [LICENCE](LICENCE) を参照してください。
