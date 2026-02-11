#!/usr/bin/env node

const fs = require('fs');
const path = require('path');

const COMMAND = process.argv[2];
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
  satori help                         Show this help message
  satori --version, -v                Show version number

Options:
  --force     Overwrite existing .github/ directory
  --dry-run   Preview what would be installed without making changes
`);
}

function showVersion() {
  const pkg = require(path.join(PACKAGE_ROOT, 'package.json'));
  console.log(pkg.version);
}

switch (COMMAND) {
  case 'init':
    init();
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
