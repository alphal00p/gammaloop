import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {spawnSync} from 'node:child_process';

const workspace = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../../../../../..');
// Source43 gives Q the presentation tag tensor-label:q. The current owner uses
// the same rank-one Q with the same component labels and arguments, without that
// display tag. These presentation-only serialization aliases has identity component/port action.
const legacyQ = 'gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q';
const currentQ = 'gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q';
export function liveModel(capture) {
  const binary = process.env.SUNSET_LIVE_RULE_PROBE ??
    path.join(workspace, 'target/debug/examples/sunset_live_rule_probe');
  const run = spawnSync(binary, [path.join(capture, 'post-import-model.json'), '--fingerprint'],
    {encoding: 'utf8', env: process.env, timeout: 120000});
  const fingerprint = run.stdout?.trim();
  const settings = JSON.parse(fs.readFileSync(path.join(capture, 'post-import-settings.json')));
  if (run.status !== 0 || !/^[0-9a-f]{64}$/.test(fingerprint) ||
      settings.model_fingerprint !== fingerprint) {
    throw new Error(`complete live-model owner differs from captured fingerprint: ${run.stderr}`);
  }
  return fingerprint;
}
export function liveRule(capture, files) {
  const binary = process.env.SUNSET_LIVE_RULE_PROBE ??
    path.join(workspace, 'target/debug/examples/sunset_live_rule_probe');
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'sunset-live-rule-'));
  try {
    let aliasedFiles = 0;
    const adapted = files.map((file, index) => {
      const source = fs.readFileSync(file, 'utf8');
      const value = source.replaceAll(legacyQ, currentQ)
        .replaceAll('gammalooprs::{spenso::index,spenso::index-label:e}::edge',
          'gammalooprs::{spenso::index}::edge')
        .replaceAll('gammalooprs::{spenso::index,spenso::index-label:v}::vertex',
          'gammalooprs::{spenso::index}::vertex');
      if (source !== value) aliasedFiles++;
      const destination = path.join(scratch, `${index}.json`);
      fs.writeFileSync(destination, value);
      return destination;
    });
    const args = [path.join(capture, 'post-import-model.json'),
      path.join(capture, 'post-import-settings.json'), ...adapted];
    const run = spawnSync(binary, args, {encoding: 'utf8', env: process.env,
      maxBuffer: 1024 * 1024, timeout: 120000});
    if (run.status !== 0) {
      throw new Error(`production live-rule owner rejected capture: ${run.error ?? ''}\n` +
        adapted.reduce((text, file, index) => text.replaceAll(file, `owner-file-${index}`), run.stderr.slice(-20000)));
    }
    return {binary, files: files.length, presentationIdentityAliasFiles: aliasedFiles,
      result: adapted.reduce((text, file, index) => text.replaceAll(file, `owner-file-${index}`), run.stdout.trim())};
  } finally {
    fs.rmSync(scratch, {recursive: true, force: true});
  }
}

export function liveCards(capture, files) {
  const binary = process.env.SUNSET_LIVE_RULE_PROBE ??
    path.join(workspace, 'target/debug/examples/sunset_live_rule_probe');
  const run = spawnSync(binary, [path.join(capture, 'post-import-model.json'), '--cards', ...files],
    {encoding: 'utf8', env: process.env, timeout: 120000});
  if (run.status !== 0) throw new Error(`live external parameter card rejected: ${run.stderr}`);
  return {files: files.length, result: run.stdout.trim()};
}

export function liveGlobals(rows) {
  const binary = process.env.SUNSET_LIVE_RULE_PROBE ??
    path.join(workspace, 'target/debug/examples/sunset_live_rule_probe');
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'sunset-global-owners-'));
  try {
    const file = path.join(scratch, 'globals.json');
    fs.writeFileSync(file, JSON.stringify(rows));
    const run = spawnSync(binary, ['--globals', file],
      {encoding: 'utf8', env: process.env, timeout: 120000});
    if (run.status !== 0) throw new Error(`global tensor owner rejected: ${run.stderr}`);
    return {owners: rows.length, result: run.stdout.trim()};
  } finally {fs.rmSync(scratch, {recursive: true, force: true});}
}

export function liveDots(files) {
  const binary = process.env.SUNSET_LIVE_RULE_PROBE ??
    path.join(workspace, 'target/debug/examples/sunset_live_rule_probe');
  const run = spawnSync(binary, ['--dot', ...files],
    {encoding: 'utf8', env: process.env, maxBuffer: 16 * 1024 * 1024, timeout: 120000});
  if (run.status !== 0) throw new Error(`production DOT parser rejected capture: ${run.stderr}`);
  const parsed = JSON.parse(run.stdout);
  if (!Array.isArray(parsed) || parsed.length !== files.length) {
    throw new Error('production DOT parser returned an incomplete graph inventory');
  }
  return parsed;
}
