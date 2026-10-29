import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [historical, current, proof, receipt] = process.argv.slice(2);
if (!historical || !current || !proof || !receipt) {
  throw new Error('usage: node verify-receipt.mjs <historical-capture> <current-capture> <proof-dir> <receipt.json>');
}
const script = path.join(path.dirname(fileURLToPath(import.meta.url)), 'receipt.mjs');
const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'sunset-receipt-'));
try {
  const generated = path.join(scratch, 'receipt.json');
  const run = spawnSync(process.execPath, [script, historical, current, proof, generated],
    {encoding: 'utf8'});
  if (run.status !== 0) throw new Error(run.stderr.slice(0, 1000));
  if (!fs.readFileSync(generated).equals(fs.readFileSync(receipt))) {
    throw new Error('registered receipt no longer binds exact source, binary, proof and capture files');
  }
} finally {
  fs.rmSync(scratch, {recursive: true, force: true});
}
console.log('registered receipt matches a fresh exact artifact inventory');
