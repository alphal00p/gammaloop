import assert from "node:assert/strict";
import crypto from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const workspace = path.resolve(here, "../../../..");
const target = path.join(workspace, "target/oracle-numerator-migration");
const cases = {
  xs_qed_eedd_L1: {
    old: "f5a21f1a84870bc7", current: "d3ad2f3a0802b6ab", diagrams: 1, topologies: 17,
    assignments: 34, content: "b467d6be7a12aefa", named: "dd30b2e206f30267", ids: "58ebcb07ba6255d3",
  },
  xs_qed_eemm_L1: {
    old: "8b7b40bba961420d", current: "414e2a7a954beae8", diagrams: 1, topologies: 17,
    assignments: 34, content: "831ae534b022fd0f", named: "f40bfdc791015105", ids: "edd739fe2aef9d55",
  },
  xs_qed_eeee_1pi_L1: {
    old: "0db0b7cc52558270", current: "1fa84d9ad4fdf0b0", diagrams: 1, topologies: 3,
    assignments: 4, content: "8686eb17d2b0f293", named: "f5ac9fd01b9eabd7", ids: "f777e735b39d93bf",
  },
};
const pools = [1, 8, 32];
const sha = (bytes) => crypto.createHash("sha256").update(bytes).digest("hex");

const manifestPath = path.join(here, "manifest.sha256");
if (fs.existsSync(manifestPath)) {
  for (const line of fs.readFileSync(manifestPath, "utf8").trim().split("\n")) {
    const [expected, relative] = line.split(/\s+/, 2);
    assert.equal(sha(fs.readFileSync(path.join(here, relative))), expected, relative);
  }
}
const inventoryPath = path.join(here, "changed-inventory.json");
if (fs.existsSync(inventoryPath)) {
  const inventory = JSON.parse(fs.readFileSync(inventoryPath));
  for (const entry of inventory.files) {
    const bytes = fs.readFileSync(path.join(workspace, entry.path));
    assert.equal(bytes.length, entry.bytes, entry.path);
    assert.equal(sha(bytes), entry.sha256, entry.path);
  }
}

const approval = fs.readFileSync(path.join(here, "review/approved-certificate-review.json"));
assert.equal(sha(approval), "1202230f36622faf1f9dbcf23c29eba8c8ba8c512ae3da14b901194a349e4b30");
assert.equal(JSON.parse(approval).verdict, "APPROVE_CERTIFICATE");

const oldSource = fs.readFileSync(path.join(here, "sources/old-main.rs"), "utf8");
const newSource = fs.readFileSync(path.join(workspace, "crates/feynkit-generator/tests/generation_oracle/main.rs"), "utf8");
assert.equal(sha(Buffer.from(oldSource)), "499b49e4d5b58f9796e79e8393399ddc5e3fee24062fe36ae590dffd7bd18aa8");
assert.equal(sha(Buffer.from(newSource)), "987e9c274c0faa2710689cd10d6b4903414b19e62b9f1833c6c5e3c463b4049b");
const oldLines = oldSource.split("\n"), newLines = newSource.split("\n");
assert.equal(oldLines.length, newLines.length);
const changedLines = [];
for (let index = 0; index < oldLines.length; index++) {
  if (oldLines[index] !== newLines[index]) changedLines.push({ line: index + 1, old: oldLines[index], current: newLines[index] });
}
assert.equal(changedLines.length, 3);
for (const [name, expected] of Object.entries(cases)) {
  const change = changedLines.find((entry) => entry.old.includes(`fast ${name}: Row`));
  assert.ok(change, name);
  assert.equal(change.old.replace(expected.old, expected.current), change.current);
  assert.equal((oldSource.match(new RegExp(expected.old, "g")) ?? []).length, 1);
  assert.equal((newSource.match(new RegExp(expected.current, "g")) ?? []).length, 1);
  assert.equal(newSource.includes(expected.old), false);
}

assert.equal(sha(fs.readFileSync(path.join(target, "old/generation_oracle"))), "c2ef6aea7f9978606ca0891a5061eafaaf75f952574fb5b3894da04a39cbbe84");
assert.equal(sha(fs.readFileSync(path.join(target, "new/generation_oracle"))), "8cde4358a163d8a024bc74bebdcc25c0a4f17f1fe68a79ae557299846316496a");
assert.equal(sha(fs.readFileSync(path.join(target, "bundle/report.pdf"))), "4d8fbf159fafe6ae1f309a8003a76f5de217e24c8c2a32ffee56991fbc311ef1");
assert.equal(sha(fs.readFileSync(path.join(target, "bundle/html/index.html"))), "189ac863098623a2b4bd61fecc16c76e423d447ec34cece8d17b1c13f9b2cb99");

const prime = 0x100000001b3n, mask = 0xffffffffffffffffn, offset = 0xcbf29ce484222325n;
class Digest {
  constructor(value = offset) { this.value = value; }
  bytes(bytes) { let value = this.value; for (const byte of bytes) value = ((value ^ BigInt(byte)) * prime) & mask; return new Digest(value); }
  text(text) { return this.bytes(Buffer.from(text)).bytes([0xff]); }
  name(name) { return name == null ? this.bytes([0xfe]) : this.text(name); }
  number(number) { const bytes = Buffer.alloc(8); bytes.writeBigUInt64LE(number); return this.bytes(bytes); }
  json(value) { return this.bytes(Buffer.from(JSON.stringify(value))); }
}
const multiset = (values) => values.sort((left, right) => left < right ? -1 : left > right ? 1 : 0).reduce((digest, value) => digest.number(value), new Digest()).value;
const hex = (value) => value.toString(16).padStart(16, "0");
function payloadDigests(payload) {
  payload = structuredClone(payload); let numerator = new Digest();
  for (const field of ["numerator", "numerator_prefactor", "projector"]) { numerator = numerator.json(payload[field]); delete payload[field]; }
  for (const field of ["num_interface", "numerator_prefactor_interface", "projector_interface"]) if (field in payload) { numerator = numerator.json(payload[field]); delete payload[field]; }
  const fragments = [];
  const take = (owner) => { const iface = owner.num_interface; delete owner.num_interface; const expression = owner.numerator; owner.numerator = null; return [expression, iface]; };
  for (const vertex of payload.vertices) fragments.push(take(vertex));
  for (const edge of payload.edges) fragments.push(take(edge[1]));
  for (const [expression, iface] of fragments) { numerator = numerator.json(expression); if (iface !== undefined) numerator = numerator.json(iface); }
  for (const field of ["model", "symmetry_factor", "overall_factor"]) delete payload[field];
  return [new Digest().json(payload).value, numerator.value];
}
function summary(data) {
  const content = [], numerators = [], ids = []; let named = new Digest();
  for (const diagram of data.diagrams) {
    const key = new Digest().json(diagram.key).value, [frame, numerator] = payloadDigests(diagram.payload);
    content.push(new Digest().number(key).number(BigInt(diagram.aut)).text(diagram.overall_factor).number(frame).value);
    numerators.push(new Digest().number(key).number(numerator).value);
    const id = BigInt(`0x${diagram.id}`), bytes = Buffer.alloc(16);
    bytes.writeBigUInt64LE(id & mask, 0); bytes.writeBigUInt64LE(id >> 64n, 8);
    ids.push(new Digest().bytes(bytes).value);
    named = named.name(diagram.name).number(key);
  }
  for (const group of data.groups) {
    named = named.number(BigInt(group.master));
    for (const member of group.members) named = named.name(member.name).text(member.ratio).text(member.overall_factor);
  }
  return { diagrams: data.diagrams.length, content: hex(multiset(content)), numerators: hex(multiset(numerators)), named: hex(named.value), ids: hex(multiset(ids)) };
}

const results = [];
let baseline = null;
for (const side of ["old", "new"]) {
  for (const pool of pools) {
    const log = fs.readFileSync(path.join(here, `logs/${side}-pool-${pool}.txt`), "utf8");
    assert.match(log, new RegExp(`pool=${pool}`));
    assert.match(log, /record=unset diagnostic=unset update=unset/);
    assert.match(log, /running 3 tests/);
    assert.match(log, side === "old" ? /exit=101/ : /exit=0/);
    assert.match(log, side === "old" ? /0 passed; 3 failed/ : /3 passed; 0 failed/);
    for (const [name, expected] of Object.entries(cases)) {
      const directory = path.join(target, `${side}/dump-p${pool}`);
      const prefix = `${name}.classes-`;
      const none = fs.readFileSync(path.join(directory, `${prefix}None-t${pool}.json`));
      const auto = fs.readFileSync(path.join(directory, `${prefix}Auto-t${pool}.json`));
      const noneSources = fs.readFileSync(path.join(directory, `${prefix}None-t${pool}.sources.json`));
      const autoSources = fs.readFileSync(path.join(directory, `${prefix}Auto-t${pool}.sources.json`));
      assert.deepEqual(none, auto);
      assert.deepEqual(none, noneSources);
      assert.deepEqual(none, autoSources);
      const found = summary(JSON.parse(none));
      assert.deepEqual(found, {
        diagrams: expected.diagrams, content: expected.content, numerators: expected.current,
        named: expected.named, ids: expected.ids,
      });
      const census = fs.readFileSync(path.join(directory, `${name}.census`));
      if (baseline == null && name === "xs_qed_eedd_L1") baseline = none;
      if (side === "new" || pool !== 1) {
        const reference = fs.readFileSync(path.join(target, `old/dump-p1/${name}.classes-None-t1.json`));
        assert.deepEqual(none, reference);
        assert.deepEqual(census, fs.readFileSync(path.join(target, `old/dump-p1/${name}.census`)));
      }
      results.push({ side, pool, case: name, summary: found, payload_sha256: sha(none), census_sha256: sha(census) });
    }
  }
}

const ownerSource = fs.readFileSync(path.join(workspace, "crates/feynkit-generator/tests/generation_oracle/classes.rs"), "utf8");
assert.match(newSource, /FEYNKIT_ORACLE_THREADS/);
assert.match(ownerSource, /for mode in \[Classes::None, Classes::Auto\]/);
assert.match(ownerSource, /threads\(\)\.into_iter\(\)\.enumerate\(\)/);

const docs = spawnSync("node", [path.join(here, "canonical-doc-check.mjs")], { encoding: "utf8" });
assert.equal(docs.status, 0, docs.stderr);

console.log(JSON.stringify({
  schema: 1,
  approval_sha256: sha(approval),
  changed_lines: changedLines.map(({ line, old, current }) => ({ line, old_sha256: sha(Buffer.from(old)), current_sha256: sha(Buffer.from(current)) })),
  old_elf_sha256: sha(fs.readFileSync(path.join(target, "old/generation_oracle"))),
  new_elf_sha256: sha(fs.readFileSync(path.join(target, "new/generation_oracle"))),
  pools,
  executed_tests: 18,
  results,
  documentation: JSON.parse(docs.stdout),
}, null, 2));
