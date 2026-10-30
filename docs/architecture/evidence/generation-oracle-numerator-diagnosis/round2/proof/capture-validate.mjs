// Immutable capture/FNV validation. Semantic equivalence is intentionally owned
// by semantic-validate.mjs and is not inferred from these hashes.
import assert from "node:assert/strict";
import crypto from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const evidence = path.resolve(here, "../..");
const cases = {
  xs_qed_eedd_L1: [0xf5a21f1a84870bc7n, 0xd3ad2f3a0802b6abn],
  xs_qed_eemm_L1: [0x8b7b40bba961420dn, 0x414e2a7a954beae8n],
  xs_qed_eeee_1pi_L1: [0x0db0b7cc52558270n, 0x1fa84d9ad4fdf0b0n],
};
const prime = 0x100000001b3n;
const mask = 0xffffffffffffffffn;

class Digest {
  constructor() { this.value = 0xcbf29ce484222325n; }
  bytes(bytes) {
    for (const byte of bytes) this.value = ((this.value ^ BigInt(byte)) * prime) & mask;
    return this;
  }
  json(value) { return this.bytes(new TextEncoder().encode(JSON.stringify(value))); }
  number(value) {
    const bytes = new Uint8Array(8);
    new DataView(bytes.buffer).setBigUint64(0, value, true);
    return this.bytes(bytes);
  }
}

function numeratorDigest(data) {
  const entries = data.diagrams.map((diagram) => {
    const payload = diagram.payload, numerator = new Digest();
    for (const field of ["numerator", "numerator_prefactor", "projector"]) numerator.json(payload[field]);
    for (const field of ["num_interface", "numerator_prefactor_interface", "projector_interface"]) {
      if (Object.hasOwn(payload, field)) numerator.json(payload[field]);
    }
    for (const fragment of [...payload.vertices, ...payload.edges.map((edge) => edge[1])]) {
      numerator.json(fragment.numerator);
      if (Object.hasOwn(fragment, "num_interface")) numerator.json(fragment.num_interface);
    }
    return new Digest().number(new Digest().json(diagram.key).value).number(numerator.value).value;
  }).sort((left, right) => left < right ? -1 : left > right ? 1 : 0);
  return entries.reduce((digest, entry) => digest.number(entry), new Digest()).value;
}

for (const line of fs.readFileSync(path.join(evidence, "manifest.sha256"), "utf8").trim().split("\n")) {
  const [expected, relative] = line.split(/\s+/, 2);
  const bytes = fs.readFileSync(path.join(evidence, relative));
  assert.equal(crypto.createHash("sha256").update(bytes).digest("hex"), expected, relative);
}

const rows = [];
for (const [name, expected] of Object.entries(cases)) {
  for (const [index, side] of ["historical", "current"].entries()) {
    const raw = fs.readFileSync(path.join(evidence, side, `${name}.json`));
    const provenance = path.join(evidence, "round2", "provenance", side);
    const none = fs.readFileSync(path.join(provenance, `${name}.classes-None-t1.json`));
    const auto = fs.readFileSync(path.join(provenance, `${name}.classes-Auto-t1.json`));
    assert.deepEqual(raw, none, `${side}/${name}: round-one payload != fresh Classes::None`);
    assert.deepEqual(none, auto, `${side}/${name}: Classes::None != Classes::Auto`);
    const found = numeratorDigest(JSON.parse(raw));
    assert.equal(found, expected[index], `${side}/${name}: numerator digest`);
    const census = fs.readFileSync(path.join(provenance, `${name}.census`), "utf8");
    assert.match(census, /completed: true/);
    assert.match(census, /D0 key=/);
    rows.push({
      side,
      case: name,
      digest: `0x${found.toString(16).padStart(16, "0")}`,
      none_auto_equal: true,
      payload_sha256: crypto.createHash("sha256").update(raw).digest("hex"),
    });
  }
}

console.log(JSON.stringify({ schema: 1, rows }, null, 2));
