import assert from "node:assert/strict";
import crypto from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.dirname(fileURLToPath(import.meta.url));
const cases = {
  xs_qed_eedd_L1: [0xf5a21f1a84870bc7n, 0xd3ad2f3a0802b6abn],
  xs_qed_eemm_L1: [0x8b7b40bba961420dn, 0x414e2a7a954beae8n],
  xs_qed_eeee_1pi_L1: [0x0db0b7cc52558270n, 0x1fa84d9ad4fdf0b0n],
};
const spellings = new Map([
  [
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q",
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q",
  ],
  [
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::K",
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:k}::K",
  ],
  [
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::P",
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:p}::P",
  ],
  [
    "gammalooprs::{spenso::index}::edge",
    "gammalooprs::{spenso::index,spenso::index-label:e}::edge",
  ],
  [
    "gammalooprs::{spenso::index}::vertex",
    "gammalooprs::{spenso::index,spenso::index-label:v}::vertex",
  ],
]);
const mask = (1n << 64n) - 1n;
const offset = 0xcbf29ce484222325n;
const prime = 0x100000001b3n;

for (const line of fs.readFileSync(path.join(root, "manifest.sha256"), "utf8").trim().split("\n")) {
  const [expected, relative] = line.split(/\s+/, 2);
  const bytes = fs.readFileSync(path.join(root, relative));
  assert.equal(crypto.createHash("sha256").update(bytes).digest("hex"), expected);
}

function read(folder, name) {
  return JSON.parse(fs.readFileSync(path.join(root, folder, `${name}.json`)));
}

function mappedText(value) {
  for (const [current, historical] of spellings) {
    value = value.replaceAll(current, historical);
  }
  return value;
}

function historicalView(value) {
  if (typeof value === "string") return mappedText(value);
  if (Array.isArray(value)) return value.map(historicalView);
  if (value && typeof value === "object") {
    return Object.fromEntries(
      Object.entries(value)
        .filter(([key]) => key !== "num_interface" && !key.endsWith("_interface"))
        .map(([key, child]) => [key, historicalView(child)]),
    );
  }
  return value;
}

function bytesOfNumber(value) {
  const bytes = Buffer.alloc(8);
  bytes.writeBigUInt64LE(value);
  return bytes;
}

function fnv(bytes, initial = offset) {
  let digest = initial;
  for (const byte of bytes) digest = ((digest ^ BigInt(byte)) * prime) & mask;
  return digest;
}

function jsonDigest(value, initial = offset) {
  return fnv(Buffer.from(JSON.stringify(value)), initial);
}

function numeratorDigest(data) {
  const entries = data.diagrams.map((diagram) => {
    const payload = diagram.payload;
    let numerator = offset;
    for (const field of ["numerator", "numerator_prefactor", "projector"]) {
      numerator = jsonDigest(payload[field], numerator);
    }
    for (const field of [
      "num_interface",
      "numerator_prefactor_interface",
      "projector_interface",
    ]) {
      if (field in payload) numerator = jsonDigest(payload[field], numerator);
    }
    const fragments = [...payload.vertices, ...payload.edges.map((edge) => edge[1])];
    for (const fragment of fragments) {
      numerator = jsonDigest(fragment.numerator, numerator);
      if ("num_interface" in fragment) {
        numerator = jsonDigest(fragment.num_interface, numerator);
      }
    }
    const key = jsonDigest(diagram.key);
    return fnv(Buffer.concat([bytesOfNumber(key), bytesOfNumber(numerator)]));
  });
  entries.sort((left, right) => (left < right ? -1 : left > right ? 1 : 0));
  return entries.reduce((digest, entry) => fnv(bytesOfNumber(entry), digest), offset);
}

function ports(metadata, owner) {
  return metadata.diagrams[0].owners.find((entry) => entry.owner === owner).ports;
}

const rows = [];
for (const [name, [expectedHistorical, expectedCurrent]] of Object.entries(cases)) {
  const historical = read("historical", name);
  const current = read("current", name);
  assert.equal(numeratorDigest(historical), expectedHistorical);
  assert.equal(numeratorDigest(current), expectedCurrent);
  assert.deepEqual(historical, historicalView(current));

  const oldMetadata = read("metadata-historical", name);
  const newMetadata = read("metadata-current", name);
  assert.deepEqual(oldMetadata.model, newMetadata.model);
  assert.deepEqual(oldMetadata.numerators, newMetadata.numerators);
  assert.deepEqual(oldMetadata.grouping, newMetadata.grouping);
  assert.equal(oldMetadata.diagrams.length, 1);
  assert.equal(newMetadata.diagrams.length, 1);
  const oldOwners = new Map(oldMetadata.diagrams[0].owners.map((entry) => [entry.owner, entry]));
  let axes = 0;
  const permutations = [];
  for (const currentOwner of newMetadata.diagrams[0].owners) {
    const historicalOwner = oldOwners.get(currentOwner.owner);
    assert.ok(historicalOwner, `${name}: missing historical owner ${currentOwner.owner}`);
    assert.equal(historicalOwner.metric, currentOwner.metric);
    assert.equal(historicalOwner.composite, currentOwner.composite);
    const currentPorts = currentOwner.ports.map(mappedText);
    assert.deepEqual(
      [...historicalOwner.ports].sort(),
      [...currentPorts].sort(),
      `${name}: typed port multiset for ${currentOwner.owner}`,
    );
    assert.equal(new Set(currentPorts).size, currentPorts.length);
    axes += currentPorts.length;
    const permutation = historicalOwner.ports.map((port) => currentPorts.indexOf(port));
    if (!permutation.every((axis, index) => axis === index)) {
      permutations.push({ owner: currentOwner.owner, historical_axis_to_current_axis: permutation });
    }
  }

  const trailingFactorMutation = structuredClone(current);
  trailingFactorMutation.diagrams[0].payload.numerator += "*2";
  assert.notDeepEqual(historical, historicalView(trailingFactorMutation));
  const dimensionMutation = structuredClone(newMetadata);
  const firstPortOwner = dimensionMutation.diagrams[0].owners.find((owner) => owner.ports.length);
  firstPortOwner.ports[0] = firstPortOwner.ports[0].replace("(4,", "(5,");
  assert.notDeepEqual(
    [...ports(oldMetadata, firstPortOwner.owner)].sort(),
    [...ports(dimensionMutation, firstPortOwner.owner).map(mappedText)].sort(),
  );
  const ratioMutation = structuredClone(current);
  ratioMutation.groups[0].members[0].ratio = "2";
  assert.notDeepEqual(historical, historicalView(ratioMutation));

  rows.push({
    case: name,
    historical_digest: `0x${expectedHistorical.toString(16).padStart(16, "0")}`,
    current_digest: `0x${expectedCurrent.toString(16).padStart(16, "0")}`,
    typed_axes: axes,
    permutations,
    normalized_complete_dump_equal: true,
    mutations_rejected: ["trailing_numerator_factor", "typed_dimension", "group_ratio"],
  });
}

console.log(JSON.stringify({ schema: 1, rows }, null, 2));
