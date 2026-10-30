// Rebuild the bounded semantic contract from the immutable round-one captures.
// Validation consumes semantic-contract.json directly and never regenerates it.
import crypto from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const evidence = path.resolve(here, "../..");
const cases = ["xs_qed_eedd_L1", "xs_qed_eemm_L1", "xs_qed_eeee_1pi_L1"];
const read = (folder, name) =>
  JSON.parse(fs.readFileSync(path.join(evidence, folder, `${name}.json`), "utf8"));
const sha = (value) => crypto.createHash("sha256").update(value).digest("hex");

function names(source) {
  const found = new Set();
  let i = 0;
  while (i < source.length) {
    const c = source[i];
    if (!/[A-Za-z_]/.test(c)) {
      i += c === "𝑖" ? "𝑖".length : 1;
      continue;
    }
    let j = i;
    while (j < source.length) {
      if (/[A-Za-z0-9_:]/.test(source[j])) {
        j++;
        continue;
      }
      if (source[j] === "{") {
        let depth = 0;
        do {
          if (source[j] === "{") depth++;
          if (source[j] === "}") depth--;
          j++;
        } while (depth && j < source.length);
        continue;
      }
      break;
    }
    found.add(source.slice(i, j));
    i = j;
  }
  return [...found].sort();
}

function payloadExpressions(payload) {
  return [
    payload.numerator,
    payload.numerator_prefactor,
    payload.projector,
    ...payload.vertices.map((vertex) => vertex.numerator),
    ...payload.edges.map((edge) => edge[1].numerator),
  ];
}

const contract = {
  schema: 1,
  authority: {
    historical_symbols_revision: "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94",
    historical_symbols_sha256: "0a13fee39574d1f8a819cdac3ecc0c9607a266b765907fcae0eaaaf146da174b",
    current_symbols_revision: "b9f788922341e7bde2e801a054423ae70ed7e9ae",
    current_symbols_sha256: "b580b7c70db4f11bf79fd02f8109a715c01b74d63a560a54d49868b8e5c4aa9d",
    constructor_map: {
      Q: [
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q",
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q",
      ],
      K: [
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:k}::K",
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::K",
      ],
      P: [
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:p}::P",
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::P",
      ],
      edge: [
        "gammalooprs::{spenso::index,spenso::index-label:e}::edge",
        "gammalooprs::{spenso::index}::edge",
      ],
      vertex: [
        "gammalooprs::{spenso::index,spenso::index-label:v}::vertex",
        "gammalooprs::{spenso::index}::vertex",
      ],
    },
  },
  cases: {},
  schema_signatures: { payload: {}, metadata: {} },
};

function observeSchema(category, value, pointer = "$", signatures = contract.schema_signatures) {
  if (Array.isArray(value)) {
    for (const child of value) observeSchema(category, child, `${pointer}/*`, signatures);
    return;
  }
  if (!value || typeof value !== "object") return;
  const keys = Object.keys(value).sort();
  const found = signatures[category][pointer] ??= [];
  if (!found.some((signature) => JSON.stringify(signature) === JSON.stringify(keys))) found.push(keys);
  for (const [key, child] of Object.entries(value)) {
    observeSchema(category, child, `${pointer}/${key}`, signatures);
  }
}

for (const name of cases) {
  const historical = read("historical", name);
  const current = read("current", name);
  const historicalMetadata = read("metadata-historical", name);
  const currentMetadata = read("metadata-current", name);
  observeSchema("payload", historical);
  observeSchema("payload", current);
  observeSchema("metadata", historicalMetadata);
  observeSchema("metadata", currentMetadata);
  const payloadHeads = { historical: new Set(), current: new Set() };
  for (const [side, data] of [
    ["historical", historical],
    ["current", current],
  ]) {
    for (const expression of payloadExpressions(data.diagrams[0].payload)) {
      for (const name of names(expression)) payloadHeads[side].add(name);
    }
  }
  const ownerContract = {};
  const currentPayload = current.diagrams[0].payload;
  const interfaces = new Map([
    ...currentPayload.vertices.map((vertex, index) => [`vertex:${index}`, vertex.num_interface]),
    ...currentPayload.edges.map((edge, index) => [`edge:${index}`, edge[1].num_interface]),
    ["numerator", currentPayload.num_interface],
    ["prefactor", currentPayload.numerator_prefactor_interface],
    ["projector", currentPayload.projector_interface],
  ]);
  const historicalOwners = new Map(
    historicalMetadata.diagrams[0].owners.map((owner) => [owner.owner, owner]),
  );
  for (const owner of currentMetadata.diagrams[0].owners) {
    const old = historicalOwners.get(owner.owner);
    ownerContract[owner.owner] = {
      historical_ports: old.ports,
      current_ports: owner.ports,
      historical_flags: [old.metric, old.composite],
      current_flags: [owner.metric, owner.composite],
      current_interface_present: interfaces.get(owner.owner) !== undefined,
    };
  }
  contract.cases[name] = {
    diagram_id: current.diagrams[0].id,
    model_sha256: sha(JSON.stringify(currentMetadata.model)),
    options_sha256: sha(JSON.stringify(currentMetadata.numerators)),
    grouping_sha256: sha(JSON.stringify(currentMetadata.grouping)),
    report: currentMetadata.report,
    payload_names: {
      historical: [...payloadHeads.historical].sort(),
      current: [...payloadHeads.current].sort(),
    },
    owners: ownerContract,
  };
}

for (const category of Object.keys(contract.schema_signatures)) {
  for (const signatures of Object.values(contract.schema_signatures[category])) {
    signatures.sort((left, right) => JSON.stringify(left).localeCompare(JSON.stringify(right)));
  }
}

process.stdout.write(`${JSON.stringify(contract, null, 2)}\n`);
