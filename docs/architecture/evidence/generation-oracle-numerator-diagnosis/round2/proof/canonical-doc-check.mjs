import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "../../../../../..");
const source = "docs/architecture/generation-oracle-numerator-diagnosis-correction.typ";
const title = "Three one-diagram QED oracle numerator diagnosis correction";
const report = fs.readFileSync(path.join(root, source), "utf8");
assert.equal(report.match(/^= .+$/gm)?.length, 1);
assert.equal(report.split("\n", 1)[0], `= ${title}`);

const registry = fs.readFileSync(path.join(root, "docs/developers.toml"), "utf8");
const entries = registry.split("[[section.note]]").filter((entry) =>
  entry.includes('id = "generation-oracle-numerator-diagnosis"'),
);
assert.equal(entries.length, 1);
assert.ok(entries[0].includes(`title = "${title}"`));
assert.ok(entries[0].includes(`source = "${source}"`));
assert.ok(entries[0].includes('status = "Round-two correction complete; awaiting new independent Opus A/B review"'));

console.log(JSON.stringify({ schema: 1, title, source, registry_entries: entries.length }, null, 2));
