import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const workspace = path.resolve(here, "../../../..");
const source = "docs/architecture/generation-oracle-numerator-migration.typ";
const title = "Three one-diagram QED numerator oracle migration";
const report = fs.readFileSync(path.join(workspace, source), "utf8");
assert.equal(report.match(/^= .+$/gm)?.length, 1);
assert.equal(report.split("\n", 1)[0], `= ${title}`);
const registry = fs.readFileSync(path.join(workspace, "docs/developers.toml"), "utf8");
const entries = registry.split("[[section.note]]").filter((entry) => entry.includes('id = "generation-oracle-numerator-migration"'));
assert.equal(entries.length, 1);
assert.ok(entries[0].includes(`title = "${title}"`));
assert.ok(entries[0].includes(`source = "${source}"`));
assert.ok(entries[0].includes('status = "Migration candidate complete; awaiting new independent Opus migration/test review"'));
console.log(JSON.stringify({ schema: 1, title, source, registry_entries: entries.length }, null, 2));
