#!/usr/bin/env node

import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";

const sha256 = (path) => createHash("sha256").update(readFileSync(path)).digest("hex");
const text = (path) => readFileSync(path, "utf8");
const json = (path) => JSON.parse(text(path));
const equal = (actual, expected, label) => {
  if (JSON.stringify(actual) !== JSON.stringify(expected)) throw new Error(`${label}: ${JSON.stringify(actual)} != ${JSON.stringify(expected)}`);
};

const paths = {
  before: ".local/review-main-core-test-inference-opus-20261029T190454Z-1747112/scopes-004.txt",
  after: "target/wp3-core-test-inference-round2/scopes-after.txt",
  invalid: "target/wp3-core-test-inference-round2/scopes-before.txt",
  oldManifest: "target/wp3-core-test-inference-round2/manifest.sha256",
  independentBefore: ".local/review-main-core-scope-round2-opus-20261029T191535Z-2142322/scopes-004d82f3.json",
  independentAfter: ".local/review-main-core-scope-round2-opus-20261029T191535Z-2142322/scopes-5126de4f.json",
};
equal(sha256(paths.before), "1baa6025891913837b2f9e34a86041622e1f691b997f9309a4952887f31fc375", "genuine before hash");
equal(sha256(paths.after), "7c65577d308b5dceaa60bc88bd01cb325b827aa3fbb734a154ff8bd9e78eefb7", "after hash");
equal(sha256(paths.invalid), "782a0115d8504e6c4c8a8bf45bf0557b2ce2ef503a6e57f79e5ae4169eac2f2b", "invalid artifact hash");
equal(sha256(paths.oldManifest), "fbae53a7b0e05dd77a00b7144089f366a3518d141bdf083d404eea73a5eebb56", "failed round2 manifest hash");
if (!text(paths.before).includes("scopes 110 bad 2") || !text(paths.before).includes("MISMATCH generation-sampling-staged") || !text(paths.before).includes("MISMATCH gammaloop-architecture")) throw new Error("genuine before content mismatch");
if (!text(paths.after).includes("scopes 110 bad 0") || text(paths.after).includes("MISMATCH")) throw new Error("after content mismatch");
if (!text(paths.invalid).endsWith("scopes 110 bad 0\n110 bad 2\n")) throw new Error("invalid artifact characterization changed");

const before = json(paths.independentBefore);
const after = json(paths.independentAfter);
equal(sha256(paths.independentBefore), "80200d9ab34e59a1cfe772fd32dba7739ad969be20d69f35a068773a74299ed0", "independent before hash");
equal(sha256(paths.independentAfter), "93a92627cb7a110d6a1a6f34db29e1ba107561b84bd96bfa47e660d8129486d1", "independent after hash");
equal([before.scopes, before.bad], [110, 2], "before counts");
equal([after.scopes, after.bad], [110, 0], "after counts");
const bad = before.bad_rows.map((row) => ({note:row.note, symbol:row.symbol, expected:row.expected, actual:row.actual})).sort((a, b) => a.note.localeCompare(b.note));
equal(bad, [
  {note:"gammaloop-architecture",symbol:"gammalooprs::integrands::process::ProcessIntegrand::warm_up",expected:"b1711442cd184b4e9776b84d92c313d85f70b7b0",actual:"2f6aea58d34fcb798c742d41feac8240dfc5ae55"},
  {note:"generation-sampling-staged",symbol:"graph-derived reflection controls",expected:"a07a79ea2a2ab2ba9e0544357f113fd7ceaf2a60",actual:"a5d00202b935534a9afd09dd37a5d7484ccbfea9"},
], "before mismatch set");
const key = (row) => [row.note,row.path,row.symbol,row.anchor,row.end_anchor];
equal(before.rows.map(key), after.rows.map(key), "ordered scope keys");
const changedExpected = before.rows.map((row, index) => row.expected === after.rows[index].expected ? null : after.rows[index].note).filter(Boolean).sort();
equal(changedExpected, ["gammaloop-architecture","generation-sampling-staged"], "changed expected digests");
equal(sha256("docs/developers.toml"), "9a086d4d7e9936f145b827e1f79109f3dfa77e16ce9ee522129f7aff8ff4bb71", "current registry bytes");

const receipt = json("docs/architecture/evidence/generation-integration/wp3/core-test-inference/receipt.json");
const correction = json("docs/architecture/evidence/generation-integration/wp3/core-test-inference/round2-scope-correction.json");
equal(receipt.static.registry.before.record, paths.before, "receipt before path");
equal(receipt.static.registry.before.result_sha256, sha256(paths.before), "receipt before hash");
equal(correction.before.record, paths.before, "correction before path");
equal(correction.before.record_sha256, sha256(paths.before), "correction before hash");

console.log(JSON.stringify({status:"PASS", before:{revision:before.rev,scopes:110,bad:2}, after:{revision:after.rev,scopes:110,bad:0}, corrected_scope_keys:changedExpected, invalid_round2_artifact:{path:paths.invalid,sha256:sha256(paths.invalid),credit:"NONE"}}, null, 2));
