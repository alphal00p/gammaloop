// Follow-up scope-B probes: tuple-arity escapes with expression-bearing objects, correct
// namespaced coefficient traps. Same harness contract as probes.mjs (semantic entrypoint only).
import fs from "node:fs";
import path from "node:path";
import crypto from "node:crypto";
import { spawnSync } from "node:child_process";
const [evidence, semantic, out] = process.argv.slice(2);
const E = "xs_qed_eedd_L1", MM = "xs_qed_eemm_L1", PI = "xs_qed_eeee_1pi_L1";
const FOLDERS = ["historical", "current", "metadata-historical", "metadata-current"];
const pv = (j) => j.diagrams[0].payload;
const own = (m, o) => m.diagrams[0].owners.find((x) => x.owner === o);
const both = (c, f) => [["current", c, f], ["historical", c, f]];
const live = (c, o, f) => ["metadata-current", "metadata-historical"].map((m) => [m, c, (j) => { own(j, o).expression = f(own(j, o).expression); }]);
const P = 0x100000001b3n, M = (1n << 64n) - 1n, OFF = 0xcbf29ce484222325n;
const fb = (h, b) => { for (const x of b) h = ((h ^ BigInt(x)) * P) & M; return h; };
const fj = (h, v) => fb(h, Buffer.from(JSON.stringify(v)));
const fu = (h, n) => { const b = Buffer.alloc(8); b.writeBigUInt64LE(n); return fb(h, b); };
function fnvDir(dir) { const r = {}; for (const s of ["historical", "current"]) for (const c of [E, MM, PI]) {
  const j = JSON.parse(fs.readFileSync(path.join(dir, s, `${c}.json`)));
  const per = j.diagrams.map((d) => { const p = d.payload; let h = OFF;
    for (const f of ["numerator", "numerator_prefactor", "projector"]) h = fj(h, p[f]);
    for (const f of ["num_interface", "numerator_prefactor_interface", "projector_interface"]) if (Object.hasOwn(p, f)) h = fj(h, p[f]);
    for (const fr of [...p.vertices, ...p.edges.map((e) => e[1])]) { h = fj(h, fr.numerator); if (Object.hasOwn(fr, "num_interface")) h = fj(h, fr.num_interface); }
    return fu(fu(OFF, fj(OFF, d.key)), h); }).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  r[`${s}/${c}`] = per.reduce((h, v) => fu(h, v), OFF).toString(16).padStart(16, "0"); } return r; }
const manifest = Object.fromEntries(fs.readFileSync(path.join(evidence, "manifest.sha256"), "utf8").trim().split("\n").map((l) => l.split(/\s+/)).map(([h, f]) => [f, h]));
const shaOk = (dir) => FOLDERS.every((f) => [E, MM, PI].every((c) => crypto.createHash("sha256").update(fs.readFileSync(path.join(dir, f, `${c}.json`))).digest("hex") === manifest[`${f}/${c}.json`]));
const rogue = (j, k, expr) => { const e = pv(j).edges[k]; const d = structuredClone(e[1]); delete d.num_interface; d.numerator = expr; e.push(d); };
const GC = "UFO::{}::GC_1";
const C = [
  ["G1m-mirrored-edge2-third-element-garbage-numerator", "reject", "both payloads: edges[2].push(edge data w/o num_interface, numerator 'garbage(')", both(E, (j) => rogue(j, 2, "garbage("))],
  ["G1m2-mirrored-edge0-third-element-zero-numerator", "reject", "both payloads: edges[0].push(edge data, numerator '0')", both(E, (j) => rogue(j, 0, "0"))],
  ["G1m3-mirrored-eemm-edge4-third-element", "reject", "eemm both payloads: edges[4].push(edge data, numerator '-1')", both(MM, (j) => rogue(j, 4, "-1"))],
  ["G1h-historical-only-third-element", "reject", "historical only: edges[2].push(rogue)", [["historical", E, (j) => rogue(j, 2, "garbage(")]]],
  ["G3b-mirrored-edge-extra-incidence-eeee", "reject", "eeee both: edges[1].push({source,target})", both(PI, (j) => { pv(j).edges[1].push({ source: 9, target: 9 }); })],
  ["G19-mirrored-edge-short-tuple", "reject", "both: edges[0] = [edges[0][1]] (drop incidence)", both(E, (j) => { pv(j).edges[0] = [pv(j).edges[0][1]]; })],
  ["T2b-gc-split-power", "accept", "aggregate GC_1^2 -> GC_1*GC_1 everywhere (equivalent)", [...both(E, (j) => { pv(j).numerator = pv(j).numerator.replace(`${GC}^2`, `${GC}*${GC}`); }), ...live(E, "numerator", (s) => s.replace(`${GC}^2`, `${GC}*${GC}`))]],
  ["T2c-gc-into-sum", "reject", "aggregate GC_1^2 -> (GC_1+GC_1)*GC_1/2? use (GC_1+0)*GC_1 everywhere (sum atomic, conservative)", [...both(E, (j) => { pv(j).numerator = pv(j).numerator.replace(`${GC}^2`, `(${GC}+0)*${GC}`); }), ...live(E, "numerator", (s) => s.replace(`${GC}^2`, `(${GC}+0)*${GC}`))]],
  ["T4b-gc-power-zero", "reject", "aggregate * GC_1^0 everywhere (conservative)", [...both(E, (j) => { pv(j).numerator += `*${GC}^0`; }), ...live(E, "numerator", (s) => s + `*${GC}^0`)]],
  ["T6-gc-power-swap-local", "reject", "aggregate GC_1^2 -> GC_1^(1+1) (exponent sum, conservative)", [...both(E, (j) => { pv(j).numerator = pv(j).numerator.replace(`${GC}^2`, `${GC}^(1+1)`); }), ...live(E, "numerator", (s) => s.replace(`${GC}^2`, `${GC}^(1+1)`))]],
  ["T7-neg-pow-precedence-aggregate", "reject", "aggregate prefixed with -𝑖^2* (= +1 via -(𝑖^2)? no: -(−1)=+1) -> accept-equivalent check", [...both(E, (j) => { pv(j).numerator = "-𝑖^2*" + pv(j).numerator; }), ...live(E, "numerator", (s) => "-𝑖^2*" + s)]],
  ["T8-neg-pow-wrong-aggregate", "reject", "aggregate prefixed with (-𝑖)^2* (= -1)", [...both(E, (j) => { pv(j).numerator = "(-𝑖)^2*" + pv(j).numerator; }), ...live(E, "numerator", (s) => "(-𝑖)^2*" + s)]],
];
C[10][1] = "accept"; // -𝑖^2 = -(𝑖^2) = +1, an equivalent form
const baseF = JSON.stringify(fnvDir(evidence));
const results = [];
for (const [id, expect, desc, edits] of C) {
  const dir = path.join(out, id);
  fs.rmSync(dir, { recursive: true, force: true });
  for (const f of FOLDERS) fs.cpSync(path.join(evidence, f), path.join(dir, f), { recursive: true });
  let editError = null;
  try { for (const [folder, name, fn] of edits) { const f = path.join(dir, folder, `${name}.json`); const before = fs.readFileSync(f, "utf8"); const j = JSON.parse(before); fn(j); const after = JSON.stringify(j, null, 2); if (JSON.stringify(JSON.parse(before)) === JSON.stringify(j)) throw new Error(`no-op ${folder}/${name}`); fs.writeFileSync(f, after); } } catch (e) { editError = e.message; }
  const r = spawnSync("node", [semantic, dir], { encoding: "utf8" });
  let failures = []; try { failures = JSON.parse(r.stdout).failures; } catch { failures = ["<unparsed>"]; }
  let fnvChanged; try { fnvChanged = JSON.stringify(fnvDir(dir)) !== baseF; } catch (e) { fnvChanged = `fnv-error ${e.message.slice(0, 40)}`; }
  results.push({ id, expect, semantic: r.status === 0 ? "accept" : "reject", capture_sha_unchanged: shaOk(dir), fnv_changed: fnvChanged, editError, desc, failures: failures.slice(0, 4) });
}
fs.writeFileSync(path.join(out, "results.json"), JSON.stringify(results, null, 1) + "\n");
for (const x of results) console.log(`${x.id.padEnd(52)} exp=${x.expect.padEnd(7)} sem=${x.semantic} shaSame=${x.capture_sha_unchanged} fnvΔ=${x.fnv_changed} ${x.editError ? "EDIT-ERR " + x.editError : ""}\n      ${(x.failures[0] ?? "").slice(0, 170)}`);
