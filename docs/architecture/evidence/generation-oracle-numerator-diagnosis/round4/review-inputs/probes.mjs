// Scope-B round-3 reviewer probes. Copies the 4 consumed raw folders, edits parsed JSON,
// runs ONLY the given semantic entrypoint (no manifest/SHA/fixed FNV). Capture identity
// (SHA vs frozen manifest, reviewer FNV over numerator fields) is recorded separately.
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
const ALL4 = (c, f) => FOLDERS.map((x) => [x, c, f]);
// independent FNV-1a-64 of the oracle numerator digest (capture identity only)
const P = 0x100000001b3n, M = (1n << 64n) - 1n, OFF = 0xcbf29ce484222325n;
const fb = (h, b) => { for (const x of b) h = ((h ^ BigInt(x)) * P) & M; return h; };
const fj = (h, v) => fb(h, Buffer.from(JSON.stringify(v)));
const fu = (h, n) => { const b = Buffer.alloc(8); b.writeBigUInt64LE(n); return fb(h, b); };
function fnvDir(dir) {
  const r = {};
  for (const s of ["historical", "current"]) for (const c of [E, MM, PI]) {
    const j = JSON.parse(fs.readFileSync(path.join(dir, s, `${c}.json`)));
    const per = j.diagrams.map((d) => { const p = d.payload; let h = OFF;
      for (const f of ["numerator", "numerator_prefactor", "projector"]) h = fj(h, p[f]);
      for (const f of ["num_interface", "numerator_prefactor_interface", "projector_interface"]) if (Object.hasOwn(p, f)) h = fj(h, p[f]);
      for (const fr of [...p.vertices, ...p.edges.map((e) => e[1])]) { h = fj(h, fr.numerator); if (Object.hasOwn(fr, "num_interface")) h = fj(h, fr.num_interface); }
      return fu(fu(OFF, fj(OFF, d.key)), h); }).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
    r[`${s}/${c}`] = per.reduce((h, v) => fu(h, v), OFF).toString(16).padStart(16, "0");
  }
  return r;
}
const manifest = Object.fromEntries(fs.readFileSync(path.join(evidence, "manifest.sha256"), "utf8").trim().split("\n").map((l) => l.split(/\s+/)).map(([h, f]) => [f, h]));
const shaOk = (dir) => FOLDERS.every((f) => [E, MM, PI].every((c) => crypto.createHash("sha256").update(fs.readFileSync(path.join(dir, f, `${c}.json`))).digest("hex") === manifest[`${f}/${c}.json`]));
const deep = (v, f) => typeof v === "string" ? f(v) : Array.isArray(v) ? v.map((x) => deep(x, f)) : v && typeof v === "object" ? (Object.keys(v).forEach((k) => (v[k] = deep(v[k], f))), v) : v;
const swapKeys = (o) => Array.isArray(o) ? o.map(swapKeys) : o && typeof o === "object" ? Object.fromEntries(Object.keys(o).reverse().map((k) => [k, swapKeys(o[k])])) : o;

const C = [
  // ---- closure gaps: array/tuple arity & position ----
  ["G1-mirrored-edge-tuple-extra-numerator", "reject", "both payloads: edges[2] gains 3rd tuple element = clone of edge data with numerator 'garbage('",
    both(E, (j) => { const e = pv(j).edges[2]; e.push({ ...structuredClone(e[1]), numerator: "garbage(" }); })],
  ["G1b-onesided-edge-tuple-extra-numerator", "reject", "current only: same extra tuple element",
    [["current", E, (j) => { const e = pv(j).edges[2]; e.push({ ...structuredClone(e[1]), numerator: "garbage(" }); }]]],
  ["G2-mirrored-edge-tuple-extra-iface", "reject", "both payloads: edges[2] 3rd element with num_interface",
    both(E, (j) => { const e = pv(j).edges[2]; e.push(structuredClone(e[1])); })],
  ["G3-mirrored-edge-tuple-extra-incidence", "reject", "both payloads: edges[2] gains extra {source,target} element",
    both(E, (j) => { const e = pv(j).edges[2]; e.push({ source: 0, target: 1 }); })],
  ["G4-mirrored-edge-tuple-position-swap", "reject", "both payloads: edges[2] = [data, incidence]",
    both(E, (j) => { const e = pv(j).edges[2]; e.reverse(); })],
  ["G5-mirrored-scalar-to-array-aut", "capture", "both: diagram aut 1 -> [1] (type change, value)",
    both(E, (j) => { j.diagrams[0].aut = [j.diagrams[0].aut]; })],
  ["G6-mirrored-scalar-to-array-expression-string", "capture", "both: symmetry_factor -> ['garbage('] (expression-like string in scalar array)",
    both(E, (j) => { pv(j).symmetry_factor = ["garbage("]; })],
  ["G7-mirrored-number-to-string-source_slot", "capture", "both: edges[2][1].source_slot number -> string",
    both(E, (j) => { const d = pv(j).edges[2][1]; d.source_slot = String(d.source_slot); })],
  ["G8-mirrored-nesting-cut-element-array", "reject", "both: cuts[0].cut[0] wrapped in extra array",
    both(E, (j) => { const c = pv(j).cuts[0].cut; c[0] = [c[0]]; })],
  ["G9-mirrored-extra-group", "reject", "both: second group",
    both(E, (j) => { j.groups.push(structuredClone(j.groups[0])); })],
  ["G10-mirrored-missing-key-half_edge_order", "reject", "both: delete payload.half_edge_order",
    both(E, (j) => { delete pv(j).half_edge_order; })],
  ["G11-mirrored-missing-edge-optional-iface-both", "reject", "both payload+contract: delete edges[2][1].num_interface on current (optional signature admitted)",
    [["current", E, (j) => { delete pv(j).edges[2][1].num_interface; }]]],
  ["G12-historical-gains-edge-iface", "reject", "historical edges[2][1] gains num_interface (signature admitted by union)",
    [["historical", E, (j) => { pv(j).edges[2][1].num_interface = []; }]]],
  ["G13-mirrored-numerator-in-vertex-extra", "reject", "both payloads: extra vertex (no owner rows)",
    both(E, (j) => { pv(j).vertices.push(structuredClone(pv(j).vertices[0])); })],
  ["G14-metadata-owner-in-extra-diagram-field", "reject", "metadata-current: owner gets ports as string",
    [["metadata-current", E, (j) => { own(j, "vertex:0").ports = own(j, "vertex:0").ports.join(","); }]]],
  ["G15-numerator-type-array", "reject", "both payloads+live: vertex:0 numerator wrapped in array",
    both(E, (j) => { pv(j).vertices[0].numerator = [pv(j).vertices[0].numerator]; })],
  ["G16-iface-type-string", "reject", "current: vertex:0 num_interface as joined string",
    [["current", E, (j) => { pv(j).vertices[0].num_interface = pv(j).vertices[0].num_interface.join(","); }]]],
  ["G17-iface-null", "reject", "current: vertex:0 num_interface null",
    [["current", E, (j) => { pv(j).vertices[0].num_interface = null; }]]],
  ["G18-key-edges-extra-numerator-array", "reject", "both: diagram key.edges[0] gains extra element with numerator",
    both(E, (j) => { j.diagrams[0].key.edges.push({ diagram: { ...structuredClone(j.diagrams[0].key.edges[0].diagram) } }); })],
  // ---- legitimate equivalences must still accept ----
  ["V1-all-files-key-reversal", "accept", "every object in all 4 files of all cases with reversed key order",
    [E, MM, PI].flatMap((c) => FOLDERS.map((f) => [f, c, (j) => swapKeys(j)]))],
  ["V2-aggregate-factor-reorder-mirrored", "accept", "eedd aggregate numerator top-level factors reversed in both payloads and live",
    null],
  ["V3-current-iface-transpose-eedd-edge2", "accept", "current edge:2 serialized interface and live ports reversed consistently",
    [["current", E, (j) => { pv(j).edges[2][1].num_interface.reverse(); }], ["metadata-current", E, (j) => { own(j, "edge:2").ports.reverse(); }]]],
  ["V4-coef-equivalent-i-squared", "accept", "eedd edge:4 (-𝑖)* -> 𝑖*𝑖*𝑖* everywhere local",
    [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = e.numerator.replace("(-𝑖)*", "𝑖*𝑖*𝑖*"); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "𝑖*𝑖*𝑖*"))]],
  ["V5-coef-double-neg", "accept", "eedd edge:4 (-𝑖)* -> -(-(-𝑖))* everywhere local",
    [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = e.numerator.replace("(-𝑖)*", "-(-(-𝑖))*"); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "-(-(-𝑖))*"))]],
  ["V6-iface-whitespace-all", "accept", "current: spaces after every comma in all serialized interfaces of eemm",
    [["current", MM, (j) => deep(j, (s) => s.includes("::bis(") || s.includes("::mink(") || s.includes("::cof(") ? (s.startsWith("spenso::{spenso::") && !s.includes("*") ? s.replaceAll(",", ", ") : s) : s)]]],
  ["V7-owner-table-reverse-both", "accept", "metadata owners arrays reversed (both)",
    ["metadata-current", "metadata-historical"].map((m) => [m, PI, (j) => { j.diagrams[0].owners.reverse(); }])],
  // ---- coefficient semantics traps ----
  ["T1-coef-i-squared-wrong", "reject", "edge:4 (-𝑖)* -> 𝑖*𝑖* (=-1) everywhere local",
    [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = e.numerator.replace("(-𝑖)*", "𝑖*𝑖*"); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "𝑖*𝑖*"))]],
  ["T2-factor-moved-into-sum", "reject", "aggregate: GC_1^2 -> GC_1*(GC_1) unchanged? no: -> (GC_1+0)*GC_1",
    [...both(E, (j) => { pv(j).numerator = pv(j).numerator.replace("GC_1^2", "(GC_1+0)*GC_1"); }), ...live(E, "numerator", (s) => s.replace("GC_1^2", "(GC_1+0)*GC_1"))]],
  ["T3-div-coefficient", "reject", "aggregate gains *2/2 (rational not reduced, conservative)",
    [...both(E, (j) => { pv(j).numerator += "*2/2"; }), ...live(E, "numerator", (s) => s + "*2/2")]],
  ["T4-power-zero", "reject", "aggregate gains *GC_1^0 (conservative)",
    [...both(E, (j) => { pv(j).numerator += "*GC_1^0"; }), ...live(E, "numerator", (s) => s + "*GC_1^0")]],
  ["T5-local-factor-moved-between-owners", "accept?", "eedd: move (-𝑖) from edge:4 to edge:3 (local redistribution; aggregate product preserved)",
    null],
];
// dynamic edits needing inspection of data
function editsFor(id) {
  if (id === "V2-aggregate-factor-reorder-mirrored") {
    const rev = (s) => { // split top-level '*' and reverse
      const parts = []; let d = 0, cur = "";
      for (const ch of s) { if (ch === "(") d++; if (ch === ")") d--; if (ch === "*" && d === 0) { parts.push(cur); cur = ""; } else cur += ch; }
      parts.push(cur); return parts.reverse().join("*");
    };
    return [...both(E, (j) => { pv(j).numerator = rev(pv(j).numerator); }), ...live(E, "numerator", rev)];
  }
  if (id === "T5-local-factor-moved-between-owners") {
    return [...both(E, (j) => { const a = pv(j).edges[4][1], b = pv(j).edges[3][1]; a.numerator = a.numerator.replace("(-𝑖)*", ""); b.numerator = "(-𝑖)*" + b.numerator; }),
      ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "")), ...live(E, "edge:3", (s) => "(-𝑖)*" + s)];
  }
}
const baseF = JSON.stringify(fnvDir(evidence));
const results = [];
for (let [id, expect, desc, edits] of C) {
  edits ??= editsFor(id);
  const dir = path.join(out, id);
  fs.rmSync(dir, { recursive: true, force: true });
  for (const f of FOLDERS) fs.cpSync(path.join(evidence, f), path.join(dir, f), { recursive: true });
  let editError = null;
  try {
    for (const [folder, name, fn] of edits) {
      const f = path.join(dir, folder, `${name}.json`);
      const before = fs.readFileSync(f, "utf8"); const j = JSON.parse(before); const r = fn(j); const v = r && typeof r === "object" && !Array.isArray(r) && (r.diagrams || r.case) ? r : j;
      const after = JSON.stringify(v, null, 2);
      if (after === before || (JSON.stringify(JSON.parse(before)) === JSON.stringify(v) && id[0] !== "V")) throw new Error(`no-op edit ${folder}/${name}`);
      fs.writeFileSync(f, after);
    }
  } catch (e) { editError = e.message; }
  const r = spawnSync("node", [semantic, dir], { encoding: "utf8" });
  let failures = []; try { failures = JSON.parse(r.stdout).failures; } catch { failures = ["<unparsed>"]; }
  let fnvChanged = null; try { fnvChanged = JSON.stringify(fnvDir(dir)) !== baseF; } catch (e) { fnvChanged = `fnv-error ${e.message}`; }
  results.push({ id, expect, semantic: r.status === 0 ? "accept" : "reject", capture_sha_unchanged: shaOk(dir), fnv_changed: fnvChanged, editError, desc, failures: failures.slice(0, 4) });
}
fs.writeFileSync(path.join(out, "results.json"), JSON.stringify(results, null, 1) + "\n");
for (const x of results) console.log(`${x.id.padEnd(48)} exp=${x.expect.padEnd(8)} sem=${x.semantic} shaSame=${x.capture_sha_unchanged} fnvΔ=${x.fnv_changed} ${x.editError ? "EDIT-ERR " + x.editError : ""}\n      ${(x.failures[0] ?? "").slice(0, 170)}`);
