// Reviewer-written FNV-1a-64 recomputation (capture identity only, not semantic proof).
// Re-derived from the oracle digest definition: per diagram fnv(key json) then fnv(numerator fields), sorted, folded.
import fs from "node:fs";
const dir = process.argv[2];
const OFF = 14695981039346656037n, PR = 1099511628211n, MASK = (1n << 64n) - 1n;
const fnv = (h, bytes) => { for (const b of bytes) h = ((h ^ BigInt(b)) * PR) & MASK; return h; };
const js = (h, v) => fnv(h, Buffer.from(JSON.stringify(v), "utf8"));
const u64 = (h, n) => { const b = Buffer.alloc(8); b.writeBigUInt64LE(n); return fnv(h, b); };
const out = {};
for (const side of ["historical", "current"]) for (const c of ["xs_qed_eedd_L1", "xs_qed_eemm_L1", "xs_qed_eeee_1pi_L1"]) {
  const j = JSON.parse(fs.readFileSync(`${dir}/${side}/${c}.json`, "utf8"));
  const per = j.diagrams.map((d) => {
    const p = d.payload; let h = OFF;
    for (const f of ["numerator", "numerator_prefactor", "projector"]) h = js(h, p[f]);
    for (const f of ["num_interface", "numerator_prefactor_interface", "projector_interface"]) if (f in p) h = js(h, p[f]);
    for (const fr of [...p.vertices, ...p.edges.map((e) => e[1])]) { h = js(h, fr.numerator); if ("num_interface" in fr) h = js(h, fr.num_interface); }
    return u64(u64(OFF, js(OFF, d.key)), h);
  }).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  out[`${side}/${c}`] = per.reduce((h, v) => u64(h, v), OFF).toString(16).padStart(16, "0");
}
console.log(JSON.stringify(out));
