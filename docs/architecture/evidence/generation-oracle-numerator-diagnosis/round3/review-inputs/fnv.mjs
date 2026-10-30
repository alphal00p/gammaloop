// Independent FNV-1a 64 recomputation of the six numerator digests (capture identity only).
import fs from "node:fs";
const dir = process.argv[2];
const P = 0x100000001b3n, M = 0xffffffffffffffffn;
class D { constructor() { this.h = 0xcbf29ce484222325n; } bytes(b) { for (const x of b) this.h = ((this.h ^ BigInt(x)) * P) & M; return this; }
  json(v) { return this.bytes(new TextEncoder().encode(JSON.stringify(v))); } num(n) { const b = new Uint8Array(8); new DataView(b.buffer).setBigUint64(0, n, true); return this.bytes(b); } }
for (const side of ["historical", "current"]) for (const c of ["xs_qed_eedd_L1", "xs_qed_eemm_L1", "xs_qed_eeee_1pi_L1"]) {
  const j = JSON.parse(fs.readFileSync(`${dir}/${side}/${c}.json`, "utf8"));
  const vals = j.diagrams.map((d) => { const p = d.payload, n = new D();
    for (const f of ["numerator", "numerator_prefactor", "projector"]) n.json(p[f]);
    for (const f of ["num_interface", "numerator_prefactor_interface", "projector_interface"]) if (Object.hasOwn(p, f)) n.json(p[f]);
    for (const fr of [...p.vertices, ...p.edges.map((e) => e[1])]) { n.json(fr.numerator); if (Object.hasOwn(fr, "num_interface")) n.json(fr.num_interface); }
    return new D().num(new D().json(d.key).h).num(n.h).h; }).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  const m = vals.reduce((d, v) => d.num(v), new D()).h;
  console.log(side, c, m.toString(16).padStart(16, "0"));
}
