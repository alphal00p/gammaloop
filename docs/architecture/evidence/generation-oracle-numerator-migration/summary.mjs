// Independent re-implementation of census.rs Summary digests (content, numerators, named, ids)
// over write_exact_dump JSON. Usage: node summary.mjs <dump.json>...
import fs from "node:fs";
const M = (1n << 64n) - 1n, OFF = 0xcbf29ce484222325n, PR = 0x100000001b3n;
class D { constructor(v = OFF) { this.v = v; }
  bytes(b) { let v = this.v; for (const x of b) v = ((v ^ BigInt(x)) * PR) & M; return new D(v); }
  text(s) { return this.bytes(Buffer.from(s, "utf8")).bytes([0xff]); }
  name(n) { return n == null ? this.bytes([0xfe]) : this.text(n); }
  num(n) { const b = Buffer.alloc(8); b.writeBigUInt64LE(n); return this.bytes(b); }
  json(v) { return this.bytes(Buffer.from(JSON.stringify(v), "utf8")); } }
const ms = (xs) => xs.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0)).reduce((d, x) => d.num(x), new D()).v;
function digestJson(p) {
  p = structuredClone(p); let n = new D();
  for (const f of ["numerator", "numerator_prefactor", "projector"]) { n = n.json(p[f]); delete p[f]; }
  for (const f of ["num_interface", "numerator_prefactor_interface", "projector_interface"]) if (f in p) { n = n.json(p[f]); delete p[f]; }
  const frags = [];
  const take = (o) => { let i; if ("num_interface" in o) { i = o.num_interface; delete o.num_interface; } const t = o.numerator; o.numerator = null; return [t, i]; };
  for (const v of p.vertices) frags.push(take(v));
  for (const e of p.edges) frags.push(take(e[1]));
  for (const [f, i] of frags) { n = n.json(f); if (i !== undefined) n = n.json(i); }
  for (const f of ["model", "symmetry_factor", "overall_factor"]) delete p[f];
  return [new D().json(p).v, n.v];
}
const hex = (v) => "0x" + v.toString(16).padStart(16, "0");
for (const file of process.argv.slice(2)) {
  const d = JSON.parse(fs.readFileSync(file));
  const content = [], nums = [], ids = []; let named = new D();
  for (const g of d.diagrams) {
    const key = new D().json(g.key).v; const [frame, num] = digestJson(g.payload);
    content.push(new D().num(key).num(BigInt(g.aut)).text(g.overall_factor).num(frame).v);
    nums.push(new D().num(key).num(num).v);
    const id = BigInt("0x" + g.id); const b = Buffer.alloc(16); b.writeBigUInt64LE(id & M, 0); b.writeBigUInt64LE(id >> 64n, 8);
    ids.push(new D().bytes(b).v);
    named = named.name(g.name).num(key);
  }
  for (const gr of d.groups) { named = named.num(BigInt(gr.master)); for (const m of gr.members) named = named.name(m.name).text(m.ratio).text(m.overall_factor); }
  console.log(JSON.stringify({ file, diagrams: d.diagrams.length, content: hex(ms(content)), numerators: hex(ms(nums)), named: hex(named.v), ids: hex(ms(ids)) }));
}
