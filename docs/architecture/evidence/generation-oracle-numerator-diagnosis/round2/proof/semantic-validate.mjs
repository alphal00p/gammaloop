// Independent scope-B checker: factorized typed-interface reconstruction without expansion.
// Usage: node independent.mjs <evidence-dir>   (exit 0 = all checks pass)
import fs from "node:fs";
import path from "node:path";
import crypto from "node:crypto";
import { fileURLToPath } from "node:url";

const root = process.argv[2];
if (!root) throw new Error("usage: node semantic-validate.mjs <evidence-dir>");
const here = path.dirname(fileURLToPath(import.meta.url));
const CONTRACT = JSON.parse(fs.readFileSync(path.join(here, "semantic-contract.json"), "utf8"));
const CASES = ["xs_qed_eedd_L1", "xs_qed_eemm_L1", "xs_qed_eeee_1pi_L1"];
const failures = [];
const note = (c, m) => failures.push(`${c}: ${m}`);
const read = (f, n) => JSON.parse(fs.readFileSync(path.join(root, f, `${n}.json`), "utf8"));

// Closed constructor map, from authoritative symbol!(...) definitions:
// historical crates/feynkit-graph/src/symbols.rs tags vs current Notation registrations.
const HIST = {
  Q: "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q",
  K: "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:k}::K",
  P: "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:p}::P",
  edge: "gammalooprs::{spenso::index,spenso::index-label:e}::edge",
  vertex: "gammalooprs::{spenso::index,spenso::index-label:v}::vertex",
};
const CUR = {
  Q: "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q",
  K: "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::K",
  P: "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::P",
  edge: "gammalooprs::{spenso::index}::edge",
  vertex: "gammalooprs::{spenso::index}::vertex",
};
for (const key of Object.keys(HIST)) {
  if (CONTRACT.authority.constructor_map[key][0] !== HIST[key]
      || CONTRACT.authority.constructor_map[key][1] !== CUR[key]) {
    throw new Error(`constructor contract mismatch for ${key}`);
  }
}
const digest = (value) => crypto.createHash("sha256").update(JSON.stringify(value)).digest("hex");

// ---------- parser (symbol-level, so the map acts on whole heads, never substrings) ----------
function lex(s) {
  const t = [];
  let i = 0;
  while (i < s.length) {
    const c = s[i];
    if (c === " ") { i++; continue; }
    if ("()+-*^,/".includes(c)) { t.push({ k: c }); i++; continue; }
    if (s.startsWith("𝑖", i)) { t.push({ k: "i" }); i += "𝑖".length; continue; }
    if (/[0-9]/.test(c)) { let j = i; while (j < s.length && /[0-9]/.test(s[j])) j++; t.push({ k: "num", v: BigInt(s.slice(i, j)) }); i = j; continue; }
    if (/[A-Za-z_]/.test(c)) {
      let j = i;
      while (j < s.length) {
        if (/[A-Za-z0-9_:]/.test(s[j])) { j++; continue; }
        if (s[j] === "{") { let d = 0; do { if (s[j] === "{") d++; if (s[j] === "}") d--; j++; } while (d > 0 && j < s.length); continue; }
        break;
      }
      t.push({ k: "name", v: s.slice(i, j) }); i = j; continue;
    }
    throw new Error(`lex: unexpected ${JSON.stringify(c)} at ${i} in ${s.slice(0, 80)}`);
  }
  return t;
}
function parse(s) {
  const t = lex(s); let p = 0;
  const peek = () => t[p]?.k, eat = (k) => { if (t[p]?.k !== k) throw new Error(`parse: want ${k} got ${t[p]?.k} @${p}`); return t[p++]; };
  function sum() { const terms = [prod()]; while (peek() === "+" || peek() === "-") { const op = eat(peek()).k; const f = prod(); terms.push(op === "-" ? { k: "neg", a: f } : f); } return terms.length === 1 ? terms[0] : { k: "sum", terms }; }
  function prod() { const fs = [unary()]; while (peek() === "*" || peek() === "/") { const op = eat(peek()).k; const f = unary(); fs.push(op === "/" ? { k: "inv", a: f } : f); } return fs.length === 1 ? fs[0] : { k: "prod", fs }; }
  function unary() { if (peek() === "-") { eat("-"); return { k: "neg", a: unary() }; } return pow(); }
  function pow() { const b = atom(); if (peek() === "^") { eat("^"); return { k: "pow", b, e: unary() }; } return b; }
  function atom() {
    if (peek() === "(") { eat("("); const e = sum(); eat(")"); return e; }
    if (peek() === "num") return { k: "num", v: eat("num").v };
    if (peek() === "i") { eat("i"); return { k: "i" }; }
    const n = eat("name").v;
    if (peek() === "(") { eat("("); const args = [sum()]; while (peek() === ",") { eat(","); args.push(sum()); } eat(")"); return { k: "call", h: n, args }; }
    return { k: "sym", h: n };
  }
  const e = sum(); if (p !== t.length) throw new Error("parse: trailing tokens"); return e;
}
const show = (e) => {
  switch (e.k) {
    case "num": return e.v.toString(); case "i": return "𝑖"; case "sym": return e.h;
    case "call": return `${e.h}(${e.args.map(show).join(",")})`;
    case "neg": return `-(${show(e.a)})`; case "inv": return `/(${show(e.a)})`;
    case "pow": return `(${show(e.b)})^(${show(e.e)})`;
    case "prod": return e.fs.map((f) => `(${show(f)})`).join("*");
    case "sum": return e.terms.map((f) => `(${show(f)})`).join("+");
  }
};
function mapHeads(e, from, to) {
  const m = new Map(Object.keys(from).map((k) => [from[k], to[k]]));
  const r = (x) => {
    if (x.k === "call" || x.k === "sym") { const h = m.get(x.h) ?? x.h; return x.k === "call" ? { ...x, h, args: x.args.map(r) } : { ...x, h }; }
    if (x.k === "sum") return { ...x, terms: x.terms.map(r) };
    if (x.k === "prod") return { ...x, fs: x.fs.map(r) };
    if (x.k === "pow") return { ...x, b: r(x.b), e: r(x.e) };
    if (x.k === "neg" || x.k === "inv") return { ...x, a: r(x.a) };
    return x;
  };
  return r(e);
}
function heads(e, out = new Set()) {
  if (e.k === "call" || e.k === "sym") out.add(e.h);
  for (const c of e.args ?? e.terms ?? e.fs ?? [e.a, e.b, e.e].filter(Boolean)) heads(c, out);
  return out;
}

// ---------- typed slots and free-index algebra ----------
const REP = { "spenso::{spenso::representation,spenso::self_dual}::bis": "sd", "spenso::{spenso::representation,spenso::self_dual}::mink": "sd", "spenso::{spenso::representation,spenso::dualizable}::cof": "dual" };
const LOWER = "spenso::{spenso::lower}::dind";
function slotOf(e) {
  let lower = false, s = e;
  if (s.k === "call" && s.h === LOWER) { if (s.args.length !== 1) throw new Error("dind arity"); lower = true; s = s.args[0]; }
  if (s.k !== "call" || !(s.h in REP)) { if (lower) throw new Error("dind of non-slot"); return null; }
  if (s.args.length !== 2 || s.args[0].k !== "num") throw new Error(`slot arity/dimension ${show(s)}`);
  if (lower && REP[s.h] !== "dual") throw new Error("lowered self-dual slot");
  return { key: `${s.h}|${s.args[0].v}|${show(s.args[1])}`, rep: REP[s.h], up: !lower, text: show(e), base: show(s) };
}
// returns Map key -> {rep, ups, downs, texts}
function free(e) {
  switch (e.k) {
    case "num": case "i": case "sym": return new Map();
    case "neg": case "inv": return free(e.a);
    case "pow": { const f = free(e.b); if (f.size) throw new Error("power of indexed factor"); if (free(e.e).size) throw new Error("indexed exponent"); return f; }
    case "call": {
      const sl = slotOf(e); if (sl) return addSlot(new Map(), sl);
      let acc = new Map();
      for (const a of e.args) acc = combine(acc, free(a));
      return acc;
    }
    case "prod": return e.fs.map(free).reduce(combine, new Map());
    case "sum": {
      const fs = e.terms.map(free), ref = sig(fs[0]);
      for (const f of fs) if (sig(f) !== ref) throw new Error(`sum terms with unequal free interfaces: ${sig(f)} vs ${ref}`);
      return fs[0];
    }
  }
}
function addSlot(m, s) { const x = m.get(s.key) ?? { rep: s.rep, ups: 0, downs: 0, base: s.base }; s.up ? x.ups++ : x.downs++; m.set(s.key, x); return m; }
function combine(a, b) {
  const m = new Map([...a].map(([k, v]) => [k, { ...v }]));
  for (const [k, v] of b) { const x = m.get(k) ?? { rep: v.rep, ups: 0, downs: 0, base: v.base }; x.ups += v.ups; x.downs += v.downs; m.set(k, x); }
  for (const [k, x] of m) {
    if (x.rep === "sd") { const n = x.ups + x.downs; if (n === 2) m.delete(k); else if (n > 2) throw new Error(`self-dual index used ${n} times: ${k}`); }
    else { if (x.ups === 1 && x.downs === 1) m.delete(k); else if (x.ups > 1 || x.downs > 1) throw new Error(`dualizable variance clash ${k} up=${x.ups} down=${x.downs}`); }
  }
  return m;
}
const portsOf = (m) => [...m.values()].map((x) => (x.rep === "dual" && x.downs ? `${LOWER}(${x.base})` : x.base)).sort();
const sig = (m) => JSON.stringify(portsOf(m));
const canonPort = (s) => { const sl = slotOf(parse(s)); if (!sl) throw new Error(`not a slot: ${s}`); return sl.text; };

// ---------- factor multiset (no distribution: sums stay atomic) ----------
function factors(e, acc = { coef: [1n, 0n], fs: [] }) {
  const mul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
  switch (e.k) {
    case "num": acc.coef = mul(acc.coef, [e.v, 0n]); break;
    case "i": acc.coef = mul(acc.coef, [0n, 1n]); break;
    case "neg": acc.coef = mul(acc.coef, [-1n, 0n]); factors(e.a, acc); break;
    case "prod": for (const f of e.fs) factors(f, acc); break;
    case "pow": if (e.e.k === "num" && e.e.v > 0n) { for (let n = 0n; n < e.e.v; n++) factors(e.b, acc); } else acc.fs.push(show(e)); break;
    default: acc.fs.push(show(e));
  }
  return acc;
}

// ---------- historical cook decoding by structural alignment ----------
function alignCooked(cooked, plain, map, c, owner) {
  // Cooking reorders symmetric arguments and commutative factors/terms by cooked
  // symbol order, so alignment searches those orderings with backtracking.
  const r = walk(cooked, plain, map);
  if (!r) throw new Error(`${owner}: historical live expression does not decode to payload fragment`);
  for (const [k, v] of r) map.set(k, v);
  function kids(x) { return x.args ?? x.terms ?? x.fs ?? [x.a, x.b, x.e].filter(Boolean); }
  function walk(a, b, m) {
    const ck = a.k === "sym" ? /^idenso::\{([^}]*)\}::cooked_[0-9a-f]+$/.exec(a.h) : null;
    if (ck) {
      if (b.k !== "call") return null;
      const hb = /^gammalooprs::\{([^}]*)\}::(hedge|edge|vertex)$/.exec(b.h);
      if (!hb || hb[1] !== ck[1]) return null;
      const v = show(b), prev = m.get(a.h);
      if (prev && prev !== v) return null;
      for (const [k, w] of m) if (w === v && k !== a.h) return null;
      const n = new Map(m); n.set(a.h, v); return n;
    }
    if (a.k !== b.k) return null;
    if ((a.k === "call" || a.k === "sym") && a.h !== b.h) return null;
    if (a.k === "num" && a.v !== b.v) return null;
    const ca = kids(a), cb = kids(b);
    if (ca.length !== cb.length) return null;
    const unordered = a.k === "prod" || a.k === "sum" || (a.k === "call" && /^[^:]*::\{symmetric[,}]/.test(a.h));
    if (!unordered) { let cur = m; for (let i = 0; i < ca.length; i++) { cur = walk(ca[i], cb[i], cur); if (!cur) return null; } return cur; }
    const used = new Array(cb.length).fill(false);
    const go = (i, cur) => { if (i === ca.length) return cur; for (let j = 0; j < cb.length; j++) { if (used[j]) continue; const n = walk(ca[i], cb[j], cur); if (n) { used[j] = true; const r = go(i + 1, n); if (r) return r; used[j] = false; } } return null; };
    return go(0, m);
  }
}

const report = [];
for (const name of CASES) {
  const H = read("historical", name), C = read("current", name);
  const MH = read("metadata-historical", name), MC = read("metadata-current", name);
  const expected = CONTRACT.cases[name];
  try {
    // 1. Constructor and symbol names are parsed as complete heads and must belong
    // to the closed inventory frozen from the authoritative producer definitions.
    const allowedNames = {
      historical: new Set(expected.payload_names.historical),
      current: new Set(expected.payload_names.current),
    };
    // 2. Owner tables: exact same owner list (no duplicates, missing or extra).
    const owners = (M) => M.diagrams.map((d) => d.owners.map((o) => o.owner));
    if (MH.diagrams.length !== 1 || MC.diagrams.length !== 1 || H.diagrams.length !== 1 || C.diagrams.length !== 1) note(name, "diagram count");
    const oh = owners(MH)[0], oc = owners(MC)[0];
    if (new Set(oh).size !== oh.length || new Set(oc).size !== oc.length) note(name, "duplicate owner");
    const P = C.diagrams[0].payload, PH = H.diagrams[0].payload;
    const expectedOwners = [...P.vertices.map((_, i) => `vertex:${i}`), ...P.edges.map((_, i) => `edge:${i}`), "numerator", "prefactor", "projector"];
    for (const [lab, o] of [["historical", oh], ["current", oc]]) if (JSON.stringify([...o].sort()) !== JSON.stringify([...expectedOwners].sort())) note(name, `${lab} owner table ${JSON.stringify(o)} != payload owners`);
    if (MH.diagrams[0].id !== C.diagrams[0].id || MC.diagrams[0].id !== C.diagrams[0].id || H.diagrams[0].id !== C.diagrams[0].id) note(name, "diagram id binding");
    const fragment = (Pl, owner) => {
      const [k, i] = owner.split(":");
      if (k === "vertex") return { expr: Pl.vertices[+i].numerator, iface: Pl.vertices[+i].num_interface, has: "num_interface" in Pl.vertices[+i] };
      if (k === "edge") return { expr: Pl.edges[+i][1].numerator, iface: Pl.edges[+i][1].num_interface, has: "num_interface" in Pl.edges[+i][1] };
      const f = { numerator: ["numerator", "num_interface"], prefactor: ["numerator_prefactor", "numerator_prefactor_interface"], projector: ["projector", "projector_interface"] }[owner];
      return { expr: Pl[f[0]], iface: Pl[f[1]], has: f[1] in Pl };
    };
    // 3. Every interface-like key must be at one exact owner location. Presence,
    // absence and ordered contents are checked before any comparison omits it.
    const currentInterfaceKeys = [], historicalInterfaceKeys = [];
    const scan = (o, p, out) => { if (o && typeof o === "object") for (const [k, v] of Object.entries(o)) { if (k.endsWith("_interface")) out.push(`${p}/${k}`); scan(v, `${p}/${k}`, out); } };
    scan(C, "", currentInterfaceKeys); scan(H, "", historicalInterfaceKeys);
    const allowed = new Set([...P.vertices.map((_, i) => `/diagrams/0/payload/vertices/${i}/num_interface`), ...P.edges.map((_, i) => `/diagrams/0/payload/edges/${i}/1/num_interface`), ...["num_interface", "numerator_prefactor_interface", "projector_interface"].map((f) => `/diagrams/0/payload/${f}`)]);
    for (const key of currentInterfaceKeys) if (!allowed.has(key)) note(name, `unexpected interface field ${key}`);
    if (historicalInterfaceKeys.length) note(name, `historical carries interface fields ${historicalInterfaceKeys.join(",")}`);

    const mOwnH = new Map(MH.diagrams[0].owners.map((o) => [o.owner, o])), mOwnC = new Map(MC.diagrams[0].owners.map((o) => [o.owner, o]));
    const cookMap = new Map();
    let axes = 0; const perms = []; const local = [];
    const half = P.half_edge_order;
    const vertexOfHalf = (n) => { const h = half[n]; return P.edges[h.edge][0][h.endpoint]; };
    for (const owner of expectedOwners) {
      const fc = fragment(P, owner), fh = fragment(PH, owner), oH = mOwnH.get(owner), oC = mOwnC.get(owner);
      if (!oH || !oC) { note(name, `owner ${owner} missing`); continue; }
      const ec = parse(fc.expr), eh = parse(fh.expr);
      for (const [side, expression] of [["current", ec], ["historical", eh]]) {
        for (const head of heads(expression)) {
          if (!allowedNames[side].has(head)) note(name, `${owner}: unsupported ${side} symbol head ${head}`);
        }
      }
      // 4. Constructor map at symbol level.
      if (show(mapHeads(ec, CUR, HIST)) !== show(eh)) note(name, `${owner}: payload expression differs beyond constructor map`);
      // 5. Payload <-> live owner expression binding (current direct; historical via cook decoding).
      if (show(parse(oC.expression)) !== show(ec)) note(name, `${owner}: current live expression != payload fragment`);
      alignCooked(parse(oH.expression), eh, cookMap, name, owner);
      // 6. Effective free interface from the factorized expression.
      const fr = free(ec), frh = free(eh);
      if (sig(fr) !== sig(frh)) note(name, `${owner}: free interfaces differ`);
      const pc = oC.ports.map(canonPort), ph = oH.ports.map(canonPort);
      if (JSON.stringify([...pc].sort()) !== sig(fr)) note(name, `${owner}: current live ports ${JSON.stringify(pc)} != free ${sig(fr)}`);
      if (JSON.stringify([...ph].sort()) !== sig(frh)) note(name, `${owner}: historical live ports != free ${sig(frh)}`);
      if (new Set(pc).size !== pc.length || new Set(ph).size !== ph.length) note(name, `${owner}: repeated port`);
      // 7. Serialized interfaces are admitted only at known owners and must equal
      // the ordered live ports. Missing interfaces use the exact producer-observed
      // inferred order frozen in the semantic contract; no sorting theorem is assumed.
      const ownerExpected = expected.owners[owner];
      if (fc.has !== ownerExpected.current_interface_present) note(name, `${owner}: interface presence differs from producer contract`);
      if (fc.has) {
        if (!Array.isArray(fc.iface)) note(name, `${owner}: malformed serialized interface`);
        else if (JSON.stringify(fc.iface.map(canonPort)) !== JSON.stringify(pc)) note(name, `${owner}: serialized interface != live ordered ports`);
      } else if (JSON.stringify(pc) !== JSON.stringify(ownerExpected.current_ports)) {
        note(name, `${owner}: missing interface does not use observed inferred order`);
      }
      if (fh.has) note(name, `${owner}: historical interface present`);
      if (JSON.stringify(ph) !== JSON.stringify(ownerExpected.historical_ports)) note(name, `${owner}: historical live order differs from observed producer order`);
      // 8. Bijective transport.
      const perm = ph.map((p) => pc.indexOf(p));
      if (perm.some((x) => x < 0) || new Set(perm).size !== perm.length || perm.length !== pc.length) note(name, `${owner}: transport not bijective`);
      if (!perm.every((x, i) => x === i)) perms.push({ owner, perm });
      axes += pc.length;
      if (JSON.stringify([oH.metric, oH.composite]) !== JSON.stringify(ownerExpected.historical_flags)
          || JSON.stringify([oC.metric, oC.composite]) !== JSON.stringify(ownerExpected.current_flags)) {
        note(name, `${owner}: producer-observed metric/composite flags changed`);
      }
      // 9. Graph owner binding of hedge ports.
      for (const p of pc) {
        const m = /hedge\((\d+),1\)\)*$/.exec(p); if (!m) { note(name, `${owner}: non-hedge port ${p}`); continue; }
        const n = +m[1], [k, i] = owner.split(":");
        if (k === "edge" && half[n].edge !== +i) note(name, `${owner}: port hedge ${n} not on edge`);
        if (k === "vertex" && vertexOfHalf(n) !== +i) note(name, `${owner}: port hedge ${n} not at vertex`);
        if (!["edge", "vertex"].includes(k) && !P.edges[half[n].edge][1].external) note(name, `${owner}: aggregate port ${n} not external`);
      }
      if (owner.startsWith("vertex") || owner.startsWith("edge")) local.push(ec);
    }
    // 10. Cook decoding injective.
    const decoded = [...cookMap.values()];
    if (new Set(decoded).size !== decoded.length) note(name, "cook decoding not injective");
    // 11. Aggregate chain: numerator == product of all local fragments (factor multiset, sums atomic).
    const agg = factors(parse(P.numerator)), prodLocal = local.reduce((acc, e) => factors(e, acc), { coef: [1n, 0n], fs: [] });
    if (JSON.stringify([...agg.fs].sort()) !== JSON.stringify([...prodLocal.fs].sort())) note(name, "aggregate numerator factors != product of local factors");
    if (agg.coef.join() !== prodLocal.coef.join()) note(name, `aggregate coefficient ${agg.coef} != local product ${prodLocal.coef}`);
    // 12. Full contraction of prefactor*numerator*projector.
    const total = combine(combine(free(parse(P.numerator_prefactor)), free(parse(P.numerator))), free(parse(P.projector)));
    if (total.size) note(name, `uncontracted indices ${sig(total)}`);
    // 13. Each internal hedge index is shared by exactly its vertex and edge owner.
    // 14. Non-numerator content: remove only the exact interface fields already
    // checked above, never arbitrary recursively named fields.
    const stripKnownInterfaces = (data) => {
      const copy = structuredClone(data), payload = copy.diagrams[0].payload;
      delete payload.num_interface;
      delete payload.numerator_prefactor_interface;
      delete payload.projector_interface;
      for (const vertex of payload.vertices) delete vertex.num_interface;
      for (const edge of payload.edges) delete edge[1].num_interface;
      return copy;
    };
    const exprKeys = new Set(["numerator", "numerator_prefactor", "projector"]);
    const noExpr = (o) => Array.isArray(o) ? o.map(noExpr) : o && typeof o === "object" ? Object.fromEntries(Object.entries(o).filter(([k]) => !exprKeys.has(k)).map(([k, v]) => [k, noExpr(v)])) : o;
    const canon = (o) => Array.isArray(o) ? o.map(canon) : o && typeof o === "object" ? Object.fromEntries(Object.keys(o).sort().map((k) => [k, canon(o[k])])) : o;
    if (JSON.stringify(canon(noExpr(stripKnownInterfaces(H)))) !== JSON.stringify(canon(noExpr(stripKnownInterfaces(C))))) note(name, "non-numerator dump content differs (exact, no text mapping)");
    // 15. Group / factor chain.
    for (const [side, data] of [["historical", H], ["current", C]]) {
      const groups = data.groups, diagram = data.diagrams[0];
      if (groups.length !== 1 || groups[0].members.length !== 1 || groups[0].master !== 0
          || groups[0].members[0].source_index !== 0 || groups[0].members[0].ratio !== "1") note(name, `${side} group shape`);
      const member = groups[0].members[0];
      try { parse(member.ratio); parse(member.overall_factor); parse(diagram.overall_factor); }
      catch (error) { note(name, `${side} unsupported scalar syntax: ${error.message}`); }
      if (member.id !== diagram.id || member.name !== diagram.name
          || member.overall_factor !== diagram.overall_factor
          || diagram.payload.overall_factor !== diagram.overall_factor) note(name, `${side} group/member/source-factor binding`);
    }
    if (P.numerator_prefactor !== "1" || PH.numerator_prefactor !== "1"
        || MC.numerators.prefactor !== P.numerator_prefactor
        || MH.numerators.prefactor !== PH.numerator_prefactor) note(name, "prefactor binding");
    for (const k of ["model", "numerators", "grouping", "report", "case"]) if (JSON.stringify(MH[k]) !== JSON.stringify(MC[k])) note(name, `metadata ${k} differs`);
    if (MC.case !== name) note(name, "metadata case name");
    if (MC.diagrams[0].id !== expected.diagram_id || MH.diagrams[0].id !== expected.diagram_id) note(name, "semantic-contract diagram id");
    if (digest(MC.model) !== expected.model_sha256 || digest(MH.model) !== expected.model_sha256) note(name, "semantic-contract model");
    if (digest(MC.numerators) !== expected.options_sha256 || digest(MH.numerators) !== expected.options_sha256) note(name, "semantic-contract options");
    if (digest(MC.grouping) !== expected.grouping_sha256 || digest(MH.grouping) !== expected.grouping_sha256) note(name, "semantic-contract grouping");
    if (MC.report !== expected.report || MH.report !== expected.report) note(name, "semantic-contract report");
    if (cookMap.size !== 14) note(name, `historical cooked-symbol map has ${cookMap.size}, expected 14`);
    report.push({ name, axes, perms: perms.map((p) => `${p.owner}:${JSON.stringify(p.perm)}`), cooked: cookMap.size, aggFactors: agg.fs.length, coef: agg.coef.join(",") });
  } catch (e) { note(name, `exception ${e.message}`); }
}
console.log(JSON.stringify({ report, failures }, null, 1));
process.exit(failures.length ? 1 : 0);
