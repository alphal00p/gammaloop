#!/usr/bin/env node
// Contract validator for site/data (see DATA_FORMAT.md). No dependency beyond Node.
//
//   node tools/tests/check_site.mjs [SITE_DIR] [--strict] [--no-mock] [--max N]
//
// SITE_DIR defaults to ../../site relative to this file. Errors violate the contract and make the exit status 1;
// warnings point at things that are allowed but probably unintended (with --strict they fail too). --no-mock fails if
// any file still carries the marker "mock": true. --max limits the number of messages printed per category.

import { readFileSync, existsSync, statSync } from "node:fs";
import { resolve, dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const argv = process.argv.slice(2);
const opts = { strict: false, noMock: false, max: 200 };
let siteArg = null;
for (let i = 0; i < argv.length; i++) {
  const a = argv[i];
  if (a === "--strict") opts.strict = true;
  else if (a === "--no-mock") opts.noMock = true;
  else if (a === "--max") opts.max = Number(argv[++i]) || 200;
  else if (a === "-h" || a === "--help") {
    console.log("usage: node check_site.mjs [SITE_DIR] [--strict] [--no-mock] [--max N]");
    process.exit(0);
  } else siteArg = a;
}
const here = dirname(fileURLToPath(import.meta.url));
const SITE = resolve(siteArg ?? join(here, "..", "..", "site"));

const errors = [];
const warnings = [];
const mockFiles = [];
const err = (where, msg) => errors.push(`${where}: ${msg}`);
const warn = (where, msg) => warnings.push(`${where}: ${msg}`);

// ---------------------------------------------------------------------------------------------- constants
const TYPES = ["tree", "pole", "evanescent", "finite", "breaking"];
const KINDS = ["gauge", "fermion", "sterile", "scalar", "ghost", "antighost", "auxiliary", "source"];
const GROUPS = ["SU(3)", "SU(2)", "U(1)"];
const TERM_ID = /^[A-Za-z0-9_-]+$/;
const CLASS_ID = /^[A-Za-z0-9_.:,-]+$/;
const CARD_KEYS = [
  "gauge-kinetic", "gauge-self", "gauge-fixing", "ghost", "fermion-kinetic", "fermion-gauge", "yukawa", "higgs-kinetic",
  "higgs-potential", "brst-sources", "breaking", "ct-fields", "ct-couplings", "ct-evanescent", "st-gauge", "st-fermion",
  "st-source", "st-scalar", "normalisation", "anomaly", "broken-shift", "broken-masses", "broken-mixing", "broken-tadpole",
  "broken-st", "broken-interactions",
];
// method cards shown in the Theory view (gdocs.theory.METHOD_KEYS); terms do not refer to them
const METHOD_KEYS = ["method-overview", "method-bmhv", "method-st", "method-beyond"];
const RICH_TAGS = new Set(["p", "em", "strong", "a", "ul", "ol", "li", "code", "sub", "sup", "span", "br", "math", "i", "b"]);
const PRIVATE = [
  [/(?:^|[\s"'(=:])\/(?:home|shared|common|Users|nix\/store|tmp|var\/folders)\/[^\s"'<>)]*/, "absolute local path"],
  [/file:\/\//, "file:// URL"],
  [/[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+\.[A-Za-z]{2,}(?:\.[A-Za-z]{2,})*/, "e-mail address"],
  [/S-\d{3}-\d{4}\.\d{2}\.\d{2}-[A-Z0-9]{16,}/, "licence key"],
  [/SYMBOLICA_LICENSE\s*=\s*\S+/, "licence assignment"],
];

// ---------------------------------------------------------------------------------------------- helpers
const isObj = (x) => x !== null && typeof x === "object" && !Array.isArray(x);
const isStr = (x) => typeof x === "string";
const isInt = (x) => Number.isInteger(x);
const absent = (x) => x === undefined || x === null;

function readJson(rel) {
  const p = join(SITE, rel);
  if (!existsSync(p)) {
    err(rel, "file not found");
    return null;
  }
  let text;
  try {
    text = readFileSync(p, "utf8");
  } catch (e) {
    err(rel, `cannot read: ${e.message}`);
    return null;
  }
  for (const [re, what] of PRIVATE) {
    const m = text.match(re);
    if (m) err(rel, `contains a private item (${what}): ${JSON.stringify(m[0].trim().slice(0, 60))}`);
  }
  try {
    const data = JSON.parse(text);
    const mock = isObj(data) ? data.mock === true : Array.isArray(data) && data.length > 0 && data.every((x) => isObj(x) && x.mock === true);
    if (mock) mockFiles.push(rel);
    return data;
  } catch (e) {
    err(rel, `invalid JSON: ${e.message}`);
    return null;
  }
}

function fileExists(rel) {
  try {
    return statSync(join(SITE, rel)).isFile();
  } catch {
    return false;
  }
}

function checkM(where, m, { required = true } = {}) {
  if (absent(m)) {
    if (required) err(where, "math object required");
    return;
  }
  if (!isObj(m)) return err(where, `math object expected, found ${typeof m}`);
  if (!isStr(m.html) || !m.html.trim()) return err(where + ".html", "non-empty string required");
  if (!/^\s*<math[\s>]/.test(m.html)) err(where + ".html", "must start with a <math> element");
  if (!/<\/math>\s*$/.test(m.html)) err(where + ".html", "must end with </math>");
  if (/<script|\son[a-z]+\s*=/i.test(m.html)) err(where + ".html", "script or event handler inside mathematics");
  if (!isStr(m.text)) warn(where + ".text", "plain-text form missing");
  if (!(m.src === undefined || m.src === null || isStr(m.src))) err(where + ".src", "string or null expected");
}

// rich-text references are collected and checked once every file is loaded
const pendingRefs = [];

function checkR(where, r, { required = false } = {}) {
  if (absent(r)) {
    if (required) err(where, "rich text required");
    return;
  }
  if (!isStr(r)) return err(where, `rich text (string) expected, found ${typeof r}`);
  if (/<script|<style|<iframe|<img|\son[a-z]+\s*=/i.test(r)) err(where, "forbidden element or event handler in rich text");
  // tags outside <math>
  const outside = r.replace(/<math[\s\S]*?<\/math>/g, "");
  for (const m of outside.matchAll(/<\/?([A-Za-z][A-Za-z0-9-]*)/g)) {
    const tag = m[1].toLowerCase();
    if (!RICH_TAGS.has(tag)) warn(where, `tag <${tag}> is not in the rich-text set of DATA_FORMAT.md`);
  }
  if ((r.match(/<math[\s>]/g) || []).length !== (r.match(/<\/math>/g) || []).length) err(where, "unbalanced <math> elements");
  for (const m of outside.matchAll(/<a\b([^>]*)>/g)) {
    const attrs = m[1];
    const get = (name) => {
      const mm = attrs.match(new RegExp(`${name}\\s*=\\s*"([^"]*)"`));
      return mm ? mm[1] : null;
    };
    const page = get("data-report-page"), term = get("data-term"), cls = get("data-class"), th = get("data-theory");
    const href = get("href");
    if (page === null && term === null && cls === null && th === null && href === null) warn(where, `anchor without target: <a${attrs}>`);
    if (href !== null && /^(https?:)?\/\//.test(href)) warn(where, `external link ${href} (the site works offline)`);
    pendingRefs.push({ where, page, term, cls, th });
  }
}

function checkRefs(where, refs) {
  if (absent(refs)) return;
  if (!Array.isArray(refs)) return err(where, "array of report references expected");
  refs.forEach((x, i) => {
    const w = `${where}[${i}]`;
    if (!isObj(x)) return err(w, "object {label, page} expected");
    if (!isStr(x.label) || !x.label) err(w + ".label", "non-empty string required");
    pendingRefs.push({ where: w, page: x.page === undefined ? null : String(x.page), strictPage: true });
  });
}

function checkStr(where, x, { required = true, nonEmpty = true } = {}) {
  if (absent(x)) {
    if (required) err(where, "string required");
    return false;
  }
  if (!isStr(x)) {
    err(where, `string expected, found ${typeof x}`);
    return false;
  }
  if (nonEmpty && !x.trim()) err(where, "empty string");
  return true;
}

function checkIdList(where, list, known, what) {
  if (absent(list)) return;
  if (!Array.isArray(list)) return err(where, "array expected");
  list.forEach((id, i) => {
    if (!isStr(id)) err(`${where}[${i}]`, "string id expected");
    else if (!known.has(id)) err(`${where}[${i}]`, `unknown ${what} "${id}"`);
  });
}

// ---------------------------------------------------------------------------------------------- manifest
const manifest = readJson("data/manifest.json");
const phases = {};
const phaseFields = {};
let reportPages = null;
if (isObj(manifest)) {
  const W = "data/manifest.json";
  checkStr(W + ".title", manifest.title);
  checkStr(W + ".subtitle", manifest.subtitle);
  checkStr(W + ".model", manifest.model, { required: false });
  if (!absent(manifest.loop_order) && !isInt(manifest.loop_order)) err(W + ".loop_order", "integer expected");
  if (!absent(manifest.generated) && !(isStr(manifest.generated) && /^\d{4}-\d{2}-\d{2}$/.test(manifest.generated)))
    err(W + ".generated", "date YYYY-MM-DD expected");
  if (checkStr(W + ".report_pdf", manifest.report_pdf) && !fileExists(manifest.report_pdf)) err(W + ".report_pdf", `file ${manifest.report_pdf} not found`);
  if (!isInt(manifest.report_pages) || manifest.report_pages < 1) err(W + ".report_pages", "positive integer required");
  else reportPages = manifest.report_pages;
  checkR(W + ".status", manifest.status, { required: true });
  for (const k of ["classes", "theory", "resources"]) {
    if (checkStr(`${W}.${k}`, manifest[k]) && !fileExists(manifest[k])) err(`${W}.${k}`, `file ${manifest[k]} not found`);
  }
  if (!Array.isArray(manifest.phases) || manifest.phases.length === 0) err(W + ".phases", "non-empty array required");
  else {
    const seen = new Set();
    manifest.phases.forEach((p, i) => {
      const w = `${W}.phases[${i}]`;
      if (!isObj(p)) return err(w, "object expected");
      if (checkStr(w + ".id", p.id)) {
        if (seen.has(p.id)) err(w + ".id", `duplicate phase id "${p.id}"`);
        seen.add(p.id);
      }
      checkStr(w + ".title", p.title);
      if (checkStr(w + ".file", p.file) && fileExists(p.file)) phases[p.id] = { file: p.file, data: readJson(p.file) };
      else if (isStr(p.file)) err(w + ".file", `file ${p.file} not found`);
    });
    if (!seen.has("unbroken")) err(W + ".phases", 'no phase with id "unbroken" (parents of broken-phase terms refer to it)');
  }
} else if (manifest !== null) err("data/manifest.json", "object expected");

// ---------------------------------------------------------------------------------------------- classes (ids first)
const classesFile = isObj(manifest) && isStr(manifest.classes) ? manifest.classes : "data/classes.json";
const classesData = readJson(classesFile);
const classIds = new Map();
if (isObj(classesData)) {
  (Array.isArray(classesData.classes) ? classesData.classes : []).forEach((c) => {
    if (isObj(c) && isStr(c.id)) classIds.set(c.id, c);
  });
}

// ---------------------------------------------------------------------------------------------- theory (keys first)
const theoryFile = isObj(manifest) && isStr(manifest.theory) ? manifest.theory : "data/theory.json";
const theoryData = readJson(theoryFile);
const cardKeys = new Set(isObj(theoryData) && isObj(theoryData.cards) ? Object.keys(theoryData.cards) : []);

// ---------------------------------------------------------------------------------------------- phases
const allFieldIds = new Set();
const termIds = {};
for (const [pid, { file, data }] of Object.entries(phases)) {
  if (data === null) continue;
  const W = file;
  if (!isObj(data)) {
    err(W, "object expected");
    continue;
  }
  if (data.phase !== pid) err(W + ".phase", `"${data.phase}" does not match the manifest id "${pid}"`);
  checkStr(W + ".title", data.title);
  checkR(W + ".intro", data.intro);
  if (!absent(data.overview)) {
    if (!Array.isArray(data.overview)) err(W + ".overview", "array expected");
    else
      data.overview.forEach((o, i) => {
        const w = `${W}.overview[${i}]`;
        if (!isObj(o)) return err(w, "object expected");
        checkStr(w + ".id", o.id);
        checkStr(w + ".title", o.title);
        checkM(w + ".math", o.math);
        checkR(w + ".text", o.text);
      });
  }
  // fields
  const fields = new Map();
  if (!Array.isArray(data.fields)) err(W + ".fields", "array required");
  else
    data.fields.forEach((f, i) => {
      const w = `${W}.fields[${i}]`;
      if (!isObj(f)) return err(w, "object expected");
      if (checkStr(w + ".id", f.id)) {
        if (fields.has(f.id)) err(w + ".id", `duplicate field id "${f.id}"`);
        fields.set(f.id, f);
        allFieldIds.add(f.id);
      }
      checkM(w + ".symbol", f.symbol);
      checkStr(w + ".name", f.name);
      if (!KINDS.includes(f.kind)) err(w + ".kind", `"${f.kind}" is not one of ${KINDS.join(", ")}`);
      if (!(f.group === null || f.group === undefined || GROUPS.includes(f.group))) err(w + ".group", `"${f.group}" is not one of ${GROUPS.join(", ")} or null`);
      if (!isInt(f.ghost_number)) err(w + ".ghost_number", "integer required");
      if (!["bose", "fermi"].includes(f.statistics)) err(w + ".statistics", '"bose" or "fermi" required');
    });
  phaseFields[pid] = fields;
  // facets
  const facets = {};
  if (!isObj(data.facets)) err(W + ".facets", "object required");
  else
    for (const key of ["sector", "type", "group", "kind"]) {
      const list = data.facets[key];
      if (absent(list)) {
        if (key !== "kind") err(`${W}.facets.${key}`, "array required");
        continue;
      }
      if (!Array.isArray(list)) {
        err(`${W}.facets.${key}`, "array expected");
        continue;
      }
      facets[key] = new Map();
      list.forEach((x, i) => {
        const w = `${W}.facets.${key}[${i}]`;
        if (!isObj(x) || !isStr(x.id) || !isStr(x.label)) return err(w, "object {id, label} with strings required");
        if (facets[key].has(x.id)) err(w, `duplicate id "${x.id}"`);
        facets[key].set(x.id, 0);
      });
      if (key === "type") for (const id of facets.type.keys()) if (!TYPES.includes(id)) err(`${W}.facets.type`, `unknown type id "${id}" (allowed: ${TYPES.join(", ")})`);
    }
  // terms
  const ids = new Set();
  termIds[pid] = ids;
  if (!Array.isArray(data.terms) || data.terms.length === 0) {
    err(W + ".terms", "non-empty array required");
    continue;
  }
  data.terms.forEach((t, i) => {
    let w = `${W}.terms[${i}]`;
    if (!isObj(t)) return err(w, "object expected");
    if (isStr(t.id)) w += ` (${t.id})`;
    if (!isStr(t.id) || !TERM_ID.test(t.id)) err(w + ".id", "string of [A-Za-z0-9_-] required");
    else if (ids.has(t.id)) err(w + ".id", "duplicate term id");
    else ids.add(t.id);
  });
}

// second pass over terms (needs every phase's ids)
for (const [pid, { file, data }] of Object.entries(phases)) {
  if (!isObj(data) || !Array.isArray(data.terms)) continue;
  const W = file;
  const fields = phaseFields[pid] || new Map();
  const facets = {};
  if (isObj(data.facets)) for (const k of ["sector", "type", "group", "kind"]) if (Array.isArray(data.facets[k])) facets[k] = new Map(data.facets[k].filter(isObj).map((x) => [x.id, 0]));
  data.terms.forEach((t, i) => {
    if (!isObj(t)) return;
    const w = `${W}.terms[${i}]${isStr(t.id) ? ` (${t.id})` : ""}`;
    checkStr(w + ".name", t.name);
    checkM(w + ".operator", t.operator);
    if (!Array.isArray(t.fields) || t.fields.length === 0) err(w + ".fields", "non-empty array required");
    else {
      t.fields.forEach((f, j) => {
        if (!fields.has(f)) err(`${w}.fields[${j}]`, `unknown field "${f}" (not in ${W}.fields)`);
      });
      if (isInt(t.legs) && t.legs !== t.fields.length) warn(w + ".legs", `${t.legs} legs but ${t.fields.length} fields`);
      const gh = t.fields.reduce((s, f) => s + (fields.get(f)?.ghost_number ?? 0), 0);
      if (isInt(t.ghost_number) && t.ghost_number !== gh) warn(w + ".ghost_number", `${t.ghost_number}, but the fields add up to ${gh}`);
    }
    if (!isInt(t.legs)) err(w + ".legs", "integer required");
    if (!isInt(t.ghost_number)) err(w + ".ghost_number", "integer required");
    if (!isStr(t.sector)) err(w + ".sector", "string required");
    else if (facets.sector && !facets.sector.has(t.sector)) err(w + ".sector", `"${t.sector}" is not in facets.sector`);
    else if (facets.sector) facets.sector.set(t.sector, facets.sector.get(t.sector) + 1);
    if (!Array.isArray(t.groups)) err(w + ".groups", "array required (may be empty)");
    else
      t.groups.forEach((g, j) => {
        if (facets.group && !facets.group.has(g)) err(`${w}.groups[${j}]`, `"${g}" is not in facets.group`);
        else if (facets.group) facets.group.set(g, facets.group.get(g) + 1);
      });
    if (typeof t.evanescent !== "boolean") err(w + ".evanescent", "boolean required");
    const types = Array.isArray(t.types) ? t.types : [];
    if (!Array.isArray(t.types) || t.types.length === 0) err(w + ".types", "non-empty array required");
    types.forEach((ty, j) => {
      if (!TYPES.includes(ty)) err(`${w}.types[${j}]`, `unknown type "${ty}"`);
      else if (facets.type && !facets.type.has(ty)) err(`${w}.types[${j}]`, `"${ty}" is not in facets.type`);
      else if (facets.type) facets.type.set(ty, facets.type.get(ty) + 1);
    });
    if (new Set(types).size !== types.length) err(w + ".types", "duplicate entries");
    if (types.includes("evanescent") && t.evanescent === false) warn(w, 'type "evanescent" on a term with evanescent: false');
    // coefficients
    const co = t.coefficients;
    if (!isObj(co)) err(w + ".coefficients", "object required");
    else {
      for (const k of Object.keys(co)) {
        if (!TYPES.includes(k)) {
          warn(`${w}.coefficients.${k}`, "unknown part (ignored by the page)");
          continue;
        }
        if (!types.includes(k)) err(`${w}.coefficients.${k}`, `coefficient given for a type not listed in types`);
        const c = co[k];
        if (!isObj(c)) {
          err(`${w}.coefficients.${k}`, "object {math, ...} expected");
          continue;
        }
        checkM(`${w}.coefficients.${k}.math`, c.math);
        checkM(`${w}.coefficients.${k}.decomposition`, c.decomposition, { required: false });
        checkM(`${w}.coefficients.${k}.symbol`, c.symbol, { required: false });
        checkR(`${w}.coefficients.${k}.note`, c.note);
      }
      for (const ty of types) {
        if (!(ty in co)) (ty === "breaking" ? warn : err)(`${w}.coefficients`, `no entry for the listed type "${ty}"`);
      }
    }
    checkM(w + ".total", t.total, { required: false });
    if (absent(t.total)) warn(w + ".total", "missing");
    // rules
    if (!absent(t.rules)) {
      if (!Array.isArray(t.rules)) err(w + ".rules", "array or null expected");
      else
        t.rules.forEach((r, j) => {
          const wr = `${w}.rules[${j}]`;
          if (!isObj(r)) return err(wr, "object expected");
          checkStr(wr + ".title", r.title, { required: false });
          if (!Array.isArray(r.legs)) err(wr + ".legs", "array required");
          else
            r.legs.forEach((l, k) => {
              if (!isObj(l)) return err(`${wr}.legs[${k}]`, "object expected");
              if (!fields.has(l.field)) err(`${wr}.legs[${k}].field`, `unknown field "${l.field}"`);
              checkM(`${wr}.legs[${k}].label`, l.label, { required: false });
            });
          if (!isObj(r.parts)) err(wr + ".parts", "object required");
          else
            for (const [k, m] of Object.entries(r.parts)) {
              if (!TYPES.includes(k)) warn(`${wr}.parts.${k}`, "unknown part (ignored by the page)");
              else if (!types.includes(k)) err(`${wr}.parts.${k}`, "part for a type not listed in types");
              checkM(`${wr}.parts.${k}`, m, { required: false });
            }
          checkR(wr + ".note", r.note);
        });
    }
    // classes
    if (!absent(t.classes)) {
      if (!isObj(t.classes)) err(w + ".classes", "object {pole, breaking} or null expected");
      else
        for (const k of ["pole", "breaking"]) {
          checkIdList(`${w}.classes.${k}`, t.classes[k], classIds, "class");
          if (Array.isArray(t.classes[k]))
            t.classes[k].forEach((cid, j) => {
              const c = classIds.get(cid);
              if (!c) return;
              if (k === "breaking" && c.insertion !== "Delta1") warn(`${w}.classes.breaking[${j}]`, `class "${cid}" has no insertion of Delta1`);
              if (k === "pole" && c.insertion === "Delta1") warn(`${w}.classes.pole[${j}]`, `class "${cid}" carries an insertion of Delta1`);
            });
        }
    }
    if (!isStr(t.theory)) err(w + ".theory", "card key required");
    else if (!cardKeys.has(t.theory)) err(w + ".theory", `no card "${t.theory}" in ${theoryFile}`);
    checkR(w + ".remarks", t.remarks);
    checkRefs(w + ".report", t.report);
    if (facets.kind) {
      if (absent(t.kind)) warn(w + ".kind", "missing although the phase has facets.kind");
      else if (!facets.kind.has(t.kind)) err(w + ".kind", `"${t.kind}" is not in facets.kind`);
      else facets.kind.set(t.kind, facets.kind.get(t.kind) + 1);
    }
    if (!absent(t.parents)) {
      if (pid === "unbroken") warn(w + ".parents", "parents in the unbroken phase are ignored");
      checkIdList(w + ".parents", t.parents, termIds.unbroken || new Set(), "unbroken term");
    } else if (pid !== "unbroken") warn(w + ".parents", "missing for a term outside the unbroken phase");
  });
  for (const [k, m] of Object.entries(facets)) for (const [id, n] of m) if (n === 0) warn(`${W}.facets.${k}`, `value "${id}" is used by no term`);
}

// ---------------------------------------------------------------------------------------------- classes.json
const graphFiles = [];
if (isObj(classesData)) {
  const W = classesFile;
  const groupKeys = new Map();
  if (!Array.isArray(classesData.groups)) err(W + ".groups", "array required");
  else
    classesData.groups.forEach((g, i) => {
      const w = `${W}.groups[${i}]`;
      if (!isObj(g)) return err(w, "object expected");
      if (checkStr(w + ".key", g.key)) {
        if (groupKeys.has(g.key)) err(w + ".key", `duplicate group key "${g.key}"`);
        groupKeys.set(g.key, g);
      }
      checkStr(w + ".title", g.title);
      if (!isInt(g.ghost_number)) err(w + ".ghost_number", "integer required");
      checkIdList(w + ".classes", g.classes, classIds, "class");
      if (!Array.isArray(g.classes)) err(w + ".classes", "array required");
      checkRefs(w + ".report", g.report);
      checkR(w + ".summary", g.summary);
    });
  const listedIn = new Map();
  for (const g of groupKeys.values()) for (const c of Array.isArray(g.classes) ? g.classes : []) {
    if (listedIn.has(c)) err(`${W}.groups`, `class "${c}" listed in groups "${listedIn.get(c)}" and "${g.key}"`);
    listedIn.set(c, g.key);
  }
  if (!Array.isArray(classesData.classes)) err(W + ".classes", "array required");
  else {
    const seen = new Set();
    classesData.classes.forEach((c, i) => {
      let w = `${W}.classes[${i}]`;
      if (!isObj(c)) return err(w, "object expected");
      if (isStr(c.id)) w += ` (${c.id})`;
      if (!isStr(c.id) || !CLASS_ID.test(c.id)) err(w + ".id", "string of [A-Za-z0-9_.:,-] required");
      else if (seen.has(c.id)) err(w + ".id", "duplicate class id");
      seen.add(c.id);
      if (!groupKeys.has(c.group)) err(w + ".group", `unknown group "${c.group}"`);
      else if (listedIn.get(c.id) !== c.group) err(w + ".group", `class is not listed in groups["${c.group}"].classes`);
      if (!isInt(c.ghost_number)) err(w + ".ghost_number", "integer required");
      else if (groupKeys.has(c.group) && groupKeys.get(c.group).ghost_number !== c.ghost_number) warn(w + ".ghost_number", "differs from the ghost number of its group");
      if (!Array.isArray(c.fields)) err(w + ".fields", "array required");
      else c.fields.forEach((f, j) => { if (!allFieldIds.has(f)) warn(`${w}.fields[${j}]`, `field "${f}" is defined in no phase`); });
      checkM(w + ".label", c.label);
      if (!(c.insertion === null || c.insertion === undefined || c.insertion === "Delta1")) err(w + ".insertion", 'null or "Delta1" required');
      if (!isInt(c.n_graphs) || c.n_graphs < 0) err(w + ".n_graphs", "non-negative integer required");
      if (!absent(c.graphs)) {
        if (!isStr(c.graphs)) err(w + ".graphs", "path string expected");
        else if (!fileExists(c.graphs)) err(w + ".graphs", `file ${c.graphs} not found`);
        else graphFiles.push({ cls: c, file: c.graphs, where: w });
      } else if (c.n_graphs > 0) warn(w + ".graphs", "no graph file although n_graphs > 0");
      if (!absent(c.result)) {
        if (!isObj(c.result)) err(w + ".result", "object expected");
        else for (const k of ["pole", "breaking"]) checkM(`${w}.result.${k}`, c.result[k], { required: false });
      }
      checkR(w + ".summary", c.summary);
      if (!listedIn.has(c.id)) err(w, "class is listed in no group");
    });
  }
} else if (classesData !== null) err(classesFile, "object expected");

// ---------------------------------------------------------------------------------------------- graph files
const graphIds = new Map();
let nGraphs = 0;
for (const { cls, file } of graphFiles) {
  const g = readJson(file);
  if (!isObj(g)) {
    if (g !== null) err(file, "object expected");
    continue;
  }
  if (g.class !== cls.id) err(file + ".class", `"${g.class}" but the file is referenced by class "${cls.id}"`);
  if (!Array.isArray(g.graphs)) {
    err(file + ".graphs", "array required");
    continue;
  }
  if (isInt(cls.n_graphs) && cls.n_graphs !== g.graphs.length) err(file, `${g.graphs.length} graphs, but n_graphs of class "${cls.id}" is ${cls.n_graphs}`);
  nGraphs += g.graphs.length;
  const local = new Set();
  g.graphs.forEach((x, i) => {
    let w = `${file}.graphs[${i}]`;
    if (!isObj(x)) return err(w, "object expected");
    if (isStr(x.id)) w += ` (${x.id})`;
    if (!isStr(x.id) || !x.id) err(w + ".id", "string required");
    else {
      if (local.has(x.id)) err(w + ".id", "duplicate graph id in the file");
      else if (graphIds.has(x.id)) warn(w + ".id", `graph id also used in ${graphIds.get(x.id)}`);
      local.add(x.id);
      graphIds.set(x.id, file);
    }
    if (!isStr(x.svg) || !/^\s*<svg[\s>]/.test(x.svg)) err(w + ".svg", "string starting with <svg required");
    else {
      if (!/viewBox\s*=/.test(x.svg.slice(0, 2000))) warn(w + ".svg", "no viewBox on the root element");
      if (/<script|\son[a-z]+\s*=/i.test(x.svg)) err(w + ".svg", "script or event handler in the drawing");
      if (/(?:href|src)\s*=\s*"(?:https?:)?\/\//i.test(x.svg) || /url\(\s*['"]?(?:https?:)?\/\//i.test(x.svg)) err(w + ".svg", "external resource in the drawing");
    }
    if (!Array.isArray(x.loop)) warn(w + ".loop", "array of field ids expected");
    else x.loop.forEach((f, j) => { if (!allFieldIds.has(f)) warn(`${w}.loop[${j}]`, `field "${f}" is defined in no phase`); });
    if (!(x.insertion === null || x.insertion === undefined || x.insertion === "Delta1")) err(w + ".insertion", 'null or "Delta1" required');
    for (const k of ["topology", "symmetry_factor", "sign"]) if (!absent(x[k]) && !isStr(x[k]) && typeof x[k] !== "number") err(`${w}.${k}`, "string expected");
    checkR(w + ".caption", x.caption);
    if (!absent(x.contributes)) {
      if (!isObj(x.contributes)) err(w + ".contributes", "object expected");
      // null: not determined for this graph (e.g. the evanescent flag, or the pole flag at ghost number one)
      else for (const [k, v] of Object.entries(x.contributes)) if (typeof v !== "boolean" && v !== null) err(`${w}.contributes.${k}`, "boolean or null expected");
    }
  });
}

// ---------------------------------------------------------------------------------------------- theory.json
if (isObj(theoryData)) {
  const W = theoryFile;
  if (!isObj(theoryData.cards)) err(W + ".cards", "object required");
  else {
    for (const [k, c] of Object.entries(theoryData.cards)) {
      const w = `${W}.cards["${k}"]`;
      if (!isObj(c)) {
        err(w, "object expected");
        continue;
      }
      checkStr(w + ".title", c.title);
      checkR(w + ".body", c.body, { required: true });
      checkRefs(w + ".report", c.report);
      if (!CARD_KEYS.includes(k) && !METHOD_KEYS.includes(k)) warn(w, "key not in the table of DATA_FORMAT.md section 7");
    }
    for (const k of CARD_KEYS) if (!(k in theoryData.cards)) warn(W + ".cards", `card "${k}" of DATA_FORMAT.md section 7 is missing`);
  }
} else if (theoryData !== null) err(theoryFile, "object expected");

// ---------------------------------------------------------------------------------------------- resources.json
const resFile = isObj(manifest) && isStr(manifest.resources) ? manifest.resources : "data/resources.json";
const resData = readJson(resFile);
if (resData !== null) {
  const list = Array.isArray(resData) ? resData : isObj(resData) && Array.isArray(resData.resources) ? resData.resources : null;
  if (!list) err(resFile, "list of resources expected");
  else
    list.forEach((r, i) => {
      const w = `${resFile}[${i}]`;
      if (!isObj(r)) return err(w, "object expected");
      if (checkStr(w + ".path", r.path)) {
        if (/^\/|^[a-z]+:/i.test(r.path) || r.path.includes("..")) err(w + ".path", "relative path inside the site required");
        else if (!fileExists(r.path)) err(w + ".path", `file ${r.path} not found`);
        else if (isInt(r.bytes) && statSync(join(SITE, r.path)).size !== r.bytes) warn(w + ".bytes", `${r.bytes}, but the file has ${statSync(join(SITE, r.path)).size} bytes`);
      }
      checkStr(w + ".title", r.title);
      checkStr(w + ".description", r.description, { required: false });
      if (!absent(r.bytes) && !isInt(r.bytes)) err(w + ".bytes", "integer expected");
    });
}

// ---------------------------------------------------------------------------------------------- cross references
for (const ref of pendingRefs) {
  if (ref.page !== null && ref.page !== undefined) {
    const n = Number(ref.page);
    if (!Number.isInteger(n) || n < 1) err(ref.where, `report page "${ref.page}" is not a positive integer`);
    else if (reportPages && n > reportPages) err(ref.where, `report page ${n} beyond the last page ${reportPages}`);
  } else if (ref.strictPage) err(ref.where + ".page", "integer page required");
  if (ref.term) {
    const [ph, ...rest] = ref.term.split(":");
    const id = rest.join(":");
    if (!id) err(ref.where, `data-term="${ref.term}" must read PHASE:ID`);
    else if (!termIds[ph]) err(ref.where, `data-term="${ref.term}": unknown phase "${ph}"`);
    else if (!termIds[ph].has(id)) err(ref.where, `data-term="${ref.term}": no term "${id}" in phase "${ph}"`);
  }
  if (ref.cls && !classIds.has(ref.cls)) err(ref.where, `data-class="${ref.cls}": unknown class`);
  if (ref.th && !cardKeys.has(ref.th)) err(ref.where, `data-theory="${ref.th}": unknown card`);
}

// ---------------------------------------------------------------------------------------------- report
if (mockFiles.length) {
  const msg = `${mockFiles.length} file(s) carry "mock": true: ${mockFiles.slice(0, 8).join(", ")}${mockFiles.length > 8 ? ", ..." : ""}`;
  if (opts.noMock) err("mock", msg);
  else console.log(`note: ${msg}`);
}
const nTerms = Object.fromEntries(Object.entries(termIds).map(([k, v]) => [k, v.size]));
console.log(
  `checked ${SITE}: terms ${JSON.stringify(nTerms)}, ${classIds.size} classes, ${nGraphs} graphs in ${graphFiles.length} files, ${cardKeys.size} theory cards`,
);
const show = (label, list) => {
  if (!list.length) return;
  console.log(`\n${list.length} ${label}:`);
  for (const m of list.slice(0, opts.max)) console.log(`  ${m}`);
  if (list.length > opts.max) console.log(`  ... ${list.length - opts.max} more`);
};
show("error(s)", errors);
show("warning(s)", warnings);
const fail = errors.length > 0 || (opts.strict && warnings.length > 0);
console.log(fail ? "\nFAILED" : "\nOK");
process.exit(fail ? 1 : 0);
