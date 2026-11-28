// Single-page reference for the renormalised Standard Model with a BMHV gamma5.
// Vanilla ES module, no dependency, no build step. Data contract: DATA_FORMAT.md (site/data/*.json).
//
// Routes (hash based, every state is a shareable link):
//   #/                                   home
//   #/phase/PHASE[?filters]              operator table of a phase
//   #/phase/PHASE/term/ID[?filters]      the same with the detail of one term
//   #/classes[/CLASS_ID]                 Green's functions by group (optionally one class opened)
//   #/theory[/KEY]                       theory cards
//   #/report[/PAGE]                      the report (PDF)
//   #/resources                          downloadable inputs and results
// Filters: q (text), f (field, repeated), fm=any|all, t (type), s (sector), g (group), k (kind), l (legs),
// gh (ghost number), sort, dir.

"use strict";

// ------------------------------------------------------------------------------------------------ constants

const TYPE_ORDER = ["tree", "pole", "evanescent", "finite", "breaking"];
const TYPE_SHORT = { tree: "Tree", pole: "Singular", evanescent: "Evanescent", finite: "Finite restoring", breaking: "Breaking" };
const TYPE_LONG = {
  tree: "Tree level",
  pole: "Singular counterterm",
  evanescent: "Evanescent counterterm",
  finite: "Finite restoring counterterm",
  breaking: "Breaking insertion",
};
// Fallback only: the build writes these definitions with MathML into the manifest ("types").
const TYPE_DEF = {
  tree: "Classical action, including gauge fixing, ghosts and the BRST source terms.",
  pole: "Singular counterterm on a non-evanescent operator, continued to d dimensions.",
  evanescent: "Singular counterterm on an evanescent operator, one that vanishes in four dimensions.",
  finite: "Finite symmetry-restoring counterterm, which restores the Slavnov–Taylor identity.",
  breaking: "Insertion of the tree-level breaking. It is not part of the action; it is listed for its Feynman rule.",
};
const KIND_ORDER = ["gauge", "fermion", "sterile", "scalar", "ghost", "antighost", "auxiliary", "source"];
const KIND_LABEL = {
  gauge: "Gauge fields", fermion: "Fermions", sterile: "Sterile partners", scalar: "Scalars", ghost: "Ghosts",
  antighost: "Antighosts", auxiliary: "Auxiliary fields", source: "BRST sources",
};
const WIDE = "(min-width: 1180px)";
const GRAPH_PAGE = 48;
const THEME_KEY = "gamma5docs.theme";

// ------------------------------------------------------------------------------------------------ small helpers

const esc = (s) =>
  String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
const enc = (s) => encodeURIComponent(String(s));
const isObj = (x) => x !== null && typeof x === "object" && !Array.isArray(x);
const arr = (x) => (Array.isArray(x) ? x : []);
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const hasMath = (m) => isObj(m) && typeof m.html === "string" && m.html.trim() !== "";
// rich name of a term (symbols as MathML, written by the build); the plain name is the fallback and the tooltip
const nameHTML = (t) => (typeof t?.name_html === "string" && t.name_html.trim() ? t.name_html : esc(t?.name || t?.id || ""));
const nameEl = (tag, t, props = {}) => h(tag, { ...props, html: nameHTML(t), title: props.title ?? null });
const smooth = () => (window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth");

function h(tag, props, ...kids) {
  const el = document.createElement(tag);
  if (props) {
    for (const [k, v] of Object.entries(props)) {
      if (v === null || v === undefined || v === false) continue;
      if (k === "class") el.className = v;
      else if (k === "html") el.innerHTML = v;
      else if (k === "text") el.textContent = v;
      else if (k === "dataset") Object.assign(el.dataset, v);
      else if (k.startsWith("on") && typeof v === "function") el.addEventListener(k.slice(2), v);
      else if (v === true) el.setAttribute(k, "");
      else el.setAttribute(k, String(v));
    }
  }
  appendKids(el, kids);
  return el;
}

function appendKids(el, kids) {
  for (const k of kids.flat(Infinity)) {
    if (k === null || k === undefined || k === false) continue;
    el.append(k instanceof Node ? k : String(k));
  }
}

function plural(n, one, many = one + "s") {
  return `${n.toLocaleString("en-GB")} ${n === 1 ? one : many}`;
}

function fmtBytes(n) {
  if (!Number.isFinite(n)) return "";
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(n < 10240 ? 1 : 0)} kB`;
  return `${(n / 1024 / 1024).toFixed(1)} MB`;
}

// ------------------------------------------------------------------------------------------------ theme

const theme = {
  get() {
    try {
      const t = localStorage.getItem(THEME_KEY);
      return t === "light" || t === "dark" ? t : "system";
    } catch {
      return "system";
    }
  },
  set(t) {
    try {
      if (t === "system") localStorage.removeItem(THEME_KEY);
      else localStorage.setItem(THEME_KEY, t);
    } catch {
      /* storage unavailable: the choice lasts for this page only */
    }
    if (t === "system") delete document.documentElement.dataset.theme;
    else document.documentElement.dataset.theme = t;
    theme.paint(t);
  },
  paint(t) {
    const btn = document.getElementById("theme-toggle");
    if (!btn) return;
    const label = { system: "Theme: follows the system", light: "Theme: light (Alt+click: follow the system)", dark: "Theme: dark (Alt+click: follow the system)" }[t];
    btn.setAttribute("aria-label", label + " (activate to change)");
    btn.title = label;
    btn.dataset.mode = t;
    const txt = btn.querySelector(".theme-label");
    if (txt) txt.textContent = { system: "Auto", light: "Light", dark: "Dark" }[t];
  },
  init() {
    const t = theme.get();
    if (t !== "system") document.documentElement.dataset.theme = t;
    theme.paint(t);
    // light <-> dark; the system setting is the starting point and is restored by Alt+click
    document.getElementById("theme-toggle")?.addEventListener("click", (e) => {
      if (e.altKey) return theme.set("system");
      const cur = theme.get();
      const shown = cur === "system" ? (window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light") : cur;
      theme.set(shown === "dark" ? "light" : "dark");
    });
  },
};

// ------------------------------------------------------------------------------------------------ data loading

const cache = new Map();
const mockFiles = new Set();

function loadJSON(path) {
  if (!cache.has(path)) {
    const p = fetch(path, { cache: "no-cache" })
      .then((r) => {
        if (!r.ok) throw new Error(`${path}: the server answered ${r.status} ${r.statusText || ""}`.trim());
        return r.json();
      })
      .then((data) => {
        const mock = isObj(data)
          ? data.mock === true
          : Array.isArray(data) && data.length > 0 && data.every((x) => isObj(x) && x.mock === true);
        if (mock) {
          mockFiles.add(path);
          paintMockNotice();
        }
        return data;
      })
      .catch((e) => {
        cache.delete(path); // a later visit may retry
        throw e instanceof SyntaxError ? new Error(`${path}: invalid JSON (${e.message})`) : e;
      });
    cache.set(path, p);
  }
  return cache.get(path);
}

const store = { manifest: null, phases: new Map(), classes: null, theory: null };

async function getManifest() {
  if (!store.manifest) store.manifest = await loadJSON("data/manifest.json");
  return store.manifest;
}

function phaseEntries(man) {
  return arr(man?.phases).filter((p) => isObj(p) && typeof p.id === "string");
}

async function getPhase(pid) {
  if (store.phases.has(pid)) return store.phases.get(pid);
  const man = await getManifest();
  const entry = phaseEntries(man).find((p) => p.id === pid);
  if (!entry) throw new Error(`There is no phase "${pid}" in the manifest.`);
  const data = await loadJSON(entry.file || `data/${pid}.json`);
  const P = indexPhase(pid, entry, data);
  store.phases.set(pid, P);
  return P;
}

async function getAllPhases() {
  const man = await getManifest();
  const out = [];
  for (const e of phaseEntries(man)) {
    try {
      out.push(await getPhase(e.id));
    } catch (err) {
      console.warn(err);
    }
  }
  return out;
}

function indexPhase(pid, entry, data) {
  data = isObj(data) ? data : {};
  const fields = arr(data.fields).filter((f) => isObj(f) && typeof f.id === "string");
  const fieldById = new Map(fields.map((f) => [f.id, f]));
  const terms = arr(data.terms).filter((t) => isObj(t) && typeof t.id === "string");
  const termById = new Map(terms.map((t) => [t.id, t]));
  const facets = isObj(data.facets) ? data.facets : {};
  const facet = (k) => arr(facets[k]).filter((x) => isObj(x) && x.id !== undefined);
  const typeLabels = { ...TYPE_SHORT };
  for (const x of facet("type")) if (typeof x.label === "string") typeLabels[x.id] = x.label;
  terms.forEach((t, i) => {
    t._order = i;
    t._types = TYPE_ORDER.filter((k) => arr(t.types).includes(k));
    const fnames = arr(t.fields).map((f) => `${f} ${fieldById.get(f)?.name ?? ""} ${fieldById.get(f)?.symbol?.text ?? ""}`);
    t._search = [t.id, t.name, t.operator?.text, t.sector, ...fnames, ...arr(t.groups)].join(" ").toLowerCase();
  });
  return {
    id: pid, title: data.title || entry.title || pid, entry, data, fields, fieldById, terms, termById,
    facets: { sector: facet("sector"), type: facet("type"), group: facet("group"), kind: facet("kind") },
    typeLabels,
  };
}

async function getClasses() {
  if (store.classes) return store.classes;
  const man = await getManifest();
  const data = await loadJSON(man.classes || "data/classes.json");
  const groups = arr(data?.groups).filter((g) => isObj(g) && typeof g.key === "string");
  const classes = arr(data?.classes).filter((c) => isObj(c) && typeof c.id === "string");
  const classById = new Map(classes.map((c) => [c.id, c]));
  const groupByKey = new Map(groups.map((g) => [g.key, g]));
  // classes not listed in any group still get a home
  const listed = new Set(groups.flatMap((g) => arr(g.classes)));
  const orphans = classes.filter((c) => !listed.has(c.id));
  if (orphans.length) {
    for (const c of orphans) {
      const g = groupByKey.get(c.group);
      if (g) g.classes = [...arr(g.classes), c.id];
    }
    const stray = orphans.filter((c) => !groupByKey.has(c.group));
    if (stray.length) {
      const g = { key: "_other", title: "Other classes", ghost_number: null, classes: stray.map((c) => c.id) };
      groups.push(g);
      groupByKey.set(g.key, g);
    }
  }
  store.classes = { data, groups, classes, classById, groupByKey };
  return store.classes;
}

async function getTheory() {
  if (store.theory) return store.theory;
  const man = await getManifest();
  const data = await loadJSON(man.theory || "data/theory.json");
  store.theory = { cards: isObj(data?.cards) ? data.cards : {} };
  return store.theory;
}

// field symbols across phases (for graph captions)
function anyField(id) {
  for (const P of store.phases.values()) {
    const f = P.fieldById.get(id);
    if (f) return f;
  }
  return null;
}

// reverse index class -> terms (built from every loaded phase)
function termsOfClass(cid) {
  const out = [];
  for (const P of store.phases.values()) {
    for (const t of P.terms) {
      const c = isObj(t.classes) ? t.classes : {};
      const roles = [];
      if (arr(c.pole).includes(cid)) roles.push("pole");
      if (arr(c.breaking).includes(cid)) roles.push("breaking");
      if (roles.length) out.push({ P, t, roles });
    }
  }
  return out;
}

// ------------------------------------------------------------------------------------------------ mathematics

let viewerPromise = null;
let viewerModule = null;

function loadViewer() {
  if (!viewerPromise) {
    viewerPromise = import("./viewer/mathview.js")
      .then((m) => {
        viewerModule = m && typeof m.renderMath === "function" ? m : null;
        return viewerModule;
      })
      .catch(() => {
        viewerModule = null;
        return null; // the static MathML stays as it is (permanent fallback)
      });
  }
  return viewerPromise;
}

async function viewerReady(timeout = 700) {
  const m = await Promise.race([loadViewer(), sleep(timeout).then(() => undefined)]);
  if (m && typeof m.ready === "function") {
    try {
      await Promise.race([m.ready(), sleep(timeout)]);
    } catch {
      /* font not loaded: render anyway */
    }
  }
}

// Fill el with the expression. Static MathML goes in first (no empty box, no layout jump); the interactive viewer,
// when available, replaces it.
function mathInto(el, M, { interactive = true } = {}) {
  if (!hasMath(M)) {
    el.textContent = "—";
    return el;
  }
  // views build their content detached and insert it later: wait until the box is in the document
  const upgrade = (mod, tries = 0) => {
    if (!mod) return;
    if (!el.isConnected) {
      if (tries < 400) (tries < 120 ? requestAnimationFrame : (f) => setTimeout(f, 100))(() => upgrade(mod, tries + 1));
      return;
    }
    try {
      const r = mod.renderMath(el, M, { interactive });
      if (r && typeof r.catch === "function") r.catch(() => (el.innerHTML = M.html));
    } catch (e) {
      console.warn("mathview failed; static MathML kept", e);
      el.innerHTML = M.html;
    }
  };
  el.innerHTML = M.html;
  if (viewerModule) queueMicrotask(() => upgrade(viewerModule));
  else if (interactive) loadViewer().then(upgrade);
  return el;
}

function mathBox(M, { cls = "", interactive = true, label = null } = {}) {
  const box = h("div", { class: `mbox ${cls}`.trim(), role: label ? "group" : null, "aria-label": label });
  mathInto(box, M, { interactive });
  return box;
}

const staticMath = (M) => (hasMath(M) ? M.html : "");

// ------------------------------------------------------------------------------------------------ rich text

const RICH_DROP = "script,style,iframe,object,embed,link,meta,form,input,button,textarea,select";

function rich(html, { tag = "div", cls = "rich" } = {}) {
  const el = document.createElement(tag);
  el.className = cls;
  if (typeof html !== "string" || !html.trim()) return el;
  const tpl = document.createElement("template");
  tpl.innerHTML = html;
  tpl.content.querySelectorAll(RICH_DROP).forEach((n) => n.remove());
  tpl.content.querySelectorAll("*").forEach((n) => {
    for (const a of [...n.attributes]) {
      if (/^on/i.test(a.name) || (/^(href|src|xlink:href)$/i.test(a.name) && /^\s*javascript:/i.test(a.value))) n.removeAttribute(a.name);
    }
  });
  tpl.content.querySelectorAll("a").forEach(linkify);
  el.append(tpl.content);
  return el;
}

function linkify(a) {
  const d = a.dataset;
  if (d.reportPage) {
    a.href = `#/report/${enc(d.reportPage)}`;
    a.classList.add("xref", "xref-report");
    a.title ||= `Report, page ${d.reportPage}`;
  } else if (d.term) {
    const i = d.term.indexOf(":");
    const [ph, id] = i > 0 ? [d.term.slice(0, i), d.term.slice(i + 1)] : ["unbroken", d.term];
    a.href = `#/phase/${enc(ph)}/term/${enc(id)}`;
    a.classList.add("xref", "xref-term");
  } else if (d.class) {
    a.href = `#/classes/${enc(d.class)}`;
    a.classList.add("xref", "xref-class");
  } else if (d.theory) {
    a.href = `#/theory/${enc(d.theory)}`;
    a.classList.add("xref", "xref-theory");
  } else if (a.getAttribute("href") && /^(https?:)?\/\//.test(a.getAttribute("href"))) {
    a.target = "_blank";
    a.rel = "noopener";
  }
}

// ------------------------------------------------------------------------------------------------ shared pieces

function typeTag(type, labels = TYPE_SHORT, { long = false } = {}) {
  const text = long ? TYPE_LONG[type] || labels[type] || type : labels[type] || TYPE_SHORT[type] || type;
  return h("span", { class: `tag t-${type}` }, h("span", { class: "tag-mark", "aria-hidden": "true" }), text);
}

function typeTagsHTML(types, labels) {
  return types
    .map((t) => `<span class="tag t-${esc(t)}"><span class="tag-mark" aria-hidden="true"></span>${esc(labels[t] || TYPE_SHORT[t] || t)}</span>`)
    .join("");
}

function fieldChipHTML(f, id) {
  const sym = f && hasMath(f.symbol) ? f.symbol.html : esc(id);
  const name = f?.name ? esc(f.name) : esc(id);
  return `<span class="fchip k-${esc(f?.kind || "unknown")}" title="${name}"><span class="fchip-sym">${sym}</span><span class="sr-only">${name}</span></span>`;
}

function fieldChip(P, id, { href = null } = {}) {
  const f = P?.fieldById.get(id) || anyField(id);
  const tag = href ? "a" : "span";
  const el = h(tag, { class: `fchip k-${f?.kind || "unknown"}`, href, title: f?.name || id });
  el.innerHTML = `<span class="fchip-sym">${f && hasMath(f.symbol) ? f.symbol.html : esc(id)}</span>${
    href ? `<span class="fchip-name">${esc(f?.name || id)}</span>` : `<span class="sr-only">${esc(f?.name || id)}</span>`
  }`;
  return el;
}

function reportRefs(refs, { cls = "refs" } = {}) {
  const list = arr(refs).filter((r) => isObj(r) && Number.isFinite(Number(r.page)));
  if (!list.length) return null;
  return h(
    "ul",
    { class: cls, "aria-label": "Report references" },
    list.map((r) =>
      h("li", null, h("a", { class: "ref", href: `#/report/${enc(r.page)}` }, h("span", { class: "ref-label" }, r.label || "Report"), h("span", { class: "ref-page" }, `p. ${r.page}`))),
    ),
  );
}

function errorBox(err, what = "This view") {
  console.warn(err);
  return h(
    "div",
    { class: "error-box", role: "alert" },
    h("p", { class: "error-title" }, `${what} could not be shown.`),
    h("p", null, err?.message || String(err)),
    h("p", { class: "muted" }, "The rest of the site is unaffected. ",
      h("button", { type: "button", class: "btn small", onclick: () => { cache.clear(); store.phases.clear(); store.manifest = null; store.classes = null; store.theory = null; current = null; route(); } }, "Retry"),
      " or ", h("a", { href: "#/" }, "go to the overview"), "."),
  );
}

function loadingBox(text = "Loading…") {
  return h("div", { class: "loading", role: "status" }, h("span", { class: "spinner", "aria-hidden": "true" }), text);
}

function paintMockNotice() {
  const el = document.getElementById("mock-notice");
  if (!el) return;
  if (!mockFiles.size) {
    el.hidden = true;
    return;
  }
  const files = [...mockFiles].map((p) => p.replace(/^data\//, ""));
  const nGraph = files.filter((f) => f.startsWith("graphs/")).length;
  const main = files.filter((f) => !f.startsWith("graphs/"));
  const list = [...main, ...(nGraph ? [plural(nGraph, "graph file")] : [])].join(", ");
  el.hidden = false;
  el.innerHTML = "";
  el.append(
    h("strong", null, "Sample data. "),
    `These files are placeholders with the shape of the real data, not its content: ${list}. Their coefficients, rules, graphs and texts are not results.`,
  );
}

function setTitle(...parts) {
  const base = "Renormalised SM, BMHV γ₅";
  document.title = [...parts.filter(Boolean), base].join(" · ");
}

function setNav(key) {
  document.querySelectorAll("#main-nav a").forEach((a) => {
    if (a.dataset.nav === key) a.setAttribute("aria-current", "page");
    else a.removeAttribute("aria-current");
  });
}

// ------------------------------------------------------------------------------------------------ graphs

const graphObserver =
  "IntersectionObserver" in window
    ? new IntersectionObserver(
        (entries) => {
          for (const e of entries) {
            if (!e.isIntersecting) continue;
            const el = e.target;
            graphObserver.unobserve(el);
            paintGraph(el);
          }
        },
        { rootMargin: "300px 0px" },
      )
    : null;
const graphOf = new WeakMap();

function paintGraph(el) {
  const g = graphOf.get(el);
  if (!g || el.dataset.painted) return;
  el.dataset.painted = "1";
  el.innerHTML = typeof g.svg === "string" && /^\s*<svg[\s>]/.test(g.svg) ? g.svg : `<span class="muted">no drawing</span>`;
  const svg = el.querySelector("svg");
  if (svg) {
    svg.setAttribute("aria-hidden", "true");
    svg.removeAttribute("width");
    svg.removeAttribute("height");
  }
}

function viewBoxRatio(svg) {
  const m = typeof svg === "string" && svg.slice(0, 600).match(/viewBox\s*=\s*"\s*[-\d.]+[\s,]+[-\d.]+[\s,]+([\d.]+)[\s,]+([\d.]+)/);
  if (!m) return "4 / 3";
  const w = Number(m[1]), hh = Number(m[2]);
  return w > 0 && hh > 0 ? `${w} / ${hh}` : "4 / 3";
}

function fieldSymHTML(id) {
  const f = anyField(id);
  return f && hasMath(f.symbol) ? f.symbol.html : esc(id);
}

const DELTA1 = '<math class="mv-inline"><msub><mi mathvariant="normal">Δ</mi><mn>1</mn></msub></math>';
const CONTRIB = [
  ["pole", "S", "contributes to the singular part"],
  ["evanescent", "E", "contributes to an evanescent operator"],
  ["breaking", "B", "contributes to the breaking"],
];

function graphCard(g, onOpen) {
  const art = h("div", { class: "graph-svg", style: `aspect-ratio: ${viewBoxRatio(g.svg)}` });
  graphOf.set(art, g);
  if (graphObserver) graphObserver.observe(art);
  else paintGraph(art);
  const loop = arr(g.loop);
  const contrib = isObj(g.contributes) ? g.contributes : {};
  const marks = CONTRIB.filter(([k]) => contrib[k]).map(([k, s, t]) => `<span class="gmark t-${k === "pole" ? "pole" : k}" title="${esc(t)}">${s}<span class="sr-only"> ${esc(t)}</span></span>`).join("");
  const cap = h("span", { class: "graph-cap" });
  cap.innerHTML =
    `<span class="graph-loop">${loop.length ? loop.map(fieldSymHTML).join('<span class="sep"> </span>') : '<span class="muted">no loop fields</span>'}</span>` +
    `<span class="graph-meta">${g.symmetry_factor !== undefined && g.symmetry_factor !== null ? `<span title="symmetry factor">sym. ${esc(g.symmetry_factor)}</span>` : ""}` +
    `${g.insertion ? `<span class="gmark t-breaking ins" title="one insertion of the breaking Δ₁">Δ₁</span>` : ""}${marks}</span>` +
    `<span class="graph-id">${esc(g.id ?? "")}</span>`;
  const btn = h(
    "button",
    {
      type: "button", class: "graph-card",
      "aria-label": `Graph ${g.id ?? ""}${loop.length ? `, internal lines ${loop.join(" ")}` : ""}. Enlarge.`,
      onclick: () => onOpen(g),
    },
    art,
    cap,
  );
  return btn;
}

// lightbox (one <dialog> for the whole page). The frame (head with the buttons, body) is built once; a step only
// replaces the title, the drawing and the side panel, so keyboard focus stays on the button that was used.
const lightbox = {
  dlg: null, list: [], index: 0, cls: null, returnFocus: null, els: null,
  init() {
    this.dlg = document.getElementById("lightbox");
    if (!this.dlg) return;
    this.dlg.addEventListener("click", (e) => {
      if (e.target === this.dlg) this.close(); // click on the backdrop
    });
    this.dlg.addEventListener("close", () => {
      this.returnFocus?.focus?.({ preventScroll: true });
      this.returnFocus = null;
    });
    // arrow keys anywhere while the dialog is open (focus may sit on the body after a click on the drawing)
    document.addEventListener("keydown", (e) => {
      if (!this.dlg.open || e.defaultPrevented || e.altKey || e.ctrlKey || e.metaKey) return;
      if (e.target?.matches?.("input, select, textarea")) return;
      if (e.key === "ArrowRight") { this.step(1); e.preventDefault(); }
      else if (e.key === "ArrowLeft") { this.step(-1); e.preventDefault(); }
    });
    const btn = (cls, label, text, fn) => h("button", { type: "button", class: `btn ${cls}`, "aria-label": label, onclick: fn }, text);
    const title = h("h2", { id: "lb-title" });
    const prev = btn("lb-prev", "Previous graph", "‹ Previous", () => this.step(-1));
    const next = btn("lb-next", "Next graph", "Next ›", () => this.step(1));
    const close = btn("lb-close", "Close", "Close", () => this.close());
    const art = h("div", { class: "lb-art" });
    const side = h("div", { class: "lb-side" });
    this.dlg.replaceChildren(
      h("div", { class: "lb-head" }, title, h("div", { class: "lb-actions" }, prev, next, close)),
      h("div", { class: "lb-body" }, art, side),
    );
    this.dlg.setAttribute("aria-labelledby", "lb-title");
    this.els = { title, prev, next, close, art, side };
  },
  open(list, index, cls) {
    if (!this.dlg) return;
    this.list = list;
    this.index = index;
    this.cls = cls;
    if (!this.dlg.open) this.returnFocus = document.activeElement;
    this.paint();
    if (!this.dlg.open) {
      if (typeof this.dlg.showModal === "function") this.dlg.showModal();
      else this.dlg.setAttribute("open", "");
    }
    this.els.close.focus();
  },
  close() {
    if (!this.dlg?.open) return;
    if (typeof this.dlg.close === "function") this.dlg.close();
    else this.dlg.removeAttribute("open");
  },
  step(d) {
    if (this.list.length < 2) return;
    this.index = (this.index + d + this.list.length) % this.list.length;
    this.paint();
  },
  paint() {
    const g = this.list[this.index];
    const E = this.els;
    const cls = this.cls;
    const contrib = isObj(g.contributes) ? g.contributes : {};
    const signText = (v) => (String(v) === "-1" ? "−1" : String(v) === "1" ? "+1" : String(v));
    const has = (v) => v !== undefined && v !== null;
    const rows = [
      ["Graph", h("code", null, g.id ?? "")],
      ["Green's function", cls ? h("span", { html: classLabelHTML(cls) }) : ""],
      ["Topology", g.topology ?? null],
      ["Internal lines", arr(g.loop).length ? h("span", { class: "lb-loop" }, arr(g.loop).map((f) => fieldChip(null, f))) : null],
      ["Symmetry factor", g.symmetry_factor ?? null],
      ["Closed-loop sign", has(g.sign) ? `${signText(g.sign)}, as in the report's tables` : null],
      ["Wick sign", has(g.wick_sign) ? `${signText(g.wick_sign)}, in the package's leg order (a different quantity)` : null],
      ["Insertion", g.insertion ? (g.insertion === "Delta1" ? h("span", { html: `one insertion of ${DELTA1}` }) : g.insertion) : "none"],
      [
        "Contributes to",
        CONTRIB.filter(([k]) => contrib[k]).length
          ? h("span", { class: "lb-contrib" }, CONTRIB.filter(([k]) => contrib[k]).map(([k, , t]) => h("span", { class: `tag t-${k}` }, h("span", { class: "tag-mark", "aria-hidden": "true" }), t.replace(/^contributes to /, ""))))
          : "—",
      ],
    ];
    E.title.textContent = `Graph ${this.index + 1} of ${this.list.length}`;
    E.prev.disabled = E.next.disabled = this.list.length < 2;
    E.art.innerHTML = typeof g.svg === "string" ? g.svg : "";
    const svg = E.art.querySelector("svg");
    if (svg) {
      svg.removeAttribute("width");
      svg.removeAttribute("height");
      svg.setAttribute("role", "img");
      svg.setAttribute("aria-label", `Drawing of graph ${g.id ?? ""}`);
    }
    const lg = graphLegend();
    E.side.replaceChildren(
      h("dl", { class: "lb-meta" }, rows.filter(([, v]) => v !== null && v !== "").map(([k, v]) => [h("dt", null, k), h("dd", null, v)])),
      g.caption ? rich(g.caption, { cls: "rich lb-caption" }) : "",
      lg ? h("details", { class: "lb-legend" }, h("summary", null, "Line styles of the drawings"), lg) : "",
      h("p", { class: "muted small" }, "Arrow keys step through the graphs of the current selection; Escape closes."),
    );
  },
};

// A class is named by its label (math object); a build without labels falls back to the plain text or the id.
function classLabelHTML(cls) {
  if (hasMath(cls?.label)) return cls.label.html;
  return esc(typeof cls?.text === "string" && cls.text ? cls.text : cls?.id ?? "");
}
function classLabelText(cls) {
  return cls?.label?.text || cls?.text || cls?.id || "";
}

// Legend of the line styles of the drawings, when the build provides one (classes.json "legend", optional).
function graphLegend() {
  const leg = arr(store.classes?.data?.legend).filter((x) => isObj(x) && typeof x.svg === "string" && /^\s*<svg[\s>]/.test(x.svg));
  if (!leg.length) return null;
  return h("ul", { class: "graph-legend", "aria-label": "Line styles of the drawings" },
    leg.map((x) => {
      const art = h("span", { class: "gl-art", "aria-hidden": "true", html: x.svg });
      return h("li", null, art, typeof x.html === "string" && x.html ? rich(x.html, { tag: "span", cls: "rich" }) : h("span", null, x.label || x.kind || ""));
    }));
}

// Collapsible block of one class of Green's functions, with lazy loading of its graphs.
function classBlock(cls, { open = false, backlinks = false, via = null, idPrefix = "cls" } = {}) {
  const det = h("details", { class: `class-block${cls.insertion ? " has-insertion" : ""}`, id: `${idPrefix}-${cls.id.replace(/[^A-Za-z0-9_-]/g, "_")}` });
  const sum = h("summary", { class: "class-sum" });
  const label = h("span", { class: "class-label" });
  label.innerHTML = classLabelHTML(cls);
  const group = store.classes?.groupByKey.get(cls.group);
  sum.append(
    h("span", { class: "class-caret", "aria-hidden": "true" }),
    label,
    h("span", { class: "class-info" },
      group && !backlinks ? h("span", { class: "class-group" }, group.title) : null,
      cls.insertion ? h("span", { class: "gmark t-breaking ins", title: "one insertion of the breaking Δ₁" }, "Δ₁") : null,
      h("span", { class: "class-count" }, Number.isFinite(cls.n_graphs) ? plural(cls.n_graphs, "graph") : "graphs")),
  );
  det.append(sum);
  const body = h("div", { class: "class-body" });
  det.append(body);
  let built = false;
  const build = () => {
    if (built) return;
    built = true;
    if (via && via.length) body.append(h("p", { class: "class-via" }, "Through ", via.map((v, i) => [i ? ", " : "", v])));
    if (cls.summary) body.append(rich(cls.summary, { cls: "rich class-summary" }));
    const res = isObj(cls.result) ? cls.result : {};
    const isZero = (m) => hasMath(m) && String(m.text ?? "").trim() === "0";
    const resRows = [];
    if (hasMath(res.pole) && !isZero(res.pole)) resRows.push(["pole", "Singular part", res.pole]);
    if (hasMath(res.breaking) && !isZero(res.breaking)) resRows.push(["breaking", "Breaking", res.breaking]);
    const zeros = [isZero(res.pole) ? "no pole" : null, isZero(res.breaking) ? "no finite breaking" : null].filter(Boolean);
    if (zeros.length) body.append(h("p", { class: "class-zero" }, h("strong", null, "Result: zero"), ` (${zeros.join(", ")}).`));
    if (resRows.length)
      body.append(
        h("div", { class: "class-results" },
          resRows.map(([k, l, m]) => h("div", { class: `coef-block t-${k}` }, h("div", { class: "coef-head" }, typeTag(k), h("span", { class: "coef-title" }, `Result: ${l.toLowerCase()}`)), mathBox(m)))),
      );
    const meta = [h("code", { class: "class-id" }, cls.id), cls.fields ? h("span", { class: "class-fields" }, arr(cls.fields).map((f) => fieldChip(null, f))) : null];
    if (group && !backlinks) {
      meta.push(h("a", { href: `#/classes/${enc(cls.id)}`, class: "small" }, "Open in the Green's-function browser"));
      if (group.report) meta.push(reportRefs(group.report, { cls: "refs inline" }));
    }
    body.append(h("div", { class: "class-meta" }, meta));
    if (backlinks) {
      const refs = termsOfClass(cls.id);
      body.append(
        h("div", { class: "backlinks" },
          h("h4", null, "Terms determined with this class"),
          refs.length
            ? h("ul", { class: "backlink-list" }, refs.map(({ P, t, roles }) =>
                h("li", null, nameEl("a", t, { href: `#/phase/${enc(P.id)}/term/${enc(t.id)}`, title: t.name }), h("span", { class: "muted small" }, ` ${P.title}, ${roles.map((r) => (r === "pole" ? "singular part" : "finite part through the breaking")).join(" and ")}`))))
            : h("p", { class: "muted" }, "No operator refers to this class.")),
      );
    }
    const gal = h("div", { class: "gallery" });
    body.append(gal);
    buildGallery(gal, cls);
  };
  det.addEventListener("toggle", () => det.open && build());
  if (open) {
    det.open = true;
    build();
  }
  return det;
}

async function buildGallery(root, cls) {
  if (!cls.graphs) {
    root.append(h("p", { class: "muted" }, "No graph file is available for this class."));
    return;
  }
  const wait = loadingBox("Loading the graphs…");
  root.append(wait);
  let data;
  try {
    data = await loadJSON(cls.graphs);
  } catch (e) {
    wait.replaceWith(errorBox(e, "The graphs of this class"));
    return;
  }
  wait.remove();
  const graphs = arr(data?.graphs).filter(isObj);
  if (!graphs.length) {
    root.append(h("p", { class: "muted" }, "The graph file of this class is empty."));
    return;
  }
  const loopCount = new Map();
  for (const g of graphs) for (const f of new Set(arr(g.loop))) loopCount.set(f, (loopCount.get(f) || 0) + 1);
  const contribCount = Object.fromEntries(CONTRIB.map(([k]) => [k, graphs.filter((g) => g.contributes?.[k]).length]));
  const sel = new Set();
  const selC = new Set();
  let page = 0;
  const bar = h("div", { class: "gallery-bar" });
  const grid = h("div", { class: "graph-grid", role: "list" });
  const pager = h("nav", { class: "pager", "aria-label": "Graph pages" });
  const status = h("p", { class: "gallery-status", "aria-live": "polite" });
  root.append(bar, status, grid, pager);

  const filtered = () =>
    graphs.filter((g) => {
      const loop = arr(g.loop);
      for (const f of sel) if (!loop.includes(f)) return false;
      for (const k of selC) if (!g.contributes?.[k]) return false;
      return true;
    });

  if (loopCount.size > 1) {
    const chips = [...loopCount.entries()].sort((a, b) => kindRank(a[0]) - kindRank(b[0]) || b[1] - a[1]).map(([f, n]) => {
      const b = h("button", { type: "button", class: "chip chip-field", "aria-pressed": "false", title: anyField(f)?.name || f });
      b.innerHTML = `<span class="chip-sym">${fieldSymHTML(f)}</span><span class="chip-n">${n}</span><span class="sr-only"> ${esc(anyField(f)?.name || f)}, ${n} graphs</span>`;
      b.addEventListener("click", () => {
        sel.has(f) ? sel.delete(f) : sel.add(f);
        b.setAttribute("aria-pressed", String(sel.has(f)));
        page = 0;
        paint();
      });
      return b;
    });
    bar.append(h("div", { class: "gallery-filter", role: "group", "aria-label": "Internal lines include" }, h("span", { class: "flabel" }, "Internal lines include"), chips));
  }
  const cchips = CONTRIB.filter(([k]) => contribCount[k] > 0 && contribCount[k] < graphs.length).map(([k, s, t]) => {
    const b = h("button", { type: "button", class: `chip chip-type t-${k}`, "aria-pressed": "false" }, h("span", { class: "tag-mark", "aria-hidden": "true" }), t.replace(/^contributes to (the |an )?/, ""), h("span", { class: "chip-n" }, String(contribCount[k])));
    b.addEventListener("click", () => {
      selC.has(k) ? selC.delete(k) : selC.add(k);
      b.setAttribute("aria-pressed", String(selC.has(k)));
      page = 0;
      paint();
    });
    return b;
  });
  if (cchips.length) bar.append(h("div", { class: "gallery-filter", role: "group", "aria-label": "Contributes to" }, h("span", { class: "flabel" }, "Contributes to"), cchips));

  function paint() {
    const list = filtered();
    const pages = Math.max(1, Math.ceil(list.length / GRAPH_PAGE));
    page = Math.min(page, pages - 1);
    const from = page * GRAPH_PAGE;
    const shown = list.slice(from, from + GRAPH_PAGE);
    grid.querySelectorAll(".graph-svg").forEach((el) => graphObserver?.unobserve(el));
    grid.replaceChildren(...shown.map((g, i) => h("div", { role: "listitem", class: "graph-item" }, graphCard(g, () => lightbox.open(list, from + i, cls)))));
    status.textContent = list.length
      ? `${list.length === graphs.length ? "" : `${plural(list.length, "graph")} of ${graphs.length} selected. `}Showing ${from + 1}–${from + shown.length}${pages > 1 ? ` (page ${page + 1} of ${pages})` : ""}.`
      : "No graph matches the selection.";
    pager.replaceChildren();
    if (pages > 1) {
      const go = (p) => {
        page = p;
        paint();
        root.scrollIntoView({ block: "nearest", behavior: smooth() });
      };
      pager.append(h("button", { type: "button", class: "btn", disabled: page === 0, onclick: () => go(page - 1) }, "‹ Previous"));
      const select = h("select", { "aria-label": "Page", class: "pager-select" });
      for (let p = 0; p < pages; p++) select.append(h("option", { value: p, selected: p === page }, `${p * GRAPH_PAGE + 1}–${Math.min(list.length, (p + 1) * GRAPH_PAGE)}`));
      select.addEventListener("change", () => go(Number(select.value)));
      pager.append(select, h("span", { class: "muted small" }, `of ${list.length}`));
      pager.append(h("button", { type: "button", class: "btn", disabled: page >= pages - 1, onclick: () => go(page + 1) }, "Next ›"));
    }
  }
  paint();
}

function kindRank(fid) {
  const k = anyField(fid)?.kind;
  const i = KIND_ORDER.indexOf(k);
  return i < 0 ? 99 : i;
}

// ------------------------------------------------------------------------------------------------ router

const main = () => document.getElementById("main");
let current = null; // { key, update(route), destroy() }

function parseRoute() {
  let raw = location.hash.replace(/^#/, "");
  if (!raw.startsWith("/")) raw = "/" + raw;
  const qi = raw.indexOf("?");
  const path = qi >= 0 ? raw.slice(0, qi) : raw;
  const query = new URLSearchParams(qi >= 0 ? raw.slice(qi + 1) : "");
  const parts = path.split("/").filter(Boolean).map((p) => {
    try {
      return decodeURIComponent(p);
    } catch {
      return p;
    }
  });
  return { parts, query, raw };
}

async function route() {
  lightbox.close();
  const r = parseRoute();
  const [head, a] = r.parts;
  let key, factory;
  if (!head) [key, factory] = ["home", () => HomeView()];
  else if (head === "phase" && a) [key, factory] = [`phase:${a}`, () => PhaseView(a)];
  else if (head === "classes") [key, factory] = ["classes", () => ClassesView()];
  else if (head === "theory") [key, factory] = ["theory", () => TheoryView()];
  else if (head === "report") [key, factory] = ["report", () => ReportView()];
  else if (head === "resources") [key, factory] = ["resources", () => ResourcesView()];
  else [key, factory] = ["missing", () => MissingView()];
  if (current && current.key === key) {
    current.update(r);
    return;
  }
  current?.destroy?.();
  const root = main();
  viewerModule?.unmount?.(root);
  root.replaceChildren();
  window.scrollTo(0, 0);
  const view = factory();
  current = { key, ...view };
  setNav(key.startsWith("phase:") ? key : key);
  try {
    await view.mount(root, r);
  } catch (e) {
    root.replaceChildren(h("div", { class: "page" }, errorBox(e)));
  }
}

// scroll positions of term details (phase:id), restored when a term is reopened through the history; a term
// opened from the table (navigate) starts at the top
const DETAIL_SCROLL = new Map();
let freshHash = null;

function navigate(hash, { replace = false } = {}) {
  if (!replace) freshHash = hash;
  if (replace) {
    history.replaceState(null, "", hash);
  } else if (location.hash !== hash) {
    location.hash = hash;
  }
}

// ------------------------------------------------------------------------------------------------ home

function HomeView() {
  return {
    async mount(root) {
      setTitle();
      const man = await getManifest();
      const wrap = h("div", { class: "page home" });
      root.append(wrap);
      // the views first, then the master formula; the statement of what is generated and what is taken from the
      // report is one disclosure below them
      // the first paragraph of the status (what the page is) stays visible; the rest folds
      const status = man.status ? rich(man.status, { cls: "rich status" }) : null;
      const statusLead = status && status.firstElementChild ? rich("", { cls: "rich status-lead" }) : null;
      if (statusLead) statusLead.append(status.firstElementChild);
      const statusBox = status && status.childElementCount
        ? h("details", { class: "status-box" },
            h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), "What is generated by the automation and what is taken from the report"),
            status)
        : null;
      wrap.append(
        h("header", { class: "hero" },
          h("p", { class: "eyebrow" }, `Renormalised Lagrangian · ${!man.loop_order || man.loop_order === 1 ? "one loop" : `${man.loop_order} loops`}`),
          h("h1", { tabindex: "-1" }, man.title || "Renormalised Standard Model"),
          man.subtitle ? h("p", { class: "lede" }, man.subtitle) : null,
          statusLead),
      );
      const tiles = h("section", { class: "tiles", "aria-label": "Views" });
      const master = h("section", { class: "master", "aria-labelledby": "master-h" });
      const legend = h("section", { class: "legend-sec", "aria-labelledby": "legend-h" });
      const conv = h("section", { class: "conventions", "aria-labelledby": "conv-h" });
      wrap.append(tiles, master, statusBox || "", legend, conv);

      const phases = phaseEntries(man);
      let first = null;
      try {
        const ub = phases.find((p) => p.id === "unbroken") || phases[0];
        first = ub ? await getPhase(ub.id) : null;
      } catch (e) {
        master.append(errorBox(e, "The master formula"));
      }
      await viewerReady();
      const ov = arr(first?.data.overview).filter((o) => isObj(o) && hasMath(o.math));
      if (ov.length) {
        const o = ov[0];
        appendKids(master, [h("h2", { id: "master-h" }, o.title || "Renormalised action"), mathBox(o.math, { cls: "mbox-hero", label: o.title }), o.text ? rich(o.text, { cls: "rich master-text" }) : null]);
        if (ov.length > 1)
          master.append(h("p", { class: "small" }, h("a", { href: `#/phase/${enc(first.id)}` }, `Its ${ov.length - 1} parts and every operator are listed in the ${first.title.toLowerCase()} ›`)));
      } else if (first) master.remove();

      // tiles
      const tile = (href, title, text, stats) =>
        h("a", { class: "tile", href }, h("span", { class: "tile-title" }, title), h("span", { class: "tile-text" }, text), stats ? h("span", { class: "tile-stats" }, stats) : null);
      for (const e of phases) {
        let P = null;
        try {
          P = await getPhase(e.id);
        } catch {
          /* tile without numbers */
        }
        const counts = P ? TYPE_ORDER.map((t) => [t, P.terms.filter((x) => x._types.includes(t)).length]).filter(([, n]) => n) : [];
        tiles.append(
          tile(`#/phase/${enc(e.id)}`, e.title || P?.title || e.id, P ? `${plural(P.terms.length, "operator")} of the renormalised action, with filters, coefficients, Feynman rules and the graphs behind each.` : "Operator table.",
            counts.length ? h("span", { class: "tile-tags" }, counts.map(([t, n]) => h("span", { class: `tag t-${t}` }, h("span", { class: "tag-mark", "aria-hidden": "true" }), `${P.typeLabels[t] || TYPE_SHORT[t]} ${n}`))) : null),
        );
      }
      let C = null;
      try {
        C = await getClasses();
      } catch {
        /* tile without numbers */
      }
      const nG = C ? C.classes.reduce((s, c) => s + (Number.isFinite(c.n_graphs) ? c.n_graphs : 0), 0) : 0;
      tiles.append(tile("#/classes", "Green's functions", C ? `${plural(C.groups.length, "group")} of 1PI functions, ${plural(C.classes.length, "class", "classes")} with definite external fields, ${plural(nG, "one-loop graph")}.` : "The classes of Green's functions and their graphs."));
      let T = null;
      try {
        T = await getTheory();
      } catch {
        /* tile without numbers */
      }
      tiles.append(tile("#/theory", "Theory", T ? `${plural(Object.keys(T.cards).length, "card")} explaining each part of the action and how its coefficients are fixed.` : "Theory cards."));
      tiles.append(tile("#/report", "Report", `The full report${man.report_pages ? `, ${man.report_pages} pages` : ""}, opened at the page of every reference followed on this site.`));
      tiles.append(tile("#/resources", "Resources", "Model files, order store, run report and final numbers of the automation, for download."));

      // legend
      const labels = first?.typeLabels || TYPE_SHORT;
      legend.append(
        h("h2", { id: "legend-h" }, "Types of terms"),
        h("dl", { class: "legend" }, TYPE_ORDER.map((t) => {
          const def = typeof man.types?.[t] === "string" ? rich(man.types[t], { tag: "dd", cls: "rich" }) : h("dd", null, TYPE_DEF[t]);
          return h("div", { class: `legend-row t-${t}` }, h("dt", null, typeTag(t, labels)), def);
        })),
      );
      if (typeof man.conventions === "string" && man.conventions.trim()) conv.append(h("h2", { id: "conv-h" }, "Conventions"), rich(man.conventions, { cls: "rich conv-list" }));
      else conv.remove();
      wrap.querySelector("h1")?.focus({ preventScroll: true });
    },
    update() {},
  };
}

// ------------------------------------------------------------------------------------------------ phase view

function readFilters(q) {
  const set = (k) => new Set(q.getAll(k).flatMap((v) => v.split("~")).filter(Boolean));
  return {
    q: q.get("q") || "",
    fields: set("f"),
    mode: q.get("fm") === "any" ? "any" : "all",
    types: set("t"),
    sectors: set("s"),
    groups: set("g"),
    kinds: set("k"),
    legs: q.get("l") || "",
    ghost: q.get("gh") || "",
    sort: q.get("sort") || "",
    dir: q.get("dir") === "desc" ? "desc" : "asc",
  };
}

function writeFilters(F) {
  const q = new URLSearchParams();
  if (F.q) q.set("q", F.q);
  for (const f of F.fields) q.append("f", f);
  if (F.fields.size > 1 && F.mode === "any") q.set("fm", "any");
  for (const t of F.types) q.append("t", t);
  for (const s of F.sectors) q.append("s", s);
  for (const g of F.groups) q.append("g", g);
  for (const k of F.kinds) q.append("k", k);
  if (F.legs) q.set("l", F.legs);
  if (F.ghost) q.set("gh", F.ghost);
  if (F.sort) {
    q.set("sort", F.sort);
    if (F.dir === "desc") q.set("dir", "desc");
  }
  return q.toString();
}

function matches(t, F) {
  if (F.q) {
    for (const tok of F.q.toLowerCase().split(/\s+/).filter(Boolean)) if (!t._search.includes(tok)) return false;
  }
  if (F.fields.size) {
    const fl = arr(t.fields);
    if (F.mode === "all") {
      for (const f of F.fields) if (!fl.includes(f)) return false;
    } else if (![...F.fields].some((f) => fl.includes(f))) return false;
  }
  if (F.types.size && !t._types.some((x) => F.types.has(x))) return false;
  if (F.sectors.size && !F.sectors.has(t.sector)) return false;
  if (F.groups.size && !arr(t.groups).some((g) => F.groups.has(g))) return false;
  if (F.kinds.size && !F.kinds.has(t.kind)) return false;
  if (F.legs && String(t.legs) !== F.legs) return false;
  if (F.ghost && String(t.ghost_number) !== F.ghost) return false;
  return true;
}

const SORTERS = {
  name: (a, b) => String(a.name || a.id).localeCompare(String(b.name || b.id), "en"),
  fields: (a, b) => (a.legs ?? 0) - (b.legs ?? 0) || arr(a.fields).join(" ").localeCompare(arr(b.fields).join(" ")),
  tree: coefSorter(["tree"]),
  pole: coefSorter(["pole", "evanescent"]),
  finite: coefSorter(["finite"]),
  types: (a, b) => TYPE_ORDER.indexOf(a._types[0]) - TYPE_ORDER.indexOf(b._types[0]) || a._types.length - b._types.length,
};

// columns whose empty cells sort last in either direction
const coefValue = (keys) => (t) => {
  const c = isObj(t.coefficients) ? t.coefficients : {};
  return !keys.some((k) => hasMath(c[k]?.math));
};
const EMPTY_LAST = { tree: coefValue(["tree"]), pole: coefValue(["pole", "evanescent"]), finite: coefValue(["finite"]) };

function coefSorter(keys) {
  const pick = (t) => {
    const c = isObj(t.coefficients) ? t.coefficients : {};
    for (const k of keys) if (hasMath(c[k]?.math)) return c[k].math.text || "~";
    return null;
  };
  return (a, b) => {
    const x = pick(a), y = pick(b);
    if (x === null && y === null) return 0;
    if (x === null) return 1;
    if (y === null) return -1;
    return x.localeCompare(y, "en", { numeric: true });
  };
}

const COLUMNS = [
  { key: "name", label: "Operator", cls: "col-op" },
  { key: "fields", label: "Fields", cls: "col-fields" },
  { key: "tree", label: "Tree level", cls: "col-coef col-tree", type: "tree" },
  { key: "pole", label: "Singular (1/ϵ)", cls: "col-coef col-pole", type: "pole" },
  { key: "finite", label: "Finite", cls: "col-coef col-finite", type: "finite" },
  { key: "types", label: "Type", cls: "col-types" },
];

function PhaseView(pid) {
  let P, F, root, rows, list = [], els = {}, detailId = null, savedScroll = 0, searchTimer = 0;
  const wide = window.matchMedia(WIDE);
  const phone = window.matchMedia("(max-width: 760px)");
  const onWide = () => layout();
  const onKey = (e) => {
    if (e.key === "Escape" && detailId && !lightbox.dlg?.open && !e.defaultPrevented) {
      const tgt = e.target;
      if (tgt?.matches?.("input, select, textarea") && !tgt.closest(".detail")) return;
      e.preventDefault();
      closeDetail();
    }
  };

  const hrefTable = () => {
    const qs = writeFilters(F);
    return `#/phase/${enc(pid)}${qs ? "?" + qs : ""}`;
  };
  const hrefTerm = (id) => {
    const qs = writeFilters(F);
    return `#/phase/${enc(pid)}/term/${enc(id)}${qs ? "?" + qs : ""}`;
  };

  function closeDetail() {
    navigate(hrefTable());
  }

  function layout() {
    const open = !!detailId;
    els.explorer.classList.toggle("has-detail", open);
    els.explorer.classList.toggle("detail-full", open && !wide.matches);
    root.classList.toggle("detail-full-page", open && !wide.matches);
  }

  async function mount(r, rr) {
    root = h("div", { class: "page phase" });
    r.append(root);
    root.append(loadingBox());
    P = await getPhase(pid);
    await viewerReady();
    root.replaceChildren();
    setTitle(P.title);
    F = readFilters(rr.query);
    const d = P.data;
    // the table comes first: the introduction shows its first paragraph, the rest and the master formulae fold
    const intro = d.intro ? rich(d.intro, { cls: "rich intro" }) : null;
    let introMore = null;
    if (intro) {
      const paras = [...intro.children];
      if (paras.length > 1) {
        const rest = rich("", { cls: "rich intro" });
        rest.append(...paras.slice(1));
        introMore = h("details", { class: "intro-more" }, h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), "More about this phase"), rest);
      }
    }
    const head = h("header", { class: "page-head" }, h("p", { class: "eyebrow" }, "Renormalised Lagrangian"), h("h1", { tabindex: "-1" }, P.title), intro, introMore);
    const over = arr(d.overview).filter(isObj);
    let ovSec = null;
    if (over.length) {
      const list = h("div", { class: "overview" },
        over.map((o) => {
          const det = h("details", { class: "ov-card" });
          const body = h("div", { class: "ov-body" });
          det.append(h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), h("span", { class: "ov-title" }, o.title || o.id || "Formula")), body);
          let done = false;
          det.addEventListener("toggle", () => {
            if (!det.open || done) return;
            done = true;
            if (hasMath(o.math)) body.append(mathBox(o.math, { cls: "mbox-ov", label: o.title }));
            if (o.text) body.append(rich(o.text, { cls: "rich ov-text" }));
          });
          return det;
        }));
      ovSec = h("details", { class: "ov-block", "aria-label": "Master formulae" },
        h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), h("span", { class: "ov-block-title" }, `Master formulae (${over.length})`),
          h("span", { class: "muted small" }, over.map((o) => o.title).filter(Boolean).slice(0, 3).join(" · ") + (over.length > 3 ? " · …" : ""))),
        list);
    }
    els.explorer = h("div", { class: "explorer" });
    els.mainCol = h("div", { class: "explorer-main" });
    els.detail = h("aside", { class: "detail", "aria-label": "Term detail", tabindex: "-1" });
    els.explorer.append(els.mainCol, els.detail);
    appendKids(root, [head, ovSec, els.explorer]);
    buildFilters();
    buildTable();
    apply(rr, { initial: true });
    wide.addEventListener("change", onWide);
    phone.addEventListener("change", onWide);
    document.addEventListener("keydown", onKey);
    if (!detailId) head.querySelector("h1")?.focus({ preventScroll: true });
  }

  // ---------------------------------------------------------------- filter bar
  function chipGroup(label, items, selected, { cls = "", onToggle, render }) {
    const box = h("div", { class: `fgroup ${cls}`.trim(), role: "group", "aria-label": label });
    box.append(h("span", { class: "flabel" }, label));
    const wrap = h("div", { class: "chips" });
    for (const it of items) {
      const b = h("button", { type: "button", class: `chip ${it.cls || ""}`.trim(), "aria-pressed": String(selected.has(it.id)), dataset: { id: it.id }, title: it.title || null });
      if (render) b.innerHTML = render(it);
      else b.textContent = it.label;
      b.addEventListener("click", () => {
        onToggle(it.id);
      });
      wrap.append(b);
    }
    box.append(wrap);
    return box;
  }

  function buildFilters() {
    const form = h("form", { class: "filters", role: "search", "aria-label": `Filter the operators of the ${P.title.toLowerCase()}`, onsubmit: (e) => e.preventDefault() });
    els.filters = form;
    // row 1: search and selects
    const q = h("input", { type: "search", id: `q-${pid}`, class: "search", placeholder: "Name, field or id", value: F.q, autocomplete: "off", spellcheck: "false" });
    q.addEventListener("input", () => {
      clearTimeout(searchTimer);
      searchTimer = setTimeout(() => {
        F.q = q.value.trim();
        refresh();
      }, 120);
    });
    q.addEventListener("keydown", (e) => {
      if (e.key === "Escape" && q.value) {
        q.value = "";
        F.q = "";
        refresh();
        e.preventDefault();
      }
    });
    els.q = q;
    const sel = (id, label, values, cur, fmt, on) => {
      const s = h("select", { id: `${id}-${pid}` }, h("option", { value: "" }, "any"), values.map((v) => h("option", { value: String(v), selected: String(v) === cur }, fmt(v))));
      s.addEventListener("change", () => on(s.value));
      return h("label", { class: "fselect", for: `${id}-${pid}` }, h("span", { class: "flabel" }, label), s);
    };
    const legsVals = [...new Set(P.terms.map((t) => t.legs).filter(Number.isFinite))].sort((a, b) => a - b);
    const ghVals = [...new Set(P.terms.map((t) => t.ghost_number).filter(Number.isFinite))].sort((a, b) => a - b);
    els.legs = sel("legs", "Legs", legsVals, F.legs, (v) => String(v), (v) => { F.legs = v; refresh(); });
    els.ghost = sel("gh", "Ghost number", ghVals, F.ghost, (v) => (v > 0 ? `+${v}` : String(v)).replace("-", "−"), (v) => { F.ghost = v; refresh(); });
    const reset = h("button", { type: "button", class: "btn btn-quiet", onclick: () => resetFilters() }, "Reset filters");
    els.skip = h("a", { class: "skip-results", href: "#", onclick: (e) => {
      e.preventDefault();
      const row = els.tbody.querySelector("tr.term-row");
      if (row) {
        row.focus({ preventScroll: true });
        row.scrollIntoView({ block: "center", behavior: smooth() });
      } else els.result.focus?.();
    } }, "Skip to the results");
    form.append(els.skip);
    form.append(
      h("div", { class: "frow frow-top" },
        h("label", { class: "fsearch", for: `q-${pid}` }, h("span", { class: "flabel" }, "Search"), q),
        els.legs, els.ghost, reset),
    );
    // rows of chips
    const toggle = (setName) => (id) => {
      F[setName].has(id) ? F[setName].delete(id) : F[setName].add(id);
      refresh();
    };
    const typeItems = P.facets.type.length ? P.facets.type.map((x) => ({ id: x.id, label: x.label })) : TYPE_ORDER.map((t) => ({ id: t, label: TYPE_SHORT[t] }));
    const typeRow = chipGroup("Type", typeItems, F.types, { cls: "fg-type", onToggle: toggle("types"), render: (it) => `<span class="tag-mark" aria-hidden="true"></span>${esc(it.label)}` });
    for (const b of typeRow.querySelectorAll(".chip")) b.classList.add("chip-type", `t-${b.dataset.id}`);
    form.append(h("div", { class: "frows frows-type" }, typeRow));
    els.chipRows = h("div", { class: "frows" });
    if (P.facets.sector.length) els.chipRows.append(chipGroup("Sector", P.facets.sector, F.sectors, { onToggle: toggle("sectors") }));
    if (P.facets.group.length) els.chipRows.append(chipGroup("Gauge group", P.facets.group, F.groups, { onToggle: toggle("groups") }));
    if (P.facets.kind.length) els.chipRows.append(chipGroup("Kind", P.facets.kind, F.kinds, { onToggle: toggle("kinds") }));
    // sector, group, kind and fields fold into one disclosure (open when one of them is in use)
    const anyFacet = F.sectors.size || F.groups.size || F.kinds.size || F.fields.size;
    els.more = h("details", { class: "fmore", open: !!anyFacet });
    els.moreSum = h("span", { class: "muted small" });
    els.more.append(h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), h("span", { class: "flabel" }, P.facets.kind.length ? "More filters: sector, group, kind and fields" : "More filters: sector, group and fields"), els.moreSum));
    const moreBody = h("div", { class: "fmore-body" });
    els.more.append(moreBody);
    form.append(els.more);
    moreBody.append(els.chipRows);
    // fields
    const used = new Set(P.terms.flatMap((t) => arr(t.fields)));
    const byKind = new Map();
    for (const f of P.fields) {
      if (!used.has(f.id)) continue;
      const k = KIND_ORDER.includes(f.kind) ? f.kind : "other";
      if (!byKind.has(k)) byKind.set(k, []);
      byKind.get(k).push(f);
    }
    const fdet = h("details", { class: "ffields", open: window.matchMedia("(min-width: 761px)").matches || F.fields.size > 0 });
    els.fieldsSummary = h("span", { class: "ffields-sum" });
    const modeBtn = (m, label) => {
      const b = h("button", { type: "button", class: "seg-btn", "aria-pressed": String(F.mode === m), dataset: { mode: m } }, label);
      b.addEventListener("click", () => {
        F.mode = m;
        refresh();
      });
      return b;
    };
    els.mode = h("div", { class: "seg", role: "group", "aria-label": "Combine the selected fields" }, modeBtn("all", "contains all selected"), modeBtn("any", "contains any"));
    fdet.append(h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), h("span", { class: "flabel" }, "Fields"), els.fieldsSummary));
    const fbody = h("div", { class: "ffields-body" }, h("div", { class: "fmode" }, els.mode));
    for (const k of [...KIND_ORDER, "other"]) {
      const fl = byKind.get(k);
      if (!fl) continue;
      fbody.append(
        chipGroup(KIND_LABEL[k] || "Other", fl.map((f) => ({ id: f.id, label: f.id, title: f.name, cls: "chip-field" })), F.fields, {
          cls: "fg-fields",
          onToggle: toggle("fields"),
          render: (it) => {
            const f = P.fieldById.get(it.id);
            return `<span class="chip-sym">${hasMath(f?.symbol) ? f.symbol.html : esc(it.id)}</span><span class="sr-only">${esc(f?.name || it.id)}</span>`;
          },
        }),
      );
    }
    fdet.append(fbody);
    moreBody.append(fdet);
    // result bar
    els.result = h("div", { class: "result-bar", "aria-live": "polite" });
    els.sortSel = h("select", { id: `sort-${pid}` },
      h("option", { value: "" }, "Table order"),
      COLUMNS.map((c) => h("option", { value: c.key }, c.key === "name" ? "Name" : c.key === "fields" ? "Number of legs" : c.label)));
    els.sortSel.addEventListener("change", () => {
      F.sort = els.sortSel.value;
      F.dir = "asc";
      refresh();
    });
    els.mainCol.append(form, h("div", { class: "result-row" }, els.result, h("label", { class: "mobile-sort", for: `sort-${pid}` }, "Sort", els.sortSel)));
  }

  function syncControls() {
    if (els.q.value.trim() !== F.q) els.q.value = F.q;
    if (els.sortSel) els.sortSel.value = SORTERS[F.sort] ? F.sort : "";
    els.legs.querySelector("select").value = F.legs;
    els.ghost.querySelector("select").value = F.ghost;
    els.filters.querySelectorAll(".fgroup").forEach((g) => {
      const label = g.getAttribute("aria-label");
      const set = g.classList.contains("fg-type") ? F.types : g.classList.contains("fg-fields") ? F.fields : label === "Sector" ? F.sectors : label === "Gauge group" ? F.groups : label === "Kind" ? F.kinds : null;
      if (!set) return;
      g.querySelectorAll(".chip").forEach((b) => b.setAttribute("aria-pressed", String(set.has(b.dataset.id))));
    });
    els.mode.querySelectorAll(".seg-btn").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.mode === F.mode)));
    els.mode.classList.toggle("is-idle", F.fields.size < 2);
    const nf = F.sectors.size + F.groups.size + F.kinds.size + F.fields.size;
    els.moreSum.textContent = nf ? `${nf} selected` : "";
    if (nf && !els.more.open) els.more.open = true;
    els.skip.textContent = `Skip to the results (${plural(list.length, "operator")})`;
    els.fieldsSummary.textContent = F.fields.size ? `${F.fields.size} selected, ${F.mode === "all" ? "all required" : "any one suffices"}` : "none selected";
  }

  function resetFilters() {
    const sort = F.sort, dir = F.dir;
    F = readFilters(new URLSearchParams());
    F.sort = sort;
    F.dir = dir;
    refresh();
    els.q.focus();
  }

  // ---------------------------------------------------------------- table
  function buildTable() {
    const table = h("table", { class: "terms" });
    const cap = h("caption", { class: "sr-only" }, `Operators of the ${P.title.toLowerCase()}. Activate a row to open its detail.`);
    const thead = h("thead");
    const tr = h("tr");
    els.sortBtns = {};
    for (const c of COLUMNS) {
      const btn = h("button", { type: "button", class: "sort-btn" }, c.type ? h("span", { class: `col-mark t-${c.type}`, "aria-hidden": "true" }) : null, h("span", null, c.label), h("span", { class: "sort-ind", "aria-hidden": "true" }));
      btn.addEventListener("click", () => {
        if (F.sort === c.key) {
          if (F.dir === "asc") F.dir = "desc";
          else {
            F.sort = "";
            F.dir = "asc";
          }
        } else {
          F.sort = c.key;
          F.dir = "asc";
        }
        refresh();
      });
      const th = h("th", { scope: "col", class: c.cls }, btn);
      els.sortBtns[c.key] = th;
      tr.append(th);
    }
    thead.append(tr);
    els.tbody = h("tbody");
    table.append(cap, thead, els.tbody);
    els.empty = h("div", { class: "empty", hidden: true }, h("p", null, "No operator matches these filters."), h("button", { type: "button", class: "btn", onclick: () => resetFilters() }, "Reset filters"));
    els.tableWrap = h("div", { class: "table-wrap" }, table, els.empty);
    els.mainCol.append(els.tableWrap);
    rows = new Map();
    const labels = P.typeLabels;
    for (const t of P.terms) {
      const row = document.createElement("tr");
      row.tabIndex = 0;
      row.dataset.id = t.id;
      row.className = "term-row";
      const c = isObj(t.coefficients) ? t.coefficients : {};
      const cell = (k, label) => {
        const m = c[k]?.math;
        if (!hasMath(m)) return `<td class="col-coef col-${k} is-empty" data-label="${esc(label)}"><span class="none" aria-hidden="true">–</span><span class="sr-only">none</span></td>`;
        return `<td class="col-coef col-${k}" data-label="${esc(label)}"><div class="cell-math">${m.html}</div></td>`;
      };
      const pole = hasMath(c.pole?.math) ? "pole" : hasMath(c.evanescent?.math) ? "evanescent" : null;
      const poleCell = pole
        ? `<td class="col-coef col-pole${pole === "evanescent" ? " is-ev" : ""}" data-label="Singular"><div class="cell-math">${c[pole].math.html}</div>${pole === "evanescent" ? `<span class="ev-note">evanescent operator</span>` : ""}</td>`
        : `<td class="col-coef col-pole is-empty" data-label="Singular"><span class="none" aria-hidden="true">–</span><span class="sr-only">none</span></td>`;
      row.innerHTML =
        `<th scope="row" class="col-op"><a class="term-name" href="${esc(hrefTermStatic(t.id))}" tabindex="-1">${nameHTML(t)}</a>` +
        `<div class="cell-math op-math">${staticMath(t.operator)}</div></th>` +
        `<td class="col-fields" data-label="Fields"><span class="fchips">${arr(t.fields).map((f) => fieldChipHTML(P.fieldById.get(f), f)).join("")}</span></td>` +
        cell("tree", "Tree") + poleCell + cell("finite", "Finite") +
        `<td class="col-types" data-label="Type"><span class="tags">${typeTagsHTML(t._types, labels)}</span></td>`;
      rows.set(t.id, row);
    }
    loadViewer().then((m) => {
      if (!m?.tidyStatic) return;
      for (const row of rows.values()) for (const cell of row.querySelectorAll(".cell-math")) m.tidyStatic(cell);
    });
    els.tbody.addEventListener("click", (e) => {
      const row = e.target.closest("tr.term-row");
      if (!row) return;
      if (e.target.closest("a") && (e.metaKey || e.ctrlKey || e.shiftKey || e.button === 1)) return; // open in new tab
      e.preventDefault();
      navigate(hrefTerm(row.dataset.id));
    });
    els.tbody.addEventListener("keydown", (e) => {
      const row = e.target.closest("tr.term-row");
      if (!row) return;
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        navigate(hrefTerm(row.dataset.id));
      } else if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        e.preventDefault();
        const sib = e.key === "ArrowDown" ? row.nextElementSibling : row.previousElementSibling;
        sib?.focus();
      } else if (e.key === "Home" || e.key === "End") {
        e.preventDefault();
        (e.key === "Home" ? els.tbody.firstElementChild : els.tbody.lastElementChild)?.focus();
      }
    });
  }

  function hrefTermStatic(id) {
    return `#/phase/${enc(pid)}/term/${enc(id)}`;
  }

  // the first filter change after the view opened adds one history entry (Back returns to the unfiltered table);
  // later changes replace it
  let filterPushed = false;
  function refresh() {
    const target = detailId ? hrefTerm(detailId) : hrefTable();
    if (!filterPushed && !detailId && location.hash !== target) {
      filterPushed = true;
      history.pushState(null, "", target);
    } else navigate(target, { replace: true });
    filterAndSort();
  }

  function filterAndSort() {
    list = P.terms.filter((t) => matches(t, F));
    syncControls();
    if (F.sort && SORTERS[F.sort]) {
      const s = SORTERS[F.sort];
      const sign = F.dir === "desc" ? -1 : 1;
      // terms without a value stay last in both directions; the direction applies to the others
      const empty = EMPTY_LAST[F.sort];
      list = [...list].sort((a, b) => {
        if (empty) {
          const ea = empty(a), eb = empty(b);
          if (ea !== eb) return ea ? 1 : -1;
        }
        return sign * s(a, b) || a._order - b._order;
      });
    }
    els.tbody.replaceChildren(...list.map((t) => rows.get(t.id)));
    els.empty.hidden = list.length > 0;
    for (const c of COLUMNS) {
      const th = els.sortBtns[c.key];
      const on = F.sort === c.key;
      th.setAttribute("aria-sort", on ? (F.dir === "asc" ? "ascending" : "descending") : "none");
      th.querySelector(".sort-ind").textContent = on ? (F.dir === "asc" ? "▲" : "▼") : "";
    }
    paintResult();
    for (const [id, row] of rows) {
      const on = id === detailId;
      row.classList.toggle("is-active", on);
      if (on) row.setAttribute("aria-current", "true");
      else row.removeAttribute("aria-current");
    }
    if (detailId) paintDetailNav();
  }

  function paintResult() {
    const r = els.result;
    r.replaceChildren();
    const active = [];
    const label = (facet, id) => P.facets[facet]?.find((x) => x.id === id)?.label || id;
    const remove = (fn) => () => {
      fn();
      refresh();
    };
    if (F.q) active.push([`Text “${F.q}”`, remove(() => (F.q = ""))]);
    if (F.types.size) active.push([`Type: ${[...F.types].map((t) => P.typeLabels[t] || t).join(" or ")}`, remove(() => F.types.clear())]);
    if (F.sectors.size) active.push([`Sector: ${[...F.sectors].map((s) => label("sector", s)).join(" or ")}`, remove(() => F.sectors.clear())]);
    if (F.groups.size) active.push([`Group: ${[...F.groups].map((s) => label("group", s)).join(" or ")}`, remove(() => F.groups.clear())]);
    if (F.kinds.size) active.push([`Kind: ${[...F.kinds].map((s) => label("kind", s)).join(" or ")}`, remove(() => F.kinds.clear())]);
    if (F.legs) active.push([`${F.legs} legs`, remove(() => (F.legs = ""))]);
    if (F.ghost) active.push([`Ghost number ${F.ghost.replace("-", "−")}`, remove(() => (F.ghost = ""))]);
    if (F.fields.size) {
      const names = [...F.fields];
      active.push([
        h("span", null, F.fields.size > 1 ? (F.mode === "all" ? "Fields, all of: " : "Fields, any of: ") : "Field: ", h("span", { class: "fchips", html: names.map((f) => fieldChipHTML(P.fieldById.get(f), f)).join("") })),
        remove(() => F.fields.clear()),
      ]);
    }
    r.append(h("p", { class: "count" }, h("strong", null, list.length.toLocaleString("en-GB")), ` of ${plural(P.terms.length, "operator")}`, active.length ? " match" : ""));
    if (active.length) {
      r.append(
        h("ul", { class: "active-filters", "aria-label": "Active filters" },
          active.map(([text, fn]) => h("li", null, h("button", { type: "button", class: "afilter", onclick: fn, "aria-label": `Remove filter: ${typeof text === "string" ? text : text.textContent}` }, text, h("span", { class: "x", "aria-hidden": "true" }, "×"))))),
        h("button", { type: "button", class: "btn btn-quiet small", onclick: () => resetFilters() }, "Clear all"),
      );
    }
  }

  // ---------------------------------------------------------------- detail
  function paintDetailNav() {
    const nav = els.detail.querySelector(".detail-nav");
    if (!nav) return;
    const seq = list.some((t) => t.id === detailId) ? list : P.terms;
    const i = seq.findIndex((t) => t.id === detailId);
    const prev = i > 0 ? seq[i - 1] : null, next = i >= 0 && i < seq.length - 1 ? seq[i + 1] : null;
    const pa = nav.querySelector(".nav-prev"), na = nav.querySelector(".nav-next"), pos = nav.querySelector(".nav-pos");
    if (!pa || !na || !pos) return;
    const set = (a, t) => {
      if (t) {
        a.href = hrefTerm(t.id);
        a.removeAttribute("aria-disabled");
        a.title = t.name || t.id;
      } else {
        a.removeAttribute("href");
        a.setAttribute("aria-disabled", "true");
        a.title = "";
      }
    };
    set(pa, prev);
    set(na, next);
    pos.textContent = i >= 0 ? `${i + 1} of ${seq.length}${seq === list ? "" : " (outside the filter)"}` : "";
  }

  async function openDetail(id, { focus = true } = {}) {
    const t = P.termById.get(id);
    const D = els.detail;
    if (!wide.matches && !detailId) savedScroll = window.scrollY;
    rememberScroll();
    detailId = id;
    layout();
    viewerModule?.unmount?.(D);
    D.replaceChildren();
    if (!t) {
      D.append(
        h("div", { class: "detail-nav" }, h("a", { class: "btn", href: hrefTable() }, "Close")),
        h("div", { class: "error-box", role: "alert" }, h("p", { class: "error-title" }, `There is no term “${id}” in the ${P.title.toLowerCase()}.`), h("p", null, "The link may refer to an older build. The table lists every term of this build.")),
      );
      setTitle(P.title);
      filterAndSort();
      return;
    }
    setTitle(t.name, P.title);
    D.append(
      h("div", { class: "detail-nav" },
        h("a", { class: "btn nav-close", href: hrefTable(), "aria-label": "Close the detail (Escape)" }, wide.matches ? "× Close" : "‹ Table"),
        h("span", { class: "nav-spacer" }),
        h("a", { class: "btn nav-prev", "aria-label": "Previous term" }, "‹ Previous"),
        h("span", { class: "nav-pos muted small" }),
        h("a", { class: "btn nav-next", "aria-label": "Next term" }, "Next ›")),
    );
    const art = await renderTerm(P, t, { fieldHref: (f) => `#/phase/${enc(pid)}?f=${enc(f)}` });
    if (detailId !== id) return; // another term was opened meanwhile
    D.append(art);
    filterAndSort();
    // a term reopened through the history returns to where the reader left it
    const key = `${pid}:${id}`;
    const back = freshHash === location.hash ? 0 : DETAIL_SCROLL.get(key) || 0;
    freshHash = null;
    const restore = () => {
      if (wide.matches) D.scrollTop = back;
      else window.scrollTo(0, back);
    };
    restore();
    if (back) setTimeout(restore, 400);
    if (focus) D.querySelector("h2")?.focus({ preventScroll: true });
    const row = rows.get(id);
    if (row && wide.matches) {
      const rb = row.getBoundingClientRect();
      if (rb.top < 60 || rb.bottom > window.innerHeight) row.scrollIntoView({ block: "center" });
    }
  }

  function rememberScroll() {
    if (!detailId) return;
    DETAIL_SCROLL.set(`${pid}:${detailId}`, wide.matches ? els.detail.scrollTop : window.scrollY);
  }

  function closeDetailView() {
    const was = detailId;
    rememberScroll();
    detailId = null;
    viewerModule?.unmount?.(els.detail);
    els.detail.replaceChildren();
    layout();
    setTitle(P.title);
    filterAndSort();
    if (!wide.matches) window.scrollTo(0, savedScroll);
    const row = was && rows.get(was);
    if (row && row.isConnected) row.focus({ preventScroll: wide.matches });
  }

  function apply(r, { initial = false } = {}) {
    const NF = readFilters(r.query);
    const changed = writeFilters(NF) !== writeFilters(F);
    F = NF;
    const [, , sub, id] = r.parts;
    const want = sub === "term" && id ? id : null;
    if (want && want !== detailId) {
      filterAndSort();
      openDetail(want, { focus: !initial });
    } else if (!want && detailId) {
      closeDetailView();
    } else if (changed || initial) {
      filterAndSort();
    }
  }

  return {
    mount,
    update: (r) => apply(r),
    destroy() {
      rememberScroll();
      wide.removeEventListener("change", onWide);
      phone.removeEventListener("change", onWide);
      document.removeEventListener("keydown", onKey);
      clearTimeout(searchTimer);
    },
  };
}

// Full description of one term (used by the phase view)
async function renderTerm(P, t, { fieldHref }) {
  await viewerReady(400);
  const art = h("article", { class: "term-detail" });
  const c = isObj(t.coefficients) ? t.coefficients : {};
  const sectorLabel = P.facets.sector.find((s) => s.id === t.sector)?.label || t.sector;
  const kindLabel = t.kind ? P.facets.kind.find((s) => s.id === t.kind)?.label || t.kind : null;
  const facts = [
    sectorLabel ? `${sectorLabel} sector` : null,
    kindLabel,
    arr(t.groups).length ? arr(t.groups).join(", ") : null,
    Number.isFinite(t.legs) ? plural(t.legs, "leg") : null,
    Number.isFinite(t.ghost_number) ? `ghost number ${t.ghost_number > 0 ? "+" : ""}${t.ghost_number}`.replace("-", "−") : null,
    t.evanescent ? "evanescent operator" : null,
  ].filter(Boolean);
  art.append(
    h("header", { class: "td-head" },
      h("p", { class: "eyebrow" }, facts.join(" · ")),
      nameEl("h2", t, { tabindex: "-1" }),
      h("div", { class: "td-tags" }, h("span", { class: "tags" }, t._types.map((x) => typeTag(x, P.typeLabels))), h("code", { class: "term-id", title: "Term id" }, t.id)),
      h("div", { class: "td-fields" }, arr(t.fields).map((f) => fieldChip(P, f, { href: fieldHref(f) })))),
  );

  // operator and total
  const sOp = h("section", { class: "td-sec", "aria-labelledby": `h-op-${t.id}` }, h("h3", { id: `h-op-${t.id}` }, "Operator and full term"));
  sOp.append(h("div", { class: "labelled" }, h("span", { class: "mlabel" }, "Operator"), mathBox(t.operator, { label: "Operator" })));
  // the breaking insertion is not part of the action: its "total" is the insertion itself
  const onlyInsertion = t._types.length > 0 && t._types.every((x) => x === "breaking");
  const totalLabel = onlyInsertion ? "Term of the insertion Δ₁ (not part of the action)" : "Term in the renormalised Lagrangian";
  if (hasMath(t.total)) sOp.append(h("div", { class: "labelled" }, h("span", { class: "mlabel", html: onlyInsertion ? `Term of the insertion ${DELTA1} (not part of the action)` : esc(totalLabel) }), mathBox(t.total, { label: "Full term" })));
  art.append(sOp);

  // coefficients
  const sCo = h("section", { class: "td-sec", "aria-labelledby": `h-co-${t.id}` }, h("h3", { id: `h-co-${t.id}` }, "Coefficients"));
  const keys = TYPE_ORDER.filter((k) => isObj(c[k]));
  if (!keys.length) sCo.append(h("p", { class: "muted" }, "No coefficient is listed for this term."));
  for (const k of keys) {
    const e = c[k];
    const blk = h("div", { class: `coef-block t-${k}` }, h("div", { class: "coef-head" }, typeTag(k, P.typeLabels), h("span", { class: "coef-title" }, TYPE_LONG[k])));
    if (k === "finite" && hasMath(e.symbol)) {
      const sym = h("span", { class: "coef-symbol" }, "Report name ", h("span", { class: "inline-math", html: e.symbol.html }));
      blk.querySelector(".coef-head").append(sym);
    }
    blk.append(mathBox(e.math, { label: TYPE_LONG[k] }));
    if (hasMath(e.decomposition)) blk.append(h("div", { class: "labelled sub" }, h("span", { class: "mlabel" }, "Through the renormalisation constants"), mathBox(e.decomposition, { label: "Decomposition" })));
    if (e.note) blk.append(rich(e.note, { cls: "rich coef-note" }));
    sCo.append(blk);
  }
  art.append(sCo);

  // Feynman rules
  const sFR = h("section", { class: "td-sec", "aria-labelledby": `h-fr-${t.id}` }, h("h3", { id: `h-fr-${t.id}` }, "Feynman rules"));
  const rules = arr(t.rules).filter(isObj);
  if (!rules.length) sFR.append(h("p", { class: "muted" }, "No vertex function is listed for this term."));
  else {
    sFR.append(h("p", { class: "fr-note" }, "The Feynman rule is i times the vertex function V given here. All momenta are incoming unless a leg label says otherwise."));
    for (const r of rules) {
      const blk = h("div", { class: "rule" });
      blk.append(h("h4", null, r.title || "Vertex function"));
      const legs = arr(r.legs).filter(isObj);
      if (legs.length)
        blk.append(
          h("ol", { class: "legs", "aria-label": "Legs" },
            legs.map((l) => h("li", { title: (P.fieldById.get(l.field) || anyField(l.field))?.name || l.field },
              hasMath(l.label) ? h("span", { class: "leg-label", html: l.label.html }) : fieldChip(P, l.field)))),
        );
      const parts = isObj(r.parts) ? r.parts : {};
      for (const k of TYPE_ORDER) {
        if (!hasMath(parts[k])) continue;
        blk.append(h("div", { class: `part t-${k}` }, h("div", { class: "part-head" }, typeTag(k, P.typeLabels), h("span", { class: "muted small" }, `V, ${TYPE_LONG[k].toLowerCase()}`)), mathBox(parts[k], { label: `Vertex function, ${TYPE_LONG[k]}` })));
      }
      if (r.note) blk.append(rich(r.note, { cls: "rich rule-note" }));
      sFR.append(blk);
    }
  }
  art.append(sFR);

  // Exposé
  const sEx = h("section", { class: "td-sec expose", "aria-labelledby": `h-ex-${t.id}` }, h("h3", { id: `h-ex-${t.id}` }, "Exposé"));
  art.append(sEx);
  let T = null;
  try {
    T = await getTheory();
  } catch (e) {
    sEx.append(errorBox(e, "The theory card"));
  }
  const card = T && t.theory ? T.cards[t.theory] : null;
  if (card) {
    sEx.append(
      h("div", { class: "theory-card" },
        h("h4", null, card.title || t.theory),
        rich(card.body, { cls: "rich prose" }),
        h("div", { class: "card-foot" }, reportRefs(card.report, { cls: "refs inline" }), h("a", { class: "small", href: `#/theory/${enc(t.theory)}` }, "This card among all theory cards ›"))),
    );
  } else if (T && t.theory) sEx.append(h("p", { class: "muted" }, `The theory card “${t.theory}” is missing from the data.`));
  if (t.remarks) sEx.append(h("div", { class: "remarks" }, h("h4", null, "Remarks on this term"), rich(t.remarks, { cls: "rich prose" })));
  const refs = reportRefs(t.report);
  if (refs) sEx.append(h("div", { class: "td-refs" }, h("h4", null, "In the report"), refs));

  // graphs
  const sGr = h("div", { class: "td-graphs" }, h("h4", { class: "graphs-h" }, "Graphs behind this term"));
  const legendSlot = h("div", { class: "legend-slot" });
  sGr.append(legendSlot);
  sEx.append(sGr);
  let C = null;
  try {
    C = await getClasses();
  } catch (e) {
    sGr.append(errorBox(e, "The Green's functions"));
  }
  if (C) {
    const lg = graphLegend();
    if (lg) legendSlot.append(h("details", { class: "lb-legend" }, h("summary", null, "Line styles of the drawings"), lg));
    const own = isObj(t.classes) ? t.classes : {};
    const parents = [...new Set(arr(t.parents).filter((x) => typeof x === "string"))];
    const hasFinite = t._types.includes("finite");
    const brkTitle = hasFinite ? "Finite part, through the breaking" : "Classes that fix the finite coefficient (to zero)";
    const groupOf = (title, k, items, openN) => {
      const g = h("div", { class: `graph-group gg-${k}` }, h("h5", null, typeTag(k === "pole" ? "pole" : "breaking", P.typeLabels), " ", title));
      if (!items.length) g.append(h("p", { class: "muted small" }, k === "pole" ? "No class is listed for the singular part." : "No class with an insertion of the breaking is listed."));
      items.forEach(({ cls, via }, i) => g.append(classBlock(cls, { via, idPrefix: `td-${k}`, open: i < openN })));
      return g;
    };
    // the first class of each part is open (all of them when there are at most three)
    const openCount = (n, total) => (total <= 3 ? n : Math.min(1, n));
    const resolve = (ids) => arr(ids).map((cid) => C.classById.get(cid)).filter(Boolean);
    const insertionOnly = t._types.length > 0 && t._types.every((x) => x === "breaking");
    const insertionText = () => h("p", { class: "muted" }, rich(`${DELTA1} is a tree-level insertion and is not determined by loop graphs; its vertex enters every graph of the ghost-number-one classes of Green's functions (`, { tag: "span", cls: "" }),
      h("a", { href: "#/classes" }, "Green's functions"), ", ", h("a", { href: "#/report/92" }, "report Section 4.2"), ").");
    if (parents.length) {
      let UP = null;
      try {
        UP = await getPhase("unbroken");
      } catch {
        /* parents shown by id */
      }
      sGr.before(
        h("div", { class: "parents" },
          h("h4", null, "Descends from"),
          h("p", { class: "small muted" }, "Terms of the unbroken phase from which this term follows under the shift of the Higgs doublet. Their graphs are listed below; graphs of the broken phase are not drawn."),
          h("ul", { class: "parent-list" },
            parents.map((p) => {
              const pt = UP?.termById.get(p);
              return h("li", null, pt ? nameEl("a", pt, { href: `#/phase/unbroken/term/${enc(p)}`, title: pt.name }) : h("a", { href: `#/phase/unbroken/term/${enc(p)}` }, p), pt ? h("span", { class: "tags" }, pt._types.map((x) => typeTag(x, UP.typeLabels))) : null);
            }))),
      );
      const collect = (k) => {
        const map = new Map();
        for (const cls of resolve(own[k])) map.set(cls.id, { cls, via: [] });
        for (const p of parents) {
          const pt = UP?.termById.get(p);
          for (const cls of resolve(pt?.classes?.[k])) {
            if (!map.has(cls.id)) map.set(cls.id, { cls, via: [] });
            map.get(cls.id).via.push(pt ? nameEl("a", pt, { href: `#/phase/unbroken/term/${enc(p)}`, title: pt.name }) : p);
          }
        }
        return [...map.values()];
      };
      const pole = collect("pole"), brk = collect("breaking");
      if (!pole.length && !brk.length) sGr.append(insertionOnly ? insertionText() : h("p", { class: "muted" }, "No loop graph enters the coefficients of this term or of its parents: they have tree-level coefficients only."));
      else {
        const total = pole.length + brk.length;
        sGr.append(groupOf("Singular part", "pole", pole, openCount(pole.length, total)), groupOf(brkTitle, "breaking", brk, openCount(brk.length, total)));
      }
    } else {
      const pole = resolve(own.pole).map((cls) => ({ cls })), brk = resolve(own.breaking).map((cls) => ({ cls }));
      const total = pole.length + brk.length;
      if (!total) {
        const treeOnly = t._types.every((x) => x === "tree");
        sGr.append(insertionOnly ? insertionText() : h("p", { class: "muted" }, treeOnly ? "This operator has a tree-level coefficient only; no loop graph enters its determination." : "No class of Green's functions is attached to this term."));
      } else {
        if (pole.length || t._types.some((x) => x === "pole" || x === "evanescent")) sGr.append(groupOf("Singular part", "pole", pole, openCount(pole.length, total)));
        if (brk.length || hasFinite) sGr.append(groupOf(brkTitle, "breaking", brk, openCount(brk.length, total)));
      }
    }
  }
  return art;
}

// Table of contents: open beside the content on wide screens, a collapsed block above it on narrow ones.
function tocBox(label, count, kids) {
  const wideToc = window.matchMedia("(min-width: 1101px)").matches;
  const nav = h("nav", { class: "toc", "aria-label": label });
  const det = h("details", { class: "toc-det", open: wideToc }, h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), `Contents`, h("span", { class: "muted small" }, count)));
  const cols = h("div", { class: "toc-cols" });
  appendKids(cols, kids);
  det.append(cols);
  nav.append(det);
  return nav;
}

// Highlight the entry of a table of contents whose section is in view.
function scrollSpy(toc, sections, linkOf) {
  if (!("IntersectionObserver" in window) || !sections.length) return null;
  const visible = new Map();
  const io = new IntersectionObserver(
    (entries) => {
      for (const e of entries) visible.set(e.target, e.isIntersecting);
      const first = sections.find((s) => visible.get(s));
      if (!first) return;
      toc.querySelectorAll("a[aria-current]").forEach((a) => a.removeAttribute("aria-current"));
      const a = linkOf(first);
      if (a) {
        a.setAttribute("aria-current", "location");
        const tb = toc.getBoundingClientRect(), ab = a.getBoundingClientRect();
        if (toc.scrollHeight > toc.clientHeight && (ab.top < tb.top || ab.bottom > tb.bottom)) toc.scrollTop += ab.top - tb.top - tb.height / 3;
      }
    },
    { rootMargin: "-80px 0px -55% 0px" },
  );
  sections.forEach((s) => io.observe(s));
  return io;
}

// ------------------------------------------------------------------------------------------------ Green's functions

function ClassesView() {
  let C, root, content, spy;
  function focusClass(r, animate) {
    const cid = r.parts[1];
    if (!cid) return;
    const el = content.querySelector(`#gf-${CSS.escape(cid.replace(/[^A-Za-z0-9_-]/g, "_"))}`);
    if (!el) {
      content.prepend(h("div", { class: "error-box", role: "alert" }, h("p", { class: "error-title" }, `There is no class “${cid}”.`)));
      return;
    }
    if (!el.open) el.open = true;
    el.dispatchEvent(new Event("toggle"));
    requestAnimationFrame(() => el.scrollIntoView({ block: "start", behavior: animate ? smooth() : "auto" }));
    el.querySelector("summary")?.focus({ preventScroll: true });
    setTitle(classLabelText(C.classById.get(cid)) || cid, "Green's functions");
  }
  return {
    async mount(r0, r) {
      setTitle("Green's functions");
      const page = h("div", { class: "page gf" });
      root = h("div", { class: "with-toc" });
      r0.append(page);
      page.append(loadingBox());
      C = await getClasses();
      await getAllPhases();
      await viewerReady();
      page.replaceChildren();
      const nG = C.classes.reduce((s, c) => s + (Number.isFinite(c.n_graphs) ? c.n_graphs : 0), 0);
      const g0 = C.groups.filter((g) => g.ghost_number === 0), g1 = C.groups.filter((g) => g.ghost_number === 1), gx = C.groups.filter((g) => g.ghost_number !== 0 && g.ghost_number !== 1);
      const tocList = (title, gs) =>
        gs.length
          ? [h("p", { class: "toc-h" }, title), h("ol", null, gs.map((g) => h("li", null, h("a", { href: `#/classes`, dataset: { target: `grp-${g.key}` } }, g.title), h("span", { class: "toc-n" }, String(arr(g.classes).length)))))]
          : null;
      const toc = tocBox("Groups of Green's functions", plural(C.groups.length, "group"), [tocList("Ghost number 0", g0), tocList("Ghost number 1, with Δ₁", g1), tocList("Other", gx)]);
      toc.addEventListener("click", (e) => {
        const a = e.target.closest("a[data-target]");
        if (!a) return;
        e.preventDefault();
        document.getElementById(a.dataset.target)?.scrollIntoView({ block: "start", behavior: smooth() });
        document.getElementById(a.dataset.target)?.querySelector("h2")?.focus({ preventScroll: true });
      });
      content = h("div", { class: "toc-content" });
      page.append(
        h("header", { class: "page-head" },
          h("p", { class: "eyebrow" }, "One-loop 1PI functions"),
          h("h1", { tabindex: "-1" }, "Green's functions"),
          h("p", { class: "lede-small" }, `${plural(C.groups.length, "group")}, ${plural(C.classes.length, "class", "classes")} with definite external fields, ${plural(nG, "graph")}. The groups of ghost number 0 fix the singular counterterm; those of ghost number 1 carry one insertion of the breaking Δ₁ and fix the finite restoring counterterm. Open a class to see its result, the terms it determines and its graphs.`),
          graphLegend()),
      );
      const sec = (title, gs) => {
        if (!gs.length) return null;
        return h("section", { class: "gf-part" }, h("h2", { class: "part-title" }, title),
          gs.map((g) => {
            const cls = arr(g.classes).map((id) => C.classById.get(id)).filter(Boolean);
            const n = cls.reduce((s, c) => s + (Number.isFinite(c.n_graphs) ? c.n_graphs : 0), 0);
            return h("section", { class: "gf-group", id: `grp-${g.key}`, "aria-labelledby": `grp-h-${g.key}` },
              h("div", { class: "gf-head" },
                h("h2", { id: `grp-h-${g.key}`, tabindex: "-1" }, g.title),
                h("p", { class: "gf-facts muted small" }, [Number.isFinite(g.ghost_number) ? `ghost number ${g.ghost_number}` : null, plural(cls.length, "class", "classes"), plural(n, "graph")].filter(Boolean).join(" · ")),
                reportRefs(g.report, { cls: "refs inline" })),
              g.summary ? rich(g.summary, { cls: "rich prose gf-summary" }) : null,
              h("div", { class: "class-list" }, cls.map((c) => classBlock(c, { backlinks: true, idPrefix: "gf" }))));
          }));
      };
      appendKids(content, [sec("Ghost number 0: the singular counterterm", g0), sec("Ghost number 1: the breaking and the finite counterterm", g1), sec("Other classes", gx)]);
      root.append(toc, content);
      page.append(root);
      spy = scrollSpy(toc, [...content.querySelectorAll(".gf-group")], (s) => toc.querySelector(`a[data-target="${s.id}"]`));
      if (r.parts[1]) focusClass(r, false);
      else page.querySelector("h1")?.focus({ preventScroll: true });
    },
    update(r) {
      if (r.parts[1]) focusClass(r, true);
      else {
        window.scrollTo(0, 0);
        setTitle("Green's functions");
      }
    },
    destroy() {
      spy?.disconnect();
    },
  };
}

// ------------------------------------------------------------------------------------------------ theory

const CARD_ORDER = [
  "method-overview", "method-bmhv", "method-st", "method-beyond",
  "gauge-kinetic", "gauge-self", "gauge-fixing", "ghost", "fermion-kinetic", "fermion-gauge", "yukawa", "higgs-kinetic",
  "higgs-potential", "brst-sources", "breaking", "ct-fields", "ct-couplings", "ct-evanescent", "st-gauge", "st-fermion",
  "st-source", "st-scalar", "normalisation", "anomaly", "broken-shift", "broken-masses", "broken-mixing", "broken-tadpole",
  "broken-st", "broken-interactions",
];
const CARD_PARTS = [
  ["Method", ["method-overview", "method-bmhv", "method-st", "method-beyond"]],
  ["Classical action", ["gauge-kinetic", "gauge-self", "gauge-fixing", "ghost", "fermion-kinetic", "fermion-gauge", "yukawa", "higgs-kinetic", "higgs-potential", "brst-sources", "breaking"]],
  ["Singular counterterm", ["ct-fields", "ct-couplings", "ct-evanescent"]],
  ["Finite restoring counterterm", ["st-gauge", "st-fermion", "st-source", "st-scalar", "normalisation", "anomaly"]],
  ["Broken phase", ["broken-shift", "broken-masses", "broken-mixing", "broken-tadpole", "broken-st", "broken-interactions"]],
];

function TheoryView() {
  let content, spy;
  let anchor = null;
  const go = (r, animate) => {
    const key = r.parts[1];
    if (!key) return;
    const el = document.getElementById(`card-${key}`);
    if (el) {
      el.scrollIntoView({ block: "start", behavior: animate ? smooth() : "auto" });
      el.querySelector("h2")?.focus({ preventScroll: true });
      setTitle(el.querySelector("h2")?.textContent, "Theory");
      // the maths above the card may still change its layout (font, viewer): keep the card at the top for a while
      anchor?.disconnect();
      if ("ResizeObserver" in window && content) {
        const until = performance.now() + 1500;
        let last = null;
        anchor = new ResizeObserver(() => {
          if (performance.now() > until) return anchor?.disconnect();
          const top = el.getBoundingClientRect().top;
          if (last !== null && Math.abs(top - last) < 1) return;
          last = top;
          el.scrollIntoView({ block: "start", behavior: "auto" });
        });
        anchor.observe(content);
        setTimeout(() => anchor?.disconnect(), 1600);
        document.fonts?.ready.then(() => performance.now() < until && el.scrollIntoView({ block: "start", behavior: "auto" }));
      }
    }
  };
  return {
    async mount(r0, r) {
      setTitle("Theory");
      const page = h("div", { class: "page theory" });
      const root = h("div", { class: "with-toc" });
      r0.append(page);
      page.append(loadingBox());
      const T = await getTheory();
      const phases = await getAllPhases();
      await viewerReady();
      page.replaceChildren();
      const keys = Object.keys(T.cards);
      const known = new Set(CARD_ORDER);
      const parts = CARD_PARTS.map(([title, ks]) => [title, ks.filter((k) => T.cards[k])]).filter(([, ks]) => ks.length);
      const extra = keys.filter((k) => !known.has(k));
      if (extra.length) parts.push(["Other", extra]);
      const toc = tocBox("Theory cards", plural(keys.length, "card"),
        parts.map(([title, ks]) => [h("p", { class: "toc-h" }, title), h("ol", null, ks.map((k) => h("li", null, h("a", { href: `#/theory/${enc(k)}`, dataset: { key: k } }, T.cards[k].title || k))))]));
      content = h("div", { class: "toc-content" });
      page.append(h("header", { class: "page-head" }, h("p", { class: "eyebrow" }, "Background"), h("h1", { tabindex: "-1" }, "Theory"),
        h("p", { class: "lede-small" }, "Each operator of the tables refers to one of these cards. A card explains the part of the action it covers and how its coefficients are determined; the terms that refer to it are listed below the text.")));
      for (const [title, ks] of parts) {
        content.append(h("h2", { class: "part-title" }, title));
        for (const k of ks) {
          const c = T.cards[k];
          const refs = phases.map((P) => [P, P.terms.filter((t) => t.theory === k)]).filter(([, ts]) => ts.length);
          const art = h("article", { class: "card", id: `card-${k}`, "aria-labelledby": `card-h-${k}` },
            h("h2", { id: `card-h-${k}`, tabindex: "-1" }, c.title || k),
            rich(c.body, { cls: "rich prose" }),
            reportRefs(c.report, { cls: "refs inline" }));
          if (refs.length) {
            const box = h("div", { class: "card-terms" }, h("h3", null, "Terms referring to this card"));
            for (const [P, ts] of refs) {
              const ul = h("ul", { class: "term-links" }, ts.map((t) => h("li", null, nameEl("a", t, { href: `#/phase/${enc(P.id)}/term/${enc(t.id)}`, title: t.name }))));
              const det = h("details", { class: "card-phase", open: ts.length <= 12 }, h("summary", null, h("span", { class: "class-caret", "aria-hidden": "true" }), `${P.title}: ${plural(ts.length, "term")}`), ul);
              box.append(det);
            }
            art.append(box);
          } else art.append(h("p", { class: "muted small" }, "No operator refers to this card."));
          content.append(art);
        }
      }
      root.append(toc, content);
      page.append(root);
      spy = scrollSpy(toc, [...content.querySelectorAll(".card")], (s) => toc.querySelector(`a[data-key="${CSS.escape(s.id.replace(/^card-/, ""))}"]`));
      if (r.parts[1]) go(r, false);
      else page.querySelector("h1")?.focus({ preventScroll: true });
    },
    update(r) {
      if (r.parts[1]) go(r, true);
      else window.scrollTo(0, 0);
    },
    destroy() {
      spy?.disconnect();
      anchor?.disconnect();
    },
  };
}

// ------------------------------------------------------------------------------------------------ report

function ReportView() {
  let man, frameBox, input, pageInfo;
  const show = (r) => {
    const pages = Number(man.report_pages) || null;
    let page = Number.parseInt(r.parts[1], 10);
    if (!Number.isFinite(page) || page < 1) page = null;
    if (page && pages && page > pages) page = pages;
    const src = `${man.report_pdf || "report/one_loop_gamma5.pdf"}${page ? `#page=${page}` : ""}`;
    // a new frame: changing only the fragment of a loaded PDF does not move every viewer
    const frame = h("iframe", { src, title: `Report${page ? `, page ${page}` : ""}`, class: "pdf-frame" });
    frameBox.replaceChildren(frame);
    input.value = page || "";
    pageInfo.textContent = page ? `Opened at page ${page}${pages ? ` of ${pages}` : ""}.` : pages ? `${pages} pages.` : "";
    setTitle(page ? `Report, page ${page}` : "Report");
  };
  return {
    async mount(r0, r) {
      man = await getManifest();
      const pdf = man.report_pdf || "report/one_loop_gamma5.pdf";
      const root = h("div", { class: "page report" });
      r0.append(root);
      input = h("input", { type: "number", min: "1", max: man.report_pages || null, id: "report-page", inputmode: "numeric" });
      pageInfo = h("span", { class: "muted small", "aria-live": "polite" });
      const form = h("form", { class: "report-go", onsubmit: (e) => {
        e.preventDefault();
        const n = Number.parseInt(input.value, 10);
        if (Number.isFinite(n)) navigate(`#/report/${n}`);
      } }, h("label", { for: "report-page" }, "Page"), input, h("button", { type: "submit", class: "btn" }, "Open"));
      frameBox = h("div", { class: "pdf-box" });
      root.append(
        h("header", { class: "page-head compact" },
          h("div", null, h("p", { class: "eyebrow" }, "Report"), h("h1", { tabindex: "-1" }, man.title || "Report")),
          h("div", { class: "report-tools" }, form, pageInfo,
            h("a", { class: "btn", href: pdf, target: "_blank", rel: "noopener" }, "Open in a new tab"),
            h("a", { class: "btn btn-primary", href: pdf, download: "" }, "Download PDF"))),
        frameBox,
        h("p", { class: "muted small pdf-fallback" }, "If the document does not appear, your browser does not display PDF files inline; use the download link."),
      );
      show(r);
    },
    update(r) {
      show(r);
    },
  };
}

// ------------------------------------------------------------------------------------------------ resources

function ResourcesView() {
  return {
    async mount(r0) {
      setTitle("Resources");
      const man = await getManifest();
      const root = h("div", { class: "page resources" });
      r0.append(root);
      root.append(h("header", { class: "page-head" }, h("p", { class: "eyebrow" }, "Downloads"), h("h1", { tabindex: "-1" }, "Resources"),
        h("p", { class: "lede-small" }, "Inputs and results behind the page: the report and its Typst sources, the model specification and the outputs of the automation run, the final numbers of the report, its per-class results and checks, and the hand-written model and toolkit notes of the report. Third-party documents are not included.")));
      let data;
      try {
        data = await loadJSON(man.resources || "data/resources.json");
      } catch (e) {
        root.append(errorBox(e, "The list of resources"));
        return;
      }
      const list = (Array.isArray(data) ? data : arr(data?.resources)).filter((x) => isObj(x) && typeof x.path === "string");
      if (!list.length) {
        root.append(h("p", { class: "muted" }, "The list of resources is empty."));
        return;
      }
      const table = h("table", { class: "res-table" },
        h("thead", null, h("tr", null, h("th", { scope: "col" }, "Resource"), h("th", { scope: "col" }, "Description"), h("th", { scope: "col", class: "num" }, "Size"), h("th", { scope: "col" }, h("span", { class: "sr-only" }, "Download")))),
        h("tbody", null, list.flatMap((x, i) => [
          // optional "group" key: a heading row before the first entry of each group
          ...(typeof x.group === "string" && x.group && x.group !== list[i - 1]?.group
            ? [h("tr", { class: "res-group" }, h("th", { scope: "colgroup", colspan: "4" }, x.group))] : []),
          h("tr", null,
            h("th", { scope: "row", "data-label": "Resource" }, h("a", { href: x.path, target: "_blank", rel: "noopener" }, x.title || x.path.split("/").pop()), h("code", { class: "res-path" }, x.path)),
            h("td", { "data-label": "Description" }, x.description || ""),
            h("td", { class: "num", "data-label": "Size" }, Number.isFinite(x.bytes) ? fmtBytes(x.bytes) : ""),
            h("td", { class: "res-dl" }, h("a", { class: "btn small", href: x.path, download: "" }, "Download", h("span", { class: "sr-only" }, ` ${x.title || x.path}`))))])));
      root.append(h("div", { class: "res-wrap" }, table));
      root.querySelector("h1")?.focus({ preventScroll: true });
    },
    update() {},
  };
}

function MissingView() {
  return {
    async mount(r0) {
      setTitle("Not found");
      r0.append(h("div", { class: "page" }, h("h1", { tabindex: "-1" }, "No such view"), h("p", null, "This link does not correspond to a view of the site. ", h("a", { href: "#/" }, "Go to the overview"), ".")));
    },
    update() {},
  };
}

// ------------------------------------------------------------------------------------------------ boot

// inline formulae wider than their paragraph get their own scrolling box (phones, narrow columns)
function fitInlineMath(root) {
  if (!root) return;
  for (const m of root.querySelectorAll(".rich math:not([display='block'])")) {
    if (m.parentElement?.classList.contains("im-wide")) continue;
    const box = m.closest("p, li, dd, div");
    if (!box) continue;
    if (m.getBoundingClientRect().width > box.clientWidth + 1) {
      const w = document.createElement("span");
      w.className = "im-wide";
      m.replaceWith(w);
      w.append(m);
    }
  }
}
let fitTimer = 0;
function scheduleFit() {
  clearTimeout(fitTimer);
  fitTimer = setTimeout(() => {
    fitInlineMath(main());
    document.fonts?.ready.then(() => fitInlineMath(main()));
  }, 250);
}

async function boot() {
  theme.init();
  lightbox.init();
  // printing: paint the drawings of the open classes, which are otherwise painted when scrolled into view
  window.addEventListener("beforeprint", () => document.querySelectorAll("details[open] .graph-svg:not([data-painted])").forEach(paintGraph));
  // inline maths that overflow (theory cards, class summaries, remarks) are boxed after each change of the view
  if ("MutationObserver" in window) new MutationObserver(scheduleFit).observe(main(), { childList: true, subtree: true });
  window.addEventListener("resize", scheduleFit);
  // sticky offsets follow the real height of the top bar (it wraps on narrow screens)
  const bar = document.querySelector(".topbar");
  if (bar && "ResizeObserver" in window) {
    new ResizeObserver(() => document.documentElement.style.setProperty("--topbar-h", `${Math.round(bar.getBoundingClientRect().height)}px`)).observe(bar);
  }
  document.getElementById("skip")?.addEventListener("click", (e) => {
    e.preventDefault();
    const m = main();
    (m.querySelector("h1") || m).focus();
  });
  loadViewer(); // start early; failures leave the static MathML in place
  try {
    const man = await getManifest();
    const nav = document.getElementById("main-nav-phases");
    if (nav) {
      nav.replaceChildren(...phaseEntries(man).map((p) => h("li", null, h("a", { href: `#/phase/${enc(p.id)}`, dataset: { nav: `phase:${p.id}` } }, p.title || p.id))));
    }
    const foot = document.getElementById("foot-meta");
    if (foot) foot.textContent = [man.generated ? `Data generated ${man.generated}` : null, man.model ? `model ${man.model}` : null, man.loop_order ? `loop order ${man.loop_order}` : null].filter(Boolean).join(" · ");
  } catch (e) {
    main().replaceChildren(errorBox(e, "The site"));
    return;
  }
  window.addEventListener("hashchange", route);
  route();
}

boot();
