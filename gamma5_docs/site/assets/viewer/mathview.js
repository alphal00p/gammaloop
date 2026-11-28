// Expression viewer of the gamma5 documentation site (DATA_FORMAT.md, section 8).
//
//   renderMath(el, M, { interactive = false })  fills `el` with the math object M = {html, text, src, view}
//   ready()                                     Promise, resolved when the math font is loaded
//
// Static mode inserts the MathML printed by Symbolica (M.html). Interactive mode mounts Symbolica's own pager
// (symbolica_pager.js, exactly MathMLViewer.javascript()) on the pre-computed payload M.view: a static adapter
// answers the pager's requests from the table computed by gdocs.mathview (answer_key there builds the same keys
// as `answerKey` here). Controls whose request has no pre-computed answer are disabled, never sent anywhere.
// No network access beyond this directory: the math font is registered from ./STIXTwoMath-Regular.woff2 before
// the pager runs, so the pager's CDN fallback is never used.

const FONT_URL = new URL('./STIXTwoMath-Regular.woff2', import.meta.url).href;
const PAGER_URL = new URL('./symbolica_pager.js', import.meta.url).href;
const UNAVAILABLE = 'Not available in the static page';

let fontPromise = null;
let pagerPromise = null;

function installStyles(doc = document) {
  if (doc.getElementById('symbolica-notebook-math-font')) return;
  // The pager inserts a @font-face with a CDN URL only when this element does not exist yet.
  const font = doc.createElement('style');
  font.id = 'symbolica-notebook-math-font';
  font.textContent = `@font-face{font-family:"STIX Two Math";src:url("${FONT_URL}") format("woff2");font-display:swap}`;
  doc.head.append(font);
  const css = doc.createElement('style');
  css.id = 'gdocs-mathview-style';
  css.textContent = `
.mv { --foreground: currentColor; --background: var(--mv-bg, Canvas);
  --border: color-mix(in srgb, currentColor 22%, transparent);
  --muted: color-mix(in srgb, currentColor 9%, transparent);
  --muted-foreground: color-mix(in srgb, currentColor 68%, transparent);
  --primary: var(--mv-accent, #16846b); }
.mv math, math.mv-inline { font-family: "STIX Two Math", math; font-size: var(--mv-scale, 1.12em);
  font-synthesis: none; color: inherit; }
.mv-static { display: block; max-width: 100%; overflow-x: auto; overflow-y: hidden;
  padding: 0.35em 0.1em; scrollbar-width: thin; }
.mv-static > math { overflow: visible; }
.mv-flow { display: flex; flex-wrap: wrap; align-items: baseline; column-gap: 0; row-gap: 0.45em;
  max-width: 100%; }
.mv-flow > math { flex: 0 0 auto; }
.mv-static mtable[data-symbolica-product] > mtr > mtd,
.mv-static mtable[data-symbolica-sum] > mtr > mtd { text-align: left; }
.mv-static [data-symbolica-continuation-separator] { display: none; }
math.mv-inline { font-size: 1.08em; }
.mv-interactive { display: block; max-width: 100%; }
.mv-zoom { display: flex; gap: 0.35em; justify-content: flex-end; margin: 0.2em 0 -0.2em;
  font: 12px/1.4 system-ui, sans-serif; }
.mv-zoom button { font: inherit; color: inherit; background: transparent; cursor: pointer;
  border: 1px solid color-mix(in srgb, currentColor 25%, transparent); border-radius: 5px; padding: 1px 8px; }
.mv-zoom button:hover { background: color-mix(in srgb, currentColor 8%, transparent); }
`;
  doc.head.append(css);
}

/** Resolves when the math font is loaded (or failed to load; the page then uses the system math font). */
export function ready() {
  if (!fontPromise) {
    installStyles();
    fontPromise = document.fonts
      ? document.fonts.load('16px "STIX Two Math"').then(f => f.length > 0, () => false)
      : Promise.resolve(false);
  }
  return fontPromise;
}

function loadPager() {
  if (!pagerPromise) pagerPromise = import(PAGER_URL);
  return pagerPromise;
}

// ------------------------------------------------------------------------------------------------ static mode
function rootFlow(container) {
  // Long top-level sums: Symbolica marks the root row with flex wrapping. Chromium does not lay out MathML rows
  // as flex boxes, so (as the pager does) the terms become separate <math> roots in an HTML flex row.
  for (const math of container.querySelectorAll(':scope > math')) {
    const row = math.children.length === 1 ? math.firstElementChild : null;
    if (!row || row.localName !== 'mrow' || row.style.flexWrap !== 'wrap') continue;
    const flow = document.createElement('div');
    flow.className = 'mv-flow';
    for (const child of [...row.children]) {
      const term = document.createElementNS(math.namespaceURI, 'math');
      for (const { name, value } of math.attributes) {
        if (name !== 'style') term.setAttribute(name, value);
      }
      if (child.localName === 'mrow') {
        for (const { name, value } of child.attributes) if (name !== 'style') term.setAttribute(name, value);
        term.append(...child.childNodes);
      } else {
        term.append(child);
      }
      flow.append(term);
    }
    math.replaceWith(flow);
  }
}

function flattenProducts(container) {
  // Symbolica's static product rows are a fallback for a fixed width (the pager refits them by measurement).
  // A static page keeps one row per product and scrolls horizontally inside the element when it is too wide.
  for (const table of container.querySelectorAll('mtable[data-symbolica-product]')) {
    const rows = [...table.children].filter(r => r.localName === 'mtr');
    if (rows.length < 2) continue;
    const first = rows[0].querySelector(':scope > mtd > mrow');
    if (!first) continue;
    for (const row of rows.slice(1)) {
      const body = row.querySelector(':scope > mtd > mrow');
      if (body) first.append(...body.childNodes);
      row.remove();
    }
    for (const n of table.querySelectorAll('[data-symbolica-continuation]')) n.remove();
    for (const n of table.querySelectorAll('[data-symbolica-continuation-separator]')) {
      n.removeAttribute('data-symbolica-continuation-separator');
    }
  }
}

function renderStatic(el, M) {
  el.classList.add('mv', 'mv-static');
  el.innerHTML = M.html;
  flattenProducts(el);
  rootFlow(el);
  if (M.text && !el.hasAttribute('title')) el.setAttribute('aria-label', M.text);
}

/**
 * Tidy static MathML that was inserted without renderMath (for example table cells written as HTML strings):
 * products become one row and long top-level sums wrap, as in renderStatic.
 */
export function tidyStatic(el) {
  flattenProducts(el);
  rootFlow(el);
}

/** Close every interactive viewer below `root` (call before the elements are removed from the document). */
export function unmount(root) {
  if (!root) return;
  const list = root._mvHandle ? [root] : [];
  for (const el of root.querySelectorAll('.mv-interactive')) list.push(el);
  for (const el of list) {
    if (el._mvHandle) {
      try { el._mvHandle.close(); } catch (e) { /* already closed */ }
      el._mvHandle = null;
    }
  }
}

// ------------------------------------------------------------------------------------------------ interactive mode
/** The key of a request; identical to gdocs.mathview.answer_key. */
function answerKey(state, msg) {
  const a = msg.action;
  const path = msg.value === undefined ? null : msg.value;
  const lay = Number(msg.layout || 0);
  const ctx = msg.context || '';
  let k;
  if (a === 'preview') {
    const popup = Boolean(msg.popup);
    k = [state, a, path, lay, ctx, popup ? 1 : 0];
    if (popup) k.push(Number(msg.start || 0), Number(msg.size || 0));
  } else if (a === 'open') k = [state, a, path, lay, ctx, Number(msg.start || 0), Number(msg.size || 0)];
  else if (a === 'expand') k = [state, a, path, lay, ctx, Number(msg.size || 0)];
  else if (a === 'expand_all' || a === 'collapse') k = [state, a, 'value' in msg ? path : null, lay, ctx];
  else if (a === 'next') k = [state, a, lay];
  else if (a === 'previous' || a === 'back') k = [state, a];
  else return null;
  return JSON.stringify(k);
}

class StaticSource {
  constructor(view) {
    this.view = view;
    this.state = '0';
    this.page = view.pages[0];
    this.preview = null;      // {msg, value} of the last preview served
    this.onChange = null;
    this.closed = false;
  }
  has(msg, state = this.state) {
    const k = answerKey(state, msg);
    return Boolean(k && this.view.answers[k]);
  }
  send(msg) {
    if (msg.action === 'attach' || msg.action === 'dispose') return null;
    const k = answerKey(this.state, msg);
    const a = k ? this.view.answers[k] : null;
    let out;
    if (!a) {
      out = msg.action === 'preview'
        ? ['preview', { request: msg.request, error: UNAVAILABLE + '.' }]
        : ['page', this.page];
    } else if (a[0] === 'v') {
      const value = { ...this.view.previews[a[1]], context: a[2], request: msg.request };
      if (msg.popup) this.preview = { msg, value };
      out = ['preview', value];
    } else {
      this.state = a[2];
      this.page = this.view.pages[a[1]];
      this.preview = null;
      out = ['page', this.page];
    }
    if (this.onChange) setTimeout(this.onChange, 0);
    return out;
  }
}

function variant(page, layout) {
  if (!layout) return page;
  const alt = (page.layouts || []).find(l => (l.layout || 0) === layout);
  return alt ? { ...page, ...alt } : page;
}

function markerTarget(page, button) {
  const holder = button.closest('[data-symbolica-marker]');
  if (!holder) return null;
  const p = variant(page, Number(button.dataset.layout || 0));
  const i = (p.marker_ids || []).indexOf(Number(holder.dataset.symbolicaMarker));
  if (i < 0) return null;
  return { path: p.markers[i][0], start: (p.marker_starts || [])[i] || 0, layout: Number(button.dataset.layout || 0) };
}

function setAvailable(button, ok) {
  if (!button) return;
  if (ok) {
    if (button.dataset.mvOff) {
      delete button.dataset.mvOff;
      button.removeAttribute('title');
    }
    return;
  }
  if (!button.disabled) button.disabled = true;
  if (!button.dataset.mvOff) {
    button.dataset.mvOff = '1';
    button.title = UNAVAILABLE;
  }
}

function guardControls(root, src) {
  const page = src.page;
  const size = Number(page.page_size || 25);
  // omission markers of the main view
  for (const b of root.querySelectorAll('.math-content button.omission')) {
    const t = markerTarget(page, b);
    const ok = t && src.has({ action: 'preview', value: t.path, layout: t.layout, popup: true, start: t.start, size });
    setAvailable(b, Boolean(ok));
  }
  // toolbar: back, previous, next, collapse all
  const [back, previous, next] = root.querySelectorAll('.toolbar > button');
  if (back && !back.hidden) setAvailable(back, src.has({ action: 'back' }));
  if (previous && !previous.hidden && !previous.disabled) setAvailable(previous, src.has({ action: 'previous' }));
  if (next && !next.hidden && !next.disabled) setAvailable(next, src.has({ action: 'next', layout: page.layout || 0 }));
  const collapse = root.querySelector('.toolbar > .collapse-all');
  if (collapse && !collapse.hidden) setAvailable(collapse, src.has({ action: 'collapse', layout: page.layout || 0 }));
  // the pinned popup
  const pv = src.preview;
  const popup = root.querySelector('.preview');
  if (!pv || !popup || popup.hidden) return;
  const m = pv.msg, v = pv.value;
  const base = { value: m.value, layout: m.layout, context: m.context };
  const start = v.start || 0, psize = Number(m.size || size);
  const q = (extra) => src.has({ ...base, action: 'preview', popup: true, size: psize, ...extra });
  setAvailable(popup.querySelector('.preview-next'), v.next === null || q({ start: v.next }));
  if (start > 0) setAvailable(popup.querySelector('.preview-previous'), q({ start: Math.max(0, start - psize) }) || q({ start: 0 }));
  const select = popup.querySelector('.preview-size');
  if (select) for (const o of select.options) o.disabled = Number(o.value) !== psize && !q({ start: 0, size: Number(o.value) });
  setAvailable(popup.querySelector('.preview-inspect'), src.has({ ...base, action: 'open', start, size: psize }));
  setAvailable(popup.querySelector('.preview-expand'), src.has({ ...base, action: 'expand', size: psize }));
  setAvailable(popup.querySelector('.preview-show-all'), src.has({ ...base, action: 'expand_all' }));
  setAvailable(popup.querySelector('.preview-collapse'), src.has({ ...base, action: 'collapse' }));
  for (const b of popup.querySelectorAll('.preview-math button.omission')) {
    const t = markerTarget(v, b);
    if (!t) { setAvailable(b, false); continue; }
    const same = JSON.stringify(t.path) === JSON.stringify(m.value);
    const ok = same ? q({ start: t.start })
      : src.has({ action: 'preview', value: t.path, layout: t.layout, context: v.context, popup: true, start: t.start, size: psize });
    setAvailable(b, Boolean(ok));
  }
}

async function renderInteractive(el, M) {
  installStyles();
  const { mount } = await loadPager();
  el.classList.remove('mv-static');
  el.classList.add('mv', 'mv-interactive');
  el.replaceChildren();
  // match the math size of the surrounding text (the pager's default is 21px)
  const px = parseFloat(getComputedStyle(el).fontSize) || 16;
  const zoom = document.createElement('div');
  zoom.className = 'mv-zoom';
  zoom.innerHTML = '<button type="button" data-k="-" title="Zoom out (key −)" aria-label="Zoom out">−</button>' +
    '<button type="button" data-k="0" title="Reset zoom (key 0)" aria-label="Reset zoom, now 100%">100%</button>' +
    '<button type="button" data-k="+" title="Zoom in (key +)" aria-label="Zoom in">+</button>';
  const host = document.createElement('div');
  el.append(zoom, host);
  const src = new StaticSource(M.view);
  const mounted = mount({ el: host, page: src.page, send: msg => (src.closed ? null : src.send(msg)) });
  // a closed viewer answers no further request (its element may already have left the document)
  const handle = { ...mounted, close() { src.closed = true; src.onChange = null; return mounted.close(); } };
  const root = host.shadowRoot;
  const style = document.createElement('style');
  style.textContent = `.math-content, .preview-math { font-size: ${Math.round(px * 1.18)}px !important; }
    math { font-size: inherit !important; }
    .preview { color: var(--foreground); }
    .zoom-feedback { display: none !important; }`;
  root.prepend(style);
  // the pager's zoom message is hidden (it covered the formula); the reset button shows the current zoom
  const feedback = root.querySelector('.zoom-feedback');
  const reset = zoom.querySelector('[data-k="0"]');
  if (feedback && reset) {
    new MutationObserver(() => {
      const m = /(\d+)%/.exec(feedback.textContent || '');
      if (m) {
        reset.textContent = `${m[1]}%`;
        reset.setAttribute('aria-label', `Reset zoom, now ${m[1]}%`);
      }
    }).observe(feedback, { childList: true, characterData: true, subtree: true });
  }
  let frame = 0;
  const guard = () => {
    cancelAnimationFrame(frame);
    frame = requestAnimationFrame(() => guardControls(root, src));
  };
  src.onChange = guard;
  new MutationObserver(guard).observe(root, { subtree: true, childList: true, attributes: true,
    attributeFilter: ['disabled', 'hidden'] });
  guard();
  zoom.addEventListener('click', e => {
    const k = e.target.closest('button')?.dataset.k;
    if (!k) return;
    const region = root.querySelector('.math');
    region.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, composed: true }));
  });
  if (M.text) host.setAttribute('aria-label', M.text);
  el._mvHandle = handle;
  return handle;
}

/**
 * Fill `el` with the math object `M`. With `interactive: true` and `M.view` present, the Symbolica pager is
 * mounted on the pre-computed payload (returns a Promise of the mount handle); otherwise the static MathML is
 * inserted (returns null).
 */
export function renderMath(el, M, { interactive = false } = {}) {
  installStyles();
  if (el._mvHandle) {
    el._mvHandle.close();
    el._mvHandle = null;
  }
  if (interactive && M && M.view) {
    renderStatic(el, M);  // shown until the pager module is loaded
    return renderInteractive(el, M).catch(err => {
      console.warn('mathview: interactive viewer unavailable', err);
      el.classList.remove('mv-interactive');
      renderStatic(el, M);
      return null;
    });
  }
  renderStatic(el, M || { html: '' });
  return null;
}
