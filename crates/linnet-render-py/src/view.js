/* DOM measurement only. Native Rust owns layout, routing and SVG construction. */
const documentOwners = new WeakMap();
const FONT_VERSIONS = 3;
const METRIC_COUNT = 2048;
const MATHML_BYTES = 2 * 1024 * 1024;
// Remote widget comms have 15 seconds to transfer the font before the separate
// five-second browser layout deadline starts. Static output survives either.
const FONT_TRANSFER_TIMEOUT_MS = 15000;
const TIMEOUT_MS = 5000;
const PROBE_VERSION = 1;
const EMPTY_FONT_KEY = 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855';

function ownerFor(doc) {
  let owner = documentOwners.get(doc);
  if (!owner) {
    owner = new DocumentMeasurements(doc);
    documentOwners.set(doc, owner);
  }
  return owner;
}

function waitFor(promise, win, signal, message, timeout = TIMEOUT_MS) {
  return new Promise((resolve, reject) => {
    let settled = false;
    const finish = (complete, value) => {
      if (settled) return;
      settled = true;
      win.clearTimeout(timer);
      signal.removeEventListener('abort', cancel);
      complete(value);
    };
    const cancel = () => finish(reject, signal.reason || new Error('Diagram view disposed'));
    const timer = win.setTimeout(() => {
      const error = new Error(message);
      error.name = 'LinnetMeasurementTimeout';
      finish(reject, error);
    }, timeout);
    signal.addEventListener('abort', cancel, { once: true });
    if (signal.aborted) cancel();
    Promise.resolve(promise).then(value => finish(resolve, value), error => finish(reject, error));
  });
}

async function fontCSS(request, model, win, signal) {
  if (typeof request.font_css === 'string') return request.font_css;
  if (!model) throw new Error('The diagram measurement request has no bundled font CSS');
  let changed;
  try {
    const available = new Promise((resolve, reject) => {
      let received = false;
      changed = () => {
        const font = model.get('font');
        if (!font || font.key !== request.font_key) return;
        received = true;
        if (typeof font.css !== 'string') reject(new Error('The bundled diagram font CSS is missing'));
        else resolve(font.css);
      };
      // A backend may answer synchronously; subscribe before reading or sending.
      model.on('change:font', changed);
      changed();
      if (!received && !signal.aborted) model.send({ action: 'font', key: request.font_key });
    });
    return await waitFor(available, win, signal, 'Timed out waiting for the bundled diagram font',
      FONT_TRANSFER_TIMEOUT_MS);
  } finally {
    if (changed) model.off('change:font', changed);
  }
}

async function fontReady(doc, owner, signal = owner.abort.signal) {
  await waitFor((async () => {
    const faces = await doc.fonts.load('9px "Linnet Offline Math"');
    await doc.fonts.ready;
    if (!faces.length || faces.some(face => face.status !== 'loaded')) {
      throw new Error('The bundled diagram math font could not be loaded');
    }
  })(), owner.win, signal, 'Timed out loading the bundled diagram font');
}

class DocumentMeasurements {
  constructor(doc) {
    this.doc = doc;
    this.win = doc.defaultView;
    this.abort = new this.win.AbortController();
    this.fonts = new Map();
    this.metrics = new Map();
    this.inflight = new Map();
    this.mathmlBytes = 0;
    const pagehide = event => {
      // A BFCache document and its Blob URLs remain live on restoration.
      if (event.persisted) return;
      this.abort.abort(new Error('The diagram document was unloaded'));
      for (const font of this.fonts.values()) {
        font.abort.abort(this.abort.signal.reason);
        this.revokeFont(font);
      }
      this.fonts.clear();
      this.metrics.clear();
      this.mathmlBytes = 0;
      this.win.removeEventListener('pagehide', pagehide);
    };
    this.win.addEventListener('pagehide', pagehide);
  }

  retainFont(request, model) {
    if (!request.font_key) throw new Error('The diagram measurement request has no font key');
    if (this.abort.signal.aborted) throw this.abort.signal.reason;
    let font = this.fonts.get(request.font_key);
    if (!font) {
      font = { refs: 0, urls: [], css: null, abort: new this.win.AbortController() };
      this.fonts.set(request.font_key, font);
      font.ready = (async () => {
        const css = await fontCSS(request, model, this.win, font.abort.signal);
        if (font.abort.signal.aborted) throw font.abort.signal.reason;
        font.css = css.replace(/data:font\/woff2;base64,([A-Za-z0-9+/=]+)/g, (_, encoded) => {
          const bytes = Uint8Array.from(this.win.atob(encoded), char => char.charCodeAt(0));
          const url = this.win.URL.createObjectURL(new this.win.Blob([bytes], { type: 'font/woff2' }));
          font.urls.push(url);
          return url;
        });
        if (!font.urls.length) throw new Error('The bundled diagram font data is missing');
        return font;
      })().catch(error => {
        if (this.fonts.get(request.font_key) === font) this.fonts.delete(request.font_key);
        this.revokeFont(font);
        throw error;
      });
    }
    this.fonts.delete(request.font_key);
    this.fonts.set(request.font_key, font);
    font.refs += 1;
    this.trimFonts();
    let released = false;
    return {
      ready: font.ready,
      release: () => {
        if (released) return;
        released = true;
        font.refs -= 1;
        if (!font.refs && font.css === null) {
          font.abort.abort(new Error('Diagram font request disposed'));
          if (this.fonts.get(request.font_key) === font) this.fonts.delete(request.font_key);
        }
        this.trimFonts();
      },
    };
  }

  revokeFont(font) {
    for (const url of font.urls) this.win.URL.revokeObjectURL(url);
    font.urls.length = 0;
  }

  trimFonts() {
    // Keep idle fonts warm across detach/remount, but never evict a live lease.
    for (const [key, font] of this.fonts) {
      if (this.fonts.size <= FONT_VERSIONS) break;
      if (font.refs) continue;
      this.fonts.delete(key);
      this.revokeFont(font);
    }
  }

  async measure(request, font) {
    for (const label of request.labels) {
      if (typeof label.mathml !== 'string') throw new Error(`Missing MathML for diagram label ${label.id}`);
    }
    const context = [PROBE_VERSION, request.font_key, request.coordinate_system,
      Object.keys(request.font_style).sort().map(key => [key, request.font_style[key]])];
    const missing = [];
    const values = request.labels.map(label => {
      const key = JSON.stringify([...context, label.mathml]);
      const cached = this.metrics.get(key);
      if (cached) {
        this.metrics.delete(key);
        this.metrics.set(key, cached);
        return cached.metric;
      }
      let pending = this.inflight.get(key);
      if (!pending) {
        let resolve, reject;
        pending = new Promise((yes, no) => { resolve = yes; reject = no; });
        this.inflight.set(key, pending);
        missing.push({ key, label, resolve, reject });
      }
      return pending;
    });
    if (missing.length) void this.measureMissing(request, font, missing);
    const metrics = await Promise.all(values);
    // Receipts belong to the current request, never to a cached source label.
    return { id: request.id, metrics: metrics.map((metric, index) => ({ id: request.labels[index].id, ...metric })) };
  }

  async measureMissing(request, font, missing) {
    font.refs += 1;
    const frame = this.doc.createElement('iframe');
    frame.setAttribute('sandbox', 'allow-same-origin');
    Object.assign(frame.style, {
      position: 'fixed', left: '-100000px', top: '0', visibility: 'hidden',
      pointerEvents: 'none', width: '1000px', height: '1000px', border: '0',
    });
    try {
      // Each miss batch owns its probe DOM; concurrent batches cannot mutate it.
      const loaded = new Promise(resolve => frame.addEventListener('load', resolve, { once: true }));
      frame.srcdoc = '<!doctype html><html><head><meta charset="utf-8"></head><body style="margin:0"></body></html>';
      this.doc.body.append(frame);
      await waitFor(loaded, this.win, this.abort.signal, 'Timed out opening the diagram measurement frame');
      if (!frame.contentDocument) throw new Error('The diagram measurement frame has no document');
      const receipt = await measureDocument({ ...request, labels: missing.map(entry => entry.label), font_css: font.css },
        frame.contentDocument, this);
      const encoder = new TextEncoder();
      missing.forEach((entry, index) => {
        const { width, height, baseline } = receipt.metrics[index];
        const metric = { width, height, baseline };
        const bytes = encoder.encode(entry.label.mathml).byteLength;
        if (bytes <= MATHML_BYTES) {
          this.metrics.set(entry.key, { metric, bytes });
          this.mathmlBytes += bytes;
          while (this.metrics.size > METRIC_COUNT || this.mathmlBytes > MATHML_BYTES) {
            const [key, oldest] = this.metrics.entries().next().value;
            this.mathmlBytes -= oldest.bytes;
            this.metrics.delete(key);
          }
        }
        entry.resolve(metric);
      });
    } catch (error) {
      for (const entry of missing) entry.reject(error);
    } finally {
      for (const entry of missing) this.inflight.delete(entry.key);
      frame.remove();
      font.refs -= 1;
      this.trimFonts();
    }
  }
}

export async function measure(request) {
  if (!request.labels.length) return { id: request.id, metrics: [] };
  const owner = ownerFor(document);
  const font = owner.retainFont(request);
  try {
    return await owner.measure(request, await font.ready);
  } finally { font.release(); }
}

async function measureDocument(request, doc, owner) {
  const style = doc.createElement('style');
  style.textContent = request.font_css;
  const svgNamespace = 'http://www.w3.org/2000/svg';
  const xhtmlNamespace = 'http://www.w3.org/1999/xhtml';
  const probes = doc.createElementNS(svgNamespace, 'svg');
  probes.setAttribute('width', '1000pt');
  probes.setAttribute('height', '1000pt');
  probes.setAttribute('viewBox', '0 0 1000 1000');
  Object.assign(probes.style, {
    position: 'fixed', left: '-100000px', top: '0', visibility: 'hidden',
    pointerEvents: 'none', margin: '0', padding: '0',
  });
  const rows = request.labels.map(label => {
    // MathML script sizing depends on its local CSS pixels. Typeset in the
    // same 9px foreignObject context as the saved SVG before its outer scale.
    const object = doc.createElementNS(svgNamespace, 'foreignObject');
    object.setAttribute('width', '1000000');
    object.setAttribute('height', '1000000');
    const row = doc.createElementNS(xhtmlNamespace, 'div');
    row.className = 'linnet-math-label';
    Object.assign(row.style, request.font_style, {
      display: 'block', width: 'max-content', margin: '0', padding: '0',
      border: '0', textIndent: '0', textTransform: 'none',
    });
    const math = doc.createElementNS('http://www.w3.org/1998/Math/MathML', 'math');
    math.innerHTML = label.mathml;
    const baseline = doc.createElementNS(xhtmlNamespace, 'span');
    Object.assign(baseline.style, {
      display: 'inline-block', verticalAlign: 'baseline', width: '0', height: '0',
      margin: '0', padding: '0', border: '0',
    });
    row.append(math, baseline);
    object.append(row);
    probes.append(object);
    return { id: label.id, row, baseline };
  });
  doc.head.append(style);
  doc.body.append(probes);
  try {
    await fontReady(doc, owner);
    const matrix = probes.getScreenCTM();
    if (!matrix || !Number.isFinite(matrix.a) || !Number.isFinite(matrix.d)
        || matrix.a <= 0 || matrix.d <= 0 || matrix.b !== 0 || matrix.c !== 0) {
      throw new Error('The diagram measurement SVG has no valid axis-aligned transform');
    }
    const metrics = rows.map(({ id, row, baseline }) => {
      const box = row.getBoundingClientRect();
      const metric = {
        width: box.width / matrix.a, height: box.height / matrix.d,
        baseline: (baseline.getBoundingClientRect().top - box.top) / matrix.d,
      };
      if (!Object.values(metric).every(Number.isFinite) || metric.width < 0 || metric.height < 0) {
        throw new Error(`The diagram label ${id} has invalid measured dimensions`);
      }
      return { id, ...metric };
    });
    return { id: request.id, metrics };
  } finally {
    probes.remove();
    style.remove();
  }
}

export default {
  async render({ model, el }) {
    const view = globalThis.crypto?.randomUUID?.() || `${Date.now()}-${Math.random()}`;
    const request = model.get('request');
    const owner = ownerFor(el.ownerDocument);
    const needsFont = request.labels.length || (request.font_key && request.font_key !== EMPTY_FONT_KEY);
    let disposed = false;
    let frame;
    let mountAbort;
    let resize;
    let controls;
    let shown;
    let font;
    let revision = 0;
    const resource = () => {
      if (!font) font = owner.retainFont(request, model);
      return font.ready;
    };
    const status = owner.doc.createElement('div');
    status.setAttribute('role', 'status');
    status.textContent = 'Measuring diagram labels…';
    el.replaceChildren(status);

    const fail = error => {
      if (disposed) return;
      status.textContent = `Diagram label measurement failed: ${error.message || error}`;
      status.setAttribute('role', 'alert');
      if (!frame) {
        font?.release();
        font = null;
        el.replaceChildren(status);
      }
    };
    const update = () => {
      const result = model.get('rendered');
      if (disposed || !result || result.id !== request.id) return;
      if (result.error) {
        fail(new Error(result.error));
        return;
      }
      if (!result.html || result.html === shown) return;
      shown = result.html;
      const current = ++revision;
      try {
        resize?.disconnect();
        controls?.disconnect();
        mountAbort?.abort();
        frame?.remove();
        const mounted = owner.doc.createElement('iframe');
        mounted.setAttribute('sandbox', 'allow-same-origin');
        mountAbort = new owner.win.AbortController();
        const signal = mountAbort.signal;
        frame = mounted;
        mounted.title = 'Graph diagram';
        Object.assign(mounted.style, {
          width: '100%', height: '120px', border: '0', display: 'block', visibility: 'hidden',
        });
        mounted.addEventListener('load', () => {
          try {
            if (disposed || current !== revision || frame !== mounted || !mounted.contentDocument) return;
            const doc = mounted.contentDocument;
            // The sandbox blocks all document scripts, including imported markup.
            // Controls come from Rust templates before label substitution. Run
            // only exact trusted sources in the parent, also for lazy children.
            const trusted = new Set(model.get('scripts') || []);
            const executed = new WeakSet();
            const activate = () => {
              for (const script of doc.querySelectorAll('script')) {
                if (executed.has(script) || !trusted.has(script.textContent)) continue;
                executed.add(script);
                new Function('document', 'window', script.textContent)(doc, mounted.contentWindow);
              }
            };
            controls = new owner.win.MutationObserver(activate);
            controls.observe(doc.body, { childList: true, subtree: true });
            activate();
            const fit = () => {
              if (mounted.isConnected) mounted.style.height = `${Math.ceil(doc.body.getBoundingClientRect().height)}px`;
            };
            fit();
            resize = new owner.win.ResizeObserver(fit);
            resize.observe(doc.body);
            // Snapshot controls and fallback-font layout are ready now. Exact
            // measurement still waits for STIX separately; displaying a graph
            // must never wait for an unanswered font comm or fonts.ready.
            if (needsFont) {
              void Promise.resolve().then(() => {
                if (disposed || signal.aborted || current !== revision || frame !== mounted) return;
                return resource();
              }).then(resolved => {
                if (!resolved || disposed || signal.aborted || current !== revision
                    || frame !== mounted || mounted.contentDocument !== doc) return;
                const style = doc.createElement('style');
                style.setAttribute('data-linnet-font', '');
                style.textContent = resolved.css;
                doc.head.append(style);
                // Saved exports use block while exact measurement loads STIX.
                // A live snapshot must keep fallback labels painted meanwhile.
                for (const face of doc.fonts) {
                  if (face.family.replace(/["']/g, '') === 'Linnet Offline Math') face.display = 'swap';
                }
                // The existing ResizeObserver fits any changed font metrics.
                // Keep this document and its selection/camera state intact.
              }).catch(() => { /* Keep the visible fallback-font snapshot. */ });
            }
            mounted.style.removeProperty('visibility');
            const detail = { id: request.id, document: doc };
            doc.body.dispatchEvent(new mounted.contentWindow.CustomEvent('linnet:rendered', { bubbles: true, detail }));
            el.dispatchEvent(new owner.win.CustomEvent('linnet:rendered', { bubbles: true, detail }));
          } catch (error) {
            if (frame === mounted) {
              resize?.disconnect();
              controls?.disconnect();
              mountAbort?.abort();
              mounted.remove();
              frame = null;
              fail(error);
            }
          }
        }, { once: true, signal });
        mounted.srcdoc = `<!doctype html><html><head><meta charset="utf-8"><style>html{overflow:hidden}body{margin:0;display:flow-root}</style></head><body data-linnet-frame-owner>${result.html}</body></html>`;
        el.replaceChildren(mounted);
      } catch (error) { fail(error); }
    };

    model.on('change:rendered', update);
    model.send({ action: 'attach', view });
    update();
    if (request.pending !== false && (request.labels.length || !shown)) {
      try {
        const measurements = request.labels.length
          // Resolve the font before starting the layout timer: a remote comm
          // transfer must retain its independent 15-second allowance.
          ? await waitFor(owner.measure(request, await resource()), owner.win,
              owner.abort.signal, 'Diagram measurement exceeded 5 seconds')
          : { id: request.id, metrics: [] };
        if (!disposed) model.send({ action: 'measure', measurements });
      } catch (error) { fail(error); }
    }
    return () => {
      disposed = true;
      resize?.disconnect();
      controls?.disconnect();
      mountAbort?.abort();
      model.off('change:rendered', update);
      model.send({ action: 'dispose', view });
      frame?.remove();
      frame = null;
      font?.release();
    };
  },
};
