/* DOM measurement only. Native Rust owns layout, routing and SVG construction. */
export async function measure(request) {
  if (!request.labels.length) return { id: request.id, metrics: [] };
  // Isolate probes from notebook CSS, matching the exported SVG document.
  const frame = document.createElement('iframe');
  Object.assign(frame.style, {
    position: 'fixed', left: '-100000px', top: '0', visibility: 'hidden',
    pointerEvents: 'none', width: '1000px', height: '1000px', border: '0',
  });
  const loaded = new Promise(resolve => frame.addEventListener('load', resolve, { once: true }));
  frame.srcdoc = '<!doctype html><html><head><meta charset="utf-8"></head><body style="margin:0"></body></html>';
  document.body.append(frame);
  try {
    await loaded;
    return await measureDocument(request, frame.contentDocument);
  } finally { frame.remove(); }
}

async function measureDocument(request, doc) {
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
    const faces = await doc.fonts.load('9px "Linnet Offline Math"');
    await doc.fonts.ready;
    if (!faces.length || faces.some(face => face.status !== 'loaded')) {
      throw new Error('The bundled diagram math font could not be loaded');
    }
    const matrix = probes.getScreenCTM();
    if (!matrix || !Number.isFinite(matrix.a) || !Number.isFinite(matrix.d)
        || matrix.a <= 0 || matrix.d <= 0 || matrix.b !== 0 || matrix.c !== 0) {
      throw new Error('The diagram measurement SVG has no valid axis-aligned transform');
    }
    const metrics = rows.map(({ id, row, baseline }) => {
      const box = row.getBoundingClientRect();
      return {
        id, width: box.width / matrix.a, height: box.height / matrix.d,
        baseline: (baseline.getBoundingClientRect().top - box.top) / matrix.d,
      };
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
    let disposed = false;
    let frame;
    let resize;
    let shown;
    const status = document.createElement('div');
    status.setAttribute('role', 'status');
    status.textContent = 'Measuring diagram labels…';
    el.replaceChildren(status);

    const update = () => {
      const result = model.get('rendered');
      if (disposed || !result || result.id !== request.id) return;
      if (result.error) {
        status.textContent = result.error;
        status.setAttribute('role', 'alert');
        if (!frame) el.replaceChildren(status);
        return;
      }
      if (!result.html || result.html === shown) return;
      shown = result.html;
      resize?.disconnect();
      frame?.remove();
      frame = document.createElement('iframe');
      frame.title = 'Graph diagram';
      Object.assign(frame.style, { width: '100%', height: '120px', border: '0', display: 'block' });
      frame.addEventListener('load', async () => {
        if (disposed || !frame.contentDocument) return;
        const doc = frame.contentDocument;
        await doc.fonts.ready;
        if (disposed || frame.contentDocument !== doc) return;
        const fit = () => {
          if (frame?.isConnected) frame.style.height = `${Math.ceil(doc.body.getBoundingClientRect().height)}px`;
        };
        fit();
        resize = new ResizeObserver(fit);
        resize.observe(doc.body);
        const detail = { id: request.id, document: doc };
        doc.body.dispatchEvent(new frame.contentWindow.CustomEvent('linnet:rendered', { bubbles: true, detail }));
        el.dispatchEvent(new CustomEvent('linnet:rendered', { bubbles: true, detail }));
      });
      // A document mount executes the existing trusted SVG/domain controls.
      // The exact captured SVG is retained inside the authored HTML template.
      // Contain figure margins so the body's measured height includes them.
      frame.srcdoc = `<!doctype html><html><head><meta charset="utf-8"><style>html{overflow:hidden}body{margin:0;display:flow-root}</style></head><body data-linnet-frame-owner>${result.html}</body></html>`;
      el.replaceChildren(frame);
    };

    model.on('change:rendered', update);
    model.send({ action: 'attach', view });
    update();
    if (!shown) {
      try {
        const measurements = await measure(request);
        if (!disposed) model.send({ action: 'measure', measurements });
      } catch (error) {
        if (!disposed) {
          status.textContent = `Diagram label measurement failed: ${error.message || error}`;
          status.setAttribute('role', 'alert');
        }
      }
    }
    return () => {
      disposed = true;
      resize?.disconnect();
      model.off('change:rendered', update);
      model.send({ action: 'dispose', view });
      frame?.remove();
    };
  },
};
