import {frameSource, initializeNotebookSources, initializePopups} from './code-blocks.js';
import {load as loadEditor, typstLanguage} from './documentation-editor.js';
// Shared presentation and direct documentation cells. Scientific kernels are
// created only by Run; editor intelligence never starts a Python interpreter.
const assets = new URL('./', import.meta.url);
let configuration;
let configLocation = new URL('documentation-cells-config.json', assets);
const views = new Map();
const workers = new Map();
const outputFrames = new Map();
let serial = 0;
const outputTheme = () => document.documentElement.dataset.theme || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
const syncOutputThemes = () => {
  document.querySelectorAll('.documentation-type-info').forEach(popup => { popup.dataset.codeTheme = outputTheme(); });
  for (const [id, frame] of outputFrames) frame.contentWindow?.postMessage({type: 'documentation-output-theme', id, theme: outputTheme()}, '*');
};
new MutationObserver(syncOutputThemes).observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', syncOutputThemes);

addEventListener('message', event => {
  if (event.data?.type !== 'documentation-output-size') return;
  const frame = outputFrames.get(event.data.id);
  if (!frame || event.source !== frame.contentWindow || !Number.isFinite(event.data.height)) return;
  frame.style.height = `${Math.max(32, Math.min(1200, Math.ceil(event.data.height)))}px`;
});
const config = () => configuration ||= fetch(configLocation).then(async response => {
  if (!response.ok) throw Error(`Documentation cell configuration: HTTP ${response.status}`);
  const value = await response.json();
  value.compilerRoot = new URL(value.compilerRoot, configLocation).href;
  value.runtime = new URL(value.runtime, configLocation).href;
  value.cells = await fetch(new URL('cells.json', value.compilerRoot)).then(r => r.json());
  return value;
});
async function api(runtime) {
  return loadEditor(runtime || (await config()).runtime);
}
async function worker(kind) {
  if (!workers.has(kind)) workers.set(kind, (async () => {
    const value = await config();
    const instance = new Worker(new URL(`documentation-${kind}-worker.js`, value.compilerRoot), {type: 'module'});
    const pending = new Map();
    const initialized = new Promise((resolve, reject) => {
      instance.onerror = event => {
        const error = Error(event.message || `${kind} worker failed`);
        reject(error);
        for (const request of pending.values()) request.reject(error);
        pending.clear();
      };
      instance.onmessage = ({data}) => {
        if (data.type === 'measure') {
          import(new URL('graph-measure.js', value.compilerRoot).href)
            .then(module => module.measure(data.request))
            .then(receipt => instance.postMessage({type: 'measurement', id: data.id, receipt}),
              error => instance.postMessage({type: 'measurement', id: data.id, error: String(error)}));
          return;
        }
        if (data.type === 'ready') { resolve(); return; }
        if (data.type === 'error' && !data.id) { reject(Error(data.message)); return; }
        const request = pending.get(data.id);
        if (!request) return;
        pending.delete(data.id);
        data.type === 'error' ? request.reject(Error(data.message)) : request.resolve(data);
      };
    });
    instance.postMessage({type: 'init', config: value});
    await initialized;
    return data => new Promise((resolve, reject) => {
      const id = ++serial;
      pending.set(id, {resolve, reject});
      instance.postMessage({...data, id});
    });
  })().catch(error => { workers.delete(kind); throw error; }));
  return workers.get(kind);
}
function sourceFor(cell, cells, seen = new Set()) {
  if (seen.has(cell.id)) return [];
  seen.add(cell.id);
  if (cell.reason) throw Error(cell.reason);
  const dependencies = (cell.depends || []).flatMap(id => {
    const prerequisite = cells.find(c => c.id === id);
    if (!prerequisite) throw Error(`Missing prerequisite ${id}`);
    return sourceFor(prerequisite, cells, seen);
  });
  return [...dependencies, ...(cell.setup ? [[cell.id + ':setup', cell.setup]] : []),
    [cell.id, views.get(cell.id)?.state.doc.toString() ?? cell.source]];
}
function notebookContext(editor, cm) {
  const notebook = editor.closest('.live-notebook');
  if (!notebook) return '';
  const sources = [];
  for (const source of notebook.querySelectorAll('marimo-island > marimo-ui-element > marimo-code-editor')) {
    // Marimo dependencies can point forwards in authored source order.
    if (source === editor) continue;
    const native = source.shadowRoot?.querySelector('.cm-editor');
    const view = native && cm.EditorView.findFromDOM(native);
    sources.push(view?.state.doc.toString() ?? JSON.parse(source.dataset.initialValue || '""'));
  }
  return sources.join('\n\n');
}
async function intelligence(view, context, status, nativeAPI) {
  initializePopups(view.dom.getRootNode().host ? view.dom.getRootNode() : document.body);
  const cm = nativeAPI || await api();
  if (view.state.field(cm.intelligenceField, false)) return;
  let generation = 0;
  let timer;
  let signature;
  const request = async (type, position) => {
    const prefix = context();
    const document = view.state.doc;
    const line = document.lineAt(position);
    const start = prefix ? prefix + '\n\n' : '';
    const send = await worker('ty');
    const response = await send({type, source: start + document.toString(),
      line: start.split('\n').length - 1 + line.number - 1, column: position - line.from});
    status('');
    return response.result;
  };
  const fail = error => status(`Type hints unavailable: ${error.message}. Run remains available.`);
  const text = (value, markdown = true) => {
    const dom = document.createElement('div');
    dom.className = 'documentation-type-info';
    dom.dataset.codeTheme = document.documentElement.dataset.theme || (view.state.facet(cm.EditorView.darkTheme) ? 'dark' : 'light');
    const source = typeof value === 'string' ? value : value?.value || JSON.stringify(value);
    if (markdown) dom.innerHTML = cm.renderMarkdown(source);
    else { const code = document.createElement('code'); code.textContent = source; dom.append(code); }
    return dom;
  };
  const complete = async context => {
    const word = context.matchBefore(/[\w]*$/);
    if (!word || (!context.explicit && word.from === word.to)) return null;
    try {
      const results = await request('complete', context.pos);
      if (context.aborted || !view.state.field(cm.intelligenceField, false)) return null;
      return {from: word.from, options: (results || []).map(item => ({
        label: item.label, apply: item.insertText || item.label,
        detail: item.detail, info: item.documentation ? () => text(item.documentation) : undefined,
      }))};
    } catch (error) { fail(error); return null; }
  };
  view.dispatch({effects: cm.StateEffect.appendConfig.of([
    cm.intelligenceField,
    cm.ViewPlugin.fromClass(class {
      constructor() {
        this.dismiss = () => signature?.remove();
        addEventListener('scroll', this.dismiss, true);
        addEventListener('resize', this.dismiss);
      }
      destroy() {
        ++generation; clearTimeout(timer); signature?.remove();
        removeEventListener('scroll', this.dismiss, true);
        removeEventListener('resize', this.dismiss);
      }
    }),
    cm.autocompletion({override: [complete]}),
    cm.hoverTooltip(async (_view, pos) => {
      try {
        const state = view.state;
        const value = await request('hover', pos);
        if (view.state !== state || !view.state.field(cm.intelligenceField, false)) return null;
        return value ? {pos, create: () => ({dom: text(value)})} : null;
      } catch (error) { fail(error); return null; }
    }),
    cm.EditorView.updateListener.of(update => {
      if (!update.docChanged && !update.selectionSet && !update.focusChanged) return;
      const epoch = ++generation;
      clearTimeout(timer);
      signature?.remove();
      if (!view.hasFocus) return;
      if (update.focusChanged) worker('ty').catch(fail);
      const position = view.state.selection.main.head;
      const before = view.state.doc.sliceString(Math.max(0, position - 200), position);
      if (!before.includes('(')) return;
      timer = setTimeout(async () => {
        try {
          const value = await request('signature', position);
          if (!value || epoch !== generation || !view.hasFocus) return;
          const coords = view.coordsAtPos(position);
          if (!coords) return;
          signature = text(value.label, false);
          signature.classList.add('documentation-signature');
          signature.style.left = `${Math.max(8, Math.min(coords.left, innerWidth - 320))}px`;
          document.body.append(signature);
          const bounds = signature.getBoundingClientRect();
          signature.style.left = `${Math.max(8, Math.min(coords.left, innerWidth - bounds.width - 8))}px`;
          signature.style.top = `${Math.max(8, Math.min(coords.bottom + 4, innerHeight - bounds.height - 8))}px`;
        } catch (error) { fail(error); }
      }, 180);
    }),
    cm.EditorView.domEventHandlers({keydown(event) {
      if (event.key !== 'Escape' || !signature?.isConnected) return false;
      ++generation;
      signature.remove();
      event.preventDefault();
      event.stopPropagation();
      return true;
    }}),
  ])});
  // Warm only ty on first editing focus, never the scientific worker.
  if (view.hasFocus) worker('ty').catch(fail);
}
async function nativeEditor(event) {
  const editor = event.detail?.editor ?? event.composedPath().find(node => node.matches?.('marimo-code-editor[data-notebook-source]'));
  if (!editor) return;
  const runtime = editor.closest('.live-notebook')?.dataset.notebookRuntime;
  if (!runtime) return;
  const cm = await api(runtime);
  const dom = editor.shadowRoot?.querySelector('.cm-editor');
  const view = dom && cm.EditorView.findFromDOM(dom);
  if (!view) return;
  let notice = editor.parentElement.querySelector('.documentation-type-status');
  if (!notice) {
    notice = document.createElement('span');
    notice.className = 'documentation-type-status';
    notice.setAttribute('role', 'status');
    editor.after(notice);
  }
  await intelligence(view, () => notebookContext(editor, cm), message => notice.textContent = message, cm);
  editor.dataset.documentationTypeReady = 'true';
}
async function renderOutput(output, bundles, append = false) {
  if (!append) {
    for (const frame of output.querySelectorAll('iframe[data-output-id]')) outputFrames.delete(Number(frame.dataset.outputId));
    output.replaceChildren();
  }
  for (const bundle of bundles) {
    if (!Object.keys(bundle).length) continue;
    if (bundle['text/html']) {
      // Use the display owner's complete document, including its real graph
      // scripts and controls, rather than stripping scripts into inert markup.
      const frame = document.createElement('iframe');
      frame.title = 'Python example result';
      frame.setAttribute('sandbox', 'allow-scripts');
      const id = ++serial;
      frame.dataset.outputId = String(id);
      outputFrames.set(id, frame);
      frame.addEventListener('load', syncOutputThemes);
      const bridge = `<script>(() => {
        addEventListener('message', event => {
          if (event.source !== parent || event.data?.type !== 'documentation-output-theme' || event.data.id !== ${id}) return;
          if (!['light', 'dark'].includes(event.data.theme)) return;
          document.documentElement.dataset.theme = event.data.theme;
          document.documentElement.style.colorScheme = event.data.theme;
        });
        const resize = () => {
          const range = document.createRange(); range.selectNodeContents(document.body);
          const height = range.getBoundingClientRect().bottom + parseFloat(getComputedStyle(document.body).marginBottom);
          parent.postMessage({type: 'documentation-output-size', id: ${id}, height}, '*');
        };
        new ResizeObserver(resize).observe(document.body);
        document.fonts.ready.then(resize); resize();
      })();</script>`;
      const html = bundle['text/html'];
      frame.srcdoc = html.includes('</body>') ? html.replace('</body>', bridge + '</body>') : html + bridge;
      output.append(frame);
    } else if (bundle['image/svg+xml']) {
      const image = document.createElement('img');
      image.alt = 'Rendered example';
      image.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(bundle['image/svg+xml']);
      output.append(image);
      await image.decode();
    } else if (bundle['image/png']) {
      const image = document.createElement('img');
      image.alt = 'Python example result';
      image.src = 'data:image/png;base64,' + bundle['image/png'];
      output.append(image);
      await image.decode();
    } else if (bundle['text/latex']) {
      const latex = bundle['text/latex'].replace(/^\${1,2}|\${1,2}$/g, '');
      const send = await worker('typst');
      const result = await send({type: 'run', source: '#set page(width: auto, height: auto, margin: 4pt)\n#import "@preview/mitex:0.2.6": mitex\n#mitex(' + JSON.stringify(latex) + ')'});
      if (!result.success) throw Error('Math rendering failed: ' + JSON.stringify(result.outputs));
      await renderOutput(output, result.outputs, true);
    } else {
      const pre = document.createElement('pre');
      pre.textContent = bundle['text/plain'] || bundle['text/latex'] || '';
      output.append(pre);
    }
  }
}
async function activate(pre, row, cell, all, c) {
  const container = document.createElement('section');
  container.className = 'documentation-cell';
  container.dataset.cellId = cell.id;
  row.before(container);
  container.append(row);
  const controls = document.createElement('div');
  controls.className = 'notebook-cell-header';
  const toggle = document.createElement('button');
  toggle.type = 'button'; toggle.className = 'notebook-source-toggle';
  toggle.setAttribute('aria-expanded', 'false');
  toggle.textContent = '</>'; toggle.title = 'Edit source'; toggle.setAttribute('aria-label', 'Edit source');
  const run = document.createElement('button');
  run.type = 'button'; run.className = 'notebook-run';
  run.textContent = '▷'; run.title = cell.language === 'typst' ? 'Compile example' : 'Run example';
  run.setAttribute('aria-label', run.title);
  const reset = document.createElement('button');
  reset.type = 'button'; reset.textContent = 'Reset';
  const status = document.createElement('span'); status.className = 'notebook-cell-state'; status.setAttribute('role', 'status');
  const timing = document.createElement('span'); timing.className = 'notebook-cell-time';
  timing.title = cell.language === 'python' ? 'Python execution time' : 'Typst compilation and rendering time';
  const hints = document.createElement('span'); hints.className = 'documentation-type-status'; hints.setAttribute('role', 'status');
  controls.append(toggle, run, reset, status, timing);
  const copy = row.querySelector('button[data-copy-target]');
  if (copy) controls.append(copy);
  const output = document.createElement('div'); output.className = 'documentation-cell-output'; output.setAttribute('aria-live', 'polite');
  container.append(controls, hints, output);
  let view;
  async function edit() {
    if (view) {
      row.hidden = !row.hidden;
      toggle.setAttribute('aria-expanded', String(!row.hidden));
      toggle.title = row.hidden ? 'Show source' : 'Hide source';
      toggle.setAttribute('aria-label', toggle.title);
      if (!row.hidden) view.focus();
      return;
    }
    const cm = await api();
    const mount = document.createElement('div'); mount.className = 'documentation-editor';
    pre.after(mount);
    view = new cm.EditorView({parent: mount, state: cm.EditorState.create({doc: cell.execution || cell.source, extensions: [
      cm.lineNumbers(), cm.highlightActiveLine(), cm.outsetGutter, cm.history(),
      ...(cell.language === 'python' ? [cm.python()] : [typstLanguage(cm)]),
      cm.syntaxHighlighting(cm.highlightStyle),
      cm.keymap.of([...cm.historyKeymap, ...cm.completionKeymap, {key: 'Mod-Enter', run: () => {run.click(); return true;}}]),
      cm.EditorView.updateListener.of(update => {
        if (update.docChanged) pre.querySelector('code').textContent = update.state.doc.toString();
      }),
    ]})});
    views.set(cell.id, view);
    toggle.setAttribute('aria-expanded', 'true');
    toggle.title = 'Hide source'; toggle.setAttribute('aria-label', toggle.title);
    pre.hidden = true;
    row.querySelector('.code-line-numbers')?.setAttribute('hidden', '');
    if (cell.language === 'python') intelligence(view,
      () => sourceFor({...cell, source: ''}, all).slice(0, -1).map(([, code]) => code).join('\n\n'),
      message => hints.textContent = message).catch(error => hints.textContent = `Type hints unavailable: ${error.message}`);
    view.focus();
  }
  toggle.addEventListener('click', async () => {
    toggle.disabled = true;
    try { await edit(); }
    catch (error) { status.textContent = error.message; }
    finally { toggle.disabled = false; }
  });
  reset.addEventListener('click', () => {
    if (view) view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: cell.execution || cell.source}});
    renderOutput(output, []); status.textContent = ''; timing.textContent = '';
    delete container.dataset.cellState;
  });
  run.addEventListener('click', async () => {
    run.disabled = true; run.setAttribute('aria-busy', 'true');
    reset.disabled = true;
    container.dataset.cellState = 'running';
    status.textContent = 'Loading…'; timing.textContent = '';
    try {
      const send = await worker(cell.language);
      status.textContent = cell.language === 'python' ? 'Running…' : 'Compiling…';
      const request = cell.language === 'python'
        ? {session: cell.file, sources: sourceFor(cell, all), packages: cell.packages || []}
        : {source: (cell.prefix || '') + (view?.state.doc.toString() ?? cell.execution ?? cell.source), path: cell.path};
      const result = await send({type: 'run', ...request});
      await renderOutput(output, result.outputs);
      container.dataset.cellState = result.success ? 'success' : 'error';
      status.textContent = result.success ? 'Done' : 'Error';
      if (Number.isFinite(result.milliseconds)) timing.textContent = `${(result.milliseconds / 1000).toFixed(3)} s`;
    } catch (error) { status.textContent = error.message; container.dataset.cellState = 'error'; }
    finally { run.disabled = false; reset.disabled = false; run.removeAttribute('aria-busy'); }
  });
}
export async function initialize(configurationURL) {
  initializeNotebookSources();
  if (configurationURL) configLocation = new URL(configurationURL);
  document.addEventListener('focusin', event => nativeEditor(event).catch(console.error));
  document.addEventListener('documentation-editor-ready', event => nativeEditor(event).catch(console.error));
  // Notebook integration above needs no standalone example inventory. In
  // particular, a shell download command must not fetch every Python example.
  const examples = [...document.querySelectorAll('.docs-article pre:has(> code[data-lang])')].filter(pre =>
    !pre.matches('.api-signature, .api-member-signature') &&
    pre.closest('.reference-copy-row') &&
    !pre.closest('.reference-copy-row').classList.contains('code-frame') &&
    ['python', 'py', 'typst', 'typ'].includes(pre.querySelector('code').dataset.lang));
  if (!examples.length) return;
  const c = await config();
  for (const pre of examples) {
    const raw = pre.querySelector('code');
    const row = pre.closest('.reference-copy-row');
    frameSource(pre);
    const source = raw.textContent.trimEnd();
    const member = pre.closest('.api-member')?.querySelector('.api-member-heading code')?.textContent;
    const candidates = c.cells.filter(cell => cell.source === source
      && (!member || !cell.owner || cell.owner.split('.').at(-1) === member));
    const cell = candidates.find(cell => (location.pathname.includes(`/products/${cell.product}/`) || cell.product === 'spenso' && location.pathname.includes('/products/idenso/'))
      && location.pathname.replace(/index\.html$/, '').endsWith('/' + cell.route))
      ?? candidates.find(cell => cell.product && (location.pathname.includes(`/products/${cell.product}/`) || cell.product === 'spenso' && location.pathname.includes('/products/idenso/')) && !cell.route)
      ?? (candidates.length === 1 ? candidates[0] : undefined);
    if (!cell) continue;
    if (cell.reason) {
      const reason = document.createElement('p'); reason.className = 'documentation-cell-reason'; reason.textContent = cell.reason;
      row.after(reason);
    } else await activate(pre, row, cell, c.cells, c);
  }
}

// The full Marimo app owns a separate pinned CodeMirror bundle. Extend its
// existing views, retaining Marimo's state, history, keymaps and reactive graph.
export async function initializeFullNotebook(configurationURL, runtimeURL, payloadURL) {
  initializeNotebookSources();
  configLocation = new URL(configurationURL);
  document.documentElement.dataset.documentationNotebookRuntime = new URL(runtimeURL).href;
  const [cm, payload] = await Promise.all([
    loadEditor(runtimeURL, 'full'), fetch(payloadURL).then(response => response.json()),
  ]);
  const sources = new Map(Object.entries(payload.full_sources));
  const editors = new Map();
  const attach = async () => {
    document.documentElement.dataset.codeTheme = document.body.classList.contains('dark') ? 'dark' : 'light';
    for (const dom of document.querySelectorAll('[data-testid="cell-editor"] .cm-editor')) {
      if (!dom.querySelector('.cm-content[data-language="python"]')) continue;
      const view = cm.EditorView.findFromDOM(dom);
      if (!view) continue;
      const cellId = dom.closest('.marimo-cell')?.dataset.cellId;
      if (!cellId) continue;
      if (!sources.has(cellId)) sources.set(cellId, view.state.doc.toString());
      editors.set(cellId, view);
      dom.closest('[data-testid="cell-editor"]').classList.add('documentation-full-code-frame');
      dom.dataset.codeTheme = view.state.facet(cm.EditorView.darkTheme) ? 'dark' : 'light';
      dom.closest('[data-testid="cell-editor"]').dataset.codeTheme = dom.dataset.codeTheme;
      if (!view.state.field(cm.intelligenceField, false)) {
        view.dispatch({effects: cm.StateEffect.appendConfig.of([cm.outsetGutter, cm.syntaxHighlighting(cm.classHighlighter)])});
        await intelligence(view, () => [...sources].map(([owner, source]) => {
          if (owner === cellId) return '';
          return editors.get(owner)?.state.doc.toString() ?? source;
        }).join('\n\n'), message => { dom.dataset.typeStatus = message; }, cm);
        dom.dataset.documentationTypeReady = 'true';
      }
    }
  };
  let scheduled = false;
  const schedule = () => {
    if (scheduled) return;
    scheduled = true;
    requestAnimationFrame(() => { scheduled = false; attach().catch(console.error); });
  };
  new MutationObserver(records => {
    if (records.some(record => record.type === 'childList' || record.target === document.body || record.target.matches?.('.cm-editor'))) schedule();
  }).observe(document.body, {childList: true, subtree: true, attributes: true, attributeFilter: ['class']});
  document.addEventListener('focusin', schedule);
  await attach();
}
