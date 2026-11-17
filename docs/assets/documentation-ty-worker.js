// ty's WASM Workspace API, independent of both Python execution workers.
let initialized;
let workspace;
let document;
let Position;
let current = '';
async function initialize(config) {
  const root = new URL(config.compilerRoot);
  const [ty, response, stubs] = await Promise.all([
    import(new URL('ty/ty_wasm.js', root).href),
    fetch(new URL('ty/ty_wasm_bg.wasm', root)),
    fetch(new URL('type-stubs.json', root)).then(r => {
      if (!r.ok) throw Error(`Type stubs: HTTP ${r.status}`);
      return r.json();
    }),
  ]);
  if (!response.ok) throw Error(`ty: HTTP ${response.status}`);
  await ty.default({module_or_path: response});
  Position = ty.Position;
  workspace = new ty.Workspace('/project', ty.PositionEncoding.Utf16, {environment: {'python-version': '3.14'}});
  for (const [path, source] of Object.entries(stubs)) workspace.openFile('/project/' + path, source).free();
  document = workspace.openFile('/project/cell.py', '');
  return ty.version();
}
self.onmessage = async ({data}) => {
  if (data.type === 'init') {
    initialized = initialize(data.config);
    initialized.then(version => postMessage({type: 'ready', version}), error => postMessage({type: 'error', message: String(error)}));
    return;
  }
  try {
    await initialized;
    if (data.source !== current) {
      workspace.closeFile(document);
      document = workspace.openFile('/project/cell.py', data.source);
      current = data.source;
    }
    const position = new Position(data.line + 1, data.column + 1);
    let result;
    try {
      if (data.type === 'complete') {
        const items = workspace.completions(document, position);
        try {
          result = items.map(item => ({label: item.name, insertText: item.insert_text,
            detail: item.detail, documentation: item.documentation, kind: item.kind}));
        } finally { items.forEach(item => item.free()); }
      } else if (data.type === 'hover') {
        const hover = workspace.hover(document, position);
        try { result = hover ? hover.markdown : null; }
        finally { hover?.free(); }
      } else if (data.type === 'signature') {
        const help = workspace.signatureHelp(document, position);
        if (help) {
          const signatures = help.signatures;
          try {
            const selected = signatures[help.active_signature || 0];
            result = selected ? {label: selected.label, documentation: selected.documentation,
              activeParameter: selected.active_parameter} : null;
          } finally { signatures.forEach(item => item.free()); help.free(); }
        }
      }
    } finally { /* Workspace methods consume Position. */ }
    postMessage({id: data.id, type: data.type, result});
  } catch (error) {
    postMessage({id: data.id, type: 'error', message: String(error)});
  }
};
