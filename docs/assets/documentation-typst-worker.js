// The compiler and vector renderer are real, separately pinned upstream WASM.
let initialized;
let compiler;
let renderer;
let chain = Promise.resolve();
async function initialize(config) {
  const root = new URL(config.compilerRoot);
  const [c, r, response] = await Promise.all([
    import(new URL('compiler/typst_ts_web_compiler.mjs', root).href),
    import(new URL('renderer/typst_ts_renderer.mjs', root).href),
    fetch(new URL('typst-files.json', root)),
  ]);
  if (!response.ok) throw Error(`Typst filesystem: HTTP ${response.status}`);
  const files = await response.json();
  await Promise.all([
    c.default({module_or_path: fetch(new URL('compiler/typst_ts_web_compiler_bg.wasm', root))}),
    r.default({module_or_path: fetch(new URL('renderer/typst_ts_renderer_bg.wasm', root))}),
  ]);
  const builder = new c.TypstCompilerBuilder();
  builder.set_dummy_access_model();
  await builder.set_package_registry({}, spec => `/project/packages/${spec.namespace}/${spec.name}/${spec.version}`);
  for (const font of config.fonts) {
    const response = await fetch(new URL(font, root));
    if (!response.ok) throw Error(`Typst font: HTTP ${response.status}`);
    await builder.add_raw_font(new Uint8Array(await response.arrayBuffer()));
  }
  compiler = await builder.build(); // consumes builder
  for (const [path, value] of Object.entries(files)) {
    compiler.map_shadow('/project/' + path, Uint8Array.from(atob(value), c => c.charCodeAt(0)));
  }
  renderer = new r.TypstRenderer();
}
self.onmessage = ({data}) => {
  if (data.type === 'init') {
    initialized = initialize(data.config);
    initialized.then(() => postMessage({type: 'ready'}), error => postMessage({type: 'error', message: String(error)}));
    return;
  }
  if (data.type !== 'run') return;
  chain = chain.catch(() => {}).then(async () => {
    try {
      await initialized;
      const start = performance.now();
      const path = '/project/' + (data.path || 'cell.typ');
      compiler.add_source(path, data.source);
      const world = compiler.snapshot('/project', path, []);
      try {
        const result = world.get_artifact(0, 1);
        if (!result.result) throw Error(JSON.stringify(result));
        const session = renderer.session_from_artifact(result.result, 'vector');
        try {
          postMessage({type: 'result', id: data.id, success: true,
            milliseconds: performance.now() - start,
            // Body, definitions and CSS form the static image. The renderer's
            // HTML-only interaction script is not valid inside an SVG image.
            outputs: [{'image/svg+xml': renderer.svg_data(session, 0b0111)}], diagnostics: result.diagnostics});
        } finally { session.free(); }
      } finally { world.free(); }
    } catch (error) {
      postMessage({type: 'result', id: data.id, success: false, outputs: [{'text/plain': String(error)}]});
    }
  });
};
