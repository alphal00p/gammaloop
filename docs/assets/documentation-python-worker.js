// Direct Pyodide execution. This worker never imports the Marimo frontend/kernel.
let initialized;
let python;
let stage = 'imports';
let chain = Promise.resolve();
let measurementId = 0;
const measurements = new Map();
self.documentationMeasure = request => new Promise((resolve, reject) => {
  const id = ++measurementId;
  measurements.set(id, {resolve, reject});
  postMessage({type: 'measure', id, request: JSON.parse(request)});
});
async function initialize(config) {
  const runtime = new URL(config.runtime, self.location.href);
  const [{loadPyodide}, {mountNative}, lock] = await Promise.all([
    import(new URL('pyodide/pyodide.mjs', config.compilerRoot).href),
    import(new URL('native-loader.js', runtime).href),
    fetch(new URL('pyodide-lock.json', runtime)).then(r => {
      if (!r.ok) throw Error(`Python package lock: HTTP ${r.status}`);
      return r.json();
    }),
  ]);
  for (const item of Object.values(lock.packages)) {
    if (item.file_name.startsWith('__LOCAL_WHEEL__/')) {
      item.file_name = new URL('../../notebook-runtime-files/' + item.file_name.slice(16), runtime).href;
    }
  }
  stage = 'Pyodide';
  python = await loadPyodide({
    indexURL: new URL('pyodide/', config.compilerRoot).href,
    stdLibURL: new URL(config.stdlib, runtime).href,
    lockFileContents: lock,
    packageBaseUrl: new URL('wheels/', runtime).href,
  });
  stage = 'packages';
  await python.loadPackage([new URL(config.wheel, runtime).href, 'ipython', 'marimo']);
  stage = 'native mount';
  await mountNative(python);
  stage = 'filesystem';
  python.FS.mkdirTree('/project');
  for (const [path, content] of Object.entries(config.files)) {
    const filename = '/project/' + path;
    python.FS.mkdirTree(filename.slice(0, filename.lastIndexOf('/')));
    python.FS.writeFile(filename, content);
  }
  python.FS.chdir('/project');
  stage = 'IPython setup';
  await python.runPythonAsync(`
import ast, contextlib, hashlib, inspect, io, json, traceback
import IPython.display
from js import documentationMeasure
from IPython.core.interactiveshell import InteractiveShell
_docs_shell = InteractiveShell.instance()
_docs_shell.display_formatter.active_types = ["text/plain", "text/html", "image/svg+xml", "text/latex", "image/png"]
async def _docs_execute(session, sources, reset):
    namespace = {"__name__": "__main__"}
    # Symbol attributes live in the native registry, beyond Python dictionaries.
    # Keep one stable public Symbolica namespace per documentation page.
    scope = "docs_" + hashlib.sha256(session.encode()).hexdigest()[:20]
    exec("from symbolica import set_namespace\\nset_namespace(" + repr(scope) + ")", namespace)
    outputs = []
    displayed = []
    def display(*objects, **kwargs):
        displayed.extend(objects)
    IPython.display.display = display
    def publish(data, metadata=None, **kwargs):
        outputs.append({key: value for key, value in data.items()
                        if key in ("text/plain", "text/html", "image/svg+xml", "text/latex", "image/png")})
    _docs_shell.display_pub.publish = publish
    stdout, stderr = io.StringIO(), io.StringIO()
    success = True
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        try:
            for filename, source in sources:
                tree = ast.parse(source, filename=filename)
                last = tree.body.pop() if tree.body and isinstance(tree.body[-1], ast.Expr) else None
                pending = eval(compile(tree, filename, "exec", flags=ast.PyCF_ALLOW_TOP_LEVEL_AWAIT), namespace)
                if inspect.isawaitable(pending):
                    await pending
                if last is not None:
                    expression = ast.Expression(last.value)
                    value = eval(compile(expression, filename, "eval", flags=ast.PyCF_ALLOW_TOP_LEVEL_AWAIT), namespace)
                    if inspect.isawaitable(value):
                        value = await value
                    if value is not None:
                        displayed.append(value)
                for value in displayed:
                    if hasattr(value, "measurement_request") and hasattr(value, "measured"):
                        receipt = await documentationMeasure(json.dumps(value.measurement_request()))
                        value.measured(measurements=json.loads(receipt))
                    data, _ = _docs_shell.display_formatter.format(value)
                    publish(data)
                displayed.clear()
        except BaseException:
            success = False
            traceback.print_exc()
    if stdout.getvalue():
        outputs.insert(0, {"text/plain": stdout.getvalue()})
    if stderr.getvalue():
        outputs.append({"text/plain": stderr.getvalue()})
    return json.dumps({"success": success, "outputs": outputs})

`);
}
self.onmessage = ({data}) => {
  if (data.type === 'measurement') {
    const pending = measurements.get(data.id);
    if (pending) {
      measurements.delete(data.id);
      data.error ? pending.reject(Error(data.error)) : pending.resolve(JSON.stringify(data.receipt));
    }
    return;
  }
  if (data.type === 'init') {
    initialized = initialize(data.config);
    initialized.then(() => postMessage({type: 'ready'}), error => postMessage({type: 'error', message: stage + ': ' + String(error) + '\n' + (error.stack || JSON.stringify(error))}));
    return;
  }
  if (data.type !== 'run') return;
  chain = chain.catch(() => {}).then(async () => {
    try {
      await initialized;
      await python.loadPackage(data.packages || []);
      await python.loadPackagesFromImports(data.sources.map(([, source]) => source).join('\n'));
      // A requested cell executes in a fresh namespace with its explicitly
      // declared prerequisites. Edits to prerequisites are supplied by the host.
      const start = performance.now();
      python.globals.set('_docs_request', JSON.stringify(data));
      const result = await python.runPythonAsync(`
_docs_data = json.loads(_docs_request)
await _docs_execute(_docs_data["session"], _docs_data["sources"], True)
`);
      postMessage({type: 'result', id: data.id, milliseconds: performance.now() - start, ...JSON.parse(result)});
    } catch (error) {
      postMessage({type: 'result', id: data.id, success: false, outputs: [{'text/plain': String(error)}]});
    }
  });
};
