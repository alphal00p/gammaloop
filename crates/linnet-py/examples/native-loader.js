// Compile while Pyodide initializes. The full wheel remains available to the
// standalone editor; this loader preserves its native module, data and RECORD.
const native = __NATIVE_DESCRIPTOR__;
const binaryURL = new URL("__NATIVE_URL__", import.meta.url);
// Preserve a stable response URL so the browser can maintain its own compiled
// code cache. The cache worker supplies gzip-decoded bytes on static hosts.
const compileStart = performance.now();
const compiled = WebAssembly.compileStreaming(fetch(binaryURL)).then(module => {
  performance.measure("notebook:native-compile", {start: compileStart});
  return {module};
}, error => ({error}));

export async function mountNative(pyodide) {
  const result = await compiled;
  if (result.error) throw result.error;
  const module = pyodide._module;
  // Pyodide's dynamic-linking interface is version-specific. Fail explicitly
  // when the pinned runtime changes instead of silently losing native features.
  if (pyodide.version !== "314.0.0" || !module?.preloadedWasm ||
      typeof module.loadWebAssemblyModule !== "function" ||
      typeof pyodide._api?.loadDynlib !== "function") {
    throw new Error("Review the Symbolica loader for this Pyodide version");
  }
  const magic = pyodide.runPython("__import__('importlib.util', fromlist=['MAGIC_NUMBER']).MAGIC_NUMBER.hex()");
  if (magic !== native.bytecode_magic) throw new Error("Runtime bytecode does not match this Python interpreter");
  const site = pyodide.runPython("__import__('sysconfig').get_path('platlib')").replace(/^\/+/, "/");
  const path = `${site}/${native.path}`;
  if (pyodide.FS.analyzePath(path).exists) throw new Error("Support wheel contains the native binary");
  module.preloadedWasm[path] = await module.loadWebAssemblyModule(
    result.module, {loadAsync: true, nodelete: true}, path, {},
  );
  await pyodide._api.loadDynlib(path);
  // Opening the binary explicitly remains supported without retaining a second
  // copy during ordinary imports. Native linking above uses the compiled module.
  const slash = path.lastIndexOf("/");
  const node = pyodide.FS.createLazyFile(path.slice(0, slash), path.slice(slash + 1), binaryURL.href, true, false);
  let bytes;
  const read = () => {
    if (!bytes) {
      const request = new XMLHttpRequest();
      request.open("GET", binaryURL.href, false);
      request.responseType = "arraybuffer";
      request.send();
      if (request.status !== 200 || request.response?.byteLength !== native.size) throw new Error("Cannot read Symbolica binary");
      bytes = new Uint8Array(request.response);
    }
    return bytes;
  };
  node.contents = {length: native.size, get: index => read()[index]};
  node.stream_ops.read = (stream, buffer, offset, length, position) => {
    const count = Math.max(0, Math.min(native.size - position, length));
    if (count) buffer.set(read().subarray(position, position + count), offset);
    return count;
  };
  const record = `${site}/${native.record}`;
  pyodide.FS.writeFile(record, pyodide.FS.readFile(record, {encoding: "utf8"}) + native.entry);
}
