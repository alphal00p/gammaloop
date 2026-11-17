// Only content-addressed runtime binaries/packages are immutable. Authored
// guides, notebook payloads and small JS chunks keep normal HTTP freshness.
const cacheName = "alphal00p-notebook-runtime-v1";
const runtimePrefix = new URL("assets/notebook-runtime/", self.location.href).pathname;
const filesPrefix = new URL("assets/notebook-runtime-files/", self.location.href).pathname;
const documentationPrefix = new URL("assets/documentation-runtime/", self.location.href).pathname;
self.addEventListener("install", () => self.skipWaiting());
self.addEventListener("activate", event => event.waitUntil(self.clients.claim()));
self.addEventListener("fetch", event => {
  const url = new URL(event.request.url);
  const documentation = url.pathname.startsWith(documentationPrefix);
  const prefix = documentation ? documentationPrefix : url.pathname.startsWith(filesPrefix) ? filesPrefix : runtimePrefix;
  if (event.request.method !== "GET" || url.origin !== self.location.origin ||
      !url.pathname.startsWith(prefix) ||
      !/^[a-f0-9]{64}\//.test(url.pathname.slice(prefix.length)) ||
      (!documentation && !/\.(whl|zip|wasm|gz)$/.test(url.pathname))) return;
  event.respondWith((async () => {
    // CacheStorage is optional: WebKit storage/upgrade failures can leave its
    // promises pending. Bound only the lookup, never the scientific session or
    // network request. A late lookup must not overwrite the network response.
    let deadline;
    const lookup = (async () => {
      const cache = await caches.open(cacheName);
      const cached = await cache.match(event.request);
      if (cached) return {cache, cached};
      // Migrate unchanged Python archives from runtime-wide URLs. Basenames
      // only select candidates: the new URL's SHA-256 must match their bytes.
      // Native responses stay untouched to retain streaming/code-cache state.
      if (prefix === filesPrefix && /\.(whl|zip)$/.test(url.pathname)) {
        const [digest, filename] = url.pathname.slice(prefix.length).split("/");
        for (const key of await cache.keys()) {
          const path = new URL(key.url).pathname;
          if (!path.startsWith(runtimePrefix) || !path.endsWith("/" + filename)) continue;
          const previous = await cache.match(key);
          if (!previous?.ok) continue;
          const bytes = await previous.clone().arrayBuffer();
          const hash = [...new Uint8Array(await crypto.subtle.digest("SHA-256", bytes))]
            .map(byte => byte.toString(16).padStart(2, "0")).join("");
          if (hash === digest) {
            event.waitUntil(cache.put(event.request, previous.clone()).catch(() => {}));
            return {cache, cached: previous};
          }
        }
      }
      return {cache};
    })().catch(() => ({}));
    const {cache, cached} = await Promise.race([lookup, new Promise(resolve => {
      deadline = setTimeout(() => resolve({}), 1000);
    })]);
    clearTimeout(deadline);
    // Preserve the CacheStorage response itself, including its code-cache
    // metadata. Rewrapping it would discard the browser's compiled-code cache.
    if (cached) return cached;
    const native = url.pathname.endsWith("/symbolica-core.wasm");
    const response = await fetch(native ? url.href + ".gz" : event.request);
    if (!response.ok) return response;
    const headers = new Headers(response.headers);
    headers.set("Cache-Control", "public, max-age=31536000, immutable");
    if (native) {
      headers.set("Content-Type", "application/wasm");
      headers.delete("Content-Length");
      headers.delete("Content-Encoding");
    }
    const body = native ? response.body.pipeThrough(new DecompressionStream("gzip")) : response.body;
    const immutable = new Response(body, {status: response.status, statusText: response.statusText, headers});
    if (cache) {
      // Do not wait for the full download and cache write before delivering
      // headers: compileStreaming must consume bytes as they arrive. Later
      // cache hits above still return the stored response and its code metadata.
      event.waitUntil(cache.put(event.request, immutable.clone()).catch(() => {}));
    }
    return immutable;
  })());
});
