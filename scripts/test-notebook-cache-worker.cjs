// Run in the Nix development shell with Node's native Fetch/Web Streams APIs.
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const fs = require('node:fs');
const vm = require('node:vm');
const source = fs.readFileSync('crates/linnet-py/examples/notebook-cache-worker.js', 'utf8');
const origin = 'https://preview.example';
const bytes = 'verified archive';
const digest = crypto.createHash('sha256').update(bytes).digest('hex');
const previousURL = `${origin}/assets/notebook-runtime/${'a'.repeat(64)}/wheels/package.whl`;
const nextURL = `${origin}/assets/notebook-runtime-files/${digest}/package.whl`;

async function trial({fault, previous, networkFailure, bodyFailure, stored} = {}) {
  const handlers = {}, writes = [], requests = [];
  const old = new Response(previous ?? bytes);
  const cache = {
    keys: async () => [new Request(previousURL)],
    match: async request => {
      if (fault === 'match') return new Promise(() => {});
      const url = request.url;
      return url === previousURL ? old : stored;
    },
    put: async (request, response) => { writes.push({url: request.url, body: await response.text()}); },
  };
  let finish;
  const body = new ReadableStream({start(controller) { finish = () => {
    if (bodyFailure) {controller.error(new Error('interrupted body')); return;}
    controller.enqueue(new TextEncoder().encode(bytes)); controller.close();
  }; }});
  const scope = {URL, Headers, Response, Request, DecompressionStream, Uint8Array,
    crypto: crypto.webcrypto, setTimeout, clearTimeout,
    caches: {open: () => fault === 'open' ? new Promise(() => {}) : Promise.resolve(cache)},
    fetch: async request => {
      requests.push(String(request?.url || request));
      if (networkFailure) throw new Error('interrupted download');
      return new Response(body);
    },
    self: {location: {href: origin + '/notebook-cache-worker.js', origin},
      addEventListener: (name, fn) => {handlers[name] = fn;}},
  };
  vm.runInNewContext(source, scope);
  let promise;
  const tasks = [];
  const start = performance.now();
  handlers.fetch({request: new Request(nextURL), respondWith: value => {promise = value;},
    waitUntil: value => tasks.push(value)});
  if (networkFailure) {
    await assert.rejects(promise, /interrupted download/);
    assert.equal(writes.length, 0);
    return;
  }
  const response = await promise;
  if (stored) assert.equal(response, stored, 'return original cached response');
  if (fault) assert.ok(performance.now() - start < 1400, 'bounded cache lookup');
  // Headers resolve before the network body completes; no startup buffering.
  finish();
  if (bodyFailure) {
    await assert.rejects(response.text(), /interrupted body/);
    await Promise.all(tasks);
    assert.equal(writes.length, 0, 'partial bodies are never cached');
    return;
  }
  assert.equal(await response.text(), stored ? 'cached' : bytes);
  await Promise.all(tasks);
  return {requests, writes};
}

(async () => {
  const migrated = await trial();
  assert.equal(migrated.requests.length, 0, 'verified old archive reused without network');
  assert.equal(migrated.writes[0].url, nextURL);
  const corrupt = await trial({previous: 'same basename, different bytes'});
  assert.deepEqual(corrupt.requests, [nextURL], 'mismatched hash must fetch');
  for (const fault of ['open', 'match']) {
    const recovered = await trial({fault});
    assert.deepEqual(recovered.requests, [nextURL]);
  }
  await trial({stored: new Response('cached')});
  await trial({previous: 'mismatch', networkFailure: true});
  await trial({previous: 'mismatch', bodyFailure: true});
  const recovery = await trial({previous: 'mismatch'});
  assert.equal(recovery.requests.length, 1);
  console.log('Verified migration, hash mismatch, bounded lookup, streaming, original cache response, failure/retry');
})().catch(error => {console.error(error); process.exitCode = 1;});
