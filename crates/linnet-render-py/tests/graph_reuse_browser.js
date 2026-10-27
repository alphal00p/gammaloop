/* Browser regressions for the production module. Native producer requests,
   bundled STIX bytes and native captured SVGs are the only geometry inputs. */
const assert = (condition, message) => {
  if (!condition) throw new Error(message);
};
const equal = (actual, expected, message) =>
  assert(JSON.stringify(actual) === JSON.stringify(expected),
    `${message}: ${JSON.stringify({actual, expected})}`);
const digest = async value => Array.from(new Uint8Array(await crypto.subtle.digest(
  'SHA-256', new TextEncoder().encode(value))), byte => byte.toString(16).padStart(2, '0')).join('');

// Observe browser-owned objects; production exposes no test-only cache API.
class Audit {
  frames = [];
  urls = [];
  fontLoads = [];
  restores = [];
  observers = [];
  hostDocuments = new Set();

  install(win) {
    const doc = win.document;
    this.hostDocuments.add(doc);
    const create = win.URL.createObjectURL;
    const revoke = win.URL.revokeObjectURL;
    win.URL.createObjectURL = blob => {
      const url = create.call(win.URL, blob);
      this.urls.push({url, type: blob.type, document: doc, revoked: 0,
        text: blob.type.includes('css') ? blob.text() : null});
      return url;
    };
    win.URL.revokeObjectURL = url => {
      const entry = this.urls.find(item => item.url === url);
      if (entry) entry.revoked += 1;
      return revoke.call(win.URL, url);
    };
    this.restores.push(() => {
      win.URL.createObjectURL = create;
      win.URL.revokeObjectURL = revoke;
    });
    const observe = target => {
      const observer = new win.MutationObserver(records => {
        for (const record of records) for (const node of record.addedNodes) {
          if (node.nodeType !== 1) continue;
          const frames = node.localName === 'iframe' ? [node] : node.querySelectorAll('iframe');
          for (const frame of frames) {
            const entry = {frame, document: doc, probe: frame.style.visibility === 'hidden', rows: []};
            this.frames.push(entry);
          }
        }
      });
      observer.observe(target, {childList: true, subtree: true});
      this.observers.push(observer);
    };
    observe(doc.documentElement);
    const originalCreate = doc.createElement;
    doc.createElement = function (...args) {
      const element = originalCreate.apply(this, args);
      if (String(args[0]).toLowerCase() === 'iframe') {
        element.addEventListener('load', () => {
          const child = element.contentDocument;
          if (!child) return;
          const originalLoad = child.fonts.load;
          child.fonts.load = (...parameters) => {
            const call = {document: child, parameters, loaded: false};
            thisAudit.fontLoads.push(call);
            return originalLoad.apply(child.fonts, parameters).then(faces => {
              call.loaded = faces.length > 0 && faces.every(face => face.status === 'loaded');
              return faces;
            });
          };
          const rows = new win.MutationObserver(records => {
            const entry = thisAudit.frames.find(item => item.frame === element);
            if (!entry?.probe) return;
            for (const record of records) for (const node of record.addedNodes) {
              if (node.nodeType !== 1) continue;
              const math = node.localName === 'math' ? [node] : node.querySelectorAll('math');
              for (const label of math) entry.rows.push(label.innerHTML);
            }
          });
          rows.observe(child.documentElement, {childList: true, subtree: true});
          thisAudit.observers.push(rows);
        });
      }
      return element;
    };
    const thisAudit = this;
    this.restores.push(() => { doc.createElement = originalCreate; });
  }

  get probes() { return this.frames.filter(entry => entry.probe); }
  get resources() { return this.urls.filter(entry => !entry.type.includes('javascript')); }
  get live() { return this.resources.filter(entry => entry.revoked === 0); }
  close() {
    for (const observer of this.observers) observer.disconnect();
    for (const restore of this.restores.reverse()) restore();
  }
}

class Model {
  callbacks = new Map();
  messages = [];
  events = [];

  constructor(request, css, fontReply = true) {
    this.css = css;
    this.fontReply = fontReply;
    this.fields = {request: structuredClone(request), rendered: {}, font: {}};
    delete this.fields.request.font_css;
    this.el = document.createElement('section');
    this.el.style.width = '350px';
    document.body.append(this.el);
    this.el.addEventListener('linnet:rendered', event => this.events.push(event.detail));
  }

  get = name => this.fields[name];
  on = (name, callback) => {
    if (!this.callbacks.has(name)) this.callbacks.set(name, new Set());
    this.callbacks.get(name).add(callback);
  };
  off = (name, callback) => this.callbacks.get(name)?.delete(callback);
  publish(name, value) {
    this.fields[name] = value;
    for (const callback of this.callbacks.get(`change:${name}`) || []) callback();
  }
  send = message => {
    this.messages.push(message);
    if (message.action === 'font' && this.fontReply) {
      assert(message.key === this.get('request').font_key, 'font request crossed resource identity');
      // Synchronous publication detects send-before-listen races.
      assert(this.callbacks.get('change:font')?.size > 0, 'font listener installed after send');
      this.publish('font', {key: message.key, css: this.css});
    }
  };
  count(action) { return this.messages.filter(message => message.action === action).length; }
  get receipt() { return this.messages.find(message => message.action === 'measure')?.measurements; }
  async mount(module) { this.dispose = await module.default.render({model: this, el: this.el}); }
  close() {
    this.dispose?.();
    this.dispose = null;
    assert([...this.callbacks.values()].every(callbacks => callbacks.size === 0), 'widget listeners leaked');
    assert(!this.el.querySelector('iframe'), 'widget iframe leaked');
    this.el.remove();
  }
}

const waitFor = async (predicate, message) => {
  const deadline = performance.now() + 16000;
  while (!predicate()) {
    assert(performance.now() < deadline, message);
    await new Promise(resolve => setTimeout(resolve, 5));
  }
};
const receipt = (request, result) => {
  equal(result.id, request.id, 'receipt request identity');
  equal(result.metrics.map(metric => metric.id), request.labels.map(label => label.id), 'receipt label order');
  for (const metric of result.metrics) {
    equal(Object.keys(metric).sort(), ['baseline', 'height', 'id', 'width'], 'numeric-only metric');
    assert(Number.isFinite(metric.width) && metric.width > 0, 'invalid actual width');
    assert(Number.isFinite(metric.height) && metric.height > 0, 'invalid actual height');
    assert(Number.isFinite(metric.baseline) && metric.baseline >= 0 && metric.baseline <= metric.height,
      'invalid actual baseline');
  }
};
const numbers = metric => ({width: metric.width, height: metric.height, baseline: metric.baseline});
const serializedMathML = markup => {
  const math = document.createElementNS('http://www.w3.org/1998/Math/MathML', 'math');
  math.innerHTML = markup;
  return math.innerHTML;
};
const variant = (request, name, labels = request.labels) => ({...structuredClone(request), id: name,
  labels: labels.map((label, index) => ({...label, id: 10000 + index}))});

export async function run({name, request, wrapper, engine, width}) {
  const audit = new Audit();
  audit.install(window);
  const wrappers = [];
  const models = [];
  const importWrapper = async suffix => {
    const url = URL.createObjectURL(new Blob([wrapper, `\n// fresh AnyWidget wrapper ${suffix}\n`],
      {type: 'text/javascript'}));
    wrappers.push(url);
    return import(url);
  };
  const module = await importWrapper('first');
  const other = await importWrapper('second');
  const probeDelta = async (count, action, message) => {
    const before = audit.probes.length;
    const result = await action();
    await Promise.resolve();
    equal(audit.probes.length - before, count, message);
    return result;
  };
  const model = (input = request, css = input.font_css, reply = true) => {
    const value = new Model(input, css, reply);
    models.push(value);
    return value;
  };
  try {
    assert(request.font_key && request.font_css.includes('data:font/woff2;base64,'), 'missing native STIX input');
    equal(await digest(request.font_css), request.font_key, 'native CSS SHA256');
    assert(module.measure === other.measure && module.default === other.default,
      'fresh wrappers did not share the canonical raw module');
    let result;
    if (name === 'reuse') {
      const empty = await probeDelta(0, () => module.measure({id: 'empty-public', labels: []}), 'empty public request must not resolve a font');
      equal(empty, {id: 'empty-public', metrics: []}, 'empty public receipt');
      equal(audit.resources.length, 0, 'empty public request allocated font URLs');
      const first = await probeDelta(1, () => module.measure(request), 'cold batch needs one probe iframe');
      receipt(request, first);
      const reordered = variant(request, 'fresh-reordered', [...request.labels].reverse());
      const second = await probeDelta(0, () => other.measure(reordered), 'fresh IDs/order must reuse numeric metrics');
      receipt(reordered, second);
      equal(second.metrics.map(numbers), [...first.metrics].reverse().map(numbers), 'reordered exact metric values');
      second.metrics[0].width = -123;
      const third = await probeDelta(0, () => module.measure(reordered), 'receipt mutation must not poison cache');
      equal(third.metrics.map(numbers), [...first.metrics].reverse().map(numbers), 'numeric cache owns no mutable receipt');
      const changed = {...request.labels[0], mathml: '<mfrac><mi>λ</mi><mrow><mi>α</mi><mo>+</mo><mi>β</mi></mrow></mfrac>'};
      const partial = variant(request, 'partially-new', [request.labels[1], changed, request.labels[0]]);
      const mixed = await probeDelta(1, () => module.measure(partial), 'partial miss needs one batch');
      receipt(partial, mixed);
      equal(audit.probes.at(-1).rows, [serializedMathML(changed.mathml)], 'cached labels must stay out of cold batch');
      equal(numbers(mixed.metrics[0]), numbers(first.metrics[1]), 'partial hit geometry');
      result = {cold: first, reordered: third, partial: mixed};
    } else if (name === 'concurrent') {
      const inputs = [request, variant(request, 'simultaneous-new-IDs'),
        variant(request, 'simultaneous-subset', [request.labels[1], request.labels[0], request.labels[1]])];
      const outputs = await probeDelta(1, () => Promise.all(inputs.map((input, index) =>
        (index ? other : module).measure(input))), 'same-key in-flight requests must share one probe batch');
      outputs.forEach((output, index) => receipt(inputs[index], output));
      equal(outputs[0].metrics.map(numbers), outputs[1].metrics.map(numbers), 'concurrent exact geometry');
      equal(outputs[2].metrics.map(numbers), [outputs[0].metrics[1], outputs[0].metrics[0], outputs[0].metrics[1]].map(numbers),
        'overlapping in-flight label geometry');
      result = {receipts: outputs, probe_rows: audit.probes[0].rows};
    } else if (name === 'isolation') {
      const original = await module.measure(request);
      const changed = {...variant(request, 'style-change'), font_style: {...request.font_style, fontSize: '13px'}};
      const larger = await probeDelta(1, () => module.measure(changed), 'complete font style must isolate cache');
      assert(larger.metrics.some((metric, index) => metric.height > original.metrics[index].height), 'style change had no actual geometry effect');
      const css = request.font_css + '\n/* exact distinct native font resource */';
      const keyed = {...variant(request, 'key-change'), font_css: css, font_key: await digest(css)};
      const keyedResult = await probeDelta(1, () => module.measure(keyed), 'different exact font key must miss');
      equal(keyedResult.metrics.map(numbers), original.metrics.map(numbers), 'equivalent STIX bytes changed geometry');
      const coordinate = {...variant(request, 'coordinate-change'), coordinate_system: 'different-test-coordinate-contract'};
      const before = audit.probes.length;
      let rejected = false;
      try { await module.measure(coordinate); } catch { rejected = true; }
      assert(rejected || audit.probes.length === before + 1, 'coordinate contract silently reused old metrics');
      const child = document.createElement('iframe');
      const loaded = new Promise(resolve => child.addEventListener('load', resolve, {once: true}));
      child.srcdoc = '<!doctype html><body></body>';
      document.body.append(child);
      await loaded;
      audit.install(child.contentWindow);
      try {
        const isolated = model(variant(request, 'same-module-another-document'));
        child.contentDocument.body.append(isolated.el);
        await probeDelta(1, () => isolated.mount(module), 'same module in another Document must probe independently');
        receipt(isolated.get('request'), isolated.receipt);
        equal(isolated.receipt.metrics.map(numbers), original.metrics.map(numbers), 'document isolation geometry');
        isolated.close();
      } finally { child.remove(); }
      result = {original, larger, keyed: keyedResult, coordinate_rejected: rejected};
    } else if (name === 'font-protocol') {
      const pending = model(request, request.font_css, false);
      const mounting = pending.mount(module);
      await waitFor(() => pending.count('font') === 1, 'font request never sent');
      equal(audit.probes.length, 0, 'font retrieval must precede probe iframe');
      pending.publish('font', {key: 'wrong-key', css: request.font_css});
      await new Promise(resolve => setTimeout(resolve, 20));
      assert(!pending.receipt && pending.callbacks.get('change:font')?.size > 0,
        'wrong font key was accepted or matching listener discarded');
      pending.publish('font', {key: request.font_key, css: request.font_css});
      await mounting;
      receipt(request, pending.receipt);
      equal(pending.callbacks.get('change:font')?.size ?? 0, 0, 'font reply listener leaked');
      const fresh = model(variant(request, 'fresh-private-widget'));
      await probeDelta(0, () => fresh.mount(other), 'fresh widget must reuse canonical module metrics');
      receipt(fresh.get('request'), fresh.receipt);
      equal(fresh.count('font'), 0, 'warm view retransmitted font CSS');
      equal(pending.count('font'), 1, 'cold view requested CSS more than once');
      const cancelledCSS = request.font_css + '\n/* disposed font request */';
      const cancelledInput = {...variant(request, 'disposed-before-font'), font_css: cancelledCSS,
        font_key: await digest(cancelledCSS)};
      const cancelled = model(cancelledInput, cancelledCSS, false);
      cancelled.fields.rendered = {id: cancelledInput.id, html: request.completed_html.replaceAll(
        `<style data-linnet-font="true"><![CDATA[${request.font_css}]]></style>`, '')};
      await cancelled.mount(other);
      await waitFor(() => cancelled.count('font') === 1, 'completed snapshot never requested its font');
      cancelled.dispose();
      cancelled.dispose = null;
      await waitFor(() => (cancelled.callbacks.get('change:font')?.size ?? 0) === 0,
        'disposed font retrieval kept its reply listener');
      cancelled.close();
      const retried = model({...cancelledInput, id: 'retry-disposed-font'});
      await retried.mount(module);
      receipt(retried.get('request'), retried.receipt);
      equal(retried.count('font'), 1, 'aborted font entry prevented fresh retrieval');
      result = {font_actions: pending.count('font') + fresh.count('font'), receipts: [pending.receipt, fresh.receipt]};
    } else if (name === 'failed-font') {
      const brokenCSS = request.font_css.replace(/data:font\/woff2;base64,[A-Za-z0-9+/=]+/, 'data:font/woff2;base64,AAAA');
      assert(brokenCSS !== request.font_css, 'font failure control did not alter real font bytes');
      const broken = {...variant(request, 'broken-font'), font_css: brokenCSS, font_key: await digest(brokenCSS)};
      let failures = 0;
      for (let attempt = 0; attempt < 2; attempt++) {
        try {
          await module.measure({...broken, id: `broken-${attempt}`});
        } catch (error) {
          assert(/font|network|load/i.test(String(error)), `unexpected font failure: ${error}`);
          failures += 1;
        }
      }
      equal(failures, 2, 'invalid supplied font unexpectedly measured');
      equal(audit.probes.length, 2, 'failed in-flight font metrics must be retried in a fresh probe');
      const recovered = await probeDelta(1, () => module.measure(request), 'valid font must recover after failed resource');
      receipt(request, recovered);
      await probeDelta(0, () => other.measure(variant(request, 'recovered-warm')), 'successful retry must populate metrics');
      assert(audit.probes.every(entry => !entry.frame.isConnected), 'failed probe iframe leaked');
      result = {failures, recovered};
    } else if (name === 'native-error-lease') {
      const failed = model(request);
      await failed.mount(module);
      receipt(request, failed.receipt);
      const failedURLs = [...audit.resources];
      assert(failedURLs.length > 0 && failedURLs.every(entry => entry.revoked === 0),
        'successful initial measurement did not retain a real font resource');
      failed.publish('rendered', {id: request.id, error: 'native render failed'});
      assert(failed.el.querySelector('[role="alert"]')?.textContent.includes('native render failed'),
        'backend failure was not shown');
      assert(!failed.el.querySelector('iframe') && failed.el.isConnected, 'failed widget must remain attached without a final frame');
      for (let index = 0; index < 4; index++) {
        const css = request.font_css + `\n/* after-native-error-${index} */`;
        const input = {...variant(request, `after-native-error-${index}`), font_css: css, font_key: await digest(css)};
        const measured = await module.measure(input);
        receipt(input, measured);
        equal(measured.metrics.map(numbers), failed.receipt.metrics.map(numbers), 'font resource churn changed real geometry');
      }
      assert(failed.el.isConnected && failed.callbacks.get('change:rendered')?.size === 1,
        'test detached failed widget instead of exercising lease release');
      assert(failedURLs.every(entry => entry.revoked === 1), 'backend error pinned an unused font lease in an attached widget');
      const retried = model(variant(request, 'after-native-error-retry'));
      await probeDelta(0, () => retried.mount(other), 'native failure must preserve valid numeric metrics');
      receipt(retried.get('request'), retried.receipt);
      equal(retried.count('font'), 1, 'evicted failed-view font did not retrieve again');
      result = {attached_failure_released: true, evicted_urls: failedURLs.length, retry_reused_metrics: true};
    } else if (name === 'metric-limit') {
      const labels = Array.from({length: 2049}, (_, id) => ({id,
        mathml: `<mrow><mi>α</mi><!--metric-entry-${id}--></mrow>`}));
      const batch = {...request, id: 'entry-capacity', labels};
      const measured = await probeDelta(1, () => module.measure(batch), 'entry-capacity cold batch');
      receipt(batch, measured);
      const recent = variant(request, 'recent-entry', [labels.at(-1)]);
      await probeDelta(0, () => module.measure(recent), 'most recent of 2049 metrics must remain');
      const oldest = variant(request, 'oldest-entry', [labels[0]]);
      const evicted = await probeDelta(1, () => module.measure(oldest), '2048-entry LRU must evict oldest');
      equal(numbers(evicted.metrics[0]), numbers(measured.metrics[0]), 'entry eviction changed actual geometry');
      result = {entries: labels.length, retained_recent: true, oldest_remeasured: true};
    } else if (name === 'metric-bytes') {
      // Unicode comments exercise UTF-8 accounting without changing ink or
      // manufacturing geometry. 3*~900KiB exceeds 2MiB but not by character count.
      const labels = Array.from({length: 3}, (_, id) => ({id,
        mathml: `<mrow><mi>α</mi><!--${'λ'.repeat(450000)}-${id}--></mrow>`}));
      const measured = await module.measure({...request, id: 'byte-capacity', labels});
      await probeDelta(0, () => module.measure(variant(request, 'byte-recent', [labels[2]])), 'recent byte-budget entry retained');
      const oldest = await probeDelta(1, () => module.measure(variant(request, 'byte-oldest', [labels[0]])), 'UTF-8 2MiB LRU eviction');
      equal(numbers(oldest.metrics[0]), numbers(measured.metrics[0]), 'byte eviction changed actual geometry');
      const oversized = {id: 77, mathml: `<mrow><mi>β</mi><!--${'λ'.repeat(1048577)}--></mrow>`};
      const input = variant(request, 'oversized', [oversized]);
      const first = await probeDelta(1, () => module.measure(input), 'oversized label must still measure');
      const second = await probeDelta(1, () => module.measure({...input, id: 'oversized-again'}), 'oversized label must not be retained');
      equal(first.metrics.map(numbers), second.metrics.map(numbers), 'oversized actual geometry');
      result = {utf8_bytes: labels.map(label => new TextEncoder().encode(label.mathml).length), oversized_measured_twice: true};
    } else if (name === 'native') {
      const {inspect} = await import('/harness.js');
      const manifest = await (await fetch('/manifest')).json();
      const records = [];
      for (const slot of ['a', 'b']) {
        const job = `${engine}-${width}-${slot}`;
        const input = manifest.requests[job];
        const record = model(input);
        record.el.style.width = `${width}px`;
        const previous = new Set(records.flatMap(item => item.input.labels.map(label => label.mathml)));
        const missing = [...new Set(input.labels.map(label => label.mathml))].filter(markup => !previous.has(markup));
        await probeDelta(missing.length ? 1 : 0, () => record.mount(slot === 'a' ? module : other),
          `native fresh ${slot} probe count`);
        equal(audit.probes.at(-1).rows, missing.map(serializedMathML), 'native partial batch must contain only uncached exact markup');
        receipt(input, record.receipt);
        records.push({job, input, record});
      }
      assert(records[0].input.id !== records[1].input.id, 'producer reused a native view');
      equal(records.map(({record}) => record.count('font')), [1, 0], 'private widget font transport reuse');
      const captures = [];
      for (const {job, input, record} of [...records].reverse()) {
        const response = await fetch(`/receipt/${job}`, {method: 'POST', headers: {'Content-Type': 'application/json'},
          body: JSON.stringify(record.receipt)});
        const captured = await response.json();
        assert(response.ok && captured.native, `native completion failed: ${JSON.stringify(captured)}`);
        const block = `<style data-linnet-font="true"><![CDATA[${input.font_css}]]></style>`;
        assert(captured.html.includes(block), 'native font-owned block missing');
        const notebook = captured.html.replaceAll(block, '');
        assert(!notebook.includes(input.font_css), 'notebook transport still embeds full font CSS');
        record.publish('rendered', {id: 'stale-view', html: notebook});
        assert(!record.el.querySelector('iframe'), 'stale view completion mounted');
        record.publish('rendered', {id: input.id, html: notebook});
        await waitFor(() => record.events.length === 1, 'native view never emitted rendered');
        const doc = record.events[0].document;
        const beforeFonts = audit.fontLoads.length;
        const geometry = await inspect(doc, input);
        assert(audit.fontLoads.slice(0, beforeFonts).some(call => call.document === doc && call.loaded),
          'rendered event preceded actual STIX font readiness in child document');
        const urls = [...doc.querySelectorAll('link[href],style')].map(node => node.href || node.textContent).join('\n');
        assert(urls.includes('blob:') && !urls.includes('data:font/woff2;base64,'), 'mounted child did not use compact shared Blob CSS');
        const svg = await (await fetch(`/saved/${job}`)).text();
        equal(svg, captured.svg, 'saved standalone SVG bytes changed');
        assert(svg.includes(input.font_css), 'standalone SVG lost selfcontained font');
        const authored = new DOMParser().parseFromString(captured.html, 'text/html');
        const allocations = target => [...target.querySelectorAll('.linnet-math-label')].map(label => {
          const node = label.closest('foreignObject');
          return [node.getAttribute('width'), node.getAttribute('height')];
        });
        equal(allocations(doc), allocations(authored),
          'notebook reuse changed native label allocation');
        captures.push({job, receipt: record.receipt, geometry, svg_sha256: await digest(svg)});
        const measures = record.count('measure');
        record.dispose();
        record.dispose = null;
        assert([...record.callbacks.values()].every(callbacks => callbacks.size === 0), 'disposed view leaked listeners');
        record.events.length = 0;
        await probeDelta(0, () => record.mount(other), 'completed remount must not probe');
        await waitFor(() => record.events.length === 1, 'completed remount never rendered');
        equal(record.count('measure'), measures, 'completed remount remeasured');
        await inspect(record.events[0].document, input);
      }
      result = {captures, completion_order: captures.map(capture => capture.job)};
    } else if (name === 'font-lifetime') {
      // A native capture is supplied separately; never fabricate a final SVG.
      assert(request.completed_html, 'font-lifetime requires a retained native HTML capture');
      const active = [];
      const versions = [];
      const versionURLs = [];
      for (let index = 0; index < 5; index++) {
        const css = request.font_css + `\n/* resource-version-${index} */`;
        const input = {...variant(request, `mounted-version-${index}`, index ? request.labels : []), font_css: css, font_key: await digest(css)};
        const record = model(input);
        const block = `<style data-linnet-font="true"><![CDATA[${request.font_css}]]></style>`;
        record.fields.rendered = {id: input.id, html: request.completed_html.replaceAll(block, '')};
        const beforeURLs = audit.resources.length;
        await record.mount(module);
        await waitFor(() => record.events.length === 1, 'active font version never mounted');
        const resources = audit.resources.slice(beforeURLs);
        assert(resources.length > 0, 'new font version did not allocate its own resource');
        versionURLs.push(resources);
        active.push(record);
        versions.push(input);
      }
      assert(audit.live.length > 3, 'active font resources unexpectedly capped or absent');
      assert(audit.resources.every(entry => entry.revoked === 0), 'active font resource revoked');
      const activeURLs = audit.live.map(entry => entry.url);
      window.dispatchEvent(new PageTransitionEvent('pagehide', {persisted: true}));
      equal(audit.live.map(entry => entry.url), activeURLs, 'BFCache pagehide revoked reusable fonts');
      active[0].close(); active[1].close(); active[2].close(); active[3].close();
      assert(audit.resources.filter(entry => entry.revoked).length > 0, 'idle version overflow never revoked');
      const mountedDoc = active[4].events[0].document;
      const references = [...mountedDoc.querySelectorAll('link[href],style')].map(node => node.href || node.textContent).join('\n');
      for (const entry of audit.resources) if (references.includes(entry.url)) equal(entry.revoked, 0, 'still-mounted CSS URL revoked');
      active[4].close();
      for (let index = 0; index < versionURLs.length; index++) {
        assert(versionURLs[index].every(entry => entry.revoked === (index < 2 ? 1 : 0)),
          `three-version LRU retained/revoked the wrong font version ${index}`);
      }
      const stale = model(versions[0]);
      stale.fields.rendered = {id: stale.get('request').id,
        html: request.completed_html.replaceAll(`<style data-linnet-font="true"><![CDATA[${request.font_css}]]></style>`, '')};
      const creations = audit.resources.length;
      await stale.mount(other);
      await waitFor(() => stale.events.length === 1, 'evicted resource failed to remount');
      assert(audit.resources.length > creations, 'evicted font version did not recreate URLs');
      stale.close();
      window.dispatchEvent(new PageTransitionEvent('pagehide', {persisted: false}));
      equal(audit.live.length, 0, 'non-BFCache pagehide must release idle resource URLs');
      assert(audit.resources.every(entry => entry.revoked === 1), 'font URL leaked or revoked twice');
      result = {font_versions: versions.length, created_urls: audit.resources.length, revoked_urls: audit.resources.filter(entry => entry.revoked).length};
    } else throw new Error(`Unknown graph reuse case: ${name}`);
    assert(audit.resources.every(entry => audit.hostDocuments.has(entry.document)), 'font Blob URL created in a disposable child window');
    return {name, result, probe_frames: audit.probes.length, probe_rows: audit.probes.map(entry => entry.rows.length),
      resources: audit.resources.map(({url, type, revoked}) => ({url, type, revoked})),
      actual_font_loads: audit.fontLoads.filter(call => call.loaded).length};
  } finally {
    for (const record of models) record.close();
    window.dispatchEvent(new PageTransitionEvent('pagehide', {persisted: false}));
    for (const url of wrappers) URL.revokeObjectURL(url);
    audit.close();
  }
}
