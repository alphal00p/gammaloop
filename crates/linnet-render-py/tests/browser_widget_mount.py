"""Check widget snapshot visibility independently of font readiness.

Run with Playwright installed and --chromium pointing to a browser executable.
These fixtures exercise mounting and controls, not native label measurements.
"""

import argparse
import base64
import json
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[3]
CHECK = r"""async ({javascript, controls, css, font, mode, theme}) => {
  const assert = (value, message) => { if (!value) throw new Error(message); };
  const tick = () => new Promise(resolve => requestAnimationFrame(resolve));
  const waitFor = async (predicate, message) => {
    const deadline = performance.now() + 2500;
    while (!predicate()) {
      assert(performance.now() < deadline, message);
      await new Promise(resolve => setTimeout(resolve, 5));
    }
  };
  document.documentElement.dataset.theme = theme;
  document.documentElement.style.colorScheme = theme;
  const host = document.createElement('section');
  host.style.width = '500px';
  document.body.append(host);
  const frames = [], events = [], sent = [], callbacks = new Map();
  const create = document.createElement;
  document.createElement = function (...args) {
    const element = create.apply(this, args);
    if (args[0] === 'iframe') {
      const entry = {frame:element, readyReads:0, loadCalls:0};
      frames.push(entry);
      // Install before the production load handler. Display must neither read
      // fonts.ready nor invoke exact-font loading, even when either rejects.
      element.addEventListener('load', () => {
        const doc = element.contentDocument;
        entry.doc = doc;
        Object.defineProperty(doc.fonts, 'ready', {configurable:true, get() {
          entry.readyReads++;
          return mode === 'rejected' ? Promise.reject(new Error('font unavailable'))
            : new Promise(() => {});
        }});
        doc.fonts.load = () => {
          entry.loadCalls++;
          return Promise.reject(new Error('display must not measure labels'));
        };
      });
    }
    return element;
  };
  const graph = label => `<figure style="margin:0"><svg data-linnet-interactive="true"
    viewBox="0 0 120 80" width="120pt" height="80pt">
    <rect width="120" height="80" fill="#ffffff"/>
    <a data-linnet-kind="node" data-linnet-id="0" data-linnet-detail='{"edges":[]}'>
      <circle cx="60" cy="40" r="8" fill="#000000"/>
    </a><style>${css}</style><script>${controls}</script></svg>
    <figcaption>${label}</figcaption></figure>`;
  const state = {
    request:{id:'snapshot', labels:[], pending:false, font_key:`font-${mode}`},
    rendered:{id:'snapshot', html:graph('first diagram')}, scripts:[controls], font:{},
  };
  const publish = (field, value) => {
    state[field] = value;
    for (const callback of callbacks.get(`change:${field}`) || []) callback();
  };
  const model = {
    get: name => state[name],
    on(name, callback) {
      if (!callbacks.has(name)) callbacks.set(name, new Set());
      callbacks.get(name).add(callback);
    },
    off: (name, callback) => callbacks.get(name)?.delete(callback),
    send(message) {
      sent.push(message);
      if (message.action === 'font' && mode === 'rejected')
        publish('font', {key:state.request.font_key, css:'invalid font resource'});
    },
  };
  host.addEventListener('linnet:rendered', event => {
    const doc = event.detail.document;
    const frame = host.querySelector('iframe');
    assert(frame.contentDocument === doc, 'stale frame emitted ready');
    assert(getComputedStyle(frame).visibility === 'visible', 'ready frame remains hidden');
    assert(doc.querySelector('.linnet-viewport'), 'ready preceded graph initialization');
    assert(parseFloat(frame.style.height) >= doc.body.getBoundingClientRect().height,
      'ready preceded frame sizing');
    events.push(event.detail);
  });
  const url = URL.createObjectURL(new Blob([javascript], {type:'text/javascript'}));
  let cleanup;
  try {
    const module = await import(url);
    const started = performance.now();
    const rendering = module.default.render({model, el:host});
    assert(host.querySelector('iframe').style.visibility === 'hidden',
      'initial iframe was inserted visibly');
    cleanup = await rendering;
    await waitFor(() => events.length === 1, 'font state prevented initial snapshot');
    const firstMs = performance.now() - started;
    const firstDoc = events[0].document;
    assert(firstDoc.querySelector('svg').dataset.theme === theme, 'initial notebook theme lost');
    assert(firstDoc.documentElement.style.colorScheme === theme, 'initial canvas scheme mismatch');
    await waitFor(() => sent.some(message => message.action === 'font'), 'font not requested in background');
    const nextTheme = theme === 'dark' ? 'light' : 'dark';
    document.documentElement.dataset.theme = nextTheme;
    document.documentElement.style.colorScheme = nextTheme;
    publish('rendered', {id:'snapshot', html:graph('replacement diagram')});
    assert(host.querySelector('iframe').style.visibility === 'hidden',
      'replacement iframe was inserted visibly');
    await waitFor(() => events.length === 2, 'font state prevented replacement snapshot');
    const doc = events[1].document, frame = host.querySelector('iframe');
    const svg = doc.querySelector('svg'), viewport = svg.querySelector('.linnet-viewport');
    assert(doc !== firstDoc && !frames[0].frame.isConnected, 'old frame remained mounted');
    assert(svg.dataset.theme === nextTheme, 'replacement notebook theme lost');
    assert(doc.documentElement.style.colorScheme === nextTheme, 'replacement canvas scheme mismatch');
    viewport.dispatchEvent(new doc.defaultView.KeyboardEvent('keydown', {key:'+', bubbles:true}));
    svg.querySelector('[data-linnet-kind="node"]').dispatchEvent(
      new doc.defaultView.MouseEvent('click', {bubbles:true, shiftKey:true}));
    assert(svg.linnetSelection.nodes.includes(0), 'selection is inert');
    const camera = viewport.getAttribute('viewBox');
    assert(svg.querySelector('output').textContent === '105%', 'zoom is inert');
    if (mode === 'delayed') {
      publish('font', {key:state.request.font_key, css:font});
      await waitFor(() => doc.querySelector('style[data-linnet-font]'), 'late font CSS never arrived');
      await tick(); await tick();
      assert([...doc.fonts].some(face => face.family.replace(/["']/g, '') === 'Linnet Offline Math'
        && face.display === 'swap'), 'late font can suppress fallback label paint');
      assert(!firstDoc.querySelector('style[data-linnet-font]'), 'late font reached a stale document');
      assert(frame === host.querySelector('iframe') && frame.contentDocument === doc,
        'late font replaced the initialized frame');
      assert(svg.linnetSelection.nodes.includes(0), 'late font cleared selection');
      assert(viewport.getAttribute('viewBox') === camera, 'late font changed zoom');
      assert(parseFloat(frame.style.height) >= doc.body.getBoundingClientRect().height,
        'late font left the frame height stale');
    }
    assert(frames.every(entry => entry.readyReads === 0 && entry.loadCalls === 0),
      'snapshot display waited for exact fonts');
    assert(getComputedStyle(frame).visibility === 'visible', 'snapshot became hidden');
    cleanup(); cleanup = null;
    await waitFor(() => [...callbacks.values()].every(value => value.size === 0),
      'disposed view retained font or render listeners');
    assert(!host.querySelector('iframe'), 'disposed view retained its iframe');
    assert(sent.filter(message => message.action === 'font').length === 1,
      'replacement duplicated the font lease');
    return {mode, theme, first_ms:firstMs, displays:events.length,
      font_ready_reads:frames.reduce((sum, entry) => sum + entry.readyReads, 0)};
  } finally {
    cleanup?.();
    document.createElement = create;
    host.remove();
    URL.revokeObjectURL(url);
  }
} """


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromium")
    args = parser.parse_args()
    font = base64.b64encode(
        (ROOT / "crates/linnet-render-py/fonts/STIXTwoMath-Regular.woff2").read_bytes()
    ).decode("ascii")
    sources = {
        "javascript": (ROOT / "crates/linnet-render-py/src/view.js").read_text(),
        "controls": (ROOT / "crates/linnest/src/svg/interactive.js").read_text(),
        "css": (ROOT / "crates/linnest/src/svg/interactive.css").read_text(),
        "font": '@font-face{font-family:"Linnet Offline Math";src:url("data:font/woff2;base64,'
        + font
        + '");font-display:block}figcaption{font-family:"Linnet Offline Math"}',
    }
    with sync_playwright() as driver:
        browser = driver.chromium.launch(executable_path=args.chromium)
        try:
            for mode in ("pending", "rejected", "delayed"):
                for theme in ("light", "dark"):
                    page = browser.new_page(
                        color_scheme="dark" if theme == "light" else "light"
                    )
                    errors = []
                    page.on(
                        "pageerror",
                        lambda error, errors=errors: errors.append(str(error)),
                    )
                    result = page.evaluate(
                        CHECK, {**sources, "mode": mode, "theme": theme}
                    )
                    assert not errors, errors
                    print(json.dumps(result), flush=True)
                    page.close()
        finally:
            browser.close()


if __name__ == "__main__":
    main()
