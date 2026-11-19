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
  const windowErrors = [];
  const onError = event => windowErrors.push(event.message);
  window.addEventListener('error', onError);
  // A parent-realm observer cannot reliably deliver child-document targets.
  const ParentResizeObserver = window.ResizeObserver;
  window.ResizeObserver = class extends ParentResizeObserver {
    observe(target, options) {
      assert(target.ownerDocument === document, 'parent observer received a child-document target');
      return super.observe(target, options);
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
  // Match the native SVG serializer: without CDATA, markup in JS comments is
  // parsed as SVG and the script no longer matches the trusted controls source.
  const graph = label => `<figure style="margin:0"><svg data-linnet-interactive="true"
    viewBox="0 0 120 80" width="120pt" height="80pt">
    <script>window.untrustedOutputExecuted = true</script>
    <rect width="120" height="80" fill="#ffffff"/>
    <a data-linnet-kind="node" data-linnet-id="0" data-linnet-detail='{"edges":[]}'>
      <circle cx="60" cy="40" r="8" fill="#000000"/>
    </a><style>${css}</style><script><![CDATA[${controls}]]></script></svg>
    <figcaption>${label}</figcaption><button hidden id="untrusted" onclick="window.untrustedHandlerExecuted=true">Untrusted</button></figure>`;
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
    assert([...doc.querySelectorAll('script')].some(script => script.textContent === controls),
      'fixture changed the trusted controls source');
    assert(doc.querySelector('.linnet-viewport'), 'ready preceded graph initialization');
    assert(frame.getAttribute('sandbox') === 'allow-same-origin allow-scripts',
      'display frame sandbox changed');
    assert(doc.querySelector('meta[http-equiv="Content-Security-Policy"]').content.includes("script-src 'none'"),
      'display frame must block markup scripts and handlers');
    doc.querySelector('#untrusted').click();
    assert(!doc.defaultView.untrustedHandlerExecuted, 'untrusted inline handler executed');
    assert(!doc.defaultView.untrustedOutputExecuted, 'untrusted output script executed');
    assert([...doc.querySelectorAll('script')].every(script => script.type === 'application/x-linnet-inert'),
      'output scripts were not inert before frame load');
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
    const graphBackground = viewport.querySelector('rect');
    const scroll = window.scrollY;
    graphBackground.dispatchEvent(new doc.defaultView.PointerEvent('pointerdown',
      {button:0, pointerId:1, bubbles:true}));
    assert(doc.activeElement === viewport, 'graph pointer did not focus keyboard viewport');
    assert(window.scrollY === scroll, 'graph focus scrolled notebook');
    graphBackground.dispatchEvent(new doc.defaultView.PointerEvent('pointerup',
      {button:0, pointerId:1, bubbles:true}));
    doc.activeElement.dispatchEvent(new doc.defaultView.KeyboardEvent('keydown', {key:'+', bubbles:true}));
    svg.querySelector('[data-linnet-kind="node"]').dispatchEvent(
      new doc.defaultView.MouseEvent('click', {bubbles:true, shiftKey:true}));
    assert(svg.linnetSelection.nodes.includes(0), 'selection is inert');
    const camera = viewport.getAttribute('viewBox');
    assert(doc.querySelector('.linnet-control-row output').textContent === '105%', 'zoom is inert');
    assert(!svg.querySelector('.linnet-toolbar'), 'navigation still overlaps the graph');
    assert(doc.querySelector('.linnet-control-row').getBoundingClientRect().bottom <= svg.getBoundingClientRect().top + 1,
      'navigation row overlaps graph labels');
    if (mode === 'delayed') {
      let completedFontFits;
      doc.fonts.addEventListener('loadingdone', () => {
        // Check after all loadingdone listeners, before a deferred resize fit.
        queueMicrotask(() => {
          completedFontFits = parseFloat(frame.style.height) >= doc.body.getBoundingClientRect().height;
        });
      }, {once:true});
      publish('font', {key:state.request.font_key, css:font});
      await waitFor(() => doc.querySelector('style[data-linnet-font]'), 'late font CSS never arrived');
      // Font decoding is asynchronous and need not finish within two frames.
      // Exercise its actual completion boundary, including the immediate fit.
      await waitFor(() => completedFontFits !== undefined, 'late font never completed');
      assert(completedFontFits, 'font completion deferred frame sizing until another resize turn');
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
    for (const width of [240, 500, 320, 500]) {
      host.style.width = `${width}px`;
      await tick(); await tick(); await tick();
      assert(parseFloat(frame.style.height) >= doc.body.getBoundingClientRect().height,
        'responsive frame height became stale');
    }
    cleanup(); cleanup = null;
    await waitFor(() => [...callbacks.values()].every(value => value.size === 0),
      'disposed view retained font or render listeners');
    assert(!host.querySelector('iframe'), 'disposed view retained its iframe');
    assert(!windowErrors.length, `window errors: ${windowErrors.join('; ')}`);
    assert(sent.filter(message => message.action === 'font').length === 1,
      'replacement duplicated the font lease');
    return {mode, theme, first_ms:firstMs, displays:events.length,
      font_ready_reads:frames.reduce((sum, entry) => sum + entry.readyReads, 0)};
  } finally {
    cleanup?.();
    document.createElement = create;
    window.ResizeObserver = ParentResizeObserver;
    window.removeEventListener('error', onError);
    host.remove();
    URL.revokeObjectURL(url);
  }
} """


def check_pointer_controls(browser, sources):
    page = browser.new_page(has_touch=True, viewport={"width": 820, "height": 600})
    try:
        page.set_content("""<div style="height:700px"></div>
          <svg data-linnet-interactive="true" viewBox="0 0 120 80" width="480" height="320">
          <rect width="120" height="80" fill="white" pointer-events="none"/>
          <a data-linnet-kind="node" data-linnet-id="0" data-linnet-detail='{"edges":[]}'>
            <circle cx="60" cy="40" r="8" fill="black"/>
          </a><foreignObject x="5" y="5" width="35" height="15">
          <input xmlns="http://www.w3.org/1999/xhtml" style="width:25px" aria-label="Graph input"/>
          </foreignObject><foreignObject x="5" y="25" width="35" height="15">
          <button xmlns="http://www.w3.org/1999/xhtml">Action</button></foreignObject></svg>""")
        page.add_style_tag(content=sources["css"])
        page.evaluate(sources["controls"])
        viewport = page.locator(".linnet-viewport")
        viewport.scroll_into_view_if_needed()
        bounds = viewport.bounding_box()
        point = {
            "x": bounds["x"] + bounds["width"] * 0.8,
            "y": bounds["y"] + bounds["height"] * 0.8,
        }
        scroll = page.evaluate("scrollY")
        page.mouse.click(**point)
        assert viewport.evaluate("x => x === document.activeElement"), (
            "mouse focus lost"
        )
        assert page.evaluate("scrollY") == scroll, "mouse focus scrolled page"
        page.keyboard.press("=")
        assert page.locator("output").inner_text() == "105%"
        page.keyboard.press("0")
        camera = viewport.get_attribute("viewBox")
        page.locator("svg[data-linnet-interactive]").evaluate(
            "x=>x.setAttribute('tabindex','0')"
        )
        page.mouse.click(**point)
        assert viewport.evaluate("x => x === document.activeElement"), (
            "outer SVG retained focus"
        )
        page.keyboard.press("=")
        assert page.locator("output").inner_text() == "105%"
        page.keyboard.press("0")
        page.mouse.move(**point)
        page.mouse.down()
        page.mouse.move(point["x"] + 30, point["y"] + 20, steps=5)
        page.mouse.up()
        assert viewport.get_attribute("viewBox") != camera, "drag did not pan"
        page.keyboard.press("0")
        assert viewport.get_attribute("viewBox") == camera, "camera reset failed"
        page.locator("circle").hover()
        assert page.get_by_role("tooltip").is_visible(), "graph tooltip missing"
        page.locator("circle").click(modifiers=["Shift"])
        assert page.locator("svg[data-linnet-interactive]").evaluate(
            "x => x.linnetSelection.nodes.includes(0)"
        ), "pointer selection failed"
        field = page.get_by_role("textbox", name="Graph input")
        field.click()
        assert field.evaluate("x => x === document.activeElement"), "input focus stolen"
        button = page.get_by_role("button", name="Action")
        button.click()
        assert button.evaluate("x => x === document.activeElement"), (
            "button focus stolen"
        )
        viewport.scroll_into_view_if_needed()
        bounds = viewport.bounding_box()
        page.touchscreen.tap(
            bounds["x"] + bounds["width"] * 0.8, bounds["y"] + bounds["height"] * 0.8
        )
        assert viewport.get_attribute("viewBox") == camera, "touch changed camera"
        print(
            json.dumps(
                {
                    "pointer_controls": "mouse, keyboard, drag, reset, selection, input, button, touch",
                    "passed": True,
                }
            ),
            flush=True,
        )
    finally:
        page.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    engines = parser.add_mutually_exclusive_group(required=True)
    engines.add_argument("--chromium")
    engines.add_argument("--webkit")
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
        browser = (
            driver.webkit.launch(executable_path=args.webkit)
            if args.webkit
            else driver.chromium.launch(executable_path=args.chromium)
        )
        try:
            check_pointer_controls(browser, sources)
            for mode in ("pending", "rejected", "delayed"):
                for theme in ("light", "dark"):
                    page = browser.new_page(
                        color_scheme="dark" if theme == "light" else "light"
                    )
                    errors = []
                    csp_errors = []
                    page.on(
                        "console",
                        lambda message, csp_errors=csp_errors: (
                            csp_errors.append(message.text)
                            if message.type == "error"
                            and "content security policy" in message.text.lower()
                            else None
                        ),
                    )
                    page.on(
                        "pageerror",
                        lambda error, errors=errors: errors.append(str(error)),
                    )
                    try:
                        result = page.evaluate(
                            CHECK, {**sources, "mode": mode, "theme": theme}
                        )
                    except Exception as error:
                        error.add_note(
                            f"mode={mode}, theme={theme}, page errors={errors}"
                        )
                        raise
                    assert not errors, errors
                    # Each mounted snapshot deliberately clicks one hostile handler.
                    # Require its refusal, preserving failures for other CSP errors.
                    assert len(csp_errors) == 2 and all(
                        "inline event handler" in message.lower()
                        for message in csp_errors
                    ), csp_errors
                    print(json.dumps(result), flush=True)
                    page.close()
        finally:
            browser.close()


if __name__ == "__main__":
    main()
