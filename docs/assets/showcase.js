/* GammaLoop showcase animation.
   Everything on the 1920x1080 stage is a pure function of time, so the same
   program drives the interactive preview and the frame-by-frame video render. */
(() => {
  'use strict';

  const W = 1920;
  const H = 1080;
  const clamp = (x, a = 0, b = 1) => Math.min(b, Math.max(a, x));
  const lerp = (a, b, p) => a + (b - a) * p;
  const ease = {
    linear: (p) => p,
    inOut: (p) => (p < 0.5 ? 2 * p * p : 1 - Math.pow(-2 * p + 2, 2) / 2),
    out: (p) => 1 - Math.pow(1 - p, 3),
    in: (p) => p * p * p,
    outBack: (p) => {
      const c1 = 1.2;
      const c3 = c1 + 1;
      return 1 + c3 * Math.pow(p - 1, 3) + c1 * Math.pow(p - 1, 2);
    },
  };

  function rng(seed) {
    let s = seed >>> 0;
    return () => {
      s = (s + 0x6d2b79f5) | 0;
      let t = Math.imul(s ^ (s >>> 15), 1 | s);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gaussian(next) {
    const u = Math.max(next(), 1e-12);
    const v = next();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

  function el(tag, attrs = {}, ...children) {
    const node = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (k === 'class') node.className = v;
      else if (k === 'style') node.style.cssText = v;
      else if (k === 'html') node.innerHTML = v;
      else node.setAttribute(k, v);
    }
    for (const c of children) node.append(c);
    return node;
  }
  const SVG = 'http://www.w3.org/2000/svg';
  function svg(tag, attrs = {}, ...children) {
    const node = document.createElementNS(SVG, tag);
    for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
    for (const c of children) node.append(c);
    return node;
  }

  // The GammaLoop mark from assets/gammalooplogo-dark.svg, recoloured through CSS classes.
  const LOGO_SVG = "<svg viewBox=\"0 0 628.582677165 311.811023622\" width=\"628.582677165pt\" height=\"311.811023622pt\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:h5=\"http://www.w3.org/1999/xhtml\"><g transform=\"translate(28.346456693 28.346456693)\"><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(148.818897638 63.779527559)\" d=\"M 0 0h 382.677165354\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"square\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(540)\" d=\"M 0 0m 0 63.779527559c 17.612230211 0 31.88976378 -14.277533568 31.88976378 -31.88976378c 0 -17.612230211 -14.277533568 -31.88976378 -31.88976378 -31.88976378\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"square\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(488.976377953)\" d=\"M 0 0m 42.519685039 0h -42.519685039\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(488.976377953)\" d=\"M 0 0v 42.519685039\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(488.976377953 85.039370079)\" d=\"M 0 0v 55.275590551\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(0 63.779527559)\" d=\"M 0 0m 123.307086614 148.818897638c 0 -29.763779528 -12.755905512 -110.551181102 -123.307086614 -148.818897638\"/><path fill=\"none\" class=\"mask\" stroke-width=\"48\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(38.267716535 63.779527559)\" d=\"M 0 0m 110.551181102 0c -36.566929134 0 -110.551181102 72.283464567 -110.551181102 148.818897638\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"round\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(38.267716535 63.779527559)\" d=\"M 0 0m 110.551181102 0c -36.566929134 0 -110.551181102 72.283464567 -110.551181102 148.818897638\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(38.267716535 212.598425197)\" d=\"M 0 0c 0 23.482973615 19.036711425 42.519685039 42.519685039 42.519685039c 23.482973615 0 42.519685039 -19.036711425 42.519685039 -42.519685039\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"square\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(148.818897638)\" d=\"M 0 0m 0 29.763779528v -29.763779528\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"square\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(148.818897638 97.795275591)\" d=\"M 0 0v 29.763779528h 59.527559055\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(233.858267717)\" d=\"M 0 0m 85.039370079 42.519685039c 0 -23.482973615 -19.036711425 -42.519685039 -42.519685039 -42.519685039c -23.482973615 0 -42.519685039 19.036711425 -42.519685039 42.519685039\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(233.858267717 85.039370079)\" d=\"M 0 0m 0 0c 0 23.482973615 19.036711425 42.519685039 42.519685039 42.519685039c 23.482973615 0 42.519685039 -19.036711425 42.519685039 -42.519685039\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(361.417322835)\" d=\"M 0 0m 85.039370079 42.519685039c 0 -23.482973615 -19.036711425 -42.519685039 -42.519685039 -42.519685039c -23.482973615 0 -42.519685039 19.036711425 -42.519685039 42.519685039\"/><path fill=\"none\" class=\"ink\" stroke-width=\"25\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\" stroke-miterlimit=\"4\" transform=\"translate(361.417322835 85.039370079)\" d=\"M 0 0m 0 0c 0 23.482973615 19.036711425 42.519685039 42.519685039 42.519685039c 23.482973615 0 42.519685039 -19.036711425 42.519685039 -42.519685039\"/></g></svg>";

  // ---------------------------------------------------------------- timeline
  class Timeline {
    constructor() {
      this.segs = [];
      this.end = 0;
    }
    add(start, dur, fn, easing = ease.inOut) {
      this.segs.push({ start, dur, fn, easing });
      this.end = Math.max(this.end, start + dur);
      return this;
    }
    seek(t) {
      for (const s of this.segs) {
        const raw = s.dur > 0 ? clamp((t - s.start) / s.dur) : t >= s.start ? 1 : 0;
        s.fn(s.easing(raw), t - s.start, t);
      }
    }
  }
  // The whole stage is a pure function of time: build() wires every scene to one timeline.
  function build(stage) {
    const tl = new Timeline();

    const scenes = [];
    function scene(id, name, start, end) {
      const node = el('section', { class: `scene ${id}` });
      stage.appendChild(node);
      const FADE = 0.55;
      tl.add(
        start - 0.01,
        end - start + 0.02,
        (_, __, t) => {
          const on = t >= start && t < end;
          node.classList.toggle('on', on);
          node.style.opacity = on ? clamp(Math.min((t - start) / FADE, (end - t) / FADE)) : 0;
        },
        ease.linear,
      );
      scenes.push({ id, name, start, end });
      return node;
    }

    function fadeIn(node, at, dur = 0.7, y = 26) {
      node.style.opacity = 0;
      tl.add(
        at,
        dur,
        (p) => {
          node.style.opacity = p;
          // SVG nodes keep their transform attribute: a CSS transform would override it.
          if (y) node.style.transform = `translateY(${(1 - p) * y}px)`;
        },
        ease.out,
      );
      return node;
    }
    function fadeOut(node, at, dur = 0.5) {
      tl.add(
        at,
        dur,
        (p) => {
          if (p > 0) node.style.opacity = 1 - p;
        },
        ease.inOut,
      );
    }
    function drawPath(path, at, dur, easing = ease.inOut) {
      let len = null;
      path.style.opacity = 0;
      tl.add(
        at,
        dur,
        (p) => {
          if (len === null) len = path.getTotalLength();
          // Drop the dash pattern once fully drawn so no seam remains at the dash boundary.
          path.style.strokeDasharray = p >= 1 ? 'none' : `${len} ${len}`;
          path.style.strokeDashoffset = p >= 1 ? '0' : String(len * (1 - p));
          path.style.opacity = p > 0 ? 1 : 0;
        },
        easing,
      );
    }
    function typeText(node, text, at, dur) {
      node.textContent = '';
      tl.add(at, dur, (p) => (node.textContent = text.slice(0, Math.round(text.length * p))), ease.linear);
    }

    // ------------------------------------------------------------ terminal
    // A terminal is a list of timed chunks re-rendered from scratch on every seek.
    const seg = (cls, text) => ({ cls, text });
    const C = (t) => seg('cmd', t);
    const F = (t) => seg('flag', t);
    const S = (t) => seg('str', t);
    const O = (t) => seg('out', t);
    const HI = (t) => seg('hi', t);
    const OK = (t) => seg('ok', t);
    const B = (t) => seg('blue', t);

    function segsHtml(segs, budget) {
      let html = '';
      let left = budget;
      for (const s of segs) {
        if (left <= 0) break;
        const text = s.text.slice(0, left);
        left -= text.length;
        html += s.cls ? `<span class="${s.cls}">${esc(text)}</span>` : esc(text);
      }
      return html;
    }
    const segsLen = (segs) => segs.reduce((n, s) => n + s.text.length, 0);
    const toSegs = (line) => (typeof line === 'string' ? [O(line)] : line);

    class Term {
      constructor(parent, { title = 'gammaloop', extraClass = '' } = {}) {
        this.root = el('div', { class: `term ${extraClass}` });
        this.root.append(
          el(
            'div',
            { class: 'bar' },
            el('span', { class: 'dot' }),
            el('span', { class: 'dot' }),
            el('span', { class: 'dot' }),
            el('span', { class: 'title' }, title),
          ),
        );
        this.body = el('div', { class: 'body' });
        this.root.append(this.body);
        parent.append(this.root);
        this.chunks = [];
        this.mounted = null;
        tl.add(0, 1e9, (_, __, t) => this.render(t), ease.linear);
      }
      type(at, dur, segs, { prompt = true } = {}) {
        const all = prompt ? [seg('prompt', '$ '), ...segs] : segs;
        this.chunks.push({ kind: 'type', start: at, end: at + dur, segs: all, prompt });
        return at + dur;
      }
      lines(at, per, lines) {
        const segLines = lines.map(toSegs);
        this.chunks.push({ kind: 'lines', start: at, end: at + per * segLines.length, lines: segLines });
        return at + per * segLines.length;
      }
      prompt(at) {
        this.chunks.push({ kind: 'prompt', start: at, end: at });
        return at;
      }
      clear(at) {
        this.chunks.push({ kind: 'clear', start: at, end: at });
        return at;
      }
      mount(at, node, until) {
        this.chunks.push({ kind: 'mount', start: at, end: until, node });
        return until;
      }
      render(t) {
        let from = 0;
        this.chunks.forEach((c, i) => {
          if (c.kind === 'clear' && c.start <= t) from = i + 1;
        });
        const visible = this.chunks.slice(from).filter((c) => c.start <= t);
        const mounted = visible.find((c) => c.kind === 'mount' && t < c.end);
        if (mounted) {
          if (this.mounted !== mounted.node) {
            this.body.replaceChildren(mounted.node);
            this.mounted = mounted.node;
            this.body.classList.add('tui-host');
          }
          return;
        }
        if (this.mounted) {
          this.mounted = null;
          this.body.classList.remove('tui-host');
        }
        let html = '';
        const blink = Math.floor(t * 2.4) % 2 === 0;
        visible.forEach((c, i) => {
          const next = this.chunks[from + visible.indexOf(c) + 1];
          const idle = !next || next.start > t;
          if (c.kind === 'type') {
            const total = segsLen(c.segs);
            const p = c.end > c.start ? clamp((t - c.start) / (c.end - c.start)) : 1;
            const n = Math.round(total * p);
            html += segsHtml(c.segs, n);
            if (t < c.end) html += '<span class="caret"></span>';
            else if (idle && blink) html += '<span class="caret"></span>';
            html += '\n';
          } else if (c.kind === 'lines') {
            const per = (c.end - c.start) / c.lines.length;
            const k = Math.min(c.lines.length, Math.floor((t - c.start) / per) + 1);
            for (let j = 0; j < k; j++) {
              const line = c.lines[j];
              const hi = line.some((s) => s.cls === 'line-hi');
              const body = segsHtml(
                line.filter((s) => s.cls !== 'line-hi'),
                1e9,
              );
              html += hi ? `<span class="line-hi">${body}</span>\n` : `${body}\n`;
            }
          } else if (c.kind === 'prompt') {
            html += `<span class="prompt">$ </span>${blink ? '<span class="caret"></span>' : ''}`;
          }
          void i;
        });
        this.body.innerHTML = html;
        this.body.scrollTop = this.body.scrollHeight;
      }
    }

    // ------------------------------------------------------ syntax highlight
    const KW = {
      python: ['from', 'import', 'with', 'as', 'def', 'return', 'print', 'for', 'in', 'if', 'else', 'None', 'True', 'False', 'class'],
      rust: ['use', 'let', 'mut', 'for', 'in', 'fn', 'pub', 'struct', 'impl', 'match', 'Some', 'None', 'Ok', 'Err', 'return', 'as', 'self', 'crate'],
    };
    function tokenize(code, lang) {
      const out = [];
      const re = new RegExp(
        [
          lang === 'python' ? '(#[^\\n]*)' : '(//[^\\n]*)',
          '("(?:[^"\\\\]|\\\\.)*")',
          '(\\b\\d+(?:\\.\\d+)?(?:e-?\\d+)?\\b)',
          `(\\b(?:${KW[lang].join('|')})\\b)`,
          '([A-Za-z_][A-Za-z0-9_]*)(?=\\s*[(!])',
          '(\\b[A-Z][A-Za-z0-9_]*\\b)',
          '([A-Za-z_][A-Za-z0-9_]*)',
        ].join('|'),
        'g',
      );
      let last = 0;
      let m;
      while ((m = re.exec(code))) {
        if (m.index > last) out.push(seg('', code.slice(last, m.index)));
        const [text] = m;
        const cls = m[1] ? 'cm' : m[2] ? 'st' : m[3] ? 'nu' : m[4] ? 'kw' : m[5] ? 'fn' : m[6] ? 'ty' : '';
        out.push(seg(cls, text));
        last = m.index + text.length;
      }
      if (last < code.length) out.push(seg('', code.slice(last)));
      return out;
    }
    function typeCode(node, code, lang, at, dur) {
      const segs = tokenize(code, lang);
      const total = segsLen(segs);
      node.innerHTML = '';
      tl.add(
        at,
        dur,
        (p, _, t) => {
          const n = Math.round(total * p);
          node.innerHTML = segsHtml(segs, n) + (t < at + dur ? '<span class="caret"></span>' : '');
          node.scrollTop = node.scrollHeight;
        },
        ease.linear,
      );
    }

    // ---------------------------------------------------- Feynman primitives
    function wavyPath(x1, y1, x2, y2, amp = 9, waves = 6) {
      const dx = x2 - x1;
      const dy = y2 - y1;
      const L = Math.hypot(dx, dy);
      const nx = -dy / L;
      const ny = dx / L;
      const N = waves * 2;
      let d = `M${x1} ${y1}`;
      for (let i = 0; i < N; i++) {
        const t0 = i / N;
        const t1 = (i + 1) / N;
        const tm = (t0 + t1) / 2;
        const s = i % 2 ? -1 : 1;
        d += ` Q${x1 + dx * tm + nx * amp * 2 * s} ${y1 + dy * tm + ny * amp * 2 * s} ${x1 + dx * t1} ${y1 + dy * t1}`;
      }
      return d;
    }
    function curlyPath(x1, y1, x2, y2, r = 10, loops = 7) {
      const dx = x2 - x1;
      const dy = y2 - y1;
      const L = Math.hypot(dx, dy);
      const ux = dx / L;
      const uy = dy / L;
      const nx = -uy;
      const ny = ux;
      const tail = 14;
      const inner = L - 2 * tail;
      const total = loops * 2 * Math.PI;
      const a = inner / total;
      const pts = [[x1, y1]];
      const steps = loops * 28;
      for (let i = 0; i <= steps; i++) {
        const th = (total * i) / steps;
        const along = tail + a * th - r * Math.sin(th);
        const perp = r * (1 - Math.cos(th)) - r;
        pts.push([x1 + ux * along + nx * perp, y1 + uy * along + ny * perp]);
      }
      pts.push([x2, y2]);
      return 'M' + pts.map(([x, y]) => `${x.toFixed(1)} ${y.toFixed(1)}`).join(' L');
    }
    function smoothPath(points) {
      let d = `M${points[0][0]} ${points[0][1]}`;
      for (let i = 0; i < points.length - 1; i++) {
        const p0 = points[Math.max(i - 1, 0)];
        const p1 = points[i];
        const p2 = points[i + 1];
        const p3 = points[Math.min(i + 2, points.length - 1)];
        const c1 = [p1[0] + (p2[0] - p0[0]) / 6, p1[1] + (p2[1] - p0[1]) / 6];
        const c2 = [p2[0] - (p3[0] - p1[0]) / 6, p2[1] - (p3[1] - p1[1]) / 6];
        d += ` C${c1[0]} ${c1[1]} ${c2[0]} ${c2[1]} ${p2[0]} ${p2[1]}`;
      }
      return d;
    }
    function arrowHead(x, y, angle, size = 13, cls = 'vertex') {
      const a = svg('path', {
        class: cls,
        d: `M${size} 0 L${-size * 0.8} ${size * 0.7} L${-size * 0.45} 0 L${-size * 0.8} ${-size * 0.7} Z`,
        transform: `translate(${x} ${y}) rotate(${(angle * 180) / Math.PI})`,
      });
      return a;
    }
    function fermion(x1, y1, x2, y2, cls = 'edge') {
      const g = svg('g');
      const line = svg('path', { class: cls, d: `M${x1} ${y1} L${x2} ${y2}` });
      const head = arrowHead((x1 + x2) / 2, (y1 + y2) / 2, Math.atan2(y2 - y1, x2 - x1));
      g.append(line, head);
      return { g, line, head };
    }

    // ---------------------------------------------------------------- schedule
    const T = {
      open: [0, 8],
      lu: [8, 24.5],
      cli: [24.5, 36.5],
      gen: [36.5, 55],
      int: [55, 73.5],
      api: [73.5, 85.5],
      eco: [85.5, 95.5],
      out: [95.5, 101.5],
    };

    // ------------------------------------------------------ background field
    const field = stage.querySelector('canvas');
    const fctx = field.getContext('2d');
    field.width = W;
    field.height = H;
    const fieldPts = (() => {
      const next = rng(20260926);
      return Array.from({ length: 260 }, () => ({
        x: next() * W,
        y: next() * H,
        r: 1 + next() * 2.2,
        ph: next() * Math.PI * 2,
        sp: 0.15 + next() * 0.35,
        a: 0.05 + next() * 0.16,
      }));
    })();
    tl.add(
      0,
      1e9,
      (_, __, t) => {
        fctx.clearRect(0, 0, W, H);
        for (const p of fieldPts) {
          const x = p.x + Math.sin(t * p.sp + p.ph) * 18;
          const y = p.y + Math.cos(t * p.sp * 0.8 + p.ph) * 12;
          const tw = 0.7 + 0.3 * Math.sin(t * 1.3 + p.ph * 3);
          fctx.globalAlpha = p.a * tw;
          fctx.fillStyle = '#b893c7';
          fctx.beginPath();
          fctx.arc(x, y, p.r, 0, Math.PI * 2);
          fctx.fill();
        }
        fctx.globalAlpha = 1;
      },
      ease.linear,
    );

    function logoNode(className = 'logo') {
      const wrapper = document.createElement('div');
      wrapper.innerHTML = LOGO_SVG;
      const node = wrapper.firstElementChild;
      node.classList.add(className);
      return node;
    }

    // ================================================================ SCENES
    // -- Cold open ----------------------------------------------------------
    {
      const [a, b] = T.open;
      const sc = scene('s-open', 'Cold open', a, b);
      const logo = logoNode();
      sc.append(logo);
      const paths = [...logo.querySelectorAll('path.ink')];
      const mask = logo.querySelector('path.mask');
      mask.style.opacity = 0;
      tl.add(a + 1.2, 0.01, (p) => (mask.style.opacity = p), ease.linear);
      paths.forEach((path, i) => drawPath(path, a + 0.25 + i * 0.14, 1.1, ease.inOut));

      const word = el('div', { class: 'wordmark' }, 'GammaLoop');
      const tag = el('div', { class: 'tagline' }, 'Local cancellation. Global precision.');
      const sub = el(
        'div',
        { class: 'sub' },
        'Differential collider cross-sections with Local Unitarity · CLI · Rust · Python',
      );
      sc.append(word, tag, sub);
      fadeIn(word, a + 2.6, 0.9, 30);
      fadeIn(tag, a + 3.7, 0.9, 20);
      fadeIn(sub, a + 4.7, 0.8, 14);
    }

    // -- Local Unitarity ----------------------------------------------------
    {
      const [a, b] = T.lu;
      const sc = scene('s-lu', 'Local Unitarity', a, b);
      const copy = el('div', { class: 'copy' });
      const kicker = el('div', { class: 'kicker' }, 'The method');
      const display = el('div', { class: 'display', html: 'Real and virtual, <em>cancelled locally.</em>' });
      const lede = el(
        'div',
        { class: 'lede' },
        'Local Unitarity rewrites one forward-scattering graph as the sum of its cuts on shared loop momenta. Infrared singularities cancel point by point, before any Monte Carlo integral is taken.',
      );
      copy.append(kicker, display, lede);
      sc.append(copy);
      fadeIn(kicker, a + 0.4);
      fadeIn(display, a + 0.7, 0.8);
      fadeIn(lede, a + 1.2, 0.8);

      // Equation
      const eq = el('div', { class: 'eq' });
      const sup = (text) => `<sup style="font-size:.6em">${text}</sup>`;
      const faint = (text) => `<span style="color:var(--gl-ink-faint)">${text}</span>`;
      eq.innerHTML =
        `σ = ∫ d${sup('3')}k d${sup('3')}l <span style="display:inline-block;transform:translateY(2px)">∑</span><span class="sub">cuts</span> ${faint('[')}` +
        [1, 2, 3, 4]
          .map((i) => `<span class="term-g" data-c="${i}"><i>G</i><span class="sub">${i}</span></span>`)
          .join(` ${faint('+')} `) +
        `${faint(']')}(k, l)`;
      const fin = el('div', { class: 'fin' });
      sc.append(eq, fin);
      fadeIn(eq, a + 2.2, 0.8);
      const terms = [...eq.querySelectorAll('.term-g')];

      // Diagram
      const box = el('div', { class: 'diagram' });
      const g = svg('svg', { viewBox: '0 0 1000 900' });
      box.append(g);
      sc.append(box);
      const Lv = [160, 450];
      const Rv = [840, 450];
      const Bv = [500, 190];
      const Cv = [500, 710];
      const ext = [
        svg('path', { class: 'edge photon', d: wavyPath(20, 330, Lv[0], Lv[1], 8, 5) }),
        svg('path', { class: 'edge photon', d: wavyPath(20, 570, Lv[0], Lv[1], 8, 5) }),
        svg('path', { class: 'edge photon', d: wavyPath(Rv[0], Rv[1], 980, 330, 8, 5) }),
        svg('path', { class: 'edge photon', d: wavyPath(Rv[0], Rv[1], 980, 570, 8, 5) }),
      ];
      const LB = fermion(Lv[0], Lv[1], Bv[0], Bv[1]);
      const LC = fermion(Cv[0], Cv[1], Lv[0], Lv[1]);
      const BR = fermion(Bv[0], Bv[1], Rv[0], Rv[1]);
      const CR = fermion(Rv[0], Rv[1], Cv[0], Cv[1]);
      const G = svg('path', { class: 'edge', d: curlyPath(Bv[0], Bv[1] + 22, Cv[0], Cv[1] - 22, 11, 8) });
      const verts = [Lv, Rv, Bv, Cv].map(([x, y]) => svg('circle', { class: 'vertex', cx: x, cy: y, r: 9 }));
      const labels = [
        svg('text', { class: 'extlabel', x: 30, y: 300 }, 'p₁'),
        svg('text', { class: 'extlabel', x: 30, y: 620 }, 'p₂'),
        svg('text', { class: 'extlabel', x: 940, y: 300 }, 'p₁'),
        svg('text', { class: 'extlabel', x: 940, y: 620 }, 'p₂'),
        svg('text', { class: 'elabel', x: 300, y: 295 }, 'q'),
        svg('text', { class: 'elabel', x: 300, y: 630 }, 'q̄'),
        svg('text', { class: 'elabel', x: 540, y: 460 }, 'g'),
      ];
      ext.forEach((p) => g.append(p));
      g.append(LB.g, LC.g, BR.g, CR.g, G, ...verts, ...labels);
      const internal = [LB.line, LC.line, BR.line, CR.line, G];
      [...ext, ...internal].forEach((p, i) => drawPath(p, a + 1.4 + i * 0.12, 0.9, ease.out));
      [LB.head, LC.head, BR.head, CR.head, ...verts, ...labels].forEach((n) => fadeIn(n, a + 2.4, 0.4, 0));

      // loop-momentum marker circulating on the outer loop
      const loopPath = svg('path', {
        d: `M${Lv[0]} ${Lv[1]} L${Bv[0]} ${Bv[1]} L${Rv[0]} ${Rv[1]} L${Cv[0]} ${Cv[1]} Z`,
        fill: 'none',
        stroke: 'none',
      });
      g.append(loopPath);
      const kdot = svg('circle', { r: 11, fill: '#b893c7' });
      const klabel = svg('text', { class: 'klabel' }, 'k');
      g.append(kdot, klabel);
      // The gluon rung carries the second loop momentum l, flowing from the top vertex down.
      const ldot = svg('circle', { r: 9, fill: '#76d3df' });
      const llabel = svg('text', { class: 'klabel', fill: '#76d3df' }, 'l');
      g.append(ldot, llabel);
      let loopLen = null;
      tl.add(
        a + 2.6,
        b - a - 2.6,
        (_, local, t) => {
          if (loopLen === null) loopLen = loopPath.getTotalLength();
          const on = t >= a + 2.6;
          [kdot, klabel, ldot, llabel].forEach((n) => (n.style.opacity = on ? 1 : 0));
          const pt = loopPath.getPointAtLength((local * 190) % loopLen);
          kdot.setAttribute('cx', pt.x);
          kdot.setAttribute('cy', pt.y);
          klabel.setAttribute('x', pt.x + 18);
          klabel.setAttribute('y', pt.y - 14);
          const ly = lerp(Bv[1] + 40, Cv[1] - 40, (local / 3.2) % 1);
          ldot.setAttribute('cx', Bv[0]);
          ldot.setAttribute('cy', ly);
          llabel.setAttribute('x', Bv[0] - 44);
          llabel.setAttribute('y', ly + 10);
        },
        ease.linear,
      );

      // Cuts
      const cuts = [
        { pts: [[325, 90], [325, 810]], hot: [LB.line, LC.line], label: 'cut 1 · virtual correction', tag: [325, 62] },
        { pts: [[675, 90], [675, 810]], hot: [BR.line, CR.line], label: 'cut 2 · virtual correction', tag: [675, 62] },
        {
          pts: [[240, 90], [325, 320], [500, 450], [675, 580], [760, 810]],
          hot: [LB.line, G, CR.line],
          label: 'cut 3 · real emission',
          tag: [240, 62],
        },
        {
          pts: [[240, 810], [325, 580], [500, 450], [675, 320], [760, 90]],
          hot: [LC.line, G, BR.line],
          label: 'cut 4 · real emission',
          tag: [760, 62],
        },
      ];
      const defs = svg('defs');
      g.append(defs);
      // One shared caption names the active cut; a small numeral stays at each cut's top end.
      const caption = svg('text', { class: 'cutlabel', x: 500, y: 20, 'text-anchor': 'middle' }, '');
      caption.style.fontSize = '30px';
      g.append(caption);
      const cutStart = a + 4.2;
      const per = 2.3;
      cuts.forEach((cut, i) => {
        const clipId = `cutclip${i}`;
        const rect = svg('rect', { x: 0, y: 0, width: 1000, height: 0 });
        defs.append(svg('clipPath', { id: clipId }, rect));
        const path = svg('path', { class: 'cutline', d: smoothPath(cut.pts), 'clip-path': `url(#${clipId})` });
        const tag = svg('text', { class: 'cutlabel', x: cut.tag[0], y: cut.tag[1], 'text-anchor': 'middle' }, String(i + 1));
        tag.style.opacity = 0;
        g.append(path, tag);
        const t0 = cutStart + i * per;
        const upward = cut.pts[0][1] > cut.pts[cut.pts.length - 1][1];
        tl.add(
          t0,
          b - t0,
          (_, local, t) => {
            const draw = clamp(local / 0.7);
            const active = t >= t0 && t < t0 + per;
            const ghostPhase = t >= cutStart + 4 * per;
            const hgt = 900 * ease.out(draw);
            rect.setAttribute('height', hgt);
            rect.setAttribute('y', upward ? 900 - hgt : 0);
            path.style.opacity = active ? 1 : ghostPhase ? 0.6 : 0.18;
            path.style.stroke = ghostPhase ? 'var(--gl-accent)' : 'var(--gl-cut)';
            tag.style.opacity = active ? clamp(local / 0.4) : ghostPhase ? 0.9 : 0.35;
            tag.style.fill = ghostPhase ? 'var(--gl-accent)' : 'var(--gl-cut)';
            if (active) {
              caption.textContent = cut.label;
              caption.style.opacity = clamp(local / 0.4);
              caption.style.fill = 'var(--gl-cut)';
            }
            if (i === cuts.length - 1 && ghostPhase) {
              caption.textContent = 'all four cuts localized in the same momenta k and l';
              caption.style.opacity = clamp((t - (cutStart + 4 * per)) / 0.5);
              caption.style.fill = 'var(--gl-accent)';
            }
            const term = terms[i];
            term.style.color = active ? 'var(--gl-cut)' : ghostPhase ? 'var(--gl-ink-strong)' : 'var(--gl-ink-muted)';
            if (active) cut.hot.forEach((e) => e.classList.add('hot'));
            else cut.hot.forEach((e) => e.classList.remove('hot'));
          },
          ease.linear,
        );
      });
      // Reset all hot edges during the final ghost phase so the loop settles.
      tl.add(cutStart + 4 * per, 0.01, () => internal.forEach((e) => e.classList.remove('hot')), ease.linear);
      const finAt = cutStart + 4 * per + 0.3;
      typeText(fin, 'finite at every k and l', finAt, 0.9);
    }

    // -- CLI help ------------------------------------------------------------
    {
      const [a, b] = T.cli;
      const sc = scene('s-cli', 'The command line', a, b);
      const head = el('div', { class: 'head', style: 'position:absolute;left:90px;top:60px;display:flex;gap:28px;align-items:baseline' });
      const kicker = el('div', { class: 'kicker' }, 'The command line');
      const display = el('div', { class: 'display', style: 'font-size:48px' }, 'One stateful application.');
      head.append(kicker, display);
      sc.append(head);
      fadeIn(kicker, a + 0.3);
      fadeIn(display, a + 0.5);

      const term = new Term(sc, { title: 'zsh — gammaloop' });
      term.root.style.cssText = 'left:90px;top:150px;width:1740px;height:800px;font-size:22px';
      let t = a + 0.9;
      t = term.type(t, 1.6, [C('nix profile add '), S('github:alphal00p/gammaloop#gammaloop')]);
      t = term.type(t + 0.5, 0.9, [C('gammaloop '), F('--help')]);
      const cmd = (name, desc) => [seg('cmd', `  ${name.padEnd(13)}`), O(desc)];
      const lines = [
        'Generate cross-section or amplitude integrands, integrate them, and manage a persistent state.',
        '',
        [O('Usage: '), C('gammaloop '), F('[OPTIONS] [COMMAND]')],
        '',
        [B('Commands:')],
        cmd('display', 'Inspect models, processes, integrands, settings, and named session data'),
        cmd('set', 'Change global, runtime, model, or process settings'),
        cmd('import', 'Import models or previously generated graphs into the active state'),
        cmd('save', 'Export graph data, standalone evaluators, schemas, or the current state'),
        cmd('run', 'Execute a stored command block or an inline command list'),
        [seg('line-hi', ''), ...cmd('generate', 'Generate cross-section or amplitude integrands from a process specification')],
        [seg('line-hi', ''), ...cmd('integrate', 'Integrate one or more generated process integrands')],
        cmd('inspect', 'Inspect integrand contributions at one phase-space point or momentum configuration'),
        cmd('approach', 'Sample an integrand while approaching a phase-space point along selected scaling axes'),
        cmd('evaluate', 'Evaluate a vacuum amplitude analytically or numerically with Vakint'),
        cmd('renormalize', 'Compute and export ultraviolet renormalization contributions for an amplitude'),
        cmd('profile', 'Profile ultraviolet or infrared scaling limits of a generated integrand'),
        cmd('3Drep', 'Validate a graph or build its diagnostic three-dimensional energy representation'),
        cmd('quit', 'End the interactive GammaLoop session without executing further commands'),
        '',
        [B('Options:')],
        [F('  -s, --state-folder <PATH>  '), O('Path to the state folder [default: ./gammaloop_state]')],
        [F('      --clean-state          '), O('Remove the resolved state folder before startup')],
        [F('      --read-only-state      '), O('Never write into the state folder')],
        [F('  -o, --save-on-exit         '), O('Save state to file after each call')],
      ];
      t = term.lines(t + 0.35, 0.13, lines);
      term.prompt(t + 0.3);
    }

    // -- Generate ------------------------------------------------------------
    {
      const [a, b] = T.gen;
      const sc = scene('s-gen', 'Process generation', a, b);
      const head = el('div', { class: 'head' });
      const kicker = el('div', { class: 'kicker' }, 'Process generation');
      const display = el('div', { class: 'display' }, 'From process to integrand.');
      head.append(kicker, display);
      sc.append(head);
      fadeIn(kicker, a + 0.3);
      fadeIn(display, a + 0.5);

      const term = new Term(sc, { title: 'zsh — gammaloop' });
      term.root.style.fontSize = '21px';
      let t = a + 0.9;
      t = term.type(t, 4.2, [
        C('gammaloop '),
        F('--clean-state -s '),
        C('./gammaloop_state/aa_aa '),
        C('run '),
        F('-c '),
        S("'\n    import model sm-default.json;\n    generate amp a a > a a | a t t~ QCD==0 QED==4 [{1}] -p aa_aa -i 1L;\n    save dot;\n    display processes;\n    quit -o'"),
      ]);
      const tModel = t + 0.5;
      term.lines(tModel, 0.2, [[O('Loaded model '), HI('sm-default.json'), O(' (Standard Model, Feynman gauge).')]]);
      const tGen = tModel + 1.4;
      term.lines(tGen, 0.35, [
        [O('Generated '), HI('6'), O(' amplitude graphs for '), C('a a > a a'), O(' (one-loop top-quark boxes).')],
        [O('Numerator grouping: '), HI('3'), O(' representatives in one group · master '), C('GL0')],
      ]);
      const tTable = tGen + 1.6;
      const cols = [
        ['integrand', 10],
        ['graph', 7],
        ['evals', 6],
        ['expr build', 11],
        ['spenso', 8],
        ['symbolica', 10],
        ['compile', 8],
      ];
      const tableLine = (cells, cls) => [
        O('  '),
        ...cells.map((cell, k) => seg(cls && k < 2 ? 'cmd' : cls ? 'out' : 'out', String(cell).padEnd(cols[k][1] + 2))),
      ];
      term.lines(tTable, 0.28, [
        [B('Integrand generation summary')],
        tableLine(cols.map((c) => c[0])),
        tableLine(['1L', 'GL0', 4, '0.41 s', '0.12 s', '2.31 s', '6.8 s'], true),
        tableLine(['1L', 'GL2', 4, '0.39 s', '0.11 s', '2.27 s', '6.5 s'], true),
        tableLine(['1L', 'GL4', 4, '0.40 s', '0.12 s', '2.30 s', '6.7 s'], true),
      ]);
      const tDot = tTable + 0.28 * 5 + 0.6;
      term.lines(tDot, 0.25, [[O('Saved DOT graphs under '), C('processes/aa_aa/1L/graphs/'), O(' (GL0, GL2, GL4).')]]);
      const tProc = tDot + 1.0;
      term.lines(tProc, 0.22, [
        [B('Processes')],
        [O('  '), HI('#0'), O('  '), C('aa_aa'), O('   amplitude   a a > a a | a t t~   QCD==0 QED==4 [{1}]   '), C('1L')],
      ]);
      const tSave = tProc + 1.2;
      term.lines(tSave, 0.25, [[OK('State saved'), O(' to '), C('./gammaloop_state/aa_aa'), O(' · run.toml replays it.')]]);
      term.prompt(tSave + 0.5);

      // Right panel: three box graphs, then the persisted state tree.
      const graphs = el('div', { class: 'graphs', style: 'left:1200px;width:660px' });
      const gs = svg('svg', { viewBox: '0 0 660 800' });
      graphs.append(gs);
      sc.append(graphs);
      const orders = [
        { name: 'GL0', master: true, labels: ['₁', '₂', '₃', '₄'] },
        { name: 'GL2', master: false, labels: ['₁', '₂', '₄', '₃'] },
        { name: 'GL4', master: false, labels: ['₁', '₃', '₂', '₄'] },
      ];
      orders.forEach((o, i) => {
        const cx = 110 + i * 220;
        const cy = 130;
        const s = 36; // half box size
        const corners = [
          [cx - s, cy - s],
          [cx + s, cy - s],
          [cx + s, cy + s],
          [cx - s, cy + s],
        ];
        const group = svg('g');
        gs.append(group);
        const edges = [];
        for (let k = 0; k < 4; k++) {
          const [x1, y1] = corners[k];
          const [x2, y2] = corners[(k + 1) % 4];
          const f = fermion(x1, y1, x2, y2);
          f.line.style.strokeWidth = '4';
          f.head.setAttribute('transform', f.head.getAttribute('transform') + ' scale(0.75)');
          group.append(f.g);
          edges.push(f.line, f.head);
        }
        const legLen = 46;
        const dirs = [
          [-1, -1],
          [1, -1],
          [1, 1],
          [-1, 1],
        ];
        corners.forEach(([x, y], k) => {
          const [dx, dy] = dirs[k];
          const ex = x + (dx * legLen) / Math.SQRT2;
          const ey = y + (dy * legLen) / Math.SQRT2;
          const leg = svg('path', { class: 'edge photon', d: wavyPath(x, y, ex, ey, 5, 3) });
          leg.style.strokeWidth = '4';
          const lab = svg('text', { class: 'label', x: ex + dx * 16, y: ey + dy * 12 + 8, 'text-anchor': 'middle' }, `γ${o.labels[k]}`);
          group.append(leg, lab);
          edges.push(leg, lab);
          group.append(svg('circle', { class: 'vertex', cx: x, cy: y, r: 5 }));
        });
        const name = svg('text', { class: `label ${o.master ? 'master' : ''}`, x: cx, y: cy + 118, 'text-anchor': 'middle' }, o.master ? `${o.name} · master` : o.name);
        group.append(name);
        group.style.opacity = 0;
        const t0 = tGen + 0.2 + i * 0.35;
        tl.add(t0, 0.8, (p) => {
          group.style.opacity = p;
          group.style.transform = `translateY(${(1 - p) * 20}px)`;
        }, ease.out);
        edges.filter((e) => e.tagName === 'path').forEach((e, j) => drawPath(e, t0 + j * 0.05, 0.6, ease.out));
      });
      // DOT excerpt
      const dot = el('div', {
        class: 'tree',
        style: 'left:1200px;top:430px;font-size:18px;line-height:1.4;color:var(--gl-ink-muted)',
      });
      sc.append(dot);
      const dotText = `digraph GL0 {\n  0:4 -> 2:5 [id=4 lmb_id="0" name="e4" particle="t"];\n  3:6 -> 0:7 [id=5 name="e5" particle="t"];\n  2:8 -> 1:9 [id=6 name="e6" particle="t"];\n  1:10 -> 3:11 [id=7 name="e7" particle="t"];\n  exte0 -> 3:0 [id=0 particle="a" pin="x:@-left"];\n  …\n}`;
      typeText(dot, dotText, tDot + 0.2, 1.3);
      const tree = el('div', { class: 'tree', style: 'left:1200px;top:640px;font-size:20px;line-height:1.5' });
      sc.append(tree);
      tree.style.opacity = 0;
      fadeIn(tree, tSave + 0.1, 0.5, 10);
      const treeLines = [
        ['dir', 'gammaloop_state/aa_aa/'],
        ['new', '├─ run.toml'],
        ['new', '├─ state_manifest.toml'],
        ['new', '├─ global_settings.toml'],
        ['new', '├─ default_runtime_settings.toml'],
        ['dir', '└─ processes/aa_aa/'],
        ['dir', '   └─ 1L/'],
        ['', '      ├─ graphs/  GL0.dot  GL2.dot  GL4.dot'],
        ['', '      └─ compiled evaluator'],
      ];
      tl.add(
        tSave + 0.3,
        1.6,
        (p) => {
          const k = Math.round(treeLines.length * p);
          tree.innerHTML = treeLines
            .slice(0, k)
            .map(([cls, text]) => (cls ? `<span class="${cls}">${esc(text)}</span>` : esc(text)))
            .join('\n');
        },
        ease.linear,
      );
    }

    // -- Integrate -----------------------------------------------------------
    {
      const [a, b] = T.int;
      const sc = scene('s-int', 'Monte Carlo integration', a, b);
      const head = el('div', { class: 'head' });
      const kicker = el('div', { class: 'kicker' }, 'Monte Carlo integration');
      const display = el('div', { class: 'display' }, 'Watch it converge.');
      head.append(kicker, display);
      sc.append(head);
      fadeIn(kicker, a + 0.3);
      fadeIn(display, a + 0.5);

      const term = new Term(sc, { title: 'zsh — gammaloop' });
      let t = a + 0.8;
      t = term.type(t, 3.0, [
        C('gammaloop '),
        F('-s '),
        C('./gammaloop_state/aa_aa '),
        C('run '),
        F('-c '),
        S("'integrate -p aa_aa -i 1L --n-cores 8 --target 3.63182591790267920e-6 0.0'"),
      ]);
      const tTui = t + 0.4;

      // Dashboard model
      const target = 3.6318259179026792e-6;
      const NIT = 12;
      const next = rng(1337);
      const iters = [];
      for (let i = 1; i <= NIT; i++) {
        // The adaptive grid trains over the first iterations, then the per-iteration error plateaus.
        const rel = lerp(0.012, 0.0024, Math.min(i - 1, 5) / 5);
        const z = gaussian(next);
        iters.push({ n: 1e6, rel, est: target * (1 + rel * z), err: target * rel });
      }
      const cumulative = [];
      let sw = 0;
      let swx = 0;
      iters.forEach((it, i) => {
        const w = 1 / (it.err * it.err);
        sw += w;
        swx += w * it.est;
        const mean = swx / sw;
        const err = Math.sqrt(1 / sw);
        let chi = 0;
        for (let j = 0; j <= i; j++) chi += Math.pow((iters[j].est - mean) / iters[j].err, 2);
        cumulative.push({ mean, err, chi: i > 0 ? chi / i : 0 });
      });
      const perIter = 0.95;
      const tuiDur = NIT * perIter + 1.0;

      const tui = el('div', { class: 'tui' });
      const tabs = el(
        'div',
        { class: 'tabs' },
        el('span', { class: 'tab on' }, 'Overview'),
        el('span', { class: 'tab' }, 'Discrete'),
        el('span', { class: 'tab' }, 'Max Weight'),
        el('span', { class: 'tab', style: 'margin-left:auto;color:var(--gl-ink-muted)' }, 'aa_aa@1L · 8 cores · ratatui'),
      );
      const grid = el('div', { class: 'grid' });
      const summary = el('div', { class: 'panel' }, el('span', { class: 'ptitle' }, 'Results summary · aa_aa@1L'));
      const summaryBody = el('div', { class: 'kv tnum' });
      summary.append(summaryBody);
      const progress = el('div', { class: 'panel' }, el('span', { class: 'ptitle' }, 'Iteration progress'));
      const progressBody = el('div', { class: 'kv tnum' });
      const bar = el('div', { class: 'bar' }, el('i'));
      const stab = el('div', { class: 'stab tnum' });
      progress.append(progressBody, bar, stab);
      const chart = el('div', { class: 'panel chart' }, el('span', { class: 'ptitle' }, 're [+] · deviation from target · ±3 % window'));
      const cs = svg('svg', { viewBox: '0 0 1660 330', preserveAspectRatio: 'none' });
      chart.append(cs);
      grid.append(summary, progress, chart);
      const hints = [
        ['p', 'phase'],
        ['g', 'history'],
        ['s', 'sort'],
        ['r / c / w', 'rel err, chi^2, mwi'],
        ['a', '|re|'],
        ['?', 'help'],
        ['x / Ctrl-C', 'training re'],
      ];
      const foot = el('div', { class: 'foot', html: hints.map(([k, v]) => `<span><b>${k}</b> ${v}</span>`).join('') });
      tui.append(tabs, grid, foot);

      // chart statics
      const X0 = 70;
      const X1 = 1630;
      const Y0 = 20;
      const Y1 = 290;
      const yOf = (dev) => lerp((Y0 + Y1) / 2, Y0, dev / 3); // dev in percent
      const xOf = (i) => lerp(X0, X1, (i - 0.5) / NIT);
      const axis = svg('g');
      [-3, -2, -1, 0, 1, 2, 3].forEach((v) => {
        axis.append(svg('line', { x1: X0, x2: X1, y1: yOf(v), y2: yOf(v), stroke: v === 0 ? '#b893c7' : '#382d3b', 'stroke-width': v === 0 ? 2 : 1, 'stroke-dasharray': v === 0 ? '8 8' : 'none' }));
        axis.append(svg('text', { x: X0 - 12, y: yOf(v) + 7, 'text-anchor': 'end', fill: '#8d818f', 'font-size': 18, 'font-family': 'IBM Plex Mono, monospace' }, `${v > 0 ? '+' : ''}${v}%`));
      });
      for (let i = 1; i <= NIT; i++) axis.append(svg('text', { x: xOf(i), y: Y1 + 28, 'text-anchor': 'middle', fill: '#8d818f', 'font-size': 18, 'font-family': 'IBM Plex Mono, monospace' }, String(i)));
      axis.append(svg('text', { x: X1, y: yOf(0) - 10, 'text-anchor': 'end', fill: '#b893c7', 'font-size': 18, 'font-family': 'IBM Plex Mono, monospace' }, 'target 3.6318e-6'));
      cs.append(axis);
      const band = svg('path', { fill: '#b893c7', opacity: 0.18 });
      const meanLine = svg('path', { fill: 'none', stroke: '#f8effa', 'stroke-width': 3 });
      const pointsG = svg('g');
      cs.append(band, meanLine, pointsG);
      const pointNodes = iters.map((it, i) => {
        const g = svg('g');
        const dev = ((it.est - target) / target) * 100;
        const e = it.rel * 100;
        g.append(svg('line', { x1: xOf(i + 1), x2: xOf(i + 1), y1: yOf(dev - e), y2: yOf(dev + e), stroke: '#d8b9e3', 'stroke-width': 2 }));
        g.append(svg('line', { x1: xOf(i + 1) - 8, x2: xOf(i + 1) + 8, y1: yOf(dev - e), y2: yOf(dev - e), stroke: '#d8b9e3', 'stroke-width': 2 }));
        g.append(svg('line', { x1: xOf(i + 1) - 8, x2: xOf(i + 1) + 8, y1: yOf(dev + e), y2: yOf(dev + e), stroke: '#d8b9e3', 'stroke-width': 2 }));
        g.append(svg('circle', { cx: xOf(i + 1), cy: yOf(dev), r: 6, fill: '#e7b66f' }));
        g.style.opacity = 0;
        pointsG.append(g);
        return g;
      });

      const fmtE = (x, d = 4) => x.toExponential(d).replace('e-', 'e-').replace('e+', 'e+');
      const kv = (rows) => rows.map(([k, v, cls = '']) => `<span class="k">${k}</span><span class="v ${cls}">${v}</span>`).join('');
      const pad = (n) => String(Math.floor(n)).padStart(2, '0');

      term.mount(tTui, tui, tTui + tuiDur);
      tl.add(
        tTui,
        tuiDur,
        (_, local) => {
          const raw = local / perIter;
          const done = Math.min(NIT, Math.floor(raw));
          const frac = done >= NIT ? 1 : raw - done;
          const cur = Math.min(NIT, done + 1);
          const samplesDone = done * 1e6 + (done >= NIT ? 0 : frac * 1e6);
          const c = done > 0 ? cumulative[done - 1] : null;
          const rel = c ? c.err / Math.abs(c.mean) : null;
          const sigma = c ? (c.mean - target) / c.err : null;
          summaryBody.innerHTML = kv([
            ['re [+]', c ? `${fmtE(c.mean)} ± ${fmtE(c.err, 1)}` : '—'],
            ['rel err', c ? `${(rel * 100).toFixed(3)} %` : '—', c && rel < 1e-3 ? 'pos' : ''],
            ['chi^2 / dof', c && done > 1 ? c.chi.toFixed(2) : '—'],
            ['target', `${fmtE(target)}  (goal rel err 1.0e-3)`],
            ['Δ target', c ? `${sigma >= 0 ? '+' : ''}${sigma.toFixed(2)} σ` : '—', c && Math.abs(sigma) < 2 ? 'pos' : 'neg'],
            ['phase', 'real · summed orientations · LMB channels'],
          ]);
          const rate = 1.85e6 + 4e4 * Math.sin(local * 3);
          const remaining = c ? Math.max(0, ((rel / 1e-3) ** 2 - 1) * (done * 1e6)) / rate : 220;
          progressBody.innerHTML = kv([
            ['Iteration', `${cur} / ∞   (n_start 1e6, n_increase 0)`],
            ['# samples total', fmtE(samplesDone, 2)],
            ['#samples/s', `${fmtE(rate, 2)}  (8 cores)`],
            ['ETA to ⟨I⟩ target', done >= NIT ? 'reached' : `00:${pad(remaining / 60)}:${pad(remaining % 60)}`, done >= NIT ? 'pos' : ''],
          ]);
          bar.firstElementChild.style.width = `${(done >= NIT ? 1 : frac) * 100}%`;
          const f128 = 0.88 + 0.05 * Math.sin(local * 0.7);
          stab.textContent = `stability: f64 ${(100 - f128).toFixed(2)}% · f128 ${f128.toFixed(2)}% · arb 0.00%`;
          // chart
          pointNodes.forEach((g, i) => (g.style.opacity = i < done ? 1 : 0));
          if (done > 0) {
            const up = [];
            const dn = [];
            let line = '';
            for (let i = 0; i < done; i++) {
              const cm = cumulative[i];
              const dev = ((cm.mean - target) / target) * 100;
              const e = (cm.err / target) * 100;
              up.push([xOf(i + 1), yOf(dev + e)]);
              dn.push([xOf(i + 1), yOf(dev - e)]);
              line += `${i ? 'L' : 'M'}${xOf(i + 1)} ${yOf(dev)} `;
            }
            meanLine.setAttribute('d', line);
            band.setAttribute('d', done > 1 ? `M${up.map((p) => p.join(' ')).join(' L')} L${dn.reverse().map((p) => p.join(' ')).join(' L')} Z` : '');
          } else {
            meanLine.setAttribute('d', '');
            band.setAttribute('d', '');
          }
        },
        ease.linear,
      );
      const last = cumulative[NIT - 1];
      const tEnd = tTui + tuiDur;
      const relLast = ((last.err / last.mean) * 100).toFixed(3);
      const sigLast = ((last.mean - target) / last.err).toFixed(2);
      term.lines(tEnd, 0.3, [
        [B('Results summary · aa_aa@1L')],
        [O('  re [+]   '), HI(`${fmtE(last.mean)} ± ${fmtE(last.err, 1)}`), O(`   rel err ${relLast} %   chi^2/dof ${last.chi.toFixed(2)}   Δ target ${sigLast} σ`)],
        [OK('Target relative accuracy 1.0e-3 reached'), O(` after ${NIT} iterations (${fmtE(NIT * 1e6, 1)} samples).`)],
        [O('Results and observables written to '), C('./gammaloop_state/aa_aa/integration_workspace/')],
      ]);
      term.prompt(tEnd + 1.5);
    }

    // -- APIs -----------------------------------------------------------------
    {
      const [a, b] = T.api;
      const sc = scene('s-api', 'Rust and Python APIs', a, b);
      const head = el('div', { class: 'head' });
      const kicker = el('div', { class: 'kicker' }, 'Rust and Python');
      const display = el('div', { class: 'display' }, 'One state. Three interfaces.');
      const lede = el('div', { class: 'lede', style: 'position:absolute;left:90px;top:972px;max-width:none' }, 'The CLI, the Rust facade, and the Python package load the same persisted state and run the same commands.');
      head.append(kicker, display);
      sc.append(head, lede);
      fadeIn(kicker, a + 0.3);
      fadeIn(display, a + 0.5);
      fadeIn(lede, a + 8.4);

      const py = el('div', { class: 'pane py' }, el('div', { class: 'bar' }, el('span', { class: 'lang' }, 'Python'), el('span', {}, 'first_gammaloop.py')));
      const rs = el('div', { class: 'pane rs' }, el('div', { class: 'bar' }, el('span', { class: 'lang' }, 'Rust'), el('span', {}, 'src/main.rs')));
      const pyCode = el('div', { class: 'code' });
      const rsCode = el('div', { class: 'code' });
      py.append(pyCode);
      rs.append(rsCode);
      sc.append(py, rs);
      fadeIn(py, a + 1.0, 0.6, 20);
      fadeIn(rs, a + 1.2, 0.6, 20);
      const pySrc = `from gammaloop import GammaLoopAPI

  api = GammaLoopAPI(state_folder="gammaloop_state/bubble")
  api.run("import model scalars-default.json")
  api.run(
      "generate amp scalar_1 > scalar_1 [{1}] "
      "--allowed-vertex-interactions V_3_SCALAR_122 "
      "-p bubble -i one_loop"
  )

  result = api.evaluate_sample(
      [0.1, 0.2, 0.3],
      process_id=0,
      integrand_name="one_loop",
  )
  print(result.integrand_result)`;
      const rsSrc = `use gammaloop_api::commands::evaluate_samples::{
      evaluate_sample, EvaluateSamples,
  };
  use gammaloop_api::{state::CommandHistory, StateLoadOption};

  let mut loaded = StateLoadOption {
      state_folder: Some("gammaloop_state/bubble".into()),
      ..StateLoadOption::default()
  }
  .load()?;

  for raw in [
      "import model scalars-default.json",
      "generate amp scalar_1 > scalar_1 [{1}] \\
       --allowed-vertex-interactions V_3_SCALAR_122 \\
       -p bubble -i one_loop",
  ] {
      let command = CommandHistory::from_raw_string(raw)?;
      loaded.cli_session().execute_command(command)?;
  }

  let result = evaluate_sample(&mut loaded.state, &EvaluateSamples {
      process_id: Some(0),
      integrand_name: Some("one_loop".into()),
      points: arr2(&[[0.1, 0.2, 0.3]]).view(),
      use_arb_prec: false,
      minimal_output: false,
      momentum_space: false,
      return_generated_events: None,
      integrator_weights: None,
      discrete_dims: None,
      graph_names: None,
      orientations: None,
  })?;
  println!("{result}");`;
      typeCode(pyCode, pySrc, 'python', a + 1.6, 5.2);
      typeCode(rsCode, rsSrc, 'rust', a + 1.9, 6.4);

      // Hub: shared state and connectors
      const hub = el('div', { class: 'chip hub', style: 'left:1100px;top:92px' }, 'gammaloop_state/bubble/');
      const cli = el('div', { class: 'chip', style: 'left:1100px;top:160px' }, '$ gammaloop -s gammaloop_state/bubble run -c "display processes"');
      sc.append(hub, cli);
      fadeIn(hub, a + 7.6, 0.6, 12);
      fadeIn(cli, a + 8.0, 0.6, 12);
      const link = svg('svg', { class: 'link', viewBox: `0 0 ${W} ${H}` });
      sc.append(link);
      const stroke = { fill: 'none', stroke: '#b893c7', 'stroke-width': 3, 'stroke-dasharray': '10 8' };
      const l1 = svg('path', { d: 'M1130 144 C1080 200, 700 190, 512 250', ...stroke });
      const l2 = svg('path', { d: 'M1290 144 C1330 200, 1390 220, 1410 250', ...stroke });
      link.append(l1, l2);
      [l1, l2].forEach((p, i) => {
        p.style.opacity = 0;
        tl.add(a + 8.2 + i * 0.2, 0.6, (q) => (p.style.opacity = q), ease.out);
      });
    }

    // -- Ecosystem -------------------------------------------------------------
    {
      const [a, b] = T.eco;
      const sc = scene('s-eco', 'The αLoop ecosystem', a, b);
      const head = el('div', { class: 'head' });
      const kicker = el('div', { class: 'kicker', style: 'text-transform:none;letter-spacing:0.06em' }, 'The αLoop ecosystem');
      const display = el('div', { class: 'display' }, 'Five crates, one loop.');
      head.append(kicker, display);
      sc.append(head);
      fadeIn(kicker, a + 0.3);
      fadeIn(display, a + 0.5);

      const g = svg('svg', { viewBox: `0 0 ${W} ${H}` });
      sc.append(g);
      const cx = 960;
      const cy = 600;
      const nodes = [
        { id: 'gammaloop', title: 'GammaLoop', sub: 'Local Unitarity cross-sections · CLI, Rust, Python', pkg: 'gammaloop-api', x: cx, y: cy, core: true },
        { id: 'linnet', title: 'Linnet', sub: 'Half-edge graphs, Linnest layout, Clinnet rendering', pkg: 'linnet', x: cx - 470, y: cy - 220 },
        { id: 'spenso', title: 'Spenso', sub: 'Typed tensors and executable tensor networks', pkg: 'spenso', x: cx + 470, y: cy - 220 },
        { id: 'idenso', title: 'Idenso', sub: 'Symbolic tensor identities and rewrites', pkg: 'idenso', x: cx + 470, y: cy + 220 },
        { id: 'vakint', title: 'Vakint', sub: 'Vacuum-integral matching and evaluation', pkg: 'vakint', x: cx - 470, y: cy + 220 },
      ];
      const NW = 360;
      const NH = 150;
      // The loop through the four libraries, drawn as a photon line.
      const ring = svg('path', {
        class: 'link',
        d: [
          wavyPath(nodes[1].x + NW / 2, nodes[1].y, nodes[2].x - NW / 2, nodes[2].y, 7, 9),
          wavyPath(nodes[2].x, nodes[2].y + NH / 2, nodes[3].x, nodes[3].y - NH / 2, 7, 4),
          wavyPath(nodes[3].x - NW / 2, nodes[3].y, nodes[4].x + NW / 2, nodes[4].y, 7, 9),
          wavyPath(nodes[4].x, nodes[4].y - NH / 2, nodes[1].x, nodes[1].y + NH / 2, 7, 4),
        ].join(' '),
      });
      g.append(ring);
      drawPath(ring, a + 0.9, 2.2, ease.inOut);
      const spokes = nodes.slice(1).map((n) => {
        const dx = n.x - cx;
        const dy = n.y - cy;
        const L = Math.hypot(dx, dy);
        const ux = dx / L;
        const uy = dy / L;
        const p = svg('path', { class: 'edge', d: `M${cx + ux * 120} ${cy + uy * 60} L${n.x - ux * 150} ${n.y - uy * 60}` });
        p.style.strokeWidth = '3';
        g.append(p);
        return p;
      });
      spokes.forEach((p, i) => drawPath(p, a + 2.2 + i * 0.15, 0.6, ease.out));
      nodes.forEach((n, i) => {
        const grp = svg('g');
        grp.append(svg('rect', { class: `node ${n.core ? 'core' : ''}`, x: n.x - NW / 2, y: n.y - NH / 2, width: NW, height: NH, rx: 18 }));
        grp.append(svg('text', { class: 'ntitle', x: n.x, y: n.y - 18, 'text-anchor': 'middle' }, n.title));
        grp.append(svg('text', { class: 'npkg', x: n.x, y: n.y + 14, 'text-anchor': 'middle' }, n.pkg));
        const sub = svg('text', { class: 'nsub', x: n.x, y: n.y + 46, 'text-anchor': 'middle' });
        const words = n.sub.split(' ');
        const half = Math.ceil(words.length / 2);
        if (n.sub.length > 34) {
          sub.append(svg('tspan', { x: n.x, dy: 0 }, words.slice(0, half).join(' ')));
          sub.append(svg('tspan', { x: n.x, dy: 24 }, words.slice(half).join(' ')));
          sub.setAttribute('y', n.y + 40);
        } else sub.textContent = n.sub;
        grp.append(sub);
        grp.style.opacity = 0;
        g.append(grp);
        const t0 = n.core ? a + 0.6 : a + 1.6 + (i - 1) * 0.45;
        tl.add(t0, 0.6, (p) => {
          grp.style.opacity = p;
          grp.style.transform = `translateY(${(1 - p) * 16}px)`;
        }, ease.out);
      });
      const legend = el('div', { class: 'legend' }, 'Written in Rust · symbolic algebra with Symbolica · diagrams rendered with Typst');
      sc.append(legend);
      fadeIn(legend, a + 4.2, 0.8, 12);
    }

    // -- Outro -----------------------------------------------------------------
    {
      const [a, b] = T.out;
      const sc = scene('s-out', 'Outro', a, b);
      const logo = logoNode();
      sc.append(logo);
      const paths = [...logo.querySelectorAll('path.ink')];
      paths.forEach((p, i) => drawPath(p, a + 0.2 + i * 0.06, 0.8, ease.out));
      const lines = el('div', { class: 'lines' });
      const big = el('div', { class: 'big' }, 'GammaLoop');
      const u1 = el('div', { class: 'url' }, 'alphal00p.github.io/gammaloop');
      const u2 = el('div', { class: 'url' }, 'github.com/alphal00p/gammaloop');
      const small = el('div', { class: 'small' }, 'The αLoop collaboration · Research associated with Local Unitarity is supported by the Swiss National Science Foundation');
      const credit = el('div', { class: 'small' }, 'Music: “Hard Boiled” by Kevin MacLeod (incompetech.com), Creative Commons Attribution 4.0');
      lines.append(big, u1, u2, small, credit);
      sc.append(lines);
      fadeIn(big, a + 1.0, 0.8, 20);
      fadeIn(u1, a + 1.6, 0.7, 14);
      fadeIn(u2, a + 1.9, 0.7, 14);
      fadeIn(small, a + 2.6, 0.8, 10);
      fadeIn(credit, a + 3.0, 0.8, 10);
    }

    return { tl, scenes, duration: T.out[1] };
  }

  // ================================================================ playback
  // Soundtrack: an <audio> element whose clock drives the animation while sound is on. Audio
  // only starts from a user gesture, so the browser never blocks it.
  class Soundtrack {
    constructor(src) {
      this.src = src;
      this.audio = null;
      this.active = false;
      this.pending = 0;
      this.waiting = false;
    }
    get available() {
      return Boolean(this.src) && typeof Audio === 'function';
    }
    get running() {
      return this.active && this.audio !== null && !this.audio.paused && this.audio.readyState >= 3;
    }
    start(time) {
      if (!this.audio) {
        this.audio = new Audio(this.src);
        this.audio.preload = 'auto';
      }
      this.active = true;
      this.pending = time;
      if (this.audio.readyState >= 1) this.begin();
      else if (!this.waiting) {
        this.waiting = true;
        this.audio.addEventListener(
          'loadedmetadata',
          () => {
            this.waiting = false;
            if (this.active) this.begin();
          },
          { once: true },
        );
      }
    }
    begin() {
      this.audio.currentTime = this.pending;
      this.audio.play().catch(() => {});
    }
    time() {
      return this.audio.currentTime;
    }
    pause() {
      if (this.audio) this.audio.pause();
    }
    stop() {
      this.active = false;
      if (this.audio) this.audio.pause();
    }
  }

  function mount(host, options = {}) {
    const render = options.render === true;
    const wrapper = el('div', { class: `gl-showcase ${render ? 'render' : 'interactive'}` });
    const frame = el('div', { class: 'gl-showcase-frame' });
    const stage = el('div', { class: 'gl-showcase-stage', 'aria-hidden': 'true' });
    stage.append(el('canvas', { class: 'gl-showcase-field' }));
    frame.append(stage);
    wrapper.append(frame);
    host.replaceChildren(wrapper);

    const { tl, scenes, duration } = build(stage);
    const clampT = (t) => clamp(t, 0, duration);
    let current = 0;
    const fmtTime = (s) => `${Math.floor(s / 60)}:${(s % 60).toFixed(1).padStart(4, '0')}`;

    function fit() {
      const width = render ? W : frame.clientWidth || W;
      stage.style.transform = `scale(${width / W})`;
    }
    fit();
    if (!render && typeof ResizeObserver === 'function') new ResizeObserver(fit).observe(frame);

    let updatePlayer = () => {};
    function seek(t) {
      current = clampT(t);
      tl.seek(current);
      updatePlayer();
    }
    const ready = (document.fonts ? document.fonts.ready : Promise.resolve()).then(() => seek(0));
    const controller = { duration, scenes, seek, ready, element: wrapper, play: () => {}, pause: () => {} };
    if (render) return controller;

    // Interactive controls
    const music = new Soundtrack(options.soundtrack);
    let playing = false;
    let origin = 0;
    let soundOn = options.autoplay !== true; // autoplay cannot carry sound, so it starts muted
    let started = false;
    let autoPaused = false;

    const cover = el(
      'button',
      { class: 'gl-showcase-cover', type: 'button', 'aria-label': 'Play the GammaLoop tour' },
      el('span', { class: 'gl-showcase-cover-button', 'aria-hidden': 'true' }, '▶'),
      el('span', { class: 'gl-showcase-cover-label' }, `Play the tour · ${fmtTime(duration)}${music.available ? ' · with sound' : ''}`),
    );
    frame.append(cover);
    const playBtn = el('button', { type: 'button', class: 'gl-showcase-play', 'aria-label': 'Play' }, 'Play');
    const range = el('input', { type: 'range', min: '0', max: String(Math.round(duration * 1000)), value: '0', step: '33', 'aria-label': 'Tour timeline' });
    const sceneEl = el('span', { class: 'gl-showcase-scene' });
    const timeEl = el('span', { class: 'gl-showcase-time' }, `0:00.0 / ${fmtTime(duration)}`);
    const soundBtn = el('button', { type: 'button', class: 'gl-showcase-sound', 'aria-pressed': String(soundOn) }, soundOn ? 'Sound on' : 'Sound off');
    const player = el('div', { class: 'gl-showcase-player' }, playBtn, range, sceneEl, timeEl);
    if (music.available) player.append(soundBtn);
    wrapper.append(player);

    updatePlayer = () => {
      range.value = String(Math.round(current * 1000));
      timeEl.textContent = `${fmtTime(current)} / ${fmtTime(duration)}`;
      const sc = scenes.find((s) => current >= s.start && current < s.end);
      sceneEl.textContent = sc ? sc.name : '';
    };

    function play() {
      if (playing) return;
      playing = true;
      started = true;
      cover.hidden = true;
      playBtn.textContent = 'Pause';
      playBtn.setAttribute('aria-label', 'Pause');
      origin = performance.now() - current * 1000;
      if (soundOn && music.available) music.start(current);
      requestAnimationFrame(tick);
    }
    function pause() {
      if (!playing) return;
      playing = false;
      playBtn.textContent = 'Play';
      playBtn.setAttribute('aria-label', 'Play');
      music.pause();
    }
    function tick(now) {
      if (!playing) return;
      let t = (now - origin) / 1000;
      if (music.running) {
        t = music.time();
        origin = now - t * 1000;
      }
      if (t >= duration) {
        t = 0;
        origin = now;
        if (music.active) music.start(0);
      }
      seek(t);
      requestAnimationFrame(tick);
    }
    function scrub(t) {
      seek(t);
      origin = performance.now() - current * 1000;
      if (playing && music.active) music.start(current);
      else if (!playing) music.stop();
    }
    function toggleSound() {
      soundOn = !soundOn;
      soundBtn.textContent = soundOn ? 'Sound on' : 'Sound off';
      soundBtn.setAttribute('aria-pressed', String(soundOn));
      if (soundOn && playing) music.start(current);
      else if (!soundOn) music.stop();
    }

    // The cover shows a still from the middle of the opening; a first play starts from the top.
    cover.addEventListener('click', () => {
      seek(0);
      play();
    });
    playBtn.addEventListener('click', () => (playing ? pause() : play()));
    soundBtn.addEventListener('click', toggleSound);
    range.addEventListener('input', () => scrub(Number(range.value) / 1000));
    wrapper.addEventListener('keydown', (e) => {
      if (e.target === range && (e.key === 'ArrowLeft' || e.key === 'ArrowRight')) return;
      if (e.key === ' ') {
        e.preventDefault();
        playing ? pause() : play();
      } else if (e.key === 'ArrowRight') scrub(current + (e.shiftKey ? 5 : 1));
      else if (e.key === 'ArrowLeft') scrub(current - (e.shiftKey ? 5 : 1));
      else if (e.key === 'Home') scrub(0);
    });
    // Do not keep animating or playing while the tour is scrolled out of view.
    if (typeof IntersectionObserver === 'function') {
      new IntersectionObserver(
        ([entry]) => {
          if (!entry.isIntersecting && playing) {
            pause();
            autoPaused = true;
          } else if (entry.isIntersecting && autoPaused) {
            autoPaused = false;
            play();
          }
        },
        { threshold: 0.1 },
      ).observe(frame);
    }

    controller.play = play;
    controller.pause = pause;
    ready.then(() => {
      const reduce = globalThis.matchMedia && globalThis.matchMedia('(prefers-reduced-motion: reduce)').matches;
      if (options.autoplay && !reduce) {
        seek(0);
        play();
      } else seek(options.poster ?? 0);
    });
    return controller;
  }

  function autoMount() {
    const render = new URLSearchParams(location.search).has('render');
    document.querySelectorAll('[data-showcase]').forEach((host, index) => {
      const controller = mount(host, {
        render,
        autoplay: host.hasAttribute('data-autoplay'),
        poster: host.dataset.poster ? Number(host.dataset.poster) : 5.6,
        soundtrack: host.dataset.soundtrack || null,
      });
      if (index === 0) globalThis.showcase = controller;
    });
  }
  globalThis.GammaLoopShowcase = { mount };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', autoMount);
  else autoMount();
})();
