(() => {
  /* Scope state to the SVG: identical graph IDs can occur in several outputs.
     Use block comments because notebook iframe wrappers may flatten newlines. */
  for (const svg of document.querySelectorAll('svg[data-linnet-interactive]')) {
    if (svg.linnetSelection !== undefined) continue;
    /* Composite channels own fixed viewports; only their parent controls size. */
    const embedded = !!svg.parentElement?.closest('svg');
    const targets = [...svg.querySelectorAll('[data-linnet-kind]')];
    const selected = { node: new Set(), edge: new Set(), halfedge: new Set() };
    const edgeHighlights = new Map();
    const selection = () => ({
      nodes: [...selected.node].sort((a, b) => a - b),
      edges: [...selected.edge].sort((a, b) => a - b),
      half_edges: [...selected.halfedge].sort((a, b) => a - b),
    });
    Object.defineProperty(svg, 'linnetSelection', { get: selection });
    const bounds = svg.viewBox.baseVal;
    const box = { x: bounds.x, y: bounds.y, width: bounds.width, height: bounds.height };
    /* SVG lengths resolve pt/mm to CSS pixels. Read the authored width before
       replacing it with the responsive viewer width. */
    const authoredScale = embedded ? 1 : svg.width.baseVal.value / box.width;
    /* Enlarge small standalone diagrams to a readable notebook baseline.
       This is 100% for zoom/reset; fixed and embedded views keep their own
       authored scale. Bound enlargement independently of the column width. */
    const previewScale = embedded || Number(svg.dataset.linnetViewportHeight) > 0 ? 1
      : Math.max(1, Math.min(2.25, 450 / (box.width * authoredScale), 300 / (box.height * authoredScale)));
    const naturalScale = authoredScale * previewScale;
    let scale = 1;
    let fitted = true;
    let pinned = false;
    let hoverDismissed = false;
    let dragging = false;
    let suppressClick = false;
    let inspected;
    let resize;
    let resizeFrame = () => {};
    let closePopup = () => {};
    let syncPopupTheme = () => {};
    const lifecycle = new AbortController();
    const frames = [];

    /* Notebook components can put their iframe inside a shadow root or slot.
       Follow that rendered ancestry to reach the actual theme/scroll owner. */
    const parent = element => {
      const container = element.assignedSlot || element.parentElement || element.getRootNode().host;
      if (container) return container;
      try { return element.ownerDocument.defaultView.frameElement; }
      catch (_) { return null; }
    };
    const theme = () => {
      /* Marimo executes scripted HTML in an iframe. Read its notebook theme
         when accessible; standalone SVGs retain the system color preference. */
      let element = parent(svg);
      while (element) {
        const explicit = element.getAttribute('data-theme');
        const jupyter = element.getAttribute('data-jp-theme-light');
        if (explicit === 'dark' || explicit === 'light') return explicit;
        if (jupyter !== null) return jupyter === 'true' ? 'light' : 'dark';
        if (element.classList.contains('dark')) return 'dark';
        if (element.classList.contains('light')) return 'light';
        element = parent(element);
      }
      return null;
    };
    const syncHost = () => {
      if (!svg.isConnected || frames.some(frame => !frame.isConnected)
          || frames[0] && frames[0].contentDocument !== document) {
        closePopup();
        lifecycle.abort();
        observer.disconnect();
        if (resize) resize.disconnect();
        return;
      }
      for (let ancestor = parent(svg); ancestor; ancestor = parent(ancestor)) {
        /* Reserve notebook scrollbars in both inline and iframe displays.
           Opening the inspector must not move the hovered target sideways:
           otherwise pointerout closes it and pointerover immediately reopens it.
           Hosts may enable scrolling after mounting, or replace inline styles. */
        const style = ancestor.ownerDocument.defaultView.getComputedStyle(ancestor);
        if (['auto', 'scroll'].includes(style.overflowY) && style.scrollbarGutter === 'auto') {
          ancestor.style.scrollbarGutter = 'stable';
        }
      }
      const value = theme();
      if (value) svg.setAttribute('data-theme', value);
      else svg.removeAttribute('data-theme');
      syncPopupTheme();
      try {
        const frame = window.frameElement;
        if (frame && frame.parentElement && document.body) {
          const style = frame.ownerDocument.defaultView.getComputedStyle(frame.parentElement);
          if (document.body.style.color !== style.color) document.body.style.color = style.color;
          if (document.body.style.fontFamily !== style.fontFamily) document.body.style.fontFamily = style.fontFamily;
          /* Browsers paint an opaque canvas when an iframe's declared scheme
             differs from its document, even when both backgrounds are transparent. */
          const scheme = frame.ownerDocument.defaultView.getComputedStyle(frame).colorScheme;
          if (document.documentElement.style.colorScheme !== scheme) document.documentElement.style.colorScheme = scheme;
        }
      } catch (_) { /* Cross-origin embeddings retain their own text styling. */ }
    };
    const observer = new MutationObserver(syncHost);
    syncHost();
    let ancestor = parent(svg);
    while (ancestor) {
      if (ancestor.localName === 'iframe') frames.push(ancestor);
      observer.observe(ancestor, { attributes: true, childList: true, attributeFilter: ['class', 'style', 'data-theme', 'data-jp-theme-light'] });
      ancestor = parent(ancestor);
    }
    try {
      const frame = window.frameElement;
      if (frame && document.body && !svg.closest('[data-linnet-frame-owner]')) {
        /* A graph-only notebook frame needs no second browser body inset.
           Leave authored page margins and surrounding notebook padding alone. */
        const roots = [...document.body.children].filter(element => !['style', 'script', 'link'].includes(element.localName));
        if (roots.length === 1 && (roots[0] === svg || roots[0].localName === 'figure' && roots[0].contains(svg))
          && getComputedStyle(document.body).margin === '8px') document.body.style.margin = '0';
        /* Measure content rather than the viewport: Marimo's own iframe
           observer grows outputs but deliberately does not shrink them. */
        document.documentElement.style.overflow = 'hidden';
        frame.style.display = 'block';
        resizeFrame = () => {
          if (!svg.isConnected) { resize.disconnect(); return; }
          const style = getComputedStyle(document.body);
          const contents = document.createRange();
          contents.selectNodeContents(document.body);
          const bottom = contents.getBoundingClientRect().bottom + window.scrollY;
          const height = Math.ceil(bottom + parseFloat(style.paddingBottom) + parseFloat(style.marginBottom));
          frame.style.height = height + 'px';
        };
        resize = new ResizeObserver(resizeFrame);
        resize.observe(document.body);
        resize.observe(svg);
      }
    } catch (_) { /* The embedding controls sizing across origins. */ }

    const html = (tag, text) => {
      const element = document.createElementNS('http://www.w3.org/1999/xhtml', tag);
      if (text !== undefined) element.textContent = text;
      return element;
    };
    const svgElement = tag => document.createElementNS('http://www.w3.org/2000/svg', tag);
    /* Keep the drawing in its own camera. Inspector and controls use CSS-pixel
       coordinates, so zooming never scales their text or changes the notebook height. */
    const viewport = svgElement('svg');
    viewport.classList.add('linnet-viewport');
    viewport.setAttribute('tabindex', '0');
    viewport.setAttribute('role', 'group');
    viewport.setAttribute('aria-label', 'Graph. Drag to pan; plus or minus to zoom; zero to fit.');
    for (const child of [...svg.children]) {
      if (!['style', 'script', 'defs', 'title', 'desc'].includes(child.localName)) viewport.append(child);
    }
    const panSurface = svgElement('rect');
    panSurface.classList.add('linnet-pan-surface');
    svg.append(panSurface, viewport);
    const camera = { ...box };
    const toolbar = svgElement('foreignObject');
    const controls = html('div');
    controls.className = 'linnet-toolbar';
    controls.setAttribute('role', 'group');
    controls.setAttribute('aria-label', 'Graph navigation');
    const zoomLabel = html('output');
    zoomLabel.setAttribute('aria-label', 'Graph zoom');
    const paintCamera = () => {
      viewport.setAttribute('viewBox', `${camera.x} ${camera.y} ${camera.width} ${camera.height}`);
      zoomLabel.textContent = `${Math.round(100 * scale)}%`;
      schedulePopupPosition();
    };
    const zoom = (factor, anchor = { x: camera.x + camera.width / 2, y: camera.y + camera.height / 2 }) => {
      const next = Math.max(Math.min(.25, fitScale()), Math.min(20, scale * factor));
      const ratio = scale / next;
      scale = next;
      fitted = false;
      camera.x = anchor.x + (camera.x - anchor.x) * ratio;
      camera.y = anchor.y + (camera.y - anchor.y) * ratio;
      camera.width *= ratio;
      camera.height *= ratio;
      paintCamera();
    };
    const fitScale = () => Math.min(1,
      viewport.width.baseVal.value / (box.width * naturalScale),
      viewport.height.baseVal.value / (box.height * naturalScale));
    const reset = () => { fitted = true; layout(); };
    controls.append(zoomLabel);
    const help = html('span');
    help.className = 'linnet-navigation-help';
    help.tabIndex = 0;
    help.setAttribute('role', 'img');
    const shortcuts = 'Drag to pan · Ctrl/⌘ + scroll or +/− to zoom · 0 to fit · 100% = default size';
    help.setAttribute('aria-label', shortcuts);
    const hint = html('span', shortcuts);
    hint.className = 'linnet-navigation-hint';
    hint.setAttribute('aria-hidden', 'true');
    help.append(hint);
    controls.append(help);
    toolbar.append(controls);
    /* Use the accessible notebook document's top layer. Neither SVG bounds
       nor a notebook scroll pane should clip details or change graph geometry.
       Opaque sandbox frames keep their popup in their own document. */
    const popupDocument = frames.at(-1)?.ownerDocument || document;
    const panel = html('div');
    panel.className = 'linnet-inspector-popup';
    panel.setAttribute('popover', 'manual');
    panel.style.cssText = 'position:fixed;inset:auto;margin:0;padding:0;border:0;'
      + 'background:transparent;width:320px;max-width:calc(100vw - 16px);overflow:visible;';
    const popupRoot = panel.attachShadow({mode:'open'});
    const popupStyle = html('style');
    /* The SVG stylesheet remains the sole owner of inspector styling. Copy
       its rules into the portal's shadow root, without notebook CSS resets. */
    const prefix = 'svg[data-linnet-interactive] ';
    const inspectorRules = new Set();
    for (const sheet of document.styleSheets) {
      let rules;
      try { rules = sheet.cssRules; } catch (_) { continue; }
      for (const rule of rules) {
        if (rule.selectorText?.startsWith(prefix + '.linnet-inspector'))
          inspectorRules.add(rule.cssText.replaceAll(prefix, ''));
      }
    }
    popupStyle.textContent = [...inspectorRules].join('\n');
    const content = html('div');
    content.className = 'linnet-inspector';
    content.style.cssText = 'max-height:calc(100vh - 16px);overflow:auto;text-align:left;';
    popupRoot.append(popupStyle, content);
    svg.append(toolbar);
    let anchor = {x:.5, y:.5};
    let positionFrame;
    const positionPopup = () => {
      positionFrame = null;
      if (!panel.isConnected || !inspected?.isConnected) return;
      const target = inspected.getBoundingClientRect();
      let x = target.left + anchor.x * target.width;
      let y = target.top + anchor.y * target.height;
      for (const frame of frames) {
        const rect = frame.getBoundingClientRect();
        const sx = rect.width / frame.offsetWidth, sy = rect.height / frame.offsetHeight;
        x = rect.left + (frame.clientLeft + x) * sx;
        y = rect.top + (frame.clientTop + y) * sy;
      }
      const view = popupDocument.defaultView;
      const width = panel.offsetWidth, height = panel.offsetHeight;
      const right = view.innerWidth - 8, bottom = view.innerHeight - 8;
      panel.style.left = Math.max(8, Math.min(right - width,
        x + 12 + width <= right ? x + 12 : x - width - 12)) + 'px';
      panel.style.top = Math.max(8, Math.min(bottom - height,
        y + 12 + height <= bottom ? y + 12 : y - height - 12)) + 'px';
    };
    const schedulePopupPosition = () => {
      if (!panel.isConnected || positionFrame) return;
      positionFrame = requestAnimationFrame(positionPopup);
    };
    syncPopupTheme = () => {
      const style = getComputedStyle(svg);
      for (const name of ['--linnet-selection', '--linnet-panel-bg', '--linnet-panel-ink'])
        panel.style.setProperty(name, style.getPropertyValue(name));
    };
    const rememberAnchor = (target, event) => {
      const rect = target.getBoundingClientRect();
      anchor = event?.type.startsWith('pointer') || event?.type === 'click' && event.detail > 0
        ? {x:(event.clientX - rect.left) / (rect.width || 1), y:(event.clientY - rect.top) / (rect.height || 1)}
        : {x:.5, y:.5};
    };
    /* An SVG document has no HTML body. A zero-sized foreignObject supplies
       the HTML mounting context; the popover still paints in the top layer. */
    const popupMount = popupDocument.body || svg.appendChild(svgElement('foreignObject'));
    const close = event => {
      /* Removing a pinned popup can expose a hit target under its Close
         button. Resume previews on pointer motion, not that synthetic enter. */
      if (event?.type === 'click' || event?.type === 'keydown') hoverDismissed = true;
      if (panel.matches(':popover-open')) panel.hidePopover();
      panel.remove();
      content.replaceChildren();
      cancelAnimationFrame(positionFrame);
      positionFrame = null;
      pinned = false;
    };
    closePopup = close;
    panel.addEventListener('linnet-close', close);
    panel.addEventListener('keydown', event => {
      if (event.key !== 'Escape') return;
      event.preventDefault();
      inspected?.focus({preventScroll:true});
      close(event);
    });
    for (const view of new Set([window, ...frames.map(frame => frame.ownerDocument.defaultView)])) {
      view.addEventListener('scroll', schedulePopupPosition, {capture:true, passive:true, signal:lifecycle.signal});
      view.addEventListener('resize', schedulePopupPosition, {passive:true, signal:lifecycle.signal});
    }
    for (const frame of frames) frame.addEventListener('load', syncHost, {signal:lifecycle.signal});
    window.addEventListener('pagehide', close, {signal:lifecycle.signal});
    window.matchMedia('(prefers-color-scheme: dark)')
      .addEventListener('change', syncPopupTheme, {signal:lifecycle.signal});
    let layoutWidth = 0;
    const layout = () => {
      const screenWidth = Math.max(1, svg.getBoundingClientRect().width);
      const width = embedded ? box.width : screenWidth;
      layoutWidth = screenWidth;
      const fixedHeight = Number(svg.dataset.linnetViewportHeight) || 0;
      const graphWidth = width;
      /* Keep the graph and caption centered without reserving empty space for
         an inspector. The camera height follows the bounded initial drawing size. */
      const fitWidth = Math.min(1, graphWidth / (box.width * naturalScale));
      const graphHeight = embedded ? box.height : fixedHeight || box.height * naturalScale * fitWidth;
      /* Navigation floats over one corner; it never reserves a canvas row. */
      const controlWidth = Math.min(84, Math.max(1, width - 8));
      toolbar.setAttribute('x', Math.max(0, width - controlWidth - 4));
      toolbar.setAttribute('y', 4);
      toolbar.setAttribute('width', controlWidth);
      toolbar.setAttribute('height', 28);
      controls.style.setProperty('--linnet-hint-width', `${Math.max(1, width - 8)}px`);
      viewport.setAttribute('x', 0);
      viewport.setAttribute('y', 0);
      viewport.setAttribute('width', graphWidth);
      viewport.setAttribute('height', graphHeight);
      const center = fitted
        ? { x: box.x + box.width / 2, y: box.y + box.height / 2 }
        : { x: camera.x + camera.width / 2, y: camera.y + camera.height / 2 };
      if (fitted) scale = fitScale();
      camera.width = graphWidth / (naturalScale * scale);
      camera.height = graphHeight / (naturalScale * scale);
      camera.x = center.x - camera.width / 2;
      camera.y = center.y - camera.height / 2;
      paintCamera();
      panSurface.setAttribute('y', 0);
      panSurface.setAttribute('width', graphWidth);
      panSurface.setAttribute('height', graphHeight);
      const height = graphHeight;
      if (!embedded) {
        svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
        svg.setAttribute('height', height);
        svg.style.height = `${height}px`;
      }
      requestAnimationFrame(resizeFrame);
    };
    if (!embedded) {
      svg.style.width = '100%';
      svg.style.maxWidth = '100%';
    }
    const layoutObserver = new ResizeObserver(() => {
      if (!svg.isConnected) { layoutObserver.disconnect(); return; }
      if (Math.abs(svg.getBoundingClientRect().width - layoutWidth) > .5) requestAnimationFrame(layout);
    });
    layoutObserver.observe(svg);
    layout();
    const show = (target, event) => {
      inspected = target;
      rememberAnchor(target, event);
      content.replaceChildren();
      const kind = target.dataset.linnetKind;
      const id = Number(target.dataset.linnetId);
      const detail = JSON.parse(target.dataset.linnetDetail);
      const dismiss = html('button', '×');
      dismiss.setAttribute('aria-label', 'Close graph details');
      dismiss.hidden = !pinned;
      dismiss.addEventListener('click', close);
      const header = html('div');
      header.className = 'linnet-inspector-header';
      const identity = `${kind === 'node' ? 'Node' : kind === 'halfedge' ? 'Half-edge' : 'Edge'} ${id}`;
      const title = detail.title ? `${detail.title} · ${identity}` : `${identity}${detail.name ? ` · ${detail.name}` : ''}`;
      header.append(html('strong', title), dismiss);
      content.append(header);
      const details = html('div');
      details.className = 'linnet-inspector-details';
      if (kind !== 'node') {
        if (kind === 'halfedge') {
          details.append(html('span', `Edge ${detail.edge} · Node ${detail.node} · ${detail.flow}`));
          if (detail.pair !== null) details.append(html('span', `Paired half-edge: ${detail.pair}`));
        }
        if (detail.particle !== undefined) {
          const particle = html('span');
          particle.append(html('span', 'Particle: '), html('strong', String(detail.particle)));
          if (detail.pdg !== undefined) particle.append(` (${detail.pdg})`);
          details.append(particle);
        }
        if ('source' in detail || 'sink' in detail) {
          const source = detail.source == null ? 'External' : `Node ${detail.source}`;
          const sink = detail.sink == null ? 'External' : detail.source == null ? `Node ${detail.sink}` : String(detail.sink);
          const connector = detail.orientation === 'undirected' ? '—' : detail.orientation === 'reversed' ? '←' : '→';
          details.append(html('span', `${source} ${connector} ${sink}`));
        }
        if (detail['external-state']) {
          const state = String(detail['external-state']);
          details.append(html('span', state.charAt(0).toUpperCase() + state.slice(1)));
        }
      } else {
        if (detail.edges) details.append(html('span', `Edges: ${detail.edges.join(', ') || 'none'}`));
      }
      content.append(details);
      if (Array.isArray(detail.properties) && detail.properties.length) {
        const properties = html('dl');
        properties.className = 'linnet-inspector-properties';
        for (const [label, value] of detail.properties) {
          properties.append(html('dt', String(label)), html('dd', String(value)));
        }
        content.append(properties);
      }
      if (selected[kind].has(id)) {
        const current = selection();
        const argumentsText = Object.entries(current)
          .filter(([, entries]) => entries.length)
          .map(([name, entries]) => `${name}=[${entries.join(', ')}]`)
          .join(', ');
        const construction = html('div');
        construction.className = 'linnet-inspector-construction';
        construction.append(html('span', 'Subgraph construction:'), html('code', `graph.subgraph(${argumentsText})`));
        content.append(construction);
      }
      const shortcuts = html('div', `${pinned ? 'Pinned · Close to resume hover' : 'Preview · Click or Enter/Space to pin'} · Shift/Ctrl/⌘-click: toggle selection`);
      shortcuts.className = 'linnet-inspector-shortcuts';
      content.append(shortcuts);
      for (const previous of popupDocument.querySelectorAll('.linnet-inspector-popup')) {
        if (previous !== panel) previous.dispatchEvent(new CustomEvent('linnet-close'));
      }
      if (!panel.isConnected) popupMount.append(panel);
      syncPopupTheme();
      content.setAttribute('role', pinned ? 'dialog' : 'tooltip');
      panel.style.pointerEvents = pinned ? 'auto' : 'none';
      if (!panel.matches(':popover-open')) panel.showPopover();
      positionPopup();
      svg.dispatchEvent(new CustomEvent('linnet-inspect', { bubbles: true, detail: { kind, id, data: detail } }));
    };
    const activate = (target, event) => {
      event.preventDefault();
      if (event.shiftKey || event.ctrlKey || event.metaKey) {
        const kind = target.dataset.linnetKind;
        const id = Number(target.dataset.linnetId);
        if (selected[kind].has(id)) selected[kind].delete(id);
        else selected[kind].add(id);
        for (const item of targets) item.setAttribute('aria-pressed', String(selected[item.dataset.linnetKind].has(Number(item.dataset.linnetId))));
        svg.dispatchEvent(new CustomEvent('linnet-selection-change', { bubbles: true, detail: selection() }));
      }
      pinned = true;
      show(target, event);
    };
    const targetOf = event => event.target.closest && event.target.closest('[data-linnet-kind]');
    /* Domain overlays and hover share the native structural half-edge paths.
       One opacity layer composites shared vertices and adjoining halves once. */
    let shadeSerial = 0;
    svg.linnetShade = ({ nodes = [], edges = [], half_edges = [] },
      { width = 8, padding = 2, outline = 0, fillOpacity = .12 } = {}) => {
      /* Read current geometry: domain viewers can replace the drawing while
         keeping this camera, and several regions may share the same carriers. */
      const carriers = [...svg.querySelectorAll('[data-linnet-carrier]')];
      const nodeShapes = new Map([...svg.querySelectorAll('[data-linnet-node]')]
        .map(node => [Number(node.dataset.linnetNode), node]));
      const included = { node: new Set(nodes), edge: new Set(edges), halfedge: new Set(half_edges) };
      const group = svgElement('g');
      group.classList.add('linnet-edge-highlight');
      if (outline > 0) {
        /* Outline the union, not each primitive: no seams where halves meet
           vertices. A faint interior and opaque rim stay distinct when nested. */
        const filter = svgElement('filter');
        filter.id = `linnet-region-${Math.random().toString(36).slice(2)}-${shadeSerial++}`;
        for (const [key, value] of Object.entries({x: '-50%', y: '-50%', width: '200%', height: '200%'}))
          filter.setAttribute(key, value);
        /* Construct nodes, not markup inside this embedded SVG script. */
        const primitive = (tag, attrs, parent = filter) => {
          const element = svgElement(tag);
          for (const [key, value] of Object.entries(attrs)) element.setAttribute(key, value);
          parent.append(element);
          return element;
        };
        primitive('feMorphology', {in: 'SourceAlpha', operator: 'dilate', radius: outline, result: 'outer'});
        primitive('feComposite', {in: 'outer', in2: 'SourceAlpha', operator: 'out', result: 'rim'});
        primitive('feFlood', {'flood-color': 'currentColor'});
        primitive('feComposite', {in2: 'rim', operator: 'in', result: 'border'});
        const tint = primitive('feComponentTransfer', {in: 'SourceGraphic', result: 'tint'});
        primitive('feFuncA', {type: 'linear', slope: fillOpacity}, tint);
        const merge = primitive('feMerge', {});
        for (const input of ['tint', 'border']) primitive('feMergeNode', {in: input}, merge);
        const defs = svgElement('defs');
        defs.append(filter);
        group.append(defs);
        group.setAttribute('filter', `url(#${filter.id})`);
        group.style.opacity = '1';
      }
      /* Measure a child of the camera: Firefox's matrix for a nested <svg>
         omits its own viewBox transform. The empty layer has no painted ink. */
      viewport.append(group);
      const inverse = group.getScreenCTM().inverse();
      group.remove();
      const place = (copy, original) => {
        const transform = inverse.multiply(original.getScreenCTM());
        copy.setAttribute('transform', `matrix(${transform.a} ${transform.b} ${transform.c} ${transform.d} ${transform.e} ${transform.f})`);
        group.append(copy);
      };
      /* Authored Typst drawings expose sampled hit regions rather than native
         carriers. Preserve their painted geometry, including custom glyphs. */
      if (svg.dataset.linnetRenderer !== 'native') {
        for (const item of targets) {
          const detail = JSON.parse(item.dataset.linnetDetail);
          const node = item.dataset.linnetKind === 'node';
          if (node ? !included.node.has(Number(item.dataset.linnetId))
            : !included.edge.has(detail.edge) && !included.halfedge.has(detail['half-edge'])) continue;
          for (const rect of item.querySelectorAll('rect')) {
            const copy = rect.cloneNode(false);
            copy.setAttribute('fill', 'currentColor');
            copy.setAttribute('rx', Math.min(rect.width.baseVal.value, rect.height.baseVal.value) / 2);
            place(copy, rect);
          }
        }
        return group;
      }
      for (const carrier of carriers) {
        const detail = JSON.parse(carrier.dataset.linnetDetail);
        if (!included.edge.has(detail.edge) && !included.halfedge.has(detail['half-edge'])) continue;
        const original = carrier.querySelector('path');
        const path = svgElement('path');
        let d = original.getAttribute('d');
        const node = nodeShapes.get(detail.node);
        if (node && included.node.has(detail.node)) {
          /* Join selected vertices to their incident halves even when labels
             reserved a larger vertex disc in the original drawing. */
          const center = new DOMPoint(node.cx.baseVal.value, node.cy.baseVal.value)
            .matrixTransform(original.getScreenCTM().inverse().multiply(node.getScreenCTM()));
          d = detail.flow === 'source' ? `M${center.x} ${center.y} L${d.slice(1)}` : `${d} L${center.x} ${center.y}`;
        }
        path.setAttribute('d', d);
        path.setAttribute('fill', 'none');
        path.setAttribute('stroke', 'currentColor');
        path.setAttribute('stroke-width', width);
        /* Half-edge membership ends exactly at the cut, including narrow gaps. */
        path.setAttribute('stroke-linecap', 'butt');
        path.setAttribute('stroke-linejoin', 'round');
        path.dataset.linnetShadeHalfedge = detail['half-edge'];
        place(path, original);
      }
      for (const id of included.node) {
        const node = nodeShapes.get(id);
        if (!node) continue;
        const copy = node.cloneNode(false);
        copy.removeAttribute('data-linnet-node');
        copy.setAttribute('r', node.r.baseVal.value + padding);
        copy.setAttribute('fill', 'currentColor');
        copy.setAttribute('stroke', 'none');
        copy.dataset.linnetShadeNode = id;
        place(copy, node);
      }
      return group;
    };
    const highlightEdge = (target, active) => {
      const kind = target.dataset.linnetKind;
      if (kind === 'node') return;
      const id = Number(target.dataset.linnetId);
      const key = `${kind}:${id}`;
      let group = edgeHighlights.get(key);
      if (!group) {
        group = svg.linnetShade(kind === 'edge' ? { edges: [id] } : { half_edges: [id] });
        viewport.append(group);
        edgeHighlights.set(key, group);
      }
      group.toggleAttribute('data-selected', selected[kind].has(id));
      group.style.display = active || selected[kind].has(id) ? '' : 'none';
    };
    const graphPoint = event => new DOMPoint(event.clientX, event.clientY).matrixTransform(viewport.getScreenCTM().inverse());
    const inGraph = event => event.target === panSurface || viewport.contains(event.target);
    svg.addEventListener('wheel', event => {
      if (!inGraph(event)) return;
      if (!event.ctrlKey && !event.metaKey) return;
      event.preventDefault();
      const unit = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? viewport.height.baseVal.value : 1;
      const delta = Math.max(-100, Math.min(100, event.deltaY * unit));
      zoom(Math.exp(-delta * .001), graphPoint(event));
    }, { passive: false });
    let pointer;
    svg.addEventListener('pointerdown', event => {
      if (event.button !== 0 || !inGraph(event)) return;
      if (event.target === panSurface) viewport.focus();
      pointer = { id: event.pointerId, x: event.clientX, y: event.clientY, camera: { ...camera }, matrix: viewport.getScreenCTM().inverse() };
      suppressClick = false;
    });
    svg.addEventListener('pointermove', event => {
      if (hoverDismissed) {
        hoverDismissed = false;
        const target = targetOf(event);
        if (!pinned && target) show(target, event);
      }
      if (!pinned && panel.isConnected && inspected) {
        rememberAnchor(inspected, event);
        schedulePopupPosition();
      }
      if (!pointer || pointer.id !== event.pointerId) return;
      if (!(event.buttons & 1)) { endPan(event); return; }
      const dx = event.clientX - pointer.x, dy = event.clientY - pointer.y;
      if (!dragging && Math.hypot(dx, dy) < 4) return;
      dragging = true;
      fitted = false;
      viewport.setPointerCapture(event.pointerId);
      viewport.classList.add('linnet-panning');
      camera.x = pointer.camera.x - dx * pointer.matrix.a - dy * pointer.matrix.c;
      camera.y = pointer.camera.y - dx * pointer.matrix.b - dy * pointer.matrix.d;
      paintCamera();
    });
    const endPan = event => {
      if (!pointer || pointer.id !== event.pointerId) return;
      suppressClick = dragging;
      dragging = false;
      pointer = null;
      viewport.classList.remove('linnet-panning');
      if (viewport.hasPointerCapture(event.pointerId)) viewport.releasePointerCapture(event.pointerId);
    };
    svg.addEventListener('pointerup', endPan);
    svg.addEventListener('pointercancel', endPan);
    viewport.addEventListener('lostpointercapture', endPan);
    viewport.addEventListener('keydown', event => {
      if (event.key === '+' || event.key === '=') zoom(1.05);
      else if (event.key === '-') zoom(1 / 1.05);
      else if (event.key === '0') reset();
      else if (['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown'].includes(event.key)) {
        fitted = false;
        camera.x += (event.key === 'ArrowLeft' ? -1 : event.key === 'ArrowRight' ? 1 : 0) * camera.width / 10;
        camera.y += (event.key === 'ArrowUp' ? -1 : event.key === 'ArrowDown' ? 1 : 0) * camera.height / 10;
        paintCamera();
      } else return;
      event.preventDefault();
    });
    /* Domain controls also live inside the viewport. Cancel a drag's click
       before it reaches their handlers, not just the built-in inspector. */
    svg.addEventListener('click', event => {
      if (suppressClick) { suppressClick = false; event.stopImmediatePropagation(); }
    }, true);
    svg.addEventListener('click', event => {
      const target = targetOf(event);
      if (target) { activate(target, event); highlightEdge(target, target.hasAttribute('data-linnet-active')); }
    });
    svg.addEventListener('keydown', event => {
      if (event.key === 'Escape') { if (inspected) inspected.focus(); close(event); }
      const target = targetOf(event);
      if (target && (event.key === 'Enter' || event.key === ' ')) {
        activate(target, event);
        highlightEdge(target, target.hasAttribute('data-linnet-active'));
      }
    });
    for (const [name, active] of [['pointerover', true], ['pointerout', false], ['focusin', true], ['focusout', false]]) {
      svg.addEventListener(name, event => {
        const target = targetOf(event);
        if (!target || dragging) return;
        const related = event.relatedTarget?.closest?.('[data-linnet-kind]');
        if (related && svg.contains(related) && related.dataset.linnetKind === target.dataset.linnetKind && related.dataset.linnetId === target.dataset.linnetId) return;
        if (!pinned) {
          if (name === 'focusin') hoverDismissed = false;
          if (active && !hoverDismissed) show(target, event);
          else if (!active) close();
        }
        highlightEdge(target, active);
        for (const item of targets) {
          if (item.dataset.linnetKind === target.dataset.linnetKind && item.dataset.linnetId === target.dataset.linnetId) {
            item.toggleAttribute('data-linnet-active', active);
          }
        }
      });
    }
  }
})();
