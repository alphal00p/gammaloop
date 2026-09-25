(() => {
  /* Scope state to the SVG: identical graph IDs can occur in several outputs.
     Use block comments because notebook iframe wrappers may flatten newlines. */
  for (const svg of document.querySelectorAll('svg[data-linnet-interactive]')) {
    if (svg.linnetSelection !== undefined) continue;
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
    let pinned = false;
    let dragging = false;
    let suppressClick = false;
    let inspected;
    let resize;
    let resizeFrame = () => {};
    const frames = [];

    const theme = () => {
      /* Marimo executes scripted HTML in an iframe. Read its notebook theme
         when accessible; standalone SVGs retain the system color preference. */
      let element = svg.parentElement;
      while (element) {
        const explicit = element.getAttribute('data-theme');
        const jupyter = element.getAttribute('data-jp-theme-light');
        if (explicit === 'dark' || explicit === 'light') return explicit;
        if (jupyter !== null) return jupyter === 'true' ? 'light' : 'dark';
        if (element.classList.contains('dark')) return 'dark';
        if (element.classList.contains('light')) return 'light';
        if (element.parentElement) element = element.parentElement;
        else {
          try { element = element.ownerDocument.defaultView.frameElement; }
          catch (_) { element = null; }
        }
      }
      return null;
    };
    const syncTheme = () => {
      if (!svg.isConnected || frames.some(frame => !frame.isConnected)) {
        observer.disconnect();
        if (resize) resize.disconnect();
        return;
      }
      const value = theme();
      if (value) svg.setAttribute('data-theme', value);
      else svg.removeAttribute('data-theme');
      try {
        const frame = window.frameElement;
        if (frame && frame.parentElement && document.body) {
          const style = frame.ownerDocument.defaultView.getComputedStyle(frame.parentElement);
          document.body.style.color = style.color;
          document.body.style.fontFamily = style.fontFamily;
          /* Browsers paint an opaque canvas when an iframe's declared scheme
             differs from its document, even when both backgrounds are transparent. */
          document.documentElement.style.colorScheme = frame.ownerDocument.defaultView.getComputedStyle(frame).colorScheme;
        }
      } catch (_) { /* Cross-origin embeddings retain their own text styling. */ }
    };
    const observer = new MutationObserver(syncTheme);
    syncTheme();
    let ancestor = svg.parentElement;
    while (ancestor) {
      if (ancestor.localName === 'iframe') frames.push(ancestor);
      observer.observe(ancestor, { attributes: true, childList: true, attributeFilter: ['class', 'data-theme', 'data-jp-theme-light'] });
      if (ancestor.parentElement) ancestor = ancestor.parentElement;
      else {
        try { ancestor = ancestor.ownerDocument.defaultView.frameElement; }
        catch (_) { ancestor = null; }
      }
    }
    try {
      const frame = window.frameElement;
      if (frame && document.body) {
        /* Measure content rather than the viewport: Marimo's own iframe
           observer grows outputs but deliberately does not shrink them. */
        document.documentElement.style.overflow = 'hidden';
        resizeFrame = () => {
          if (!svg.isConnected) { resize.disconnect(); return; }
          const style = getComputedStyle(document.body);
          const contents = document.createRange();
          contents.selectNodeContents(document.body);
          const bottom = contents.getBoundingClientRect().bottom + window.scrollY;
          const height = Math.ceil(bottom + parseFloat(style.paddingBottom) + parseFloat(style.marginBottom) + 6);
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
       coordinates, so zooming never scales text or changes the notebook height. */
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
      zoomLabel.textContent = `${Math.round(100 * box.width / camera.width)}%`;
    };
    const zoom = (factor, anchor = { x: camera.x + camera.width / 2, y: camera.y + camera.height / 2 }) => {
      const scale = Math.max(.25, Math.min(20, box.width / camera.width * factor));
      const ratio = box.width / scale / camera.width;
      camera.x = anchor.x + (camera.x - anchor.x) * ratio;
      camera.y = anchor.y + (camera.y - anchor.y) * ratio;
      camera.width *= ratio;
      camera.height *= ratio;
      paintCamera();
    };
    const reset = () => { Object.assign(camera, box); paintCamera(); };
    controls.append(zoomLabel);
    const hint = html('span', 'Drag to pan · Ctrl/⌘ + scroll or +/− to zoom · 0 to fit');
    hint.className = 'linnet-navigation-hint';
    controls.append(hint);
    toolbar.append(controls);
    svg.append(toolbar);
    const panel = svgElement('foreignObject');
    const content = html('div');
    content.className = 'linnet-inspector';
    content.setAttribute('role', 'status');
    panel.append(content);
    svg.append(panel);
    let layoutWidth = 0;
    const layout = () => {
      const width = Math.max(1, svg.getBoundingClientRect().width);
      layoutWidth = width;
      const wide = width >= 600;
      const panelWidth = wide ? 260 : width;
      const graphWidth = wide ? width - panelWidth - 16 : width;
      const graphHeight = Math.max(240, Math.min(520, graphWidth * box.height / box.width));
      const toolbarHeight = 40;
      toolbar.setAttribute('width', width);
      toolbar.setAttribute('height', toolbarHeight);
      viewport.setAttribute('x', 0);
      viewport.setAttribute('y', toolbarHeight);
      viewport.setAttribute('width', graphWidth);
      viewport.setAttribute('height', graphHeight);
      panSurface.setAttribute('y', toolbarHeight);
      panSurface.setAttribute('width', graphWidth);
      panSurface.setAttribute('height', graphHeight);
      panel.setAttribute('x', wide ? graphWidth + 16 : 0);
      panel.setAttribute('y', wide ? toolbarHeight : toolbarHeight + graphHeight + 12);
      panel.setAttribute('width', panelWidth);
      const panelHeight = panel.style.display === 'none' ? 0 : content.offsetHeight;
      panel.setAttribute('height', panelHeight);
      const height = toolbarHeight + (wide ? Math.max(graphHeight, panelHeight) : graphHeight + (panelHeight ? panelHeight + 12 : 0));
      svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
      svg.setAttribute('height', height);
      svg.style.height = `${height}px`;
      requestAnimationFrame(resizeFrame);
    };
    const close = () => {
      pinned = false;
      content.replaceChildren();
      panel.style.display = 'none';
      layout();
    };
    svg.style.width = '100%';
    svg.style.maxWidth = '100%';
    const layoutObserver = new ResizeObserver(() => {
      if (!svg.isConnected) { layoutObserver.disconnect(); return; }
      if (Math.abs(svg.getBoundingClientRect().width - layoutWidth) > .5) requestAnimationFrame(layout);
    });
    layoutObserver.observe(svg);
    reset();
    close();
    const show = target => {
      inspected = target;
      content.replaceChildren();
      panel.style.display = '';
      const kind = target.dataset.linnetKind;
      const id = Number(target.dataset.linnetId);
      const detail = JSON.parse(target.dataset.linnetDetail);
      const dismiss = html('button', '×');
      dismiss.setAttribute('aria-label', 'Close graph details');
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
        const construction = html('div');
        construction.className = 'linnet-inspector-construction';
        construction.append(html('span', 'Subgraph construction:'), html('code', `graph.subgraph(nodes=[${current.nodes.join(', ')}], edges=[${current.edges.join(', ')}], half_edges=[${current.half_edges.join(', ')}])`));
        content.append(construction);
      }
      const shortcuts = html('div', `${pinned ? 'Pinned · Close to resume hover' : 'Preview · Click or Enter/Space to pin'} · Shift/Ctrl/⌘-click: toggle selection`);
      shortcuts.className = 'linnet-inspector-shortcuts';
      content.append(shortcuts);
      layout();
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
      show(target);
    };
    const targetOf = event => event.target.closest && event.target.closest('[data-linnet-kind]');
    const highlightEdge = (target, active) => {
      const kind = target.dataset.linnetKind;
      if (kind === 'node') return;
      const id = Number(target.dataset.linnetId);
      const key = `${kind}:${id}`;
      let group = edgeHighlights.get(key);
      if (!group) {
        /* Composite the overlapping hit regions once, so selection adds a
           smooth translucent halo instead of darkening each sampled box. */
        group = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        group.classList.add('linnet-edge-highlight');
        const inverse = viewport.getCTM().inverse();
        for (const item of targets) {
          if (item.dataset.linnetKind === 'node') continue;
          const detail = JSON.parse(item.dataset.linnetDetail);
          if ((kind === 'edge' ? detail.edge : detail['half-edge']) !== id) continue;
          for (const rect of item.querySelectorAll('rect')) {
            const copy = rect.cloneNode(false);
            const transform = inverse.multiply(rect.getCTM());
            copy.setAttribute('transform', `matrix(${transform.a} ${transform.b} ${transform.c} ${transform.d} ${transform.e} ${transform.f})`);
            copy.removeAttribute('fill');
            copy.setAttribute('rx', Math.min(rect.width.baseVal.value, rect.height.baseVal.value) / 2);
            group.append(copy);
          }
        }
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
      if (!pointer || pointer.id !== event.pointerId) return;
      if (!(event.buttons & 1)) { endPan(event); return; }
      const dx = event.clientX - pointer.x, dy = event.clientY - pointer.y;
      if (!dragging && Math.hypot(dx, dy) < 4) return;
      dragging = true;
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
        camera.x += (event.key === 'ArrowLeft' ? -1 : event.key === 'ArrowRight' ? 1 : 0) * camera.width / 10;
        camera.y += (event.key === 'ArrowUp' ? -1 : event.key === 'ArrowDown' ? 1 : 0) * camera.height / 10;
        paintCamera();
      } else return;
      event.preventDefault();
    });
    svg.addEventListener('click', event => {
      if (suppressClick) { suppressClick = false; return; }
      const target = targetOf(event);
      if (target) { activate(target, event); highlightEdge(target, false); }
    });
    svg.addEventListener('keydown', event => {
      if (event.key === 'Escape') { if (inspected) inspected.focus(); close(); }
      const target = targetOf(event);
      if (target && (event.key === 'Enter' || event.key === ' ')) {
        activate(target, event);
        highlightEdge(target, false);
      }
    });
    for (const [name, active] of [['pointerover', true], ['pointerout', false], ['focusin', true], ['focusout', false]]) {
      svg.addEventListener(name, event => {
        const target = targetOf(event);
        if (!target || dragging) return;
        const related = event.relatedTarget?.closest?.('[data-linnet-kind]');
        if (related && svg.contains(related) && related.dataset.linnetKind === target.dataset.linnetKind && related.dataset.linnetId === target.dataset.linnetId) return;
        if (!pinned) {
          if (active) show(target);
          else close();
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
