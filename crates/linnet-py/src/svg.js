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
    const original = Object.fromEntries(['viewBox', 'width', 'height'].map(key => [key, svg.getAttribute(key)]));
    const bounds = svg.viewBox.baseVal;
    const box = { x: bounds.x, y: bounds.y, width: bounds.width, height: bounds.height };
    let panel;
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

    const close = () => {
      if (panel) panel.remove();
      panel = null;
      for (const [key, value] of Object.entries(original)) {
        if (value === null) svg.removeAttribute(key);
        else svg.setAttribute(key, value);
      }
      requestAnimationFrame(resizeFrame);
    };
    const html = (tag, text) => {
      const element = document.createElementNS('http://www.w3.org/1999/xhtml', tag);
      if (text !== undefined) element.textContent = text;
      return element;
    };
    const show = target => {
      inspected = target;
      close();
      const kind = target.dataset.linnetKind;
      const id = Number(target.dataset.linnetId);
      const detail = JSON.parse(target.dataset.linnetDetail);
      const width = Math.max(box.width, 310);
      const height = 1;
      panel = document.createElementNS('http://www.w3.org/2000/svg', 'foreignObject');
      panel.setAttribute('x', box.x);
      panel.setAttribute('y', box.y + box.height + 8);
      panel.setAttribute('width', width);
      panel.setAttribute('height', height);
      const content = html('div');
      content.className = 'linnet-inspector';
      content.setAttribute('role', 'status');
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
      const shortcuts = html('div', 'Click or Enter/Space: details · Shift/Ctrl/⌘-click: toggle selection');
      shortcuts.className = 'linnet-inspector-shortcuts';
      content.append(shortcuts);
      panel.append(content);
      svg.append(panel);
      svg.setAttribute('viewBox', `${box.x} ${box.y} ${width} ${box.height + height + 8}`);
      svg.setAttribute('width', width + 'pt');
      svg.setAttribute('height', (box.height + height + 8) + 'pt');
      /* Fit the complete panel after layout, including wrapped construction
         expressions. Neither long selections nor narrow displays need a scrollbar. */
      const fittedHeight = content.offsetHeight;
      panel.setAttribute('height', fittedHeight);
      svg.setAttribute('viewBox', `${box.x} ${box.y} ${width} ${box.height + fittedHeight + 8}`);
      svg.setAttribute('height', (box.height + fittedHeight + 8) + 'pt');
      requestAnimationFrame(resizeFrame);
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
        const inverse = svg.getCTM().inverse();
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
        svg.append(group);
        edgeHighlights.set(key, group);
      }
      group.toggleAttribute('data-selected', selected[kind].has(id));
      group.style.display = active || selected[kind].has(id) ? '' : 'none';
    };
    svg.addEventListener('click', event => {
      const target = targetOf(event);
      if (target) { activate(target, event); highlightEdge(target, false); }
    });
    svg.addEventListener('keydown', event => {
      if (event.key === 'Escape') { close(); if (inspected) inspected.focus(); }
      const target = targetOf(event);
      if (target && (event.key === 'Enter' || event.key === ' ')) {
        activate(target, event);
        highlightEdge(target, false);
      }
    });
    for (const [name, active] of [['pointerover', true], ['pointerout', false], ['focusin', true], ['focusout', false]]) {
      svg.addEventListener(name, event => {
        const target = targetOf(event);
        if (!target) return;
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
