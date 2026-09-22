(() => {
  // Scope state to the SVG: identical graph IDs can occur in several outputs.
  for (const svg of document.querySelectorAll('svg[data-linnet-interactive]')) {
    if (svg.linnetSelection !== undefined) continue;
    const targets = [...svg.querySelectorAll('[data-linnet-kind]')];
    const selected = { node: new Set(), edge: new Set() };
    const selection = () => ({
      nodes: [...selected.node].sort((a, b) => a - b),
      edges: [...selected.edge].sort((a, b) => a - b),
    });
    Object.defineProperty(svg, 'linnetSelection', { get: selection });
    const original = Object.fromEntries(['viewBox', 'width', 'height'].map(key => [key, svg.getAttribute(key)]));
    const bounds = svg.viewBox.baseVal;
    const box = { x: bounds.x, y: bounds.y, width: bounds.width, height: bounds.height };
    let panel;
    let inspected;
    let resize;
    const frames = [];

    const theme = () => {
      // Marimo executes scripted HTML in an iframe. Read its notebook theme
      // when accessible; standalone SVGs retain the system color preference.
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
        resize = new ResizeObserver(() => {
          if (!svg.isConnected) { resize.disconnect(); return; }
          const style = getComputedStyle(document.body);
          frame.style.height = Math.ceil(document.body.getBoundingClientRect().height + parseFloat(style.marginTop) + parseFloat(style.marginBottom)) + 'px';
        });
        resize.observe(document.body);
      }
    } catch (_) { /* The embedding controls sizing across origins. */ }

    const close = () => {
      if (panel) panel.remove();
      panel = null;
      for (const [key, value] of Object.entries(original)) {
        if (value === null) svg.removeAttribute(key);
        else svg.setAttribute(key, value);
      }
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
      const height = 116;
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
      content.append(dismiss, html('strong', `${kind === 'node' ? 'Node' : 'Edge'} ${id}`));
      content.append(html('div', Object.entries(detail).map(([key, value]) => `${key}: ${typeof value === 'object' ? JSON.stringify(value) : value}`).join(' · ')));
      const current = selection();
      content.append(html('code', `subgraph(nodes=[${current.nodes.join(', ')}], edges=[${current.edges.join(', ')}])`));
      content.append(html('div', 'Shift/Ctrl/⌘-click to toggle selection. Esc to close.'));
      panel.append(content);
      svg.append(panel);
      svg.setAttribute('viewBox', `${box.x} ${box.y} ${width} ${box.height + height + 8}`);
      svg.setAttribute('width', width + 'pt');
      svg.setAttribute('height', (box.height + height + 8) + 'pt');
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
    svg.addEventListener('click', event => {
      const target = targetOf(event);
      if (target) activate(target, event);
    });
    svg.addEventListener('keydown', event => {
      if (event.key === 'Escape') { close(); if (inspected) inspected.focus(); }
      const target = targetOf(event);
      if (target && (event.key === 'Enter' || event.key === ' ')) activate(target, event);
    });
    for (const [name, active] of [['pointerover', true], ['pointerout', false]]) {
      svg.addEventListener(name, event => {
        const target = targetOf(event);
        if (!target) return;
        for (const item of targets) {
          if (item.dataset.linnetKind === target.dataset.linnetKind && item.dataset.linnetId === target.dataset.linnetId) {
            item.toggleAttribute('data-linnet-active', active);
          }
        }
      });
    }
  }
})();
