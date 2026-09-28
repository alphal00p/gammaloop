(() => {
  for (const root of document.querySelectorAll('.feynkit-collection')) {
    if (root.dataset.ready) continue;
    root.dataset.ready = 'true';
    const entries = [...root.querySelectorAll('template')];
    const mount = (index, stage) => {
      const fragment = entries[index].content.cloneNode(true);
      const svg = fragment.querySelector('svg');
      svg.dataset.linnetViewportHeight = '360';
      const sources = [...fragment.querySelectorAll('script')];
      sources.forEach(source => source.remove());
      stage.replaceChildren(fragment);
      /* The templates contain the native Linnet renderer. Activate its own
         script after insertion; graph navigation and inspection stay shared. */
      for (const source of sources) {
        const script = document.createElement('script');
        script.textContent = source.textContent;
        stage.append(script);
      }
    };
    const thumbnails = entries.map(entry => {
      const thumbnail = entry.content.querySelector('svg').cloneNode(true);
      thumbnail.querySelectorAll('script').forEach(script => script.remove());
      return thumbnail;
    });
    /* A collapsed amplitude has no live graph to bridge the notebook theme.
       Observe accessible frame ancestors so its static thumbnails match too. */
    const ancestors = [];
    let ancestor = root.parentElement;
    while (ancestor) {
      ancestors.push(ancestor);
      try { ancestor = ancestor.parentElement || ancestor.ownerDocument.defaultView.frameElement; }
      catch (_) { ancestor = null; }
    }
    let previousTheme;
    const syncTheme = () => {
      if (!root.isConnected) { observer.disconnect(); return; }
      let theme = '';
      for (const element of ancestors) {
        const explicit = element.dataset.theme;
        const jupyter = element.getAttribute('data-jp-theme-light');
        if (explicit === 'dark' || explicit === 'light') { theme = explicit; break; }
        if (jupyter !== null) { theme = jupyter === 'true' ? 'light' : 'dark'; break; }
        if (element.classList.contains('dark')) { theme = 'dark'; break; }
        if (element.classList.contains('light')) { theme = 'light'; break; }
      }
      try {
        const frame = window.frameElement;
        if (frame?.parentElement) {
          const style = frame.ownerDocument.defaultView.getComputedStyle(frame.parentElement);
          document.body.style.color = style.color;
          document.body.style.fontFamily = style.fontFamily;
          document.documentElement.style.colorScheme = style.colorScheme;
        }
      } catch (_) { /* Cross-origin hosts supply their own frame theme. */ }
      if (theme === previousTheme) return;
      previousTheme = theme;
      root.dataset.theme = theme;
      root.style.colorScheme = theme || 'normal';
      const images = root.querySelectorAll('.fk-thumbnail');
      thumbnails.forEach((thumbnail, index) => {
        thumbnail.setAttribute('data-theme', theme);
        images[index].src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(new XMLSerializer().serializeToString(thumbnail));
      });
    };
    const observer = new MutationObserver(syncTheme);
    ancestors.forEach(element => observer.observe(element, {attributes:true, childList:true, attributeFilter:['class', 'data-theme', 'data-jp-theme-light']}));
    syncTheme();
    const rows = root.querySelectorAll('.fk-row');
    for (const [index, row] of [...rows].entries()) {
      row.addEventListener('toggle', () => {
        const stage = row.querySelector('.fk-stage');
        if (row.open && !stage.firstChild) mount(index, stage);
      });
    }
    const buttons = [...root.querySelectorAll('.fk-strip button')];
    for (const [index, button] of buttons.entries()) {
      button.addEventListener('click', () => {
        buttons.forEach(other => other.setAttribute('aria-pressed', String(other === button)));
        mount(index, root.querySelector('.fk-stage'));
        root.querySelector('.fk-caption').textContent = button.dataset.caption;
      });
    }
    if (buttons.length) buttons[0].click();
  }
})();
