(() => {
  for (const root of document.querySelectorAll('.feynkit-collection')) {
    if (root.dataset.ready) continue;
    root.dataset.ready = 'true';
    const entries = [...root.querySelectorAll('template')];
    const mount = (index, stage) => {
      const fragment = entries[index].content.cloneNode(true);
      for (const svg of fragment.querySelectorAll('svg')) {
        svg.dataset.linnetViewportHeight = '360';
      }
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
    root.addEventListener('feynkit-theme', () => {
      const images = root.querySelectorAll('.fk-thumbnail');
      thumbnails.forEach((thumbnail, index) => {
        thumbnail.setAttribute('data-theme', root.dataset.theme);
        images[index].src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(new XMLSerializer().serializeToString(thumbnail));
      });
    });
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
