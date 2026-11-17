// Static source frames never install an editor or execution control.
export function frameSource(pre, {compact = false} = {}) {
  if (!compact && pre.closest('.cm-tooltip, .documentation-type-info')) return;
  let row = pre.closest('.reference-copy-row, .example-wrap, .codehilite, .code-frame');
  if (!row) {
    row = document.createElement('div');
    pre.before(row);
    row.append(pre);
  }
  if (row.classList.contains('code-frame')) return row;
  row.classList.add('code-frame');
  if (compact) { row.classList.add('documentation-example'); return row; }
  const gutter = document.createElement('span');
  gutter.className = 'code-line-numbers';
  gutter.setAttribute('aria-hidden', 'true');
  // Marimo copies the container's textContent. Generated gutter text keeps
  // visual line numbers out of its native clipboard value and text selection.
  gutter.dataset.lineNumbers = pre.querySelector('code').textContent.replace(/\n$/, '').split('\n').map((_, index) => index + 1).join('\n');
  row.prepend(gutter);
  return row;
}
export function initializeRustdoc() {
  import(new URL('navigation.js', import.meta.url)).then(module => module.initializeNavigation());
  const header = document.querySelector('.alphal00p-rustdoc-bar');
  const search = header?.querySelector('[data-search-open]');
  if (search) {
    search.setAttribute('aria-label', 'Search Rust API');
    search.title = 'Search Rust API (S or /)';
    search.addEventListener('click', () => window.searchState.focus());
  }
  header?.querySelector('[data-theme-toggle]')?.addEventListener('click', () => {
    const theme = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark';
    // Native Settings may later choose independently; nothing writes on page load.
    window.switchTheme(theme, true);
    try {
      localStorage.setItem('rustdoc-use-system-theme', 'false');
      localStorage.setItem('alphal00p-docs-theme', theme);
      localStorage.setItem('alphal00p-rustdoc-theme-synced', theme);
    } catch { /* Native controls still work without persistent storage. */ }
  });
  // Native initial source highlighting can run before the source font settles.
  // Re-anchor the same requested range once; keep subsequent native selections intact.
  const sourceRange = document.body.classList.contains('src') && /^#\d+(?:-\d+)?$/.test(location.hash) ? location.hash : null;
  if (sourceRange) document.fonts.ready.then(() => requestAnimationFrame(() => {
    if (location.hash === sourceRange) document.getElementById(sourceRange.slice(1).split('-')[0])?.scrollIntoView({block: 'start'});
  }));
  document.querySelectorAll('.docblock .example-wrap pre:has(> code)').forEach(frameSource);
}

// Markdown outputs are inserted reactively by both native Marimo frontends.
export function initializeNotebookSources(scope = document.body) {
  initializePopups(scope);
  const frame = () => scope.querySelectorAll('.codehilite pre:has(> code)').forEach(frameSource);
  frame();
  new MutationObserver(records => {
    if (records.some(record => [...record.addedNodes].some(node => node.nodeType === 1 &&
      (node.matches('pre, .codehilite') || node.querySelector('pre'))))) frame();
  }).observe(scope, {childList: true, subtree: true});
}

// Preserve native placement/flipping, then keep open popups above the toolbar.
const popupScopes = new WeakSet();
export function initializePopups(scope = document.body) {
  if (popupScopes.has(scope)) return;
  popupScopes.add(scope);
  const selector = '.cm-tooltip, .documentation-signature, [data-radix-popper-content-wrapper]:has([cmdk-root])';
  const popups = new Set();
  const narrowInfos = new WeakSet();
  let pending;
  const schedule = () => { pending ||= requestAnimationFrame(place); };
  const resize = new ResizeObserver(schedule);
  const place = () => {
    pending = 0;
    for (const popup of popups) if (!popup.isConnected) {
      resize.unobserve(popup); popups.delete(popup);
    }
    const viewport = window.visualViewport;
    const left = (viewport?.offsetLeft || 0) + 8;
    const top = (viewport?.offsetTop || 0) + 8;
    const right = left + (viewport?.width || innerWidth) - 16;
    const toolbar = document.querySelector('.notebook-toolbar, [data-testid="chrome-footer"]')?.getBoundingClientRect();
    const bottom = Math.min(top + (viewport?.height || innerHeight) - 16, toolbar?.top ?? innerHeight) - 8;
    for (const popup of scope.querySelectorAll(selector)) {
      if (!popups.has(popup)) { popups.add(popup); resize.observe(popup); }
      let safeLeft = left, safeRight = right;
      if (document.documentElement.dataset.documentationNotebookRuntime) {
        // Native full-app panes clip at the navigation rail. Viewport bounds
        // alone can leave readable tooltip text painted underneath that rail.
        const editor = popup.closest('.cm-editor') || document.querySelector('.documentation-full-code-frame .cm-editor.cm-focused');
        for (let pane = editor?.parentElement; pane; pane = pane.parentElement) {
          if (!/^(auto|scroll|hidden|clip)$/.test(getComputedStyle(pane).overflowX)) continue;
          const box = pane.getBoundingClientRect();
          if (!box.width) continue;
          safeLeft = Math.max(left, box.left + pane.clientLeft + 8);
          safeRight = Math.min(right, box.left + pane.clientLeft + pane.clientWidth - 8);
          break;
        }
        const width = `${Math.max(0, safeRight - safeLeft)}px`;
        if (popup.style.getPropertyValue('--documentation-popup-width') !== width) popup.style.setProperty('--documentation-popup-width', width);
      }
      const completionInfo = popup.matches('.cm-completionInfo');
      const info = completionInfo && getComputedStyle(popup).position !== 'static';
      if (info) {
        // Native side selection assumes bounded documentation. A viewport-wide
        // info box can otherwise be clamped back across its own suggestions.
        const width = `min(36rem, ${Math.max(0, safeRight - safeLeft)}px)`;
        if (popup.style.width !== width) popup.style.setProperty('width', width, 'important');
      } else if (completionInfo && popup.style.getPropertyPriority('width') === 'important') {
        popup.style.removeProperty('width');
      }
      let rect = popup.getBoundingClientRect();
      const list = popup.closest('.cm-tooltip-autocomplete');
      const options = info ? list?.querySelector('ul')?.getBoundingClientRect() : null;
      const beside = options && (popup.classList.contains('cm-completionInfo-right')
        ? safeRight - options.right - 4 >= rect.width
        : popup.classList.contains('cm-completionInfo-left') && options.left - safeLeft - 4 >= rect.width);
      const narrowInfo = info && (!beside || popup.matches('.cm-completionInfo-right-narrow, .cm-completionInfo-left-narrow'));
      if (!narrowInfo && narrowInfos.has(popup)) {
        for (const name of ['top', 'left', 'max-height']) popup.style.removeProperty(name);
        narrowInfos.delete(popup);
        rect = popup.getBoundingClientRect();
      }
      // CodeMirror parks unmeasured/hidden documentation offscreen.
      if (!rect.width || !rect.height || rect.top < -100000) continue;
      // Native narrow info panes use a viewport-sized width and the option's
      // top. Place them in the free vertical space, never over the option list.
      if (options && narrowInfo) {
        narrowInfos.add(popup);
        const below = bottom - options.bottom - 4;
        const above = options.top - top - 4;
        const down = below >= Math.min(rect.height, 160) || below >= above;
        const height = `${Math.max(0, down ? below : above)}px`;
        if (popup.style.maxHeight !== height) popup.style.maxHeight = height;
        rect = popup.getBoundingClientRect();
        const desiredTop = down ? options.bottom + 4 : options.top - rect.height - 4;
        const desiredLeft = Math.max(safeLeft, Math.min(options.left, safeRight - rect.width));
        const style = getComputedStyle(popup);
        const x = parseFloat(style.left), y = parseFloat(style.top);
        if (Math.abs(desiredLeft - rect.left) > .5) {
          popup.style.setProperty('left', `${(Number.isFinite(x) ? x : popup.offsetLeft) + desiredLeft - rect.left}px`, 'important');
          popup.style.right = 'auto';
        }
        if (Math.abs(desiredTop - rect.top) > .5) {
          popup.style.setProperty('top', `${(Number.isFinite(y) ? y : popup.offsetTop) + desiredTop - rect.top}px`, 'important');
          popup.style.bottom = 'auto';
        }
        continue;
      }
      // Documentation is positioned relative to the completion list. Measure
      // after correcting its parent so each panel stays independently visible.
      const dx = Math.max(safeLeft, Math.min(rect.left, safeRight - rect.width)) - rect.left;
      const dy = Math.max(top, Math.min(rect.top, bottom - rect.height)) - rect.top;
      // Correct layout coordinates, rather than translating: WebKit includes
      // the original untransformed box in the page's scrollable overflow.
      // offsetLeft/Top round fractional native positions to integers. Preserve
      // the resolved coordinates so repeated native measurement can settle.
      if (Math.abs(dx) > .5) {
        const current = parseFloat(getComputedStyle(popup).left);
        popup.style.setProperty('left', `${(Number.isFinite(current) ? current : popup.offsetLeft) + dx}px`, info ? 'important' : ''); popup.style.right = 'auto';
      }
      if (Math.abs(dy) > .5) {
        const current = parseFloat(getComputedStyle(popup).top);
        popup.style.setProperty('top', `${(Number.isFinite(current) ? current : popup.offsetTop) + dy}px`, info ? 'important' : ''); popup.style.bottom = 'auto';
      }
    }
  };
  new MutationObserver(records => {
    if (records.some(({type, target, addedNodes, removedNodes}) => type === 'attributes'
      ? target.matches(selector)
      : [...addedNodes, ...removedNodes].some(node => node.nodeType === 1 && (node.matches(selector) || node.querySelector(selector))))) schedule();
  }).observe(scope, {childList: true, subtree: true, attributes: true, attributeFilter: ['style', 'class']});
  addEventListener('resize', schedule);
  addEventListener('scroll', schedule, true);
  window.visualViewport?.addEventListener('resize', schedule);
  window.visualViewport?.addEventListener('scroll', schedule);
  schedule();
}
