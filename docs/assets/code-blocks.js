// Static source frames never install an editor or execution control.
export function frameSource(pre) {
  if (pre.closest('.cm-tooltip, .documentation-type-info')) return;
  let row = pre.closest('.reference-copy-row, .example-wrap, .codehilite, .code-frame');
  if (!row) {
    row = document.createElement('div');
    pre.before(row);
    row.append(pre);
  }
  if (row.classList.contains('code-frame')) return row;
  row.classList.add('code-frame');
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
    const toolbar = document.querySelector('.notebook-toolbar')?.getBoundingClientRect();
    const bottom = Math.min(top + (viewport?.height || innerHeight) - 16, toolbar?.top ?? innerHeight) - 8;
    for (const popup of scope.querySelectorAll(selector)) {
      if (!popups.has(popup)) { popups.add(popup); resize.observe(popup); }
      const rect = popup.getBoundingClientRect();
      // CodeMirror parks unmeasured/hidden documentation offscreen.
      if (!rect.width || !rect.height || rect.top < -100000) continue;
      // Documentation is positioned relative to the completion list. Measure
      // after correcting its parent so each panel stays independently visible.
      const dx = Math.max(left, Math.min(rect.left, right - rect.width)) - rect.left;
      const dy = Math.max(top, Math.min(rect.top, bottom - rect.height)) - rect.top;
      // Correct layout coordinates, rather than translating: WebKit includes
      // the original untransformed box in the page's scrollable overflow.
      if (Math.abs(dx) > .5) { popup.style.left = `${popup.offsetLeft + dx}px`; popup.style.right = 'auto'; }
      if (Math.abs(dy) > .5) { popup.style.top = `${popup.offsetTop + dy}px`; popup.style.bottom = 'auto'; }
    }
  };
  new MutationObserver(records => {
    if (records.some(({type, target, addedNodes, removedNodes}) => type === 'attributes'
      ? target.matches(selector)
      : [...addedNodes, ...removedNodes].some(node => node.nodeType === 1 && (node.matches(selector) || node.querySelector(selector))))) schedule();
  }).observe(scope, {childList: true, subtree: true, attributes: true, attributeFilter: ['style']});
  addEventListener('resize', schedule);
  addEventListener('scroll', schedule, true);
  window.visualViewport?.addEventListener('resize', schedule);
  window.visualViewport?.addEventListener('scroll', schedule);
  schedule();
}
