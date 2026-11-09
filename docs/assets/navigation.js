// The site bar and local contents drawer share one responsive state owner.
export function initializeNavigation() {
  const body = document.body;
  document.documentElement.classList.add('js');
  const narrow = matchMedia('(max-width: 52rem)');
  let focusedToggle;
  const globalToggle = document.querySelector('[data-portal-menu-toggle]');
  const globalNav = document.querySelector('#portal-navigation');
  const localToggle = document.querySelector('[data-menu-toggle]');
  const sidebar = document.querySelector('.docs-sidebar');
  const backdrop = document.querySelector('[data-sidebar-backdrop]');
  if (backdrop) backdrop.tabIndex = -1;
  const background = document.querySelectorAll('.docs-main, .docs-toc, .skip-link');
  const setGlobal = (requested, restore = false) => {
    const open = requested && narrow.matches;
    if (open) setLocal(false);
    const focused = globalNav?.contains(document.activeElement);
    globalNav?.classList.toggle('portal-nav-open', open);
    globalNav?.toggleAttribute('inert', narrow.matches && !open);
    globalToggle?.setAttribute('aria-expanded', String(open));
    globalToggle?.setAttribute('aria-label', open ? 'Close site navigation' : 'Open site navigation');
    if (!open && restore && focused && narrow.matches) globalToggle?.focus();
  };
  const setLocal = (requested, focus = false) => {
    const open = requested && narrow.matches;
    if (open) setGlobal(false);
    body.classList.toggle('sidebar-open', open);
    sidebar?.toggleAttribute('inert', narrow.matches && !open);
    localToggle?.setAttribute('aria-expanded', String(open));
    localToggle?.setAttribute('aria-label', open ? 'Close documentation contents' : 'Open documentation contents');
    background.forEach(element => element.toggleAttribute('inert', open));
    if (open && focus) {
      const target = sidebar?.querySelector('[aria-current="page"]') || sidebar?.querySelector('a, button, select');
      target?.focus({preventScroll: true});
      if (target) {
        const box = target.getBoundingClientRect();
        const pane = sidebar.getBoundingClientRect();
        if (box.top < pane.top || box.bottom > pane.bottom) sidebar.scrollTop += box.top - pane.top - pane.height / 2;
      }
    }
  };
  // Enhance the existing groups only after JS loads: every link remains available without JS.
  sidebar?.querySelectorAll('.sidebar-group').forEach((group, index) => {
    const title = group.querySelector(':scope > .sidebar-group-title');
    if (!title) return;
    const items = document.createElement('div');
    items.id = `contents-group-${index}`;
    while (title.nextSibling) items.append(title.nextSibling);
    const button = document.createElement('button');
    button.type = 'button';
    button.append(...title.childNodes);
    button.setAttribute('aria-controls', items.id);
    const expanded = index === 0 || !!group.querySelector('[aria-current]') || !!items.querySelector('[aria-current]');
    button.setAttribute('aria-expanded', String(expanded));
    items.hidden = !expanded;
    button.addEventListener('click', () => {
      items.hidden = !items.hidden;
      button.setAttribute('aria-expanded', String(!items.hidden));
    });
    title.append(button);
    group.append(items);
  });
  requestAnimationFrame(() => {
    const current = sidebar?.querySelector('[aria-current="page"]');
    if (!current || narrow.matches) return;
    const box = current.getBoundingClientRect();
    const pane = sidebar.getBoundingClientRect();
    if (box.top < pane.top || box.bottom > pane.bottom) sidebar.scrollTop += box.top - pane.top - pane.height / 2;
  });
  const close = () => {
    const localFocused = sidebar?.contains(document.activeElement);
    setLocal(false);
    setGlobal(false, true);
    if (localFocused && narrow.matches) localToggle?.focus();
  };
  close();
  narrow.addEventListener('change', () => {
    const globalFocused = document.activeElement === globalToggle || focusedToggle === globalToggle;
    const localFocused = document.activeElement === localToggle || focusedToggle === localToggle;
    close();
    if (!narrow.matches && globalFocused) (globalNav?.querySelector('[aria-current]') || globalNav?.querySelector('a'))?.focus();
    if (!narrow.matches && localFocused) (sidebar?.querySelector('[aria-current="page"]') || sidebar?.querySelector('a, button, select'))?.focus();
  });
  globalToggle?.addEventListener('click', () => {
    const open = globalToggle.getAttribute('aria-expanded') !== 'true';
    setGlobal(open);
    if (open) (globalNav?.querySelector('[aria-current]') || globalNav?.querySelector('a'))?.focus();
  });
  localToggle?.addEventListener('click', () => setLocal(!body.classList.contains('sidebar-open'), true));
  backdrop?.addEventListener('click', () => { setLocal(false); localToggle?.focus(); });
  sidebar?.addEventListener('click', event => { if (event.target.closest('a')) setLocal(false); });
  globalNav?.addEventListener('click', event => { if (event.target.closest('a')) setGlobal(false); });
  document.addEventListener('click', event => {
    if (!event.target.closest('.portal-header')) setGlobal(false, true);
  });
  document.addEventListener('focusin', event => {
    focusedToggle = event.target === globalToggle || event.target === localToggle ? event.target : null;
    if (!event.target.closest('.portal-header')) setGlobal(false);
    if (body.classList.contains('sidebar-open') && !event.target.closest('.docs-sidebar, .docs-local-bar, [data-sidebar-backdrop]')) setLocal(false);
  });
  document.addEventListener('keydown', event => {
    if (event.key !== 'Escape') return;
    if (body.classList.contains('sidebar-open')) { setLocal(false); localToggle?.focus(); event.preventDefault(); }
    else if (globalToggle?.getAttribute('aria-expanded') === 'true') { setGlobal(false); globalToggle.focus(); event.preventDefault(); }
  });
  return {close};
}
