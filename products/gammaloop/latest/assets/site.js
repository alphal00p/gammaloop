(() => {
  const root = document.documentElement;
  const body = document.body;
  const docsRoot = body.dataset.docsRoot || "";
  const menu = document.querySelector("[data-menu-toggle]");
  const backdrop = document.querySelector("[data-sidebar-backdrop]");
  const searchDialog = document.querySelector("[data-search-dialog]");
  const searchInput = document.querySelector("[data-search-input]");
  const searchResults = document.querySelector("[data-search-results]");
  const searchButtons = document.querySelectorAll("[data-search-open]");
  const themeButton = document.querySelector("[data-theme-toggle]");
  const productSelect = document.querySelector("[data-product-select]");

  const setMenuOpen = (open) => {
    body.classList.toggle("sidebar-open", open);
    menu?.setAttribute("aria-expanded", String(open));
  };
  const closeMenu = () => setMenuOpen(false);
  menu?.addEventListener("click", () => setMenuOpen(!body.classList.contains("sidebar-open")));
  backdrop?.addEventListener("click", closeMenu);
  document.querySelectorAll(".docs-sidebar a").forEach((link) => link.addEventListener("click", closeMenu));

  const storedTheme = localStorage.getItem("alphal00p-docs-theme");
  const preferredTheme = matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  root.dataset.theme = storedTheme || preferredTheme;
  themeButton?.addEventListener("click", () => {
    root.dataset.theme = root.dataset.theme === "dark" ? "light" : "dark";
    localStorage.setItem("alphal00p-docs-theme", root.dataset.theme);
  });

  productSelect?.addEventListener("change", (event) => {
    if (event.target.value) window.location.assign(event.target.value);
  });

  let indexPromise;
  const getIndex = () => indexPromise ||= fetch(`${docsRoot}search-index.json`)
    .then((response) => response.ok ? response.json() : [])
    .catch(() => []);

  const renderSearch = async (query) => {
    if (!searchResults) return;
    const normalized = query.trim().toLowerCase();
    searchResults.replaceChildren();
    if (!normalized) {
      const hint = document.createElement("li");
      hint.className = "search-empty";
      hint.textContent = "Search tutorials, manual headings, and API entries.";
      searchResults.append(hint);
      return;
    }
    const terms = normalized.split(/\s+/);
    const matches = (await getIndex()).filter((entry) => {
      const haystack = `${entry.title} ${entry.summary} ${entry.kind}`.toLowerCase();
      return terms.every((term) => haystack.includes(term));
    }).slice(0, 30);
    for (const entry of matches) {
      const item = document.createElement("li");
      item.className = "search-result";
      const link = document.createElement("a");
      link.href = `${docsRoot}${entry.href}`;
      const title = document.createElement("span");
      title.className = "search-result-title";
      title.textContent = entry.title;
      const summary = document.createElement("span");
      summary.className = "search-result-summary";
      summary.textContent = entry.summary || entry.kind;
      link.append(title, summary);
      item.append(link);
      searchResults.append(item);
    }
    if (!matches.length) {
      const empty = document.createElement("li");
      empty.className = "search-empty";
      empty.textContent = "No matching documentation.";
      searchResults.append(empty);
    }
  };

  const openSearch = () => {
    if (!searchDialog || !searchInput || !searchResults) return false;
    searchDialog.showModal();
    searchInput.focus();
    renderSearch(searchInput.value);
    return true;
  };
  searchButtons.forEach((button) => button.addEventListener("click", openSearch));
  searchInput?.addEventListener("input", (event) => renderSearch(event.target.value));
  searchDialog?.addEventListener("click", (event) => {
    if (event.target === searchDialog) searchDialog.close();
  });
  document.addEventListener("keydown", (event) => {
    const shortcut = (event.key === "/" && !/input|textarea/i.test(document.activeElement?.tagName)) ||
      (event.key.toLowerCase() === "k" && (event.metaKey || event.ctrlKey));
    if (shortcut && openSearch()) {
      event.preventDefault();
    }
  });

  const tocLinks = [...document.querySelectorAll(".toc-link")];
  const headings = tocLinks.map((link) => document.getElementById(link.hash.slice(1))).filter(Boolean);
  if (headings.length && "IntersectionObserver" in window) {
    const observer = new IntersectionObserver((entries) => {
      const visible = entries.find((entry) => entry.isIntersecting);
      if (!visible) return;
      tocLinks.forEach((link) => link.classList.toggle("active", link.hash === `#${visible.target.id}`));
    }, { rootMargin: "-20% 0px -70%" });
    headings.forEach((heading) => observer.observe(heading));
  }
})();
