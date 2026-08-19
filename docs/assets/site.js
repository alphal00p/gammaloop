(() => {
  const root = document.documentElement;
  const body = document.body;
  root.classList.add("js");
  const docsRoot = body.dataset.docsRoot || "";
  const searchIndex = body.dataset.searchIndex || `${docsRoot}search-index.json`;
  const searchRoot = body.dataset.searchRoot || docsRoot;
  const menu = document.querySelector("[data-menu-toggle]");
  const sidebar = document.querySelector(".docs-sidebar");
  const backdrop = document.querySelector("[data-sidebar-backdrop]");
  const searchDialog = document.querySelector("[data-search-dialog]");
  const searchInput = document.querySelector("[data-search-input]");
  const searchResults = document.querySelector("[data-search-results]");
  const searchButtons = document.querySelectorAll("[data-search-open]");
  const themeButton = document.querySelector("[data-theme-toggle]");
  const productSelect = document.querySelector("[data-product-select]");
  const themeColor = document.querySelector('meta[name="theme-color"]');
  const mobileNavigation = matchMedia("(max-width: 52rem)");

  const setMenuOpen = (requestedOpen, focusNavigation = false) => {
    const open = requestedOpen && mobileNavigation.matches;
    const hidden = mobileNavigation.matches && !open;
    body.classList.toggle("sidebar-open", open);
    sidebar?.toggleAttribute("inert", hidden);
    if (hidden) sidebar?.setAttribute("aria-hidden", "true");
    else sidebar?.removeAttribute("aria-hidden");
    menu?.setAttribute("aria-expanded", String(open));
    menu?.setAttribute("aria-label", open ? "Close navigation" : "Open navigation");
    if (open && focusNavigation) {
      requestAnimationFrame(() => {
        const target = sidebar?.querySelector('[aria-current="page"]')
          || sidebar?.querySelector("[aria-current]")
          || sidebar?.querySelector("a");
        target?.focus();
      });
    }
  };
  const closeMenu = () => setMenuOpen(false);
  setMenuOpen(false);
  mobileNavigation.addEventListener("change", () => setMenuOpen(false));
  menu?.addEventListener("click", () => setMenuOpen(!body.classList.contains("sidebar-open"), true));
  backdrop?.addEventListener("click", () => {
    closeMenu();
    menu?.focus();
  });
  document.querySelectorAll(".docs-sidebar a").forEach((link) => link.addEventListener("click", closeMenu));
  const currentSidebarLink = sidebar?.querySelector('[aria-current="page"]');
  requestAnimationFrame(() => {
    if (!sidebar || !currentSidebarLink) return;
    const linkTop = currentSidebarLink.offsetTop;
    const linkBottom = linkTop + currentSidebarLink.offsetHeight;
    if (linkTop < sidebar.scrollTop || linkBottom > sidebar.scrollTop + sidebar.clientHeight) {
      sidebar.scrollTop = Math.max(0, linkTop - sidebar.clientHeight / 2);
    }
  });

  const storedTheme = localStorage.getItem("alphal00p-docs-theme");
  const preferredTheme = matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  const applyTheme = (theme) => {
    root.dataset.theme = theme;
    const nextTheme = theme === "dark" ? "light" : "dark";
    themeButton?.setAttribute("aria-label", `Switch to ${nextTheme} theme`);
    themeButton?.setAttribute("aria-pressed", String(theme === "dark"));
    themeColor?.setAttribute("content", theme === "dark" ? "#211a23" : "#f9f6f0");
  };
  applyTheme(storedTheme || preferredTheme);
  themeButton?.addEventListener("click", () => {
    const theme = root.dataset.theme === "dark" ? "light" : "dark";
    applyTheme(theme);
    localStorage.setItem("alphal00p-docs-theme", theme);
  });

  productSelect?.addEventListener("change", (event) => {
    if (event.target.value) window.location.assign(event.target.value);
  });

  let indexPromise;
  const getIndex = () => indexPromise ||= fetch(searchIndex)
    .then((response) => response.ok ? response.json() : [])
    .catch(() => []);

  const normalizeSearch = (value) => String(value || "").normalize("NFKD").toLowerCase();
  const searchScore = (entry, query, terms) => {
    const title = normalizeSearch(entry.title);
    const summary = normalizeSearch(entry.summary);
    const kind = normalizeSearch(entry.kind);
    const text = normalizeSearch(entry.text);
    const haystack = `${title} ${summary} ${kind} ${text}`;
    if (!terms.every((term) => haystack.includes(term))) return null;

    const words = title.split(/[^a-z0-9_.:-]+/).filter(Boolean);
    let score = title === query ? 1200 : title.startsWith(query) ? 700 : title.includes(query) ? 420 : 0;
    for (const term of terms) {
      if (words.includes(term)) score += 150;
      else if (title.includes(term)) score += 90;
      if (kind.includes(term)) score += 35;
      if (summary.includes(term)) score += 15;
      if (text.includes(term)) score += 3;
    }
    if (/command|setting/.test(kind)) score += 120;
    else if (/product|tutorial|manual|guide|reference|version history|python-api|rust-api/.test(kind)) score += 80;
    else if (!/developer/.test(kind)) score += 30;
    if (/developer/.test(kind)) {
      score -= / heading/.test(kind) ? 100 : 60;
      if (!/developer current /.test(kind)) score -= 120;
    }
    return score;
  };

  const renderSearch = async (query) => {
    if (!searchResults) return;
    const normalized = normalizeSearch(query.trim());
    searchResults.replaceChildren();
    if (!normalized) {
      const hint = document.createElement("li");
      hint.className = "search-empty";
      hint.textContent = body.dataset.searchIndex
        ? "Search all five projects and the developer architecture notes."
        : "Search tutorials, manual headings, and API entries.";
      searchResults.append(hint);
      return;
    }
    const terms = normalized.split(/\s+/);
    const matches = (await getIndex())
      .map((entry, index) => ({ entry, index, score: searchScore(entry, normalized, terms) }))
      .filter((match) => match.score !== null)
      .sort((left, right) => right.score - left.score ||
        left.entry.title.length - right.entry.title.length ||
        left.entry.title.localeCompare(right.entry.title) || left.index - right.index)
      .slice(0, 30)
      .map((match) => match.entry);
    for (const entry of matches) {
      const item = document.createElement("li");
      item.className = "search-result";
      const link = document.createElement("a");
      link.href = `${searchRoot}${entry.href}`;
      const title = document.createElement("span");
      title.className = "search-result-title";
      title.textContent = entry.title;
      const kind = document.createElement("span");
      kind.className = "search-result-kind";
      kind.textContent = entry.kind;
      const summary = document.createElement("span");
      summary.className = "search-result-summary";
      summary.textContent = entry.summary || entry.kind;
      link.append(title, kind, summary);
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
    if (event.key === "Escape" && body.classList.contains("sidebar-open")) {
      closeMenu();
      menu?.focus();
      return;
    }
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

  const publicationFilters = document.querySelector("[data-publication-filters]");
  if (publicationFilters) {
    const search = publicationFilters.querySelector("[data-publication-search]");
    const author = publicationFilters.querySelector("[data-publication-author]");
    const year = publicationFilters.querySelector("[data-publication-year]");
    const type = publicationFilters.querySelector("[data-publication-type]");
    const sort = publicationFilters.querySelector("[data-publication-sort]");
    const count = publicationFilters.querySelector("[data-publication-count]");
    const list = document.querySelector("[data-publication-list]");
    const cards = [...document.querySelectorAll("[data-publication]")];
    const params = new URLSearchParams(location.search);
    search.value = params.get("q") || "";
    author.value = params.get("author") || "";
    year.value = params.get("year") || "";
    type.value = params.get("type") || "";
    sort.value = params.get("sort") || "newest";

    const updatePublications = () => {
      const query = search.value.trim().toLowerCase();
      const matches = cards.filter((card) => {
        const people = card.dataset.people.split(/\s+/);
        const types = card.dataset.types.split("|");
        return (!query || card.textContent.toLowerCase().includes(query)) &&
          (!author.value || people.includes(author.value)) &&
          (!year.value || card.dataset.year === year.value) &&
          (!type.value || types.includes(type.value));
      });
      const visible = new Set(matches);
      cards.forEach((card) => card.hidden = !visible.has(card));
      const ordered = [...cards].sort((left, right) => sort.value === "cited"
        ? Number(right.dataset.citations) - Number(left.dataset.citations) || right.dataset.date.localeCompare(left.dataset.date)
        : right.dataset.date.localeCompare(left.dataset.date));
      ordered.forEach((card) => list.append(card));
      count.textContent = `${matches.length} publication${matches.length === 1 ? "" : "s"}`;

      const next = new URLSearchParams();
      if (query) next.set("q", search.value.trim());
      if (author.value) next.set("author", author.value);
      if (year.value) next.set("year", year.value);
      if (type.value) next.set("type", type.value);
      if (sort.value !== "newest") next.set("sort", sort.value);
      history.replaceState(null, "", `${location.pathname}${next.size ? `?${next}` : ""}${location.hash}`);
    };
    publicationFilters.addEventListener("input", updatePublications);
    publicationFilters.addEventListener("change", updatePublications);
    updatePublications();
  }

  document.querySelectorAll("[data-copy-target]").forEach((button) => {
    button.addEventListener("click", async () => {
      const target = document.getElementById(button.dataset.copyTarget);
      if (!target) return;
      try {
        await navigator.clipboard.writeText(target.textContent.trim());
        const label = button.textContent;
        button.textContent = "Copied";
        setTimeout(() => button.textContent = label, 1500);
      } catch {
        target.focus?.();
      }
    });
  });

  document.querySelectorAll("[data-reference-filter]").forEach((input) => {
    const root = input.closest("[data-reference-filter-root]");
    const scope = root?.querySelector("[data-reference-filter-scope]") ||
      document.querySelector("[data-reference-filter-scope]");
    if (!scope) return;
    const entries = [...scope.querySelectorAll("[data-reference-entry]")];
    const groups = [...scope.querySelectorAll("[data-reference-group]")];
    const count = input.parentElement?.querySelector("[data-reference-filter-count]") ||
      document.querySelector("[data-reference-filter-count]");
    const initialQuery = new URLSearchParams(location.search).get("q");
    if (initialQuery && !input.value) input.value = initialQuery;
    const updateReference = () => {
      const terms = input.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
      let visible = 0;
      entries.forEach((entry) => {
        const haystack = (entry.dataset.referenceSearch || entry.textContent).toLowerCase();
        const match = terms.every((term) => haystack.includes(term));
        entry.hidden = !match;
        if (match) visible += 1;
      });
      groups.forEach((group) => {
        group.hidden = !group.querySelector("[data-reference-entry]:not([hidden])");
      });
      if (count) {
        count.textContent = `${visible} of ${entries.length} entr${entries.length === 1 ? "y" : "ies"}`;
      }
    };
    input.addEventListener("input", () => {
      const params = new URLSearchParams(location.search);
      if (input.value.trim()) params.set("q", input.value.trim());
      else params.delete("q");
      history.replaceState(null, "", `${location.pathname}${params.size ? `?${params}` : ""}${location.hash}`);
      updateReference();
    });
    updateReference();
  });

  const revealReferenceFragment = () => {
    if (!location.hash) return;
    const target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
    if (!target) return;
    if (target.matches("a[data-reference-redirect]") && target.href) {
      location.replace(target.href);
      return;
    }
    if (target.matches("details")) target.open = true;
    let parent = target.parentElement?.closest("details");
    while (parent) {
      parent.open = true;
      parent = parent.parentElement?.closest("details");
    }
  };
  addEventListener("hashchange", revealReferenceFragment);
  revealReferenceFragment();
})();
