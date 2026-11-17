(() => {
  const root = document.documentElement;
  const body = document.body;
  root.classList.add("js");
  const docsRoot = body.dataset.docsRoot || "";
  const searchIndex = body.dataset.searchIndex || `${docsRoot}search-index.json`;
  const searchRoot = body.dataset.searchRoot || docsRoot;

  // Marimo copies same-origin sheets into editor shadows but strips @import.
  // Expose the shared source owner directly, and load it before mounting views.
  const sourceStylesReady = Promise.all([...document.querySelectorAll('link[rel="stylesheet"]')]
    .filter(link => new URL(link.href).pathname.endsWith('/assets/site.css'))
    .map(link => {
      // One preferred set keeps every sheet active, including in WebKit.
      // The marimo prefix also lets the runtime discover versioned CSS URLs.
      link.title = 'marimo';
      const source = document.createElement('link');
      source.rel = 'stylesheet'; source.title = 'marimo';
      const sourceURL = new URL('code.css', link.href);
      sourceURL.search = new URL(link.href).search;
      source.href = sourceURL;
      const loaded = new Promise(resolve => {
        source.onload = () => resolve(true);
        source.onerror = () => resolve(false);
      });
      link.before(source);
      return loaded;
    }));

  // The front-page tour loads its stage and styles only where a host asks for them.
  if (document.querySelector("[data-showcase]")) {
    const stylesheet = document.createElement("link");
    stylesheet.rel = "stylesheet";
    stylesheet.href = `${docsRoot}assets/showcase.css`;
    document.head.append(stylesheet);
    const script = document.createElement("script");
    script.src = `${docsRoot}assets/showcase.js`;
    document.head.append(script);
  }
  const searchDialog = document.querySelector("[data-search-dialog]");
  const searchInput = document.querySelector("[data-search-input]");
  const searchResults = document.querySelector("[data-search-results]");
  const searchStatus = searchResults ? document.createElement("span") : null;
  if (searchStatus) {
    searchStatus.className = "search-status";
    searchStatus.setAttribute("role", "status");
    searchResults.before(searchStatus);
  }
  const searchButtons = document.querySelectorAll("[data-search-open]");
  const themeButton = document.querySelector("[data-theme-toggle]");
  const productSelect = document.querySelector("[data-product-select]");
  const themeColor = document.querySelector('meta[name="theme-color"]');
  let navigation;
  const siteScript = [...document.scripts].find(script => new URL(script.src || location.href).pathname.endsWith('/assets/site.js'));
  import(new URL('navigation.js', siteScript.src)).then(module => {
    navigation = module.initializeNavigation();
  });

  // Islands cache their initial theme. Keep widget shadow-root boundaries in
  // sync so an initial `.light` class cannot override the documentation theme.
  const notebookThemeRoots = new WeakSet();
  const notebookTables = new WeakMap();
  let notebookIntegration;
  const syncNotebookTheme = () => {
    const theme = root.dataset.theme;
    document.querySelectorAll(".live-notebook marimo-code-editor").forEach((editor) => {
      const value = JSON.stringify(theme);
      if (editor.dataset.theme !== value) editor.dataset.theme = value;
    });
    document.querySelectorAll(".live-notebook *").forEach((widget) => {
      const shadow = widget.shadowRoot;
      if (!shadow) return;
      if (!notebookThemeRoots.has(shadow)) {
        notebookThemeRoots.add(shadow);
        notebookThemeObserver.observe(shadow, { childList: true, subtree: true });
      }
      shadow.querySelectorAll(".marimo > .contents").forEach((boundary) => {
        boundary.dataset.theme = theme;
        boundary.style.colorScheme = theme;
        boundary.classList.toggle("dark", theme === "dark");
        boundary.classList.toggle("light", theme === "light");
      });
      shadow.querySelectorAll(".documentation-type-info").forEach((popup) => {
        if (popup.dataset.codeTheme !== theme) popup.dataset.codeTheme = theme;
      });
      if (widget.matches("marimo-table")) {
        const table = shadow.querySelector("table");
        const scroll = table?.parentElement;
        if (table && scroll && notebookTables.get(widget)?.table !== table) {
          notebookTables.get(widget)?.observer.disconnect();
          const hint = shadow.querySelector(".notebook-table-scroll-hint") || document.createElement("p");
          hint.className = "notebook-table-scroll-hint";
          hint.textContent = "Scroll table horizontally →";
          shadow.querySelector('[part="table-wrapper"]').after(hint);
          const update = () => {
            if (!widget.isConnected) { observer.disconnect(); return; }
            const overflow = scroll.scrollWidth > scroll.clientWidth + 1;
            widget.toggleAttribute("data-notebook-table-overflow", overflow);
            widget.toggleAttribute("data-notebook-table-quiet", !overflow && widget.dataset.totalRows === "1" &&
              widget.dataset.pagination === "false" && !widget.hasAttribute("data-selection"));
            hint.hidden = !overflow;
            if (overflow) {
              scroll.tabIndex = 0;
              scroll.setAttribute("role", "region");
              scroll.setAttribute("aria-label", "Table; scroll horizontally to see all columns");
            } else {
              scroll.removeAttribute("tabindex");
              scroll.removeAttribute("role");
              scroll.removeAttribute("aria-label");
            }
          };
          let resizeFrame;
          const observer = new ResizeObserver(() => {
            cancelAnimationFrame(resizeFrame);
            resizeFrame = requestAnimationFrame(update);
          });
          observer.observe(scroll); observer.observe(table);
          notebookTables.set(widget, {table, observer});
          update();
        }
      }
      if (widget.matches("marimo-code-editor")) {
        const source = widget.matches("marimo-island > marimo-ui-element > marimo-code-editor");
        if (source) widget.setAttribute("data-notebook-source", "");
        notebookIntegration?.styleEditor(shadow, source);
      }
    });
  };
  const notebookThemeObserver = new MutationObserver(syncNotebookTheme);
  notebookThemeObserver.observe(root, { attributes: true, attributeFilter: ["data-theme"] });
  document.querySelectorAll(".live-notebook").forEach((notebook) => {
    notebookThemeObserver.observe(notebook, { childList: true, subtree: true });
  });

  let explicitTheme;
  try {
    const stored = localStorage.getItem("alphal00p-docs-theme");
    if (stored === "dark" || stored === "light") explicitTheme = stored;
  } catch { /* Theme selection remains usable when storage is unavailable. */ }
  const syncExplicitTheme = (force = false) => {
    if (!explicitTheme) return;
    try {
      const nativeTheme = localStorage.getItem("rustdoc-theme");
      const lastSynced = localStorage.getItem("alphal00p-rustdoc-theme-synced");
      // A later native Settings choice belongs to Rustdoc until the next site toggle.
      const followsSite = nativeTheme === lastSynced
        && localStorage.getItem("rustdoc-use-system-theme") === "false";
      if (force || nativeTheme === null || followsSite) {
        localStorage.setItem("rustdoc-use-system-theme", "false");
        localStorage.setItem("rustdoc-theme", explicitTheme);
        localStorage.setItem("alphal00p-rustdoc-theme-synced", explicitTheme);
      }
    } catch { /* Native references retain their defaults without storage. */ }
    document.querySelectorAll("a.notebook-full-link").forEach((link) => {
      const url = new URL(link.href);
      url.searchParams.set("theme", explicitTheme);
      link.href = url.href;
    });
  };
  const preferredTheme = matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  const applyTheme = (theme, forceSync = false) => {
    root.dataset.theme = theme;
    syncNotebookTheme();
    syncExplicitTheme(forceSync);
    const nextTheme = theme === "dark" ? "light" : "dark";
    themeButton?.setAttribute("aria-label", `Switch to ${nextTheme} theme`);
    themeButton?.setAttribute("aria-pressed", String(theme === "dark"));
    themeColor?.setAttribute("content", theme === "dark" ? "#211a23" : "#f9f6f0");
  };
  applyTheme(explicitTheme || preferredTheme);
  themeButton?.addEventListener("click", () => {
    const theme = root.dataset.theme === "dark" ? "light" : "dark";
    explicitTheme = theme;
    applyTheme(theme, true);
    try { localStorage.setItem("alphal00p-docs-theme", theme); }
    catch { /* Keep the chosen appearance for this page. */ }
  });

  productSelect?.addEventListener("change", (event) => {
    if (event.target.value) window.location.assign(event.target.value);
  });

  let documentationConfiguration;
  const loadDocumentationConfiguration = () => documentationConfiguration ||= (async () => {
    const script = document.querySelector('script[src*="/assets/site.js"]');
    const url = new URL('documentation-cells-config.json', script.src);
    // This mutable pointer selects immutable runtime assets after publication.
    const response = await fetch(url, {cache: 'no-cache'});
    if (!response.ok) throw Error(`Documentation cell configuration: HTTP ${response.status}`);
    return {url, value: await response.json()};
  })();
  // Observe the pinned scientific worker without rewriting immutable runtime
  // assets or rewrapping fetch responses (which loses WASM code-cache metadata).
  const notebookStartup = {phase: "Starting notebook worker…", listeners: new Set()};
  let rejectNotebookStartup;
  let resolveNotebookExecution;
  notebookStartup.failure = new Promise((_, reject) => { rejectNotebookStartup = reject; });
  // Failure may arrive before the dynamic import finishes.
  notebookStartup.failure.catch(() => {});
  const reportNotebookStartup = (phase) => {
    notebookStartup.phase = phase;
    for (const update of notebookStartup.listeners) update(phase);
  };
  const observeNotebookWorker = (runtimeURL) => {
    const workerURL = new URL("assets/worker-D7hWi02M.js", runtimeURL).href;
    const OriginalWorker = window.Worker;
    // This wrapper executes before the original module's native-loader import.
    // Preserve incoming RPC messages until that module has installed handlers.
    const observe = async (url) => {
      const send = (value) => self.postMessage({notebookStartup: value});
      const queued = [];
      const hold = (event) => { event.stopImmediatePropagation(); queued.push(event.data); };
      self.addEventListener("message", hold);
      let phase = "Loading notebook worker module…";
      const pending = new Map();
      let nextRequest = 0;
      new PerformanceObserver(list => {
        for (const resource of list.getEntries()) {
          for (const [id, request] of pending) {
            if (resource.name === request.url && resource.startTime >= request.start - 1) pending.delete(id);
          }
        }
      }).observe({entryTypes: ["resource"]});
      const fetchOriginal = self.fetch;
      self.fetch = async (...args) => {
        const url = String(args[0]?.url || args[0]);
        const id = ++nextRequest;
        pending.set(id, {url, start: performance.now()});
        try {
          const response = await fetchOriginal(...args);
          if (!response.ok) {
            phase = `Runtime download returned HTTP ${response.status}: ${new URL(url, self.location.href).pathname}`;
            send({phase});
          }
          if (pending.has(id)) pending.get(id).headers = true;
          return response;
        } catch (error) {
          pending.delete(id);
          send({error: `Runtime download failed: ${error.message}`});
          throw error;
        }
      };
      self.addEventListener("unhandledrejection", (event) => send({error: String(event.reason?.message || event.reason)}));
      const errorOriginal = console.error;
      console.error = (...args) => {
        // These pinned Marimo catches log and swallow failed dependency promises,
        // leaving the session mounted but never instantiating its cells.
        if (["Failed to load notebook dependencies", "Failed to load packages from micropip"].includes(args[0])) {
          send({error: args.map(value => value?.message || String(value)).join(": ")});
        }
        errorOriginal.apply(console, args);
      };
      const compileOriginal = WebAssembly.compileStreaming;
      WebAssembly.compileStreaming = async (...args) => {
        phase = "Downloading and compiling scientific WebAssembly…";
        send({phase});
        try { return await compileOriginal.apply(WebAssembly, args); }
        finally { phase = "Initializing Python and installing runtime packages…"; }
      };
      const timer = setInterval(() => {
        const span = self.t?.getSpans?.().findLast(span => !span.endTime);
        const stages = {loadPyodide: phase, "startSession.runPython": "Importing Marimo and creating the Python session…",
          "pyodide.loadPackage": "Installing notebook Python packages…", "micropip.install": "Installing notebook dependencies…"};
        const downloads = [...pending.values()].map(({url, start, headers}) =>
          `${new URL(url, self.location.href).pathname.split("/").pop()} (${Math.floor((performance.now() - start) / 1000)} s${headers ? " headers received; awaiting completion" : " awaiting response headers"})`);
        send({phase: downloads.length ? `Downloading ${downloads.join(", ")}…` : stages[span?.name] || phase});
      }, 1000);
      self.addEventListener("message", (event) => {
        if (event.data?.notebookStartupStop) { event.stopImmediatePropagation(); clearInterval(timer); }
      });
      try {
        await import(url);
        send({loaded: true});
        self.removeEventListener("message", hold);
        for (const data of queued) self.dispatchEvent(new MessageEvent("message", {data}));
      } catch (error) {
        clearInterval(timer);
        send({error: `Notebook worker failed to load: ${error.message}`});
      }
    };
    window.Worker = class extends OriginalWorker {
      postMessage(message, ...options) {
        // Observe the exact code handed to the pinned Marimo worker, before
        // later edits or queued/running UI updates can change editor state.
        const request = message?.data;
        const command = request?.method === "bridge" && request.params?.functionName === "put_control_request"
          ? request.params.payload : null;
        if (command?.type === "execute-cells" && Array.isArray(command.cellIds) && Array.isArray(command.codes)) {
          document.dispatchEvent(new CustomEvent("notebook-execution-request", {
            detail: {appId: request.params.appId, cellIds: command.cellIds, codes: command.codes},
          }));
        }
        return super.postMessage(message, ...options);
      }
      constructor(url, options) {
        if (String(url) !== workerURL) { super(url, options); return; }
        notebookStartup.executing = new Promise(resolve => { resolveNotebookExecution = resolve; });
        const wrapper = URL.createObjectURL(new Blob([`(${observe.toString()})(${JSON.stringify(workerURL)})`], {type: "text/javascript"}));
        super(wrapper, {...options, type: "module"});
        let released = false;
        const release = () => { if (!released) { released = true; URL.revokeObjectURL(wrapper); } };
        const fail = (message) => { release(); rejectNotebookStartup(Error(message)); };
        this.addEventListener("error", (event) => fail(`Notebook worker failed: ${event.message || "module could not be loaded"}`));
        this.addEventListener("messageerror", () => fail("Notebook worker sent an unreadable message"));
        this.addEventListener("message", ({data}) => {
          if (data?.notebookStartup) {
            const value = data.notebookStartup;
            if (value.loaded) release();
            if (value.error) fail(value.error);
            if (value.phase) reportNotebookStartup(value.phase);
          }
          const rpc = data?.data;
          if (rpc?.id === "initializedError") fail(`Python initialization failed: ${rpc.payload?.error}`);
          if (rpc?.id === "kernelMessage") {
            const message = JSON.parse(rpc.payload.message);
            if (message.op === "cell-op" && ["queued", "running"].includes(message.data.status)) {
              resolveNotebookExecution();
              this.postMessage({notebookStartupStop: true});
              release();
            }
          }
        });
      }
    };
  };
  let notebookRuntime;
  let notebookRuntimeURL;
  let notebookIntegrationPromise;
  // Marimo 0.25.1 islands expose neither cell timings nor syntax extensions.
  // Keep this pinned bundle boundary together; review it on runtime upgrades.
  const loadNotebookIntegration = (runtimeURL) => notebookIntegrationPromise ||= Promise.all([
    import(new URL("useTheme-DzXmYZ5V.js", runtimeURL).href),
    import(new URL("cells-Cv93oav_.js", runtimeURL).href),
    loadDocumentationConfiguration().then(async ({url, value}) => {
      const module = await import(new URL('documentation-editor.js', new URL(value.compilerRoot, url)));
      return module.load(runtimeURL);
    }),
  ]).then(([{ _: store }, { d: cells }, editor]) => {
    const {EditorView, StateEffect, StateField, syntaxHighlighting, classHighlighter, highlightActiveLine, highlightActiveLineGutter, outsetGutter} = editor;
    if (!store?.get(cells)?.cellRuntime || !store.sub || !EditorView?.findFromDOM ||
        !StateEffect?.appendConfig || !StateField?.define || !syntaxHighlighting || !classHighlighter ||
        !highlightActiveLine || !highlightActiveLineGutter || !outsetGutter) {
      throw new Error("Review the pinned Marimo editor and timing integration");
    }
    // React can reconfigure an existing view on theme changes, dropping appended
    // extensions. Track the extension in editor state, rather than view identity.
    const styled = StateField.define({ create: () => true, update: () => true });
    notebookIntegration = {
      cellRuntime: (id) => store.get(cells).cellRuntime[id],
      editorSource: (element) => element && EditorView.findFromDOM(element)?.state.doc.toString(),
      focusEditor: (element) => {
        const editor = element.closest?.(".cm-editor");
        const view = editor && EditorView.findFromDOM(editor);
        if (!view) return false;
        view.focus();
        return true;
      },
      subscribe: (update) => store.sub(cells, update),
      styleEditor: (shadow, source) => {
        const element = shadow.querySelector(".cm-editor");
        const view = element && EditorView.findFromDOM(element);
        if (!view || view.state.field(styled, false)) return;
        view.dispatch({ effects: StateEffect.appendConfig.of([styled, syntaxHighlighting(classHighlighter),
          ...(source ? [highlightActiveLine(), highlightActiveLineGutter(), outsetGutter,
            EditorView.updateListener.of(update => {
              if (update.docChanged) shadow.host.dispatchEvent(new CustomEvent("notebook-source-change", {
                bubbles: true, composed: true, detail: update.state.doc.toString(),
              }));
            })] : [])]) });
        if (source) document.dispatchEvent(new CustomEvent("documentation-editor-ready", {detail: {editor: shadow.host}}));
      },
    };
    syncNotebookTheme();
    return notebookIntegration;
  }).catch((error) => {
    // An unavailable integration must not stop Python or invent durations.
    console.warn("Notebook editor styling and cell timings unavailable", error);
  });
  let notebookElementId = 0;
  const loadNotebook = async (container) => {
    container.dataset.notebookState = "loading";
    const status = document.createElement("p");
    status.setAttribute("role", "status");
    status.setAttribute("aria-atomic", "true");
    status.dataset.state = "busy";
    status.textContent = "Downloading notebook…";
    const notebook = document.createElement("div");
    notebook.className = "live-notebook-cells";
    notebook.setAttribute("aria-busy", "true");
    let started = false;
    let restarting = false;
    let runGeneration = 0;
    let showSetup = false;
    let unsubscribeTimings;
    const editedCode = new Map();
    const lastRunCode = new WeakMap();
    const currentCode = (island) => {
      const editor = island.querySelector(sourceSelector);
      const dom = editor?.shadowRoot?.querySelector(".cm-editor");
      return notebookIntegration?.editorSource(dom) ?? editedCode.get(Number(island.dataset.cellIdx))
        ?? JSON.parse(editor?.dataset.initialValue || '""');
    };
    const markEdited = (island, code) => {
      if (code === lastRunCode.get(island)) delete island.dataset.notebookEdited;
      else island.dataset.notebookEdited = "true";
    };
    const sourceSelector = ":scope > marimo-ui-element > marimo-code-editor";
    const runSelector = ':scope > .relative > .absolute button:has(svg.lucide-play)';
    const toolbar = document.createElement("div");
    toolbar.className = "notebook-toolbar";
    toolbar.setAttribute("role", "region");
    toolbar.setAttribute("aria-label", "Notebook runtime controls");
    const more = document.createElement("details");
    more.className = "notebook-more";
    const summary = document.createElement("summary");
    summary.textContent = "More";
    summary.title = "Notebook options and help";
    const options = document.createElement("div");
    options.className = "notebook-options";
    more.append(summary, options);
    const closeOptions = (restoreFocus = false) => {
      if (!more.open) return;
      more.open = false;
      if (restoreFocus) summary.focus({ preventScroll: true });
    };
    document.addEventListener("keydown", (event) => {
      if (event.key === "Escape" && more.open && !searchDialog?.open) {
        event.preventDefault();
        closeOptions(true);
      }
    });
    document.addEventListener("click", (event) => {
      if (!searchDialog?.open && !searchDialog?.contains(event.target) && !more.contains(event.target)) closeOptions();
    });
    more.addEventListener("focusout", (event) => {
      // Modal search restores this control when it closes.
      if (!searchDialog?.open && event.relatedTarget && !more.contains(event.relatedTarget)) closeOptions();
    });
    // Reserve the viewport-wide bar's real height, including safe-area padding.
    const placeToolbar = () => {
      let bottom = 0;
      document.querySelectorAll(".notebook-toolbar").forEach((bar) => {
        bar.style.bottom = `${bottom}px`;
        bottom += bar.getBoundingClientRect().height;
      });
      document.documentElement.style.setProperty("--notebook-bar-height", `${bottom}px`);
    };
    let toolbarFrame;
    const scheduleToolbar = () => {
      if (toolbarFrame) return;
      toolbarFrame = requestAnimationFrame(() => { toolbarFrame = 0; placeToolbar(); });
    };
    const toolbarResize = new ResizeObserver(scheduleToolbar);
    toolbarResize.observe(container);
    toolbarResize.observe(toolbar);
    window.addEventListener("resize", scheduleToolbar);
    const runAll = document.createElement("button");
    runAll.type = "button";
    runAll.textContent = "Run all";
    runAll.title = "Start a fresh notebook session and run every cell, retaining code edits";
    runAll.disabled = true;
    const restart = document.createElement("button");
    restart.type = "button";
    restart.textContent = "Restart original";
    restart.title = "Reload the page to restart Python and run the original example, discarding edits";
    restart.disabled = false;
    restart.addEventListener("click", () => {
      // A pinned reset can be invoked far below a lazily initialized notebook.
      // Restore its viewport before reload so the fresh loader is observed.
      container.scrollIntoView({ block: "start", behavior: "instant" });
      window.location.reload();
    });
    runAll.className = "notebook-run-all";
    toolbar.append(status, runAll, more);
    options.append(restart);
    const setupToggle = document.createElement("button");
    setupToggle.type = "button";
    setupToggle.textContent = "Show setup cells";
    setupToggle.setAttribute("aria-expanded", "false");
    setupToggle.disabled = true;
    setupToggle.addEventListener("click", () => {
      showSetup = !showSetup;
      notebook.querySelectorAll("[data-notebook-setup]").forEach((cell) => cell.hidden = !showSetup);
      setupToggle.textContent = showSetup ? "Hide setup cells" : "Show setup cells";
      setupToggle.setAttribute("aria-expanded", String(showSetup));
    });
    options.prepend(setupToggle);
    container.append(toolbar);
    const help = document.createElement("p");
    help.className = "notebook-help";
    help.textContent = "Edit source, then Run cell to update its output. Hide source keeps the result visible. Show setup cells reveals the runtime bootstrap. Run all starts a new session and keeps edits; Restart original discards them.";
    options.append(help);
    const updateCells = () => {
      notebook.querySelectorAll("marimo-island").forEach((island) => {
        const label = island.parentElement.querySelector(":scope > .notebook-cell-header > .notebook-cell-state");
        const state = island.dataset.status;
        const runtime = notebookIntegration?.cellRuntime(island.cellId);
        const text = ["running", "queued"].includes(state) ? state
          : island.querySelector(".text-error") ? "error"
          : island.hasAttribute("data-notebook-edited") ? "edited"
          : runtime?.stopped ? (runtime.runElapsedTimeMs === null ? "not run" : "stopped")
          : started ? "finished" : "waiting for Python";
        if (label && label.textContent !== text) label.textContent = text;
        if (label) label.hidden = false;
        const busy = ["running", "queued"].includes(state);
        const timing = island.parentElement.querySelector(":scope > .notebook-cell-header > .notebook-cell-time");
        const elapsed = runtime?.runElapsedTimeMs;
        if (timing) {
          const available = started && !busy && !["running", "queued"].includes(runtime?.status) &&
            Number.isFinite(elapsed) && elapsed >= 0;
          timing.hidden = !available;
          if (available) {
            const duration = elapsed < 1 ? "<1 ms" : elapsed < 1000
              ? `${Math.round(elapsed)} ms` : `${(elapsed / 1000).toFixed(2)} s`;
            if (timing.textContent !== duration) timing.textContent = duration;
            timing.dataset.elapsedMs = String(elapsed);
            timing.setAttribute("aria-label", `Last cell execution: ${duration}`);
          } else {
            timing.removeAttribute("data-elapsed-ms");
          }
        }
        island.parentElement.setAttribute("aria-busy", String(busy));
        const run = island.parentElement.querySelector(":scope > .notebook-cell-header > .notebook-run");
        const upstreamRun = island.querySelector(runSelector);
        const disabled = busy || !started || !upstreamRun || upstreamRun.disabled;
        if (run && run.disabled !== disabled) run.disabled = disabled;
        run?.setAttribute("aria-busy", String(busy));
        const sourceToggle = island.parentElement.querySelector(":scope > .notebook-cell-header > .notebook-source-toggle");
        const source = island.querySelector(sourceSelector)?.closest("marimo-ui-element");
        if (source && sourceToggle) {
          // Marimo replaces the payload's source host during hydration.
          source.id = sourceToggle.getAttribute("aria-controls");
          source.hidden = island.dataset.notebookSourceHidden === "true";
          const action = source.hidden ? "Show source" : "Hide source";
          sourceToggle.setAttribute("aria-label", action);
          sourceToggle.title = action;
          sourceToggle.setAttribute("aria-expanded", String(!source.hidden));
        }
        // Label upstream controls; execution remains Marimo's own click handler.
        for (const [icon, name] of [["play", "Run cell"], ["copy", "Copy code"]]) {
          island.querySelectorAll(`:scope > .relative > .absolute button:has(svg.lucide-${icon})`).forEach((button) => {
            button.setAttribute("aria-label", name);
            button.title = name;
            button.dataset.notebookAction = name;
          });
        }
      });
    };
    // Marimo mirrors actual cell execution to each island's data-status.
    // Initial idle cells precede Python startup, so wait for the first run.
    const updateNotebookState = (mutations = []) => {
      if (restarting) return;
      for (const mutation of mutations) {
        if (mutation.attributeName === "data-status" &&
            ["queued", "running"].includes(mutation.target.dataset.status)) {
          markEdited(mutation.target, currentCode(mutation.target));
        }
      }
      const runtimes = [...notebook.querySelectorAll("marimo-island")]
        .map((island) => notebookIntegration?.cellRuntime(island.cellId));
      const busy = runtimes.some((runtime) => ["queued", "running"].includes(runtime?.status)) ||
        Boolean(notebook.querySelector('marimo-island:is([data-status="queued"], [data-status="running"])'));
      const failed = Boolean(notebook.querySelector("marimo-island .text-error"));
      const edited = Boolean(notebook.querySelector("[data-notebook-edited]"));
      // React may render a completed cell directly as idle. Kernel measurements
      // establish execution even when no transient busy DOM state was observed.
      started ||= busy || runtimes.some((runtime) => Number.isFinite(runtime?.runElapsedTimeMs));
      updateCells();
      runAll.disabled = busy || !started;
      if (!started && !failed) return;
      status.dataset.state = busy ? "busy" : failed ? "error" : edited ? "stale" : "ready";
      const running = [...notebook.querySelectorAll('marimo-island[data-status="running"]')];
      let message = busy ? "Running notebook…" : failed ? "A cell failed. See its error below."
        : edited ? "Code changed. Run the cell to update its results." : "Ready — all queued cells finished.";
      const setupError = notebook.querySelector('[data-notebook-setup] .text-error');
      if (setupError && !busy) message = `Notebook setup failed: ${setupError.textContent}. Setup cells contain the error details. Restart original to retry.`;
      if (running.some((island) => island.dataset.notebookPhase === "install")) {
        message = "Downloading and installing Python packages…";
      } else if (running.some((island) => island.dataset.notebookPhase === "activate")) {
        message = "Activating Symbolica and loading its modules…";
      }
      if (status.textContent !== message) status.textContent = message;
      notebook.setAttribute("aria-busy", String(busy));
    };
    const executionObserver = new MutationObserver(updateNotebookState);
    executionObserver.observe(notebook, {
      subtree: true, childList: true, attributes: true,
      attributeFilter: ["data-status", "data-notebook-edited", "class", "disabled"],
    });
    document.addEventListener("notebook-execution-request", ({detail}) => {
      detail.cellIds.forEach((id, index) => {
        const island = [...notebook.querySelectorAll("marimo-island")]
          .find(island => island.cellId === id && island.dataset.appId === detail.appId);
        if (!island || typeof detail.codes[index] !== "string") return;
        lastRunCode.set(island, detail.codes[index]);
        markEdited(island, currentCode(island));
      });
    });
    notebook.addEventListener("notebook-source-change", (event) => {
      const island = event.target.closest("marimo-island");
      if (!island || typeof event.detail !== "string") return;
      editedCode.set(Number(island.dataset.cellIdx), event.detail);
      markEdited(island, event.detail);
    });
    // Marimo UI elements publish their complete value with this event; native
    // input contains only DOM edits, which can be virtualized by CodeMirror.
    notebook.addEventListener("marimo-value-input", (event) => {
      const editor = event.detail?.element;
      if (!editor?.matches("marimo-code-editor")) return;
      const island = editor.closest("marimo-island");
      if (island && island.querySelector(sourceSelector) === editor && typeof event.detail.value === "string") {
        // A debounced widget value can trail a newer CodeMirror transaction.
        editedCode.set(Number(island.dataset.cellIdx),
          notebookIntegration?.editorSource(editor.shadowRoot?.querySelector(".cm-editor")) ?? event.detail.value);
      }
    }, { capture: true });
    try {
      const url = new URL(`${docsRoot}assets/notebooks/${container.dataset.notebook}.json`, document.baseURI);
      // Payloads contain authored source and must revalidate across site updates.
      // The runtime worker separately caches only content-addressed binaries.
      const response = await fetch(url, { cache: "no-cache" });
      if (!response.ok) throw new Error("This documentation build does not include the browser notebook assets.");
      const payload = await response.json();
      // Islands run Python in a worker, whose URL can belong to the CDN.
      // Resolve our wheel against the version-local asset URL on this page.
      const wheel = new URL(payload.runtime_wheel || payload.wheel, url).href;
      let cells = payload.body.replaceAll("__NOTEBOOK_WHEEL_URL__", wheel);
      if (payload.dependency_wheel) {
        cells = cells.replaceAll("__DEPENDENCY_WHEEL_URL__", new URL(payload.dependency_wheel, url).href);
      }
      cells = cells.replaceAll("__SPENSO_STUB_URL__", new URL("spynso3.pyi", url).href);
      notebook.innerHTML = cells;
      const originalBody = cells;
      const replaceSession = async () => {
        restarting = true;
        runAll.disabled = true;
        status.dataset.state = "busy";
        status.textContent = "Restarting and running all cells…";
        try {
          container.dataset.notebookRuntime = notebookRuntimeURL;
      const runtime = await notebookRuntime;
          // Dispose widgets while their kernel session still accepts cleanup.
          notebook.replaceChildren();
          await new Promise(requestAnimationFrame);
          const template = document.createElement("template");
          template.innerHTML = originalBody;
          const script = template.content.querySelector('script[type="application/vnd.marimo.islands+json"]');
          const app = JSON.parse(script.textContent);
          for (const [index, code] of editedCode) app.cells[index].code = code;
          // A changed hidden bootstrap makes initialize replace the active
          // session even when the user's cells have no edits. Keeping the
          // bridge active lets late widget cleanup use its session transition.
          app.cells[0].code += `\n# Run all generation ${++runGeneration}`;
          script.textContent = JSON.stringify(app);
          started = false;
          notebook.replaceChildren(template.content);
          prepareCells();
          restarting = false;
          await runtime.initialize();
        } catch (error) {
          restarting = false;
          status.dataset.state = "error";
          status.textContent = `${error.message} Reload this page to retry.`;
          restart.disabled = false;
        }
      };
      runAll.addEventListener("click", () => void replaceSession());
      if (payload.full_notebook) {
        const full = document.createElement("a");
        full.href = new URL(payload.full_notebook, url).href;
        full.textContent = "Open full notebook";
        full.className = "notebook-full-link";
        const narrow = window.matchMedia("(max-width: 40rem)");
        const placeFullLink = () => {
          if (narrow.matches) options.insertBefore(full, help);
          else toolbar.insertBefore(full, more);
        };
        narrow.addEventListener("change", placeFullLink);
        placeFullLink();
        syncExplicitTheme();
      }
      const prepareCells = () => {
        // Capture hidden bootstrap code before Marimo consumes its source elements.
        notebook.querySelectorAll("marimo-island").forEach((island, index) => {
          const wrapper = document.createElement("div");
          wrapper.className = "notebook-cell";
          island.before(wrapper);
          const label = document.createElement("p");
          label.className = "notebook-cell-state";
          label.id = `notebook-state-${++notebookElementId}`;
          label.textContent = "waiting for Python";
          label.hidden = false;
          const header = document.createElement("div");
          header.className = "notebook-cell-header";
          header.append(label);
          const prose = payload.prose_cells?.includes(index);
          if (prose) {
            wrapper.dataset.notebookProse = "";
            // The surrounding guide already supplies the page title. Retain the
            // first prose cell's introduction and every later notebook heading.
            if (index === payload.prose_cells[0]) wrapper.dataset.notebookIntro = "";
          }
          else {
            const timing = document.createElement("span");
            timing.className = "notebook-cell-time";
            timing.title = "Last cell execution time, measured by the Python kernel; excludes queueing and startup";
            timing.hidden = true;
            header.append(timing);
            const run = document.createElement("button");
            run.type = "button";
            run.className = "notebook-run";
            run.setAttribute("aria-label", "Run cell");
            run.title = "Run cell";
            run.innerHTML = '<svg viewBox="0 0 20 20" aria-hidden="true"><path d="M6 3.5 16 10 6 16.5Z" fill="currentColor"/></svg>';
            run.setAttribute("aria-describedby", label.id);
            run.disabled = true;
            run.addEventListener("click", () => {
              island.querySelector(runSelector)?.click();
            });
            header.append(run);
          }
          wrapper.append(header, island);
          const code = decodeURIComponent(island.querySelector("marimo-cell-code")?.textContent || "")
            || JSON.parse(island.querySelector("marimo-code-editor")?.dataset.initialValue || '""');
          lastRunCode.set(island, editedCode.get(index) ?? code);
          if (code.includes("_micropip.install(")) island.dataset.notebookPhase = "install";
          else if (/^\s*(?:from|import) symbolica\b/m.test(code)) island.dataset.notebookPhase = "activate";
          // Collapse generated installation work as a whole. Source-only
          // toggles below always retain the authored cell's output.
          if (island.dataset.notebookPhase === "install" || payload.setup_cells?.includes(index)) {
            wrapper.dataset.notebookSetup = "";
            wrapper.id = `notebook-setup-${++notebookElementId}`;
            wrapper.hidden = !showSetup;
            label.hidden = false;
            const editor = island.querySelector(sourceSelector);
            if (editor) editor.dataset.disabled = "true";
          } else if (!prose) {
            const source = island.querySelector(sourceSelector)?.closest("marimo-ui-element");
            if (source && payload.hidden_source_cells?.includes(index)) {
              source.hidden = true;
              island.dataset.notebookSourceHidden = "true";
            }
            if (source) {
              source.id = `notebook-source-${++notebookElementId}`;
              const toggle = document.createElement("button");
              toggle.type = "button";
              toggle.className = "notebook-source-toggle";
              toggle.setAttribute("aria-label", source.hidden ? "Show source" : "Hide source");
              toggle.title = source.hidden ? "Show source" : "Hide source";
              toggle.innerHTML = '<svg viewBox="0 0 20 20" fill="none" stroke="currentColor" stroke-width="1.4" aria-hidden="true"><path d="m6 5-5 5 5 5M14 5l5 5-5 5M12 3l-4 14"/></svg>';
              toggle.setAttribute("aria-controls", source.id);
              toggle.setAttribute("aria-expanded", String(!source.hidden));
              toggle.addEventListener("click", () => {
                island.dataset.notebookSourceHidden = String(island.dataset.notebookSourceHidden !== "true");
                updateCells();
              });
              header.append(toggle);
              if (code.split("\n").length > 15) {
                const expand = document.createElement("button");
                expand.type = "button";
                expand.className = "notebook-source-expand";
                expand.textContent = "Expand code";
                expand.setAttribute("aria-expanded", "false");
                expand.setAttribute("aria-controls", source.id);
                expand.addEventListener("click", () => {
                  const expanded = expand.getAttribute("aria-expanded") !== "true";
                  island.querySelector(sourceSelector)?.toggleAttribute("data-notebook-expanded", expanded);
                  expand.setAttribute("aria-expanded", String(expanded));
                  expand.textContent = expanded ? "Collapse code" : "Expand code";
                });
                header.append(expand);
              }
            }
          }
        });
        const setup = [...notebook.querySelectorAll("[data-notebook-setup]")];
        setupToggle.disabled = setup.length === 0;
        setupToggle.setAttribute("aria-controls", setup.map((cell) => cell.id).join(" "));
        notebook.querySelectorAll("marimo-island > marimo-ui-element > marimo-code-editor").forEach((editor) => {
          editor.dataset.theme = JSON.stringify(root.dataset.theme === "dark" ? "dark" : "light");
          editor.dataset.showCopyButton = "true";
          if (container.dataset.notebook === "rendering_api") editor.dataset.maxHeight = "360";
          // Readable source while the editor's JavaScript downloads. The wrapper
          // is retained by the islands runtime; :defined hides this once mounted.
          const preview = document.createElement("pre");
          preview.className = "live-notebook-code-preview";
          const code = document.createElement("code");
          code.textContent = JSON.parse(editor.dataset.initialValue || '""');
          preview.append(code);
          editor.after(preview);
        });
      };
      prepareCells();
      status.textContent = "Downloading notebook runtime…";
      container.append(notebook);
      if (!notebookRuntime) {
        if (!(await sourceStylesReady).every(Boolean)) throw new Error("Source stylesheet failed to load");
        if (payload.cache_worker && "serviceWorker" in navigator) {
          let cacheDeadline;
          const cacheEvents = new AbortController();
          try {
            const cacheReady = (async () => {
              const workerURL = new URL(payload.cache_worker, url);
              const registration = await navigator.serviceWorker.register(workerURL, {updateViaCache: "none"});
              const installing = registration.installing || registration.waiting;
              if (installing && installing.state !== "activated") {
                await new Promise((resolve, reject) => installing.addEventListener("statechange", () => {
                  if (installing.state === "activated") resolve();
                  else if (installing.state === "redundant") reject(new Error("Cache worker update failed"));
                }, {signal: cacheEvents.signal}));
              }
              // Activation claims this page; controllerchange may already have
              // fired by the time the activation promise resolves.
            })();
            // Caching is optional. A restricted or stalled worker must not
            // prevent the scientific runtime from starting over the network.
            await Promise.race([cacheReady, new Promise(resolve => {
              cacheDeadline = setTimeout(resolve, 2000);
            })]);
          } catch {
            // Storage and worker restrictions affect caching, not correctness.
          } finally {
            clearTimeout(cacheDeadline);
            cacheEvents.abort();
          }
        }
        const head = new DOMParser().parseFromString(payload.head, "text/html");
        // Use the documentation fonts instead of downloading Marimo's unused
        // Fira Mono, Lora and PT Sans. Keep the runtime and math stylesheets.
        head.querySelectorAll("link").forEach((link) => {
          if (!["fonts.googleapis.com", "fonts.gstatic.com"].includes(new URL(link.href).hostname)) {
            if (link.title.startsWith('marimo')) link.title = 'marimo';
            document.head.append(link);
          }
        });
        const runtimeURL = payload.runtime ? new URL(payload.runtime, url).href : head.querySelector('script[type="module"]').src;
        notebookRuntimeURL = runtimeURL;
        observeNotebookWorker(runtimeURL);
        notebookRuntime = import(runtimeURL).catch((error) => {
          notebookRuntime = undefined;
          throw error;
        });
        void loadNotebookIntegration(runtimeURL);
      }
      container.dataset.notebookRuntime = notebookRuntimeURL;
      const runtime = await notebookRuntime;
      void notebookIntegrationPromise.then((integration) => {
        if (!integration || container.dataset.notebookState === "error") return;
        unsubscribeTimings = integration.subscribe(() => updateNotebookState());
        updateNotebookState();
      });
      addEventListener("pagehide", (event) => {
        if (!event.persisted) unsubscribeTimings?.();
      });
      const startupBegan = performance.now();
      const updateStartup = (phase) => {
        container.dataset.notebookStartup = phase;
        const elapsed = Math.floor((performance.now() - startupBegan) / 1000);
        if (!started) status.textContent = phase + (elapsed >= 10 ? ` (${elapsed} s since worker startup; waiting for this step.)` : "");
      };
      notebookStartup.listeners.add(updateStartup);
      updateStartup(notebookStartup.phase);
      try {
        await Promise.race([(async () => {
          await runtime.initialized;
          // A second notebook may have arrived after the module's initial scan.
          await runtime.initialize();
          await notebookStartup.executing;
        })(), notebookStartup.failure]);
      } finally { notebookStartup.listeners.delete(updateStartup); }
      container.dataset.notebookState = "active";
    } catch (error) {
      unsubscribeTimings?.();
      executionObserver.disconnect();
      notebook.setAttribute("aria-busy", "false");
      container.dataset.notebookState = "error";
      status.dataset.state = "error";
      status.textContent = `${error.message} Reload this page to retry.`;
      restart.disabled = false;
    }
  };
  const notebookObserver = new IntersectionObserver((entries, observer) => {
    for (const { target: container, isIntersecting } of entries) {
      if (!isIntersecting) continue;
      observer.unobserve(container);
      // Independent notebooks fetch together and share the same runtime import.
      void loadNotebook(container);
    }
  }, { rootMargin: "200px" });
  document.querySelectorAll("[data-notebook]").forEach((container) => notebookObserver.observe(container));

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

  let searchGeneration = 0;
  let activeSearchResult = -1;
  const setActiveSearchResult = (index, scroll = false) => {
    const links = [...(searchResults?.querySelectorAll("a") || [])];
    activeSearchResult = links.length ? (index + links.length) % links.length : -1;
    links.forEach((link, i) => link.setAttribute("aria-selected", String(i === activeSearchResult)));
    const active = links[activeSearchResult];
    if (active) {
      searchInput.setAttribute("aria-activedescendant", active.id);
      if (scroll) active.scrollIntoView({ block: "nearest" });
    } else searchInput?.removeAttribute("aria-activedescendant");
  };
  const searchKindLabel = (kind) => {
    const developer = String(kind).match(/^Developers · developer (current|proposal|archived|investigation)\b/);
    if (developer) return `Developer note · ${developer[1][0].toUpperCase()}${developer[1].slice(1)}`;
    const labels = {
      "python-api": "Python API", "rust-api": "Rust API", "typst api": "Typst API",
      "guides": "Guide", "start here": "Overview", "getting started": "Getting started",
      "version history": "Version history", "companion tools": "Companion tools",
      "settings namespace": "Settings namespace", "setting": "Setting", "command": "Command",
      "advanced": "Advanced", "product": "Overview", "reference": "Reference",
      "showcases": "Showcases", "tutorial": "Tutorial", "examples": "Examples",
      "dependency": "Dependency", "topology": "Topology", "collaboration": "Collaboration", "talk": "Talk",
    };
    return String(kind).split(" · ").map((part) => labels[part] || part).join(" · ");
  };
  const renderSearch = async (query) => {
    if (!searchResults) return;
    const generation = ++searchGeneration;
    const normalized = normalizeSearch(query.trim());
    searchResults.replaceChildren();
    searchStatus.textContent = "";
    setActiveSearchResult(-1);
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
    const index = await getIndex();
    if (generation !== searchGeneration) return;
    searchResults.replaceChildren();
    const seen = new Set();
    const matches = index
      .map((entry, index) => ({ entry, index, score: searchScore(entry, normalized, terms) }))
      .filter((match) => match.score !== null)
      .sort((left, right) => right.score - left.score ||
        left.entry.title.length - right.entry.title.length ||
        left.entry.title.localeCompare(right.entry.title) || left.index - right.index)
      .filter(({ entry }) => {
        // Shared Python modules have identical reference pages in multiple
        // product manuals. Keep one destination for the same documented object.
        const target = /python-api/.test(entry.kind)
          ? entry.href.replace(/^products\/[^/]+\/(snapshots\/[^/]+|[^/]+)\/reference\/python\//, "$1/python/")
          : entry.href;
        const key = `${target}\n${entry.title}`;
        if (seen.has(key)) return false;
        seen.add(key);
        return true;
      })
      .slice(0, 30)
      .map((match) => match.entry);
    for (const entry of matches) {
      const item = document.createElement("li");
      item.className = "search-result";
      item.setAttribute("role", "presentation");
      const link = document.createElement("a");
      link.href = `${searchRoot}${entry.href}`;
      link.id = `search-result-${generation}-${searchResults.children.length}`;
      link.setAttribute("role", "option");
      link.addEventListener("pointermove", () => setActiveSearchResult([...searchResults.querySelectorAll("a")].indexOf(link)));
      link.addEventListener("focus", () => setActiveSearchResult([...searchResults.querySelectorAll("a")].indexOf(link)));
      const title = document.createElement("span");
      title.className = "search-result-title";
      title.textContent = entry.title;
      const kind = document.createElement("span");
      kind.className = "search-result-kind";
      kind.textContent = searchKindLabel(entry.kind);
      const summary = document.createElement("span");
      summary.className = "search-result-summary";
      summary.textContent = entry.summary || searchKindLabel(entry.kind);
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
    searchStatus.textContent = matches.length ? `${matches.length} ${matches.length === 1 ? "result" : "results"}` : "No results";
    setActiveSearchResult(0);
  };

  let searchFocus;
  const openSearch = () => {
    if (!searchDialog || !searchInput || !searchResults) return false;
    navigation?.close();
    if (searchDialog.open) return true;
    searchFocus = document.activeElement;
    while (searchFocus?.shadowRoot?.activeElement) searchFocus = searchFocus.shadowRoot.activeElement;
    searchDialog.showModal();
    searchInput.setAttribute("aria-expanded", "true");
    searchInput.focus();
    searchInput.select();
    renderSearch(searchInput.value);
    return true;
  };
  searchDialog?.addEventListener("close", () => {
    searchInput?.setAttribute("aria-expanded", "false");
    if (searchFocus?.isConnected && !notebookIntegration?.focusEditor(searchFocus)) {
      searchFocus.focus({ preventScroll: true });
    }
  });
  // Capture before editor keymaps and native search-input Escape handling.
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && searchDialog?.open) {
      event.preventDefault();
      event.stopPropagation();
      searchDialog.close();
      return;
    }
    if (event.key.toLowerCase() === "k" && (event.metaKey || event.ctrlKey) && openSearch()) {
      event.preventDefault();
      event.stopPropagation();
    }
  }, { capture: true });
  searchButtons.forEach((button) => button.addEventListener("click", openSearch));
  if (searchInput && searchResults) {
    searchResults.id ||= "documentation-search-results";
    searchResults.setAttribute("role", "listbox");
    searchResults.setAttribute("aria-label", "Documentation results");
    searchInput.setAttribute("role", "combobox");
    searchInput.setAttribute("aria-autocomplete", "list");
    searchInput.setAttribute("aria-controls", searchResults.id);
    searchInput.setAttribute("aria-expanded", "false");
    searchInput.addEventListener("keydown", (event) => {
      if (event.key === "ArrowDown" || event.key === "ArrowUp") {
        event.preventDefault();
        setActiveSearchResult(activeSearchResult + (event.key === "ArrowDown" ? 1 : -1), true);
      }
    });
    searchDialog.querySelector("form")?.addEventListener("submit", (event) => {
      event.preventDefault();
      searchResults.querySelectorAll("a")[activeSearchResult]?.click();
    });
    searchDialog.querySelectorAll('button[value="close"], [data-search-close]').forEach((button) => {
      button.type = "button";
      button.addEventListener("click", () => searchDialog.close());
    });
  }
  searchInput?.addEventListener("input", (event) => renderSearch(event.target.value));
  searchDialog?.addEventListener("click", (event) => {
    if (event.target === searchDialog) searchDialog.close();
  });
  document.addEventListener("keydown", (event) => {
    const shortcut = !event.composedPath().some((target) => target.isContentEditable) && (
      event.key === "/" && !event.metaKey && !event.ctrlKey && !event.altKey &&
      !/input|textarea|select/i.test(document.activeElement?.tagName)
    );
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
    const restorePublicationFilters = () => {
      const params = new URLSearchParams(location.search);
      search.value = params.get("q") || "";
      author.value = params.get("author") || "";
      year.value = params.get("year") || "";
      type.value = params.get("type") || "";
      sort.value = params.get("sort") || "newest";
    };
    restorePublicationFilters();

    const updatePublications = (historyMode = "replace") => {
      const query = search.value.trim().toLowerCase();
      const matches = cards.filter((card) => {
        const people = card.dataset.people.split(/\s+/);
        const types = card.dataset.types.split("|");
        return (!query || `${card.textContent} ${card.dataset.arxiv}`.toLowerCase().includes(query)) &&
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
      const url = `${location.pathname}${next.size ? `?${next}` : ""}${location.hash}`;
      if (historyMode && url !== `${location.pathname}${location.search}${location.hash}`) {
        history[`${historyMode}State`](null, "", url);
      }
    };
    publicationFilters.addEventListener("input", (event) => {
      if (event.target === search) updatePublications();
    });
    publicationFilters.addEventListener("change", () => updatePublications("push"));
    addEventListener("popstate", () => { restorePublicationFilters(); updatePublications(false); });
    publicationFilters.addEventListener("submit", (event) => event.preventDefault());
    publicationFilters.addEventListener("reset", () => {
      // Native values reset after this event returns.
      setTimeout(() => { updatePublications("push"); search.focus(); }, 0);
    });
    updatePublications();
  }

  const talkFilters = document.querySelector("[data-talk-filters]");
  if (talkFilters) {
    const search = talkFilters.querySelector("[data-talk-search]");
    const year = talkFilters.querySelector("[data-talk-year]");
    const count = talkFilters.querySelector("[data-talk-count]");
    const cards = [...document.querySelectorAll("[data-talk]")];
    const restoreTalkFilters = () => {
      const params = new URLSearchParams(location.search);
      search.value = params.get("q") || "";
      year.value = params.get("year") || "";
    };
    restoreTalkFilters();
    const updateTalks = (historyMode = "replace") => {
      const query = search.value.trim().toLowerCase();
      let matches = 0;
      for (const card of cards) {
        card.hidden = Boolean((query && !card.textContent.toLowerCase().includes(query)) ||
          (year.value && card.dataset.year !== year.value));
        if (!card.hidden) matches++;
      }
      document.querySelectorAll(".talk-year").forEach((section) => {
        section.hidden = !section.querySelector("[data-talk]:not([hidden])");
      });
      count.textContent = `${matches} talk${matches === 1 ? "" : "s"}`;
      const next = new URLSearchParams();
      if (query) next.set("q", search.value.trim());
      if (year.value) next.set("year", year.value);
      const url = `${location.pathname}${next.size ? `?${next}` : ""}${location.hash}`;
      if (historyMode && url !== `${location.pathname}${location.search}${location.hash}`) {
        history[`${historyMode}State`](null, "", url);
      }
    };
    talkFilters.addEventListener("input", (event) => {
      if (event.target === search) updateTalks();
    });
    talkFilters.addEventListener("change", () => updateTalks("push"));
    addEventListener("popstate", () => { restoreTalkFilters(); updateTalks(false); });
    talkFilters.addEventListener("submit", (event) => event.preventDefault());
    talkFilters.addEventListener("reset", () => {
      setTimeout(() => { updateTalks("push"); search.focus(); }, 0);
    });
    updateTalks();
  }

  const pythonTokens = /(?<comment>#[^\n]*)|(?<string>[rRuUbBfF]{0,2}(?:"""[\s\S]*?"""|'''[\s\S]*?'''|"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'))|(?<keyword>\b(?:False|None|True|and|as|assert|async|await|break|case|class|continue|def|del|elif|else|except|finally|for|from|global|if|import|in|is|lambda|match|nonlocal|not|or|pass|raise|return|try|while|with|yield)\b)|(?<number>\b(?:0[xX][\dA-Fa-f](?:_?[\dA-Fa-f])*|\d(?:_?\d)*(?:\.\d(?:_?\d)*)?(?:[eE][+-]?\d(?:_?\d)*)?j?)\b)|(?<call>\b[A-Za-z_]\w*(?=\s*\())|(?<operator>->|:=|==|!=|<=|>=|\/\/|\*\*|<<|>>|[-+*/%@&|^~<>:=])/g;
  const syntaxSelector = ".syntax-call, .syntax-comment, .syntax-keyword, .syntax-number, .syntax-operator, .syntax-string";
  document.querySelectorAll("code[data-lang='python']").forEach((code) => {
    const walker = document.createTreeWalker(code, NodeFilter.SHOW_TEXT, {
      acceptNode: (node) => node.parentElement?.closest(syntaxSelector)
        ? NodeFilter.FILTER_REJECT
        : NodeFilter.FILTER_ACCEPT,
    });
    const textNodes = [];
    while (walker.nextNode()) textNodes.push(walker.currentNode);
    textNodes.forEach((node) => {
      const source = node.data;
      const highlighted = document.createDocumentFragment();
      let cursor = 0;
      let found = false;
      pythonTokens.lastIndex = 0;
      for (const match of source.matchAll(pythonTokens)) {
        found = true;
        highlighted.append(source.slice(cursor, match.index));
        const token = document.createElement("span");
        token.className = `syntax-${Object.keys(match.groups).find((name) => match.groups[name])}`;
        token.textContent = match[0];
        highlighted.append(token);
        cursor = match.index + match[0].length;
      }
      if (!found) return;
      highlighted.append(source.slice(cursor));
      node.replaceWith(highlighted);
    });
  });

  document.querySelectorAll(".reference-table-wrap").forEach((region) => {
    const hint = document.createElement("p");
    hint.className = "reference-table-scroll-hint";
    hint.textContent = "Scroll table horizontally →";
    region.before(hint);
    const update = () => {
      const overflowing = region.scrollWidth > region.clientWidth + 1;
      hint.hidden = !overflowing;
      if (overflowing) {
        region.tabIndex = 0;
        region.setAttribute("role", "region");
        region.setAttribute("aria-label", "Reference table; scroll horizontally for more columns");
      } else {
        region.removeAttribute("tabindex");
        region.removeAttribute("role");
        region.removeAttribute("aria-label");
      }
    };
    const observer = new ResizeObserver(update);
    observer.observe(region);
    if (region.firstElementChild) observer.observe(region.firstElementChild);
    update();
  });

  const copyButtons = new Map([...document.querySelectorAll("[data-copy-target]")]
    .map((button) => [button.dataset.copyTarget, button]));
  let copySerial = 0;
  document.querySelectorAll(".docs-article pre, .citation-card pre").forEach((pre) => {
    let button = pre.id ? copyButtons.get(pre.id) : undefined;
    let row = pre.parentElement?.classList.contains("reference-copy-row")
      ? pre.parentElement
      : undefined;
    if (!button) {
      if (!pre.id) {
        let target;
        do {
          target = `copy-block-${++copySerial}`;
        } while (document.getElementById(target));
        pre.id = target;
      }
      button = document.createElement("button");
      button.type = "button";
      button.dataset.copyTarget = pre.id;
      button.textContent = pre.closest(".api-doc-examples")
        ? "Copy example"
        : pre.matches(".api-signature, .api-member-signature")
          ? "Copy signature"
          : "Copy code";
      copyButtons.set(pre.id, button);
    }
    if (!row) {
      row = document.createElement("div");
      row.className = "reference-copy-row api-example-copy-row";
      pre.before(row);
      row.append(pre);
    }
    if (button.parentElement !== row) row.append(button);
  });

  const writeClipboard = async (value) => {
    if (navigator.clipboard?.writeText) {
      try {
        await navigator.clipboard.writeText(value);
        return;
      } catch {
        // Fall through for browsers that expose Clipboard without granting it.
      }
    }
    const input = document.createElement("textarea");
    input.value = value;
    input.readOnly = true;
    input.setAttribute("aria-hidden", "true");
    Object.assign(input.style, { position: "fixed", inset: "0 auto auto -100vw", opacity: "0" });
    body.append(input);
    let copied = false;
    try {
      input.select();
      copied = document.execCommand("copy");
    } finally {
      input.remove();
    }
    if (!copied) throw new Error("clipboard copy unavailable");
  };

  document.querySelectorAll("[data-copy-target]").forEach((button) => {
    const label = button.textContent.trim() || "Copy code";
    button.setAttribute("aria-label", label);
    button.title = label;
    const status = document.createElement("span");
    status.className = "copy-status";
    status.setAttribute("role", "status");
    button.parentElement?.append(status);
    let resetTimer;
    const resetCopyState = () => {
      delete button.dataset.copyState;
      delete button.dataset.copyFeedback;
      button.setAttribute("aria-label", label);
      button.title = label;
      status.textContent = "";
    };
    button.addEventListener("click", async () => {
      const target = document.getElementById(button.dataset.copyTarget);
      if (!target) return;
      clearTimeout(resetTimer);
      try {
        await writeClipboard(target.textContent.trim());
        button.dataset.copyState = "copied";
        button.dataset.copyFeedback = "Copied";
        button.setAttribute("aria-label", "Copied");
        button.title = "Copied";
        status.textContent = `${label.replace(/^Copy /, "")} copied.`;
      } catch {
        button.dataset.copyState = "error";
        button.dataset.copyFeedback = "Copy unavailable";
        button.setAttribute("aria-label", "Copy unavailable");
        button.title = "Copy unavailable";
        status.textContent = "Copy unavailable. Select the code and copy it manually.";
      }
      button.focus({ preventScroll: true });
      resetTimer = setTimeout(resetCopyState, 1800);
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
  if (document.querySelector('.docs-article pre code[data-lang], .live-notebook')) {
    const script = document.querySelector('script[src*="/assets/site.js"]');
    if (script) {
      loadDocumentationConfiguration().then(async ({url, value}) => {
        if (value.cacheWorker && 'serviceWorker' in navigator) {
          navigator.serviceWorker.register(new URL(value.cacheWorker, url), {updateViaCache: 'none'}).catch(console.warn);
        }
        const module = await import(new URL(value.entry, url));
        await module.initialize(url);
      }).catch(error => console.error('Documentation cells:', error));
    }
  }
})();
