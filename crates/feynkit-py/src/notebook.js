/* Shared by rich displays whose content, rather than a single graph, sizes the
   notebook frame. Ancestor themes also work in accessible Marimo iframes. */
(() => {
  for (const root of document.querySelectorAll("[data-feynkit-notebook]")) {
    if (root.dataset.notebookReady) continue;
    root.dataset.notebookReady = "true";
    root.dataset.linnetFrameOwner = "true";
    const ancestors = [];
    let ancestor = root.parentElement,
      frameObserver;
    while (ancestor) {
      ancestors.push(ancestor);
      try {
        ancestor =
          ancestor.parentElement ||
          ancestor.ownerDocument.defaultView.frameElement;
      } catch (_) {
        ancestor = null;
      }
    }
    let previousTheme;
    const syncTheme = () => {
      if (!root.isConnected) {
        observer.disconnect();
        frameObserver?.disconnect();
        return;
      }
      let theme = "";
      for (const element of ancestors) {
        const explicit = element.dataset.theme,
          jupyter = element.getAttribute("data-jp-theme-light");
        if (explicit === "dark" || explicit === "light") {
          theme = explicit;
          break;
        }
        if (jupyter !== null) {
          theme = jupyter === "true" ? "light" : "dark";
          break;
        }
        if (element.classList.contains("dark")) {
          theme = "dark";
          break;
        }
        if (element.classList.contains("light")) {
          theme = "light";
          break;
        }
      }
      try {
        const frame = window.frameElement;
        if (frame?.parentElement) {
          const style = frame.ownerDocument.defaultView.getComputedStyle(
            frame.parentElement,
          );
          document.body.style.color = style.color;
          document.body.style.fontFamily = style.fontFamily;
          document.documentElement.style.colorScheme = style.colorScheme;
        }
      } catch (_) {
        /* Cross-origin hosts supply their own frame theme. */
      }
      if (theme === previousTheme) return;
      previousTheme = theme;
      root.dataset.theme = theme;
      root.style.colorScheme = theme || "light dark";
      root.dispatchEvent(new Event("feynkit-theme"));
    };
    const observer = new MutationObserver(syncTheme);
    ancestors.forEach((element) =>
      observer.observe(element, {
        attributes: true,
        childList: true,
        attributeFilter: ["class", "data-theme", "data-jp-theme-light"],
      }),
    );
    syncTheme();
    try {
      const frame = window.frameElement;
      if (frame && document.body) {
        document.documentElement.style.overflow = "hidden";
        // Measure content instead of scrollHeight, which includes the old
        // viewport height and prevents a collapsed notebook frame from shrinking.
        const fitFrame = () => {
          if (!root.isConnected) {
            frameObserver.disconnect();
            return;
          }
          const style = getComputedStyle(document.body);
          const bottom =
            Math.max(
              ...[...document.body.children].map(
                (element) => element.getBoundingClientRect().bottom,
              ),
            ) + window.scrollY;
          frame.style.height =
            Math.ceil(
              bottom +
                parseFloat(style.paddingBottom) +
                parseFloat(style.marginBottom) +
                6,
            ) + "px";
        };
        frameObserver = new ResizeObserver(() =>
          requestAnimationFrame(fitFrame),
        );
        frameObserver.observe(root);
      }
    } catch (_) {
      /* Cross-origin embeddings control their own frame sizing. */
    }
  }
})();
