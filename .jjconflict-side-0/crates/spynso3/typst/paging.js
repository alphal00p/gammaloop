export default {
  render({ model, el }) {
    const root = el.shadowRoot || el.attachShadow({ mode: "open" });
    root.innerHTML = `<style>
      :host { display:block; color:var(--foreground,var(--jp-ui-font-color1,#202624)); font:13px/1.5 system-ui,sans-serif; }
      * { box-sizing:border-box; }
      .toolbar,.parts { display:flex; flex-wrap:wrap; align-items:center; gap:8px; }
      .toolbar { padding:8px 0; }
      button,select { font:inherit; color:inherit; background:var(--background,var(--jp-layout-color1,#fff)); border:1px solid var(--border,var(--jp-border-color2,#dce2df)); border-radius:6px; min-height:30px; padding:4px 10px; transition:background .12s,border-color .12s; }
      button { cursor:pointer; font-weight:500; display:inline-flex; align-items:center; justify-content:center; }
      button.page-arrow { width:30px; padding:5px; }
      button:hover:not(:disabled),select:hover:not(:disabled) { background:var(--muted,var(--jp-layout-color2,#f3f5f4)); }
      button:disabled,select:disabled { opacity:.4; cursor:default; }
      button:focus-visible,select:focus-visible,input:focus-visible { outline:2px solid var(--primary,var(--jp-brand-color1,#16846b)); outline-offset:2px; }
      label { display:inline-flex; align-items:center; gap:7px; color:var(--muted-foreground,var(--jp-ui-font-color2,#68716c)); }
      .size { margin-left:4px; }
      .scroll-option { margin-left:4px; cursor:pointer; white-space:nowrap; min-height:30px; padding:4px 9px; border:1px solid var(--border,var(--jp-border-color2,#dce2df)); border-radius:6px; }
      .scroll-option:hover { background:var(--muted,var(--jp-layout-color2,#f3f5f4)); }
      input { accent-color:var(--primary,var(--jp-brand-color1,#16846b)); margin:0; }
      .status { color:var(--muted-foreground,var(--jp-ui-font-color2,#68716c)); font-size:12px; padding:0 0 8px; }
      .math { overflow-x:auto; border-block:1px solid var(--border,var(--jp-border-color2,#e2e7e4)); padding:12px 0; }
      .math [data-spenso-math] { overflow:visible; padding:0; }
      .spenso-page-heading { font:12px/1.5 system-ui,sans-serif; color:var(--muted-foreground,var(--jp-ui-font-color2,#68716c)); margin:8px 0; }
      .spenso-page-heading:first-of-type { margin-top:0; }
      .math.horizontal [data-spenso-sum] { flex-wrap:nowrap!important; max-width:none!important; }
      .parts { margin-top:8px; }
      .parts:empty { display:none; }
      .parts button { font-size:12px; }
    </style><nav class="toolbar" aria-label="Expression pages"></nav><div class="status" role="status" aria-live="polite"></div><div class="math"></div><nav class="parts" aria-label="Subexpressions"></nav>`;
    const nav = root.querySelector(".toolbar");
    const status = root.querySelector(".status");
    const math = root.querySelector(".math");
    const parts = root.querySelector(".parts");
    let serial = 0, pending = false, live = false, horizontal = true, timer;
    const view = crypto.randomUUID?.() || `${Date.now()}-${Math.random()}`;

    function send(action, value) {
      if (pending || !live) return;
      pending = true;
      serial++;
      root.querySelectorAll("button,select").forEach(e => e.disabled = true);
      status.textContent = "Rendering page…";
      model.send({ action, value, request: `${view}:${serial}` });
      timer = setTimeout(() => {
        pending = false;
        status.textContent = "No response from the kernel. Re-run the cell to reconnect.";
      }, 12000);
    }

    function button(parent, label, action, value, disabled = false) {
      const b = document.createElement("button");
      b.textContent = label;
      if (action === "previous" || action === "next") {
        b.className = "page-arrow";
        b.title = label;
        b.setAttribute("aria-label", label);
        const path = action === "previous" ? "M15 18l-6-6 6-6" : "M9 6l6 6-6 6";
        b.innerHTML = `<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="${path}"/></svg>`;
      }
      b.disabled = disabled || !live;
      b.onclick = () => send(action, value);
      parent.append(b);
    }

    function update() {
      const p = model.get("page");
      if (p.connected === view) live = true;
      if (pending && p.request !== `${view}:${serial}`) return;
      if (typeof p.request === "string" && p.request.startsWith(view + ":") && Number(p.request.split(":").at(-1)) < serial) return;
      clearTimeout(timer);
      pending = false;
      nav.replaceChildren();
      parts.replaceChildren();
      if (p.breadcrumbs.length) button(nav, "← Parent", "back");
      button(nav, "Previous", "previous", null, !p.previous);
      button(nav, "Next", "next", null, p.next == null);
      const label = document.createElement("label");
      label.className = "size";
      const select = document.createElement("select");
      select.setAttribute("aria-label", "Maximum terms per page");
      select.title = "Maximum terms per page";
      for (const n of [25, 100, 250, 500]) {
        const o = document.createElement("option");
        o.value = n;
        o.textContent = `${n} terms`;
        o.selected = n === p.page_size;
        select.append(o);
      }
      select.disabled = !live;
      select.onchange = () => send("size", Number(select.value));
      label.append(select);
      nav.append(label);
      const scrollLabel = document.createElement("label");
      scrollLabel.className = "scroll-option";
      const checkbox = document.createElement("input");
      checkbox.type = "checkbox";
      checkbox.checked = horizontal;
      math.classList.toggle("horizontal", horizontal);
      checkbox.onchange = () => {
        horizontal = checkbox.checked;
        math.classList.toggle("horizontal", horizontal);
      };
      scrollLabel.append(checkbox, "Horizontal scroll");
      nav.append(scrollLabel);
      const location = p.breadcrumbs.length ? p.breadcrumbs.join(" › ") + " · " : "";
      const range = p.total ? `${p.unit} ${p.start + 1}–${p.end} of ${p.total}` : "Expression preview";
      const budget = p.budget ? ` · Page reduced by ${p.budget} limit` : "";
      status.textContent = (!live ? "Navigation requires a live notebook kernel. " : "") + (p.error || location + range + budget);
      // Browsers do not register @font-face declarations inside a shadow root.
      // Reuse the standard notebook math font once per document; all other
      // viewer styles stay scoped to this shadow root.
      if (!document.getElementById("spenso-notebook-math-font")) {
        const font = p.html.match(/@font-face\s*\{[^}]*\}/)?.[0];
        if (font) {
          const style = document.createElement("style");
          style.id = "spenso-notebook-math-font";
          style.textContent = font;
          document.head.append(style);
        }
      }
      math.innerHTML = p.html;
      math.scrollTop = 0;
      math.scrollLeft = 0;
      for (const [id, text] of p.holes) button(parts, text, "open", id);
    }

    model.on("change:page", update);
    update();
    model.send({ action: "attach", view });
    return () => {
      model.send({ action: "dispose", view });
      clearTimeout(timer);
      model.off("change:page", update);
      root.replaceChildren();
    };
  }
};
