(() => {
  for (const root of document.querySelectorAll(".feynkit-ltd-result")) {
    if (root.dataset.ready) continue;
    root.dataset.ready = "true";
    const data = JSON.parse(
      root.querySelector("[data-ltd]").textContent,
      (key, value) => key.endsWith("_html") ? root.expressionDisplay.math(value) : value,
    );
    const native = root
      .querySelector("[data-drawing]")
      .content.querySelector("svg");
    data.graph = {
      paths: [...native.querySelectorAll("[data-linnet-carrier]")].map((a) => ({
        ...JSON.parse(a.dataset.linnetDetail),
        d: a.querySelector("path").getAttribute("d"),
      })),
      nodes: [...native.querySelectorAll("[data-linnet-node]")].map((a) => ({
        id: +a.dataset.linnetNode,
        x: +a.getAttribute("cx"),
        y: +a.getAttribute("cy"),
      })),
    };
    const momentumArrows = new Map(
      [...native.querySelectorAll("[data-linnet-momentum]")].map((arrow) => [
        JSON.parse(arrow.dataset.linnetDetail).edge,
        arrow,
      ]),
    );
    const labelPositions = new Map(
      [...native.querySelectorAll('a[data-linnet-kind="edge"][transform]')]
        // Sampled edge hit regions also have transforms, but carry a half-edge
        // identity. Only the label anchor records the searched label position.
        .filter((label) => !("half-edge" in JSON.parse(label.dataset.linnetDetail)))
        .map((label) => {
          const rect = label.querySelector("rect");
          return [
            JSON.parse(label.dataset.linnetDetail).edge,
            `${label.getAttribute("transform")} translate(${+rect.getAttribute("width") / 2} ${+rect.getAttribute("height") / 2})`,
          ];
        }),
    );
    const navigation = native.querySelector("script").textContent;
    const physicalScale =
      native.width.baseVal.value / native.viewBox.baseVal.width;
    const instance = "ltd-" + Math.random().toString(36).slice(2);
    const ns = "http://www.w3.org/2000/svg";
    const make = (tag, attrs, parent) => {
      const el = document.createElementNS(ns, tag);
      for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
      parent.append(el);
      return el;
    };
    const sub = (s, i) => `${s}<sub>${i}</sub>`;
    const surfaceName = (s) => sub(s.kind === "E" ? "η" : "H", s.id);
    const sign = (c) => (String(c).startsWith("-") ? -1 : 1);
    function sumTerms(terms) {
      return (
        terms
          .map(
            ([c, t], i) =>
              `${sign(c) < 0 ? " − " : i ? " + " : ""}${String(c).replace("-", "") === "1" && t ? "" : String(c).replace("-", "")}${t}`,
          )
          .join("")
          .trim() || "0"
      );
    }
    function linear(expr, g, multiplier = 1) {
      const terms = [
        ...expr.internal.map(([i, c]) => [c, g.energy_html[g.edges[i]]]),
        ...expr.external.map(([i, c]) => [c, g.external_energy_html[i]]),
      ];
      if (expr.constant !== "0") terms.push([expr.constant, ""]);
      return sumTerms(
        terms.map(([c, t]) => [
          multiplier < 0 ? (sign(c) < 0 ? String(c).slice(1) : `-${c}`) : c,
          t,
        ]),
      );
    }
    function routing(edge, g) {
      const s = g.routing[edge].signature;
      return sumTerms([
        ...s.loop_signature.flatMap((c, i) =>
          c ? [[String(c), g.loop_energy_html[i]]] : [],
        ),
        ...s.external_signature.flatMap((c, i) =>
          c ? [[String(c), g.external_energy_html[i]]] : [],
        ),
      ]);
    }
    function selected(st) {
      const { g, r } = current(st),
        pair = r.propagator_factors[st.edge];
      return {
        g,
        r,
        pair,
        s: g.surfaces[st.surface ?? pair?.[st.branch < 0 ? 0 : 1]],
      };
    }
    function treeEdges(g, r) {
      return g.routing.map((_, i) => i).filter((i) => !r.cuts.includes(i));
    }
    function partition(g, r, edge) {
      // The inspected q enters this component: its boundary equation is A_e - q_e.
      if (!g.routing[edge]) return new Set();
      const a = new Set([g.routing[edge].head]);
      const tree = treeEdges(g, r).filter((i) => i !== edge);
      let changed = true;
      while (changed) {
        changed = false;
        for (const i of tree) {
          const { tail, head } = g.routing[i];
          if (a.has(tail) !== a.has(head)) {
            a.add(tail);
            a.add(head);
            changed = true;
          }
        }
      }
      return a;
    }
    function poleSign(st, index) {
      const { r } = current(st);
      if (r.cuts.includes(index)) return r.signs[index] === "Reversed" ? -1 : 1;
      return index === st.edge ? -st.branch : 0;
    }
    function cutCycle(st) {
      const { g, r } = current(st),
        cut = st.traceCut;
      if (!r.cuts.includes(cut)) return null;
      const e = g.routing[cut],
        queue = [{ node: e.head, path: [] }],
        seen = new Set([e.head]);
      let path;
      while (queue.length) {
        const p = queue.shift();
        if (p.node === e.tail) {
          path = p.path;
          break;
        }
        for (const i of treeEdges(g, r)) {
          const t = g.routing[i],
            next =
              t.tail === p.node ? t.head : t.head === p.node ? t.tail : null;
          if (next === null || seen.has(next)) continue;
          seen.add(next);
          queue.push({
            node: next,
            path: [...p.path, [i, t.tail === p.node ? 1 : -1]],
          });
        }
      }
      if (!path) return null;
      const directions = new Map([[cut, 1], ...path]);
      // Use the diagram's ordered reference basis, not a guessed momentum signature.
      const anchors = g.reference_chords;
      const latest = anchors.findLast((i) => directions.has(i)),
        orientation = directions.get(latest);
      if (!orientation || orientation !== poleSign(st, cut)) return null;
      for (const [i, d] of directions) directions.set(i, d * orientation);
      return { directions, latest, loop: anchors.indexOf(latest) };
    }
    function factor(s, st, edge, power = 1) {
      return `<button type="button" class="lp-surface hs-factored-factor" style="--hs-color:var(--lp-${s.kind.toLowerCase()})" data-surface="${s.id}" data-factor-id="surface-${s.id}" data-factor-power="${power}" data-factor-edge="${edge}" data-hover-edge="${edge}" data-kind="${s.kind}" aria-label="Inspect ${s.kind} surface ${s.id} of edge ${selected(st).g.edges[edge] ?? "—"}" aria-pressed="${s.id === selected(st).s?.id && edge === st.edge}">${surfaceName(s)}${power > 1 ? `<sup data-power-label>${power}</sup>` : ""}<span class="hs-swatch"></span></button>`;
    }
    const factorTrees = new Map();
    function equation(st) {
      const { g, r } = selected(st);
      if (!factorTrees.has(r.id)) {
        const entries = r.terms
          .flatMap((term) =>
            term.chains.map((chain) => ({
              coefficient: term.coefficient,
              factors: [
                ...term.energies.map((i) => `energy-${i}`),
                ...chain.map((i) => `surface-${i}`),
              ],
            })),
          )
          .map((entry, id) => ({ ...entry, id }));
        factorTrees.set(r.id, {
          entries,
          tree: entries.length
            ? root.expressionDisplay.factorTerms(entries)
            : null,
        });
      }
      const { entries, tree } = factorTrees.get(r.id);
      function denominator(factors) {
        const remaining = [...factors],
          parts = [];
        // Keep complete propagator pairs together after extracting common factors.
        for (const [key, power] of root.expressionDisplay.powers(factors.filter((key) => key.startsWith("energy-")))) {
          const edge = +key.slice(7);
          for (let i = 0; i < power; i++) remaining.splice(remaining.indexOf(key), 1);
          parts.push(
            root.expressionDisplay.energyFactor(g.energy_html[g.edges[edge]], edge, power),
          );
        }
        const counts = new Map(root.expressionDisplay.powers(remaining));
        for (const [edge, pair] of Object.entries(r.propagator_factors)) {
          const ids = [...new Set(pair)].filter(id => counts.has(`surface-${id}`));
          if (!ids.length) continue;
          parts.push(`<span class="lp-factor-pair" data-factor-pair="${edge}" data-active="${+edge === st.edge}" aria-label="Surface factors for edge ${g.edges[edge]}">${ids.map(id => {
            const power = counts.get(`surface-${id}`);
            counts.delete(`surface-${id}`);
            return factor(g.surfaces[id], st, +edge, power);
          }).join("")}</span>`);
        }
        for (const [key, power] of counts) parts.push(factor(g.surfaces[+key.slice(8)], st, null, power));
        return parts.join("");
      }
      function math(node) {
        const leaf = !node.children.length,
          coefficient = leaf ? entries[node.terms[0]].coefficient : "1",
          numerator = `<span data-coefficient="${coefficient}">${coefficient.replace("-", "−")}</span>`,
          fraction = node.factors.length
            ? root.expressionDisplay.fraction(numerator, denominator(node.factors), `data-path-active="true" data-shared-terms="${node.terms.length}"`)
            : "";
        if (leaf)
          return `<span class="hs-factor-leaf" data-ltd-term="${node.terms[0]}">${fraction || numerator}</span>`;
        const branches = node.children
          .map(
            (child, index) =>
              `<span class="hs-factor-summand">${index ? '<span aria-hidden="true">+</span>' : ""}${math(child)}</span>`,
          )
          .join("");
        const bracket = `<span class="hs-factored-bracket"><span class="hs-bracket-edge" aria-hidden="true"></span><span class="hs-factored-sum">${branches}</span><span class="hs-bracket-edge" aria-hidden="true"></span></span>`;
        return `<span class="hs-factor-node">${fraction}${fraction ? '<span aria-hidden="true">·</span>' : ""}${bracket}</span>`;
      }
      return `<div class="lp-equation hs-factored-equation" aria-label="Factored expression for residue ${r.id}"><span class="hs-factored-lhs">${sub("R", r.id)} =</span>${tree ? math(tree) : "0"}</div>`;
    }
    function globalRouting(st) {
      const { g } = selected(st);
      if (!g.loops) return "";
      const anchors = g.reference_chords
        .filter((i) => i !== null)
        .map(
          (i) =>
            `<span class="lp-math">${g.momentum_html[g.edges[i]]} = ${routing(i, g)}</span>`,
        );
      return `<div class="lp-routing"><div class="lp-caption">Global routing · integrate ${g.loop_energy_html.join(" → ")} · close below</div><div class="lp-routing-content">${anchors.join("")}</div></div>`;
    }
    function cutStructure(st) {
      const { g, r } = selected(st);
      return `<div class="lp-cut-structure" aria-live="polite"><span class="lp-caption">Cut structure</span><span class="lp-math">${sub("Σ", r.id)} =</span><div class="lp-cut-vector" style="--cut-edges:${g.edges.length}" role="group" aria-label="Signed cut vector for tree ${r.id}">${g.edges
        .map((e, i) => {
          const sigma = r.cuts.includes(i)
            ? r.signs[i] === "Reversed"
              ? -1
              : 1
            : 0;
          const value = sigma === 0 ? "0" : sigma < 0 ? "−1" : "+1";
          const label =
            sigma === 0
              ? `e${e}: uncut tree edge`
              : `e${e}: ${sigma < 0 ? "negative" : "positive"} energy pole`;
          const content = `<span class="lp-caption">${sub("e", e)}</span><span class="lp-cut-value">${value}</span>`;
          return sigma
            ? `<button type="button" class="lp-cut-entry cursor-interaction" data-trace-cut="${i}" data-hover-edge="${i}" data-cut-component="${e}" data-sigma="${sigma}" aria-label="${label}. Focus to explain sign; ${g.scope === "representation" ? `activate for next tree cutting e${e}.` : "fixed cut assignment."}">${content}</button>`
            : `<span class="lp-cut-entry" data-cut-component="${e}" data-sigma="0" aria-label="${label}">${content}</span>`;
        })
        .join("")}</div></div>`;
    }
    function ledger(st) {
      const { g, r } = selected(st);
      return `<table class="lp-ledger" aria-label="Routing and surfaces for the selected tree"><thead><tr><th>Edge</th><th>Global routing</th><th>On this cut</th><th>Factors</th></tr></thead><tbody>${treeEdges(
        g,
        r,
      )
        .map(
          (i) =>
            `<tr data-active="${i === st.edge}"><td data-label="Edge"><button class="cursor-interaction" data-tree-edge="${i}" data-hover-edge="${i}" aria-pressed="${i === st.edge}">${sub("e", g.edges[i])}</button></td><td data-label="Routing"><span>${routing(i, g)}</span></td><td data-label="On cut"><span>${linear(r.edge_map[i], g)}</span></td><td data-label="Factors"><span>${(r.propagator_factors[i] || []).map((id) => factor(g.surfaces[id], st, i)).join(" · ")}</span></td></tr>`,
        )
        .join("")}</tbody></table>`;
    }
    function graph() {
      return '<div class="lp-graph-wrap"><div class="lp-region-choice"><span class="lp-caption"></span></div><div data-graph-host></div></div>';
    }
    function legend() {
      return '<div class="lp-graph-legend"><span>Pole +: outward · −: inward</span><span><svg class="lp-arrow-key" aria-hidden="true" viewBox="0 0 30 12"><path class="lp-cut-mark" d="M1 1L7 6L1 11 M23 1L29 6L23 11"/></svg>Cut</span><span><svg class="lp-arrow-key" aria-hidden="true" viewBox="0 0 30 12"><path class="lp-cut-mark" d="M8 1L14 6L8 11 M16 1L22 6L16 11"/></svg>Selected edge</span></div>';
    }
    function boundarySigns(st) {
      const { g, r } = current(st),
        part = partition(g, r, st.edge);
      if (st.edge === null) return "";
      return `<div class="lp-boundary-signs lp-caption" aria-label="Momentum boundary sign times pole sign"><span>Momentum boundary × pole → surface term</span>${[
        ...r.cuts,
        st.edge,
      ]
        .flatMap((i) => {
          const e = g.routing[i],
            incidence = +part.has(e.tail) - +part.has(e.head),
            sigma = poleSign(st, i);
          if (!incidence) return [];
          return [
            `<span data-boundary-edge="${g.edges[i]}" data-incidence="${incidence}" data-pole="${sigma}" data-coefficient="${incidence * sigma}">${sub("e", g.edges[i])}: q ${incidence > 0 ? "out (+)" : "in (−)"} × (${sigma > 0 ? "+" : "−"}${g.energy_html[g.edges[i]]}) → <strong>${incidence * sigma > 0 ? "+" : "−"}${g.energy_html[g.edges[i]]}</strong></span>`,
          ];
        })
        .join("")}</div>`;
    }
    function mountGraph(st) {
      if (st.svg) return;
      const svg = native.cloneNode(true);
      svg.querySelectorAll("script").forEach((el) => el.remove());
      scopeIds(svg, instance);
      const labels = root
        .querySelector("[data-drawing]")
        .content.querySelector("[data-ltd-labels]")
        .cloneNode(true);
      scopeIds(labels, `${instance}-labels`);
      // Glyphs live in SVG use shadow trees: inherit the label's theme color
      // directly rather than relying on a descendant CSS selector crossing it.
      labels.querySelectorAll("[fill]").forEach((glyph) => {
        if (glyph.getAttribute("fill") !== "none")
          glyph.setAttribute("fill", "currentColor");
      });
      st.labels = new Map(
        [...labels.querySelectorAll("[data-label-key]")].map((label) => [
          label.dataset.labelKey,
          label,
        ]),
      );
      svg.append(labels.querySelector("defs"));
      // Retain native geometry for Linnet's half-edge shader. Only the energy
      // annotations are painted; particle labels and inspector targets stay hidden.
      const reference = make(
        "g",
        {
          visibility: "hidden",
          "pointer-events": "none",
          "aria-hidden": "true",
        },
        svg,
      );
      for (const child of [...svg.children])
        if (child !== reference && !["style", "defs"].includes(child.localName))
          reference.append(child);
      reference
        .querySelectorAll("[data-linnet-node]")
        .forEach((node) => node.setAttribute("r", 3.5 / physicalScale));
      svg.classList.add("lp-graph");
      svg.dataset.linnetViewportHeight = "180";
      svg.setAttribute("role", "group");
      st.el.querySelector("[data-graph-host]").append(svg);
      st.svg = svg;
      const script = document.createElement("script");
      script.textContent = navigation;
      st.el.querySelector("[data-graph-host]").append(script);
    }
    function scopeIds(svg, prefix) {
      const ids = new Map();
      svg.querySelectorAll("[id]").forEach((el) => {
        ids.set(el.id, `${prefix}-${el.id}`);
        el.id = ids.get(el.id);
      });
      svg.querySelectorAll("*").forEach((el) =>
        [...el.attributes].forEach((attr) => {
          let value = attr.value;
          if (value.startsWith("#") && ids.has(value.slice(1)))
            value = "#" + ids.get(value.slice(1));
          value = value.replace(/url\(#([^)]+)\)/g, (match, id) =>
            ids.has(id) ? `url(#${ids.get(id)})` : match,
          );
          if (value !== attr.value)
            el.setAttributeNS(attr.namespaceURI, attr.name, value);
        }),
      );
    }
    function setLabel(st, label, key) {
      const page = st.labels.get(key);
      const use = label.querySelector("use") || make("use", {}, label);
      use.setAttribute("href", `#${page.id}`);
      use.setAttribute("x", -page.dataset.width / 2);
      use.setAttribute("y", -page.dataset.height / 2);
      label.dataset.labelKey = key;
    }
    function drawGraph(st) {
      if (!st.graphOpen) return;
      mountGraph(st);
      const svg = st.svg;
      const { g, r } = selected(st),
        part = partition(g, r, st.edge),
        scale = physicalScale;
      // Linnet owns the camera. Changing a residue replaces only this layer.
      svg.querySelector("[data-ltd-drawing]")?.remove();
      const viewport = svg.querySelector(".linnet-viewport");
      const layer = make("g", { "data-ltd-drawing": "" }, viewport);
      const bounds = native.viewBox.baseVal,
        extent = {
          x: bounds.x - 16,
          y: bounds.y - 16,
          width: bounds.width + 32,
          height: bounds.height + 32,
        },
        maskId = `${instance}-cut-gaps`,
        defs = make("defs", {}, layer),
        mask = make(
          "mask",
          {
            id: maskId,
            maskUnits: "userSpaceOnUse",
            ...extent,
            style: "mask-type:luminance",
          },
          defs,
        );
      make("rect", { ...extent, fill: "white" }, mask);
      const clipped = make("g", { mask: `url(#${maskId})` }, layer);
      // Native shading already joins the half-edge curves and node geometry.
      // Composite it once, so a vertex does not darken once per incident edge.
      const regions = svg.linnetShade({
        nodes: [...part],
        half_edges: g.graph.paths
          .filter((p) => part.has(p.node))
          .map((p) => p["half-edge"]),
      });
      regions.classList.add("lp-region");
      regions.dataset.regionLayer = "";
      clipped.append(regions);
      make("g", { "data-cycle-layer": "", "pointer-events": "none" }, clipped);
      const paths = make("g", {}, clipped),
        nodes = make("g", {}, layer),
        marks = make("g", {}, layer),
        edges = new Map();
      make("g", { "data-cycle-arrows": "", "pointer-events": "none" }, layer);
      for (const p of g.graph.paths) {
        const node = g.graph.nodes.find((n) => n.id === p.node);
        if (!node) continue;
        const d =
          p.flow === "source"
            ? `M${node.x} ${node.y} L${p.d.slice(1)}`
            : `${p.d} L${node.x} ${node.y}`;
        const index = g.edges.indexOf(p.edge),
          cut = r.cuts.includes(index);
        const path = make(
          "path",
          {
            d,
            class: `lp-path${cut ? " is-cut" : ""}${index === st.edge ? " is-selected-edge" : ""}`,
            "data-edge": p.edge,
            "data-flow": p.flow,
          },
          paths,
        );
        if (!edges.has(p.edge) || p.flow === "source")
          edges.set(p.edge, { path, index });
        if (index >= 0) {
          const hit = make(
            "path",
            {
              d,
              fill: "none",
              stroke: "transparent",
              "stroke-width": 18,
              class: "cursor-interaction",
              "data-graph-edge": index,
              "data-hover-edge": index,
            },
            marks,
          );
          hit.addEventListener("click", () => selectGraphEdge(st, index));
        }
      }
      for (const n of g.graph.nodes)
        make(
          "circle",
          { cx: n.x, cy: n.y, r: 3.5 / scale, class: "lp-node" },
          nodes,
        );
      for (const [e, { path, index }] of edges) {
        const len = path.getTotalLength();
        // Reuse the very same curve and chevron painted for Feynman diagrams.
        // The momentum direction is fixed; residue pole signs never flip it.
        const source =
          index < 0 ? g.external_routing[e].source : g.routing[index].tail;
        const target =
          index < 0 ? g.external_routing[e].destination : g.routing[index].head;
        const arrow = make(
          "g",
          {
            class: "lp-routing-arrow",
            "data-momentum-edge": e,
            "data-inspected": index === st.edge,
            "data-momentum-source": source ?? "external",
            "data-momentum-target": target ?? "external",
          },
          marks,
        );
        for (const child of momentumArrows.get(e).querySelectorAll("path"))
          arrow.append(child.cloneNode(true));
        if (index < 0) {
          const label = make(
            "g",
            {
              transform: labelPositions.get(e),
              class: "lp-math-label",
              "aria-label": g.external_routing[e].label,
            },
            marks,
          );
          setLabel(st, label, `p-${e}`);
          continue;
        }
        const at = path.getPointAtLength(len),
          before = path.getPointAtLength(Math.max(0, len - 1));
        const angle = Math.atan2(at.y - before.y, at.x - before.x),
          sigma = poleSign(st, index);
        const label = make(
          "g",
          {
            transform: labelPositions.get(e),
            "data-edge-label": index,
            "data-hover-edge": index,
            class: "lp-math-label cursor-interaction",
            "aria-label": `q${e}`,
          },
          marks,
        );
        // Keep the hit area fixed when the hover equation grows around it.
        make(
          "rect",
          { x: -8, y: -7, width: 16, height: 14, fill: "transparent" },
          label,
        );
        setLabel(st, label, `q-${e}`);
        label.addEventListener("click", () => selectGraphEdge(st, index));
        if (sigma) {
          label.setAttribute(
            "aria-label",
            `q${e}⁰ = ${sigma < 0 ? "−" : "+"}Eᵒˢ${e}`,
          );
          // Region-relative sign glyph: outward means positive pole, not positive energy flux.
          // A cut away from the boundary has no inward/outward side, so use flat faces there.
          const cut = r.cuts.includes(index),
            gap = cut ? 22 : 8,
            height = 4,
            incidence = +part.has(source) - +part.has(target);
          const direction = incidence * sigma,
            side = incidence ? (sigma > 0 ? "out" : "in") : "none",
            reach = incidence ? 2.4 : 0;
          const mark = make(
            "g",
            {
              transform: `translate(${at.x} ${at.y}) rotate(${(angle * 180) / Math.PI + (direction < 0 ? 180 : 0)}) scale(${1 / scale})`,
              "data-pole-edge": e,
              "data-pole-sign": sigma,
              "data-pole-role": cut ? "cut" : "surface",
              "data-boundary-incidence": incidence,
              "data-boundary-coefficient": incidence * sigma,
              "data-region-direction": side,
              "data-gap-width": gap,
              "data-hover-edge": index,
              "aria-label": `${cut ? "Cut" : "Surface pole"} e${e}: q${e}⁰ = ${sigma < 0 ? "−" : "+"}Eᵒˢ${e}. ${incidence ? `Pole marker points ${sigma > 0 ? "out of" : "into"} the shaded region.` : "Off the region boundary; pole sign is shown in the label."}`,
              class: "lp-pole-marker cursor-interaction",
              role: "button",
              tabindex: "0",
            },
            marks,
          );
          if (cut) {
            mark.setAttribute("data-cut-edge", e);
            mark.setAttribute("data-cut-sign", sigma);
          } else mark.setAttribute("data-surface-edge", e);
          make(
            "path",
            {
              d: incidence
                ? `M${-gap / 2 - reach} -12V${-height}L${-gap / 2 + reach} 0L${-gap / 2 - reach} ${height}V12H${gap / 2 - reach}V${height}L${gap / 2 + reach} 0L${gap / 2 - reach} ${-height}V-12Z`
                : `M${-gap / 2} -12H${gap / 2}V12H${-gap / 2}Z`,
              fill: "black",
              transform: mark.getAttribute("transform"),
              "data-cutout-edge": e,
            },
            mask,
          );
          // Cut faces and the shading share one transparent aperture. The
          // smaller aperture inspects a tree-edge pole without making it a cut.
          make(
            "path",
            {
              d: [-gap / 2, gap / 2]
                .map((x) =>
                  incidence
                    ? `M${x - reach} ${-height}L${x + reach} 0L${x - reach} ${height}`
                    : `M${x} ${-height}V${height}`,
                )
                .join(" "),
              class: "lp-cut-mark",
            },
            mark,
          );
          make(
            "rect",
            { x: -17, y: -14, width: 34, height: 28, fill: "transparent" },
            mark,
          );
          const order = r.pole_orders?.[index] || 1;
          if (cut && order > 1) {
            const badge = make("text", { x: at.x, y: at.y - 12 / scale, "text-anchor": "middle", "font-size": 10 / scale, class: "lp-pole-order", "data-pole-order": order, "pointer-events": "none" }, marks);
            badge.textContent = `order ${order}`;
          }
          mark.addEventListener("click", () => selectGraphEdge(st, index));
          mark.addEventListener("keydown", (event) => {
            if (event.key === "Enter" || event.key === " ") {
              event.preventDefault();
              selectGraphEdge(st, index);
            }
          });
        }
      }
      updatePreview(st);
    }
    function updatePreview(st) {
      if (!st.graphOpen) return;
      const svg = st.el.querySelector(".lp-graph"),
        { g } = current(st),
        cycle = cutCycle(st);
      if (!svg) return;
      const back = svg.querySelector("[data-cycle-layer]"),
        arrows = svg.querySelector("[data-cycle-arrows]"),
        scale = physicalScale;
      back.replaceChildren();
      arrows.replaceChildren();
      for (const label of svg.querySelectorAll("[data-edge-label]")) {
        const index = +label.dataset.edgeLabel,
          e = g.edges[index],
          sigma = poleSign(st, index);
        if (!sigma) continue;
        const expanded = st.previewEdge === index;
        const incidence = +svg.querySelector(`[data-pole-edge="${e}"]`).dataset
          .boundaryIncidence;
        const sign = sigma < 0 ? "minus" : "plus";
        setLabel(
          st,
          label,
          expanded
            ? `pole-${e}-${sign}`
            : `q-${e}${incidence ? "" : `-${sign}`}`,
        );
        label.classList.toggle("lp-pole-sign", expanded || !incidence);
      }
      // Keep the marker's reference region visible while tracing its sign cycle.
      svg.querySelector("[data-region-layer]").style.opacity = cycle
        ? ".08"
        : "";
      st.el.querySelector(".lp-region-choice .lp-caption").innerHTML = cycle
        ? `Cycle along ${sub("e", g.edges[cycle.latest])} (${sub("k", cycle.loop)}, last) · ${sub("σ", g.edges[st.traceCut])} = ${poleSign(st, st.traceCut) > 0 ? "+" : "−"}1`
        : "";
      svg.setAttribute(
        "aria-label",
        cycle
          ? "Restored cycle explaining the hovered cut pole sign"
          : "Fixed half-edge region entered by the inspected momentum arrow",
      );
      if (!cycle) return;
      // Keep hit targets intact during hover so the cycle cannot interrupt clicks.
      for (const path of svg.querySelectorAll(".lp-path")) {
        const physical = +path.dataset.edge,
          index = g.edges.indexOf(physical),
          direction = cycle.directions.get(index);
        if (!direction) continue;
        make(
          "path",
          {
            d: path.getAttribute("d"),
            class: "lp-cycle",
            "data-cycle-edge": physical,
          },
          back,
        );
        if (path.dataset.flow !== "sink") continue;
        const len = path.getTotalLength(),
          p = path.getPointAtLength(len * 0.5),
          a = path.getPointAtLength(Math.max(0, len * 0.5 - 0.2)),
          b = path.getPointAtLength(Math.min(len, len * 0.5 + 0.2));
        const angle =
            Math.atan2(b.y - a.y, b.x - a.x) + (direction < 0 ? Math.PI : 0),
          size = 5 / scale;
        make(
          "path",
          {
            d: `M${-size} ${-size * 0.7}L${size} 0L${-size} ${size * 0.7}`,
            class: "lp-cycle-arrow",
            transform: `translate(${p.x} ${p.y}) rotate(${(angle * 180) / Math.PI})`,
            "stroke-width": 2 / scale,
            "data-cycle-arrow": physical,
            "data-cycle-direction": direction,
          },
          arrows,
        );
      }
    }
    function previewEdge(st, index, source) {
      if (st.rendering) return;
      const { g, r } = current(st);
      st[source] = Number.isInteger(index) && g.routing[index] ? index : null;
      const edge = st.hoverEdge ?? st.focusEdge;
      if (edge === st.previewEdge) return;
      st.previewEdge = edge;
      st.traceCut = r.cuts.includes(edge) ? edge : null;
      updatePreview(st);
    }
    function chooseResidue(st, id) {
      if (!Number.isInteger(id) || !data.residues[id]) {
        st.el.querySelector("[data-sum-input]").value = st.residue;
        return;
      }
      st.residue = id;
      st.message = "";
      st.traceCut = null;
      const { r } = current(st);
      if (!r.propagator_factors[st.edge])
        st.edge = Object.keys(r.propagator_factors).map(Number)[0] ?? null;
      st.surface = null;
      render(st);
    }
    function current(st) {
      const g = data;
      return { g, r: g.residues[st.residue] };
    }
    function chooseSurface(st, id, edge) {
      const { r } = current(st),
        branch = r.propagator_factors[edge]?.indexOf(id);
      st.edge = branch === 0 || branch === 1 ? edge : null;
      st.surface = id;
      st.branch = branch === 0 ? -1 : 1;
      st.traceCut = null;
      render(st);
    }
    function selectGraphEdge(st, index) {
      const { g, r } = selected(st);
      if (r.cuts.includes(index)) {
        if (data.scope === "residue") return;
        const choices = g.residues.filter((r) => r.cuts.includes(index)),
          at = choices.findIndex((r) => r.id === st.residue);
        chooseResidue(st, choices[(at + 1) % choices.length].id);
        return;
      }
      if (!r.propagator_factors[index]) return;
      st.surface = null;
      if (st.edge === index) st.branch = -st.branch;
      st.edge = index;
      st.traceCut = null;
      render(st);
    }
    function render(st) {
      // Replacing focused controls can emit focus events against the previous graph.
      st.rendering = true;
      const { g, r, s } = selected(st);
      for (const key of ["hoverEdge", "focusEdge"])
        if (!g.routing[st[key]]) st[key] = null;
      st.previewEdge = st.hoverEdge ?? st.focusEdge;
      st.traceCut = r.cuts.includes(st.previewEdge) ? st.previewEdge : null;
      const graphWrap = st.el.querySelector(".lp-graph-wrap"),
        focusedPole =
          document.activeElement?.closest?.("[data-pole-edge]")?.dataset
            .poleEdge;
      const header = `<header class="hs-meta"><code>${g.scope === "residue" ? "LtdResidue" : "LtdRepresentation"}</code><span>${g.scope === "residue" ? `Residue ${r.id} · ${r.cuts.length} cuts` : `${g.loops} ${g.loops === 1 ? "loop" : "loops"} · ${g.residues.length} residues`}</span></header>`;
      const sum = g.scope === "representation" ? root.expressionDisplay.sum({
        ids: g.residues.map((r) => r.id),
        current: st.residue,
        symbol: "R",
        total: "L",
        label: "Tree residue",
        onSelect(index) {
          chooseResidue(st, index);
        },
        onInvalid() {
          st.message = "Enter a residue ID retained in this result.";
          render(st);
        },
      }) : "";
      const drawing = `<details class="lp-details" data-residue-graph ${st.graphOpen ? "open" : ""}><summary>Explore cut graph</summary>${graph()}${legend()}</details>`;
      st.el.innerHTML = `${header}<div class="lp-body">${sum}<div class="hs-status" role="status" aria-live="polite" data-error="${!!st.message}">${st.message}</div>${equation(st)}${drawing}<div class="lp-definition hs-definition" aria-live="polite">${s ? `${surfaceName(s)} = ${linear(s.expression, g)}` : "No uncut propagator factors"}</div><details class="lp-details" data-residue-details ${st.detailsOpen ? "open" : ""}><summary class="cursor-interaction">Details</summary>${globalRouting(st)}${cutStructure(st)}${boundarySigns(st)}${ledger(st)}<div class="lp-footer">Scalar numerator 1 · dq⁰/(2πi), closed below · spatial measure separate. Parallel arrows: momentum routing. Gaps: pole sign relative to the shaded component, not energy flux.</div></details></div>`;
      if (graphWrap)
        st.el.querySelector(".lp-graph-wrap").replaceWith(graphWrap);
      st.el.querySelectorAll("[data-tree-edge]").forEach((b) =>
        b.addEventListener("click", () => {
          const focused = document.activeElement === b,
            index = +b.dataset.treeEdge;
          selectGraphEdge(st, index);
          if (focused)
            st.el
              .querySelector(`.lp-ledger [data-tree-edge="${index}"]`)
              ?.focus({ preventScroll: true });
        }),
      );
      st.el.querySelectorAll("[data-surface]").forEach((b) =>
        b.addEventListener("click", () => {
          const focused = document.activeElement === b,
            id = +b.dataset.surface,
            edge =
              b.dataset.factorEdge === "null" ? null : +b.dataset.factorEdge,
            container = b.closest(".lp-ledger") ? ".lp-ledger" : ".lp-equation";
          chooseSurface(st, id, edge);
          if (focused)
            st.el
              .querySelector(
                `${container} [data-factor-edge="${edge}"][data-surface="${id}"]`,
              )
              ?.focus({ preventScroll: true });
        }),
      );
      st.el.querySelectorAll("[data-trace-cut]").forEach((b) =>
        b.addEventListener("click", () => {
          const focused = document.activeElement === b,
            index = +b.dataset.traceCut;
          selectGraphEdge(st, index);
          if (focused)
            st.el
              .querySelector(`[data-trace-cut="${index}"]`)
              ?.focus({ preventScroll: true });
        }),
      );
      const graphDetails = st.el.querySelector("[data-residue-graph]");
      if (graphDetails) graphDetails.addEventListener("toggle", () => {
        if (graphDetails.isConnected && st.graphOpen !== graphDetails.open) {
          st.graphOpen = graphDetails.open;
          if (st.graphOpen) drawGraph(st);
        }
      });
      const details = st.el.querySelector("[data-residue-details]");
      details.querySelector("summary").addEventListener("click", (e) => {
        // Native toggle notifications are queued; retain the choice before another action renders.
        e.preventDefault();
        details.open = !details.open;
        st.detailsOpen = details.open;
      });
      details.addEventListener("toggle", (e) => {
        if (e.target.isConnected && st.detailsOpen !== e.target.open) {
          st.detailsOpen = e.target.open;
        }
      });
      drawGraph(st);
      if (focusedPole !== undefined && st.svg)
        st.svg
          .querySelector(`[data-pole-edge="${focusedPole}"]`)
          ?.focus({ preventScroll: true });
      st.rendering = false;
    }
    const el = root.querySelector(".ltd-product");
    if (!data.residues.length) {
      el.textContent = "LtdRepresentation: no retained residues";
      continue;
    }
    const st = {
      el,
      residue: 0,
      message: "",
      edge:
        Object.keys(data.residues[0].propagator_factors).map(Number)[0] ?? null,
      surface: null,
      branch: -1,
      detailsOpen: false,
      graphOpen: false,
      traceCut: null,
      hoverEdge: null,
      focusEdge: null,
      previewEdge: null,
      svg: null,
    };
    const hovered = (target) => {
      const hit = target?.closest?.("[data-hover-edge]");
      return hit && hit.dataset.hoverEdge !== "null"
        ? +hit.dataset.hoverEdge
        : null;
    };
    el.addEventListener("pointerover", (event) => {
      if (event.pointerType !== "touch")
        previewEdge(st, hovered(event.target), "hoverEdge");
    });
    el.addEventListener("pointerout", (event) => {
      if (event.pointerType !== "touch")
        previewEdge(st, hovered(event.relatedTarget), "hoverEdge");
    });
    el.addEventListener("focusin", (event) =>
      previewEdge(st, hovered(event.target), "focusEdge"),
    );
    el.addEventListener("focusout", (event) =>
      previewEdge(st, hovered(event.relatedTarget), "focusEdge"),
    );
    render(st);
    let width = el.clientWidth;
    const observer = new ResizeObserver(() => {
      if (!el.isConnected) {
        observer.disconnect();
        return;
      }
      if (el.clientWidth !== width) {
        width = el.clientWidth;
        drawGraph(st);
      }
    });
    observer.observe(el);
  }
})();
