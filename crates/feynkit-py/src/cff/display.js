(() => {
  for (const root of document.querySelectorAll(".feynkit-cff-result")) {
    if (root.dataset.ready) continue;
    root.dataset.ready = "true";
    const data = JSON.parse(
      root.querySelector("[data-cff]").textContent,
      (key, value) => key.endsWith("_html") ? root.expressionDisplay.math(value) : value,
    );
    const state = {
      orientation: 0,
      family: 0,
      surfaces: [...new Set(data.orientations[0]?.terms[0] || [])].slice(0, data.scope === "family" ? undefined : 2),
      contours: true,
      multiple: data.scope === "family",
      familyPage: 0,
      palette: {},
    };
    let message = "",
      messageError = false;
    data.familyCount = data.orientations.reduce(
      (count, o) => count + o.terms.length,
      0,
    );
    const familyId = (index) => data.orientations[state.orientation].family_ids[index];
    root.querySelector("[data-metadata]").textContent =
      (data.pole_order ? `Pole coefficient · order ${data.pole_order} · ` : "") +
      (data.scope === "representation"
        ? `${data.loops} ${data.loops === 1 ? "loop" : "loops"} · ${data.orientations.length} orientations · ${data.familyCount.toLocaleString("en-US")} terms`
        : `Orientation ${data.orientations[0].id}` + (data.scope === "family" ? ` · family ${familyId(0)}` : ""));
    const explorer = root.querySelector("[data-explorer]");
    root.querySelector("[data-navigation-help]").hidden = data.scope !== "representation";
    root.querySelector("[data-scope-help]").textContent =
      data.generalized ? "Includes the supplied numerator and on-shell factors; the spatial measure is separate." : "Includes the on-shell energy prefactor; numerators and the spatial measure are separate.";
    if (data.scope === "family") explorer.querySelector("summary").textContent = "Explore graph";
    const internalEdges = data.edges.filter((edge) =>
      data.orientations.some((o) =>
        ["default", "reversed"].includes(o.directions[edge.id]),
      ),
    );
    const directionKey = (directions) =>
      internalEdges.map((edge) => directions[edge.id] || "absent").join(",");
    const orientationByDirections = new Map(
      data.orientations.map((o, index) => [directionKey(o.directions), index]),
    );
    function flipTarget(edgeId) {
      const directions = { ...data.orientations[state.orientation].directions };
      if (!["default", "reversed"].includes(directions[edgeId]))
        return undefined;
      directions[edgeId] =
        directions[edgeId] === "reversed" ? "default" : "reversed";
      return orientationByDirections.get(directionKey(directions));
    }
    function paletteFor(term) {
      const palette = {},
        used = new Set();
      for (const id of new Set(term)) {
        const slot = state.palette[id];
        if (Number.isInteger(slot) && slot > 0 && !used.has(slot)) {
          palette[id] = slot;
          used.add(slot);
        }
      }
      for (const id of new Set(term))
        if (!palette[id]) {
          let slot = 1;
          while (used.has(slot)) slot++;
          palette[id] = slot;
          used.add(slot);
        }
      return palette;
    }
    const ns = "http://www.w3.org/2000/svg";
    const instance = "cff-" + Math.random().toString(36).slice(2);
    let serial = 0;
    const eta = (id) =>
      `${data.surfaces[id].kind === "h" ? "H" : "η"}<sub>${data.surfaces[id].index}</sub>`;
    const factorTrees = new Map();
    function factorTree(o) {
      if (!factorTrees.has(o.id))
        factorTrees.set(
          o.id,
          root.expressionDisplay.factorTerms(
            o.terms.map((factors, id) => ({ factors: [...factors, ...o.contributions[id].energies.map(e => `energy-${e}`)], id })),
          ),
        );
      return factorTrees.get(o.id);
    }
    function factoredMath(node, term, selected, key = "root") {
      const active = node.terms.includes(state.family);
      const factors = root.expressionDisplay.powers(node.factors)
        .map(([id, power], index) => {
          const exponent = power > 1 ? `<sup data-power-label>${power}</sup>` : "";
          if (typeof id === "string") return root.expressionDisplay.energyFactor(data.energy_html[id.slice(7)], id.slice(7), power);
          if (!active)
            return `<span class="hs-factored-factor" data-factor-id="${id}" data-factor-power="${power}" data-path-active="false">${eta(id)}${exponent}<span class="hs-factor-spacer"></span></span>`;
          return `<button type="button" class="hs-factored-factor" data-factor-id="${id}" data-factor-power="${power}" data-path-active="true" data-surface="${id}" data-factor-key="${key}-${index}" style="--hs-color:${color(id, term)}" aria-label="Select ${data.surfaces[id].kind} surface ${data.surfaces[id].index} on family ${familyId(state.family)}" aria-pressed="${selected.includes(id)}">${eta(id)}${exponent}<span class="hs-swatch"></span></button>`;
        })
        .join("");
      const leaf = !node.children.length,
        contribution = leaf ? data.orientations[state.orientation].contributions[node.terms[0]] : null,
        coefficient = contribution?.coefficient || "1",
        numerator = contribution?.numerator_html || "1",
        value = `${coefficient === "1" ? "" : coefficient === "-1" ? "−" : coefficient.replaceAll("-", "−") + (numerator !== "1" ? " · " : "")}${numerator === "1" && coefficient !== "1" && coefficient !== "-1" ? "" : numerator}`;
      const explanation = contribution ? Object.entries(contribution.energy_map).map(([edge, value]) => `q${edge}⁰ = ${value}`).join("; ") : "";
      const fraction = factors
        ? root.expressionDisplay.fraction(`<span class="hs-contribution-numerator">${value}</span>`, factors, `data-path-active="${active}" data-shared-families="${node.terms.length}"`)
        : "";
      if (!node.children.length)
        return `<span class="hs-factor-leaf" data-family-leaf="${node.terms[0]}" data-coefficient="${coefficient}" title="${explanation.replaceAll("&", "&amp;").replaceAll('"', "&quot;").replaceAll("<", "&lt;")}">${fraction || `<span class="hs-contribution-numerator">${value}</span>`}<span class="hs-leaf-label" data-path-active="${active}">F${familyId(node.terms[0])}</span></span>`;
      const branches = node.children
        .map(
          (child, index) =>
            `<span class="hs-factor-summand">${index ? '<span aria-hidden="true">+</span>' : ""}${factoredMath(child, term, selected, `${key}-${index}`)}</span>`,
        )
        .join("");
      const bracket = `<span class="hs-factored-bracket"><span class="hs-bracket-edge" aria-hidden="true"></span><span class="hs-factored-sum">${branches}</span><span class="hs-bracket-edge" aria-hidden="true"></span></span>`;
      return `<span class="hs-factor-node">${fraction}${fraction ? '<span aria-hidden="true">·</span>' : ""}${bracket}</span>`;
    }
    const definition = (s) => s.expression_html;
    const color = (id, term) => {
      const slot = paletteFor(term)[id];
      return slot <= 7
        ? `var(--hs-${slot})`
        : `light-dark(hsl(${(slot * 137.508) % 360} 55% 38%),hsl(${(slot * 137.508) % 360} 65% 73%))`;
    };
    function element(name, attrs) {
      const el = document.createElementNS(ns, name);
      Object.entries(attrs).forEach(([k, v]) => el.setAttribute(k, v));
      return el;
    }
    function offsetPath(path, offset) {
      const length = path.getTotalLength(),
        points = [];
      for (let i = 0; i <= 24; i++) {
        const t = (length * i) / 24,
          a = path.getPointAtLength(Math.max(0, t - 0.2)),
          b = path.getPointAtLength(Math.min(length, t + 0.2));
        const p = path.getPointAtLength(t),
          dx = b.x - a.x,
          dy = b.y - a.y,
          n = Math.hypot(dx, dy) || 1;
        points.push(
          `${(p.x - (dy * offset) / n).toFixed(2)},${(p.y + (dx * offset) / n).toFixed(2)}`,
        );
      }
      return "M" + points.join(" L");
    }
    // Linnet supplies structural carriers separately from the painted particle
    // lines. Their metadata identifies the physical edge and the owning vertex.
    const template = root.querySelector("[data-drawing]").content;
    const native = template.querySelector("svg");
    native.dataset.linnetViewportHeight = "145";
    const navigation = native.querySelector("script").textContent;
    native.querySelectorAll("script").forEach((el) => el.remove());
    const carriers = [...native.querySelectorAll("[data-linnet-carrier]")].map(
      (anchor) => ({
        ...JSON.parse(anchor.dataset.linnetDetail),
        d: anchor.querySelector("path").getAttribute("d"),
      }),
    );
    native
      .querySelectorAll("a:not([data-linnet-carrier]),:scope > path")
      .forEach((el) => el.remove());
    native.querySelectorAll("[data-linnet-carrier]").forEach((el) => {
      el.removeAttribute("data-linnet-kind");
    });
    native.querySelectorAll("[fill]").forEach((el) => {
      if (!["none", "transparent"].includes(el.getAttribute("fill")))
        el.setAttribute("fill", "var(--hs-fg)");
    });
    native.querySelectorAll("[stroke]").forEach((el) => {
      if (!["none", "transparent"].includes(el.getAttribute("stroke")))
        el.setAttribute("stroke", "var(--hs-muted)");
    });
    for (const c of carriers)
      native.insertBefore(
        element("path", {
          d: c.d,
          fill: "none",
          stroke: "var(--hs-muted)",
          "stroke-width": 1,
          "data-edge": c.edge,
          "data-owner": c.node,
          "data-flow": c.flow,
        }),
        native.firstChild,
      );
    function drawing(host, term, selected, { mini = false } = {}) {
      const previous = mini
        ? null
        : host.querySelector("svg[data-linnet-interactive]");
      const fragment = template.cloneNode(true);
      const svg = fragment.querySelector("svg"),
        prefix = `${instance}-${serial++}-`,
        ids = new Map();
      // The native renderer puts each typeset label in a top-level group.
      // Keep its edge paths and vertex positions for family previews.
      if (mini)
        svg.querySelectorAll(":scope > g").forEach((label) => label.remove());
      svg.querySelectorAll("[id]").forEach((el) => {
        const old = el.id;
        el.id = prefix + old;
        ids.set(old, el.id);
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
      svg.setAttribute("role", mini ? "img" : "group");
      svg.setAttribute(
        "aria-label",
        `Energy-flow orientation ${data.orientations[state.orientation].id}; ${selected.length} selected surfaces. Rings mark vertex membership.`,
      );
      const bands = element("g", { "data-bands": "", "pointer-events": "none" });
      const paths = [...svg.querySelectorAll("path[data-edge]")],
        nodes = [...svg.querySelectorAll("[data-linnet-node]")];
      if (mini)
        for (const path of paths) {
          // Native edges stop at the label-sized vertex circles. Extend them
          // to the small preview dots while retaining the routed curve.
          const node = nodes.find(
              (node) => node.dataset.linnetNode === path.dataset.owner,
            ),
            center = `${node.getAttribute("cx")} ${node.getAttribute("cy")}`,
            d = path.getAttribute("d");
          const start = path.getPointAtLength(0);
          path.setAttribute(
            "d",
            path.dataset.flow === "source"
              ? `M ${center} L ${start.x} ${start.y} ${d}`
              : `${d} L ${center}`,
          );
        }
      const displayed = [...new Set(mini ? term : selected)];
      for (const edge of data.edges) {
        const halves = paths.filter((p) => Number(p.dataset.edge) === edge.id);
        const boundaryIds = displayed.filter(
          (id) =>
            data.surfaces[id].support.includes(edge.id),
        );
        halves.forEach((path) => {
          const owner = Number(path.dataset.owner);
          boundaryIds.forEach((id, slot) => {
            if (!data.surfaces[id].v.includes(owner)) return;
            bands.append(
              element("path", {
                d: offsetPath(
                  path,
                  (slot - (boundaryIds.length - 1) / 2) * (mini ? 4 : 3),
                ),
                fill: "none",
                stroke: color(id, term),
                "stroke-width": 2.4,
                "stroke-linecap": "round",
                "data-boundary-half": edge.id,
                "data-owner": owner,
                "data-surface": id,
              }),
            );
          });
        });
      }
      svg.append(bands);
      if (mini)
        for (const node of nodes) {
          node.setAttribute("r", "4.2");
          svg.append(node);
        }
      else
        for (const node of nodes) {
          const members = displayed.filter((id) =>
            data.surfaces[id].v.includes(Number(node.dataset.linnetNode)),
          );
          const x = Number(node.getAttribute("cx")),
            y = Number(node.getAttribute("cy")),
            r = Number(node.getAttribute("r"));
          members.forEach((id, index) => {
            const attrs = {
              fill: "none",
              stroke: color(id, term),
              "stroke-width": 2.5,
              "data-node-membership": id,
              "data-vertex": node.dataset.linnetNode,
              "pointer-events": "none",
            };
            if (members.length === 1)
              svg.append(element("circle", { ...attrs, cx: x, cy: y, r }));
            else {
              const a =
                  -Math.PI / 2 + (index * 2 * Math.PI) / members.length + 0.035,
                b =
                  -Math.PI / 2 +
                  ((index + 1) * 2 * Math.PI) / members.length -
                  0.035;
              svg.append(
                element("path", {
                  ...attrs,
                  d: `M ${x + r * Math.cos(a)} ${y + r * Math.sin(a)} A ${r} ${r} 0 0 1 ${x + r * Math.cos(b)} ${y + r * Math.sin(b)}`,
                }),
              );
            }
          });
        }
      const o = data.orientations[state.orientation];
      for (const edge of internalEdges) {
        const halves = paths.filter(
          (path) => Number(path.dataset.edge) === edge.id,
        );
        const source = halves.find((p) => p.dataset.flow === "source");
        if (!source) continue;
        const length = source.getTotalLength(),
          a = source.getPointAtLength(Math.max(0, length - 0.2)),
          point = source.getPointAtLength(length);
        const reversed = o.directions[edge.id] === "reversed",
          angle =
            (Math.atan2(point.y - a.y, point.x - a.x) * 180) / Math.PI +
            (reversed ? 180 : 0);
        const arrow = element("g", {
          transform: `translate(${point.x} ${point.y}) rotate(${angle})`,
        });
        const shape = mini
          ? "M -7 -4.5 L 0 0 L -7 4.5"
          : "M -3 -2.1 L 0 0 L -3 2.1";
        arrow.append(
          element("path", {
            d: shape,
            fill: "none",
            stroke: "var(--hs-fg)",
            "stroke-width": mini ? 1.8 : 1,
            "stroke-linejoin": "round",
            "pointer-events": "none",
          }),
        );
        if (!mini && data.scope === "representation") {
          arrow.classList.add("hs-arrowhead");
          arrow.dataset.flipEdge = edge.id;
          arrow.setAttribute("role", "button");
          arrow.setAttribute("tabindex", "0");
          arrow.setAttribute(
            "aria-disabled",
            String(flipTarget(edge.id) === undefined),
          );
          const label = `Reverse e${edge.id}: v${reversed ? edge.target : edge.source} → v${reversed ? edge.source : edge.target}`;
          arrow.setAttribute("aria-label", label);
          const title = element("title", {});
          title.textContent = label;
          arrow.prepend(title);
          arrow.append(element("path", { d: shape, class: "hs-arrow-hit" }));
        }
        svg.append(arrow);
      }
      if (mini) {
        svg.removeAttribute("data-linnet-interactive");
        svg.querySelectorAll("style").forEach((el) => el.remove());
        host.replaceChildren(svg);
        // Center the vertices rather than letting asymmetric external legs
        // shift the graph. Keep enough room on both sides for the full drawing.
        const box = svg.getBBox(),
          xs = nodes.map((node) => Number(node.getAttribute("cx"))),
          ys = nodes.map((node) => Number(node.getAttribute("cy"))),
          cx = nodes.length
            ? (Math.min(...xs) + Math.max(...xs)) / 2
            : box.x + box.width / 2,
          cy = nodes.length
            ? (Math.min(...ys) + Math.max(...ys)) / 2
            : box.y + box.height / 2,
          halfWidth = Math.max(cx - box.x, box.x + box.width - cx) + 6,
          halfHeight = Math.max(cy - box.y, box.y + box.height - cy) + 6;
        svg.setAttribute(
          "viewBox",
          `${cx - halfWidth} ${cy - halfHeight} ${2 * halfWidth} ${2 * halfHeight}`,
        );
        svg.setAttribute("preserveAspectRatio", "xMidYMid meet");
        return;
      }
      // Keep the shared camera alive while replacing only its drawing. Panning,
      // zoom and fit state survive orientation, family and surface changes.
      const content = element("g", { "data-cff-drawing": "" });
      for (const child of [...svg.children])
        if (child.localName !== "style") content.append(child);
      if (previous) {
        previous.querySelector("[data-cff-drawing]").replaceWith(content);
        previous.setAttribute("aria-label", svg.getAttribute("aria-label"));
      } else {
        svg.append(content);
        host.replaceChildren(svg);
        // As in diagram collections, cloned template scripts need activation.
        const script = document.createElement("script");
        script.textContent = navigation;
        host.append(script);
      }
      const viewer = previous || svg;
      // Nested circlings keep the same spacing when another surface is toggled.
      // Linnet outlines the union of the native vertices and internal paths.
      const levels = new Map();
      [...new Set(term)]
        .sort(
          (a, b) =>
            data.surfaces[a].v.length - data.surfaces[b].v.length || a - b,
        )
        .forEach((id) => {
          const members = data.surfaces[id].v;
          const contained = [...levels].filter(([other]) =>
            data.surfaces[other].v.every((vertex) => members.includes(vertex)),
          );
          levels.set(id, Math.max(0, ...contained.map(([, level]) => level + 1)));
        });
      for (const id of displayed) {
        if (!data.surfaces[id].v.length) continue;
        const members = data.surfaces[id].v,
          inside = nodes.filter((node) =>
            members.includes(Number(node.dataset.linnetNode)),
          ),
          padding = 4 + 4 * levels.get(id),
          radius =
            Math.max(0, ...inside.map((node) => Number(node.getAttribute("r")))) +
            padding;
        const region = viewer.linnetShade(
          {
            nodes: members,
            edges: data.edges
              .filter(
                (edge) =>
                  members.includes(edge.source) && members.includes(edge.target),
              )
              .map((edge) => edge.id),
          },
          state.contours
            ? { width: 2 * radius, padding, outline: 1.6, fillOpacity: 0.09 }
            : {},
        );
        region.style.color = color(id, term);
        region.dataset.shadingSurface = id;
        if (state.contours) {
          region.dataset.outlineSurface = id;
          region.dataset.outlineLevel = levels.get(id);
        }
        content.prepend(region);
      }
    }
    function renderInputs(product, o) {
      const singleFamily = data.scope === "family";
      product.querySelector(".hs-family-bar").hidden = singleFamily || !o.terms.length;
      product.querySelector("[data-family-previews]").hidden = singleFamily;
      if (singleFamily) return;
      const start = state.familyPage * 3,
        end = Math.min(start + 3, o.terms.length);
      product.querySelector("[data-page-range]").textContent =
        `${start + 1}–${end} of ${o.terms.length}`;
      product.querySelector('[data-page-step="-1"]').disabled = start === 0;
      product.querySelector('[data-page-step="1"]').disabled =
        end === o.terms.length;
      const previews = product.querySelector("[data-family-previews]");
      previews.innerHTML = o.terms
        .slice(start, end)
        .map(
          (term, j) =>
            `<button type="button" class="hs-family-choice" data-choose-family="${start + j}" aria-label="Choose family ${familyId(start + j)}" aria-pressed="${state.family === start + j}"><span class="hs-family-label">F${familyId(start + j)}</span><div class="hs-mini" data-preview-family="${start + j}"></div></button>`,
        )
        .join("");
      previews.querySelectorAll("[data-preview-family]").forEach((host) => {
        const term = o.terms[Number(host.dataset.previewFamily)];
        drawing(host, term, term, { mini: true });
      });
    }

    function render() {
      const o = data.orientations[state.orientation],
        term = o.terms[state.family] || [],
        selected = [...new Set(term)].filter((id) =>
          state.surfaces.includes(id),
        );
      state.palette = paletteFor(term);
      const product = root;
      {
        product.querySelector("[data-status]").textContent = message;
        product.querySelector("[data-status]").dataset.error =
          String(messageError);
        // Hidden SVGs have no useful viewport to measure. Initialize and refresh
        // the graph only when expanded, retaining its camera across toggles.
        if (explorer.open) {
          renderInputs(product, o);
          drawing(product.querySelector("[data-graph]"), term, selected);
        }
        product.querySelector("[data-total-sum]").hidden = data.scope !== "representation";
        product.querySelector("[data-total-sum]").innerHTML = data.scope === "representation"
          ? root.expressionDisplay.sum({
            ids: data.orientations.map((o) => o.id),
            current: state.orientation,
            symbol: "C",
            total: "C",
            label: "Orientation",
            onSelect(index) {
              chooseOrientation(index);
              render();
            },
            onInvalid() {
              message = "Enter an orientation ID retained in this result.";
              messageError = true;
              render();
            },
          }) : "";
        const tree = o.terms.length ? factorTree(o) : null,
          shared = tree?.factors.length || 0;
        product.querySelector("[data-orientation-sum]").innerHTML =
          `<div class="hs-sum-label" ${data.scope === "family" ? "hidden" : ""}>${o.terms.length} ${o.terms.length === 1 ? "family" : "families"}${o.terms.length > 1 && shared ? ` · ${shared} shared ${shared === 1 ? "factor" : "factors"}` : ""}</div><div class="hs-factored-equation" aria-label="Factored denominator expression for orientation ${o.id}, with family ${familyId(state.family)} highlighted"><span class="hs-factored-lhs">${data.scope === "family" ? `F<sub>${o.id},${familyId(state.family)}</sub>` : `C<sub>${o.id}</sub>`} =</span>${data.generalized ? (data.normalization === "1" ? "" : `<span>${data.normalization_html}</span><span>·</span>`) : data.energy_edges.length ? root.expressionDisplay.fraction(data.energy_edges.length % 2 ? "−1" : "1", `<span class="hs-energy">∏<sub>e</sub> 2${data.energy_html.e}</span>`, `title="On-shell factors for edges ${data.energy_edges.join(", ")}"`) + "<span>·</span>" : ""}${tree ? factoredMath(tree, term, selected) : "0"}</div>`;
        product.querySelector("[data-contours]").checked = state.contours;
        product.querySelector("[data-multiple]").checked = state.multiple;
        product.querySelector("[data-selection-count]").textContent =
          `${selected.length} surface${selected.length === 1 ? "" : "s"} selected`;
        product.querySelector("[data-detail]").innerHTML = selected.length
          ? selected
              .map((id) => {
                const surface = data.surfaces[id];
                return `<div class="hs-definition-row" style="--hs-color:${color(id, term)}" data-definition="${id}"><div class="hs-definition">${eta(id)} = ${definition(surface)}</div><div class="hs-sets">${surface.origin === "helper" ? "Algebraic helper · no graph region" : surface.numerator_only ? "Numerator factor" : surface.v.length ? `Inside {${surface.v.join(", ")}}` : "No graph region"}</div></div>`;
              })
              .join("")
          : '<span class="hs-muted">Select a surface to inspect its definition.</span>';
      }
    }
    explorer.addEventListener("toggle", () => {
      if (explorer.open && data.orientations.length) render();
    });
    function chooseFamily(index) {
      const term = data.orientations[state.orientation].terms[index] || [];
      const old = state.surfaces,
        dropped = old.filter((id) => !term.includes(id));
      state.family = index;
      state.familyPage = Math.floor(index / 3);
      state.surfaces = [...new Set(term)].filter((id) => old.includes(id));
      message = dropped.length
        ? `Unavailable in this family: ${dropped.map((id) => (data.surfaces[id].kind === "h" ? "H" : "η") + data.surfaces[id].index).join(", ")}.`
        : "";
      messageError = false;
    }
    function chooseOrientation(index) {
      const previous =
        data.orientations[state.orientation].terms[state.family] || [];
      const o = data.orientations[index];
      const scores = o.terms.map(
        (term) =>
          (previous.length + 1) *
            term.filter((id) => state.surfaces.includes(id)).length +
          term.filter((id) => previous.includes(id)).length,
      );
      state.orientation = index;
      chooseFamily(Math.max(0, scores.indexOf(Math.max(...scores))));
    }
    root.addEventListener("keydown", (event) => {
      if (
        (event.key === "Enter" || event.key === " ") &&
        event.target.matches("[data-flip-edge]")
      ) {
        event.preventDefault();
        event.target.dispatchEvent(new MouseEvent("click", { bubbles: true }));
        return;
      }
    });
    root.addEventListener("change", (event) => {
      const product = event.target.closest(".feynkit-cff-result");
      if (!product) return;
      if (event.target.matches("[data-contours]"))
        state.contours = event.target.checked;
      else if (event.target.matches("[data-multiple]"))
        state.multiple = event.target.checked;
      else return;
      const term =
        data.orientations[state.orientation].terms[state.family] || [];
      state.surfaces = [...new Set(term)].filter((id) =>
        state.surfaces.includes(id),
      );
      render();
    });
    root.addEventListener("click", (event) => {
      const product = event.target.closest(".feynkit-cff-result");
      if (!product) return;
      const action = event.target.closest(
        "[data-flip-edge],[data-choose-family],[data-page-step]",
      );
      if (action) {
        let focusSelector = "";
        if (action.hasAttribute("data-flip-edge")) {
          const edge = Number(action.dataset.flipEdge),
            next = flipTarget(edge);
          if (next === undefined) {
            message = `Reversing e${edge} is unavailable (a directed cycle or excluded orientation).`;
            messageError = true;
          } else chooseOrientation(next);
          focusSelector = `[data-graph] [data-flip-edge="${edge}"]`;
        } else if (action.hasAttribute("data-choose-family")) {
          chooseFamily(Number(action.dataset.chooseFamily));
          focusSelector = `[data-family-previews] [data-choose-family="${state.family}"]`;
        } else if (action.hasAttribute("data-page-step")) {
          state.familyPage += Number(action.dataset.pageStep);
          focusSelector = `[data-page-step="${action.dataset.pageStep}"]`;
        }
        render();
        product.querySelector(focusSelector)?.focus({ preventScroll: true });
        return;
      }
      const target = event.target.closest("[data-surface]");
      if (!target) return;
      const id = Number(target.dataset.surface);
      if (event.shiftKey || state.multiple)
        state.surfaces = state.surfaces.includes(id)
          ? state.surfaces.filter((other) => other !== id)
          : [...state.surfaces, id];
      else state.surfaces = [id];
      render();
      if (target.dataset.factorKey)
        product
          .querySelector(`[data-factor-key="${target.dataset.factorKey}"]`)
          ?.focus({ preventScroll: true });
    });
    if (data.orientations.length) render();
    else {
      for (const child of [...root.children])
        if (
          !["HEADER", "SCRIPT", "TEMPLATE", "NOSCRIPT"].includes(child.tagName)
        )
          child.remove();
      const empty = document.createElement("p");
      empty.textContent = "C = 0 · No orientations retained.";
      root.append(empty);
    }
  }
})();
