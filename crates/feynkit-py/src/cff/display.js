(() => {
  for (const root of document.querySelectorAll(".feynkit-cff-result")) {
    if (root.dataset.ready) continue;
    root.dataset.ready = "true";
    const data = JSON.parse(root.querySelector("[data-cff]").textContent);
    const state = {
      orientation: 0,
      family: 0,
      surfaces: [...new Set(data.orientations[0]?.terms[0] || [])].slice(0, 2),
      contours: false,
      multiple: false,
      familyPage: 0,
      palette: {},
    };
    let message = "",
      messageError = false;
    data.familyCount = data.orientations.reduce(
      (count, o) => count + o.terms.length,
      0,
    );
    root.querySelector("[data-metadata]").textContent =
      `${data.loops} ${data.loops === 1 ? "loop" : "loops"} · ${data.orientations.length} orientations · ${data.familyCount.toLocaleString("en-US")} terms`;
    const explorer = root.querySelector("[data-explorer]");
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
    function factorFamilies(entries) {
      const remaining = entries.map((entry) => ({
          ...entry,
          factors: [...entry.factors],
        })),
        factors = [];
      // Extract the multiset intersection; each leaf retains its original family ID.
      for (const id of entries[0].factors)
        if (remaining.every((entry) => entry.factors.includes(id))) {
          factors.push(id);
          remaining.forEach((entry) =>
            entry.factors.splice(entry.factors.indexOf(id), 1),
          );
        }
      const node = {
        factors,
        families: entries.map((entry) => entry.family),
        children: [],
      };
      if (entries.length === 1) return node;
      const groups = new Map();
      remaining.forEach((entry) => {
        const key = entry.factors.length
          ? `surface-${entry.factors[0]}`
          : `family-${entry.family}`;
        if (!groups.has(key)) groups.set(key, []);
        groups.get(key).push(entry);
      });
      node.children = [...groups.values()].map(factorFamilies);
      return node;
    }
    const factorTrees = new Map();
    function factorTree(o) {
      if (!factorTrees.has(o.id))
        factorTrees.set(
          o.id,
          factorFamilies(
            o.terms.map((factors, family) => ({ factors, family })),
          ),
        );
      return factorTrees.get(o.id);
    }
    function factoredMath(node, term, selected, key = "root") {
      const active = node.families.includes(state.family);
      const factors = node.factors
        .map((id, index) => {
          if (!active)
            return `<span class="hs-factored-factor" data-factor-id="${id}" data-path-active="false">${eta(id)}<span class="hs-factor-spacer"></span></span>`;
          return `<button type="button" class="hs-factored-factor" data-factor-id="${id}" data-path-active="true" data-surface="${id}" data-factor-key="${key}-${index}" style="--hs-color:${color(id, term)}" aria-label="Select ${data.surfaces[id].kind} surface ${data.surfaces[id].index} on family ${state.family + 1}" aria-pressed="${selected.includes(id)}">${eta(id)}<span class="hs-swatch"></span></button>`;
        })
        .join("");
      const fraction = factors
        ? `<span class="hs-factored-fraction" data-path-active="${active}" data-shared-families="${node.families.length}"><span>1</span><span class="hs-factored-denominator">${factors}</span></span>`
        : "";
      if (!node.children.length)
        return `<span class="hs-factor-leaf" data-family-leaf="${node.families[0]}">${fraction || "<span>1</span>"}<span class="hs-leaf-label" data-path-active="${active}">F${node.families[0] + 1}</span></span>`;
      const branches = node.children
        .map(
          (child, index) =>
            `<span class="hs-factor-summand">${index ? '<span aria-hidden="true">+</span>' : ""}${factoredMath(child, term, selected, `${key}-${index}`)}</span>`,
        )
        .join("");
      const bracket = `<span class="hs-factored-bracket"><span class="hs-bracket-edge" aria-hidden="true"></span><span class="hs-factored-sum">${branches}</span><span class="hs-bracket-edge" aria-hidden="true"></span></span>`;
      return `<span class="hs-factor-node">${fraction}${fraction ? '<span aria-hidden="true">·</span>' : ""}${bracket}</span>`;
    }
    function sumSequence(indices, label) {
      const entries = [...new Set(indices)].sort((a, b) => a - b),
        parts = [];
      entries.forEach((index, j) => {
        if (j) {
          parts.push("<span>+</span>");
          if (index > entries[j - 1] + 1)
            parts.push("<span>⋯</span><span>+</span>");
        }
        parts.push(label(index));
      });
      return parts.join("");
    }
    const definition = (s) => {
      const parts = [
        ...s.e.map((e) => [1, `E<sup>os</sup><sub>${e}</sub>`]),
        ...s.negative.map((e) => [-1, `E<sup>os</sup><sub>${e}</sub>`]),
        ...s.q.map(([q, n]) => [n, `Q<sup>0</sup><sub>${q}</sub>`]),
      ];
      return (
        parts
          .map(
            ([n, label], i) =>
              `${n < 0 ? " − " : i ? " + " : ""}${Math.abs(n) === 1 ? "" : Math.abs(n)}${label}`,
          )
          .join("") || "0"
      );
    };
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
    native.querySelectorAll("a,:scope > path").forEach((el) => el.remove());
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
      const underlay = element("g", { "data-shading": "" }),
        bands = element("g", { "data-bands": "", "pointer-events": "none" });
      svg.insertBefore(underlay, svg.firstChild);
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
            data.surfaces[id].e.includes(edge.id) ||
            data.surfaces[id].negative.includes(edge.id),
        );
        const interiorIds = mini
          ? []
          : displayed.filter(
              (id) =>
                data.surfaces[id].v.includes(edge.source) &&
                data.surfaces[id].v.includes(edge.target),
            );
        halves.forEach((path) => {
          const owner = Number(path.dataset.owner);
          interiorIds.forEach((id, slot) =>
            underlay.append(
              element("path", {
                d: offsetPath(
                  path,
                  (slot - (interiorIds.length - 1) / 2) * 3.5,
                ),
                fill: "none",
                stroke: color(id, term),
                "stroke-width": 3.2,
                opacity: 0.16,
                "stroke-linecap": "round",
                "data-interior": edge.id,
                "data-interior-surface": id,
              }),
            ),
          );
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
      if (state.contours && !mini) {
        // Give nested circlings distinct offsets, using the whole family so
        // selecting or deselecting a surface does not move the other outlines.
        const levels = new Map();
        [...new Set(term)]
          .sort(
            (a, b) =>
              data.surfaces[a].v.length - data.surfaces[b].v.length || a - b,
          )
          .forEach((id) => {
            const members = data.surfaces[id].v;
            const contained = [...levels].filter(([other]) =>
              data.surfaces[other].v.every((vertex) =>
                members.includes(vertex),
              ),
            );
            levels.set(
              id,
              Math.max(0, ...contained.map(([, level]) => level + 1)),
            );
          });
        displayed.forEach((id) => {
          const members = data.surfaces[id].v,
            inside = nodes.filter((node) =>
              members.includes(Number(node.dataset.linnetNode)),
            ),
            padding = 4 + 3 * levels.get(id),
            radius =
              Math.max(
                0,
                ...inside.map((node) => Number(node.getAttribute("r"))),
              ) + padding;
          const filter = element("filter", {
            id: prefix + "outline-" + id,
            x: "-50%",
            y: "-50%",
            width: "200%",
            height: "200%",
          });
          filter.innerHTML =
            '<feMorphology in="SourceAlpha" operator="dilate" radius="0.8" result="outer"/><feComposite in="outer" in2="SourceAlpha" operator="out" result="boundary"/>';
          filter.append(element("feFlood", { "flood-color": color(id, term) }));
          filter.append(
            element("feComposite", { in2: "boundary", operator: "in" }),
          );
          svg.append(filter);
          const region = element("g", {
            "data-outline-surface": id,
            "data-outline-level": levels.get(id),
            filter: `url(#${filter.id})`,
            fill: color(id, term),
            stroke: color(id, term),
            "stroke-width": 2 * radius,
            opacity: 0.8,
            "pointer-events": "none",
          });
          inside.forEach((node) =>
            region.append(
              element("circle", {
                cx: node.getAttribute("cx"),
                cy: node.getAttribute("cy"),
                r: radius,
                stroke: "none",
                "data-outline-vertex": node.dataset.linnetNode,
              }),
            ),
          );
          paths.forEach((path) => {
            const edge = data.edges.find(
              (edge) => edge.id === Number(path.dataset.edge),
            );
            if (members.includes(edge.source) && members.includes(edge.target))
              region.append(
                element("path", {
                  d: path.getAttribute("d"),
                  fill: "none",
                  "stroke-linecap": "round",
                  "data-outline-edge": edge.id,
                }),
              );
          });
          svg.insertBefore(region, svg.firstChild);
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
        if (!mini) {
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
    }
    function renderInputs(product, o) {
      product.querySelector(".hs-family-bar").hidden = !o.terms.length;
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
            `<button type="button" class="hs-family-choice" data-choose-family="${start + j}" aria-label="Choose family ${start + j + 1}" aria-pressed="${state.family === start + j}"><span class="hs-family-label">F${start + j + 1}</span><div class="hs-mini" data-preview-family="${start + j}"></div></button>`,
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
        const neighbors = [
          0,
          state.orientation - 1,
          state.orientation,
          state.orientation + 1,
          data.orientations.length - 1,
        ].filter((index) => index >= 0 && index < data.orientations.length);
        // Navigate retained orientations by position, but display their native IDs.
        const sum = sumSequence(neighbors, (index) => {
          const id = data.orientations[index].id;
          return index === state.orientation
            ? `<label class="hs-sum-current" data-total-orientation="${id}" title="Go to orientation ID · Enter to apply, Escape to cancel">C<sub><input type="text" inputmode="numeric" data-orientation-input aria-label="Orientation ID" value="${id}" style="width:${Math.max(2, String(id).length)}ch" /></sub></label>`
            : `<button type="button" class="hs-sum-term" data-choose-orientation="${index}" data-total-orientation="${id}" aria-label="Show orientation ${id}">C<sub>${id}</sub></button>`;
        });
        product.querySelector("[data-total-sum]").innerHTML =
          `<div class="hs-sum-math" role="group" aria-label="Orientation sum"><span>C =</span>${sum}</div>`;
        const tree = o.terms.length ? factorTree(o) : null,
          shared = tree?.factors.length || 0;
        product.querySelector("[data-orientation-sum]").innerHTML =
          `<div class="hs-sum-label">${o.terms.length} ${o.terms.length === 1 ? "family" : "families"}${o.terms.length > 1 && shared ? ` · ${shared} shared ${shared === 1 ? "factor" : "factors"}` : ""}</div><div class="hs-factored-equation" aria-label="Factored denominator expression for orientation ${o.id}, with family ${state.family + 1} highlighted"><span class="hs-factored-lhs">C<sub>${o.id}</sub> =</span>${tree ? factoredMath(tree, term, selected) : "0"}</div>`;
        product.querySelector("[data-contours]").checked = state.contours;
        product.querySelector("[data-multiple]").checked = state.multiple;
        product.querySelector("[data-selection-count]").textContent =
          `${selected.length} surface${selected.length === 1 ? "" : "s"} selected`;
        product.querySelector("[data-detail]").innerHTML = selected.length
          ? selected
              .map((id) => {
                const surface = data.surfaces[id];
                return `<div class="hs-definition-row" style="--hs-color:${color(id, term)}" data-definition="${id}"><div class="hs-definition">${eta(id)} = ${definition(surface)}</div><div class="hs-sets">Inside {${surface.v.join(", ")}}</div></div>`;
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
      if (!event.target.matches("[data-orientation-input]")) return;
      if (event.key === "Escape") {
        event.preventDefault();
        event.target.value = String(data.orientations[state.orientation].id);
        event.target.blur();
      } else if (event.key === "Enter") {
        event.preventDefault();
        const value = Number(event.target.value),
          index = data.orientations.findIndex((o) => o.id === value);
        if (
          event.target.value.trim() === "" ||
          !Number.isInteger(value) ||
          index < 0
        ) {
          message = "Enter an orientation ID retained in this result.";
          messageError = true;
        } else chooseOrientation(index);
        render();
        const input = root.querySelector("[data-orientation-input]");
        input.focus({ preventScroll: true });
        input.select();
      }
    });
    root.addEventListener("focusout", (event) => {
      // Leaving an unsubmitted edit must not rerender and swallow a neighbor click.
      if (event.target.matches("[data-orientation-input]"))
        event.target.value = String(data.orientations[state.orientation].id);
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
        "[data-choose-orientation],[data-flip-edge],[data-choose-family],[data-page-step]",
      );
      if (action) {
        let focusSelector = "";
        if (action.hasAttribute("data-choose-orientation")) {
          chooseOrientation(Number(action.dataset.chooseOrientation));
          focusSelector = "[data-orientation-input]";
        } else if (action.hasAttribute("data-flip-edge")) {
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
