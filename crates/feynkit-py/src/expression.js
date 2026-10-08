(() => {
  function factorTerms(entries) {
    const remaining = entries.map((entry) => ({
        ...entry,
        factors: [...entry.factors],
      })),
      factors = [];
    // Extract a multiset intersection, preserving powers and original term IDs.
    for (const id of entries[0].factors)
      if (remaining.every((entry) => entry.factors.includes(id))) {
        factors.push(id);
        remaining.forEach((entry) =>
          entry.factors.splice(entry.factors.indexOf(id), 1),
        );
      }
    const node = {
      factors,
      terms: entries.map((entry) => entry.id),
      children: [],
    };
    if (entries.length === 1) return node;
    const groups = new Map();
    remaining.forEach((entry) => {
      const key = entry.factors.length
        ? `factor-${entry.factors[0]}`
        : `term-${entry.id}`;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(entry);
    });
    node.children = [...groups.values()].map(factorTerms);
    return node;
  }

  for (const root of document.querySelectorAll(".feynkit-expression")) {
    if (root.expressionDisplay) continue;
    let sum;
    const mathCache = new Map(),
      mathStyles = new Set();
    function math(value) {
      if (Array.isArray(value)) return value.map(math);
      if (value && typeof value === "object")
        return Object.fromEntries(
          Object.entries(value).map(([key, html]) => [key, math(html)]),
        );
      if (typeof value !== "string") return value;
      if (!mathCache.has(value)) {
        const template = document.createElement("template");
        template.innerHTML = value;
        // Native fragments are standalone. Retain their styles once per display,
        // outside the expressions replaced during orientation/surface selection.
        for (const style of template.content.querySelectorAll("style")) {
          style.remove();
          if (!mathStyles.has(style.textContent)) {
            mathStyles.add(style.textContent);
            root.append(style);
          }
        }
        mathCache.set(value, template.innerHTML);
      }
      return mathCache.get(value);
    }
    const input = () => root.querySelector("[data-sum-input]");
    const resize = (field) => {
      field.style.width = `${Math.max(2, field.value.length)}ch`;
    };
    const restore = (field) => {
      field.value = sum.ids[sum.current];
      resize(field);
    };
    const focus = () => input()?.focus({ preventScroll: true });
    root.expressionDisplay = {
      factorTerms,
      math,
      powers(factors) {
        const counts = new Map();
        for (const factor of factors) counts.set(factor, (counts.get(factor) || 0) + 1);
        return [...counts];
      },
      energyFactor(html, id, power) {
        return `<span class="hs-factored-factor hs-energy" data-factor-id="energy-${id}" data-factor-power="${power}">${power > 1 ? "(" : ""}2${html}${power > 1 ? `)<sup data-power-label>${power}</sup>` : ""}<span class="hs-factor-spacer"></span></span>`;
      },
      fraction(numerator, denominator, attributes = "") {
        // mtext retains interactive HTML; mfrac supplies the mathematical axis,
        // independent of numerator height, surface swatches and family captions.
        return `<math class="hs-factored-fraction" ${attributes}><mfrac><mtext>${numerator}</mtext><mtext><span class="hs-factored-denominator">${denominator}</span></mtext></mfrac></math>`;
      },
      sum(options) {
        sum = options;
        // Navigate retained positions, but display and accept their native IDs.
        const { ids, current, symbol, total, label } = sum,
          indices = [
            ...new Set([0, current - 1, current, current + 1, ids.length - 1]),
          ]
            .filter((i) => i >= 0 && i < ids.length)
            .sort((a, b) => a - b),
          parts = [];
        indices.forEach((index, j) => {
          if (j) {
            parts.push("<span>+</span>");
            if (index > indices[j - 1] + 1)
              parts.push("<span>⋯</span><span>+</span>");
          }
          const id = ids[index];
          parts.push(
            index === current
              ? `<label class="hs-sum-current" data-sum-id="${id}" title="Go to ${label.toLowerCase()} ID · Enter to apply, Escape to cancel">${symbol}<sub><input type="text" inputmode="numeric" data-sum-input aria-label="${label} ID" value="${id}" style="width:${Math.max(2, String(id).length)}ch" /></sub></label>`
              : `<button type="button" class="hs-sum-term" data-sum-index="${index}" data-sum-id="${id}" aria-label="Show ${label.toLowerCase()} ${id}">${symbol}<sub>${id}</sub></button>`,
          );
        });
        return `<div class="hs-sum-math" role="group" aria-label="${label} sum"><span>${total} =</span>${parts.join("")}</div>`;
      },
    };
    root.addEventListener("click", (event) => {
      const term = event.target.closest("[data-sum-index]");
      if (!term) return;
      sum.onSelect(Number(term.dataset.sumIndex));
      focus();
    });
    root.addEventListener("keydown", (event) => {
      if (!event.target.matches("[data-sum-input]")) return;
      if (event.key === "Escape") {
        event.preventDefault();
        restore(event.target);
        event.target.blur();
      } else if (event.key === "Enter") {
        event.preventDefault();
        const value = event.target.value.trim(),
          id = Number(value),
          index = sum.ids.indexOf(id);
        if (!value || !Number.isInteger(id) || index < 0) {
          restore(event.target);
          sum.onInvalid();
        } else sum.onSelect(index);
        focus();
        input()?.select();
      }
    });
    root.addEventListener("input", (event) => {
      if (event.target.matches("[data-sum-input]")) resize(event.target);
    });
    root.addEventListener("focusout", (event) => {
      // Restore without rerendering: a pending neighboring click must survive.
      if (event.target.matches("[data-sum-input]")) restore(event.target);
    });
  }
})();
