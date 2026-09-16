// Tydenso's pure-Typst notation layer for portable Atom render trees.
//
// The canonical generic renderer lives beside this file as `render.typ`.
// Rust supplies algebraic shape and opaque portable attachments; every Spenso
// layout decision below is ordinary, document-customizable Typst code.

#import "render.typ" as core

#let _merge(base, overlay) = {
  if type(overlay) != dictionary { panic("notation overrides must be dictionaries") }
  let result = base
  for (key, value) in overlay { result.insert(key, value) }
  result
}

#let _settings(value) = _merge((
  tensor-layout: "ports",
  with-dim: false,
  parens: true,
  commas: none,
  symbol-scripts: true,
  index-gap: 0.08em,
  factor-gap: 0.12em,
), value)

#let _short-name(name) = name.split("::").last()
#let _symbol(node) = node.at("symbol", default: (:))
#let _name(node) = _symbol(node).at("name", default: "?")
#let _node-short-name(node) = _symbol(node).at(
  "short-name",
  default: _short-name(_name(node)),
)
#let _tags(node) = _symbol(node).at("tags", default: ())
#let _has-tag(node, wanted) = {
  let short = _short-name(wanted)
  _tags(node).any(tag => tag == wanted or _short-name(tag) == short)
}
#let _kind(node) = if node.at("kind", default: none) in ("variable", "symbol", "var") {
  "variable"
} else {
  node.at("kind", default: none)
}
#let _is-function(node, name: none) = (
  type(node) == dictionary
    and _kind(node) == "function"
    and (name == none or _node-short-name(node) == name)
)
#let _is-variable(node, name: none) = (
  type(node) == dictionary
    and _kind(node) == "variable"
    and (name == none or _node-short-name(node) == name)
)

#let _record(ctx, schema, identity) = {
  if type(identity) != str { return none }
  let records = ctx.attachments.filter(attachment => (
    attachment.at("schema", default: none) == schema
      and attachment.at("identity", default: none) == bytes(identity)
  ))
  if records.len() == 0 { none } else { records.first() }
}
#let _attachment(ctx, schema, identity, versions) = {
  let record = _record(ctx, schema, identity)
  if record == none or record.at("version", default: none) not in versions {
    none
  } else {
    cbor(record.data)
  }
}

#let _builtin-symbols = (
  alpha: $alpha$, beta: $beta$, gamma: $gamma$, delta: $delta$,
  epsilon: $epsilon$, zeta: $zeta$, eta: $eta$, theta: $theta$,
  iota: $iota$, kappa: $kappa$, lambda: $lambda$, mu: $mu$,
  nu: $nu$, xi: $xi$, omicron: $omicron$, pi: $pi$, rho: $rho$,
  sigma: $sigma$, tau: $tau$, upsilon: $upsilon$, phi: $phi$,
  chi: $chi$, psi: $psi$, omega: $omega$,
  Alpha: $Alpha$, Beta: $Beta$, Gamma: $Gamma$, Delta: $Delta$,
  Epsilon: $Epsilon$, Zeta: $Zeta$, Eta: $Eta$, Theta: $Theta$,
  Iota: $Iota$, Kappa: $Kappa$, Lambda: $Lambda$, Mu: $Mu$,
  Nu: $Nu$, Xi: $Xi$, Omicron: $Omicron$, Pi: $Pi$, Rho: $Rho$,
  Sigma: $Sigma$, Tau: $Tau$, Upsilon: $Upsilon$, Phi: $Phi$,
  Chi: $Chi$, Psi: $Psi$, Omega: $Omega$,
)
#let _display-symbol(name) = {
  let builtin = _builtin-symbols.at(name, default: none)
  if builtin != none { builtin } else { math.italic(name) }
}
#let _display-literal(node) = {
  if type(node) == array and node.len() == 2 and node.first() in ("symbol", "text") {
    node.at(1)
  } else {
    none
  }
}
#let _display-node(display) = {
  if type(display) != array or display.len() == 0 { return math.upright("?") }
  let kind = display.first()
  if kind == "symbol" and display.len() == 2 {
    return _display-symbol(display.at(1))
  }
  if kind == "text" and display.len() == 2 {
    return math.upright(display.at(1))
  }
  if kind == "number" and display.len() == 2 {
    return $ #display.at(1) $
  }
  if kind == "sequence" and display.len() == 2 {
    return display.at(1).map(_display-node).join()
  }
  if kind == "list" and display.len() == 2 {
    return display.at(1).map(_display-node).join[$,$]
  }
  if kind == "attach" and display.len() == 4 {
    let base = _display-node(display.at(1))
    let top = display.at(2)
    let bottom = display.at(3)
    if top != none and bottom != none {
      return math.attach(base, t: _display-node(top), b: _display-node(bottom))
    }
    if top != none { return math.attach(base, t: _display-node(top)) }
    if bottom != none { return math.attach(base, b: _display-node(bottom)) }
    return base
  }
  if kind != "math" or display.len() != 3 { return math.upright("?") }

  let head = display.at(1)
  let raw = display.at(2)
  let args = raw.map(_display-node)
  let first = () => if args.len() == 0 { math.upright("?") } else { args.first() }

  if head == "arg" { return args.join[$,$] }
  if head == "add" { return args.join[$ + $] }
  if head == "sub" { return args.join[$ - $] }
  if head == "plus" { return $ + #first() $ }
  if head == "neg" { return $ - #first() $ }
  if head == "times" { return args.join[$ times $] }
  if head == "dot" { return args.join[$ dot $] }
  if head == "factorial" { return $ #first() ! $ }
  if head == "mul" { return args.join(h(0.12em)) }
  if head in ("()", "group", "lr") { return math.lr($ (#first()) $) }
  if head == "pow" and args.len() == 2 {
    return math.attach(args.at(0), t: args.at(1))
  }
  if head == "frac" and args.len() == 2 {
    return math.frac(args.at(0), args.at(1))
  }
  if head == "root" and args.len() == 1 { return math.root(args.first()) }
  if head == "root" and args.len() == 2 {
    return math.root(args.at(0), args.at(1))
  }
  if head in ("call", "op-call") and args.len() > 0 {
    return (args.first(), math.lr($ (#args.slice(1).join[$,$]) $)).join()
  }
  if head == "op" and raw.len() == 1 {
    let literal = _display-literal(raw.first())
    return math.op(if literal == none { "?" } else { literal })
  }
  if head == "accent" and args.len() == 2 {
    let accent = _display-literal(raw.at(1))
    return math.accent(args.first(), if accent == none { "?" } else { accent })
  }
  if head == "underline" and args.len() == 1 { return math.underline(first()) }
  if head == "overline" and args.len() == 1 { return math.overline(first()) }
  if head == "cancel" and args.len() == 1 { return math.cancel(first()) }
  if head == "scripts" and args.len() == 1 { return math.scripts(first()) }
  if head == "limits" and args.len() == 1 { return math.limits(first()) }
  if head == "class" and args.len() == 2 {
    let class = _display-literal(raw.at(0))
    return math.class(if class == none { "normal" } else { class }, args.at(1))
  }
  if head == "underbrace" { return math.underbrace(..args) }
  if head == "overbrace" { return math.overbrace(..args) }
  if head == "underbracket" { return math.underbracket(..args) }
  if head == "overbracket" { return math.overbracket(..args) }
  if head == "underparen" { return math.underparen(..args) }
  if head == "overparen" { return math.overparen(..args) }
  if head == "undershell" { return math.undershell(..args) }
  if head == "overshell" { return math.overshell(..args) }
  if head == "vec" { return math.vec(..args) }
  if head == "mat" { return math.mat(..args) }
  if head == "cases" { return math.cases(..args) }
  if head == "mid" and args.len() == 1 { return math.mid(first()) }
  if head == "stretch" and args.len() == 1 { return math.stretch(first()) }
  math.upright("?")
}

#let _payload-head(ctx) = {
  let display = _attachment(
    ctx,
    "spenso.math-display",
    ctx.identity,
    (1,),
  )
  if display == none { (ctx.default)() } else { _display-node(display) }
}

#let _visual(ctx, node) = (ctx.render-visual)(node)
#let _head(ctx, node) = (ctx.head-for)(node)

#let _unwrap-display(node) = {
  let value = node
  while _is-function(value, name: "tensor_display") and value.arguments.len() == 1 {
    value = value.arguments.first()
  }
  value
}
#let _is-marker(node, name) = _is-variable(_unwrap-display(node), name: name)
#let _chain-markers(arguments) = {
  let inputs = ()
  let outputs = ()
  for (position, argument) in arguments.enumerate() {
    if _is-marker(argument, "in") { inputs.push(position) }
    if _is-marker(argument, "out") { outputs.push(position) }
  }
  if inputs.len() != 1 or outputs.len() != 1 { return none }
  (
    input: inputs.first(),
    output: outputs.first(),
    transposed: inputs.first() > outputs.first(),
  )
}

#let _representation(node, ctx) = {
  let value = _unwrap-display(node)
  let dual = false
  if _is-function(value, name: "dind") and value.arguments.len() == 1 {
    dual = true
    value = _unwrap-display(value.arguments.first())
  }
  if not _is-function(value) or not _has-tag(value, "representation") {
    return none
  }

  let data = _attachment(ctx, "spenso.representation", _name(value), (1,))
  if type(data) != array or data.len() != 3 { return none }
  let (class, palette, row) = data
  if (
    class not in ("self-dual", "dualizable", "inline-metric")
      or row not in ("top", "bottom")
  ) {
    return none
  }
  if dual and class == "dualizable" {
    row = if row == "top" { "bottom" } else { "top" }
  }
  (
    node: value,
    identity: _name(value),
    class: class,
    palette: palette,
    row: row,
    polarity: if class == "dualizable" and not dual { "ket" } else { "bra" },
    dual: dual,
  )
}
#let _representation-value(node, ctx) = {
  let rep = _representation(node, ctx)
  if rep == none or rep.node.arguments.len() != 1 {
    none
  } else {
    (..rep, dimension: rep.node.arguments.first())
  }
}
#let _slot(node, ctx) = {
  let rep = _representation(node, ctx)
  if rep == none or rep.node.arguments.len() != 2 {
    none
  } else {
    (
      ..rep,
      dimension: rep.node.arguments.at(0),
      index: rep.node.arguments.at(1),
    )
  }
}

#let _natural-index(node) = {
  if _kind(node) != "number" { return none }
  let text = node.at("text", default: none)
  if type(text) != str or text.len() == 0 { return none }
  if not text.clusters().all(cluster => cluster in (
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
  )) { return none }
  int(text)
}
#let _palette-index(palette, index) = {
  if type(palette) != array or palette.len() == 0 or palette.first() == "numeric" {
    return none
  }
  if palette.len() != 3 or palette.first() != "cyclic" { return none }
  let start = palette.at(1)
  let labels = palette.at(2)
  if type(start) != int or type(labels) != array or labels.len() == 0 or index < start {
    return none
  }
  let offset = index - start
  let cycle = calc.floor(offset / labels.len())
  let position = offset - cycle * labels.len()
  let label = _display-node(labels.at(position))
  if cycle == 0 { label } else { math.attach(label, b: $ #cycle $) }
}
#let _slot-index(slot, ctx) = {
  let index = _natural-index(slot.index)
  if index != none {
    let display = _palette-index(slot.palette, index)
    if display != none { return display }
  }
  // Exact document heads are resolved by the common renderer before its
  // payload-head fallback, so they beat spenso.math-display overlays.
  _visual(ctx, slot.index)
}

#let _phantom(value) = context if target() == "html" {
  html.elem("mphantom", value)
} else {
  hide(value)
}

#let _row-attachment(base, columns, settings) = {
  if columns.len() == 0 { return base }
  let top = ()
  let bottom = ()
  for column in columns {
    if column.row == "bottom" {
      top.push(_phantom(column.source))
      bottom.push(column.source)
    } else {
      top.push(column.source)
      bottom.push(_phantom(column.source))
    }
  }
  math.attach(
    (base, _phantom($zws$)).join(),
    t: top.join(h(settings.index-gap)),
    b: bottom.join(h(settings.index-gap)),
  )
}
#let _qualified-index(slot, source, ctx, settings) = {
  if not settings.with-dim { return source }
  math.attach(
    source,
    t: math.attach(_head(ctx, slot.node), b: _visual(ctx, slot.dimension)),
  )
}
#let _qualified-port(rep, ctx, settings) = {
  let port = if rep.class == "inline-metric" {
    sym.square.stroked.small
  } else if rep.class == "self-dual" {
    sym.circle.stroked.small
  } else if rep.class == "dualizable" {
    if rep.dual {
      sym.triangle.r.small.stroked
    } else {
      sym.triangle.l.small.stroked
    }
  } else {
    rep.class
  }
  if not settings.with-dim { return port }
  math.attach(
    port,
    t: math.attach(_head(ctx, rep.node), b: _visual(ctx, rep.dimension)),
  )
}
#let _qualified-compact-vector(compact, ctx, settings) = {
  let source = math.bold(compact.label)
  let rep = compact.representation
  if not settings.with-dim { return source }
  math.attach(
    source,
    t: math.attach(_head(ctx, rep.node), b: _visual(ctx, rep.dimension)),
  )
}

#let _parentheses(body) = math.lr($ (#body) $)
#let _brackets(body) = math.lr($ [#body] $)
#let _tight(parts) = parts.join()
#let _bra(label) = _tight((math.upright("⟨"), label, math.upright("|")))
#let _ket(label) = _tight((math.upright("|"), label, math.upright("⟩")))
#let _ports(body, bras, kets) = _tight(
  (if bras.len() > 0 { (_bra(bras.join[$,$]),) } else { () })
    + (body,)
    + (if kets.len() > 0 { (_ket(kets.join[$,$]),) } else { () }),
)
#let _call-separator(settings) = if settings.commas {
  _tight(($,$, h(settings.index-gap)))
} else {
  h(settings.index-gap)
}
#let _call-layout(base, columns, settings) = {
  if columns.len() == 0 { return base }
  let top = columns.filter(column => column.row == "top").map(column => column.source)
  let bottom = columns.filter(column => column.row == "bottom").map(column => column.source)
  let separator = _call-separator(settings)
  let body = if bottom.len() == 0 {
    top.join(separator)
  } else if top.len() == 0 {
    _tight(($;$, h(settings.index-gap), bottom.join(separator)))
  } else {
    _tight((
      top.join(separator),
      h(settings.index-gap),
      $;$,
      h(settings.index-gap),
      bottom.join(separator),
    ))
  }
  _tight((base, _parentheses(body)))
}

#let _compact-vector(node, ctx, settings) = {
  let node = _unwrap-display(node)
  if (
    not _is-function(node)
      or not _has-tag(node, "tensor")
      or not _has-tag(node, "rank1")
  ) {
    return none
  }
  let rep = none
  let labels = ()
  for argument in node.arguments {
    if _slot(argument, ctx) != none { return none }
    let candidate = _representation-value(argument, ctx)
    if candidate == none {
      labels.push(argument)
    } else if rep != none {
      return none
    } else {
      rep = candidate
    }
  }
  if rep == none { return none }

  let label = _head(ctx, node)
  let rendered = labels.map(argument => _visual(ctx, argument))
  if settings.symbol-scripts and rendered.len() > 0 {
    label = _row-attachment(
      label,
      rendered.map(source => (source: source, row: "bottom")),
      settings,
    )
  } else if rendered.len() > 0 {
    let separator = if settings.commas { $,$ } else { h(0.15em) }
    label = _tight((label, _parentheses(rendered.join(separator))))
  }
  (label: label, representation: rep)
}

#let _render-tensor(ctx, settings) = {
  let node = ctx.node
  if _has-tag(node, "rank1") {
    let compact = _compact-vector(node, ctx, settings)
    if compact != none { return compact.label }
  }

  let columns = ()
  let ordinary = ()
  let bras = ()
  let kets = ()
  let markers = _chain-markers(node.arguments)
  for (position, argument) in node.arguments.enumerate() {
    if markers != none and position in (markers.input, markers.output) { continue }
    let slot = _slot(argument, ctx)
    if slot != none {
      columns.push((
        source: _qualified-index(slot, _slot-index(slot, ctx), ctx, settings),
        row: slot.row,
      ))
      continue
    }
    let compact = _compact-vector(argument, ctx, settings)
    if compact != none {
      let rep = compact.representation
      let source = if settings.tensor-layout == "ports" {
        _qualified-port(rep, ctx, settings)
      } else {
        _qualified-compact-vector(compact, ctx, settings)
      }
      columns.push((source: source, row: rep.row))
      if settings.tensor-layout == "ports" {
        if rep.polarity == "ket" { kets.push(compact.label) } else { bras.push(compact.label) }
      }
      continue
    }
    let source = _visual(ctx, argument)
    if settings.symbol-scripts {
      columns.push((source: source, row: "bottom"))
    } else {
      ordinary.push(source)
    }
  }

  let base = _head(ctx, node)
  if ordinary.len() > 0 {
    let separator = if settings.commas { $,$ } else { h(0.15em) }
    base = _tight((base, _parentheses(ordinary.join(separator))))
  }
  base = if settings.tensor-layout == "call" {
    _call-layout(base, columns, settings)
  } else {
    _row-attachment(base, columns, settings)
  }
  if markers != none and markers.transposed {
    base = math.attach(base, t: math.upright("T"))
  }
  if settings.tensor-layout == "ports" { _ports(base, bras, kets) } else { base }
}

// A document tag or class attached to a tensor should see the same visual
// arguments that Tydenso itself uses, rather than generic Atom calls such as
// `mink(4, 1)`. The structured `arguments` remain untouched for callers that
// need the exact tree.
#let _tensor-document-context(ctx, settings) = {
  let visual-arguments = ctx.arguments.map(argument => {
    let slot = _slot(argument, ctx)
    if slot != none {
      return _qualified-index(slot, _slot-index(slot, ctx), ctx, settings)
    }
    let compact = _compact-vector(argument, ctx, settings)
    if compact != none {
      if settings.tensor-layout == "ports" { return compact.label }
      return _qualified-compact-vector(compact, ctx, settings)
    }
    _visual(ctx, argument)
  })
  let result = ctx
  result.insert("visual-arguments", visual-arguments)
  result.insert("default", () => _render-tensor(ctx, settings))
  result
}

#let _render-gamma(ctx, settings) = {
  let node = ctx.node
  if node.arguments.len() != 3 { return _render-tensor(ctx, settings) }
  let first = node.arguments.at(0)
  let second = node.arguments.at(1)
  let lorentz = node.arguments.at(2)
  let forward = _is-marker(first, "in") and _is-marker(second, "out")
  let transposed = _is-marker(first, "out") and _is-marker(second, "in")

  if forward or transposed {
    let compact = _compact-vector(lorentz, ctx, settings)
    if compact != none {
      let body = math.cancel(compact.label)
      return if transposed { math.attach(body, t: math.upright("T")) } else { body }
    }
  }

  let selected = if forward or transposed { (lorentz,) } else { node.arguments }
  let columns = ()
  for argument in selected {
    let slot = _slot(argument, ctx)
    if slot == none {
      columns.push((source: _visual(ctx, argument), row: "bottom"))
    } else {
      columns.push((
        source: _qualified-index(slot, _slot-index(slot, ctx), ctx, settings),
        row: slot.row,
      ))
    }
  }
  let body = _row-attachment(_head(ctx, node), columns, settings)
  if transposed { math.attach(body, t: math.upright("T")) } else { body }
}

#let _same-representation(left, right) = (
  left.identity == right.identity
    and left.dimension == right.dimension
    and left.class in ("self-dual", "inline-metric")
)
#let _render-dot(ctx, settings) = {
  if ctx.arguments.len() != 2 { return (ctx.default)() }
  let left = _compact-vector(ctx.arguments.at(0), ctx, settings)
  let right = _compact-vector(ctx.arguments.at(1), ctx, settings)
  if (
    left == none
      or right == none
      or not _same-representation(left.representation, right.representation)
  ) {
    return (ctx.default)()
  }
  _tight((left.label, math.class("normal", $dot$), right.label))
}

#let _render-chain(ctx, settings) = {
  if ctx.arguments.len() < 2 { return (ctx.default)() }
  let start = ctx.arguments.at(0)
  let end = ctx.arguments.at(1)
  let body = ctx.arguments.slice(2).map(argument => _visual(ctx, argument))
    .join(h(settings.factor-gap))
  if settings.parens { body = _brackets(body) }

  let columns = ()
  let prefix = none
  let suffix = none
  let compact = _compact-vector(start, ctx, settings)
  let slot = _slot(start, ctx)
  if compact != none {
    prefix = compact.label
  } else if slot != none {
    columns.push((
      source: _qualified-index(slot, _slot-index(slot, ctx), ctx, settings),
      row: slot.row,
    ))
  } else {
    prefix = _visual(ctx, start)
  }

  compact = _compact-vector(end, ctx, settings)
  slot = _slot(end, ctx)
  if compact != none {
    suffix = compact.label
  } else if slot != none {
    columns.push((
      source: _qualified-index(slot, _slot-index(slot, ctx), ctx, settings),
      row: slot.row,
    ))
  } else {
    suffix = _visual(ctx, end)
  }

  body = _row-attachment(body, columns, settings)
  _tight(
    (if prefix == none { () } else { (_bra(prefix),) })
      + (body,)
      + (if suffix == none { () } else { (_ket(suffix),) }),
  )
}

#let _render-trace(ctx, settings) = {
  if ctx.arguments.len() == 0 { return (ctx.default)() }
  let representation = ctx.arguments.first()
  let payload = ctx.arguments.slice(1)
  let factors = none
  if payload.len() == 0 {
    factors = ()
  } else if payload.len() == 1 and _is-function(payload.first(), name: "cyclic") {
    factors = payload.first().arguments
  } else if payload.len() == 1 and _is-function(payload.first(), name: "sym") {
    factors = payload
  } else {
    return (ctx.default)()
  }

  let head = math.op("Tr")
  if settings.with-dim { head = math.attach(head, b: _visual(ctx, representation)) }
  if factors.len() == 0 { return head }
  let body = factors.map(argument => _visual(ctx, argument)).join(h(settings.factor-gap))
  _tight((head, _parentheses(body)))
}

#let _lookup-document-renderer(renderers, keys) = {
  if type(renderers) != dictionary { return none }
  for key in keys {
    if type(key) == str and key in renderers { return renderers.at(key) }
  }
  none
}
#let _invoke-document-renderer(renderer, ctx) = {
  let visual = if type(renderer) == function { renderer(ctx) } else { renderer }
  if type(visual) != content {
    panic("a Tydenso notation renderer must return Typst content")
  }
  if repr(visual.func()) == "equation" and "body" in visual.fields() {
    visual.fields().body
  } else {
    visual
  }
}

#let _default-heads = (
  "spenso::gamma": $gamma$,
  "spenso::gamma0": $gamma_0$,
  "spenso::gamma5": $gamma_5$,
  "spenso::projp": $ℙ^+$,
  "spenso::projm": $ℙ^-$,
)

/// Build Tydenso's document-side notation.
///
/// Exact `heads` and `calls` win over package defaults. Call, tag, and class
/// closures receive the common renderer's full context dictionary, including
/// `arguments`, `visual-arguments`, `attachments`, and `render-visual`.
/// Tensor tag and class renderers receive slot-aware visual arguments such as
/// $mu$ rather than the underlying representation call `mink(4, 1)`.
#let notation(
  settings: (:),
  tensor-layout: none,
  with-dim: none,
  parens: none,
  commas: none,
  symbol-scripts: none,
  index-gap: none,
  factor-gap: none,
  heads: (:),
  calls: (:),
  tags: (:),
  classes: (:),
) = {
  let settings = _settings(settings)
  for (key, value) in (
    tensor-layout: tensor-layout,
    with-dim: with-dim,
    parens: parens,
    commas: commas,
    symbol-scripts: symbol-scripts,
    index-gap: index-gap,
    factor-gap: factor-gap,
  ) {
    if value != none { settings.insert(key, value) }
  }
  if settings.tensor-layout not in ("ports", "schoonschip", "call") {
    panic("tensor-layout must be \"ports\", \"schoonschip\", or \"call\"")
  }
  if settings.commas == none {
    settings.insert("commas", settings.tensor-layout == "call")
  }
  let tensor = ctx => if ctx.kind == "function" {
    // The common dispatcher may encounter the broad package tensor tag before
    // a later document tag or class on the same symbol. Re-check the complete
    // document layer here so either kind of document override still wins over
    // this package default, independent of Symbolica's tag ordering.
    let renderer = _lookup-document-renderer(tags, ctx.tags)
    if renderer == none {
      renderer = _lookup-document-renderer(classes, ctx.classes)
    }
    if renderer == none {
      _render-tensor(ctx, settings)
    } else {
      _invoke-document-renderer(renderer, _tensor-document-context(ctx, settings))
    }
  } else {
    (ctx.default)()
  }
  let defaults = core.notation(
    heads: _default-heads,
    calls: (
      "spenso::gamma": ctx => _render-gamma(ctx, settings),
      "spenso::dot": ctx => _render-dot(ctx, settings),
      "spenso::chain": ctx => _render-chain(ctx, settings),
      "spenso::trace": ctx => _render-trace(ctx, settings),
    ),
    tags: (
      tensor: tensor,
      "spenso::tensor": tensor,
    ),
    fallback-head: _payload-head,
  )
  let document = core.notation(
    heads: heads,
    calls: calls,
    tags: tags,
    classes: classes,
  )
  core.merge-notation(defaults, document)
}

/// Construct the built-in Tydenso notation with optional layout settings.
#let default-notation(settings: (:)) = notation(settings: settings)

/// Merge notation layers from least to most specific.
#let merge-notation(..layers) = core.merge-notation(..layers)

/// Render a portable Atom render tree with Tydenso's notation.
///
/// The integration callback is invoked as `annotate(atom, visual, node)`.
/// Custom renderers return only visual content; the common renderer applies the
/// exact full-call metadata afterward.
#let render(
  tree,
  notation: default-notation(),
  annotate: none,
  block: false,
) = {
  let body = core.render-tree(tree, notation: notation, annotate: annotate)
  if block { math.equation(body, block: true) } else { body }
}
