// Synced from symbolica-typst-plugin cf7b9fe59ba3fc7fe9ea70b875f5660ee2addd7a.
// Local extensions preserve tensor power grouping, complex signs, and rational coefficients.
// Generic, document-side rendering for a `symbolica` Atom render tree.
// Rust owns algebra and exact payloads; Typst owns presentation.

#let _empty-notation = (
  heads: (:),
  calls: (:),
  tags: (:),
  classes: (:),
  fallback-head: none,
  fallback-call: none,
  fallback-node: none,
  fallback: none,
)

/// Build one presentation layer for generic Atom render trees.
///
/// `heads` changes an exact namespaced variable or function head. `calls`
/// changes a complete exact function call. `tags` and `classes` are fallback
/// maps. Each renderer receives one context dictionary containing the
/// structured `node`, `symbol`, `arguments`, document `attachments`, and
/// rendering helpers. `power-base` marks calls used as the base of a power,
/// so renderers that display them as products can preserve grouping.
/// In a whole-call renderer, prefer `visual-arguments` and
/// `render-visual`; the framework then adds one exact annotation around the
/// finished visual. Use `(ctx.render)(node)` explicitly when a custom layout
/// should retain nested exact metadata. Typst calls functions stored in
/// dictionaries with parentheses around the field access, for example
/// `(ctx.render-visual)(node)` or `(ctx.default)()`.
///
/// ```typ
/// #let display = notation(
///   heads: ("model::x": ctx => $xi$),
///   calls: ("model::f": ctx => {
///     let (argument,) = ctx.visual-arguments
///     $cal(F)(#argument)$
///   }),
/// )
/// ```
#let notation(
  heads: (:),
  calls: (:),
  tags: (:),
  classes: (:),
  fallback-head: none,
  fallback-call: none,
  fallback-node: none,
  fallback: none,
) = (
  heads: heads,
  calls: calls,
  tags: tags,
  classes: classes,
  fallback-head: fallback-head,
  fallback-call: fallback-call,
  fallback-node: fallback-node,
  fallback: fallback,
)

#let _field(value, names, default: none) = {
  if type(value) != dictionary { return default }
  for name in names {
    if name in value { return value.at(name) }
  }
  default
}

#let _merge-dictionaries(base, overlay) = {
  let result = if type(base) == dictionary { base } else { (:) }
  if type(overlay) == dictionary {
    for (key, value) in overlay { result.insert(key, value) }
  }
  result
}

/// Merge presentation layers from least to most specific.
///
/// This implements package defaults, payload-local overlays, and explicit
/// document overrides without putting renderer functions in Atom metadata.
#let merge-notation(..layers) = {
  if layers.named().len() > 0 {
    panic("merge-notation accepts only positional notation layers")
  }
  let result = _empty-notation
  for layer in layers.pos() {
    if layer == none { continue }
    if type(layer) != dictionary { panic("a notation layer must be a dictionary") }
    for key in ("heads", "calls", "tags", "classes") {
      result.insert(key, _merge-dictionaries(
        result.at(key),
        _field(layer, (key,), default: (:)),
      ))
    }
    for key in ("fallback-head", "fallback-call", "fallback-node", "fallback") {
      let value = _field(layer, (key,), default: none)
      if value != none { result.insert(key, value) }
    }
  }
  result
}

/// Construct the conventional exact notation key for a namespaced head.
#let head-key(namespace, name) = {
  if namespace == none or namespace == "" { name } else { namespace + "::" + name }
}

#let _as-array(value) = {
  if value == none { () }
  else if type(value) == array { value }
  else { (value,) }
}

#let _kind(node) = {
  let kind = _field(node, ("kind", "type", "head"), default: none)
  if kind == "var" or kind == "symbol" { "variable" }
  else if kind == "call" { "function" }
  else if kind == "pow" { "power" }
  else if kind == "mul" { "product" }
  else if kind == "add" { "sum" }
  else { kind }
}

#let _symbol(node) = {
  let value = _field(node, ("symbol", "descriptor", "head-symbol", "head_symbol"), default: (:))
  if type(value) == dictionary { value }
  else if type(value) == str { (name: value, short-name: value) }
  else { (:) }
}

#let _identity(node) = {
  let symbol = _symbol(node)
  let exact = _field(
    symbol,
    ("name", "identity", "qualified-name", "qualified_name"),
    default: _field(node, ("identity", "qualified-name", "qualified_name"), default: none),
  )
  if exact != none { return exact }
  let namespace = _field(symbol, ("namespace",), default: _field(node, ("namespace",), default: none))
  let name = _field(
    symbol,
    ("short-name", "short_name"),
    default: _field(node, ("name", "short-name", "short_name"), default: none),
  )
  if name == none { none } else { head-key(namespace, name) }
}

#let _short-name(node) = _field(
  _symbol(node),
  ("short-name", "short_name"),
  default: _field(_symbol(node), ("name",), default: _field(node, ("name",), default: none)),
)

#let _namespace(node) = _field(
  _symbol(node),
  ("namespace",),
  default: _field(node, ("namespace",), default: none),
)

#let _attributes(node) = {
  let value = _field(
    _symbol(node),
    ("attributes", "attrs"),
    default: _field(node, ("attributes", "attrs"), default: ()),
  )
  if type(value) in (array, dictionary) { value } else { () }
}

#let _tags(node) = _as-array(_field(
  _symbol(node),
  ("tags", "tag"),
  default: _field(node, ("tags", "tag"), default: ()),
))

#let _classes(node) = {
  let attributes = _attributes(node)
  let direct = _field(
    _symbol(node),
    ("classes", "class"),
    default: _field(node, ("classes", "class"), default: none),
  )
  let result = _as-array(direct)
  let attribute-classes = if type(attributes) == dictionary {
    _as-array(_field(attributes, ("classes", "class"), default: none))
  } else {
    attributes
  }
  for class in attribute-classes {
    if class not in result { result.push(class) }
  }
  result
}

#let _arguments(node) = _as-array(_field(
  node,
  ("arguments", "args", "children"),
  default: (),
))

#let _atom(node) = _field(
  node,
  ("atom", "exact-atom", "exact_atom", "payload", "exact"),
  default: none,
)

#let _head-atom(node) = {
  let atom = _field(node, ("head-atom", "head_atom", "symbol-atom", "symbol_atom"), default: none)
  if atom != none { return atom }
  _field(_symbol(node), ("atom", "exact-atom", "exact_atom"), default: none)
}

#let _source(node) = {
  let source = _field(node, ("source", "typst-source", "typst_source"), default: none)
  if source != none { return source }
  _field(_symbol(node), ("source", "typst-source", "typst_source"), default: none)
}

#let _math-body(value) = {
  if type(value) == content and repr(value.func()) == "equation" and "body" in value.fields() {
    value.fields().body
  } else {
    value
  }
}

// Source is emitted by the plugin from built-in Typst math syntax only. The
// explicitly empty scope makes it impossible to depend on caller bindings.
#let _source-content(source) = {
  if type(source) != str { panic("an Atom render-tree leaf is missing Typst source") }
  _math-body(eval(source, mode: "math", scope: (:)))
}

#let _default-source(node) = {
  let source = _source(node)
  if source != none { return _source-content(source) }
  let name = _short-name(node)
  if type(name) != str { panic("an Atom render-tree symbol is missing a name") }
  if name.len() == 1 { _source-content(name) } else { _source-content(repr(name)) }
}

#let _parenthesize(value) = _math-body($ (#value) $)

#let _lookup(map, keys) = {
  if type(map) != dictionary { return none }
  for key in keys {
    if type(key) == str and key in map { return map.at(key) }
  }
  none
}

#let _invoke(renderer, ctx) = {
  let value = if type(renderer) == function { renderer(ctx) } else { renderer }
  if type(value) != content { panic("a notation renderer must return Typst content") }
  _math-body(value)
}

#let _annotate(callback, atom, visual, node) = {
  if callback == none or atom == none { return visual }
  if type(callback) != function { panic("annotate must be a function") }
  let value = callback(atom, visual, node)
  if type(value) != content { panic("annotate must return Typst content") }
  _math-body(value)
}

#let _tag-or-class-renderer(config, node) = {
  let renderer = _lookup(config.tags, _tags(node))
  if renderer != none { return renderer }
  _lookup(config.classes, _classes(node))
}

#let _negative-number(node) = {
  if _kind(node) != "number" { return none }
  let source = _source(node)
  if type(source) != str or not source.starts-with("-") { return none }
  // A complex coefficient's leading minus belongs only to its real part.
  if source.contains("𝑖") { return none }
  let positive = node
  positive.insert("source", source.slice(1))
  let text = _field(positive, ("text", "value"), default: none)
  if type(text) == str and text.starts-with("-") {
    if "text" in positive { positive.insert("text", text.slice(1)) }
    if "value" in positive { positive.insert("value", text.slice(1)) }
  }
  positive
}

// Classify exact real rational exponents without converting their digits to
// Typst integers (which would truncate Symbolica's arbitrary precision range).
#let _rational-exponent(node) = {
  if _kind(node) != "number" { return none }
  let source = _source(node)
  if type(source) != str { return none }
  let negative = source.starts-with("-")
  let magnitude = if negative { source.slice(1) } else { source }
  let parts = magnitude.split("/")
  if parts.len() not in (1, 2) or not parts.all(part => part.match(regex("^[0-9]+$")) != none) {
    return none
  }
  (
    negative: negative,
    numerator: parts.first(),
    denominator: if parts.len() == 2 { parts.last() } else { "1" },
    magnitude: (kind: "number", source: magnitude),
  )
}

// Division is stored as multiplication by negative rational powers. Only
// change presentation: keep the original base nodes and their exact metadata.
#let _denominator-factor(node) = {
  if _kind(node) != "power" { return none }
  let exponent = _rational-exponent(_field(node, ("exponent", "exp", "rhs", "right")))
  if exponent == none or not exponent.negative { return none }
  let base = _field(node, ("base", "lhs", "left"))
  if exponent.numerator == "1" and exponent.denominator == "1" { return base }
  let positive = (kind: "power", base: base, exponent: exponent.magnitude)
  // The displayed denominator is the reciprocal of the original factor.
  // Its metadata must describe the positive power, not the negative one.
  let reciprocal-atom = _field(node, ("reciprocal-atom",), default: none)
  if reciprocal-atom != none { positive.insert("atom", reciprocal-atom) }
  positive
}

#let _is-fraction-product(node) = {
  _kind(node) == "product" and _as-array(_field(
    node, ("factors", "arguments", "args", "children"), default: (),
  )).any(factor => _denominator-factor(factor) != none)
}

#let _split-sign(node) = {
  let positive = _negative-number(node)
  if positive != none { return (negative: true, node: positive) }
  if _kind(node) == "product" {
    let factors = _as-array(_field(node, ("factors", "arguments", "args", "children"), default: ()))
    if factors.len() > 0 {
      let first = _negative-number(factors.first())
      if first != none {
        let absolute = node
        let absolute-factors = factors
        let first-source = _source(first)
        let first-text = _field(first, ("text", "value"), default: none)
        if first-source == "1" or first-text == "1" {
          absolute-factors = absolute-factors.slice(1)
        } else {
          absolute-factors.at(0) = first
        }
        absolute.insert("factors", absolute-factors)
        return (negative: true, node: absolute)
      }
    }
  }
  (negative: false, node: node)
}

/// Render a generic Atom tree as Typst math content.
///
/// `annotate` is called as `annotate(atom, visual, node)`. It runs after every
/// custom renderer, so a renderer cannot accidentally discard the exact Atom.
/// Numbers remain ordinary, unannotated Typst math. Generic calls annotate only
/// their head and retain visible recursive structure; a customized whole call
/// receives one exact full-call annotation instead.
#let render-tree(tree, notation: _empty-notation, annotate: none) = {
  let config = merge-notation(notation)
  let envelope-kind = _field(tree, ("kind",), default: none)
  let root = if envelope-kind == "atom-render-tree" {
    _field(tree, ("root", "tree", "node"), default: none)
  } else {
    tree
  }
  if root == none { panic("the Atom render tree is missing its root node") }
  let attachments = if envelope-kind == "atom-render-tree" {
    _as-array(_field(tree, ("attachments",), default: ()))
  } else {
    ()
  }

  let render-node(node, exact: true, roots: true, power-base: false) = {
    if type(node) != dictionary { panic("an Atom render-tree node must be a dictionary") }
    let kind = _kind(node)

    // Resolve just a symbol/function head, including exact document overrides,
    // without rendering or annotating its arguments.
    let render-head-visual(target) = {
      let default = () => _default-source(target)
      let ctx = (
        kind: "head",
        node: target,
        symbol: _symbol(target),
        identity: _identity(target),
        namespace: _namespace(target),
        name: _short-name(target),
        tags: _tags(target),
        classes: _classes(target),
        attributes: _attributes(target),
        arguments: _arguments(target),
        visual-arguments: (),
        attachments: attachments,
        head: none,
        render: child => render-node(child),
        render-visual: child => render-node(child, exact: false),
        head-for: render-head-visual,
        default: default,
        parenthesize: _parenthesize,
      )
      let renderer = _lookup(config.heads, (_identity(target),))
      if renderer == none { renderer = config.fallback-head }
      if renderer == none { default() } else { _invoke(renderer, ctx) }
    }

    let make-context(
      semantic-kind,
      arguments: (),
      visual-arguments: (),
      default: none,
      head: none,
    ) = (
      kind: semantic-kind,
      power-base: power-base,
      node: node,
      symbol: _symbol(node),
      identity: _identity(node),
      namespace: _namespace(node),
      name: _short-name(node),
      tags: _tags(node),
      classes: _classes(node),
      attributes: _attributes(node),
      arguments: arguments,
      visual-arguments: visual-arguments,
      attachments: attachments,
      head: head,
      render: child => render-node(child),
      render-visual: child => render-node(child, exact: false),
      head-for: render-head-visual,
      default: default,
      parenthesize: _parenthesize,
    )

    let render-head(target, attach: true) = {
      let visual = render-head-visual(target)
      if attach and exact {
        // The annotation payload is the variable Atom for the function head,
        // not the complete call. Give callbacks a matching semantic node so
        // metadata never describes a function call while carrying head bytes.
        let head-atom = _head-atom(target)
        let head-node = (
          kind: "function-head",
          symbol: _symbol(target),
          source: _source(target),
          atom: head-atom,
        )
        _annotate(annotate, head-atom, visual, head-node)
      } else {
        visual
      }
    }

    if kind == "number" {
      return _source-content(_field(node, ("source", "text", "value"), default: none))
    }

    if kind == "variable" {
      let default = () => _default-source(node)
      let ctx = make-context("variable", default: default)
      let renderer = _lookup(config.heads, (_identity(node),))
      if renderer == none { renderer = _tag-or-class-renderer(config, node) }
      if renderer == none { renderer = config.fallback-head }
      if renderer == none { renderer = config.fallback }
      let visual = if renderer == none { default() } else { _invoke(renderer, ctx) }
      return if exact { _annotate(annotate, _atom(node), visual, node) } else { visual }
    }

    if kind == "function" {
      let arguments = _arguments(node)
      let visual-arguments = arguments.map(argument => render-node(argument, exact: false))
      let head = render-head(node, attach: false)
      let default = () => math.op(head) + _parenthesize(visual-arguments.join(_math-body($, $)))
      let ctx = make-context(
        "function",
        arguments: arguments,
        visual-arguments: visual-arguments,
        default: default,
        head: head,
      )
      let renderer = _lookup(config.calls, (_identity(node),))
      if renderer == none { renderer = _tag-or-class-renderer(config, node) }
      if renderer == none { renderer = config.fallback-call }
      if renderer == none { renderer = config.fallback }
      if renderer == none {
        let structural-arguments = if exact {
          arguments.map(argument => render-node(argument))
        } else {
          visual-arguments
        }
        return math.op(render-head(node, attach: exact)) + _parenthesize(structural-arguments.join(_math-body($, $)))
      }
      let visual = _invoke(renderer, ctx)
      return if exact { _annotate(annotate, _atom(node), visual, node) } else { visual }
    }

    if kind == "power" {
      let base-node = _field(node, ("base", "lhs", "left"), default: none)
      let exponent-node = _field(node, ("exponent", "exp", "rhs", "right"), default: none)
      if base-node == none or exponent-node == none { panic("a power node needs a base and exponent") }
      let exponent = _rational-exponent(exponent-node)
      let reciprocal = exponent != none and exponent.negative
      let radical = roots and exponent != none and exponent.denominator != "1"
      let unit = exponent != none and exponent.numerator == "1" and exponent.denominator == "1"
      // Metadata inside an attachment base prevents Typst from aligning a new
      // superscript with existing subscripts. In this case render the complete
      // power first and put its exact metadata outside the attachment.
      let detached-base = if exact and _atom(node) != none and _kind(base-node) == "variable" {
        render-node(base-node, exact: false, power-base: not radical and not unit)
      } else { none }
      // Source syntax puts a superscript on the last item of a sequence,
      // such as the closing argument group in a literal label `h'(c)`.
      // Keep custom sequence notation's existing whole-base convention.
      let sequence-power = detached-base != none and repr(detached-base.func()) == "sequence" and (
        detached-base == _default-source(base-node)
      )
      let merge-attachments = detached-base != none and (
        repr(detached-base.func()) == "attach" or sequence-power
      )
      let base = if detached-base == none {
        render-node(base-node, exact: exact, power-base: not radical and not unit)
      } else if merge-attachments {
        detached-base
      } else {
        _annotate(annotate, _atom(base-node), detached-base, base-node)
      }
      if _source(exponent-node) == "1" { return base }
      let visual = if radical {
        // Symbolica prints standalone fractional powers as roots. Powers moved
        // into a product's denominator retain their positive exponent instead.
        let radical = if exponent.denominator == "2" {
          math.sqrt(base)
        } else {
          math.root(_source-content(exponent.denominator), base)
        }
        if exponent.numerator == "1" { radical }
        else { math.attach(radical, t: _source-content(exponent.numerator)) }
      } else if reciprocal and exponent.numerator == "1" and exponent.denominator == "1" {
        base
      } else {
        if _kind(base-node) in ("sum", "product", "power") or (
          _kind(base-node) == "number" and (
            _negative-number(base-node) != none or (
              type(_source(base-node)) == str and _source(base-node).contains("/")
            )
          )
        ) { base = _parenthesize(base) }
        let power = if reciprocal { exponent.magnitude } else { exponent-node }
        let top = render-node(power, exact: exact)
        if sequence-power and base.children.len() > 0 {
          let parts = base.children
          parts.at(parts.len() - 1) = math.attach(parts.last(), t: top)
          parts.join()
        } else {
          math.attach(base, t: top)
        }
      }
      let visual = if reciprocal { math.frac(_source-content("1"), visual) } else { visual }
      return if merge-attachments { _annotate(annotate, _atom(node), visual, node) } else { visual }
    }

    if kind == "product" {
      let factors = _as-array(_field(node, ("factors", "arguments", "args", "children"), default: ()))
      let numerator = ()
      let denominator = ()
      for factor in factors {
        if _kind(factor) == "number" {
          let source = _source(factor)
          let rational = if type(source) == str { source.match(regex("^(-?[0-9]+)/([0-9]+)$")) }
          if rational != none {
            let (num, denom) = rational.captures
            if num != "1" { numerator.push((kind: "number", source: num)) }
            if denom != "1" { denominator.push((kind: "number", source: denom)) }
            continue
          }
        }
        let below = _denominator-factor(factor)
        if below == none { numerator.push(factor) } else { denominator.push(below) }
      }
      let render-factors(factors, roots: true) = {
        if factors.len() == 0 { return _source-content("1") }
        let minus = factors.len() > 1 and _kind(factors.first()) == "number" and _source(factors.first()) == "-1"
        let visible = if minus { factors.slice(1) } else { factors }
        let body = visible.map(factor => {
          let visual = render-node(factor, exact: exact, roots: roots)
          if _kind(factor) == "sum" and factors.len() > 1 { _parenthesize(visual) } else { visual }
        }).join([ ])
        if minus { _math-body($-#body$) } else { body }
      }
      let above = render-factors(numerator)
      if denominator.len() == 0 { return above }
      return math.frac(above, render-factors(denominator, roots: false))
    }

    if kind == "sum" {
      let terms = _as-array(_field(node, ("terms", "arguments", "args", "children"), default: ()))
      if terms.len() == 0 { return _source-content("0") }
      let result = none
      for (index, term) in terms.enumerate() {
        // The printer keeps a fractional term's sign inside its numerator.
        let signed = if _is-fraction-product(term) {
          (negative: false, node: term)
        } else { _split-sign(term) }
        let visual = render-node(signed.node, exact: exact)
        if index == 0 {
          result = if signed.negative { _math-body($-#visual$) } else { visual }
        } else if signed.negative {
          result += _math-body($-#visual$)
        } else {
          result += _math-body($+#visual$)
        }
      }
      return result
    }

    let default = () => panic("unsupported Atom render-tree node kind: " + repr(kind))
    let ctx = make-context(kind, default: default)
    let renderer = config.fallback-node
    if renderer == none { renderer = config.fallback }
    if renderer == none { return default() }
    let visual = _invoke(renderer, ctx)
    if exact { _annotate(annotate, _atom(node), visual, node) } else { visual }
  }

  render-node(root)
}
