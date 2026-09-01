#import "@preview/parsely:0.1.0"
#import "notation.typ" as tensor-notation

#let _default-grammar = (
  arg: (infix: $,$, assoc: true, prec: 4),
  add: (infix: $+$, prec: 1, assoc: true),
  sub: (infix: $-$, prec: 1),
  plus: (prefix: $+$, prec: 2),
  neg: (prefix: $-$, prec: 2),
  times: (infix: $times$, prec: 2),
  dot: (infix: $dot$, prec: 2),
  factorial: (postfix: $#parsely.tight !$, prec: 3),
  semantic-metadata: (postfix: metadata, prec: 5),
  mul: (infix: $$, prec: 2.5, assoc: true),
  "()": (match: $(#parsely.slot("expr*"))$),
  pow: (match: $#parsely.slot("base")^#parsely.slot("exp")$),
  attach: math.attach,
  accent: math.accent,
  cancel: math.cancel,
  cases: math.cases,
  class: math.class,
  frac: math.frac,
  limits: math.limits,
  lr: math.lr,
  mat: math.mat,
  mid: math.mid,
  overbrace: math.overbrace,
  overbracket: math.overbracket,
  overline: math.overline,
  overparen: math.overparen,
  overshell: math.overshell,
  root: math.root,
  scripts: math.scripts,
  stretch: math.stretch,
  underbrace: math.underbrace,
  underbracket: math.underbracket,
  underline: math.underline,
  underparen: math.underparen,
  undershell: math.undershell,
  vec: math.vec,
  op-call: (match: $op(#parsely.slot("op"))(#parsely.slot("args*"))$),
  call: (match: $#parsely.slot("fn") #parsely.tight (#parsely.slot("body*"))$),
  op: math.op,
)
#let _typst-math = math

#let _leaf(value) = {
  if type(value) == dictionary {
    value
  } else if type(value) == array {
    value.map(_leaf)
  } else if type(value) == str {
    value
  } else if type(value) == content {
    let fields = value.fields()
    if "text" in fields {
      fields.text
    } else if "children" in fields {
      fields.children.map(_leaf).join("")
    } else {
      repr(value)
    }
  } else {
    repr(value)
  }
}

#let _node-to-ast(node) = (
  head: node.head,
  args: node.args,
  slots: node.slots,
)

#let _with-semantic-metadata(grammar) = {
  let enhanced = (
    semantic-metadata: (postfix: metadata, prec: 5),
  )
  for (name, rule) in grammar {
    if name != "semantic-metadata" { enhanced.insert(name, rule) }
  }
  enhanced
}

#let _is-space(value) = {
  if type(value) == str { return value.trim() == "" }
  if type(value) == content {
    if repr(value.func()) == "space" { return true }
    if repr(value.func()) == "symbol" and "text" in value.fields() {
      return value.fields().text.trim() == ""
    }
  }
  false
}

#let _content-positional-fields = (
  attach: ("base",),
  equation: ("body",),
  frac: ("num", "denom"),
  lr: ("body",),
  root: ("index", "radicand"),
)

#let _call-content(fn, fields) = {
  let kind = repr(fn)
  if kind == "sequence" and "children" in fields {
    return fields.children.join()
  }

  let positional = ()
  for field in _content-positional-fields.at(kind, default: ()) {
    if field in fields { positional.push(fields.remove(field)) }
  }
  fn(..positional, ..fields)
}

#let _trim-math(value) = {
  let trim-array(values) = {
    let values = values.map(_trim-math)
    while values.len() > 0 and _is-space(values.first()) {
      values = values.slice(1)
    }
    while values.len() > 0 and _is-space(values.last()) {
      values = values.slice(0, values.len() - 1)
    }
    values
  }

  if type(value) == array { return trim-array(value) }
  if type(value) != content { return value }

  let kind = repr(value.func())
  if kind != "sequence" and kind not in _content-positional-fields {
    return value
  }

  let fields = (:)
  for (key, field) in value.fields() {
    fields.insert(key, _trim-math(field))
  }
  _call-content(value.func(), fields)
}

#let _ast-bytes(equation, grammar) = {
  let parsed = parsely.parse(_trim-math(equation), _with-semantic-metadata(grammar))
  let tree = parsely.walk(parsed.tree, post: _node-to-ast, leaf: _leaf)
  cbor.encode(tree)
}

#let _display-leaf(value) = {
  if type(value) == content and repr(value.func()) == "text" {
    return (
      head: "text",
      args: (value.fields().at("text", default: ""),),
      slots: (:),
    )
  }
  _leaf(value)
}

#let _display-ast-bytes(equation, grammar) = {
  let parsed = parsely.parse(_trim-math(equation), _with-semantic-metadata(grammar))
  let tree = parsely.walk(parsed.tree, post: _node-to-ast, leaf: _display-leaf)
  cbor.encode(tree)
}

#let _from-math(engine, equation, grammar: none, namespace: none) = {
  let grammar = if grammar == none { engine.grammar } else { grammar }
  let namespace = if namespace == none { engine.namespace } else { namespace }
  engine.plugin.from_ast(_ast-bytes(equation, grammar), cbor.encode(namespace))
}
#let _namespace(engine, namespace) = if namespace == none { engine.namespace } else { namespace }

#let _decompress-bundled(path) = plugin("tydenso-inflate.wasm").decompress(
  read(path, encoding: none),
)
#let _bundled-plugin() = plugin(_decompress-bundled("tydenso.wasm.zlib"))

#let _default-representation-index-row(name, namespace) = if (
  name == "spenso::bis" or (name == "bis" and namespace == "spenso")
) { "bottom" } else { "top" }

#let _builtin-representation-index-palettes = (
  "spenso::mink": ("mu", "nu", "rho", "sigma"),
  "spenso::euc": ("i", "j", "k", "l"),
  "spenso::lor": ("mu", "nu", "rho", "sigma"),
  "spenso::bis": ("a", "b", "c", "d"),
  "spenso::spf": ("alpha", "beta", "gamma", "delta"),
  "spenso::cof": ("i", "j", "k", "l"),
  "spenso::coad": ("a", "b", "c", "d"),
  "spenso::cos": ("I", "J", "K", "L"),
)

#let _default-representation-indices(name, namespace) = {
  let canonical = if name.contains("::") { name } else { namespace + "::" + name }
  _builtin-representation-index-palettes.at(canonical, default: none)
}

#let _plain-representation(value) = (
  kind: "representation",
  name: value.name,
  namespace: value.namespace,
  dimension: value.dimension,
  self-dual: value.self-dual,
  is-dual: value.is-dual,
  indices: value.at("indices", default: none),
  index-start: value.at("index-start", default: 1),
  index-row: value.at(
    "index-row",
    default: _default-representation-index-row(value.name, value.namespace),
  ),
)

#let _contains-raw-math(value) = {
  if type(value) == array {
    return value.any(_contains-raw-math)
  }
  if type(value) != content { return false }

  let kind = repr(value.func())
  if kind == "raw" { return true }
  if kind == "equation" {
    return _contains-raw-math(value.fields().at("body", default: none))
  }
  if kind == "sequence" {
    return _contains-raw-math(value.fields().at("children", default: ()))
  }
  false
}

#let _portable-index(engine, value) = {
  if type(value) == content {
    if _contains-raw-math(value) {
      panic(
        "manual index math resolved to a non-math Typst value; " +
        "a local binding may shadow the symbol (use a string such as \"mu\", " +
        "or qualify literal math as $std.sym.mu$)",
      )
    }
    return (
      kind: "display-index",
      version: 1,
      ast: _display-ast-bytes(value, engine.grammar),
    )
  }
  value
}

#let _portable-indices(engine, values) = {
  if values == none { none } else { values.map(value => _portable-index(engine, value)) }
}

#let _portable(engine, value) = {
  if type(value) == dictionary and value.at("kind", default: none) == "representation" {
    return (
      kind: "representation",
      name: value.name,
      namespace: value.namespace,
      dimension: _portable(engine, value.dimension),
      self-dual: value.self-dual,
      is-dual: value.is-dual,
      indices: _portable-indices(engine, value.at("indices", default: none)),
      index-start: value.at("index-start", default: 1),
      index-row: value.at(
        "index-row",
        default: _default-representation-index-row(value.name, value.namespace),
      ),
    )
  }
  if type(value) == dictionary and value.at("kind", default: none) == "slot" {
    return (
      kind: "slot",
      representation: _portable(engine, value.representation),
      index: _portable-index(engine, value.index),
      dual: value.dual,
    )
  }
  if type(value) == array {
    return value.map(item => _portable(engine, item))
  }
  if type(value) == dictionary {
    let result = (:)
    for (key, item) in value {
      if type(item) != function {
        result.insert(key, _portable(engine, item))
      }
    }
    return result
  }
  if type(value) == content {
    return _from-math(engine, value)
  }
  value
}

#let _construct(engine, value) = engine.plugin.construct(cbor.encode(_portable(engine, value)))

#let _portable-tensor-argument(engine, value) = {
  if type(value) == content {
    return (
      kind: "display-math",
      version: 1,
      ast: _display-ast-bytes(value, engine.grammar),
    )
  }
  _portable(engine, value)
}

#let _payload(engine, value, label: "expression") = {
  if type(value) == bytes { return value }
  if type(value) == content { return _from-math(engine, value) }
  if type(value) in (int, float, str, dictionary) { return _construct(engine, value) }
  panic(label + " must be an Atom payload, annotated math, or a Tydenso constructor value")
}

#let _slot(representation, index, dual: none) = (
  kind: "slot",
  representation: _plain-representation(representation),
  index: index,
  dual: if dual == none { representation.is-dual } else { dual },
)

#let _as-slot(representation, value, dual: none) = {
  if type(value) == dictionary and value.at("kind", default: none) == "slot" {
    value
  } else {
    _slot(representation, value, dual: dual)
  }
}

#let _atom-envelope(atom, semantic) = (
  protocol: "tymbolica",
  version: 1,
  kind: "atom",
  atom: atom,
  semantic: semantic,
)

#let _annotated-visual(atom, visual, semantic) = {
  _typst-math.attach(visual) + metadata(_atom-envelope(atom, semantic))
}

#let _annotated(engine, atom, semantic) = {
  let tree = cbor(engine.plugin.render_tree(atom))
  let visual = tensor-notation.render(
    tree,
    notation: engine.notation,
  )
  _annotated-visual(atom, visual, _portable(engine, semantic))
}

#let _validate-tags(tags) = {
  if type(tags) != array or not tags.all(tag => type(tag) == str) {
    panic("tags must be an array of strings")
  }
  for (index, tag) in tags.enumerate() {
    let parts = tag.split("::")
    if parts.len() < 2 or parts.any(part => part == "") {
      panic("tags must use canonical namespaced names such as model::positive")
    }
    if tag in tags.slice(0, index) {
      panic("tags must not contain duplicates")
    }
  }
}

#let _call(
  engine,
  name,
  arguments,
  namespace: "spenso",
  symmetric: false,
  antisymmetric: false,
  cycle-symmetric: false,
  linear: false,
  semantic-kind: "tensor",
  tags: (),
) = {
  _validate-tags(tags)
  let constructor = (
    kind: "call",
    name: name,
    namespace: namespace,
    arguments: arguments,
    symmetric: symmetric,
    antisymmetric: antisymmetric,
    cycle-symmetric: cycle-symmetric,
    linear: linear,
    tags: tags,
  )
  let atom = _construct(engine, constructor)
  _annotated(engine, atom, (
    kind: semantic-kind,
    name: name,
    namespace: namespace,
    arguments: arguments,
    symmetric: symmetric,
    antisymmetric: antisymmetric,
    cycle-symmetric: cycle-symmetric,
    linear: linear,
    tags: tags,
  ))
}

#let _symbol(engine, name, namespace: "spenso", tags: ()) = {
  if type(name) != str { panic("symbol name must be a string") }
  _validate-tags(tags)
  let constructor = (
    kind: "symbol",
    name: name,
    namespace: namespace,
    tags: tags,
  )
  let atom = _construct(engine, constructor)
  _annotated(engine, atom, (
    kind: "symbol",
    name: name,
    namespace: namespace,
    tags: tags,
  ))
}

#let _function(engine, name, namespace: "spenso", tags: ()) = {
  if type(name) != str { panic("function name must be a string") }
  _validate-tags(tags)
  (..arguments) => {
    if arguments.named().len() > 0 {
      panic("symbolic function calls accept only positional arguments")
    }
    _call(
      engine,
      name,
      arguments.pos(),
      namespace: namespace,
      semantic-kind: "function-call",
      tags: tags,
    )
  }
}

#let _tensor(
  engine,
  name,
  rank-one: false,
  namespace: "spenso",
  symmetric: false,
  antisymmetric: false,
  cycle-symmetric: false,
  linear: false,
  tags: (),
) = {
  if type(name) != str { panic("tensor name must be a string") }
  _validate-tags(tags)
  (..arguments) => {
    if arguments.named().len() > 0 {
      panic("tensor calls accept only positional arguments")
    }
    let kind = if rank-one { "vector" } else { "tensor" }
    let arguments = arguments.pos().map(
      argument => _portable-tensor-argument(engine, argument),
    )
    let constructor = (
      kind: kind,
      name: name,
      namespace: namespace,
      arguments: arguments,
      symmetric: symmetric,
      antisymmetric: antisymmetric,
      cycle-symmetric: cycle-symmetric,
      linear: linear,
      tags: tags,
    )
    let atom = _construct(engine, constructor)
    _annotated(engine, atom, constructor)
  }
}

#let _dot(engine, left, right) = _call(
  engine,
  "dot",
  (left, right),
  semantic-kind: "dot-product",
)

#let _gamma(engine, lorentz, endpoints) = {
  if endpoints.len() not in (0, 2) {
    panic("gamma needs either one Lorentz argument or two spinor endpoints")
  }
  let arguments = if endpoints.len() == 0 {
    // Spenso chain factors carry literal `in` and `out` placeholder symbols.
    ("in", "out", lorentz)
  } else {
    // The Atom stores spinor endpoints first and the Lorentz argument last.
    (endpoints.at(0), endpoints.at(1), lorentz)
  }
  // Idenso has already registered this head with its tensor tag, linearity,
  // and custom printer. Parsing the existing symbol preserves that definition.
  // `gamma` is also a Symbolica builtin. Qualifying the Idenso head avoids
  // Symbolica's parser resolving this constructor to `symbolica::gamma`.
  _call(engine, "spenso::gamma", arguments, semantic-kind: "gamma")
}

#let _spinor-factor(engine, name, endpoints) = {
  if endpoints.len() not in (0, 2) {
    panic(name + " needs either no arguments or two spinor endpoints")
  }
  let arguments = if endpoints.len() == 0 { ("in", "out") } else { endpoints }
  _call(engine, name, arguments, semantic-kind: name)
}

#let _chain(engine, start, end, factors) = _call(
  engine,
  "chain",
  (start, end) + factors,
  semantic-kind: "chain",
)

#let _cyclic(engine, factors) = _call(
  engine,
  "cyclic",
  factors,
  semantic-kind: "cyclic-projector",
)

#let _trace(engine, representation, factors) = {
  if factors.len() > 1 {
    panic("trace accepts one cyclic(...) payload, not a raw factor list")
  }
  if factors.len() == 1 {
    let projector = cbor(engine.plugin.inspect(_payload(engine, factors.first())))
    if projector.at("name", default: none) != "spenso::cyclic" {
      panic("a non-empty trace must wrap its factors with cyclic(...)")
    }
  }
  _call(
    engine,
    "trace",
    (representation,) + factors,
    semantic-kind: "trace",
  )
}

#let _representation(
  engine,
  name,
  dimension,
  namespace: "spenso",
  self-dual: false,
  is-dual: false,
  indices: none,
  index-start: 1,
  index-row: none,
) = {
  let resolved-indices = if indices == none {
    _default-representation-indices(name, namespace)
  } else {
    indices
  }
  if resolved-indices != none {
    if type(resolved-indices) != array or resolved-indices.len() == 0 {
      panic("indices must be a non-empty array, or none for the representation default")
    }
    if not resolved-indices.all(index => type(index) in (content, str, int)) {
      panic("indices must contain only Typst math content, strings, or integers")
    }
  }
  if type(index-start) != int or index-start < 0 {
    panic("index-start must be a non-negative integer")
  }
  if resolved-indices == none and index-start != 1 {
    panic("index-start requires an indices palette")
  }
  let index-row = if index-row == none {
    _default-representation-index-row(name, namespace)
  } else {
    index-row
  }
  if index-row not in ("top", "bottom") {
    panic("index-row must be \"top\" or \"bottom\"")
  }
  let descriptor = (
    kind: "representation",
    name: name,
    namespace: namespace,
    dimension: dimension,
    self-dual: self-dual,
    is-dual: is-dual,
    // Both orientations share one canonical representation head. The slot's
    // `dual` flag is the only algebraic variance marker.
    indices: resolved-indices,
    index-start: index-start,
    index-row: index-row,
  )
  (
    ..descriptor,
    slot: (index, dual: none) => _slot(descriptor, index, dual: dual),
    dual: () => _representation(
      engine,
      descriptor.name,
      dimension,
      namespace: namespace,
      self-dual: self-dual,
      is-dual: if self-dual { false } else { not is-dual },
      indices: resolved-indices,
      index-start: index-start,
      index-row: index-row,
    ),
    metric: (first, second) => _call(engine, "g", (
      _as-slot(descriptor, first),
      _as-slot(descriptor, second),
    )),
    identity: (first, second) => _call(engine, "id", (
      _as-slot(descriptor, first),
      _as-slot(descriptor, second),
    )),
    flat: (first, second) => _call(engine, "flat", (
      _as-slot(descriptor, first),
      _as-slot(descriptor, second),
    )),
  )
}

#let _print-settings(
  preset: "typst",
  with-dim: none,
  parens: none,
  commas: none,
  index-subscripts: none,
  symbol-scripts: none,
) = {
  assert(preset in ("typst", "compact"), message: "preset must be 'typst' or 'compact'")
  let defaults = if preset == "typst" {
    (with-dim: false, parens: true, commas: false, index-subscripts: true, symbol-scripts: true)
  } else {
    (with-dim: false, parens: true, commas: false, index-subscripts: false, symbol-scripts: false)
  }
  (
    preset: preset,
    with-dim: if with-dim == none { defaults.with-dim } else { with-dim },
    parens: if parens == none { defaults.parens } else { parens },
    commas: if commas == none { defaults.commas } else { commas },
    index-subscripts: if index-subscripts == none { defaults.index-subscripts } else { index-subscripts },
    symbol-scripts: if symbol-scripts == none { defaults.symbol-scripts } else { symbol-scripts },
  )
}

#let _render-semantic(node) = (
  kind: "render-node",
  node-kind: node.at("kind", default: none),
  symbol: node.at("symbol", default: none),
)

#let _resolve-notation(engine, settings: none, notation: none) = {
  if settings != none and notation != none {
    panic("pass either settings or notation to to-typst, not both")
  }
  if notation != none { return notation }
  if settings != none { return tensor-notation.default-notation(settings: settings) }
  engine.notation
}

#let _to-typst(
  engine,
  expression,
  settings: none,
  notation: none,
  block: false,
) = {
  let atom = _payload(engine, expression)
  let tree = cbor(engine.plugin.render_tree(atom))
  tensor-notation.render(
    tree,
    notation: _resolve-notation(
      engine,
      settings: settings,
      notation: notation,
    ),
    annotate: (exact, visual, node) => _annotated-visual(
      exact,
      visual,
      _render-semantic(node),
    ),
    block: block,
  )
}

#let _to-string(engine, expression, settings: _print-settings(preset: "compact")) = str(
  engine.plugin.to_string(cbor.encode((
    expr: _payload(engine, expression),
    settings: settings,
  ))),
)

#let _inspect(engine, expression) = cbor(engine.plugin.inspect(_payload(engine, expression)))

/// Create an independent Tydenso API.
///
/// The returned dictionary contains tensor constructors, document-side math
/// rendering, CBOR inspection, and Idenso transformations. Use a custom `source` only when
/// developing another compatible plugin; ordinary documents can use the
/// imported top-level functions.
///
/// ```example
/// #let tensors = init()
/// #let V = (tensors.mink)(4)
/// #let p = (tensors.vector)("p")
/// #(tensors.to-typst)(p((tensors.slot)(V, 1)))
/// ```
///
/// -> dictionary
#let init(
  /// Default namespace for symbols parsed from Typst math.
  /// -> str
  namespace: "spenso",
  /// WebAssembly plugin path or uncompressed module bytes. `none` selects the
  /// bundled compressed Tydenso engine.
  /// -> str | bytes | none
  source: none,
  /// Parser grammar used by `math` unless it receives an explicit override.
  /// -> dictionary
  grammar: _default-grammar,
  /// Default document-side tensor notation used by constructors and
  /// `to-typst`.
  /// -> dictionary
  notation: tensor-notation.default-notation(),
) = {
  let plugin-module = if source == none { _bundled-plugin() } else { plugin(source) }
  let engine = (
    plugin: plugin-module,
    grammar: grammar,
    namespace: namespace,
    notation: notation,
  )
  (
    plugin: plugin-module,
    math: (equation, grammar: none, namespace: none) => _from-math(
      engine,
      equation,
      grammar: grammar,
      namespace: namespace,
    ),
    atom: value => _payload(engine, value),
    construct: value => _construct(engine, value),
    symbol: (name, namespace: none, tags: ()) => _symbol(
      engine, name, namespace: _namespace(engine, namespace), tags: tags,
    ),
    function: (name, namespace: none, tags: ()) => _function(
      engine, name, namespace: _namespace(engine, namespace), tags: tags,
    ),
    tensor: (name, namespace: none, symmetric: false, antisymmetric: false, cycle-symmetric: false, linear: false, tags: ()) => _tensor(
      engine, name,
      namespace: _namespace(engine, namespace),
      symmetric: symmetric,
      antisymmetric: antisymmetric,
      cycle-symmetric: cycle-symmetric,
      linear: linear,
      tags: tags,
    ),
    vector: (name, namespace: none, linear: false, tags: ()) => _tensor(
      engine, name,
      rank-one: true,
      namespace: _namespace(engine, namespace),
      linear: linear,
      tags: tags,
    ),
    dot: (left, right) => _dot(engine, left, right),
    gamma: (lorentz, ..endpoints) => _gamma(
      engine, lorentz, endpoints.pos(),
    ),
    gamma0: (..endpoints) => _spinor-factor(engine, "gamma0", endpoints.pos()),
    gamma5: (..endpoints) => _spinor-factor(engine, "gamma5", endpoints.pos()),
    projp: (..endpoints) => _spinor-factor(engine, "projp", endpoints.pos()),
    projm: (..endpoints) => _spinor-factor(engine, "projm", endpoints.pos()),
    chain: (start, end, ..factors) => _chain(
      engine, start, end, factors.pos(),
    ),
    cyclic: (..factors) => _cyclic(engine, factors.pos()),
    trace: (representation, ..factors) => _trace(
      engine, representation, factors.pos(),
    ),
    representation: (name, dimension, namespace: none, self-dual: false, is-dual: false, indices: none, index-start: 1, index-row: none) => _representation(
      engine, name, dimension,
      namespace: _namespace(engine, namespace),
      self-dual: self-dual,
      is-dual: is-dual,
      indices: indices,
      index-start: index-start,
      index-row: index-row,
    ),
    mink: dimension => _representation(engine, "mink", dimension, self-dual: true),
    euc: dimension => _representation(engine, "euc", dimension, self-dual: true),
    lor: dimension => _representation(engine, "lor", dimension),
    bis: dimension => _representation(engine, "bis", dimension, self-dual: true),
    spf: dimension => _representation(engine, "spf", dimension),
    cof: dimension => _representation(engine, "cof", dimension),
    coad: dimension => _representation(engine, "coad", dimension, self-dual: true),
    cos: dimension => _representation(engine, "cos", dimension),
    slot: (representation, index, dual: none) => _slot(representation, index, dual: dual),
    metric: (representation, first, second) => _call(engine, "g", (
      _as-slot(representation, first), _as-slot(representation, second),
    )),
    identity-tensor: (representation, first, second) => _call(engine, "id", (
      _as-slot(representation, first), _as-slot(representation, second),
    )),
    flat-tensor: (representation, first, second) => _call(engine, "flat", (
      _as-slot(representation, first), _as-slot(representation, second),
    )),
    dual-representation: representation => (representation.dual)(),
    add: (..terms) => _construct(engine, (kind: "sum", terms: terms.pos())),
    mul: (..factors) => _construct(engine, (kind: "product", factors: factors.pos())),
    neg: expression => _construct(engine, (kind: "negative", expression: expression)),
    sub: (left, right) => _construct(engine, (
      kind: "sum",
      terms: (left, (kind: "negative", expression: right)),
    )),
    div: (numerator, denominator) => _construct(engine, (
      kind: "product",
      factors: (
        numerator,
        (kind: "power", base: denominator, exponent: -1),
      ),
    )),
    pow: (base, exponent) => _construct(engine, (
      kind: "power", base: base, exponent: exponent,
    )),
    print-settings: (preset: "typst", with-dim: none, parens: none, commas: none, index-subscripts: none, symbol-scripts: none) => _print-settings(
      preset: preset,
      with-dim: with-dim,
      parens: parens,
      commas: commas,
      index-subscripts: index-subscripts,
      symbol-scripts: symbol-scripts,
    ),
    to-typst: (expression, settings: none, notation: none, block: false) => _to-typst(
      engine,
      expression,
      settings: settings,
      notation: notation,
      block: block,
    ),
    to-string: (expression, settings: _print-settings(preset: "compact")) => _to-string(engine, expression, settings: settings),
    inspect: expression => _inspect(engine, expression),
    cook-function: expression => plugin-module.cook_function(_payload(engine, expression)),
    cook-indices: expression => plugin-module.cook_indices(_payload(engine, expression)),
    dirac-adjoint: expression => plugin-module.dirac_adjoint(_payload(engine, expression)),
    expand-bis: expression => plugin-module.expand_bis(_payload(engine, expression)),
    expand-color: expression => plugin-module.expand_color(_payload(engine, expression)),
    expand-metrics: expression => plugin-module.expand_metrics(_payload(engine, expression)),
    expand-mink: expression => plugin-module.expand_mink(_payload(engine, expression)),
    expand-mink-bis: expression => plugin-module.expand_mink_bis(_payload(engine, expression)),
    list-dangling: expression => cbor(plugin-module.list_dangling(_payload(engine, expression))),
    simplify-color: expression => plugin-module.simplify_color(_payload(engine, expression)),
    simplify-gamma: expression => plugin-module.simplify_gamma(_payload(engine, expression)),
    simplify-metrics: expression => plugin-module.simplify_metrics(_payload(engine, expression)),
    to-dots: expression => plugin-module.to_dots(_payload(engine, expression)),
    wrap-dummies: (expression, header) => plugin-module.wrap_dummies(
      _payload(engine, expression), _payload(engine, header, label: "header"),
    ),
    wrap-indices: (expression, header) => plugin-module.wrap_indices(
      _payload(engine, expression), _payload(engine, header, label: "header"),
    ),
  )
}

#let _default-engine() = init()

/// Parse Typst math into an exact Symbolica Atom payload.
///
/// Annotated values created by `symbol`, `function`, `tensor`, `vector`, and
/// the named tensor constructors are imported from their exact Atom metadata.
/// Use interpolation for callable bindings, for example `#F(mu, nu)`.
///
/// -> bytes
#let math(
  /// Math content, normally written as `$...$`.
  /// -> content
  equation,
  /// Parser grammar override. `none` uses the engine grammar.
  /// -> dictionary | none
  grammar: none,
  /// Namespace for unannotated parsed symbols.
  /// -> str | none
  namespace: none,
) = (_default-engine().math)(equation, grammar: grammar, namespace: namespace)

/// Convert a supported value or annotated math expression to Atom bytes.
///
/// -> bytes
#let atom(value) = (_default-engine().atom)(value)

/// Construct a named tensor function.
///
/// The result is callable with any mixture of slots, representations, scalar
/// values, display math, and compatible Atom payloads. Raw Typst math in a
/// non-slot argument is retained as one atomic display label; pass `math($...$)`
/// when that argument should instead be parsed as Symbolica algebra. Symmetry
/// flags follow Spenso's `TensorName` model and are registered on the underlying
/// Symbolica symbol.
///
/// ```example
/// #let V = mink(4)
/// #let Z = tensor("Z")
/// #to-typst(Z(slot(V, 1), slot(V, 2)))
/// ```
///
/// -> function
#let tensor(
  /// Tensor name.
  /// -> str
  name,
  /// Symbol namespace.
  /// -> str
  namespace: "spenso",
  /// Make the tensor symmetric under argument permutation.
  /// -> bool
  symmetric: false,
  /// Make the tensor antisymmetric under argument permutation.
  /// -> bool
  antisymmetric: false,
  /// Make the tensor symmetric under cyclic argument permutation.
  /// -> bool
  cycle-symmetric: false,
  /// Mark the tensor function as linear.
  /// -> bool
  linear: false,
  /// Portable semantic labels retained in attached metadata.
  /// -> array
  tags: (),
) = (_default-engine().tensor)(
  name,
  namespace: namespace,
  symmetric: symmetric,
  antisymmetric: antisymmetric,
  cycle-symmetric: cycle-symmetric,
  linear: linear,
  tags: tags,
)

/// Construct a named rank-one tensor function.
///
/// Non-slot arguments are rendered as symbol subscripts in Spenso's Typst
/// mode, while the vector slot is rendered as an abstract index. Raw Typst math
/// is preserved as the notation for one atomic label, so accents, fractions,
/// attachments, and arithmetic can be written directly. Use `math($...$)` when
/// the argument should remain an algebraic expression rather than a label.
///
/// ```example
/// #let V = mink(4)
/// #let p = vector("p")
/// #to-typst(p($arrow(x + y)$, V))
/// ```
///
/// -> function
#let vector(
  /// Vector name.
  /// -> str
  name,
  /// Symbol namespace.
  /// -> str
  namespace: "spenso",
  /// Mark the vector function as linear.
  /// -> bool
  linear: false,
  /// Portable semantic labels retained in attached metadata.
  /// -> array
  tags: (),
) = (_default-engine().vector)(
  name,
  namespace: namespace,
  linear: linear,
  tags: tags,
)

/// Construct Spenso's compact scalar product `dot(left, right)`.
///
/// Rank-one tensor calls with a representation and no explicit slot are the
/// compact Schoonschip form. Thus `dot(p(M), q(M))` emits exactly
/// `dot(p(mink(4)),q(mink(4)))` when `M = mink(4)`.
///
/// -> content
#let dot(left, right) = (_default-engine().dot)(left, right)

/// Construct an Idenso gamma tensor or a gamma chain factor.
///
/// With one argument this emits the actual Spenso factor
/// `gamma(in,out,lorentz)`. With all three arguments it emits
/// `gamma(first,second,lorentz)`, which is the explicit tensor order stored in
/// the Atom even though the Typst API places the Lorentz argument first.
///
/// ```example
/// #let M = mink(4)
/// #let B = bis(4)
/// #let mu = slot(M, 1)
/// #let a = slot(B, 1)
/// #let b = slot(B, 2)
/// #let factor = gamma(mu)
/// #let explicit = gamma(mu, a, b)
/// #to-typst(chain(a, b, factor))
/// ```
///
/// -> content
#let gamma(lorentz, ..endpoints) = (
  _default-engine().gamma
)(lorentz, ..endpoints)

/// Construct the built-in $gamma_0$ spinor factor.
///
/// With no arguments, the hidden endpoints are Spenso's `in` and `out`
/// markers. Pass two explicit spinor slots to display its indices.
///
/// -> content
#let gamma0(..endpoints) = (_default-engine().gamma0)(..endpoints)

/// Construct the built-in $gamma_5$ spinor factor.
///
/// -> content
#let gamma5(..endpoints) = (_default-engine().gamma5)(..endpoints)

/// Construct the positive-chirality spinor projector.
///
/// -> content
#let projp(..endpoints) = (_default-engine().projp)(..endpoints)

/// Construct the negative-chirality spinor projector.
///
/// -> content
#let projm(..endpoints) = (_default-engine().projm)(..endpoints)

/// Construct an open Spenso chain.
///
/// `start` and `end` may be explicit slots or compact rank-one tensors. For
/// example, if `u = vector("u")` and `B = bis(4)`, then `u(B)` is the actual
/// compact endpoint Atom `u(bis(4))`; no `bra` or `ket` head is introduced.
/// Factors such as `gamma(mu)` already contain Spenso's `in` and `out`
/// placeholders.
///
/// -> content
#let chain(start, end, ..factors) = (
  _default-engine().chain
)(start, end, ..factors)

/// Construct Spenso's inert cycle-symmetric projector.
///
/// This emits the actual Atom `cyclic(factors...)`.
///
/// -> content
#let cyclic(..factors) = (_default-engine().cyclic)(..factors)

/// Construct a canonical closed Spenso chain.
///
/// A non-empty call accepts one explicit `cyclic(...)` payload and therefore
/// emits exactly `trace(rep,cyclic(factors...))`. Passing a raw factor list is
/// rejected instead of being silently rewritten. An empty call emits
/// `trace(rep)`.
///
/// ```example
/// #let M = mink(4)
/// #let B = bis(4)
/// #let mu = slot(M, 1)
/// #let nu = slot(M, 2)
/// #to-typst(trace(B, cyclic(gamma(mu), gamma(nu))))
/// ```
///
/// -> content
#let trace(representation, ..factors) = (
  _default-engine().trace
)(representation, ..factors)

/// Construct a symbolic scalar with exact Atom metadata.
///
/// -> content
#let symbol(
  /// Symbol name.
  /// -> str
  name,
  /// Symbol namespace.
  /// -> str
  namespace: "spenso",
  /// Portable semantic labels retained in attached metadata.
  /// -> array
  tags: (),
) = (_default-engine().symbol)(name, namespace: namespace, tags: tags)

/// Construct a callable symbolic function with exact Atom metadata.
///
/// This is distinct from `symbol`: a Typst content value is not callable.
/// Calling the returned function annotates the complete function call.
///
/// -> function
#let function(
  /// Function-head name.
  /// -> str
  name,
  /// Symbol namespace.
  /// -> str
  namespace: "spenso",
  /// Portable semantic labels retained in attached metadata.
  /// -> array
  tags: (),
) = (_default-engine().function)(name, namespace: namespace, tags: tags)

/// Describe a representation and attach its tensor-building methods.
///
/// The returned dictionary is inspectable Typst data. Its `slot`, `metric`,
/// `identity`, `flat`, and `dual` fields are functions.
///
/// -> dictionary
#let representation(
  /// Symbolic representation name.
  /// -> str
  name,
  /// Dimension, which may itself be symbolic.
  /// -> int | str | bytes
  dimension,
  /// Symbol namespace.
  /// -> str
  namespace: "spenso",
  /// Whether the representation is its own dual.
  /// -> bool
  self-dual: false,
  /// Whether this value denotes the dual orientation of a non-self-dual
  /// representation.
  /// -> bool
  is-dual: false,
  /// Fixed cyclic display palette for automatic integer indices. Entries may
  /// be Typst math content, strings, or integers. `none` selects the canonical
  /// built-in palette when one exists, and numeric display otherwise.
  /// -> array | none
  indices: none,
  /// Integer represented by the first palette entry. Each wrap adds a numeric
  /// subscript, so `($mu$, $nu$)` maps 1, 2, 3 to $mu$, $nu$, $mu_1$.
  /// -> int
  index-start: 1,
  /// Preferred Typst script row for the representation's base orientation.
  /// `none` selects bottom for the built-in `spenso::bis` representation and
  /// top for every other representation. A dualizable representation's dual
  /// orientation uses the opposite row.
  /// -> str | none
  index-row: none,
) = (_default-engine().representation)(
  name, dimension,
  namespace: namespace,
  self-dual: self-dual,
  is-dual: is-dual,
  indices: indices,
  index-start: index-start,
  index-row: index-row,
)

/// Construct a Minkowski representation.
///
/// Integer slots cycle through `mu`, `nu`, `rho`, and `sigma`. This is a
/// self-dual inline-metric representation whose indices print on the top row.
///
/// -> dictionary
#let mink(dimension) = (_default-engine().mink)(dimension)

/// Construct a Euclidean representation.
///
/// Integer slots cycle through `i`, `j`, `k`, and `l`. This representation is
/// self-dual and prints its indices on the top row.
///
/// -> dictionary
#let euc(dimension) = (_default-engine().euc)(dimension)

/// Construct a Lorentz representation.
///
/// Integer slots cycle through `mu`, `nu`, `rho`, and `sigma`. Base indices
/// print on the top row and dual indices on the bottom row.
///
/// -> dictionary
#let lor(dimension) = (_default-engine().lor)(dimension)

/// Construct a bispinor representation.
///
/// Integer slots cycle through `a`, `b`, `c`, and `d`. The representation is
/// self-dual and conventionally prints all of its indices on the bottom row.
///
/// -> dictionary
#let bis(dimension) = (_default-engine().bis)(dimension)

/// Construct a spin-fundamental representation.
///
/// Integer slots cycle through `alpha`, `beta`, `gamma`, and `delta`. Base
/// indices print on the top row and dual indices on the bottom row.
///
/// -> dictionary
#let spf(dimension) = (_default-engine().spf)(dimension)

/// Construct a color-fundamental representation.
///
/// Integer slots cycle through `i`, `j`, `k`, and `l`. Base indices print on
/// the top row and dual indices on the bottom row.
///
/// -> dictionary
#let cof(dimension) = (_default-engine().cof)(dimension)

/// Construct a color-adjoint representation.
///
/// Integer slots cycle through `a`, `b`, `c`, and `d`. This representation is
/// self-dual and prints its indices on the top row.
///
/// -> dictionary
#let coad(dimension) = (_default-engine().coad)(dimension)

/// Construct a color-sextet representation.
///
/// Integer slots cycle through `I`, `J`, `K`, and `L`. Base indices print on
/// the top row and dual indices on the bottom row.
///
/// -> dictionary
#let cos(dimension) = (_default-engine().cos)(dimension)

/// Construct an indexed slot from a representation.
///
/// The result is an ordinary dictionary, so its representation, index, and
/// variance remain visible to Typst before an expression is constructed.
///
/// -> dictionary
#let slot(
  /// Representation dictionary.
  /// -> dictionary
  representation,
  /// Abstract index label. Typst math content such as `$mu_1$` is retained as
  /// safe display metadata on a symbolic index. A string instead names an
  /// ordinary Symbolica identifier and uses Symbolica's native quoted Typst
  /// output. If a local variable shadows a math name, qualify the literal
  /// symbol, for example `$std.sym.mu$`.
  /// -> int | str | bytes | content
  index,
  /// Wrap the slot with Spenso's dual-index marker. `none` inherits the
  /// representation's `is-dual` field.
  /// -> bool | none
  dual: none,
) = (_default-engine().slot)(representation, index, dual: dual)

/// Construct the metric tensor for a representation.
///
/// Index labels are converted to slots automatically; existing slots pass
/// through unchanged.
///
/// -> content
#let metric(representation, first, second) = (
  _default-engine().metric
)(representation, first, second)

/// Construct the identity tensor for a representation.
///
/// -> content
#let identity-tensor(representation, first, second) = (
  _default-engine().identity-tensor
)(representation, first, second)

/// Construct Spenso's musical isomorphism tensor for a representation.
///
/// -> content
#let flat-tensor(representation, first, second) = (
  _default-engine().flat-tensor
)(representation, first, second)

/// Return the representation paired with this representation.
///
/// -> dictionary
#let dual-representation(representation) = (
  _default-engine().dual-representation
)(representation)

/// Add symbolic expressions or constructor values.
///
/// -> bytes
#let add(..terms) = (_default-engine().add)(..terms)

/// Multiply symbolic expressions or constructor values.
///
/// -> bytes
#let mul(..factors) = (_default-engine().mul)(..factors)

/// Negate a symbolic expression or constructor value.
///
/// -> bytes
#let neg(expression) = (_default-engine().neg)(expression)

/// Subtract two symbolic expressions or constructor values.
///
/// -> bytes
#let sub(left, right) = (_default-engine().sub)(left, right)

/// Divide two symbolic expressions or constructor values.
///
/// -> bytes
#let div(numerator, denominator) = (_default-engine().div)(numerator, denominator)

/// Raise a symbolic expression or constructor value to a power.
///
/// -> bytes
#let pow(base, exponent) = (_default-engine().pow)(base, exponent)

/// Configure Spenso's string printer.
///
/// These settings are primarily for `to-string`. Typst math layout is more
/// directly configured through `notation`.
///
/// -> dictionary
#let print-settings(
  /// Base preset: `"typst"` or `"compact"`.
  /// -> str
  preset: "typst",
  /// Include representation dimensions in printed slots.
  /// -> bool | none
  with-dim: none,
  /// Use parentheses around tensor arguments.
  /// -> bool | none
  parens: none,
  /// Separate tensor arguments with commas.
  /// -> bool | none
  commas: none,
  /// Render indices as subscripts where supported.
  /// -> bool | none
  index-subscripts: none,
  /// Render non-slot tensor arguments as symbol scripts where supported.
  /// -> bool | none
  symbol-scripts: none,
) = (_default-engine().print-settings)(
  preset: preset,
  with-dim: with-dim,
  parens: parens,
  commas: commas,
  index-subscripts: index-subscripts,
  symbol-scripts: symbol-scripts,
)

/// Build Tydenso's document-side notation.
///
/// Exact `heads` and `calls` customize one namespaced Symbolica head or full
/// call. `tags` and `classes` customize a family. Full-call renderers receive
/// the complete arguments through their context; exact Atom metadata is added
/// after they return.
///
/// -> dictionary
#let notation(
  /// Reusable notation settings. Directly named options override this layer.
  /// -> dictionary
  settings: (:),
  /// Tensor layout: `"ports"` keeps the bra/ket port notation,
  /// `"schoonschip"` places bold compact vectors in aligned script columns,
  /// and `"call"` groups top and bottom rows around a semicolon.
  /// -> str | none
  tensor-layout: none,
  /// Show each representation and dimension on its indices.
  /// -> bool | none
  with-dim: none,
  /// Enclose chain bodies in brackets.
  /// -> bool | none
  parens: none,
  /// Separate ordinary tensor arguments with commas.
  /// -> bool | none
  commas: none,
  /// Place non-index arguments in tensor scripts when possible.
  /// -> bool | none
  symbol-scripts: none,
  /// Space between adjacent index columns.
  /// -> length | none
  index-gap: none,
  /// Space between factors in a chain or trace.
  /// -> length | none
  factor-gap: none,
  /// Exact namespaced head renderers or content.
  /// -> dictionary
  heads: (:),
  /// Exact namespaced full-call renderers.
  /// -> dictionary
  calls: (:),
  /// Symbolica-tag renderer fallbacks.
  /// -> dictionary
  tags: (:),
  /// Attribute/class renderer fallbacks.
  /// -> dictionary
  classes: (:),
) = tensor-notation.notation(
  settings: settings,
  tensor-layout: tensor-layout,
  with-dim: with-dim,
  parens: parens,
  commas: commas,
  symbol-scripts: symbol-scripts,
  index-gap: index-gap,
  factor-gap: factor-gap,
  heads: heads,
  calls: calls,
  tags: tags,
  classes: classes,
)

/// Merge notation layers from least to most specific.
///
/// -> dictionary
#let merge-notation(..layers) = tensor-notation.merge-notation(..layers)

/// Print and evaluate an expression as Typst math content.
///
/// The plugin exports a generic Atom tree. Tydenso's pure Typst notation lays
/// out tensors, indices, gamma matrices, dots, chains, and traces. Generic
/// algebra stays structurally visible; a customized full call receives one
/// exact metadata annotation after its renderer returns.
///
/// -> content
#let to-typst(
  /// Atom payload or constructor value.
  expression,
  /// Compatibility shorthand for a default notation with these settings.
  /// Cannot be combined with `notation`.
  /// -> dictionary | none
  settings: none,
  /// Per-call notation. `none` uses the engine default.
  /// -> dictionary | none
  notation: none,
  /// Display the result as a block equation.
  /// -> bool
  block: false,
) = (_default-engine().to-typst)(
  expression,
  settings: settings,
  notation: notation,
  block: block,
)

/// Print an expression in Spenso's compact Symbolica notation.
///
/// -> str
#let to-string(
  /// Atom payload or constructor value.
  expression,
  /// Settings created by `print-settings`.
  /// -> dictionary
  settings: _print-settings(preset: "compact"),
) = (_default-engine().to-string)(expression, settings: settings)

/// Decode an Atom payload into a recursive CBOR expression tree.
///
/// The result uses `kind` values such as `symbol`, `number`, `function`,
/// `sum`, `product`, and `power`. Function nodes include their name, arguments,
/// and symmetry flags. This view is ordinary Typst data; transforming the
/// expression still uses its lossless Atom payload.
///
/// -> dictionary
#let inspect(expression) = (_default-engine().inspect)(expression)

/// Cook a tensor function into Idenso's canonical structure.
///
/// -> bytes
#let cook-function(expression) = (_default-engine().cook-function)(expression)

/// Cook indices into Idenso's canonical structure.
///
/// -> bytes
#let cook-indices(expression) = (_default-engine().cook-indices)(expression)

/// Take the Dirac adjoint of a tensor expression.
///
/// -> bytes
#let dirac-adjoint(expression) = (_default-engine().dirac-adjoint)(expression)

/// Selectively expand bispinor structures.
///
/// -> bytes
#let expand-bis(expression) = (_default-engine().expand-bis)(expression)

/// Selectively expand color structures.
///
/// -> bytes
#let expand-color(expression) = (_default-engine().expand-color)(expression)

/// Expand all supported metric structures.
///
/// -> bytes
#let expand-metrics(expression) = (_default-engine().expand-metrics)(expression)

/// Selectively expand Minkowski structures.
///
/// -> bytes
#let expand-mink(expression) = (_default-engine().expand-mink)(expression)

/// Selectively expand combined Minkowski and bispinor structures.
///
/// -> bytes
#let expand-mink-bis(expression) = (_default-engine().expand-mink-bis)(expression)

/// Return the dangling indices as compatible Atom payloads.
///
/// -> array
#let list-dangling(expression) = (_default-engine().list-dangling)(expression)

/// Simplify color tensors.
///
/// -> bytes
#let simplify-color(expression) = (_default-engine().simplify-color)(expression)

/// Simplify gamma-matrix expressions.
///
/// -> bytes
#let simplify-gamma(expression) = (_default-engine().simplify-gamma)(expression)

/// Contract and simplify metric tensors.
///
/// -> bytes
#let simplify-metrics(expression) = (_default-engine().simplify-metrics)(expression)

/// Rewrite supported contractions as dot products.
///
/// -> bytes
#let to-dots(expression) = (_default-engine().to-dots)(expression)

/// Wrap dummy indices under the given header symbol.
///
/// -> bytes
#let wrap-dummies(expression, header) = (_default-engine().wrap-dummies)(expression, header)

/// Wrap all indices under the given header symbol.
///
/// -> bytes
#let wrap-indices(expression, header) = (_default-engine().wrap-indices)(expression, header)
