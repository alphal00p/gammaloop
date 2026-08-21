use crate::{
    shadowing::{
        CYCLIC, SYM,
        symbolica_utils::{SpensoPrintBackend, SpensoPrintSettings},
    },
    structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{IndexDisplay, IndexRow, RepresentationClass, RepresentationMetadata},
    },
};
use symbolica::{
    atom::{
        Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, NamespacedSymbol, Symbol,
        SymbolAttribute, SymbolBuilder,
    },
    coefficient::CoefficientView,
    printer::{PrintOptions, PrintState},
    symbol, tag,
};
pub struct SpensoTags {
    pub broadcast: String,
    /// Marks rank-one tensor symbols whose final argument is the tensor slot.
    ///
    /// Symbols carrying this tag must not use representation slots in earlier
    /// arguments; rank-one shorthand rewrites assume that contract.
    pub rank1: String,
    pub rank1_: Symbol,
    /// Reserved chain wiring markers registered in Spenso's Symbolica namespace.
    pub chain_in: Symbol,
    pub chain_out: Symbol,
    pub chain: Symbol,
    pub trace: Symbol,
    pub upper: String,
    pub lower: String,
    pub bracket: Symbol,
    pub pure_scalar: Symbol,
    pub scalar: Symbol,
    pub tensor: String,
    pub tensor_: Symbol,
    /// Internal wrapper used to restore tensor printing after importing an
    /// Atom whose dynamically registered Rust print callback was not exported.
    pub tensor_display: Symbol,
    pub index: String,
    pub representation: String,
    pub i_: Symbol,
    pub dot: Symbol,
    pub rep_: Symbol,
    pub self_dual: String,
    pub self_dual_: Symbol,
    pub dualizable: String,
    pub dualizable_: Symbol,
}

crate::symbolica_init_lazy_static! {
    pub static SPENSO_TAG, SPENSO_TAG_INNER: SpensoTags = SpensoTags::new;
}

pub fn scalar_store_alias(index: usize) -> Atom {
    FunctionBuilder::new(SPENSO_TAG.scalar)
        .add_arg(Atom::num(
            i64::try_from(index).expect("scalar alias index must fit in i64"),
        ))
        .finish()
}

pub fn scalar_store_alias_index(value: AtomView<'_>) -> Option<usize> {
    let AtomView::Fun(fun) = value else {
        return None;
    };
    if fun.get_symbol() != SPENSO_TAG.scalar || fun.get_nargs() != 1 {
        return None;
    }

    let AtomView::Num(index) = fun.iter().next()? else {
        return None;
    };
    match index.get_coeff_view() {
        CoefficientView::Natural(index, 1, 0, 1) => usize::try_from(index).ok(),
        _ => None,
    }
}

fn typst_builtin_name(name: &str) -> bool {
    matches!(
        name,
        "alpha"
            | "beta"
            | "gamma"
            | "delta"
            | "epsilon"
            | "zeta"
            | "eta"
            | "theta"
            | "iota"
            | "kappa"
            | "lambda"
            | "mu"
            | "nu"
            | "xi"
            | "omicron"
            | "pi"
            | "rho"
            | "sigma"
            | "tau"
            | "upsilon"
            | "phi"
            | "chi"
            | "psi"
            | "omega"
            | "Alpha"
            | "Beta"
            | "Gamma"
            | "Delta"
            | "Epsilon"
            | "Zeta"
            | "Eta"
            | "Theta"
            | "Iota"
            | "Kappa"
            | "Lambda"
            | "Mu"
            | "Nu"
            | "Xi"
            | "Omicron"
            | "Pi"
            | "Rho"
            | "Sigma"
            | "Tau"
            | "Upsilon"
            | "Phi"
            | "Chi"
            | "Psi"
            | "Omega"
    )
}

fn escape_typst_string(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

fn typst_tensor_head(symbol: Symbol) -> String {
    match symbol.get_name() {
        "spenso::projp" => return "ℙ_p".to_owned(),
        "spenso::projm" => return "ℙ_m".to_owned(),
        "spenso::gamma5" => return "gamma_5".to_owned(),
        "spenso::gamma0" => return "gamma_0".to_owned(),
        _ => {}
    }

    let name = symbol.get_stripped_name();
    if (name.chars().count() == 1
        && name
            .chars()
            .next()
            .is_some_and(|character| character.is_alphanumeric()))
        || typst_builtin_name(name)
    {
        name.to_owned()
    } else {
        format!(r#"italic("{}")"#, escape_typst_string(name))
    }
}

fn natural_index(index: AtomView<'_>) -> Option<usize> {
    let AtomView::Num(index) = index else {
        return None;
    };
    match index.get_coeff_view() {
        CoefficientView::Natural(index, 1, 0, 1) => usize::try_from(index).ok(),
        _ => None,
    }
}

fn typst_index_source(
    representation: Symbol,
    index: AtomView<'_>,
    options: &PrintOptions,
) -> Option<String> {
    if let Some(index) = natural_index(index)
        && let Some(display) = RepresentationMetadata::from_symbol(representation)
            .and_then(|metadata| metadata.index_palette.resolve(index))
    {
        return Some(display.to_typst_source());
    }

    if let AtomView::Var(variable) = index {
        let symbol = variable.get_symbol();
        if symbol.has_tag(&SPENSO_TAG.index)
            && let Some(display) = IndexDisplay::from_symbol(symbol)
        {
            return Some(display.to_typst_source());
        }
        let name = symbol.get_stripped_name();
        if typst_builtin_name(name) {
            return Some(name.to_owned());
        }
    }

    let mut output = String::new();
    index.format(&mut output, options, PrintState::new()).ok()?;
    Some(output)
}

fn typst_source(atom: AtomView<'_>, options: &PrintOptions) -> Option<String> {
    let mut output = String::new();
    atom.format(&mut output, options, PrintState::new()).ok()?;
    Some(output)
}

fn unwrap_tensor_display(mut atom: AtomView<'_>) -> AtomView<'_> {
    loop {
        let AtomView::Fun(wrapper) = atom else {
            return atom;
        };
        if wrapper.get_symbol() != SPENSO_TAG.tensor_display || wrapper.get_nargs() != 1 {
            return atom;
        }
        atom = wrapper.iter().next().expect("one-argument tensor wrapper");
    }
}

fn is_chain_marker(atom: AtomView<'_>, marker: Symbol) -> bool {
    matches!(unwrap_tensor_display(atom), AtomView::Var(variable) if variable.get_symbol() == marker)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ChainOrientation {
    Forward,
    Transposed,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ChainMarkers {
    input: usize,
    output: usize,
    orientation: ChainOrientation,
}

/// Locate the two direct chain placeholders without imposing an argument order.
///
/// Spenso's chain materializer permits ordinary tensor arguments before,
/// between, and after the placeholders. Ambiguous factors with a missing or
/// repeated placeholder are deliberately left alone.
fn chain_markers(arguments: &[AtomView<'_>]) -> Option<ChainMarkers> {
    let mut input = None;
    let mut output = None;

    for (position, argument) in arguments.iter().enumerate() {
        if is_chain_marker(*argument, SPENSO_TAG.chain_in) {
            if input.replace(position).is_some() {
                return None;
            }
        } else if is_chain_marker(*argument, SPENSO_TAG.chain_out)
            && output.replace(position).is_some()
        {
            return None;
        }
    }

    let input = input?;
    let output = output?;
    Some(ChainMarkers {
        input,
        output,
        orientation: if input < output {
            ChainOrientation::Forward
        } else {
            ChainOrientation::Transposed
        },
    })
}

fn typst_attachment(base: String, columns: Vec<(String, IndexRow)>) -> String {
    if columns.is_empty() {
        return base;
    }

    let mut top = Vec::with_capacity(columns.len());
    let mut bottom = Vec::with_capacity(columns.len());
    for (source, row) in columns {
        let hidden = format!("std.hide({source})");
        match row {
            IndexRow::Top => {
                top.push(source);
                bottom.push(hidden);
            }
            IndexRow::Bottom => {
                top.push(hidden);
                bottom.push(source);
            }
        }
    }

    format!(
        "attach(#(${base}$,std.hide($zws$)).join(),t:{},b:{})",
        top.join(" "),
        bottom.join(" ")
    )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CompactPolarity {
    Bra,
    Ket,
}

#[derive(Clone, Copy)]
struct RepresentationValue<'a> {
    symbol: Symbol,
    dimension: AtomView<'a>,
    class: RepresentationClass,
    row: IndexRow,
    polarity: CompactPolarity,
}

fn representation_value(argument: AtomView<'_>) -> Option<RepresentationValue<'_>> {
    let argument = unwrap_tensor_display(argument);
    let (representation, dual) = if let AtomView::Fun(dual) = argument
        && dual.get_symbol() == AIND_SYMBOLS.dind
        && dual.get_nargs() == 1
    {
        (unwrap_tensor_display(dual.iter().next()?), true)
    } else {
        (argument, false)
    };

    let AtomView::Fun(representation) = representation else {
        return None;
    };
    if representation.get_nargs() != 1
        || !representation
            .get_symbol()
            .has_tag(&SPENSO_TAG.representation)
    {
        return None;
    }

    let metadata = RepresentationMetadata::from_symbol(representation.get_symbol())?;
    let row = if dual && metadata.class == RepresentationClass::Dualizable {
        metadata.index_row.opposite()
    } else {
        metadata.index_row
    };
    let polarity = match metadata.class {
        RepresentationClass::Dualizable if !dual => CompactPolarity::Ket,
        RepresentationClass::SelfDual
        | RepresentationClass::InlineMetric
        | RepresentationClass::Dualizable => CompactPolarity::Bra,
    };

    Some(RepresentationValue {
        symbol: representation.get_symbol(),
        dimension: representation.iter().next()?,
        class: metadata.class,
        row,
        polarity,
    })
}

fn qualified_typst_port(
    representation: RepresentationValue<'_>,
    settings: SpensoPrintSettings,
    options: &PrintOptions,
) -> Option<String> {
    let port = "○";
    if !settings.with_dim {
        return Some(port.to_owned());
    }

    let metadata = RepresentationMetadata::from_symbol(representation.symbol)?;
    let dimension = typst_source(representation.dimension, options)?;
    Some(format!(
        "attach({port},t:attach({},b:{dimension}))",
        metadata.label.to_typst_source()
    ))
}

struct CompactVector<'a> {
    label: String,
    representation: RepresentationValue<'a>,
}

fn compact_vector<'a>(
    atom: AtomView<'a>,
    settings: SpensoPrintSettings,
    options: &PrintOptions,
) -> Option<CompactVector<'a>> {
    let AtomView::Fun(function) = unwrap_tensor_display(atom) else {
        return None;
    };
    let symbol = function.get_symbol();
    if !symbol.has_tag(&SPENSO_TAG.tensor) || !symbol.has_tag(&SPENSO_TAG.rank1) {
        return None;
    }

    let arguments = function.iter().collect::<Vec<_>>();
    let mut representation = None;
    let mut label_arguments = Vec::with_capacity(arguments.len().saturating_sub(1));
    for argument in arguments {
        if tensor_slot(argument).is_some() {
            return None;
        }
        if let Some(candidate) = representation_value(argument) {
            if representation.replace(candidate).is_some() {
                return None;
            }
        } else {
            label_arguments.push(argument);
        }
    }
    let representation = representation?;

    let mut label = typst_tensor_head(symbol);
    if settings.symbol_scripts && !label_arguments.is_empty() {
        let columns = label_arguments
            .iter()
            .map(|argument| Some((typst_source(*argument, options)?, IndexRow::Bottom)))
            .collect::<Option<Vec<_>>>()?;
        label = typst_attachment(label, columns);
    } else if !label_arguments.is_empty() {
        let separator = if settings.commas { "," } else { " " };
        label.push('(');
        label.push_str(
            &label_arguments
                .iter()
                .map(|argument| typst_source(*argument, options))
                .collect::<Option<Vec<_>>>()?
                .join(separator),
        );
        label.push(')');
    }

    Some(CompactVector {
        label,
        representation,
    })
}

#[derive(Clone, Copy)]
struct TensorSlot<'a> {
    representation: Symbol,
    dimension: AtomView<'a>,
    index: AtomView<'a>,
    row: IndexRow,
}

fn tensor_slot(index: AtomView<'_>) -> Option<TensorSlot<'_>> {
    let index = unwrap_tensor_display(index);
    let (slot, dual) = if let AtomView::Fun(dual) = index
        && dual.get_symbol() == AIND_SYMBOLS.dind
        && dual.get_nargs() == 1
    {
        (dual.iter().next()?, true)
    } else {
        (index, false)
    };

    let AtomView::Fun(representation) = slot else {
        return None;
    };
    if !representation
        .get_symbol()
        .has_tag(&SPENSO_TAG.representation)
        || representation.get_nargs() != 2
    {
        return None;
    }

    let mut arguments = representation.iter();
    let metadata = RepresentationMetadata::from_symbol(representation.get_symbol())?;
    let row = if dual && metadata.class == RepresentationClass::Dualizable {
        metadata.index_row.opposite()
    } else {
        metadata.index_row
    };

    Some(TensorSlot {
        representation: representation.get_symbol(),
        dimension: arguments.next()?,
        index: arguments.next()?,
        row,
    })
}

fn qualified_typst_index(
    slot: TensorSlot<'_>,
    index_source: String,
    options: &PrintOptions,
) -> Option<String> {
    let representation_source = RepresentationMetadata::from_symbol(slot.representation)
        .map(|metadata| metadata.label.to_typst_source())
        .unwrap_or_else(|| typst_tensor_head(slot.representation));

    let mut dimension_source = String::new();
    slot.dimension
        .format(&mut dimension_source, options, PrintState::new())
        .ok()?;

    Some(format!(
        "attach({index_source},t:attach({representation_source},b:{dimension_source}))"
    ))
}

/// Print a tagged tensor using native Typst attachments in Spenso's Typst mode.
///
/// Every script occupies the same horizontal column in the top and bottom
/// rows. The opposite row receives a hidden copy, following Physica's tensor
/// layout technique. Each representation owns its preferred row; only the
/// dual orientation of a dualizable representation flips that row.
pub fn tensor_print(
    atom: AtomView<'_>,
    options: &PrintOptions,
    _state: &PrintState,
) -> Option<String> {
    if !options.mode.is_typst() {
        return None;
    }

    let resolved = SpensoPrintSettings::resolve(options)?;
    if !matches!(resolved.backend, SpensoPrintBackend::Typst) {
        return None;
    }
    let settings = resolved.presentation;

    let AtomView::Fun(function) = atom else {
        return None;
    };
    if !function.get_symbol().has_tag(&SPENSO_TAG.tensor) {
        return None;
    }

    let function_symbol = function.get_symbol();
    let function_name = function_symbol.get_name();
    if function_name == "spenso::gamma" {
        return gamma_typst_print(function.as_view(), options, settings);
    }

    let mut columns = Vec::new();
    let mut ordinary_arguments = Vec::new();
    let mut bras = Vec::new();
    let mut kets = Vec::new();
    let arguments = function.iter().collect::<Vec<_>>();
    let markers = chain_markers(&arguments);

    for (position, argument) in arguments.into_iter().enumerate() {
        if markers.is_some_and(|markers| position == markers.input || position == markers.output) {
            continue;
        }
        if let Some(slot) = tensor_slot(argument) {
            let source = typst_index_source(slot.representation, slot.index, options)?;
            let source = if settings.with_dim {
                qualified_typst_index(slot, source, options)?
            } else {
                source
            };
            columns.push((source, slot.row));
        } else if let Some(compact) = compact_vector(argument, settings, options) {
            let port = qualified_typst_port(compact.representation, settings, options)?;
            columns.push((port, compact.representation.row));
            match compact.representation.polarity {
                CompactPolarity::Bra => bras.push(compact.label),
                CompactPolarity::Ket => kets.push(compact.label),
            }
        } else if settings.symbol_scripts {
            columns.push((typst_source(argument, options)?, IndexRow::Bottom));
        } else {
            ordinary_arguments.push(typst_source(argument, options)?);
        }
    }

    let mut base = typst_tensor_head(function.get_symbol());
    if !ordinary_arguments.is_empty() {
        let separator = if settings.commas { "," } else { " " };
        base.push('(');
        base.push_str(&ordinary_arguments.join(separator));
        base.push(')');
    }

    base = typst_attachment(base, columns);

    if matches!(
        markers,
        Some(ChainMarkers {
            orientation: ChainOrientation::Transposed,
            ..
        })
    ) {
        base = format!("attach({base},t:upright(\"T\"))");
    }

    if !bras.is_empty() {
        base = format!(r#"upright("⟨") {} upright("|") {base}"#, bras.join(","));
    }
    if !kets.is_empty() {
        base = format!(r#"{base} upright("|") {} upright("⟩")"#, kets.join(","));
    }
    if markers.is_some() && (!bras.is_empty() || !kets.is_empty()) {
        base = format!("lr(({base}))");
    }
    Some(base)
}

fn gamma_typst_print(
    atom: AtomView<'_>,
    options: &PrintOptions,
    settings: SpensoPrintSettings,
) -> Option<String> {
    let AtomView::Fun(function) = unwrap_tensor_display(atom) else {
        return None;
    };
    if function.get_nargs() != 3 {
        return None;
    }
    let arguments = function.iter().collect::<Vec<_>>();
    let [first, second, lorentz] = arguments.as_slice() else {
        return None;
    };

    let forward = is_chain_marker(*first, SPENSO_TAG.chain_in)
        && is_chain_marker(*second, SPENSO_TAG.chain_out);
    let transposed = is_chain_marker(*first, SPENSO_TAG.chain_out)
        && is_chain_marker(*second, SPENSO_TAG.chain_in);

    if (forward || transposed)
        && let Some(compact) = compact_vector(*lorentz, settings, options)
    {
        let mut output = format!("cancel({})", compact.label);
        if transposed {
            output = format!("attach({output},t:upright(\"T\"))");
        }
        return Some(output);
    }

    let selected = if forward || transposed {
        vec![*lorentz]
    } else {
        arguments
    };
    let mut columns = Vec::new();
    for argument in selected {
        if let Some(slot) = tensor_slot(argument) {
            let source = typst_index_source(slot.representation, slot.index, options)?;
            let source = if settings.with_dim {
                qualified_typst_index(slot, source, options)?
            } else {
                source
            };
            columns.push((source, slot.row));
        } else {
            columns.push((typst_source(argument, options)?, IndexRow::Bottom));
        }
    }
    let mut output = typst_attachment("gamma".to_owned(), columns);
    if transposed {
        output = format!("attach({output},t:upright(\"T\"))");
    }
    Some(output)
}

/// Route every tagged tensor through the portable Typst renderer.
///
/// Symbolica exports tensor tags but intentionally does not export custom Rust
/// functions. Locally registered Idenso tensors can also carry callbacks whose
/// legacy Typst branches do not understand representation-owned rows. The
/// temporary wrapper gives both cases one consistent renderer without changing
/// the algebraic Atom.
pub fn prepare_tensor_print(atom: &Atom) -> Atom {
    atom.replace_map_bottom_up(|view, _, output| {
        let AtomView::Fun(function) = view else {
            return;
        };
        let symbol = function.get_symbol();
        if symbol == SPENSO_TAG.tensor_display && function.get_nargs() == 1 {
            let argument = function.iter().next().expect("one-argument tensor wrapper");
            **output = FunctionBuilder::new(SPENSO_TAG.tensor_display)
                .add_arg(unwrap_tensor_display(argument))
                .finish();
            return;
        }
        if symbol.has_tag(&SPENSO_TAG.tensor) {
            **output = FunctionBuilder::new(SPENSO_TAG.tensor_display)
                .add_arg(view)
                .finish();
        }
    })
}

/// Register a generic tensor or vector head, returning an equivalent existing
/// symbol when it has already been declared.
///
/// Reusing an existing symbol is essential for dynamic frontends: Symbolica
/// deliberately rejects attempts to register the same Rust callback twice.
pub fn register_tensor_symbol(
    name: NamespacedSymbol,
    attributes: Vec<SymbolAttribute>,
    rank_one: bool,
) -> Result<Symbol, String> {
    if let Some(existing) = Symbol::get_symbol(name.clone()) {
        if !existing.has_tag(&SPENSO_TAG.tensor)
            || existing.has_tag(&SPENSO_TAG.rank1) != rank_one
            || existing.get_attributes() != attributes
        {
            return Err(format!(
                "symbol {} already exists with a different tensor declaration",
                existing.get_name()
            ));
        }
        return Ok(existing);
    }

    let tags = if rank_one {
        vec![SPENSO_TAG.tensor.clone(), SPENSO_TAG.rank1.clone()]
    } else {
        vec![SPENSO_TAG.tensor.clone()]
    };
    SymbolBuilder::new(name)
        .with_attributes(attributes)
        .with_tags(tags)
        .build()
        .map_err(|error| error.to_string())
}

/// Builds Symbolica atoms from a symbol and an optional argument list.
///
/// Spenso wildcard heads use the convention that a bare head with no arguments
/// is a variable atom, while a head with arguments is a function atom. This
/// helper keeps that convention in one place.
pub trait SymbolAtomExt {
    fn atom_with_args<'a, A>(self, args: impl IntoIterator<Item = A>) -> Atom
    where
        A: Into<AtomOrView<'a>>;
}

impl SymbolAtomExt for Symbol {
    fn atom_with_args<'a, A>(self, args: impl IntoIterator<Item = A>) -> Atom
    where
        A: Into<AtomOrView<'a>>,
    {
        let mut function = FunctionBuilder::new(self);
        let mut has_args = false;
        for arg in args {
            has_args = true;
            function = function.add_arg(arg);
        }

        if has_args {
            function.finish()
        } else {
            Atom::var(self)
        }
    }
}

macro_rules! define_numbered_tag_family_methods {
    ($(
        $(#[$meta:meta])*
        $vis:vis fn $method:ident => $base_field:ident, $prefix:literal, $symbol_method:ident;
    )*) => {
        $(
            $(#[$meta])*
            $vis fn $method<'a, const N: usize, A: Into<AtomOrView<'a>>>(
                &self,
                args: impl IntoIterator<Item = A>,
            ) -> Atom {
                let symbol = if N == 0 {
                    self.$base_field
                } else {
                    self.$symbol_method(&format!("{}{}_", $prefix, N))
                };
                symbol.atom_with_args(args)
            }
        )*
    };
}

macro_rules! define_numbered_tag_macros {
    ($d:tt; $($macro_name:ident => $method:ident;)*) => {
        $(
            #[macro_export]
            macro_rules! $macro_name {
                ($d n:literal; $d($d arg:expr),* $d(,)?) => {
                    $crate::network::tags::SPENSO_TAG.$method::<$d n, _>(
                        vec![$d($crate::shadowing::IntoAtom::into_atom($d arg)),*],
                    )
                };
                ($d n:literal $d(;)?) => {
                    $crate::network::tags::SPENSO_TAG.$method::<$d n, symbolica::atom::Atom>(
                        std::iter::empty::<symbolica::atom::Atom>(),
                    )
                };
            }
        )*
    };
}

define_numbered_tag_macros!($;
    rank1_ => rank1_;
    tensor_ => tensor_;
    rep_ => rep_;
    self_dual_ => self_dual_;
    dualizable_ => dualizable_;
    dualizable_dual_ => dualizable_dual_;
);

/// Creates a tensor-head symbol tagged with Spenso's generic tensor tag.
///
/// This expands `symbolica::symbol!` at the call site, so the symbol keeps the
/// caller's crate namespace while automatically receiving the Spenso tag. Any
/// Symbolica attributes and settings such as `print = ...` are forwarded, while
/// `tag`/`tags` remain owned by this macro so the tensor tag cannot be skipped.
#[macro_export]
macro_rules! tensor_symbol {
    ($name:ident) => {
        $crate::tensor_symbol!(stringify!($name))
    };
    ($name:ident; $($attr:ident),*) => {
        $crate::tensor_symbol!(stringify!($name); $($attr),*)
    };
    ($name:ident, $($setting:ident = $value:expr),*) => {
        $crate::tensor_symbol!(stringify!($name), $($setting = $value),*)
    };
    ($name:ident; $($attr:ident),+; $($setting:ident = $value:expr),*) => {
        $crate::tensor_symbol!(stringify!($name); $($attr),+; $($setting = $value),*)
    };
    ($id:expr) => {
        $crate::network::tags::register_tensor_symbol(
            symbolica::wrap_symbol!($id),
            Vec::new(),
            false,
        )
        .unwrap_or_else(|error| panic!("{error}"))
    };
    ($id:expr, tag = $tag:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tag = ...")
    };
    ($id:expr, tags = $tags:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tags = ...")
    };
    ($id:expr, $($setting:ident = $value:expr),*) => {
        symbolica::symbol!(
            $id,
            tag = &$crate::network::tags::SPENSO_TAG.tensor,
            $($setting = $value),*
        )
    };
    ($id:expr; $($attr:ident),*) => {
        $crate::network::tags::register_tensor_symbol(
            symbolica::wrap_symbol!($id),
            vec![$(symbolica::atom::SymbolAttribute::$attr),*],
            false,
        )
        .unwrap_or_else(|error| panic!("{error}"))
    };
    ($id:expr; $($attr:ident),+; tag = $tag:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tag = ...")
    };
    ($id:expr; $($attr:ident),+; tags = $tags:expr $(, $($rest:tt)*)?) => {
        compile_error!("tensor_symbol! owns the Spenso tensor tag; do not pass tags = ...")
    };
    ($id:expr; $($attr:ident),+; $($setting:ident = $value:expr),*) => {
        symbolica::symbol!(
            $id;
            $($attr),+;
            tag = &$crate::network::tags::SPENSO_TAG.tensor,
            $($setting = $value),*
        )
    };
}

/// Creates a tensor-head symbol tagged as a Spenso vector.
///
/// This is the rank-one tensor head constructor. It expands
/// `symbolica::symbol!` at the call site, so `vector_symbol!(p)` gets the
/// caller's namespace plus the `tensor` and `rank1` tags.
#[macro_export]
macro_rules! vector_symbol {
    ($name:ident) => {
        $crate::network::tags::register_tensor_symbol(
            symbolica::wrap_symbol!(stringify!($name)),
            Vec::new(),
            true,
        )
        .unwrap_or_else(|error| panic!("{error}"))
    };
    ($name:ident, $($setting:ident = $value:expr),* $(,)?) => {
        symbolica::symbol!(
            stringify!($name),
            tags = [
                &$crate::network::tags::SPENSO_TAG.tensor,
                &$crate::network::tags::SPENSO_TAG.rank1
            ],
            $($setting = $value),*
        )
    };
    ($name:literal) => {
        $crate::network::tags::register_tensor_symbol(
            symbolica::wrap_symbol!($name),
            Vec::new(),
            true,
        )
        .unwrap_or_else(|error| panic!("{error}"))
    };
    ($name:literal, $($setting:ident = $value:expr),* $(,)?) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.tensor,
                &$crate::network::tags::SPENSO_TAG.rank1
            ],
            $($setting = $value),*
        )
    };
}

/// Creates a representation symbol tagged with Spenso's representation tag.
///
/// This expands `symbolica::symbol!` at the call site and adds only the generic
/// representation tag.
#[macro_export]
macro_rules! representation_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tag = &$crate::network::tags::SPENSO_TAG.representation
        )
    };
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tag = &$crate::network::tags::SPENSO_TAG.representation
        )
    };
}

/// Creates a self-dual representation symbol.
///
/// This expands `symbolica::symbol!` at the call site and adds the
/// `representation` and `self_dual` tags.
#[macro_export]
macro_rules! self_dual_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.self_dual
            ]
        )
    };
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.self_dual
            ]
        )
    };
}

/// Creates a dualizable representation symbol.
///
/// This expands `symbolica::symbol!` at the call site and adds the
/// `representation` and `dualizable` tags.
#[macro_export]
macro_rules! dualizable_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.dualizable
            ]
        )
    };
    ($name:literal) => {
        symbolica::symbol!(
            $name,
            tags = [
                &$crate::network::tags::SPENSO_TAG.representation,
                &$crate::network::tags::SPENSO_TAG.dualizable
            ]
        )
    };
}

/// Creates an abstract-index symbol tagged with Spenso's index tag.
///
/// This expands `symbolica::symbol!` at the call site and adds the Spenso
/// index tag.
#[macro_export]
macro_rules! index_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tag = &$crate::network::tags::SPENSO_TAG.index
        )
    };
    ($name:literal) => {
        symbolica::symbol!($name, tag = &$crate::network::tags::SPENSO_TAG.index)
    };
}

/// Creates a function symbol tagged with Spenso's broadcast tag.
///
/// This expands `symbolica::symbol!` at the call site and adds the Spenso
/// broadcast tag.
#[macro_export]
macro_rules! broadcast_symbol {
    ($name:ident) => {
        symbolica::symbol!(
            stringify!($name),
            tag = &$crate::network::tags::SPENSO_TAG.broadcast
        )
    };
    ($name:literal) => {
        symbolica::symbol!($name, tag = &$crate::network::tags::SPENSO_TAG.broadcast)
    };
}

impl SpensoTags {
    fn print_chain_factor(value: AtomView<'_>, opt: &PrintOptions) -> Option<String> {
        fn without_visible_placeholders(value: AtomView<'_>, opt: &PrintOptions) -> Option<Atom> {
            let mut output = String::new();
            value.format(&mut output, opt, PrintState::new()).ok()?;
            let leaks_placeholder =
                [SPENSO_TAG.chain_in, SPENSO_TAG.chain_out]
                    .into_iter()
                    .any(|placeholder| {
                        let mut needle = String::new();
                        Atom::var(placeholder)
                            .as_view()
                            .format(&mut needle, opt, PrintState::new())
                            .is_ok()
                            && output.match_indices(&needle).any(|(start, _)| {
                                let end = start + needle.len();
                                let boundary = |character: char| {
                                    !character.is_alphanumeric() && character != '_'
                                };
                                output[..start].chars().next_back().is_none_or(boundary)
                                    && output[end..].chars().next().is_none_or(boundary)
                            })
                    });
            if !leaks_placeholder {
                return Some(value.to_owned());
            }

            match value {
                AtomView::Add(add) => add.iter().try_fold(Atom::Zero, |sum, term| {
                    Some(sum + without_visible_placeholders(term, opt)?)
                }),
                AtomView::Mul(mul) => mul.iter().try_fold(Atom::num(1), |product, factor| {
                    Some(product * without_visible_placeholders(factor, opt)?)
                }),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    Some(without_visible_placeholders(base, opt)?.pow(exponent.to_owned()))
                }
                AtomView::Fun(function) => {
                    let arguments = function
                        .iter()
                        .filter(|argument| {
                            !matches!(argument, AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in || var.get_symbol() == SPENSO_TAG.chain_out)
                        })
                        .map(|argument| without_visible_placeholders(argument, opt))
                        .collect::<Option<Vec<_>>>()?;
                    Some(if arguments.is_empty() {
                        Atom::var(function.get_symbol())
                    } else {
                        FunctionBuilder::new(function.get_symbol())
                            .add_args(arguments)
                            .finish()
                    })
                }
                _ => None,
            }
        }

        let sanitized = without_visible_placeholders(value, opt)?;
        let mut output = String::new();
        sanitized
            .as_view()
            .format(&mut output, opt, PrintState::new())
            .ok()?;
        Some(output)
    }

    fn print_dot(a: AtomView<'_>, opt: &PrintOptions, _state: &PrintState) -> Option<String> {
        let resolved = SpensoPrintSettings::resolve(opt)?;
        let settings = resolved.presentation;
        let parens = settings.parens;
        let AtomView::Fun(f) = a else {
            return None;
        };

        if f.get_nargs() != 2 {
            return None;
        }
        let mut argitem = f.iter();
        let a = argitem.next()?;
        let b = argitem.next()?;

        if matches!(resolved.backend, SpensoPrintBackend::Typst) {
            let a = compact_vector(a, settings, opt)?;
            let b = compact_vector(b, settings, opt)?;
            if a.representation.symbol != b.representation.symbol
                || a.representation.dimension != b.representation.dimension
                || !matches!(
                    a.representation.class,
                    RepresentationClass::SelfDual | RepresentationClass::InlineMetric
                )
            {
                return None;
            }
            return Some(format!("{} dot {}", a.label, b.label));
        }

        let mut out = String::new();
        if parens {
            out.push('(');
        }
        a.format(&mut out, opt, PrintState::new()).ok()?;
        out.push_str(match resolved.backend {
            SpensoPrintBackend::Plain => ".",
            SpensoPrintBackend::Typst => " dot ",
            SpensoPrintBackend::Latex => r"\cdot ",
        });
        b.format(&mut out, opt, PrintState::new()).ok()?;
        if parens {
            out.push(')');
        }
        Some(out)
    }

    fn print_chain(a: AtomView<'_>, opt: &PrintOptions, _state: &PrintState) -> Option<String> {
        let resolved = SpensoPrintSettings::resolve(opt)?;
        let settings = resolved.presentation;
        let AtomView::Fun(function) = a else {
            return None;
        };
        if function.get_nargs() < 2 {
            return None;
        }
        let arguments = function.iter().collect::<Vec<_>>();
        let [start, end, factors @ ..] = arguments.as_slice() else {
            return None;
        };
        let is_reserved_marker = |value| {
            is_chain_marker(value, SPENSO_TAG.chain_in)
                || is_chain_marker(value, SPENSO_TAG.chain_out)
        };

        if matches!(resolved.backend, SpensoPrintBackend::Typst) {
            let factor_source = factors
                .iter()
                .map(|factor| Self::print_chain_factor(*factor, opt))
                .collect::<Option<Vec<_>>>()?
                .join(" ");
            let mut body = if settings.parens {
                format!("lr([{factor_source}])")
            } else {
                factor_source
            };

            let mut columns = Vec::new();
            let mut prefix = None;
            let mut suffix = None;
            if is_reserved_marker(*start) {
                // Internal wiring markers are never part of the rendered endpoint.
            } else if let Some(compact) = compact_vector(*start, settings, opt) {
                prefix = Some(compact.label);
            } else if let Some(slot) = tensor_slot(*start) {
                let source = typst_index_source(slot.representation, slot.index, opt)?;
                let source = if settings.with_dim {
                    qualified_typst_index(slot, source, opt)?
                } else {
                    source
                };
                columns.push((source, slot.row));
            } else {
                prefix = Some(typst_source(*start, opt)?);
            }

            if is_reserved_marker(*end) {
                // Internal wiring markers are never part of the rendered endpoint.
            } else if let Some(compact) = compact_vector(*end, settings, opt) {
                suffix = Some(compact.label);
            } else if let Some(slot) = tensor_slot(*end) {
                let source = typst_index_source(slot.representation, slot.index, opt)?;
                let source = if settings.with_dim {
                    qualified_typst_index(slot, source, opt)?
                } else {
                    source
                };
                columns.push((source, slot.row));
            } else {
                suffix = Some(typst_source(*end, opt)?);
            }

            body = typst_attachment(body, columns);
            if let Some(prefix) = prefix {
                body = format!(r#"upright("⟨") {prefix} upright("|") {body}"#);
            }
            if let Some(suffix) = suffix {
                body = format!(r#"{body} upright("|") {suffix} upright("⟩")"#);
            }
            return Some(body);
        }

        let mut output = String::new();
        if !is_reserved_marker(*start) {
            start.format(&mut output, opt, PrintState::new()).ok()?;
        }
        if settings.parens {
            output.push('[');
        }
        for factor in factors {
            output.push_str(&Self::print_chain_factor(*factor, opt)?);
        }
        if settings.parens {
            output.push(']');
        }
        if !is_reserved_marker(*end) {
            end.format(&mut output, opt, PrintState::new()).ok()?;
        }
        Some(output)
    }

    fn print_trace(a: AtomView<'_>, opt: &PrintOptions, _state: &PrintState) -> Option<String> {
        let resolved = SpensoPrintSettings::resolve(opt)?;
        let settings = resolved.presentation;
        let AtomView::Fun(function) = a else {
            return None;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let [representation, payload @ ..] = arguments.as_slice() else {
            return None;
        };

        let factors = match payload {
            [] => Vec::new(),
            [AtomView::Fun(projector)] if projector.get_symbol() == *CYCLIC => {
                projector.iter().collect::<Vec<_>>()
            }
            // Symmetric projector traces are canonical too; preserve their
            // projector rather than presenting a raw factor list as cyclic.
            [projector @ AtomView::Fun(function)] if function.get_symbol() == *SYM => {
                vec![*projector]
            }
            _ => return None,
        };

        if matches!(resolved.backend, SpensoPrintBackend::Typst) {
            let mut head = r#"op("Tr")"#.to_owned();
            if settings.with_dim {
                head = format!("attach({head},b:{})", typst_source(*representation, opt)?);
            }
            if factors.is_empty() {
                return Some(head);
            }
            let factors = factors
                .into_iter()
                .map(|factor| Self::print_chain_factor(factor, opt))
                .collect::<Option<Vec<_>>>()?
                .join(" ");
            return Some(format!("{head} lr(({factors}))"));
        }

        let mut output = match resolved.backend {
            SpensoPrintBackend::Latex => r#"\operatorname{Tr}"#,
            SpensoPrintBackend::Plain | SpensoPrintBackend::Typst => "Tr",
        }
        .to_owned();
        if settings.with_dim {
            representation
                .format(&mut output, opt, PrintState::new())
                .ok()?;
        }
        if settings.parens {
            output.push('(');
        }
        if factors.is_empty() {
            output.push('1');
        } else {
            for factor in factors {
                output.push_str(&Self::print_chain_factor(factor, opt)?);
            }
        }
        if settings.parens {
            output.push(')');
        }
        Some(output)
    }

    fn new() -> Self {
        let broadcast = tag!("broadcast");
        let upper = tag!("upper");
        let lower = tag!("lower");
        let rank1 = tag!("rank1");
        let tensor = tag!("tensor");
        let index = tag!("index");
        let representation = tag!("representation");
        let self_dual = tag!("self_dual");
        let dualizable = tag!("dualizable");
        Self {
            chain_in: symbol!("in"),
            chain_out: symbol!("out"),
            chain: symbol!(
                "chain";Linear;
                print = Self::print_chain
            ),
            trace: symbol!(
                "trace";Linear;
                print = Self::print_trace
            ),
            rank1_: symbol!("rank1_", tags = [&tensor, &rank1], print = tensor_print),
            bracket: symbol!(
                "bracket",
                print = |a, opt, state| {
                    SpensoPrintSettings::resolve(opt)?;
                    let AtomView::Fun(f) = a else {
                        return None;
                    };
                    if f.get_nargs() == 0 {
                        return None;
                    }

                    let parenthesize = state.in_exp || state.in_exp_base;
                    let mut product_state = *state;
                    product_state.in_product = true;
                    product_state.in_sum = false;
                    product_state.level += 1;

                    let mut output = String::new();
                    if parenthesize {
                        output.push('(');
                    }
                    for (position, factor) in f.iter().enumerate() {
                        if position > 0 {
                            if opt.mode.is_latex() {
                                output.push(' ');
                            } else {
                                output.push(opt.multiplication_operator);
                            }
                        }
                        factor.format(&mut output, opt, product_state).ok()?;
                    }
                    if parenthesize {
                        output.push(')');
                    }
                    Some(output)
                }
            ),
            pure_scalar: symbol!("pure_scalar"),
            scalar: symbol!("scalar"),
            dot: symbol!("dot";Symmetric,Linear; print = Self::print_dot),
            tensor_: symbol!("tensor_", tag = tensor, print = tensor_print),
            tensor_display: symbol!(
                "tensor_display",
                print = |atom, options, state| {
                    let AtomView::Fun(wrapper) = atom else {
                        return None;
                    };
                    if wrapper.get_nargs() != 1 {
                        return None;
                    }
                    tensor_print(
                        unwrap_tensor_display(wrapper.iter().next()?),
                        options,
                        state,
                    )
                }
            ),
            i_: symbol!("i_", tag = &index),
            rep_: symbol!("rep_", tag = &representation),
            self_dual_: symbol!("self_dual_", tags = [&representation, &self_dual]),
            dualizable_: symbol!("dualizable_", tags = [&representation, &dualizable]),
            broadcast,
            upper,
            lower,
            rank1,
            tensor,
            index,
            representation,
            self_dual,
            dualizable,
        }
    }

    define_numbered_tag_family_methods! {
        pub fn rank1_ => rank1_, "rank1", rank_one_tensor_symbol;
        pub fn rep_ => rep_, "rep", representation_symbol;
    }

    pub fn tensor_symbol(&self, name: &str) -> Symbol {
        register_tensor_symbol(symbolica::wrap_symbol!(name), Vec::new(), false)
            .unwrap_or_else(|error| panic!("{error}"))
    }

    pub fn representation_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tag = &self.representation)
    }

    pub fn self_dual_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.representation, &self.self_dual])
    }

    pub fn dualizable_symbol(&self, name: &str) -> Symbol {
        symbol!(name, tags = [&self.representation, &self.dualizable])
    }

    pub fn rank_one_tensor_symbol(&self, name: &str) -> Symbol {
        register_tensor_symbol(symbolica::wrap_symbol!(name), Vec::new(), true)
            .unwrap_or_else(|error| panic!("{error}"))
    }

    define_numbered_tag_family_methods! {
        pub fn tensor_ => tensor_, "tensor", tensor_symbol;
    }

    pub fn chain<'a, 'b, 'c, A, B, F>(
        &self,
        start: A,
        end: B,
        factors: impl IntoIterator<Item = F>,
    ) -> Atom
    where
        A: Into<AtomOrView<'a>>,
        B: Into<AtomOrView<'b>>,
        F: Into<AtomOrView<'c>>,
    {
        let mut f = FunctionBuilder::new(self.chain).add_arg(start).add_arg(end);
        for factor in factors {
            f = f.add_arg(factor);
        }
        f.finish()
    }

    pub fn trace<'a, 'b, R, F>(&self, rep: R, factors: impl IntoIterator<Item = F>) -> Atom
    where
        R: Into<AtomOrView<'a>>,
        F: Into<AtomOrView<'b>>,
    {
        let mut f = FunctionBuilder::new(self.trace).add_arg(rep);
        for factor in factors {
            f = f.add_arg(factor);
        }
        f.finish()
    }

    pub fn reverse_flip_factor(&self, factor: AtomView<'_>) -> Atom {
        let tmp = symbol!("spenso::chain_flip_tmp");
        factor
            .to_owned()
            .replace(self.chain_in)
            .with(tmp)
            .replace(self.chain_out)
            .with(self.chain_in)
            .replace(tmp)
            .with(self.chain_out)
    }

    pub fn reverse_flip_factors(&self, factors: impl IntoIterator<Item = Atom>) -> Vec<Atom> {
        let mut factors = factors.into_iter().collect::<Vec<_>>();
        factors.reverse();
        factors
            .into_iter()
            .map(|factor| self.reverse_flip_factor(factor.as_view()))
            .collect()
    }

    define_numbered_tag_family_methods! {
        pub fn self_dual_ => self_dual_, "self_dual", self_dual_symbol;
        pub fn dualizable_ => dualizable_, "dualizable", dualizable_symbol;
    }

    pub fn dualizable_dual_<'a, const N: usize, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        AIND_SYMBOLS.dual(self.dualizable_::<N, A>(args))
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{
        atom::{
            Atom, AtomCore, AtomView, FunctionBuilder, NamespacedSymbol, Symbol, SymbolAttribute,
            SymbolBuilder, UserData,
        },
        function,
        printer::{PrintOptions, PrintState},
        symbol, wrap_symbol,
    };

    use crate::{
        cyclic, dind, lor, mink,
        shadowing::symbolica_utils::SpensoPrintSettings,
        structure::representation::{IndexDisplay, IndexPalette, IndexRow, LibraryRep},
    };

    use super::{
        ChainMarkers, ChainOrientation, SPENSO_TAG, SpensoTags, SymbolAtomExt, chain_markers,
        prepare_tensor_print, register_tensor_symbol, typst_tensor_head,
    };

    fn typst_options(settings: SpensoPrintSettings) -> PrintOptions {
        PrintOptions {
            custom_print_mode: settings.into(),
            ..PrintOptions::typst()
        }
    }

    fn palette_metadata() -> UserData {
        let label = UserData::List(vec![
            UserData::String("symbol".to_owned()),
            UserData::String("M".to_owned()),
        ]);
        let palette_label = |name: &str| {
            UserData::List(vec![
                UserData::String("symbol".to_owned()),
                UserData::String(name.to_owned()),
            ])
        };
        UserData::List(vec![
            UserData::String("spenso::representation-v1".to_owned()),
            UserData::String("self-dual".to_owned()),
            label,
            UserData::List(vec![
                UserData::String("cyclic".to_owned()),
                UserData::Integer(1),
                UserData::List(vec![palette_label("mu"), palette_label("nu")]),
            ]),
        ])
    }

    fn palette_representation(index: i64) -> Atom {
        let representation =
            SymbolBuilder::new(wrap_symbol!("spenso_typst_tests::palette_representation"))
                .with_tags([
                    SPENSO_TAG.representation.clone(),
                    SPENSO_TAG.self_dual.clone(),
                ])
                .with_user_data(palette_metadata())
                .build()
                .unwrap();
        FunctionBuilder::new(representation)
            .add_arg(Atom::num(4))
            .add_arg(Atom::num(index))
            .finish()
    }

    fn bottom_representation(index: Option<Symbol>) -> Atom {
        let representation = LibraryRep::new_self_dual_with_index_palette_and_row(
            "spenso_typst_tests::B",
            IndexPalette::Numeric,
            IndexRow::Bottom,
        )
        .unwrap();
        let mut atom = FunctionBuilder::new(representation.symbol()).add_arg(Atom::num(4));
        if let Some(index) = index {
            atom = atom.add_arg(Atom::var(index));
        }
        atom.finish()
    }

    fn gamma_symbol() -> Symbol {
        register_tensor_symbol(
            NamespacedSymbol::parse("spenso::gamma"),
            vec![SymbolAttribute::Linear],
            false,
        )
        .unwrap()
    }

    fn compact_vector_atom(symbol: Symbol, label: i64, representation: Atom) -> Atom {
        FunctionBuilder::new(symbol)
            .add_arg(Atom::num(label))
            .add_arg(representation)
            .finish()
    }

    fn compact_vector_atom_rep_first(symbol: Symbol, label: i64, representation: Atom) -> Atom {
        FunctionBuilder::new(symbol)
            .add_arg(representation)
            .add_arg(Atom::num(label))
            .finish()
    }

    #[test]
    fn numbered_wildcard_macros_build_variables_without_args() {
        let expr = rank1_!(0);
        let AtomView::Var(var) = expr.as_view() else {
            panic!("empty wildcard head should be a variable");
        };

        assert_eq!(var.get_symbol(), SPENSO_TAG.rank1_);
    }

    #[test]
    fn numbered_wildcard_macros_build_functions_with_args() {
        let expr = rank1_!(
            1;
            Atom::var(symbol!("a___")),
            rep_!(2; Atom::var(symbol!("d_")))
        );

        let AtomView::Fun(fun) = expr.as_view() else {
            panic!("wildcard head with args should be a function");
        };

        assert_eq!(
            fun.get_symbol(),
            SPENSO_TAG.rank_one_tensor_symbol("rank11_")
        );
        assert_eq!(fun.get_nargs(), 2);
    }

    #[test]
    fn numbered_representation_families_use_their_own_prefixes() {
        let self_dual = self_dual_!(1; Atom::var(symbol!("d_")));
        let dualizable = dualizable_!(1; Atom::var(symbol!("d_")));

        let AtomView::Fun(self_dual) = self_dual.as_view() else {
            panic!("self-dual wildcard should be a function");
        };
        let AtomView::Fun(dualizable) = dualizable.as_view() else {
            panic!("dualizable wildcard should be a function");
        };

        assert_eq!(
            self_dual.get_symbol(),
            SPENSO_TAG.self_dual_symbol("self_dual1_")
        );
        assert_eq!(
            dualizable.get_symbol(),
            SPENSO_TAG.dualizable_symbol("dualizable1_")
        );
    }

    #[test]
    fn symbol_atom_ext_uses_variable_for_empty_args() {
        let symbol = SPENSO_TAG.tensor_symbol("empty_tensor_pattern");

        assert_eq!(
            symbol.atom_with_args(std::iter::empty::<Atom>()),
            Atom::var(symbol)
        );
        assert!(matches!(
            tensor_!(0; Atom::var(symbol!("a___"))).as_view(),
            AtomView::Fun(_)
        ));
    }

    #[test]
    fn trace_uses_typst_operator_only_in_typst_print_mode() {
        let trace = SPENSO_TAG.trace(Atom::var(symbol!("rep")), [cyclic!(Atom::num(1))]);

        assert_eq!(
            trace
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            r#"op("Tr") lr((1))"#
        );
        assert_eq!(
            trace
                .printer(SpensoPrintSettings::compact().nice_symbolica())
                .to_string(),
            "Tr(1)"
        );
        assert_eq!(
            trace
                .printer(SpensoPrintSettings::typst().nice_symbolica())
                .to_string(),
            "Tr(1)"
        );
    }
    #[test]
    fn chain_hides_namespaced_endpoint_placeholders() {
        assert_eq!(SPENSO_TAG.chain_in.get_namespace(), "spenso");
        assert_eq!(SPENSO_TAG.chain_out.get_namespace(), "spenso");

        let chain = SPENSO_TAG.chain(
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            [Atom::var(SPENSO_TAG.tensor_symbol("factor"))],
        );

        let compact = chain
            .printer(SpensoPrintSettings::compact().nice_symbolica())
            .to_string();
        let typst = chain.printer(PrintOptions::typst()).to_string();

        assert_eq!(compact, "[factor]");
        assert_eq!(typst, "lr([\"factor\"])");
        assert!(!compact.contains("in"));
        assert!(!compact.contains("out"));
    }

    #[test]
    fn reverse_flip_swaps_namespaced_endpoint_placeholders() {
        let factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("flip_factor"))
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();
        let flipped = SPENSO_TAG.reverse_flip_factor(factor.as_view());
        let AtomView::Fun(flipped) = flipped.as_view() else {
            panic!("expected a tensor factor")
        };
        let arguments = flipped.iter().collect::<Vec<_>>();

        assert!(
            matches!(arguments[0], AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_out)
        );
        assert!(
            matches!(arguments[1], AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in)
        );
    }

    #[test]
    fn bracket_prints_as_an_ordered_product() {
        let bracket = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(Atom::var(symbol!("z")))
            .add_arg(Atom::var(symbol!("a")))
            .finish();
        let mut compact = SpensoPrintSettings::compact().nice_symbolica();
        compact.color_builtin_symbols = false;

        assert_eq!(bracket.printer(compact).to_string(), "z·a");
        assert_eq!(bracket.printer(PrintOptions::typst()).to_string(), "z a");
    }

    #[test]
    fn tensors_use_physica_style_typst_attachment_columns() {
        let head = crate::tensor_symbol!("spenso_typst_tests::T");
        let tensor = function!(
            head,
            Atom::num(1),
            mink!(4, symbol!("mu")),
            dind!(lor!(4, symbol!("nu")))
        );

        assert_eq!(
            prepare_tensor_print(&tensor)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($T$,std.hide($zws$)).join(),t:std.hide(1) mu std.hide(nu),b:1 std.hide(mu) nu)"
        );
    }

    #[test]
    fn vectors_default_self_dual_slots_to_the_top_row() {
        let head = crate::vector_symbol!("spenso_typst_tests::p");
        let vector = function!(head, mink!(4, symbol!("rho")));

        assert_eq!(
            prepare_tensor_print(&vector)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($p$,std.hide($zws$)).join(),t:rho,b:std.hide(rho))"
        );
    }

    #[test]
    fn metric_uses_tensor_attachments_in_spenso_typst_mode() {
        let metric = crate::network::library::symbolic::ETS
            .metric(mink!(4, symbol!("mu")), mink!(4, symbol!("nu")));

        assert_eq!(
            metric
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($g$,std.hide($zws$)).join(),t:mu nu,b:std.hide(mu) std.hide(nu))"
        );
    }

    #[test]
    fn tagged_imports_can_recover_the_generic_tensor_printer() {
        let head = SymbolBuilder::new(wrap_symbol!("spenso_typst_tests::R"))
            .with_tags([SPENSO_TAG.tensor.clone()])
            .build()
            .unwrap();
        let tensor = FunctionBuilder::new(head)
            .add_arg(mink!(4, symbol!("sigma")))
            .finish();
        assert!(head.get_print_function().is_none());

        assert_eq!(
            prepare_tensor_print(&tensor)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($R$,std.hide($zws$)).join(),t:sigma,b:std.hide(sigma))"
        );
    }

    #[test]
    fn typst_preparation_wraps_nested_tensors_once() {
        let inner_head = SymbolBuilder::new(wrap_symbol!("spenso_typst_tests::NestedInner"))
            .with_tags([SPENSO_TAG.tensor.clone()])
            .build()
            .unwrap();
        let outer_head = SymbolBuilder::new(wrap_symbol!("spenso_typst_tests::NestedOuter"))
            .with_tags([SPENSO_TAG.tensor.clone()])
            .build()
            .unwrap();
        let inner = function!(inner_head, mink!(4, symbol!("mu")));
        let outer = function!(outer_head, inner.clone(), mink!(4, symbol!("nu")));

        let wrapped_inner = FunctionBuilder::new(SPENSO_TAG.tensor_display)
            .add_arg(inner)
            .finish();
        let wrapped_outer = FunctionBuilder::new(SPENSO_TAG.tensor_display)
            .add_arg(function!(
                outer_head,
                wrapped_inner,
                mink!(4, symbol!("nu"))
            ))
            .finish();
        let prepared = prepare_tensor_print(&outer);

        assert_eq!(prepared, wrapped_outer);
        assert_eq!(prepare_tensor_print(&prepared), prepared);
        let printed = prepared
            .printer(SpensoPrintSettings::typst_options())
            .to_string();
        assert!(
            printed.starts_with(r#"attach(#($italic("NestedOuter")$"#),
            "{printed}"
        );
        assert!(
            printed.contains(r#"attach(#($italic("NestedInner")$"#),
            "{printed}"
        );
    }

    #[test]
    fn punctuation_tensor_heads_are_quoted_as_typst_strings() {
        let head = SymbolBuilder::new(wrap_symbol!("spenso_typst_tests::#"))
            .with_tags([SPENSO_TAG.tensor.clone()])
            .build()
            .unwrap();

        assert_eq!(typst_tensor_head(head), r##"italic("#")"##);
    }

    #[test]
    fn tensor_indices_resolve_the_representations_fixed_palette() {
        let head = crate::tensor_symbol!("spenso_typst_tests::PaletteTensor");
        let tensor = function!(head, palette_representation(3));

        assert_eq!(
            prepare_tensor_print(&tensor)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($italic(\"PaletteTensor\")$,std.hide($zws$)).join(),t:attach(mu,b:1),b:std.hide(attach(mu,b:1)))"
        );
    }

    #[test]
    fn tensor_dimension_qualification_decorates_the_index_inside_its_variance() {
        let head = crate::tensor_symbol!("spenso_typst_tests::QualifiedTensor");
        let tensor = function!(head, palette_representation(3));
        let mut settings = SpensoPrintSettings::typst();
        settings.with_dim = true;

        assert_eq!(
            prepare_tensor_print(&tensor)
                .printer(typst_options(settings))
                .to_string(),
            "attach(#($italic(\"QualifiedTensor\")$,std.hide($zws$)).join(),t:attach(attach(mu,b:1),t:attach(M,b:4)),b:std.hide(attach(attach(mu,b:1),t:attach(M,b:4))))"
        );
    }

    #[test]
    fn tensor_indices_use_portable_manual_display_metadata() {
        let index_display = IndexDisplay::symbol("mu")
            .unwrap()
            .with_bottom(IndexDisplay::Number(1));
        let index = SymbolBuilder::new(wrap_symbol!("spenso_typst_tests::manual_index"))
            .with_tags([SPENSO_TAG.index.clone()])
            .with_user_data(index_display.symbol_user_data())
            .build()
            .unwrap();
        let head = crate::vector_symbol!("spenso_typst_tests::manual_vector");
        let vector = function!(head, mink!(4, Atom::var(index)));

        assert_eq!(
            prepare_tensor_print(&vector)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($italic(\"manual_vector\")$,std.hide($zws$)).join(),t:attach(mu,b:1),b:std.hide(attach(mu,b:1)))"
        );
    }

    #[test]
    fn representation_owned_rows_interleave_in_original_argument_order() {
        let head = crate::tensor_symbol!("spenso_typst_tests::Interleaved");
        let tensor = function!(
            head,
            mink!(4, symbol!("mu")),
            bottom_representation(Some(symbol!("a"))),
            mink!(4, symbol!("nu")),
            bottom_representation(Some(symbol!("b")))
        );

        assert_eq!(
            prepare_tensor_print(&tensor)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($italic(\"Interleaved\")$,std.hide($zws$)).join(),t:mu std.hide(a) nu std.hide(b),b:std.hide(mu) a std.hide(nu) b)"
        );
    }

    #[test]
    fn chain_markers_are_found_anywhere_and_reversal_is_transposed() {
        let head = crate::tensor_symbol!("spenso_typst_tests::MarkerTensor");
        let forward = function!(
            head,
            mink!(4, symbol!("mu")),
            Atom::var(SPENSO_TAG.chain_in),
            bottom_representation(Some(symbol!("a"))),
            Atom::var(SPENSO_TAG.chain_out),
            mink!(4, symbol!("nu"))
        );
        let reversed = SPENSO_TAG.reverse_flip_factor(forward.as_view());
        let body = "attach(#($italic(\"MarkerTensor\")$,std.hide($zws$)).join(),t:mu std.hide(a) nu,b:std.hide(mu) a std.hide(nu))";

        assert_eq!(
            prepare_tensor_print(&forward)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            body
        );
        assert_eq!(
            prepare_tensor_print(&reversed)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            format!("attach({body},t:upright(\"T\"))")
        );

        let duplicate_input = [
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            Atom::var(SPENSO_TAG.chain_in),
        ];
        let duplicate_views = duplicate_input
            .iter()
            .map(Atom::as_view)
            .collect::<Vec<_>>();
        assert_eq!(chain_markers(&duplicate_views), None);

        let valid_arguments = if let AtomView::Fun(function) = forward.as_view() {
            function.iter().collect::<Vec<_>>()
        } else {
            unreachable!()
        };
        assert_eq!(
            chain_markers(&valid_arguments),
            Some(ChainMarkers {
                input: 1,
                output: 3,
                orientation: ChainOrientation::Forward,
            })
        );
    }

    #[test]
    fn gamma_uses_bispinor_rows_and_schoonschip_slash_notation() {
        let gamma = gamma_symbol();
        let explicit = function!(
            gamma,
            bottom_representation(Some(symbol!("a"))),
            bottom_representation(Some(symbol!("b"))),
            mink!(4, symbol!("mu"))
        );
        assert_eq!(
            prepare_tensor_print(&explicit)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($gamma$,std.hide($zws$)).join(),t:std.hide(a) std.hide(b) mu,b:a b std.hide(mu))"
        );

        let p = crate::vector_symbol!("spenso_typst_tests::p");
        let factor = function!(
            gamma,
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            compact_vector_atom_rep_first(p, 1, mink!(4))
        );
        assert_eq!(
            prepare_tensor_print(&factor)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "cancel(attach(#($p$,std.hide($zws$)).join(),t:std.hide(1),b:1))"
        );
    }

    #[test]
    fn compact_vectors_render_as_dot_products_and_chain_endpoints() {
        let p = crate::vector_symbol!("spenso_typst_tests::dot_p");
        let q = crate::vector_symbol!("spenso_typst_tests::dot_q");
        let left = compact_vector_atom(p, 1, mink!(4));
        let right = compact_vector_atom_rep_first(q, 2, mink!(4));
        let dot = function!(SPENSO_TAG.dot, left, right);

        assert_eq!(
            prepare_tensor_print(&dot)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($italic(\"dot_p\")$,std.hide($zws$)).join(),t:std.hide(1),b:1) dot attach(#($italic(\"dot_q\")$,std.hide($zws$)).join(),t:std.hide(2),b:2)"
        );

        let u = crate::vector_symbol!("spenso_typst_tests::u");
        let v = crate::vector_symbol!("spenso_typst_tests::v");
        let gamma = gamma_symbol();
        let factor = function!(
            gamma,
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            mink!(4, symbol!("mu"))
        );
        let chain = SPENSO_TAG.chain(
            compact_vector_atom(u, 1, bottom_representation(None)),
            compact_vector_atom_rep_first(v, 2, bottom_representation(None)),
            [factor],
        );
        assert_eq!(
            prepare_tensor_print(&chain)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "upright(\"⟨\") attach(#($u$,std.hide($zws$)).join(),t:std.hide(1),b:1) upright(\"|\") lr([attach(#($gamma$,std.hide($zws$)).join(),t:mu,b:std.hide(mu))]) upright(\"|\") attach(#($v$,std.hide($zws$)).join(),t:std.hide(2),b:2) upright(\"⟩\")"
        );
    }

    #[test]
    fn dot_notation_requires_equal_self_dual_representations() {
        let p = crate::vector_symbol!("spenso_typst_tests::checked_dot_p");
        let q = crate::vector_symbol!("spenso_typst_tests::checked_dot_q");
        let options = SpensoPrintSettings::typst_options();
        let state = PrintState::new();
        let invalid_pairs = [
            (
                compact_vector_atom(p, 1, mink!(4)),
                compact_vector_atom(q, 2, mink!(5)),
            ),
            (
                compact_vector_atom(p, 1, mink!(4)),
                compact_vector_atom(q, 2, bottom_representation(None)),
            ),
            (
                compact_vector_atom(p, 1, lor!(4)),
                compact_vector_atom(q, 2, lor!(4)),
            ),
        ];

        for (left, right) in invalid_pairs {
            let dot = function!(SPENSO_TAG.dot, left, right);
            assert!(SpensoTags::print_dot(dot.as_view(), &options, &state).is_none());
        }
    }

    #[test]
    fn typst_preparation_overrides_callbacks_and_preserves_idenso_heads() {
        let callback_head = SymbolBuilder::new(wrap_symbol!("spenso_typst_tests::CallbackTensor"))
            .with_tags([SPENSO_TAG.tensor.clone()])
            .with_print_function(|_, _, _| Some("wrong-callback".to_owned()))
            .build()
            .unwrap();
        let tensor = function!(callback_head, mink!(4, symbol!("mu")));
        assert_eq!(
            prepare_tensor_print(&tensor)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "attach(#($italic(\"CallbackTensor\")$,std.hide($zws$)).join(),t:mu,b:std.hide(mu))"
        );

        for (name, expected) in [
            ("spenso::projp", "ℙ_p"),
            ("spenso::projm", "ℙ_m"),
            ("spenso::gamma5", "gamma_5"),
            ("spenso::gamma0", "gamma_0"),
        ] {
            let symbol = SymbolBuilder::new(NamespacedSymbol::parse(name))
                .with_tags([SPENSO_TAG.tensor.clone()])
                .build()
                .unwrap();
            assert_eq!(typst_tensor_head(symbol), expected);
        }
    }

    #[test]
    fn canonical_trace_flattens_only_the_cyclic_factor_container() {
        let gamma = gamma_symbol();
        let gamma_mu = function!(
            gamma,
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            mink!(4, symbol!("mu"))
        );
        let gamma_nu = function!(
            gamma,
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            mink!(4, symbol!("nu"))
        );
        let trace = SPENSO_TAG.trace(mink!(4), [cyclic!(gamma_mu, gamma_nu)]);

        assert_eq!(
            prepare_tensor_print(&trace)
                .printer(SpensoPrintSettings::typst_options())
                .to_string(),
            "op(\"Tr\") lr((attach(#($gamma$,std.hide($zws$)).join(),t:mu,b:std.hide(mu)) attach(#($gamma$,std.hide($zws$)).join(),t:nu,b:std.hide(nu))))"
        );
    }
}
