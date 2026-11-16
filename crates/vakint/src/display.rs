//! Display forms of vakint's symbols (LaTeX, MathML, Typst).
//!
//! One table, `display_forms.json`, lists every symbol vakint defines or emits:
//! `forms` (a dedicated display form per print mode), `native` (Symbolica's
//! default printing is already the conventional form, with the reason) and
//! `internal` (never part of a user-visible result). The same table drives the
//! Python twin `tools/display/vakint_display.py`, which a Python session (the
//! benchmark dashboard) calls before parsing vakint output, so both render alike.
//!
//! Only the LaTeX, MathML and Typst modes are customised. Plain (Symbolica),
//! Mathematica and Sympy printing, and therefore every canonical string, file
//! and cache key, are unchanged: the print callbacks return `None` there.
//!
//! Symbolica fixes a symbol's print function when the symbol is first created,
//! and orders symbols (hence factors in every printed string) by creation order.
//! So a decorated symbol is created exactly where and when vakint created it
//! before: `vk_symbol!` calls [`create_decorated`] instead of parsing the name
//! when the name has a Rust form and does not exist yet. Nothing is registered
//! ahead of time, so the creation order of all symbols, and every plain and
//! canonical string, is the same as without display forms (tested by
//! `display::tests::display_forms_keep_the_symbol_creation_order`). A symbol
//! first created by parsing a longer text stays undecorated. The environment
//! variable `VAKINT_NO_DISPLAY_FORMS` disables the decoration.
//!
//! Scalar products: the reduction keeps `vakint::dot(p(i), q(j))`; it displays as
//! `p_i·q_j` (or `p_i^2`). [`to_spenso_dot_notation`] converts them to spenso's
//! `spenso::dot(p(i, mink(D)), q(j, mink(D)))`, the form gammaloop uses, and
//! [`from_spenso_dot_notation`] converts back (round trip tested).

use std::sync::LazyLock;

use serde::Deserialize;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::dimension::Dimension;
use spenso::structure::representation::{BaseRepName, Minkowski, RepName};
use symbolica::atom::{
    Atom, AtomCore, AtomView, FunctionBuilder, NamespacedSymbol, Symbol, SymbolAttribute,
    SymbolBuilder,
};
use symbolica::printer::{PrintMode, PrintOptions, PrintState};

/// The display table (schema `vakint.display_forms.v1`).
pub const DISPLAY_FORMS_JSON: &str = include_str!("display_forms.json");

#[derive(Debug, Clone, Deserialize)]
pub struct DisplayTable {
    pub schema: String,
    pub forms: Vec<DisplayForm>,
    pub native: Vec<NativeEntry>,
    pub internal: Vec<InternalEntry>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct DisplayForm {
    pub name: String,
    pub owner: String,
    /// Registered by vakint itself (vakint-owned symbols only).
    pub rust: bool,
    /// Number of arguments, `0` for a variable, `"*"` for any number.
    pub arity: serde_json::Value,
    #[serde(default)]
    pub attributes: Vec<String>,
    pub meaning: String,
    pub sample: String,
    pub latex: String,
    pub mathml: String,
    pub typst: String,
    #[serde(default)]
    pub when_args_equal: Option<ModeTemplates>,
    /// Wrap the form in parentheses when it is the base of a power.
    #[serde(default)]
    pub parenthesise_as_power_base: bool,
    /// Wrap a sign-led form in parentheses when it is a factor of a product.
    #[serde(default)]
    pub parenthesise_in_product: bool,
}

#[derive(Debug, Clone, Deserialize)]
pub struct ModeTemplates {
    pub latex: String,
    pub mathml: String,
    pub typst: String,
}

#[derive(Debug, Clone, Deserialize)]
pub struct NativeEntry {
    pub name: String,
    pub meaning: String,
    pub reason: String,
}

#[derive(Debug, Clone, Deserialize)]
pub struct InternalEntry {
    pub name: String,
    pub reason: String,
}

/// Above this stored size an atom keeps Symbolica's native (pageable) MathML.
pub const MATHML_MAX_ATOM_BYTES: usize = 1024;
/// Above this length a MathML form is dropped for the native printing (the
/// explorer's limit for one custom fragment is 4096 bytes including wrappers).
pub const MATHML_MAX_FRAGMENT_BYTES: usize = 3500;

pub static DISPLAY_TABLE: LazyLock<DisplayTable> = LazyLock::new(|| {
    serde_json::from_str(DISPLAY_FORMS_JSON).expect("display_forms.json is valid")
});

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Latex,
    MathMl,
    Typst,
}

impl Mode {
    fn of(options: &PrintOptions) -> Option<Self> {
        match options.mode {
            PrintMode::Latex => Some(Self::Latex),
            PrintMode::MathML => Some(Self::MathMl),
            PrintMode::Typst => Some(Self::Typst),
            _ => None,
        }
    }

    fn plus(self) -> &'static str {
        match self {
            Self::MathMl => "<mo>+</mo>",
            _ => "+",
        }
    }

    fn minus(self) -> &'static str {
        match self {
            Self::MathMl => "<mo>−</mo>",
            _ => "-",
        }
    }

    fn comma(self) -> &'static str {
        match self {
            Self::MathMl => "<mo>,</mo>",
            _ => ",",
        }
    }

    fn parenthesise(self, inner: &str) -> String {
        match self {
            Self::Latex => format!(r"\left({inner}\right)"),
            Self::MathMl => format!("<mrow><mo>(</mo>{inner}<mo>)</mo></mrow>"),
            Self::Typst => format!("({inner})"),
        }
    }

    fn subscript(self, base: &str, sub: &str) -> String {
        match self {
            Self::Latex => format!("{base}_{{{sub}}}"),
            Self::MathMl => format!("<msub><mi>{base}</mi>{sub}</msub>"),
            Self::Typst => format!("{base}_({sub})"),
        }
    }
}

fn render(atom: AtomView<'_>, options: &PrintOptions, mode: Mode) -> Option<String> {
    let mut out = String::new();
    atom.format(&mut out, options, PrintState::new()).ok()?;
    Some(if mode == Mode::MathMl {
        format!("<mrow>{out}</mrow>")
    } else {
        out
    })
}

fn is_negative_number(atom: AtomView<'_>) -> bool {
    matches!(atom, AtomView::Num(_)) && atom.to_canonical_string().starts_with('-')
}

fn is_composite(atom: AtomView<'_>) -> bool {
    match atom {
        AtomView::Add(_) => true,
        AtomView::Num(_) => is_negative_number(atom),
        AtomView::Mul(product) => product.iter().any(is_negative_number),
        _ => false,
    }
}

/// `+`, `-`, `+n`, `-n` or `+(x)` for a coefficient written before a factor.
fn signed(atom: AtomView<'_>, options: &PrintOptions, mode: Mode) -> Option<String> {
    if let AtomView::Num(_) = atom {
        let negative = is_negative_number(atom);
        let magnitude = if negative {
            -atom.to_owned()
        } else {
            atom.to_owned()
        };
        let sign = if negative { mode.minus() } else { mode.plus() };
        if magnitude == Atom::num(1) {
            return Some(sign.into());
        }
        return Some(format!(
            "{sign}{}",
            render(magnitude.as_view(), options, mode)?
        ));
    }
    let inner = render(atom, options, mode)?;
    Some(format!(
        "{}{}",
        mode.plus(),
        if is_composite(atom) || matches!(atom, AtomView::Mul(_)) {
            mode.parenthesise(&inner)
        } else {
            inner
        }
    ))
}

fn parenthesised(atom: AtomView<'_>, options: &PrintOptions, mode: Mode) -> Option<String> {
    let inner = render(atom, options, mode)?;
    Some(if is_composite(atom) {
        mode.parenthesise(&inner)
    } else {
        inner
    })
}

/// A one-argument vector `p(0)` as `p_0`; anything else as [`parenthesised`].
fn vector(atom: AtomView<'_>, options: &PrintOptions, mode: Mode) -> Option<String> {
    if let AtomView::Fun(function) = atom
        && function.get_nargs() == 1
    {
        let label = render(function.iter().next()?, options, mode)?;
        return Some(mode.subscript(function.get_symbol().get_stripped_name(), &label));
    }
    parenthesised(atom, options, mode)
}

fn join(arguments: &[AtomView<'_>], options: &PrintOptions, mode: Mode) -> Option<String> {
    let parts = arguments
        .iter()
        .map(|argument| render(*argument, options, mode))
        .collect::<Option<Vec<_>>>()?;
    Some(parts.join(mode.comma()))
}

/// Expand one template (placeholders documented in `display_forms.json`).
fn expand(
    template: &str,
    arguments: &[AtomView<'_>],
    options: &PrintOptions,
    mode: Mode,
) -> Option<String> {
    let mut out = String::new();
    let mut rest = template;
    while let Some(start) = rest.find('⟨') {
        out.push_str(&rest[..start]);
        let after = &rest[start + '⟨'.len_utf8()..];
        let end = after.find('⟩')?;
        let key = &after[..end];
        rest = &after[end + '⟩'.len_utf8()..];
        let piece = match key {
            "*" => join(arguments, options, mode)?,
            "init" => join(&arguments[..arguments.len().checked_sub(1)?], options, mode)?,
            "last" => render(*arguments.last()?, options, mode)?,
            _ if key.ends_with('*') => {
                let from = key[..key.len() - 1].parse::<usize>().ok()?;
                join(arguments.get(from..)?, options, mode)?
            }
            _ => {
                let (index, modifier) = key.split_once(':').unwrap_or((key, ""));
                let argument = *arguments.get(index.parse::<usize>().ok()?)?;
                match modifier {
                    "" => render(argument, options, mode)?,
                    "+" => signed(argument, options, mode)?,
                    "p" => parenthesised(argument, options, mode)?,
                    "vec" => vector(argument, options, mode)?,
                    _ => return None,
                }
            }
        };
        out.push_str(&piece);
    }
    out.push_str(rest);
    Some(out)
}

impl DisplayForm {
    fn template(&self, mode: Mode, equal_arguments: bool) -> &str {
        let templates = match (&self.when_args_equal, equal_arguments) {
            (Some(equal), true) => (&equal.latex, &equal.mathml, &equal.typst),
            _ => (&self.latex, &self.mathml, &self.typst),
        };
        match mode {
            Mode::Latex => templates.0,
            Mode::MathMl => templates.1,
            Mode::Typst => templates.2,
        }
    }

    fn accepts(&self, arguments: Option<usize>) -> bool {
        match (&self.arity, arguments) {
            (serde_json::Value::Number(n), None) => n.as_u64() == Some(0),
            (serde_json::Value::Number(n), Some(count)) => n.as_u64() == Some(count as u64),
            (serde_json::Value::String(s), Some(count)) => s == "*" && count > 0,
            _ => false,
        }
    }

    /// The display string of `atom` (whose head is this form's symbol), or
    /// `None` for Symbolica's default printing.
    pub fn print(
        &self,
        atom: AtomView<'_>,
        options: &PrintOptions,
        state: &PrintState,
    ) -> Option<String> {
        let mode = Mode::of(options)?;
        // Symbolica's MathML explorer shows a custom fragment as one atomic
        // block (a hole above 4096 bytes): a large argument (a long coef(...))
        // keeps the native, pageable printing.
        if mode == Mode::MathMl && atom.get_byte_size() > MATHML_MAX_ATOM_BYTES {
            return None;
        }
        let out = self.print_bare(atom, options, mode)?;
        if mode == Mode::MathMl && out.len() > MATHML_MAX_FRAGMENT_BYTES {
            return None;
        }
        // A form never carries enclosing parentheses itself: they come only from
        // the printer state (power base, or a sign-led form inside a product).
        Some(
            if (self.parenthesise_as_power_base && state.in_exp_base)
                || (self.parenthesise_in_product && state.in_product)
            {
                mode.parenthesise(&out)
            } else {
                out
            },
        )
    }

    fn print_bare(&self, atom: AtomView<'_>, options: &PrintOptions, mode: Mode) -> Option<String> {
        match atom {
            AtomView::Var(_) if self.accepts(None) => {
                expand(self.template(mode, false), &[], options, mode)
            }
            AtomView::Fun(function) if self.accepts(Some(function.get_nargs())) => {
                let arguments = function.iter().collect::<Vec<_>>();
                let equal = arguments.len() > 1 && arguments.windows(2).all(|w| w[0] == w[1]);
                expand(self.template(mode, equal), &arguments, options, mode)
            }
            _ => None,
        }
    }

    fn symbol_attributes(&self) -> Vec<SymbolAttribute> {
        self.attributes
            .iter()
            .map(|attribute| match attribute.as_str() {
                "Symmetric" => SymbolAttribute::Symmetric,
                "Antisymmetric" => SymbolAttribute::Antisymmetric,
                "Linear" => SymbolAttribute::Linear,
                "Scalar" => SymbolAttribute::Scalar,
                "Real" => SymbolAttribute::Real,
                "Integer" => SymbolAttribute::Integer,
                "Positive" => SymbolAttribute::Positive,
                other => panic!("display_forms.json: unknown attribute {other}"),
            })
            .collect()
    }

    /// A builder for this symbol with its attributes and display form; foreign
    /// owners may use it to register the form at their own definition site.
    pub fn symbol_builder(&'static self) -> SymbolBuilder {
        let builder = SymbolBuilder::new(NamespacedSymbol::parse(&self.name))
            .with_print_function(move |atom, options, state| self.print(atom, options, state));
        if self.attributes.is_empty() {
            builder
        } else {
            builder.with_attributes(self.symbol_attributes())
        }
    }
}

/// The display form of a symbol by its full name (`namespace::name`).
pub fn display_form(name: &str) -> Option<&'static DisplayForm> {
    DISPLAY_TABLE.forms.iter().find(|form| form.name == name)
}

static RUST_FORMS: LazyLock<std::collections::HashMap<&'static str, &'static DisplayForm>> =
    LazyLock::new(|| {
        if std::env::var_os("VAKINT_NO_DISPLAY_FORMS").is_some() {
            return Default::default();
        }
        DISPLAY_TABLE
            .forms
            .iter()
            .filter(|form| form.rust)
            .map(|form| (form.name.as_str(), form))
            .collect()
    });

/// Create the symbol `full_name` (`namespace::name`) with its display form, if
/// it has a Rust form and does not exist yet; `None` otherwise (the caller then
/// creates or looks it up as before). Called by `vk_symbol!` at the point where
/// it would have created the symbol, so the creation order is unchanged.
pub fn create_decorated(full_name: &str) -> Option<Symbol> {
    let form = RUST_FORMS.get(full_name)?;
    if Symbol::get_symbol(NamespacedSymbol::parse(full_name)).is_some() {
        return None;
    }
    form.symbol_builder().build().ok()
}

/// Whether `symbol` carries a display form (it was created by vakint first).
pub fn is_decorated(symbol: Symbol) -> bool {
    symbol.get_print_function().is_some()
}

fn minkowski(dimension: impl Into<Dimension>) -> Atom {
    Minkowski {}.new_rep(dimension).to_symbolic([])
}

/// A one-label rank-one vector `p(i)` (spenso rank1 tensor head, one argument).
fn bare_vector(atom: AtomView<'_>) -> bool {
    matches!(atom, AtomView::Fun(f) if f.get_nargs() == 1 && f.get_symbol().has_tag(&SPENSO_TAG.rank1))
}

/// Convert vakint's scalar products `vakint::dot(p(i), q(j))` of one-label
/// vectors into spenso's `spenso::dot(p(i, mink(D)), q(j, mink(D)))` (the form
/// gammaloop uses and that spenso's printers show as `p_i·q_j`). Display and
/// interchange only: vakint's reduction keeps the internal `vakint::dot` form;
/// [`from_spenso_dot_notation`] is the inverse.
pub fn to_spenso_dot_notation(expression: AtomView<'_>, dimension: impl Into<Dimension>) -> Atom {
    let mink = minkowski(dimension);
    let vakint_dot = crate::symbols::S.dot;
    expression.replace_map(|term, _context, out| {
        let AtomView::Fun(function) = term else {
            return;
        };
        if function.get_symbol() != vakint_dot || function.get_nargs() != 2 {
            return;
        }
        let arguments = function.iter().collect::<Vec<_>>();
        if !arguments.iter().all(|argument| bare_vector(*argument)) {
            return;
        }
        let mut dot = FunctionBuilder::new(SPENSO_TAG.dot);
        for argument in arguments {
            let AtomView::Fun(vector) = argument else {
                unreachable!()
            };
            dot = dot.add_arg(
                FunctionBuilder::new(vector.get_symbol())
                    .add_arg(vector.iter().next().unwrap())
                    .add_arg(&mink)
                    .finish(),
            );
        }
        **out = dot.finish();
    })
}

/// [`to_spenso_dot_notation`] with vakint's dimension symbol `vakint::d`.
pub fn to_spenso_dots(expression: AtomView<'_>) -> Atom {
    to_spenso_dot_notation(expression, crate::symbols::S.form_dimension)
}

/// Inverse of [`to_spenso_dot_notation`]: `spenso::dot(p(i, mink(D)), q(j, mink(D)))`
/// with one-label vectors in the same Minkowski space becomes `vakint::dot(p(i), q(j))`.
pub fn from_spenso_dot_notation(expression: AtomView<'_>) -> Atom {
    let mink_symbol = Minkowski::selfless_symbol();
    let compact = |atom: AtomView<'_>| -> Option<(Symbol, Atom, Atom)> {
        let AtomView::Fun(vector) = atom else {
            return None;
        };
        if vector.get_nargs() != 2 || !vector.get_symbol().has_tag(&SPENSO_TAG.rank1) {
            return None;
        }
        let mut arguments = vector.iter();
        let label = arguments.next()?;
        let representation = arguments.next()?;
        matches!(representation, AtomView::Fun(r) if r.get_symbol() == mink_symbol && r.get_nargs() == 1)
            .then(|| (vector.get_symbol(), label.to_owned(), representation.to_owned()))
    };
    expression.replace_map(|term, _context, out| {
        let AtomView::Fun(function) = term else {
            return;
        };
        if function.get_symbol() != SPENSO_TAG.dot || function.get_nargs() != 2 {
            return;
        }
        let mut arguments = function.iter();
        let (Some(left), Some(right)) = (
            compact(arguments.next().unwrap()),
            compact(arguments.next().unwrap()),
        ) else {
            return;
        };
        if left.2 != right.2 {
            return;
        }
        let vector = |(head, label, _): (Symbol, Atom, Atom)| {
            FunctionBuilder::new(head).add_arg(&label).finish()
        };
        **out = crate::symbols::S.dot(vector(left), vector(right));
    })
}
