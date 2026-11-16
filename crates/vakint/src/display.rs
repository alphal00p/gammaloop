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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;
    use std::collections::BTreeSet;
    use std::fmt::Write;

    const GOLDEN: &str = include_str!("display_forms.golden");

    fn modes() -> [(&'static str, PrintOptions); 4] {
        [
            ("plain", PrintOptions::file()),
            ("latex", PrintOptions::latex()),
            ("mathml", PrintOptions::mathml()),
            ("typst", PrintOptions::typst()),
        ]
    }

    /// The rendered samples of every vakint-registered form, one line per mode.
    fn rendered() -> String {
        crate::Vakint::initialize_vakint_symbols();
        // Create the decorated heads as vakint does (vk_symbol!) before the
        // samples are parsed (a parse would create them undecorated).
        for form in DISPLAY_TABLE.forms.iter().filter(|form| form.rust) {
            let _ = crate::utils::vakint_macros::vk_symbol!(form.name.as_str());
        }
        let mut out = String::new();
        for form in DISPLAY_TABLE.forms.iter().filter(|form| form.rust) {
            let sample = vk_parse!(form.sample.as_str()).unwrap();
            for (mode, options) in modes() {
                let mut text = String::new();
                sample
                    .as_view()
                    .format(&mut text, &options, PrintState::new())
                    .unwrap();
                writeln!(out, "{}\t{}\t{}", form.name, mode, text.replace('\n', " ")).unwrap();
            }
        }
        // Cross-implementation check of the spenso conversion (Python twin).
        let dots = vk_parse!(SPENSO_SAMPLE).unwrap();
        writeln!(
            out,
            "to_spenso_dots\t{}",
            to_spenso_dots(dots.as_view()).to_canonical_string()
        )
        .unwrap();
        out
    }

    /// Sample of the spenso-dot cross check (also used by the Python twin's test).
    const SPENSO_SAMPLE: &str =
        "dot(p(0),p(0))^2*coef(1/(d-4))+3*dot(p(1),k(2))+x*dot(vec(1),p(1))";

    /// The audit made permanent: every vakint form renders as stored, and a new
    /// symbol name in vakint's sources must be classified in display_forms.json.
    /// Regenerate the golden file with VAKINT_DISPLAY_BLESS=1.
    #[test]
    fn display_forms_render_as_stored() {
        let actual = rendered();
        if std::env::var_os("VAKINT_DISPLAY_BLESS").is_some() {
            std::fs::write(
                concat!(env!("CARGO_MANIFEST_DIR"), "/src/display_forms.golden"),
                &actual,
            )
            .unwrap();
            return;
        }
        for form in DISPLAY_TABLE.forms.iter().filter(|form| form.rust) {
            assert!(
                is_decorated(crate::utils::vakint_macros::vk_symbol!(form.name.as_str())),
                "{} has no display form",
                form.name
            );
        }
        for (expected, got) in GOLDEN.lines().zip(actual.lines()) {
            assert_eq!(
                expected, got,
                "display form changed (VAKINT_DISPLAY_BLESS=1 to accept)"
            );
        }
        assert_eq!(GOLDEN.lines().count(), actual.lines().count());
    }

    /// User rule (2026-11-16): a symbol's form never carries its own enclosing
    /// parentheses; they come only from the printer state (power base, product).
    #[test]
    fn forms_carry_no_enclosing_parentheses() {
        let opens = ["(", r"\left(", "<mo>(</mo>", "<mrow><mo>(</mo>"];
        for form in &DISPLAY_TABLE.forms {
            let mut templates = vec![&form.latex, &form.mathml, &form.typst];
            if let Some(equal) = &form.when_args_equal {
                templates.extend([&equal.latex, &equal.mathml, &equal.typst]);
            }
            for template in templates {
                assert!(
                    !opens.iter().any(|open| template.starts_with(open)),
                    "{}: the form starts with a parenthesis: {template}",
                    form.name
                );
            }
        }
        let _ = rendered();
        let doubled = [
            "((",
            r"\left(\left(",
            "<mo>(</mo><mo>(</mo>",
            "<mo>(</mo><mrow><mo>(</mo>",
        ];
        for form in DISPLAY_TABLE.forms.iter().filter(|form| form.rust) {
            let atom = vk_parse!(form.sample.as_str()).unwrap();
            // the first factor of the sample that has this form's head
            let mut single = None;
            atom.visitor(&mut |part| {
                if single.is_none()
                    && matches!(part, AtomView::Fun(f) if f.get_symbol().get_name() == form.name)
                {
                    single = Some(part.to_owned());
                }
                single.is_none()
            });
            let single = single.unwrap_or(atom);
            let x = vk_parse!("x").unwrap();
            let contexts = [
                single.clone(),
                &single * &x,
                symbolica::function!(symbolica::atom::Symbol::LOG, &single),
                single.pow(Atom::num(2)),
            ];
            for (mode, options) in modes() {
                for context in &contexts {
                    let mut text = String::new();
                    context
                        .as_view()
                        .format(&mut text, &options, PrintState::new())
                        .unwrap();
                    assert!(
                        !doubled.iter().any(|d| text.contains(d)),
                        "{} ({mode}): doubled parentheses in {text}",
                        form.name
                    );
                }
            }
        }
    }

    #[test]
    fn plain_printing_and_canonical_strings_are_unchanged() {
        let _ = rendered();
        for form in &DISPLAY_TABLE.forms {
            if !form.rust {
                continue;
            }
            let sample = vk_parse!(form.sample.as_str()).unwrap();
            // Plain output re-parses to the same expression and carries no
            // display markup.
            let plain = sample.to_canonical_string();
            assert_eq!(vk_parse!(plain.as_str()).unwrap(), sample, "{}", form.name);
            let mut file = String::new();
            sample
                .as_view()
                .format(&mut file, &PrintOptions::file(), PrintState::new())
                .unwrap();
            assert_eq!(vk_parse!(file.as_str()).unwrap(), sample, "{}", form.name);
            for marker in ['⟨', '<', '\\'] {
                assert!(!file.contains(marker), "{}: {file}", form.name);
            }
        }
    }

    #[test]
    fn every_symbol_in_the_sources_is_classified() {
        let table = &*DISPLAY_TABLE;
        let known = table
            .forms
            .iter()
            .map(|f| f.name.clone())
            .chain(table.native.iter().map(|n| n.name.clone()))
            .chain(table.internal.iter().map(|i| i.name.clone()))
            .collect::<BTreeSet<_>>();
        let patterns = [
            regex::Regex::new(r#"vk_symbol!\(\s*"([^"]+)""#).unwrap(),
            regex::Regex::new(r#"(?:^|[^a-z_])symbol!\(\s*"vakint::([^"]+)""#).unwrap(),
            regex::Regex::new(
                r#"(?:^|[^a-z_])symbol!\(\s*format!\(\s*"\{\}::([^"{}]+)"\s*,\s*crate::NAMESPACE"#,
            )
            .unwrap(),
        ];
        let mut unclassified = BTreeSet::new();
        let mut stack = vec![std::path::PathBuf::from(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src"
        ))];
        while let Some(path) = stack.pop() {
            for entry in std::fs::read_dir(&path).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    stack.push(path);
                    continue;
                }
                if path.extension().is_none_or(|e| e != "rs") {
                    continue;
                }
                let text = std::fs::read_to_string(&path).unwrap();
                for pattern in &patterns {
                    for capture in pattern.captures_iter(&text) {
                        let name = &capture[1];
                        // Wildcards (trailing `_`) are pattern variables.
                        if name.ends_with('_') {
                            continue;
                        }
                        let full = if name.contains("::") {
                            name.to_string()
                        } else {
                            format!("{}::{name}", crate::NAMESPACE)
                        };
                        if !known.contains(&full) {
                            unclassified.insert(format!("{full} ({})", path.display()));
                        }
                    }
                }
            }
        }
        assert!(
            unclassified.is_empty(),
            "new symbols without a display classification: add them to src/display_forms.json \
             (forms / native / internal): {unclassified:#?}"
        );
        assert_eq!(table.schema, "vakint.display_forms.v1");
    }

    #[test]
    fn large_arguments_keep_the_pageable_native_mathml() {
        let _ = rendered();
        let big = vk_parse!(
            (0..200)
                .map(|i| format!("{i}*d^{i}"))
                .collect::<Vec<_>>()
                .join("+")
                .as_str()
        )
        .unwrap();
        let coef = vk_parse!("coef(x)")
            .unwrap()
            .replace(vk_parse!("x").unwrap().to_pattern())
            .with(big.to_pattern());
        let mut text = String::new();
        coef.as_view()
            .format(&mut text, &PrintOptions::mathml(), PrintState::new())
            .unwrap();
        assert!(
            !text.contains("opacity"),
            "large coef must not use the bracket form"
        );
        let mut small = String::new();
        vk_parse!("coef(1/(d-4))")
            .unwrap()
            .as_view()
            .format(&mut small, &PrintOptions::mathml(), PrintState::new())
            .unwrap();
        assert!(small.contains("opacity"), "{small}");
    }

    /// The display forms must not change the order in which symbols are created
    /// (Symbolica orders factors by it): the symbol registry after vakint's
    /// vocabulary, a Vakint instance and a typical parse is identical with and
    /// without the forms (VAKINT_NO_DISPLAY_FORMS), name by name and in order.
    #[test]
    fn display_forms_keep_the_symbol_creation_order() {
        const CHILD: &str = "VAKINT_DISPLAY_ORDER_CHILD";
        if std::env::var_os(CHILD).is_some() {
            crate::Vakint::initialize_vakint_symbols();
            let _vakint = crate::Vakint::new().unwrap();
            let x = vk_parse!(
                "PGamma(2)*PGamma(-2)^2/(dot(p(1),p(1))*PMaster(no1))+coef(1/(d-4))*Master(1)+ep^2*Oep(2,PMaster(t105))"
            )
            .unwrap();
            println!("PLAIN {}", x.to_canonical_string());
            println!("DECORATED {}", is_decorated(crate::symbols::S.p_gamma));
            for (_, name) in symbolica::state::State::symbol_iter() {
                println!("SYMBOL {name}");
            }
            return;
        }
        let run = |disable: bool| {
            let mut command = std::process::Command::new(std::env::current_exe().unwrap());
            command
                .args([
                    "--exact",
                    "display::tests::display_forms_keep_the_symbol_creation_order",
                    "--nocapture",
                    "--test-threads",
                    "1",
                ])
                .env(CHILD, "1");
            if disable {
                command.env("VAKINT_NO_DISPLAY_FORMS", "1");
            } else {
                command.env_remove("VAKINT_NO_DISPLAY_FORMS");
            }
            let output = command.output().unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            String::from_utf8_lossy(&output.stdout)
                .lines()
                .filter(|l| {
                    l.starts_with("SYMBOL ")
                        || l.starts_with("PLAIN ")
                        || l.starts_with("DECORATED ")
                })
                .map(str::to_owned)
                .collect::<Vec<_>>()
        };
        let with_forms = run(false);
        let without = run(true);
        assert!(
            with_forms.contains(&"DECORATED true".to_owned()),
            "forms not active"
        );
        assert!(
            without.contains(&"DECORATED false".to_owned()),
            "switch not active"
        );
        let strip = |lines: &[String]| {
            lines
                .iter()
                .filter(|l| !l.starts_with("DECORATED "))
                .cloned()
                .collect::<Vec<_>>()
        };
        assert!(with_forms.len() > 100, "{}", with_forms.len());
        assert_eq!(strip(&with_forms), strip(&without));
    }

    #[test]
    fn spenso_dot_round_trip() {
        crate::Vakint::initialize_vakint_symbols();
        let source = vk_parse!(
            "dot(p(0),p(0))^2*coef(1/(d-4))+3*dot(p(1),k(2))+dot(p(1),p(2)+p(3))+x*dot(vec(1),p(1))"
        )
        .unwrap();
        let spenso = to_spenso_dots(source.as_view());
        assert_ne!(spenso, source);
        assert!(spenso.contains_symbol(SPENSO_TAG.dot));
        assert_eq!(from_spenso_dot_notation(spenso.as_view()), source);
        let text = spenso.to_canonical_string();
        assert!(text.contains("spenso::"), "{text}");
        assert!(text.contains("mink(vakint::{}::d)"), "{text}");
        // A gammaloop-style input with a symbolic dimension converts back too.
        let gl = to_spenso_dot_notation(
            vk_parse!("dot(p(1),k(1))").unwrap().as_view(),
            symbolica::symbol!("gammalooprs::dim"),
        );
        assert_eq!(
            from_spenso_dot_notation(gl.as_view()),
            vk_parse!("dot(p(1),k(1))").unwrap()
        );
    }
}
