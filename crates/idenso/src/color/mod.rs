#[cfg(test)]
use std::sync::LazyLock;

use spenso::{
    network::{library::symbolic::ExplicitKey, tags::SPENSO_TAG as T},
    shadowing::{Collectable, IntoAtom, TensorCollectExt, symbolica_utils::SpensoPrintSettings},
    structure::{
        Canonicalized, TensorStructure,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        dimension::Dimension,
        representation::RepName,
        slot::AbsInd,
    },
    tensor_symbol,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, EvaluationInfo, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
    domains::rational::Rational,
    function,
    printer::{PrintOptions, PrintState},
    symbol,
};

use crate::{
    color::{casimir::CofDimensionInvariantRewriter, simplify::ColorAlgebraSimplifier},
    representations::{ColorAntiFundamental, ColorSextet},
    selective_expand::SelectiveExpand,
    shorthands::metric::PermuteWithMetric,
};

use super::rep_symbols::RS;
use super::representations::{ColorAdjoint, ColorFundamental};

mod casimir;
mod conjugate;
mod macros;
mod simplify;

pub use conjugate::color_conj_impl;

#[derive(Debug)]
pub enum ColorError {
    NotFully(Atom),
}

pub struct ColorSymbols {
    pub nc_: Symbol,
    pub adj_: Symbol,
    /// Symbol backing the color fundamental representation function.
    pub fundamental_rep: Symbol,
    /// Symbol backing the color adjoint representation function.
    pub adjoint_rep: Symbol,
    /// The generator symbol
    pub t: Symbol,
    /// The structure constant symbol i.e. [T^a, T^b] = i f^{abc} T^c
    pub f: Symbol,
    /// The symmetric color invariant symbol.
    pub d: Symbol,
    /// The degree-k Gram invariant symbol for two symmetric traces.
    pub gram: Symbol,
    /// The degree-k Casimir eigenvalue symbol.
    pub cas: Symbol,
    /// The degree-k Dynkin index symbol.
    pub idx: Symbol,
    /// The number of colors symbol (i.e. the dimension of the fundamental representation) usually Nc=3
    pub nc: Symbol,
}

impl ColorSymbols {
    pub fn chain_t<'a>(&self, adjoint_index: impl Into<AtomOrView<'a>>) -> Atom {
        FunctionBuilder::new(self.t)
            .add_arg(adjoint_index)
            .add_arg(Atom::var(T.chain_in))
            .add_arg(Atom::var(T.chain_out))
            .finish()
    }

    pub fn explicit_t<'a, 'b, 'c>(
        &self,
        adjoint_index: impl Into<AtomOrView<'a>>,
        left_fundamental: impl Into<AtomOrView<'b>>,
        right_fundamental: impl Into<AtomOrView<'c>>,
    ) -> Atom {
        FunctionBuilder::new(self.t)
            .add_arg(adjoint_index)
            .add_arg(left_fundamental)
            .add_arg(right_fundamental)
            .finish()
    }

    pub fn structure_f<'a, 'b, 'c>(
        &self,
        a: impl Into<AtomOrView<'a>>,
        b: impl Into<AtomOrView<'b>>,
        c: impl Into<AtomOrView<'c>>,
    ) -> Atom {
        FunctionBuilder::new(self.f)
            .add_arg(a)
            .add_arg(b)
            .add_arg(c)
            .finish()
    }

    pub fn symmetric_d<'a>(&self, rep: impl Into<AtomOrView<'a>>, args: Vec<Atom>) -> Atom {
        args.into_iter()
            .fold(FunctionBuilder::new(self.d).add_arg(rep), |d, arg| {
                d.add_arg(arg)
            })
            .finish()
    }

    pub fn gram<D: IntoAtom, L: IntoAtom, R: IntoAtom>(
        &self,
        degree: D,
        left: L,
        right: R,
    ) -> Atom {
        FunctionBuilder::new(self.gram)
            .add_arg(degree.into_atom())
            .add_arg(left.into_atom())
            .add_arg(right.into_atom())
            .finish()
    }

    pub fn cas<D: IntoAtom, R: IntoAtom>(&self, degree: D, rep: R) -> Atom {
        FunctionBuilder::new(self.cas)
            .add_arg(degree.into_atom())
            .add_arg(rep.into_atom())
            .finish()
    }

    pub fn idx<D: IntoAtom, R: IntoAtom>(&self, degree: D, rep: R) -> Atom {
        FunctionBuilder::new(self.idx)
            .add_arg(degree.into_atom())
            .add_arg(rep.into_atom())
            .finish()
    }

    pub fn symmetric_generator_trace<R: IntoAtom, A: IntoAtom>(
        &self,
        rep: R,
        adjoint_indices: impl IntoIterator<Item = A>,
    ) -> Atom {
        let factors = adjoint_indices
            .into_iter()
            .map(|adjoint| self.chain_t(adjoint.into_atom()))
            .collect::<Vec<_>>();
        spenso::shadowing::trace_sym(rep.into_atom(), factors)
    }

    #[cfg(test)]
    pub(crate) fn initialize_tensor_symbols(&self) {
        let _ = self.t;
        let _ = self.f;
        let _ = self.d;

        let _ = self.gram;
        let _ = self.cas;
        let _ = self.idx;
    }

    // Generator for the adjoint representation of SU(N)
    pub fn t_strct<Aind: AbsInd>(
        &self,
        fundimd: impl Into<Dimension>,
        adim: impl Into<Dimension>,
    ) -> Canonicalized<ExplicitKey<Aind>> {
        let nc = fundimd.into();
        ExplicitKey::from_iter(
            [
                ColorAdjoint {}.new_rep(adim).cast(),
                ColorFundamental {}.new_rep(nc).to_lib(),
                ColorAntiFundamental {}.new_rep(nc).cast(),
            ],
            self.t,
            None,
        )
    }
    pub fn t_pattern(
        &self,
        fundimd: impl Into<Dimension>,
        adim: impl Into<Dimension>,
        a: impl Into<AbstractIndex>,
        i: impl Into<AbstractIndex>,
        j: impl Into<AbstractIndex>,
    ) -> Atom {
        let structure = self.t_strct::<AbstractIndex>(fundimd, adim);
        let logical_indices = [a.into(), i.into(), j.into()];
        let storage_indices = structure.layout().logical_to_canonical(&logical_indices);
        structure
            .into_canonical()
            .reindex_storage(&storage_indices)
            .unwrap()
            .permute_with_metric()
    }

    pub fn f_strct<Aind: AbsInd>(
        &self,
        adim: impl Into<Dimension>,
    ) -> Canonicalized<ExplicitKey<Aind>> {
        let adim = adim.into();
        ExplicitKey::from_iter(
            [
                ColorAdjoint {}.new_rep(adim),
                ColorAdjoint {}.new_rep(adim),
                ColorAdjoint {}.new_rep(adim),
            ],
            self.f,
            None,
        )
    }

    pub fn f_pattern(
        &self,
        adim: impl Into<Dimension>,
        a: impl Into<AbstractIndex>,
        b: impl Into<AbstractIndex>,
        c: impl Into<AbstractIndex>,
    ) -> Atom {
        let structure = self.f_strct::<AbstractIndex>(adim);
        let logical_indices = [a.into(), b.into(), c.into()];
        let storage_indices = structure.layout().logical_to_canonical(&logical_indices);
        structure
            .into_canonical()
            .reindex_storage(&storage_indices)
            .unwrap()
            .permute_with_metric()
    }
}

#[derive(Clone, Copy)]
enum ColorInvariantPrintKind {
    Gram,
    Casimir,
    Index,
}

impl ColorInvariantPrintKind {
    fn print(self, atom: AtomView<'_>, opt: &PrintOptions, state: &PrintState) -> Option<String> {
        if opt.mode.is_mathml() {
            return self.mathml(atom, opt, state);
        }
        let resolved = SpensoPrintSettings::resolve(opt)?;
        let settings = resolved.presentation;
        let AtomView::Fun(f) = atom else { return None };
        let args = f.iter().collect::<Vec<_>>();
        let arity = if matches!(self, Self::Gram) { 3 } else { 2 };
        if args.len() != arity {
            return None;
        }
        if !settings.explicit_invariants && !settings.with_dim && small_integer(args[0]) == Some(2)
        {
            let alias = match self {
                Self::Casimir if is_color_rep(args[1], "spenso::cof") => Some(("C_F", "CF")),
                Self::Casimir if is_color_rep(args[1], "spenso::coad") => Some(("C_A", "CA")),
                Self::Index if is_color_rep(args[1], "spenso::cof") => Some(("T_R", "TR")),
                _ => None,
            };
            if let Some((script, plain)) = alias {
                let label = if settings.symbol_scripts {
                    script.to_owned()
                } else {
                    Self::function_head(plain, opt)
                };
                return Some(Self::colorize(&label, opt));
            }
        }
        let (symbol, function) = match self {
            Self::Gram => ("G", "gram"),
            Self::Casimir => ("C", "cas"),
            Self::Index => ("I", "idx"),
        };
        let mut degree = String::new();
        args[0].format(&mut degree, opt, PrintState::new()).unwrap();
        let head = if settings.symbol_scripts {
            let (open, close) = resolved.script_delimiters();
            // Group symbolic or multi-digit degrees so following powers apply to
            // the complete invariant, never to part of its degree label.
            if degree.chars().count() == 1 {
                format!("{symbol}_{degree}")
            } else {
                format!("{symbol}_{open}{degree}{close}")
            }
        } else {
            Self::function_head(function, opt)
        };
        let mut arguments = args[1..]
            .iter()
            .map(|arg| Self::representation(*arg, opt))
            .collect::<Vec<_>>();
        if !settings.symbol_scripts {
            arguments.insert(0, degree);
        }
        // These are scalar function arguments, not index lists: always delimit
        // them, independently of the tensor parentheses/comma preferences.
        Some(format!(
            "{}({})",
            Self::colorize(&head, opt),
            arguments.join(",")
        ))
    }

    fn mathml_settings(opt: &PrintOptions) -> SpensoPrintSettings {
        match opt.custom_print_mode.get("spenso") {
            Some(PrintUserData::Integer(encoded)) => SpensoPrintSettings::from(*encoded as usize),
            _ => SpensoPrintSettings::typst(),
        }
    }

    fn mathml_fits(atom: AtomView<'_>, depth: usize, remaining: &mut usize) -> bool {
        if depth > 12
            || *remaining == 0
            || atom.get_byte_size() > 2048
            || atom
                .get_symbol()
                .is_some_and(|symbol| symbol.get_name().len() > 1024)
        {
            return false;
        }
        *remaining -= 1;
        atom.children()
            .all(|child| Self::mathml_fits(child, depth + 1, remaining))
    }

fn invariant_print_special(
    kind: &ColorInvariantPrintKind,
    args: &[AtomView<'_>],
    symbol_scripts: bool,
) -> Option<&'static str> {
    if small_integer(args[0])? != 2 {
        return None;
    }

    match kind {
        ColorInvariantPrintKind::Casimir if is_color_rep(args[1], "cof") => {
            Some(if symbol_scripts { "C_F" } else { "CF" })
        }
        ColorInvariantPrintKind::Casimir if is_color_rep(args[1], "coad") => {
            Some(if symbol_scripts { "C_A" } else { "CA" })
        }
        ColorInvariantPrintKind::Index if is_color_rep(args[1], "cof") => {
            Some(if symbol_scripts { "T_R" } else { "TR" })
        }
        _ => None,
    }
}

fn print_invariant_args(
    out: &mut String,
    args: &[AtomView<'_>],
    opt: &PrintOptions,
    parens: bool,
    commas: bool,
) {
    if parens {
        out.push('(');
    } else if !args.is_empty() {
        out.push(' ');
    }

    for (position, arg) in args.iter().enumerate() {
        if position > 0 {
            if commas {
                out.push(',');
            } else {
                out.push(' ');
            }
        }
        arg.format(out, opt, PrintState::new()).unwrap();
    }

    if parens {
        out.push(')');
    }
}

    fn colorize(head: &str, opt: &PrintOptions) -> String {
        if opt.color_builtin_symbols {
            nu_ansi_term::Color::Magenta.paint(head).to_string()
        } else {
            head.to_owned()
        }
    }

    fn function_head(name: &str, opt: &PrintOptions) -> String {
        match SpensoPrintSettings::resolve(opt).unwrap().backend {
            SpensoPrintBackend::Typst => format!("op({name:?})"),
            SpensoPrintBackend::Latex => format!("\\operatorname{{{name}}}"),
            SpensoPrintBackend::Plain => name.to_owned(),
        }
    }

    fn color_count(atom: AtomView<'_>, opt: &PrintOptions, state: &PrintState) -> Option<String> {
        if !matches!(atom, AtomView::Var(_)) {
            return None;
        }
        if opt.mode.is_mathml() {
            let source = if Self::mathml_settings(opt).symbol_scripts {
                "<msub><mi>N</mi><mi>c</mi></msub>"
            } else {
                "<mi mathvariant=\"normal\">Nc</mi>"
            };
            return Some(if state.in_exp_base {
                Self::mathml_group(source)
            } else {
                source.to_owned()
            });
        }
        let settings = SpensoPrintSettings::resolve(opt)?.presentation;
        let out = if settings.symbol_scripts {
            "N_c".to_owned()
        } else {
            Self::function_head("Nc", opt)
        };
        Some(Self::colorize(&out, opt))
    }
}

fn small_integer(expr: AtomView<'_>) -> Option<i64> {
    let AtomView::Num(number) = expr else {
        return None;
    };
    let CoefficientView::Natural(value, 1, 0, 1) = number.get_coeff_view() else {
        return None;
    };
    Some(value)
}

fn is_color_rep(expr: AtomView<'_>, name: &str) -> bool {
    matches!(
        expr,
        AtomView::Fun(rep) if rep.get_symbol().get_stripped_name() == name && rep.get_nargs() == 1
    )
}

spenso::symbolica_init_lazy_static! {
pub static CS, CS_INNER: ColorSymbols = || {
    fn representation_symbol(rep: Atom) -> Symbol {
        let AtomView::Fun(f) = rep.as_view() else {
            unreachable!("Color representations are symbolic functions")
        };
        f.get_symbol()
    }

    ColorSymbols {
        t: tensor_symbol!("spenso::t"; Real; print = |a, opt, _state| {

            match SpensoPrintSettings::resolve(opt) {
                Some(resolved)=>{
                    let (script_open, script_close) = resolved.script_delimiters();
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = resolved.presentation;


                    let AtomView::Fun(f)=a else {
                        return None;
                    };
                    if f.get_nargs()!=3 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let a = argitem.next().unwrap();
                    let b = argitem.next().unwrap();
                    let mut c = argitem.next().unwrap();

                    let mut out = "t".to_string();
                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens{
                        out.push(script_open);
                    }
                    a.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas{
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    b.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens{
                        out.push(script_close);
                    }
                    if symbol_scripts{
                        out.push('_');
                    }


                    if parens{
                        out.push(script_open);
                    }else if !symbol_scripts{
                        out.push(' ');
                    }

                    let AtomView::Fun(f)=c else {
                        return None;
                    };
                    if f.get_nargs()!=1 {
                        return None;
                    }
                    if f.get_symbol()!=AIND_SYMBOLS.dind{
                        return None;
                    }
                    c = f.iter().next().unwrap();
                    c.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens{
                        out.push(script_close);
                    }
                    Some(out)
                }
                _=>None}

        }),
        f: tensor_symbol!("spenso::f"; Real, Antisymmetric; print = |a, opt, _state| {

            match SpensoPrintSettings::resolve(opt) {
                Some(resolved)=>{
                    let (script_open, script_close) = resolved.script_delimiters();
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = resolved.presentation;


                    let AtomView::Fun(f)=a else {
                        return None;
                    };
                    if f.get_nargs()!=3 {
                        return None;
                    }
                    let mut argitem = f.iter();
                    let a = argitem.next().unwrap();
                    let b = argitem.next().unwrap();
                    let c = argitem.next().unwrap();

                    let mut out = "f".to_string();
                    if symbol_scripts {
                        out.push('^');
                    }
                    if opt.color_builtin_symbols {
                        out = nu_ansi_term::Color::Magenta.paint(out).to_string();
                    }

                    if parens{
                        out.push(script_open);
                    }
                    a.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas{
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    b.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas{
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    c.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens{
                        out.push(script_close);
                    }
                    Some(out)
                }
                _=>None}

        }),
        ca: symbol!("spenso::CA";Real;eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(3,1).to_multi_prec_float(prec).into()))),
        cf: symbol!("spenso::CF";Real; eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(4,3).to_multi_prec_float(prec).into()))),
        d: tensor_symbol!("spenso::d"),
        gram: symbol!("spenso::gram"; Real; print = |a, opt, _state| {
            print_color_invariant(a, opt, ColorInvariantPrintKind::Gram)
        }),
        cas: symbol!("spenso::cas"; Real; print = |a, opt, _state| {
            print_color_invariant(a, opt, ColorInvariantPrintKind::Casimir)
        }),
        idx: symbol!("spenso::idx"; Real; print = |a, opt, _state| {
            print_color_invariant(a, opt, ColorInvariantPrintKind::Index)
        }),
        trace_dummy: symbol!("x"),
        fundamental_rep: representation_symbol(
            ColorFundamental {}.to_symbolic(std::iter::empty::<Atom>()),
        ),
        adjoint_rep: representation_symbol(ColorAdjoint {}.to_symbolic(std::iter::empty::<Atom>())),
        adj_: symbol!("adj_"),
        nc_: symbol!("nc_"),
        tr: symbol!("spenso::TR";Real;eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(1,2).to_multi_prec_float(prec).into()))),
        nc: symbol!("spenso::Nc";Real;eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(3,1).to_multi_prec_float(prec).into()))),
    }
};
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColorSimplifySettings {
    /// Normalize non-color chains and metric contractions outside collected
    /// color payloads. Disable this to preserve their factorized coefficients.
    pub simplify_non_color: bool,
    /// Whether closed color chains should be evaluated as traces.
    pub evaluate_traces: bool,
    /// Whether contractions between generators on different open chains should
    /// be expanded with the fundamental Fierz identity.
    pub expand_cross_chain_fierz: bool,
    /// Whether invariant factors for `cof(N)` should be written directly in
    /// terms of the fundamental dimension.
    pub substitute_cof_dimension_invariants: bool,
}

impl Default for ColorSimplifySettings {
    fn default() -> Self {
        Self {
            simplify_non_color: true,
            evaluate_traces: true,
            expand_cross_chain_fierz: true,
            substitute_cof_dimension_invariants: false,
        }
    }
}

impl ColorSimplifySettings {
    /// Leaves collected `trace(...)` nodes unevaluated; Fierz contractions can still join them.
    pub fn without_trace_evaluation(mut self) -> Self {
        self.evaluate_traces = false;
        self
    }

    /// Keeps separate open chains and traces instead of applying cross-chain Fierz
    /// expansion.
    pub fn without_cross_chain_fierz_expansion(mut self) -> Self {
        self.expand_cross_chain_fierz = false;
        self
    }

    /// Rewrites supported `cof(N)` invariants to explicit dimension formulas.
    pub fn with_cof_dimension_invariants(mut self) -> Self {
        self.substitute_cof_dimension_invariants = true;
        self
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ColorCasimirSettings {
    /// Rewrite the explicit fundamental dimension using `d_F = C_A`.
    pub rewrite_fundamental_dimension: bool,
    /// Substitute the common fundamental normalization `T_F = 1/2`.
    pub substitute_fundamental_index: bool,
}

impl Default for ColorCasimirSettings {
    fn default() -> Self {
        Self {
            rewrite_fundamental_dimension: true,
            substitute_fundamental_index: false,
        }
    }
}

impl ColorCasimirSettings {
    /// Keep `d_F` explicit instead of applying the SU(N) relation `d_F = C_A`.
    pub fn without_fundamental_dimension_rewrite(mut self) -> Self {
        self.rewrite_fundamental_dimension = false;
        self
    }

    /// Apply the standard fundamental index normalization `T_F = 1/2`.
    pub fn with_fundamental_index_normalization(mut self) -> Self {
        self.substitute_fundamental_index = true;
        self
    }
}

/// Rewrite color dimensions and scalar invariants between SU(N) conventions.
/// Indexed generator, structure-constant and trace reduction belongs to
/// `simplify_algebra` method on [`SymbolicTensor`](crate::tensor::SymbolicTensor).
pub trait ColorSimplifier {
    /// Replace concrete QCD fundamental and adjoint dimensions by their
    /// parametric SU(Nc) expressions while preserving all color indices.
    ///
    /// This is useful for exact color-algebra comparisons: public expressions
    /// may use the physical `cof(3)` and `coad(8)` representations, while an
    /// intermediate symbolic calculation should retain its `Nc` dependence.
    fn to_parametric_color(&self) -> Atom;

    /// Rewrites the explicit representation dimensions into a Casimir basis.
    fn to_color_casimir(&self, fundamental_rep: AtomView<'_>, adjoint_rep: AtomView<'_>) -> Atom;

    /// Rewrites explicit dimensions with control over SU(N) normalization choices.
    fn to_color_casimir_with(
        &self,
        fundamental_rep: AtomView<'_>,
        adjoint_rep: AtomView<'_>,
        settings: ColorCasimirSettings,
    ) -> Atom;

    /// Rewrites supported `cof(N)` invariant factors into explicit dimension formulas.
    fn to_cof_dimension_invariants(&self) -> Atom;
}
impl ColorSimplifier for Atom {
    fn to_parametric_color(&self) -> Atom {
        self.as_view().to_parametric_color()
    }

    fn to_color_casimir(&self, fundamental_rep: AtomView<'_>, adjoint_rep: AtomView<'_>) -> Atom {
        self.to_color_casimir_with(
            fundamental_rep,
            adjoint_rep,
            ColorCasimirSettings::default(),
        )
    }

    fn to_color_casimir_with(
        &self,
        fundamental_rep: AtomView<'_>,
        adjoint_rep: AtomView<'_>,
        settings: ColorCasimirSettings,
    ) -> Atom {
        casimir::color_casimir_basis_impl(
            self.as_atom_view(),
            fundamental_rep,
            adjoint_rep,
            settings,
        )
    }

    fn to_cof_dimension_invariants(&self) -> Atom {
        CofDimensionInvariantRewriter.run(self.as_atom_view())
    }
}

impl ColorSimplifier for AtomView<'_> {
    fn to_parametric_color(&self) -> Atom {
        let adjoint = ColorAdjoint {};
        let fundamental = ColorFundamental {};
        let nc = Atom::var(CS.nc);
        self.replace(adjoint.to_symbolic([RS.d_, RS.a_]))
            .with(
                adjoint.to_symbolic([nc.clone().pow(Atom::num(2)) - Atom::one(), Atom::var(RS.a_)]),
            )
            .replace(fundamental.to_symbolic([RS.d_, RS.a_]))
            .with(fundamental.to_symbolic([nc, Atom::var(RS.a_)]))
    }

    fn to_color_casimir(&self, fundamental_rep: AtomView<'_>, adjoint_rep: AtomView<'_>) -> Atom {
        self.to_color_casimir_with(
            fundamental_rep,
            adjoint_rep,
            ColorCasimirSettings::default(),
        )
    }

    fn to_color_casimir_with(
        &self,
        fundamental_rep: AtomView<'_>,
        adjoint_rep: AtomView<'_>,
        settings: ColorCasimirSettings,
    ) -> Atom {
        casimir::color_casimir_basis_impl(
            self.as_atom_view(),
            fundamental_rep,
            adjoint_rep,
            settings,
        )
    }

    fn to_cof_dimension_invariants(&self) -> Atom {
        CofDimensionInvariantRewriter.run(self.as_atom_view())
    }

    fn expand_color(&self) -> Vec<(Atom, Atom)> {
        let cof = ColorFundamental {};
        let coaf = ColorFundamental {}.dual();
        let coad = ColorAdjoint {};

        let color_trace_pat = function!(T.trace, cof.to_symbolic([RS.b__]), RS.a___).to_pattern();
        let color_chain_pat = function!(
            T.chain,
            cof.to_symbolic([RS.b__]),
            coaf.to_symbolic([RS.c__]),
            RS.a___
        )
        .to_pattern();
        let color_d_pat = function!(CS.d, RS.a___).to_pattern();
        let color_gram_pat = function!(CS.gram, RS.a___).to_pattern();
        let color_cas_pat = function!(CS.cas, RS.a___).to_pattern();
        let color_idx_pat = function!(CS.idx, RS.a___).to_pattern();
        let cof_pat = function!(RS.f_, RS.a___, cof.to_symbolic([RS.b__]), RS.c___).to_pattern();
        let coaf_pat = function!(RS.f_, RS.a___, coaf.to_symbolic([RS.b__]), RS.c___).to_pattern();
        let coad_pat = function!(RS.f_, RS.a___, coad.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[
            color_trace_pat,
            color_chain_pat,
            color_d_pat,
            color_gram_pat,
            color_cas_pat,
            color_idx_pat,
            cof_pat,
            coad_pat,
            coaf_pat,
        ])
    }

    // fn canonize_color(&self) -> Atom {
    //     self..canonize_color()
    // }

    fn wrap_color(&self, symbol: Symbol) -> Atom {
        self.expand_color()
            .into_iter()
            .fold(Atom::Zero, |a, (c, s)| a + function!(symbol, c) * s)
    }
}

#[cfg(test)]
mod test;
