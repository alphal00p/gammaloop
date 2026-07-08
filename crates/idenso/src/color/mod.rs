#[cfg(test)]
use std::sync::LazyLock;

use spenso::{
    network::{library::symbolic::ExplicitKey, tags::SPENSO_TAG as T},
    shadowing::{Collectable, IntoAtom, TensorCollectExt, symbolica_utils::SpensoPrintSettings},
    structure::{
        TensorStructure,
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
    printer::{PrintOptions, PrintState, PrintUserData},
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
    /// The adjoint representation dimension symbol, i.e. NA = Nc^2 - 1.
    pub na: Symbol,
    /// The adjoint Casimir symbol, i.e. CA = Nc
    pub ca: Symbol,
    /// The fundamental Casimir symbol. T^a_ij T^a_jk = CF delta_ik -> CF = TR (na/nc)
    pub cf: Symbol,
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
    /// Dummy adjoint-index symbol used in color trace decompositions.
    pub trace_dummy: Symbol,
    /// The trace constant symbol i.e. Tr(T^a T^b) = TR delta^{ab}. Usually TR=1/2
    pub tr: Symbol,
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
    ) -> ExplicitKey<Aind> {
        let nc = fundimd.into();
        let res = ExplicitKey::from_iter(
            [
                ColorAdjoint {}.new_rep(adim).cast(),
                ColorFundamental {}.new_rep(nc).to_lib(),
                ColorAntiFundamental {}.new_rep(nc).cast(),
            ],
            self.t,
            None,
        );
        debug_assert!(res.rep_permutation.is_identity());
        res.structure
    }
    pub fn t_pattern(
        &self,
        fundimd: impl Into<Dimension>,
        adim: impl Into<Dimension>,
        a: impl Into<AbstractIndex>,
        i: impl Into<AbstractIndex>,
        j: impl Into<AbstractIndex>,
    ) -> Atom {
        self.t_strct(fundimd, adim)
            .reindex(&[a.into(), i.into(), j.into()])
            .unwrap()
            .permute_with_metric()
    }

    pub fn f_strct<Aind: AbsInd>(&self, adim: impl Into<Dimension>) -> ExplicitKey<Aind> {
        let adim = adim.into();
        let res = ExplicitKey::from_iter(
            [
                ColorAdjoint {}.new_rep(adim),
                ColorAdjoint {}.new_rep(adim),
                ColorAdjoint {}.new_rep(adim),
            ],
            self.f,
            None,
        );
        debug_assert!(res.rep_permutation.is_identity());
        res.structure
    }

    pub fn f_pattern(
        &self,
        adim: impl Into<Dimension>,
        a: impl Into<AbstractIndex>,
        b: impl Into<AbstractIndex>,
        c: impl Into<AbstractIndex>,
    ) -> Atom {
        self.f_strct(adim)
            .reindex(&[a.into(), b.into(), c.into()])
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

fn print_color_invariant(
    atom: AtomView<'_>,
    opt: &PrintOptions,
    kind: ColorInvariantPrintKind,
) -> Option<String> {
    match opt.custom_print_mode.get("spenso") {
        Some(PrintUserData::Integer(i)) => {
            let SpensoPrintSettings {
                parens,
                symbol_scripts,
                commas,
                ..
            } = SpensoPrintSettings::from(*i as usize);

            let AtomView::Fun(f) = atom else {
                return None;
            };
            let args = f.iter().collect::<Vec<_>>();
            let expected_nargs = match kind {
                ColorInvariantPrintKind::Gram => 3,
                ColorInvariantPrintKind::Casimir | ColorInvariantPrintKind::Index => 2,
            };
            if args.len() != expected_nargs {
                return None;
            }

            if let Some(special) = invariant_print_special(&kind, &args, symbol_scripts) {
                return Some(colorize_invariant_head(special, opt));
            }

            let head = if symbol_scripts {
                let mut rank = String::new();
                args[0].format(&mut rank, opt, PrintState::new()).unwrap();
                match kind {
                    ColorInvariantPrintKind::Gram => format!("G_{rank}"),
                    ColorInvariantPrintKind::Casimir => format!("C_{rank}"),
                    ColorInvariantPrintKind::Index => format!("I_{rank}"),
                }
            } else {
                match kind {
                    ColorInvariantPrintKind::Gram => "gram".to_string(),
                    ColorInvariantPrintKind::Casimir => "cas".to_string(),
                    ColorInvariantPrintKind::Index => "idx".to_string(),
                }
            };

            let mut out = colorize_invariant_head(&head, opt);
            let printed_args = if symbol_scripts {
                &args[1..]
            } else {
                &args[..]
            };
            print_invariant_args(&mut out, printed_args, opt, parens, commas);
            Some(out)
        }
        _ => None,
    }
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

fn colorize_invariant_head(head: &str, opt: &PrintOptions) -> String {
    if opt.color_builtin_symbols {
        nu_ansi_term::Color::Magenta.paint(head).to_string()
    } else {
        head.to_string()
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

            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(i))=>{
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = SpensoPrintSettings::from(*i as usize);


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
                        out.push('(');
                    }
                    a.format(&mut out, opt, PrintState::new()).unwrap();
                    if commas{
                        out.push(',');
                    } else {
                        out.push(' ');
                    }
                    b.format(&mut out, opt, PrintState::new()).unwrap();
                    if parens{
                        out.push(')');
                    }
                    if symbol_scripts{
                        out.push('_');
                    }


                    if parens{
                        out.push('(');
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
                        out.push(')');
                    }
                    Some(out)
                }
                _=>None}

        }),
        f: tensor_symbol!("spenso::f"; Real, Antisymmetric; print = |a, opt, _state| {

            match opt.custom_print_mode.get("spenso") {
                Some(PrintUserData::Integer(i))=>{
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = SpensoPrintSettings::from(*i as usize);


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
                        out.push('(');
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
                        out.push(')');
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
        na: symbol!("NA";Real;eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(3,1).to_multi_prec_float(prec).into()))),
        tr: symbol!("spenso::TR";Real;eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(1,2).to_multi_prec_float(prec).into()))),
        nc: symbol!("spenso::Nc";Real;eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(3,1).to_multi_prec_float(prec).into()))),
    }
};
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColorSimplifySettings {
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
            evaluate_traces: true,
            expand_cross_chain_fierz: true,
            substitute_cof_dimension_invariants: false,
        }
    }
}

impl ColorSimplifySettings {
    /// Leaves collected `trace(...)` nodes inert after chain collection.
    pub fn without_trace_evaluation(mut self) -> Self {
        self.evaluate_traces = false;
        self
    }

    /// Keeps separate open chains instead of applying cross-chain Fierz
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
    /// Rewrite the fundamental dimension `Nc` using `CA = Nc`.
    pub rewrite_nc: bool,
    /// Rewrite the adjoint dimension `NA` using `NA = Nc^2 - 1 = 2 CA CF`.
    pub rewrite_na: bool,
    /// Substitute the common fundamental normalization `TR = 1/2`.
    pub substitute_tr: bool,
}

impl Default for ColorCasimirSettings {
    fn default() -> Self {
        Self {
            rewrite_nc: true,
            rewrite_na: true,
            substitute_tr: false,
        }
    }
}

impl ColorCasimirSettings {
    /// Also apply the standard fundamental trace normalization `TR = 1/2`.
    pub fn with_trace_normalization(mut self) -> Self {
        self.substitute_tr = true;
        self
    }
}

/// Trait for applying SU(N) color algebra simplification rules to a symbolic expression.
///
/// Implementors provide a method to simplify expressions containing color factors
/// like structure constants (`f_abc`), generators (`T^a`), traces (`TR`), and the
/// number of colors (`Nc`).
pub trait ColorSimplifier {
    /// Attempts to simplify the color structure of the expression.
    ///
    /// Applies various identities of SU(N) algebra, including Fierz identities,
    /// Casimir relations, and contractions involving `f_abc` and `T^a`.
    ///
    /// # Returns
    /// - `Ok(Atom)`: The simplified expression, ideally with only color-scalar factors remaining.
    /// - `Err(ColorError::NotFully(Atom))`: If the simplification could not fully remove all
    ///   explicit color index structures (like `cof(...)`, `coad(...)`). The partially
    ///   simplified `Atom` is included in the error.
    fn simplify_color(&self) -> Atom;

    /// Simplifies color structures with explicit chain/trace settings.
    fn simplify_color_with(&self, settings: ColorSimplifySettings) -> Atom;

    /// Rewrites `Nc`/`NA` scalar factors into the `CA`, `CF` Casimir basis.
    fn to_color_casimir(&self) -> Atom;

    /// Rewrites color scalar factors with explicit control over normalization choices.
    fn to_color_casimir_with(&self, settings: ColorCasimirSettings) -> Atom;

    /// Rewrites supported `cof(N)` invariant factors into explicit dimension formulas.
    fn to_cof_dimension_invariants(&self) -> Atom;

    /// Expands factorized terms around color representation factors.
    fn expand_color(&self) -> Vec<(Atom, Atom)>;

    /// Collects factorized terms around color representation factors.
    fn collect_color(&self) -> Atom;

    fn collect_color_constants(&self) -> Atom;

    // fn canonize_color(&self) -> Atom;

    fn wrap_color(&self, symbol: Symbol) -> Atom;
}
impl ColorSimplifier for Atom {
    fn simplify_color(&self) -> Atom {
        self.simplify_color_with(ColorSimplifySettings::default())
    }

    fn simplify_color_with(&self, settings: ColorSimplifySettings) -> Atom {
        self.as_view().simplify_color_with(settings)
    }

    fn to_color_casimir(&self) -> Atom {
        self.to_color_casimir_with(ColorCasimirSettings::default())
    }

    fn to_color_casimir_with(&self, settings: ColorCasimirSettings) -> Atom {
        casimir::color_casimir_basis_impl(self.as_atom_view(), settings)
    }

    fn to_cof_dimension_invariants(&self) -> Atom {
        CofDimensionInvariantRewriter.run(self.as_atom_view())
    }

    fn expand_color(&self) -> Vec<(Atom, Atom)> {
        self.as_view().expand_color()
    }

    fn collect_color(&self) -> Atom {
        self.as_view().collect_color()
    }

    fn collect_color_constants(&self) -> Atom {
        self.as_view().collect_color_constants()
    }

    // fn canonize_color(&self) -> Atom {
    //     self.as_view().canonize_color()
    // }

    fn wrap_color(&self, symbol: Symbol) -> Atom {
        self.as_view().wrap_color(symbol)
    }
}

impl ColorSimplifier for AtomView<'_> {
    fn simplify_color(&self) -> Atom {
        self.simplify_color_with(ColorSimplifySettings::default())
    }

    fn simplify_color_with(&self, settings: ColorSimplifySettings) -> Atom {
        ColorAlgebraSimplifier { settings }.run(*self)
    }

    fn collect_color_constants(&self) -> Atom {
        self.collect_with_map(
            |a| {
                matches!(a,AtomView::Var(a) if a.get_symbol()==CS.tr || a.get_symbol()==CS.ca|| a.get_symbol()==CS.cf )
                    || matches!(a, AtomView::Fun(f) if f.get_symbol() == CS.cas || f.get_symbol() == CS.idx || f.get_symbol() == CS.gram)
            },
        )
        .unwrap_collect()
    }

    fn collect_color(&self) -> Atom {
        self.collect_reps([
            ColorAdjoint {}.into(),
            ColorFundamental {}.into(),
            ColorSextet {}.into(),
        ])
    }

    fn to_color_casimir(&self) -> Atom {
        self.to_color_casimir_with(ColorCasimirSettings::default())
    }

    fn to_color_casimir_with(&self, settings: ColorCasimirSettings) -> Atom {
        casimir::color_casimir_basis_impl(self.as_atom_view(), settings)
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
