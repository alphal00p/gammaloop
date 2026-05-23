use std::sync::LazyLock;

use spenso::{
    network::{library::symbolic::ExplicitKey, tags::SPENSO_TAG as T},
    shadowing::{Collectable, TensorCollectExt, symbolica_utils::SpensoPrintSettings},
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
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    function,
    printer::PrintState,
    symbol,
};

use crate::{
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
pub use simplify::{color_simplify_impl, color_simplify_with_impl};

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
    /// The fundamental Casimir symbol.
    pub cf: Symbol,
    /// The generator symbol
    pub t: Symbol,
    /// The structure constant symbol i.e. [T^a, T^b] = i f^{abc} T^c
    pub f: Symbol,
    /// The symmetric color invariant symbol.
    pub d: Symbol,
    /// The contracted rank-three symmetric invariant symbol.
    pub d33: Symbol,
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

    pub fn d33<'a, 'b>(
        &self,
        left: impl Into<AtomOrView<'a>>,
        right: impl Into<AtomOrView<'b>>,
    ) -> Atom {
        FunctionBuilder::new(self.d33)
            .add_arg(left)
            .add_arg(right)
            .finish()
    }

    #[cfg(test)]
    pub(crate) fn initialize_tensor_symbols(&self) {
        let _ = self.t;
        let _ = self.f;
        let _ = self.d;
        let _ = self.d33;
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

pub static CS: LazyLock<ColorSymbols> = LazyLock::new(|| {
    fn representation_symbol(rep: Atom) -> Symbol {
        let AtomView::Fun(f) = rep.as_view() else {
            unreachable!("Color representations are symbolic functions")
        };
        f.get_symbol()
    }

    ColorSymbols {
        t: tensor_symbol!("spenso::t"; Real; print = |a, opt, _state| {

            match opt.custom_print_mode {
                Some(("spenso",i))=>{
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = SpensoPrintSettings::from(i);


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

            match opt.custom_print_mode {
                Some(("spenso",i))=>{
                    let SpensoPrintSettings{
                        parens,
                        symbol_scripts,
                        commas,..
                    } = SpensoPrintSettings::from(i);


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
        ca: symbol!("spenso::CA";Real),
        cf: symbol!("spenso::CF";Real),
        d: tensor_symbol!("spenso::d"),
        d33: tensor_symbol!("spenso::d33"),
        trace_dummy: symbol!("x"),
        fundamental_rep: representation_symbol(
            ColorFundamental {}.to_symbolic(std::iter::empty::<Atom>()),
        ),
        adjoint_rep: representation_symbol(ColorAdjoint {}.to_symbolic(std::iter::empty::<Atom>())),
        adj_: symbol!("adj_"),
        nc_: symbol!("nc_"),
        na: symbol!("NA";Real),
        tr: symbol!("spenso::TR";Real),
        nc: symbol!("spenso::Nc";Real),
    }
});

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColorSimplifySettings {
    /// Whether closed color chains should be evaluated as traces.
    pub evaluate_traces: bool,
    /// Whether contractions between generators on different open chains should
    /// be expanded with the fundamental Fierz identity.
    pub expand_cross_chain_fierz: bool,
}

impl Default for ColorSimplifySettings {
    fn default() -> Self {
        Self {
            evaluate_traces: true,
            expand_cross_chain_fierz: true,
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
        color_simplify_with_impl(self.as_atom_view(), settings)
    }

    fn collect_color_constants(&self) -> Atom {
        self.collect_with_map(
            |a| matches!(a,AtomView::Var(a) if a.get_symbol()==CS.tr || a.get_symbol()==CS.ca|| a.get_symbol()==CS.cf ),
        )
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

    fn expand_color(&self) -> Vec<(Atom, Atom)> {
        let cof = ColorFundamental {};
        let coaf = ColorFundamental {}.dual();
        let coad = ColorAdjoint {};

        let cof_pat = function!(RS.f_, RS.a___, cof.to_symbolic([RS.b__]), RS.c___).to_pattern();
        let coaf_pat = function!(RS.f_, RS.a___, coaf.to_symbolic([RS.b__]), RS.c___).to_pattern();
        let coad_pat = function!(RS.f_, RS.a___, coad.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[cof_pat, coad_pat, coaf_pat])
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
