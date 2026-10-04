#[cfg(test)]
use std::sync::LazyLock;

use spenso::{
    network::{library::symbolic::ExplicitKey, tags::SPENSO_TAG as T},
    shadowing::{
        IntoAtom,
        symbolica_utils::{SpensoPrintBackend, SpensoPrintSettings},
    },
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
    printer::{PrintOptions, PrintState},
    symbol,
};

use crate::{
    color::casimir::CofDimensionInvariantRewriter, representations::ColorAntiFundamental,
    shorthands::metric::PermuteWithMetric,
};

use super::rep_symbols::RS;
use super::representations::{ColorAdjoint, ColorFundamental};

mod casimir;
mod conjugate;
mod macros;
pub(crate) mod simplify;

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
    fn print(self, atom: AtomView<'_>, opt: &PrintOptions) -> Option<String> {
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

    fn representation(mut arg: AtomView<'_>, opt: &PrintOptions) -> String {
        let resolved = SpensoPrintSettings::resolve(opt).unwrap();
        let original = arg;
        let mut dual = false;
        if let AtomView::Fun(f) = arg
            && f.get_symbol() == AIND_SYMBOLS.dind
            && f.get_nargs() == 1
        {
            dual = true;
            arg = f.iter().next().unwrap();
        }
        let mut out = String::new();
        if let AtomView::Fun(rep) = arg
            && rep.get_symbol().has_tag(&T.representation)
            && rep.get_nargs() == 1
        {
            out = match rep.get_symbol().get_name() {
                "spenso::cof" => "F".to_owned(),
                "spenso::coad" => "A".to_owned(),
                _ => {
                    Atom::var(rep.get_symbol())
                        .format(&mut out, opt, PrintState::new())
                        .unwrap();
                    out
                }
            };
            if dual {
                out = match resolved.backend {
                    SpensoPrintBackend::Typst => format!("accent({out},macron)"),
                    SpensoPrintBackend::Latex => format!("\\overline{{{out}}}"),
                    SpensoPrintBackend::Plain => format!("bar({out})"),
                };
            }
            if resolved.presentation.with_dim {
                let mut dimension = String::new();
                rep.iter()
                    .next()
                    .unwrap()
                    .format(&mut dimension, opt, PrintState::new())
                    .unwrap();
                if resolved.presentation.symbol_scripts {
                    let (open, close) = resolved.script_delimiters();
                    out.push_str(&format!("_{open}{dimension}{close}"));
                } else {
                    out.push_str(&format!("[{dimension}]"));
                }
            }
        } else {
            original.format(&mut out, opt, PrintState::new()).unwrap();
        }
        out
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

    fn color_count(atom: AtomView<'_>, opt: &PrintOptions) -> Option<String> {
        if !matches!(atom, AtomView::Var(_)) {
            return None;
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
        AtomView::Fun(rep) if rep.get_symbol().get_name() == name && rep.get_nargs() == 1
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
        d: tensor_symbol!("spenso::d"),
        gram: symbol!("spenso::gram"; Real, Scalar; print = |a, opt, _state| {
            ColorInvariantPrintKind::Gram.print(a, opt)
        }),
        cas: symbol!("spenso::cas"; Real, Scalar; print = |a, opt, _state| {
            ColorInvariantPrintKind::Casimir.print(a, opt)
        }),
        idx: symbol!("spenso::idx"; Real, Scalar; print = |a, opt, _state| {
            ColorInvariantPrintKind::Index.print(a, opt)
        }),
        fundamental_rep: representation_symbol(
            ColorFundamental {}.to_symbolic(std::iter::empty::<Atom>()),
        ),
        adjoint_rep: representation_symbol(ColorAdjoint {}.to_symbolic(std::iter::empty::<Atom>())),
        adj_: symbol!("adj_"),
        nc_: symbol!("nc_"),
        nc: symbol!("spenso::Nc";Real; print = |a, opt, _| ColorInvariantPrintKind::color_count(a, opt),eval = EvaluationInfo::constant(|_tags, prec| Ok(Rational::new(3,1).to_multi_prec_float(prec).into()))),
    }
};
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColorSimplifySettings {
    /// Whether closed color chains should be evaluated as traces.
    pub evaluate_traces: bool,
    /// Whether contractions between generators on different open chains or traces should
    /// be expanded with the fundamental Fierz identity.
    pub expand_cross_chain_fierz: bool,
    /// Whether invariant factors for `cof(N)` should be written directly in
    /// terms of the fundamental dimension.
    pub substitute_cof_dimension_invariants: bool,
    /// Whether a trace line that nothing else in its term can reach is
    /// decomposed completely in one kernel call. The kernel then also
    /// certifies colour fixed points, so the planner needs no confirming
    /// colour round, and distributes a colour sum with ports only when a
    /// rule can span it. Without it, every insertion returns to the planner.
    pub one_shot_traces: bool,
}

impl Default for ColorSimplifySettings {
    fn default() -> Self {
        Self {
            evaluate_traces: true,
            expand_cross_chain_fierz: true,
            substitute_cof_dimension_invariants: false,
            one_shot_traces: true,
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

    /// Returns to the planner after each trace insertion and confirms every
    /// colour fixed point with a no-op round.
    pub fn without_one_shot_traces(mut self) -> Self {
        self.one_shot_traces = false;
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
}

#[cfg(test)]
mod test;
