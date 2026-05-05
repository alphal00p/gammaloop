use std::{collections::HashSet, sync::LazyLock};

use itertools::Itertools;
use spenso::{
    chain,
    network::{
        library::symbolic::{ETS, ExplicitKey},
        tags::SPENSO_TAG as T,
    },
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        TensorStructure,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        dimension::Dimension,
        representation::{Minkowski, RepName},
        slot::{AbsInd, IsAbstractSlot},
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
    function,
    id::{Context, MatchSettings, Pattern, Replacement},
    printer::PrintState,
    symbol,
    utils::Settable,
};

use crate::{chain::Chain, metric::PermuteWithMetric, representations::ColorAntiFundamental};

use super::rep_symbols::RS;
use super::{
    metric::MetricSimplifier,
    representations::{Bispinor, ColorAdjoint, ColorFundamental},
};

#[derive(Debug)]
pub enum ColorError {
    NotFully(Atom),
}

pub struct ColorSymbols {
    pub nc_: Symbol,
    pub adj_: Symbol,
    /// The adjoint Casimir symbol, i.e. CA = Nc
    pub ca: Symbol,
    /// The fundamental Casimir symbol.
    pub cf: Symbol,
    /// The generator symbol
    pub t: Symbol,
    /// The structure constant symbol i.e. [T^a, T^b] = i f^{abc} T^c
    pub f: Symbol,
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

pub static CS: LazyLock<ColorSymbols> = LazyLock::new(|| ColorSymbols {
    t: symbol!("spenso::t";Real;print = |a, opt| {

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
    f: symbol!("spenso::f";Real;print = |a, opt| {

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
    adj_: symbol!("adj_"),
    nc_: symbol!("nc_"),
    tr: symbol!("spenso::TR";Real),
    nc: symbol!("spenso::Nc";Real),
});

/// Builds a fundamental color-generator factor for use inside `chain!` or
/// `trace!`.
///
/// The adjoint-index argument is converted through
/// `spenso::symbolica_atom::IntoAtom`, so it can be a typed color-adjoint slot,
/// an atom, or an atom view. The factor uses the chain placeholder indices `in`
/// and `out`; the surrounding chain or trace owns the physical endpoints.
///
/// # Examples
///
/// ```ignore
/// use idenso::color_t;
/// use spenso::{slot, trace};
///
/// let factor = color_t!(slot!(coad_na, a));
/// let expr = trace!(&cof_nc, factor);
/// ```
#[macro_export]
macro_rules! color_t {
    ($a:expr) => {
        $crate::color::CS.chain_t(spenso::symbolica_atom::IntoAtom::into_atom($a))
    };
}

/// Builds an adjoint color structure-constant atom `f(a,b,c)`.
///
/// Arguments are converted through `spenso::symbolica_atom::IntoAtom`, so typed
/// color-adjoint slots can be used directly. This avoids accidentally mixing
/// parsed symbols that print alike but are not identical atoms.
///
/// # Examples
///
/// ```ignore
/// use idenso::color_f;
/// use spenso::slot;
///
/// let expr = color_f!(slot!(coad_na, a), slot!(coad_na, b), slot!(coad_na, c));
/// ```
#[macro_export]
macro_rules! color_f {
    ($a:expr, $b:expr, $c:expr $(,)?) => {
        $crate::color::CS.structure_f(
            spenso::symbolica_atom::IntoAtom::into_atom($a),
            spenso::symbolica_atom::IntoAtom::into_atom($b),
            spenso::symbolica_atom::IntoAtom::into_atom($c),
        )
    };
}

pub fn color_conj_impl(expression: AtomView) -> Atom {
    let expr = expression.to_owned();
    let cof = ColorFundamental {};
    let coaf = ColorFundamental {}.dual();

    let expr = expr
        .replace(
            function!(
                CS.t,
                RS.i_,
                cof.to_symbolic([RS.d_, RS.a_]),
                coaf.to_symbolic([RS.d_, RS.b_])
            )
            .to_pattern(),
        )
        .with(function!(
            CS.t,
            RS.i_,
            coaf.to_symbolic([RS.d_, RS.b_]),
            cof.to_symbolic([RS.d_, RS.a_])
        ));

    expr.replace_multiple(&[
        Replacement::new(
            coaf.to_symbolic([RS.a__]).to_pattern(),
            cof.to_symbolic([RS.a__]),
        ),
        Replacement::new(
            cof.to_symbolic([RS.a__]).to_pattern(),
            coaf.to_symbolic([RS.a__]),
        ),
    ])
}

pub trait SelectiveExpand {
    fn expand_in_patterns(&self, pats: &[Pattern]) -> Vec<(Atom, Atom)>;
    fn expand_metrics(&self) -> Vec<(Atom, Atom)> {
        let metric_pat = function!(ETS.metric, RS.a__).to_pattern();
        let id_pat = function!(ETS.metric, RS.a__).to_pattern();

        self.expand_in_patterns(&[metric_pat, id_pat])
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

    fn expand_bis(&self) -> Vec<(Atom, Atom)> {
        let bis = Bispinor {};

        let bis_pat = function!(RS.f_, RS.a___, bis.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[bis_pat])
    }

    fn expand_mink(&self) -> Vec<(Atom, Atom)> {
        let mink = Minkowski {};

        let mink_pat = function!(RS.f_, RS.a___, mink.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[mink_pat])
    }

    fn expand_mink_bis(&self) -> Vec<(Atom, Atom)> {
        let mink = Minkowski {};

        let mink_pat = function!(RS.f_, RS.a___, mink.to_symbolic([RS.b__]), RS.c___).to_pattern();

        let bis = Bispinor {};

        let bis_pat = function!(RS.f_, RS.a___, bis.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[mink_pat, bis_pat])
    }
}
impl SelectiveExpand for Atom {
    fn expand_in_patterns(&self, pats: &[Pattern]) -> Vec<(Atom, Atom)> {
        self.as_view().expand_in_patterns(pats)
    }
}

impl SelectiveExpand for AtomView<'_> {
    fn expand_in_patterns(&self, pats: &[Pattern]) -> Vec<(Atom, Atom)> {
        let mut coefs = HashSet::new();

        //A (x+y)(z*B+x*C)=> A(x*x*C+y*x*C+y*z*B+y*z*B)

        for p in pats {
            for m in self.pattern_match(p, None, None) {
                coefs.insert(p.replace_wildcards(&m));
            }
        }

        let coefs = coefs.into_iter().collect_vec();

        self.coefficient_list::<i8>(&coefs)
        // .coll
    }
}

fn simplify_raw_color_tensors(expression: AtomView) -> Atom {
    let tr = Atom::var(CS.tr);

    fn t(
        a: impl Into<AbstractIndex>,
        i: impl Into<AbstractIndex>,
        j: impl Into<AbstractIndex>,
    ) -> Atom {
        CS.t_pattern(CS.nc_, CS.adj_, a, i, j)
    }

    fn f(
        a: impl Into<AbstractIndex>,
        b: impl Into<AbstractIndex>,
        c: impl Into<AbstractIndex>,
    ) -> Atom {
        CS.f_pattern(CS.adj_, a, b, c)
    }

    let coad = ColorAdjoint {}.new_rep(CS.adj_);
    let cof = ColorFundamental {}.new_rep(CS.nc_);
    let coaf = cof.dual();

    let tpat = CS
        .t_pattern(CS.nc_, CS.adj_, RS.a_, RS.i_, RS.j_)
        .to_pattern();
    let mut ncs = None;
    for m in expression.pattern_match(&tpat, None, None) {
        if let Some(ncs) = &ncs
            && ncs != &m[&CS.nc_]
        {
            panic!("Mismatched Nc values in expression")
        } else {
            ncs = Some(m[&CS.nc_].clone());
        }
    }

    let expression = if let Some(ncs) = &ncs {
        expression
            .replace(ColorFundamental {}.to_symbolic([ncs.as_view(), Atom::var(RS.a_).as_view()]))
            .with(ColorFundamental {}.to_symbolic([CS.nc, RS.a_]))
    } else {
        expression.to_owned()
    };

    let reps = vec![
        (t(RS.a_, RS.b_, RS.b_), Atom::num(0)),
        (
            t(RS.a_, RS.i_, RS.j_) * t(RS.b_, RS.j_, RS.i_),
            &tr * coad.g(RS.a_, RS.b_),
        ),
        (
            t(RS.a_, RS.i_, RS.j_).pow(Atom::num(2)),
            &tr * coad.g(RS.a_, RS.a_),
        ),
        // Prefer the separated open-line Casimir shortcut over expanding the
        // outer generators with Fierz; the expanded form is much harder to
        // clean back into a single generator line.
        (
            t(RS.i_, RS.a_, RS.b_) * t(RS.e_, RS.b_, RS.c_) * t(RS.i_, RS.c_, RS.d_),
            -(&tr / Atom::var(CS.nc_)) * t(RS.e_, RS.a_, RS.d_),
        ),
        (
            t(RS.e_, RS.a_, RS.b_) * t(RS.e_, RS.c_, RS.d_),
            &tr * (coaf.id(RS.a_, RS.d_) * coaf.id(RS.c_, RS.b_)
                - (coaf.id(RS.a_, RS.b_) * coaf.id(RS.c_, RS.d_) / CS.nc_)),
        ),
    ];

    let i = symbol!("i");
    let j = symbol!("j");
    let k = symbol!("k");

    fn ta<'a>(
        a: impl Into<AbstractIndex>,
        i: impl Into<AtomOrView<'a>>,
        j: impl Into<AtomOrView<'a>>,
    ) -> Atom {
        function!(
            CS.t,
            ColorAdjoint {}.new_rep(CS.adj_).slot(a).to_atom(),
            ColorFundamental {}.new_rep(CS.nc).pattern(i),
            ColorAntiFundamental {}.new_rep(CS.nc).pattern(j)
        )
    }

    let frep = [
        Replacement::new(
            f(RS.a_, RS.b_, RS.c_).pow(Atom::num(2)).to_pattern(),
            CS.ca * CS.adj_,
        ),
        Replacement::new(
            f(RS.a_, RS.b_, RS.c_).to_pattern(),
            (((ta(
                RS.a_,
                function!(i, RS.a_, RS.b_, RS.c_),
                function!(j, RS.a_, RS.b_, RS.c_),
            ) * ta(
                RS.b_,
                function!(j, RS.a_, RS.b_, RS.c_),
                function!(k, RS.a_, RS.b_, RS.c_),
            ) * ta(
                RS.c_,
                function!(k, RS.a_, RS.b_, RS.c_),
                function!(i, RS.a_, RS.b_, RS.c_),
            ) - ta(
                RS.a_,
                function!(i, RS.a_, RS.b_, RS.c_),
                function!(j, RS.a_, RS.b_, RS.c_),
            ) * ta(
                RS.c_,
                function!(j, RS.a_, RS.b_, RS.c_),
                function!(k, RS.a_, RS.b_, RS.c_),
            ) * ta(
                RS.b_,
                function!(k, RS.a_, RS.b_, RS.c_),
                function!(i, RS.a_, RS.b_, RS.c_),
            )) / &tr)
                * -Atom::i())
            .to_pattern(),
        ),
    ];

    let settings = MatchSettings {
        rhs_cache_size: 0,
        ..Default::default()
    };
    let replacements: Vec<Replacement> = reps
        .into_iter()
        .map(|(a, b)| {
            Replacement::new(a.to_pattern(), b.to_pattern()).with_settings(settings.clone())
        })
        .collect();

    // for r in &replacements {
    //     println!("{r}")
    // }
    // for f in &frep {
    //     println!("{f}")
    // }
    let mut expression = expression.expand_color();

    for (e, _) in &mut expression {
        let mut atom = Atom::num(0);
        let mut first = true;
        while first || e.replace_multiple_into(&replacements, &mut atom) {
            if !first {
                std::mem::swap(e, &mut atom)
            };
            first = false;
            *e = e.replace_multiple(&frep);
            *e = e.expand();
            *e = e.simplify_metrics();
        }
    }

    // let pats: Vec<LibraryRep> = vec![ColorAdjoint {}.into()];
    // let dualizablepats: Vec<LibraryRep> = vec![ColorFundamental {}.into(), ColorSextet {}.into()];

    // let mut fully_simplified = true;
    // for p in pats.iter().chain(&dualizablepats) {
    //     if expression
    //         .pattern_match(&p.to_symbolic([RS.a__]).to_pattern(), None, None)
    //         .next()
    //         .is_some()
    //     {
    //         fully_simplified = false;
    //     }
    // }
    //
    let out = expression.iter().fold(Atom::Zero, |a, (c, s)| a + c * s);
    if let Some(ncs) = ncs {
        out.replace(CS.nc).with(ncs)
    } else {
        out
    }
}

static COLOR_FUNDAMENTAL_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
    let rep = ColorFundamental {}.to_symbolic(std::iter::empty::<Atom>());
    let AtomView::Fun(f) = rep.as_view() else {
        unreachable!("Color fundamental representations are symbolic functions")
    };

    f.get_symbol()
});

static COLOR_ADJOINT_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
    let rep = ColorAdjoint {}.to_symbolic(std::iter::empty::<Atom>());
    let AtomView::Fun(f) = rep.as_view() else {
        unreachable!("Color adjoint representations are symbolic functions")
    };

    f.get_symbol()
});

static COLOR_D_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| symbol!("spenso::d"));
static COLOR_DR_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| symbol!("spenso::dR"));
static COLOR_D33_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| symbol!("spenso::d33"));
static COLOR_TRACE_DUMMY_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| symbol!("x"));

/// Applies SU(N) simplification rules to raw color tensors and chain/trace form.
pub fn color_simplify_impl(expression: AtomView) -> Atom {
    let use_raw_tensor_simplifier =
        !contains_chain_or_trace(expression) && contains_raw_color_generator(expression);

    let mut expr = expression.to_owned().expand().simplify_metrics();
    // Keep legacy raw `t(a,i,j)` expressions on the old pattern path; chain
    // collection is only unambiguous once generators use chain endpoints.
    if use_raw_tensor_simplifier {
        return simplify_raw_color_tensors(expr.as_view())
            .expand()
            .simplify_metrics();
    }

    loop {
        let collected = collect_color_lines(expr.as_view());
        let next = rewrite_color_algebra_terms(collected.as_view())
            .expand()
            .simplify_metrics();

        if next == expr {
            return next;
        }

        expr = next;
    }
}

fn rewrite_color_algebra_terms(expr: AtomView) -> Atom {
    // Terminal trace rules can create sums; product rules such as f*f -> CA*g
    // then need to run on each generated term instead of on the whole Add.
    if let AtomView::Add(add) = expr {
        return add
            .iter()
            .map(rewrite_color_algebra_terms)
            .fold(Atom::Zero, |sum, term| sum + term);
    }

    if let Some(rewritten) = rewrite_color_node(expr) {
        return rewritten;
    }

    expr.to_owned().replace_map(&color_algebra_rewrite)
}

fn collect_color_lines(expr: AtomView) -> Atom {
    let rep = ColorFundamental {}.into();
    expr.to_owned()
        .chainify(rep)
        .collect_chains(ColorFundamental {}.into())
        .replace_map(&close_color_trace_chain)
}

fn close_color_trace_chain(arg: AtomView, _context: &Context, out: &mut Settable<'_, Atom>) {
    let Some((start, end, factors)) = chain_parts(arg) else {
        return;
    };
    let Some((dimension, start_index, false)) = color_fundamental_slot(start.as_view()) else {
        return;
    };
    let Some((end_dimension, end_index, true)) = color_fundamental_slot(end.as_view()) else {
        return;
    };
    if dimension != end_dimension || start_index != end_index {
        return;
    }

    **out = trace!(ColorFundamental {}.to_symbolic([dimension]); factors);
}

// Try product-level rewrites first so trace*f contractions can fire before the
// trace terminal expands into dR and f terms.
fn color_algebra_rewrite(arg: AtomView, _context: &Context, out: &mut Settable<'_, Atom>) {
    if let Some(rewritten) = rewrite_color_node(arg) {
        **out = rewritten;
    }
}

fn rewrite_color_node(arg: AtomView) -> Option<Atom> {
    simplify_color_product(arg)
        .or_else(|| simplify_color_chain_node(arg))
        .or_else(|| simplify_color_trace_node(arg))
        .or_else(|| simplify_color_power(arg))
}

fn simplify_color_chain_node(chain: AtomView) -> Option<Atom> {
    let (start, end, factors) = chain_parts(chain)?;
    if factors.is_empty()
        || factors
            .iter()
            .all(|factor| is_chain_identity_factor(factor.as_view()))
    {
        return Some(color_metric(start.clone(), end.clone()));
    }

    for i in 0..factors.len().saturating_sub(1) {
        let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
            continue;
        };
        let Some(right) = color_generator_adjoint(factors[i + 1].as_view()) else {
            continue;
        };
        if left != right {
            continue;
        }

        return Some(Atom::var(CS.cf) * chain_with_removed_range(&start, &end, &factors, i, i + 2));
    }

    None
}

fn simplify_color_trace_node(trace: AtomView) -> Option<Atom> {
    let (rep, factors) = trace_parts(trace)?;

    if factors.is_empty()
        || factors
            .iter()
            .all(|factor| is_chain_identity_factor(factor.as_view()))
    {
        return trace_terminal_dimension(rep.as_view());
    }

    if factors.len() > 2
        && let Some(rewritten) = simplify_adjacent_trace_casimir(&rep, &factors)
    {
        return Some(rewritten);
    }

    if factors.len() > 3
        && let Some(rewritten) = simplify_separated_trace_casimir(&rep, &factors)
    {
        return Some(rewritten);
    }

    let generators = factors
        .iter()
        .map(|factor| color_generator_adjoint(factor.as_view()))
        .collect::<Option<Vec<_>>>()?;

    match generators.as_slice() {
        [_] => Some(Atom::Zero),
        [a, b] => Some(Atom::var(CS.tr) * color_metric(a.clone(), b.clone())),
        [a, b, c] => Some(
            color_d_r([a.clone(), b.clone(), c.clone()])
                + Atom::i() * Atom::num(1) / Atom::num(2)
                    * Atom::var(CS.tr)
                    * color_f([a.clone(), b.clone(), c.clone()]),
        ),
        [a, b, c, d] => simplify_four_generator_trace_terminal(a, b, c, d),
        _ => None,
    }
}

fn simplify_adjacent_trace_casimir(rep: &Atom, factors: &[Atom]) -> Option<Atom> {
    for i in 0..factors.len().saturating_sub(1) {
        let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
            continue;
        };
        let Some(right) = color_generator_adjoint(factors[i + 1].as_view()) else {
            continue;
        };
        if left != right {
            continue;
        }

        return Some(
            Atom::var(CS.cf)
                * trace_with_factors(rep.clone(), factors_excluding_range(factors, i, i + 2)),
        );
    }

    None
}

fn simplify_separated_trace_casimir(rep: &Atom, factors: &[Atom]) -> Option<Atom> {
    for i in 0..factors.len().saturating_sub(2) {
        let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
            continue;
        };
        let Some(right) = color_generator_adjoint(factors[i + 2].as_view()) else {
            continue;
        };
        if left != right {
            continue;
        }

        return Some(
            (Atom::var(CS.cf) - Atom::var(CS.ca) / Atom::num(2))
                * trace_with_factors(rep.clone(), factors_excluding_indices(factors, &[i, i + 2])),
        );
    }

    None
}

fn simplify_four_generator_trace_terminal(a: &Atom, b: &Atom, c: &Atom, d: &Atom) -> Option<Atom> {
    let x = color_adjoint_dummy_like(a)?;

    Some(
        color_d_r([a.clone(), b.clone(), c.clone(), d.clone()])
            + Atom::i() / Atom::num(2)
                * color_d_r([a.clone(), b.clone(), x.clone()])
                * color_f([c.clone(), d.clone(), x.clone()])
            + Atom::i() / Atom::num(2)
                * color_d_r([c.clone(), d.clone(), x.clone()])
                * color_f([a.clone(), b.clone(), x.clone()])
            - Atom::var(CS.tr) / Atom::num(6)
                * color_f([a.clone(), c.clone(), x.clone()])
                * color_f([b.clone(), d.clone(), x.clone()])
            + Atom::var(CS.tr) / Atom::num(3)
                * color_f([a.clone(), d.clone(), x.clone()])
                * color_f([b.clone(), c.clone(), x]),
    )
}

fn simplify_color_product(product: AtomView) -> Option<Atom> {
    let factors = multiplicative_factors(product);
    if factors.len() < 2 {
        return None;
    }

    join_color_chain_product(&factors)
        .or_else(|| simplify_trace_structure_product(&factors))
        .or_else(|| simplify_symmetric_structure_product(&factors))
        .or_else(|| simplify_two_f_loop_product(&factors))
        .or_else(|| simplify_three_f_loop_product(&factors))
        .or_else(|| simplify_d33_product(&factors))
}

fn join_color_chain_product(factors: &[Atom]) -> Option<Atom> {
    for (left_index, left_factor) in factors.iter().enumerate() {
        let Some((left_start, left_end, left_factors)) = chain_parts(left_factor.as_view()) else {
            continue;
        };
        let Some((left_end_dim, left_end_index, true)) = color_fundamental_slot(left_end.as_view())
        else {
            continue;
        };

        for (right_index, right_factor) in factors.iter().enumerate() {
            if right_index == left_index {
                continue;
            }
            let Some((right_start, right_end, right_factors)) = chain_parts(right_factor.as_view())
            else {
                continue;
            };
            let Some((right_start_dim, right_start_index, false)) =
                color_fundamental_slot(right_start.as_view())
            else {
                continue;
            };
            if left_end_dim != right_start_dim || left_end_index != right_start_index {
                continue;
            }

            let replacement = chain!(
                left_start,
                right_end;
                left_factors.into_iter().chain(right_factors.into_iter())
            );
            return Some(product_replacing_pair(
                factors,
                left_index,
                right_index,
                replacement,
            ));
        }
    }

    None
}

fn simplify_color_power(power: AtomView) -> Option<Atom> {
    let AtomView::Pow(pow) = power else {
        return None;
    };
    let (base, exponent) = pow.get_base_exp();
    if positive_integer(exponent)? != 2 {
        return None;
    }

    let args = structure_constant_args(base)?;
    let dimension = color_adjoint_dimension(&args[0])?;
    Some(Atom::var(CS.ca) * dimension)
}

fn simplify_trace_structure_product(factors: &[Atom]) -> Option<Atom> {
    for (trace_index, trace_factor) in factors.iter().enumerate() {
        let Some((rep, trace_factors)) = trace_parts(trace_factor.as_view()) else {
            continue;
        };
        let [first, second, rest @ ..] = trace_factors.as_slice() else {
            continue;
        };
        let Some(a) = color_generator_adjoint(first.as_view()) else {
            continue;
        };
        let Some(b) = color_generator_adjoint(second.as_view()) else {
            continue;
        };

        for (f_index, f_factor) in factors.iter().enumerate() {
            if f_index == trace_index {
                continue;
            }
            let Some([fa, fb, fc]) = structure_constant_args(f_factor.as_view()) else {
                continue;
            };
            if fa != a || fb != b {
                continue;
            }

            let replacement = Atom::i() * Atom::var(CS.ca) / Atom::num(2)
                * trace!(rep.clone(); std::iter::once(CS.chain_t(fc)).chain(rest.iter().cloned()));
            return Some(product_replacing_pair(
                factors,
                trace_index,
                f_index,
                replacement,
            ));
        }
    }

    None
}

fn simplify_symmetric_structure_product(factors: &[Atom]) -> Option<Atom> {
    for (symmetric_index, symmetric_factor) in factors.iter().enumerate() {
        let Some(symmetric_args) = color_symmetric_args(symmetric_factor.as_view()) else {
            continue;
        };

        for (f_index, f_factor) in factors.iter().enumerate() {
            if f_index == symmetric_index {
                continue;
            }
            let Some(structure_args) = structure_constant_args(f_factor.as_view()) else {
                continue;
            };
            let common_count = symmetric_args
                .iter()
                .filter(|arg| structure_args.iter().any(|candidate| candidate == *arg))
                .count();
            if common_count >= 2 {
                return Some(Atom::Zero);
            }
        }
    }

    None
}

fn simplify_two_f_loop_product(factors: &[Atom]) -> Option<Atom> {
    for (left_index, left_factor) in factors.iter().enumerate() {
        let Some(left) = structure_constant_args(left_factor.as_view()) else {
            continue;
        };

        for (right_index, right_factor) in factors.iter().enumerate().skip(left_index + 1) {
            let Some(right) = structure_constant_args(right_factor.as_view()) else {
                continue;
            };
            let common = common_atoms(&left, &right);
            if common.len() != 2 {
                continue;
            }
            let Some(left_open) = first_not_in(&left, &common) else {
                continue;
            };
            let Some(right_open) = first_not_in(&right, &common) else {
                continue;
            };

            let replacement =
                Atom::var(CS.ca) * color_metric(left_open.clone(), right_open.clone());
            return Some(product_replacing_pair(
                factors,
                left_index,
                right_index,
                replacement,
            ));
        }
    }

    None
}

fn simplify_three_f_loop_product(factors: &[Atom]) -> Option<Atom> {
    for indices in (0..factors.len()).combinations(3) {
        let Some(f_args) = indices
            .iter()
            .map(|index| structure_constant_args(factors[*index].as_view()))
            .collect::<Option<Vec<_>>>()
        else {
            continue;
        };

        let all_args = f_args.iter().flatten().cloned().collect::<Vec<_>>();
        let externals = all_args
            .iter()
            .filter(|arg| {
                all_args
                    .iter()
                    .filter(|candidate| *candidate == *arg)
                    .count()
                    == 1
            })
            .cloned()
            .collect::<Vec<_>>();
        if externals.len() != 3 {
            continue;
        }

        let replacement = Atom::var(CS.ca) / Atom::num(2) * color_f(externals);
        let mut excluded = vec![false; factors.len()];
        for index in indices {
            excluded[index] = true;
        }
        return Some(product_excluding(factors, &excluded) * replacement);
    }

    None
}

fn simplify_d33_product(factors: &[Atom]) -> Option<Atom> {
    for (left_index, left_factor) in factors.iter().enumerate() {
        let Some((left_rep, left_args)) = symmetric_invariant_args(left_factor.as_view()) else {
            continue;
        };
        for (right_index, right_factor) in factors.iter().enumerate().skip(left_index + 1) {
            let Some((right_rep, right_args)) = symmetric_invariant_args(right_factor.as_view())
            else {
                continue;
            };
            let [la, lb, lc] = left_args.as_slice() else {
                continue;
            };
            let [ra, rb, rc] = right_args.as_slice() else {
                continue;
            };
            if la != ra || lb != rb {
                continue;
            }
            let dimension = color_adjoint_dimension(la)?;
            let replacement =
                color_d33(left_rep, right_rep) * color_metric(lc.clone(), rc.clone()) / dimension;
            return Some(product_replacing_pair(
                factors,
                left_index,
                right_index,
                replacement,
            ));
        }
    }

    None
}

fn chain_parts(chain: AtomView) -> Option<(Atom, Atom, Vec<Atom>)> {
    let AtomView::Fun(f) = chain else {
        return None;
    };
    if f.get_symbol() != T.chain {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    let [start, end, factors @ ..] = args.as_slice() else {
        return None;
    };
    Some((start.clone(), end.clone(), factors.to_vec()))
}

fn trace_parts(trace: AtomView) -> Option<(Atom, Vec<Atom>)> {
    let AtomView::Fun(f) = trace else {
        return None;
    };
    if f.get_symbol() != T.trace {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    let [rep, factors @ ..] = args.as_slice() else {
        return None;
    };
    Some((rep.clone(), factors.to_vec()))
}

fn color_generator_adjoint(factor: AtomView) -> Option<Atom> {
    let AtomView::Fun(f) = factor else {
        return None;
    };
    if f.get_symbol() != CS.t || f.get_nargs() != 3 {
        return None;
    }

    let args = f.iter().collect::<Vec<_>>();
    if !has_chain_endpoints(args[1], args[2]) {
        return None;
    }

    Some(args[0].to_owned())
}

fn structure_constant_args(factor: AtomView) -> Option<[Atom; 3]> {
    let AtomView::Fun(f) = factor else {
        return None;
    };
    if f.get_symbol() != CS.f || f.get_nargs() != 3 {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    Some([args[0].clone(), args[1].clone(), args[2].clone()])
}

fn symmetric_invariant_args(invariant: AtomView) -> Option<(Atom, Vec<Atom>)> {
    let AtomView::Fun(f) = invariant else {
        return None;
    };
    if f.get_symbol() != *COLOR_D_SYMBOL || f.get_nargs() != 4 {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    Some((args[0].clone(), args[1..].to_vec()))
}

fn color_symmetric_args(invariant: AtomView) -> Option<Vec<Atom>> {
    let AtomView::Fun(f) = invariant else {
        return None;
    };

    if f.get_symbol() == *COLOR_DR_SYMBOL && f.get_nargs() >= 2 {
        return Some(f.iter().map(|arg| arg.to_owned()).collect());
    }
    if f.get_symbol() == *COLOR_D_SYMBOL && f.get_nargs() >= 3 {
        let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
        return Some(args[1..].to_vec());
    }

    None
}

fn color_fundamental_slot(slot: AtomView) -> Option<(Atom, Atom, bool)> {
    if let Some((dimension, index)) = representation_slot(slot, *COLOR_FUNDAMENTAL_SYMBOL) {
        return Some((dimension, index, false));
    }

    let AtomView::Fun(f) = slot else {
        return None;
    };
    if f.get_symbol() != AIND_SYMBOLS.dind || f.get_nargs() != 1 {
        return None;
    }

    representation_slot(f.iter().next()?, *COLOR_FUNDAMENTAL_SYMBOL)
        .map(|(dimension, index)| (dimension, index, true))
}

fn color_adjoint_dimension(slot: &Atom) -> Option<Atom> {
    representation_slot(slot.as_view(), *COLOR_ADJOINT_SYMBOL).map(|(dimension, _)| dimension)
}

fn color_adjoint_dummy_like(slot: &Atom) -> Option<Atom> {
    let dimension = color_adjoint_dimension(slot)?;
    Some(ColorAdjoint {}.to_symbolic([dimension, Atom::var(*COLOR_TRACE_DUMMY_SYMBOL)]))
}

fn representation_slot(slot: AtomView, symbol: Symbol) -> Option<(Atom, Atom)> {
    let AtomView::Fun(f) = slot else {
        return None;
    };
    if f.get_symbol() != symbol || f.get_nargs() != 2 {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    Some((args[0].clone(), args[1].clone()))
}

fn trace_terminal_dimension(rep: AtomView) -> Option<Atom> {
    let AtomView::Fun(f) = rep else {
        return None;
    };
    if !f.get_symbol().has_tag(&T.representation) || f.get_nargs() == 0 {
        return None;
    }

    f.iter().next().map(|dimension| dimension.to_owned())
}

fn has_chain_endpoints(left: AtomView, right: AtomView) -> bool {
    is_chain_endpoint(left, T.chain_in) && is_chain_endpoint(right, T.chain_out)
}

fn is_chain_identity_factor(factor: AtomView) -> bool {
    // `chain!` represents an empty line as a metric over the placeholder
    // endpoints; collapse it to the physical endpoint metric/dimension.
    let AtomView::Fun(f) = factor else {
        return false;
    };
    if f.get_symbol() != ETS.metric || f.get_nargs() != 2 {
        return false;
    }

    let args = f.iter().collect::<Vec<_>>();
    has_chain_endpoints(args[0], args[1])
}

fn is_chain_endpoint(arg: AtomView, expected: Symbol) -> bool {
    matches!(arg, AtomView::Var(symbol) if symbol.get_symbol() == expected)
}

fn color_metric(left: Atom, right: Atom) -> Atom {
    function!(ETS.metric, left, right)
}

fn color_f(factors: impl IntoIterator<Item = Atom>) -> Atom {
    factors
        .into_iter()
        .fold(FunctionBuilder::new(CS.f), |builder, factor| {
            builder.add_arg(factor)
        })
        .finish()
}

fn color_d_r(factors: impl IntoIterator<Item = Atom>) -> Atom {
    factors
        .into_iter()
        .fold(FunctionBuilder::new(*COLOR_DR_SYMBOL), |builder, factor| {
            builder.add_arg(factor)
        })
        .finish()
}

fn color_d33(left_rep: Atom, right_rep: Atom) -> Atom {
    function!(*COLOR_D33_SYMBOL, left_rep, right_rep)
}

fn chain_with_removed_range(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    from: usize,
    to: usize,
) -> Atom {
    let remaining = factors_excluding_range(factors, from, to);
    if remaining.is_empty() {
        color_metric(start.clone(), end.clone())
    } else {
        chain!(start.clone(), end.clone(); remaining)
    }
}

fn trace_with_factors(rep: Atom, factors: Vec<Atom>) -> Atom {
    if factors.is_empty() {
        trace_terminal_dimension(rep.as_view()).unwrap_or(rep)
    } else {
        trace!(rep; factors)
    }
}

fn factors_excluding_range(factors: &[Atom], from: usize, to: usize) -> Vec<Atom> {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| *index < from || *index >= to)
        .map(|(_, factor)| factor.clone())
        .collect()
}

fn factors_excluding_indices(factors: &[Atom], excluded: &[usize]) -> Vec<Atom> {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| !excluded.contains(index))
        .map(|(_, factor)| factor.clone())
        .collect()
}

fn multiplicative_factors(expr: AtomView) -> Vec<Atom> {
    match expr {
        AtomView::Mul(mul) => mul.iter().map(|factor| factor.to_owned()).collect(),
        _ => vec![expr.to_owned()],
    }
}

fn product_replacing_pair(
    factors: &[Atom],
    left_index: usize,
    right_index: usize,
    replacement: Atom,
) -> Atom {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != left_index && *index != right_index)
        .fold(replacement, |product, (_, factor)| product * factor)
}

fn product_excluding(factors: &[Atom], excluded: &[bool]) -> Atom {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| !excluded[*index])
        .fold(Atom::num(1), |product, (_, factor)| product * factor)
}

fn common_atoms(left: &[Atom; 3], right: &[Atom; 3]) -> Vec<Atom> {
    left.iter()
        .filter(|arg| right.iter().any(|candidate| candidate == *arg))
        .cloned()
        .collect()
}

fn first_not_in<'a>(args: &'a [Atom; 3], excluded: &[Atom]) -> Option<&'a Atom> {
    args.iter()
        .find(|arg| !excluded.iter().any(|candidate| candidate == *arg))
}

fn positive_integer(expr: AtomView) -> Option<i64> {
    let AtomView::Num(number) = expr else {
        return None;
    };
    let CoefficientView::Natural(value, 1, 0, 1) = number.get_coeff_view() else {
        return None;
    };

    (value > 0).then_some(value)
}

fn contains_raw_color_generator(expr: AtomView) -> bool {
    match expr {
        AtomView::Fun(f) => {
            if f.get_symbol() == CS.t {
                let args = f.iter().collect::<Vec<_>>();
                return args.len() != 3 || !has_chain_endpoints(args[1], args[2]);
            }

            f.iter().any(contains_raw_color_generator)
        }
        AtomView::Add(add) => add.iter().any(contains_raw_color_generator),
        AtomView::Mul(mul) => mul.iter().any(contains_raw_color_generator),
        AtomView::Pow(pow) => pow.iter().any(contains_raw_color_generator),
        AtomView::Num(_) | AtomView::Var(_) => false,
    }
}

fn contains_chain_or_trace(expr: AtomView) -> bool {
    match expr {
        AtomView::Fun(f) => {
            f.get_symbol() == T.chain
                || f.get_symbol() == T.trace
                || f.iter().any(contains_chain_or_trace)
        }
        AtomView::Add(add) => add.iter().any(contains_chain_or_trace),
        AtomView::Mul(mul) => mul.iter().any(contains_chain_or_trace),
        AtomView::Pow(pow) => pow.iter().any(contains_chain_or_trace),
        AtomView::Num(_) | AtomView::Var(_) => false,
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

    // fn canonize_color(&self) -> Atom;

    fn wrap_color(&self, symbol: Symbol) -> Atom;
}
impl ColorSimplifier for Atom {
    fn simplify_color(&self) -> Atom {
        color_simplify_impl(self.as_atom_view())
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
        color_simplify_impl(self.as_atom_view())
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
