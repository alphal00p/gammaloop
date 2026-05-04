use std::{collections::HashSet, sync::LazyLock};

use itertools::Itertools;
use spenso::{
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
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    function,
    id::{MatchSettings, Pattern, Replacement},
    printer::PrintState,
    symbol,
};

use crate::{metric::PermuteWithMetric, representations::ColorAntiFundamental};

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
    ca: symbol!("ca"),
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

pub fn color_simplify_impl(expression: AtomView) -> Atom {
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
        (
            t(RS.e_, RS.a_, RS.b_) * t(RS.e_, RS.c_, RS.d_),
            &tr * (coaf.id(RS.a_, RS.d_) * coaf.id(RS.c_, RS.b_)
                - (coaf.id(RS.a_, RS.b_) * coaf.id(RS.c_, RS.d_) / CS.nc_)),
        ),
        (
            t(RS.i_, RS.a_, RS.b_) * t(RS.e_, RS.b_, RS.c_) * t(RS.i_, RS.c_, RS.d_),
            -(&tr / Atom::var(CS.nc_)) * t(RS.e_, RS.a_, RS.d_),
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
