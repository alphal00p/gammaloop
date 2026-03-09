use std::{collections::HashSet, sync::LazyLock};

use itertools::Itertools;
use spenso::{
    network::library::symbolic::{ETS, ExplicitKey},
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
    atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol},
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
mod test {

    use insta::assert_snapshot;
    use spenso::shadowing::symbolica_utils::AtomCoreExt;
    use spenso::structure::IndexlessNamedStructure;
    use spenso::structure::PermutedStructure;

    static _CF: LazyLock<PermutedStructure<IndexlessNamedStructure<Symbol, ()>>> =
        LazyLock::new(|| {
            IndexlessNamedStructure::from_iter(
                [
                    ColorAdjoint {}.new_rep(8),
                    ColorAdjoint {}.new_rep(2),
                    ColorAdjoint {}.new_rep(2),
                    ColorAdjoint {}.new_rep(4),
                    ColorAdjoint {}.new_rep(2),
                ],
                CS.f,
                None,
            )
        });

    static CT: LazyLock<PermutedStructure<IndexlessNamedStructure<Symbol, ()>>> =
        LazyLock::new(|| {
            IndexlessNamedStructure::from_iter(
                [
                    ColorAntiFundamental {}.new_rep(3).to_lib(),
                    ColorFundamental {}.new_rep(3).to_lib(),
                    ColorAdjoint {}.new_rep(8).to_lib(),
                ],
                CS.t,
                None,
            )
        });
    use spenso::{
        network::parsing::ShadowedStructure, structure::permuted::Perm,
        tensors::symbolic::SymbolicTensor,
    };
    use symbolica::{parse, parse_lit};

    use crate::gamma::PS;
    use crate::test::test_initialize;
    use crate::{
        IndexTooling, gamma::GammaSimplifier, metric::MetricSimplifier, representations::initialize,
    };

    use super::*;

    #[test]
    fn test_color_structures() {
        let f = IndexlessNamedStructure::<Symbol, ()>::from_iter(
            [
                ColorAdjoint {}.new_rep(8),
                ColorAdjoint {}.new_rep(2),
                ColorAdjoint {}.new_rep(2),
                ColorAdjoint {}.new_rep(4),
                ColorAdjoint {}.new_rep(2),
                ColorAdjoint {}.new_rep(7),
            ],
            symbol!("test"),
            None,
        )
        .clone()
        .reindex([5, 4, 2, 3, 1, 0])
        .unwrap()
        .map_structure(|a| SymbolicTensor::from_named(&a).unwrap());

        let f_s = f.structure.structure.clone();

        let f_p = f.clone().permute_inds();

        let f_parsed = PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(
            &f_p.expression.simplify_metrics(),
        )
        .unwrap();

        assert_eq!(f.index_permutation, f_parsed.index_permutation);
        assert!(f_parsed.rep_permutation.is_identity());

        let f_p = f.clone().permute_reps_wrapped().permute_inds();

        let f_parsed = PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(
            &f_p.expression.simplify_metrics(),
        )
        .unwrap();

        assert_eq!(f.index_permutation, f_parsed.index_permutation);
        assert_eq!(f.rep_permutation, f_parsed.rep_permutation);

        println!("Parsed: {}", f_parsed.index_permutation);
        println!("OG: {}", f.index_permutation);

        println!(
            "Structure:{}\nPermuted:{}\nPermuted Structure{}\nMetric simplified{}",
            f_s,
            f_p,
            f_p.structure,
            f_p.expression.simplify_metrics()
        );

        let t = CT
            .clone()
            .reindex([4, 2, 3])
            .unwrap()
            .map_structure(|a| SymbolicTensor::from_named(&a).unwrap())
            .permute_inds();

        println!("{t}")
    }

    #[test]
    fn test_color_simplification() {
        initialize();
        let atom = parse_lit!(
            f(
                coad(Nc ^ 2 - 1, 1),
                coad(Nc ^ 2 - 1, 2),
                coad(Nc ^ 2 - 1, 3)
            ) ^ 2,
            default_namespace = "spenso"
        );
        let simplified = atom.simplify_color();

        assert_snapshot!(simplified.to_bare_ordered_string(), @"-1*ca+Nc^2*ca");
    }

    fn colored_matrix_element() -> (Atom, Atom) {
        (
            parse!(
                "-G^2
            *(
                -g(mink(D,5),mink(D,6))*Q(2,mink(D,7))
                +g(mink(D,5),mink(D,6))*Q(3,mink(D,7))
                +g(mink(D,5),mink(D,7))*Q(2,mink(D,6))
                +g(mink(D,5),mink(D,7))*Q(4,mink(D,6))
                -g(mink(D,6),mink(D,7))*Q(3,mink(D,5))
                -g(mink(D,6),mink(D,7))*Q(4,mink(D,5))
            )
            *g(mink(D,4),mink(D,7))
            *t(coad(Nc^2-1,6),cof(Nc,5),dind(cof(Nc,4)))
            *f(coad(Nc^2-1,7),coad(Nc^2-1,8),coad(Nc^2-1,9))
            *g(bis(D,0),bis(D,5))
            *g(bis(D,1),bis(D,4))
            *g(mink(D,2),mink(D,5))
            *g(mink(D,3),mink(D,6))
            *g(coad(Nc^2-1,2),coad(Nc^2-1,7))
            *g(coad(Nc^2-1,3),coad(Nc^2-1,8))
            *g(coad(Nc^2-1,6),coad(Nc^2-1,9))
            *g(cof(Nc,0),dind(cof(Nc,5)))
            *g(cof(Nc,4),dind(cof(Nc,1)))
            *gamma(bis(D,5),bis(D,4),mink(D,4))
            *vbar(1,bis(D,1))
            *u(0,bis(D,0))
            *ϵbar(2,mink(D,2))
            *ϵbar(3,mink(D,3))",
                default_namespace = "spenso"
            ),
            parse_lit!(
                -4 * TR
                    ^ 2 * Nc * G
                    ^ 4 * (Nc - 1) * (Nc + 1) * (D - 2)
                    ^ -2 * (-2 * dot(Q(0), Q(1)) * dot(Q(2), Q(2))
                        + dot(Q(0), Q(1)) * dot(Q(2), Q(3))
                        - 3 * dot(Q(0), Q(1)) * dot(Q(2), Q(4))
                        - 2 * dot(Q(0), Q(1)) * dot(Q(3), Q(3))
                        - 3 * dot(Q(0), Q(1)) * dot(Q(3), Q(4))
                        - 3 * dot(Q(0), Q(1)) * dot(Q(4), Q(4))
                        + 2 * dot(Q(0), Q(2)) * dot(Q(1), Q(2))
                        - dot(Q(0), Q(2)) * dot(Q(1), Q(3))
                        + dot(Q(0), Q(2)) * dot(Q(1), Q(4))
                        - dot(Q(0), Q(3)) * dot(Q(1), Q(2))
                        + 2 * dot(Q(0), Q(3)) * dot(Q(1), Q(3))
                        + dot(Q(0), Q(3)) * dot(Q(1), Q(4))
                        + dot(Q(0), Q(4)) * dot(Q(1), Q(2))
                        + dot(Q(0), Q(4)) * dot(Q(1), Q(3))
                        + 2 * dot(Q(0), Q(4)) * dot(Q(1), Q(4))
                        + D * dot(Q(0), Q(1)) * dot(Q(2), Q(2))
                        - D * dot(Q(0), Q(1)) * dot(Q(2), Q(3))
                        + D * dot(Q(0), Q(1)) * dot(Q(2), Q(4))
                        + D * dot(Q(0), Q(1)) * dot(Q(3), Q(3))
                        + D * dot(Q(0), Q(1)) * dot(Q(3), Q(4))
                        + D * dot(Q(0), Q(1)) * dot(Q(4), Q(4))
                        - D * dot(Q(0), Q(2)) * dot(Q(1), Q(2))
                        + D * dot(Q(0), Q(2)) * dot(Q(1), Q(3))
                        + D * dot(Q(0), Q(3)) * dot(Q(1), Q(2))
                        - D * dot(Q(0), Q(3)) * dot(Q(1), Q(3))),
                default_namespace = "spenso"
            ),
        )
    }

    #[test]
    fn compact_printing() {
        initialize();

        let (atom, _) = colored_matrix_element();
        println!(
            "{}",
            atom.printer(SpensoPrintSettings::compact().nice_symbolica())
        );

        println!(
            "{}",
            atom.printer(
                SpensoPrintSettings {
                    parens: true,
                    with_dim: false,
                    commas: false,
                    index_subscripts: false,
                    symbol_scripts: false,
                }
                .nice_symbolica()
            )
        );

        println!(
            "{}",
            atom.printer(SpensoPrintSettings::typst().nice_symbolica())
        )
    }

    #[test]
    fn t_structure() {
        println!("{}", CS.t_strct::<AbstractIndex>(3, 8));

        let _ = Atom::Zero.simplify_metrics();
    }

    #[test]
    #[should_panic]
    fn test_color_matrix_element() {
        initialize();
        let _q = symbol!("spenso::Q";Real);
        let spin_sum_rule = parse!(
            "
            g(coad(Nc^2-1, left(3)), coad(Nc^2-1, right(3)))
                * g(coad(Nc^2-1, left(2)), coad(Nc^2-1, right(2)))
                * g(cof(Nc, right(0)), dind(cof(Nc, left(0))))
                * g(cof(Nc, left(1)), dind(cof(Nc, right(1))))",
            default_namespace = "spenso"
        );

        let amplitude_color = parse!(
            "
            t(coad(Nc^2-1, 6), cof(Nc, 5), dind(cof(Nc, 4)))
                * f(coad(Nc^2-1, 7), coad(Nc^2-1, 8), coad(Nc^2-1, 9))
                * g(coad(Nc^2-1, 2), coad(Nc^2-1, 7))
                * g(coad(Nc^2-1, 3), coad(Nc^2-1, 8))
                * g(coad(Nc^2-1, 6), coad(Nc^2-1, 9))
                * g(cof(Nc, 0), dind(cof(Nc, 5)))
                * g(cof(Nc, 4), dind(cof(Nc, 1)))",
            default_namespace = "spenso"
        );
        let amplitude_color_left = amplitude_color.wrap_indices(symbol!("spenso::left"));

        // return;
        let amplitude_color_right = amplitude_color
            .dirac_adjoint::<AbstractIndex>()
            .unwrap()
            .wrap_indices(symbol!("spenso::right"));
        println!("left{amplitude_color_left}");

        println!("right{amplitude_color_right}");
        let amp_squared_color = amplitude_color_left * spin_sum_rule * amplitude_color_right;
        let simplified_color = amp_squared_color.simplify_metrics().simplify_color();
        println!("simplified_color={}", simplified_color);

        let spin_sum_rule_src = parse_lit!(
            spenso::vbar(1, bis(D, left(1)))
                * spenso::v(1, bis(D, right(1)))
                * spenso::u(0, bis(D, left(0)))
                * spenso::ubar(0, bis(D, right(0)))
                * spenso::ϵbar(2, mink(D, left(2)))
                * spenso::ϵ(2, mink(D, right(2)))
                * spenso::ϵbar(3, mink(D, left(3)))
                * spenso::ϵ(3, mink(D, right(3))),
            default_namespace = "spenso"
        );

        let spin_sum_rule_trg = parse!(
            "
            1/4*1/(D-2)^2*
            (
                (-1) * gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))*Q(1,mink(D,1337))
                * gamma(bis(D,right(0)),bis(D,left(0)),mink(D,1338))*Q(0,mink(D,1338))
                * (-1) * g(mink(D,left(2)),mink(D,right(2)))
                * (-1) * g(mink(D,left(3)),mink(D,right(3)))
                )
                * (
                    g(coad(Nc^2-1, left(3)), coad(Nc^2-1, right(3)))
                    * g(coad(Nc^2-1, left(2)), coad(Nc^2-1, right(2)))
                    * g(cof(Nc, right(0)), dind(cof(Nc, left(0))))
                    * g(cof(Nc, left(1)), dind(cof(Nc, right(1))))
                )",
            default_namespace = "spenso"
        );

        let (amplitude, tgt) = colored_matrix_element();

        let amplitude_left = amplitude.wrap_indices(symbol!("spenso::left"));

        println!("Amplitude left:\n{}", amplitude_left.collect_factors());

        println!(
            "Amplitude left cooked:\n{}",
            amplitude_left.collect_factors().cook_indices()
        );
        let amplitude_right = amplitude.wrap_indices(symbol!("spenso::right"));

        println!("Amplitude right:\n{}", amplitude_right.factor());

        println!(
            "Amplitude right conj:\n{}",
            amplitude_right
                .dirac_adjoint::<AbstractIndex>()
                .unwrap()
                .factor()
        );

        let mut amp_squared =
            amplitude_left * amplitude_right.dirac_adjoint::<AbstractIndex>().unwrap();

        println!("Amplitude squared:\n{}", amp_squared.factor());

        // let dangling_atoms = amp_squared.list_dangling::<AbstractIndex>();

        // assert_eq!(dangling_atoms.len(), 8);

        amp_squared = amp_squared
            .expand()
            .replace(spin_sum_rule_src.to_pattern())
            .with(spin_sum_rule_trg.to_pattern());

        println!("Amplitude squared spin-summed:\n{}", amp_squared);

        let mut simplified_amp_squared = amp_squared.clone();

        simplified_amp_squared = simplified_amp_squared.simplify_color();

        println!(
            "Color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.simplify_gamma();

        println!(
            "Gamma+color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.to_dots();

        assert_eq!(
            tgt,
            simplified_amp_squared.factor(),
            "{}\nnot equal to\n{}",
            tgt,
            simplified_amp_squared.factor()
        );
    }

    #[test]
    #[should_panic]
    fn test_color_matrix_element_two() {
        initialize();

        let spin_sum_rule_src = parse_lit!(
            vbar(1, bis(D, left(1)))
                * v(1, bis(D, right(1)))
                * u(0, bis(D, left(0)))
                * ubar(0, bis(D, right(0)))
                * ϵbar(2, mink(D, left(2)))
                * ϵ(2, mink(D, right(2)))
                * ϵbar(3, mink(D, left(3)))
                * ϵ(3, mink(D, right(3))),
            default_namespace = "spenso"
        );

        let spin_sum_rule_trg = parse!(
            "
            1/4*1/(D-2)^2*
            (
                (-1) * gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))*Q(1,mink(D,1337))
                * gamma(bis(D,right(0)),bis(D,left(0)),mink(D,1338))*Q(0,mink(D,1338))
                * (-1) * g(mink(D,left(2)),mink(D,right(2)))
                * (-1) * g(mink(D,left(3)),mink(D,right(3)))
                )
                ",
            default_namespace = "spenso"
        );

        let (amplitude, tgt) = colored_matrix_element();

        let amplitude_left = amplitude.wrap_dummies::<AbstractIndex>(symbol!("spenso::left"));

        println!("Amplitude left:\n{}", amplitude_left.collect_factors());

        let amplitude_right = amplitude.wrap_dummies::<AbstractIndex>(symbol!("spenso::right"));

        println!("Amplitude right:\n{}", amplitude_right.conj().factor());

        let mut amp_squared = amplitude_left * amplitude_right.conj();

        println!("Amplitude squared:\n{}", amp_squared.factor());

        let _spin_sum_pat = parse!(
            "gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))",
            default_namespace = "spenso"
        )
        .to_pattern();

        let pols_pats: Vec<Pattern> = vec![
            function!(PS.vbar, RS.x__).into(),
            function!(PS.v, RS.x__).into(),
            function!(PS.u, RS.x__).into(),
            function!(PS.ubar, RS.x__).into(),
            function!(PS.eps, RS.x__).into(),
            function!(PS.ebar, RS.x__).into(),
        ];

        let mut pols_coefs = amp_squared.expand_in_patterns(&pols_pats);

        for (c, _) in &mut pols_coefs {
            println!("c:{c}");
            let r = c
                .replace(spin_sum_rule_src.to_pattern())
                .with(spin_sum_rule_trg.to_pattern());
            println!("r:{r}");
            *c = r;
        }

        amp_squared = pols_coefs.iter().fold(Atom::Zero, |a, (c, s)| a + c * s);

        println!("Amplitude squared spin-summed:\n{}", amp_squared.factor());

        let mut simplified_amp_squared = amp_squared.clone();

        simplified_amp_squared = simplified_amp_squared.simplify_color();

        println!(
            "Color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.simplify_gamma().simplify_metrics();

        println!(
            "Gamma+color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.to_dots();

        assert_eq!(
            tgt,
            simplified_amp_squared.factor(),
            "{:#}\nnot equal to\n{:#}",
            tgt,
            simplified_amp_squared.factor()
        );
    }

    #[test]
    fn test_val() {
        initialize();
        let expr = parse_lit!(
            (G ^ 3
                * (g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                    - g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9))))
                * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
                * g(mink(4, l(0)), mink(4, l(6)))
                * g(mink(4, l(1)), mink(4, l(7)))
                * g(mink(4, l(4)), mink(4, l(8)))
                * g(mink(4, l(5)), mink(4, l(9)))
                * g(bis(4, l(2)), bis(4, l(5)))
                * g(bis(4, l(3)), bis(4, l(6)))
                * g(dind(cof(Nc, l(6))), cof(Nc, 3))
                * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
                * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
                * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
                * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
                * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
                * f(
                    coad(Nc ^ 2 - 1, l(8)),
                    coad(Nc ^ 2 - 1, l(11)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * f(
                    coad(Nc ^ 2 - 1, l(9)),
                    coad(Nc ^ 2 - 1, l(10)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * ubar(2, bis(4, l(2)))
                * v(3, bis(4, l(3)))
                * ϵ(0, mink(4, l(0)))
                * ϵ(1, mink(4, l(1)))
                * ϵbar(4, mink(4, l(4)))
                + G
                ^ 3 * (g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                    - g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                    * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(dind(cof(Nc, l(6))), cof(Nc, 3))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
                    * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
                    * f(
                        coad(Nc ^ 2 - 1, l(8)),
                        coad(Nc ^ 2 - 1, l(10)),
                        coad(Nc ^ 2 - 1, l(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, l(9)),
                        coad(Nc ^ 2 - 1, l(11)),
                        coad(Nc ^ 2 - 1, l(12))
                    )
                    * ubar(2, bis(4, l(2)))
                    * v(3, bis(4, l(3)))
                    * ϵ(0, mink(4, l(0)))
                    * ϵ(1, mink(4, l(1)))
                    * ϵbar(4, mink(4, l(4)))
                    + G
                ^ 3 * (g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9)))
                    - g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                    * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(dind(cof(Nc, l(6))), cof(Nc, 3))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
                    * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
                    * f(
                        coad(Nc ^ 2 - 1, l(8)),
                        coad(Nc ^ 2 - 1, l(9)),
                        coad(Nc ^ 2 - 1, l(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, l(10)),
                        coad(Nc ^ 2 - 1, l(11)),
                        coad(Nc ^ 2 - 1, l(12))
                    )
                    * ubar(2, bis(4, l(2)))
                    * v(3, bis(4, l(3)))
                    * ϵ(0, mink(4, l(0)))
                    * ϵ(1, mink(4, l(1)))
                    * ϵbar(4, mink(4, l(4))))
                * (-G
                    ^ 3 * (g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                        - g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9))))
                        * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                        * g(mink(4, r(0)), mink(4, r(6)))
                        * g(mink(4, r(1)), mink(4, r(7)))
                        * g(mink(4, r(4)), mink(4, r(8)))
                        * g(mink(4, r(5)), mink(4, r(9)))
                        * g(bis(4, r(2)), bis(4, r(5)))
                        * g(bis(4, r(3)), bis(4, r(6)))
                        * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                        * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                        * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                        * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                        * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                        * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                        * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                        * f(
                            coad(Nc ^ 2 - 1, r(8)),
                            coad(Nc ^ 2 - 1, r(11)),
                            coad(Nc ^ 2 - 1, r(12))
                        )
                        * f(
                            coad(Nc ^ 2 - 1, r(9)),
                            coad(Nc ^ 2 - 1, r(10)),
                            coad(Nc ^ 2 - 1, r(12))
                        )
                        * u(2, bis(4, r(2)))
                        * vbar(3, bis(4, r(3)))
                        * ϵ(4, mink(4, r(4)))
                        * ϵbar(0, mink(4, r(0)))
                        * ϵbar(1, mink(4, r(1)))
                        - G
                    ^ 3 * (g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                        - g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                        * g(mink(4, r(0)), mink(4, r(6)))
                        * g(mink(4, r(1)), mink(4, r(7)))
                        * g(mink(4, r(4)), mink(4, r(8)))
                        * g(mink(4, r(5)), mink(4, r(9)))
                        * g(bis(4, r(2)), bis(4, r(5)))
                        * g(bis(4, r(3)), bis(4, r(6)))
                        * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                        * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                        * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                        * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                        * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                        * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                        * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                        * f(
                            coad(Nc ^ 2 - 1, r(8)),
                            coad(Nc ^ 2 - 1, r(10)),
                            coad(Nc ^ 2 - 1, r(12))
                        )
                        * f(
                            coad(Nc ^ 2 - 1, r(9)),
                            coad(Nc ^ 2 - 1, r(11)),
                            coad(Nc ^ 2 - 1, r(12))
                        )
                        * u(2, bis(4, r(2)))
                        * vbar(3, bis(4, r(3)))
                        * ϵ(4, mink(4, r(4)))
                        * ϵbar(0, mink(4, r(0)))
                        * ϵbar(1, mink(4, r(1)))
                        - G
                    ^ 3 * (g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9)))
                        - g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                        * g(mink(4, r(0)), mink(4, r(6)))
                        * g(mink(4, r(1)), mink(4, r(7)))
                        * g(mink(4, r(4)), mink(4, r(8)))
                        * g(mink(4, r(5)), mink(4, r(9)))
                        * g(bis(4, r(2)), bis(4, r(5)))
                        * g(bis(4, r(3)), bis(4, r(6)))
                        * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                        * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                        * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                        * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                        * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                        * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                        * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                        * f(
                            coad(Nc ^ 2 - 1, r(8)),
                            coad(Nc ^ 2 - 1, r(9)),
                            coad(Nc ^ 2 - 1, r(12))
                        )
                        * f(
                            coad(Nc ^ 2 - 1, r(10)),
                            coad(Nc ^ 2 - 1, r(11)),
                            coad(Nc ^ 2 - 1, r(12))
                        )
                        * u(2, bis(4, r(2)))
                        * vbar(3, bis(4, r(3)))
                        * ϵ(4, mink(4, r(4)))
                        * ϵbar(0, mink(4, r(0)))
                        * ϵbar(1, mink(4, r(1)))),
            default_namespace = "spenso"
        );

        println!("{expr}");
        println!("Simplify_metrics");
        println!("{}", expr.clone().simplify_metrics());
        println!("Simplify_color");
        println!(
            "{:>}",
            expr.simplify_gamma()
                .simplify_color()
                .expand()
                .simplify_metrics()
                .to_dots()
        );
    }

    #[test]
    fn ratio_simplify() {
        test_initialize();
        let expr = parse_lit!(
            G ^ 4 * ee
                ^ 2 * f(
                    coad(ohoho, dummy(0)),
                    coad(ohoho, dummy(1)),
                    coad(ohoho, dummy(2))
                ) * t(
                    coad(ohoho, dummy(0)),
                    cof(ahaha, dummy(3)),
                    dind(cof(ahaha, dummy(4)))
                ) * t(
                    coad(ohoho, dummy(1)),
                    cof(ahaha, dummy(5)),
                    dind(cof(ahaha, dummy(3)))
                ) * t(
                    coad(ohoho, dummy(2)),
                    cof(ahaha, dummy(4)),
                    dind(cof(ahaha, dummy(5)))
                ),
            default_namespace = "spenso"
        );

        insta::assert_snapshot!(
            expr.cook_indices()
                .simplify_color()
                .to_bare_ordered_string(),@"-1𝑖*G^4*TR^2*ahaha*ee^2*ohoho"
        );
    }

    #[test]
    fn simple() {
        test_initialize();
        // tc₄ г₁₀ г₈·tc₆ г₈ г₁₀·fc₄ c₀ c₂·fc₆ c₂ c₁
        let expr = parse_lit!(
            f(
                coad(Nc ^ 2 - 1, c4),
                coad(Nc ^ 2 - 1, c0),
                coad(Nc ^ 2 - 1, c2)
            ) * f(
                coad(Nc ^ 2 - 1, c6),
                coad(Nc ^ 2 - 1, c2),
                coad(Nc ^ 2 - 1, c0)
            ) * t(coad(Nc ^ 2 - 1, c4), cof(Nc, i0), dind(cof(Nc, i2)))
                * t(coad(Nc ^ 2 - 1, c6), cof(Nc, i2), dind(cof(Nc, i0))),
            default_namespace = "spenso"
        );

        println!("{}", expr.simplify_color());
        println!(
            "{}",
            expr.simplify_color()
                .replace(parse_lit!(spenso::Nc))
                .with(3)
        );
    }

    #[test]
    fn minus_sign() {
        test_initialize();
        let expr1 = parse_lit!(
            f(coad(8, hedge(5)), coad(8, hedge(9)), coad(8, hedge(13)))
                * g(coad(8, hedge(12)), coad(8, hedge(13)))
                * g(coad(8, hedge(4)), coad(8, hedge(5)))
                * g(coad(8, hedge(8)), coad(8, hedge(9)))
                * g(cof(3, hedge(10)), dind(cof(3, hedge(11))))
                * g(cof(3, hedge(11)), dind(cof(3, hedge(15))))
                * g(cof(3, hedge(15)), dind(cof(3, hedge(14))))
                * g(cof(3, hedge(16)), dind(cof(3, hedge(17))))
                * g(cof(3, hedge(17)), dind(cof(3, hedge(7))))
                * g(cof(3, hedge(2)), dind(cof(3, hedge(3))))
                * g(cof(3, hedge(7)), dind(cof(3, hedge(6))))
                * t(
                    coad(8, hedge(12)),
                    cof(3, hedge(14)),
                    dind(cof(3, hedge(16)))
                )
                * t(coad(8, hedge(4)), cof(3, hedge(6)), dind(cof(3, hedge(2))))
                * t(coad(8, hedge(8)), cof(3, hedge(3)), dind(cof(3, hedge(10)))),
            default_namespace = "spenso"
        );

        let expr2 = parse_lit!(
            f(coad(8, hedge(5)), coad(8, hedge(9)), coad(8, hedge(13)))
                * g(coad(8, hedge(12)), coad(8, hedge(13)))
                * g(coad(8, hedge(4)), coad(8, hedge(5)))
                * g(coad(8, hedge(8)), coad(8, hedge(9)))
                * g(cof(3, hedge(11)), dind(cof(3, hedge(10))))
                * g(cof(3, hedge(14)), dind(cof(3, hedge(15))))
                * g(cof(3, hedge(15)), dind(cof(3, hedge(11))))
                * g(cof(3, hedge(17)), dind(cof(3, hedge(16))))
                * g(cof(3, hedge(3)), dind(cof(3, hedge(2))))
                * g(cof(3, hedge(6)), dind(cof(3, hedge(7))))
                * g(cof(3, hedge(7)), dind(cof(3, hedge(17))))
                * t(
                    coad(8, hedge(12)),
                    cof(3, hedge(16)),
                    dind(cof(3, hedge(14)))
                )
                * t(coad(8, hedge(4)), cof(3, hedge(2)), dind(cof(3, hedge(6))))
                * t(coad(8, hedge(8)), cof(3, hedge(10)), dind(cof(3, hedge(3)))),
            default_namespace = "spenso"
        );

        // println!(
        //     "{}",
        //     (expr1.cook_indices().canonize(AbstractIndex::Dummy)
        //         / (expr2.cook_indices().canonize(AbstractIndex::Dummy)))
        //     .cancel()
        // );
        println!(
            "{}\n",
            expr1
                .simplify_metrics()
                .cook_indices()
                .canonize(AbstractIndex::Dummy)
        );
        println!(
            "{}\n",
            expr2
                .simplify_metrics()
                .cook_indices()
                .canonize(AbstractIndex::Dummy)
        );
        assert!((expr1.simplify_color() + expr2.simplify_color()).is_zero());
    }
}
