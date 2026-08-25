use std::ops::Deref;

use color_eyre::eyre::{bail, ensure};
use idenso::{
    color::{CS, ColorSimplifier},
    representations::{ColorAdjoint, ColorFundamental},
};
use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::representation::{Minkowski, RepName},
};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol},
    function,
    poly::series::SeriesDepth,
};

use crate::utils::{GS, TENSORLIB, W_};

use super::ParsingNet;
pub type ParsingNetError = spenso::network::TensorNetworkError<
    spenso::structure::IndexlessNamedStructure<
        Symbol,
        Vec<Atom>,
        spenso::structure::representation::LibraryRep,
        super::aind::Aind,
    >,
    Symbol,
>;

pub trait NumeratorAtomExt {
    /// Truncate through an absolute integer order, preserving an existing root
    /// product's variable-independent factors outside the coefficient sum.
    /// Other roots retain native series behavior, including additive zeros.
    /// Regroup the result in the caller's owned, Taylor-independent numerator
    /// coefficient keys. These must be symbols or functions occurring linearly,
    /// outside other functions and inverses. Their definitions and all dependent
    /// energy arguments remain the caller's responsibility; this method never
    /// makes a Taylor-dependent numerator opaque.
    fn series_preserving_factors(
        &self,
        variable: Symbol,
        expansion_point: AtomView<'_>,
        depth: i64,
        numerator_family_keys: &[Atom],
    ) -> color_eyre::Result<Atom>;

    fn to_param_color(&self) -> Atom;
    // fn wrap_color(&self, symbol: Symbol) -> Atom;
    fn kill_color(&self) -> Atom;

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom;

    fn unwrap_function(&self, symbol: Symbol) -> Atom;

    #[allow(clippy::result_large_err)]
    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError>;
}

impl NumeratorAtomExt for Atom {
    fn series_preserving_factors(
        &self,
        variable: Symbol,
        expansion_point: AtomView<'_>,
        depth: i64,
        numerator_family_keys: &[Atom],
    ) -> color_eyre::Result<Atom> {
        self.as_view().series_preserving_factors(
            variable,
            expansion_point,
            depth,
            numerator_family_keys,
        )
    }

    fn to_param_color(&self) -> Atom {
        self.as_view().to_param_color()
    }
    fn kill_color(&self) -> Atom {
        self.wrap_color(GS.killing_func)
    }

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.as_view().map_mink_dim(dim)
    }
    // fn wrap_color(&self, symbol: Symbol) -> Atom {
    //     self.as_view().wrap_color(symbol)
    // }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.as_view().unwrap_function(symbol)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        self.as_view().parse_into_net()
    }
}

impl NumeratorAtomExt for AtomView<'_> {
    fn series_preserving_factors(
        &self,
        variable: Symbol,
        expansion_point: AtomView<'_>,
        depth: i64,
        numerator_family_keys: &[Atom],
    ) -> color_eyre::Result<Atom> {
        for (index, key) in numerator_family_keys.iter().enumerate() {
            ensure!(
                matches!(key.as_view(), AtomView::Var(_) | AtomView::Fun(_))
                    && !key.contains_symbol(variable),
                "numerator-family key must be a Taylor-independent symbol or function: {key}"
            );
            ensure!(
                numerator_family_keys[..index]
                    .iter()
                    .all(|other| !key.contains(other) && !other.contains(key)),
                "numerator-family keys must be distinct and cannot contain one another: {key}"
            );
        }
        let contains_family =
            |atom: AtomView<'_>| numerator_family_keys.iter().any(|key| atom.contains(key));
        let mut independent = Atom::num(1);
        let mut dependent = self.to_owned();
        if let AtomView::Mul(product) = self {
            dependent = Atom::num(1);
            for factor in product.iter() {
                if factor.contains_symbol(variable) || contains_family(factor) {
                    dependent *= factor;
                } else {
                    independent *= factor;
                }
            }
        }
        // Keep native whole-product precision lifting and cross-arm zero
        // detection before regrouping the surviving coefficient families.
        let series = dependent
            .series(variable, expansion_point, SeriesDepth::absolute(depth))?
            .to_atom();
        if numerator_family_keys.is_empty() {
            return Ok(independent * series);
        }
        // Reject nonlinear or hidden occurrences before polynomial collection
        // can distribute powers or products of coefficient families.
        let mut invalid = None;
        series.visitor(&mut |atom| {
            if invalid.is_some()
                || numerator_family_keys
                    .iter()
                    .any(|key| key.as_view() == atom)
            {
                return false;
            }
            let valid = match atom {
                AtomView::Mul(product) => {
                    product
                        .iter()
                        .filter(|factor| contains_family(*factor))
                        .count()
                        <= 1
                }
                AtomView::Fun(_) | AtomView::Pow(_) => !contains_family(atom),
                _ => true,
            };
            if !valid {
                invalid = Some(atom.to_owned());
            }
            valid
        });
        if let Some(invalid) = invalid {
            bail!(
                "numerator-family keys must occur linearly outside functions and inverses: {invalid}"
            );
        }
        let mut grouped = Atom::Zero;
        for (key, coefficient) in series.coefficient_list::<u32>(numerator_family_keys) {
            ensure!(
                key.is_one() || numerator_family_keys.contains(&key),
                "series must be linear in its numerator-family keys, found {key}"
            );
            ensure!(
                !contains_family(coefficient.as_view()),
                "numerator-family key remains hidden in series coefficient {coefficient}"
            );
            grouped += key * coefficient;
        }
        Ok(independent * grouped)
    }

    fn kill_color(&self) -> Atom {
        self.wrap_color(GS.killing_func)
    }

    fn to_param_color(&self) -> Atom {
        let adj = ColorAdjoint {};
        let fund = ColorFundamental {};
        self.replace(adj.to_symbolic([W_.d_, W_.a_]))
            .with(adj.to_symbolic([CS.nc * CS.nc - 1, Atom::var(W_.a_)]))
            .replace(fund.to_symbolic([W_.d_, W_.a_]))
            .with(fund.to_symbolic([CS.nc, W_.a_]))
    }
    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.replace(Minkowski {}.to_symbolic([W_.d_, W_.a___]))
            .with(Minkowski {}.to_symbolic([dim.into().into_owned(), Atom::var(W_.a___)]))
    }
    // fn wrap_color(&self, symbol: Symbol) -> Atom {
    //     self.expand_color()
    //         .into_iter()
    //         .fold(Atom::Zero, |a, (c, s)| a + function!(symbol, c) * s)
    // }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.replace(function!(symbol, W_.a___)).with(W_.a___)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        ParsingNet::try_from_view(
            *self,
            TENSORLIB.read().unwrap().deref(),
            &ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: false,
                        expand_schoonship: true,
                        expand_inside_chains: true,
                    },
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            },
        )
    }

    // fn parse_into_only_lib_net<T: TensorLibraryData + Clone + Default>(
    //     &self,
    //     one: T,
    //     zero: T,
    // ) -> Result<ParsingNet, ParsingNetError> {
    //     let mut lib = hep_lib(one, zero);

    //     ParsingNet::try_from_view(*self, &lib);
    // }
}

#[cfg(test)]
mod tests {
    use idenso::IndexTooling;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{
            representation::{Minkowski, RepName},
            slot::{DummyAind, IsAbstractSlot, Slot},
        },
    };
    use symbolica::{
        atom::{Atom, AtomCore},
        function, parse_lit, symbol,
    };

    use crate::{
        dot,
        graph::{FeynmanGraph, Graph, parse::IntoGraph},
        initialisation::test_initialise,
        numerator::aind::Aind,
        utils::GS,
        uv::UltravioletGraph,
    };

    use super::NumeratorAtomExt;

    #[test]
    fn series_preserves_independent_tensor_factors() {
        test_initialise().unwrap();
        let t = symbol!("series_t");
        let slot: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());
        let tensor = GS.emr_vec_index(EdgeIndex(0), slot.to_atom())
            * GS.emr_vec_index(EdgeIndex(1), slot.to_atom());
        for spectator in [parse_lit!(g), tensor * parse_lit!(u + v)] {
            let source = &spectator * parse_lit!((a + b * series_t) * (c + d * series_t));
            let result = source
                .series_preserving_factors(t, Atom::Zero.as_view(), 1, &[])
                .unwrap();
            assert_eq!(
                result,
                &spectator * parse_lit!(a * c + (a * d + b * c) * series_t)
            );
            assert_eq!(
                result
                    .series_preserving_factors(t, Atom::Zero.as_view(), 1, &[])
                    .unwrap(),
                result
            );
        }
    }

    #[test]
    fn series_preserves_laurent_order_with_mapped_numerator() {
        let t = symbol!("series_t");
        let source = parse_lit!(g * (a + b * series_t) * (c + d * series_t) / series_t ^ 2);
        for (depth, expected) in [
            (
                -1,
                parse_lit!(g * (a * c / series_t ^ 2 + (a * d + b * c) / series_t)),
            ),
            (
                0,
                parse_lit!(g * (a * c / series_t ^ 2 + (a * d + b * c) / series_t + b * d)),
            ),
        ] {
            assert_eq!(
                source
                    .series_preserving_factors(t, Atom::Zero.as_view(), depth, &[])
                    .unwrap(),
                expected
            );
        }
    }

    #[test]
    fn series_preserves_native_additive_zero() {
        let t = symbol!("series_t");
        for source in [
            Atom::Zero,
            parse_lit!(g * (a + b * series_t) / series_t - g * a / series_t - g * b),
        ] {
            let result = source
                .series_preserving_factors(t, Atom::Zero.as_view(), 0, &[])
                .unwrap();
            assert!(result.is_zero());
            assert_eq!(result, source.series(t, Atom::Zero, 0).unwrap().to_atom());
            assert!(result.replace(t).with(Atom::num(1)).is_zero());
        }
    }

    #[test]
    fn series_preserves_nonzero_expansion_point() {
        let t = symbol!("series_t");
        let source = parse_lit!(g * (a + b * series_t) * (c + d * series_t));
        let result = source
            .series_preserving_factors(t, Atom::num(2).as_view(), 1, &[])
            .unwrap();
        assert_eq!(
            result,
            parse_lit!(
                g * ((a + 2 * b) * (c + 2 * d)
                    + (series_t - 2) * (b * (c + 2 * d) + d * (a + 2 * b)))
            )
        );
    }

    #[test]
    fn series_keeps_native_additive_denominators() {
        let t = symbol!("series_t");
        let source =
            parse_lit!(s1 * g * (a + b * series_t) / D1 + s2 * g * (c + d * series_t) / D2);
        for depth in [0, 1] {
            assert_eq!(
                source
                    .series_preserving_factors(t, Atom::Zero.as_view(), depth, &[])
                    .unwrap(),
                source.series(t, Atom::Zero, depth).unwrap().to_atom()
            );
        }
    }

    #[test]
    fn series_regroups_owned_numerator_coefficient_families() {
        let t = symbol!("series_t");
        let family = symbol!("series_owned_coefficient"; Scalar);
        let keys = [0, 1, 2].map(|order| function!(family, order));
        let [n0, n1, n2] = &keys;
        let spectator = parse_lit!(g * (u + v));
        let source = &spectator
            * (n0 + n1 * Atom::var(t) + n2 * Atom::var(t).pow(2))
            * parse_lit!((a + b * series_t) / series_t ^ 2);
        for (depth, expected) in [
            (
                -1,
                n0 * parse_lit!(a / series_t ^ 2 + b / series_t) + n1 * parse_lit!(a / series_t),
            ),
            (
                0,
                n0 * parse_lit!(a / series_t ^ 2 + b / series_t)
                    + n1 * parse_lit!(a / series_t + b)
                    + n2 * parse_lit!(a),
            ),
        ] {
            let result = source
                .series_preserving_factors(t, Atom::Zero.as_view(), depth, &keys)
                .unwrap();
            assert_eq!(result, &spectator * expected);
            assert_eq!(
                result
                    .series_preserving_factors(t, Atom::Zero.as_view(), depth, &keys)
                    .unwrap(),
                result
            );
            let native = source.series(t, Atom::Zero, depth).unwrap().to_atom();
            assert!((result - native).expand().is_zero());
        }
    }

    #[test]
    fn series_family_regrouping_preserves_zero_and_nonzero_center() {
        let t = symbol!("series_t");
        let family = symbol!("series_owned_coefficient"; Scalar);
        let keys = [0, 1, 2].map(|order| function!(family, order));
        let [n0, n1, n2] = &keys;
        let source = n0 * parse_lit!((a + b * series_t) / series_t)
            - n0 * parse_lit!(a / series_t)
            - n0 * parse_lit!(b);
        for zero in [Atom::Zero, source] {
            let result = zero
                .series_preserving_factors(t, Atom::Zero.as_view(), 0, &keys)
                .unwrap();
            assert!(result.is_zero());
            assert!(result.replace(t).with(Atom::num(1)).is_zero());
        }

        let shift = parse_lit!(series_t - 2);
        let source = (n0 + n1 * &shift + n2 * shift.pow(2))
            * (parse_lit!(a) + parse_lit!(b) * &shift)
            / shift.pow(2);
        let result = source
            .series_preserving_factors(t, Atom::num(2).as_view(), 0, &keys)
            .unwrap();
        assert_eq!(
            result,
            n0 * (parse_lit!(a) / shift.pow(2) + parse_lit!(b) / &shift)
                + n1 * (parse_lit!(a) / &shift + parse_lit!(b))
                + n2 * parse_lit!(a)
        );
        assert!(
            (result - source.series(t, Atom::num(2), 0).unwrap().to_atom())
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn series_family_specialization_preserves_ose_derivatives() {
        let t = symbol!("series_t");
        let family = symbol!("series_ose_coefficient"; Scalar);
        let keys = [0, 1, 2].map(|order| {
            function!(
                family,
                order,
                parse_lit!(sigma),
                parse_lit!(z),
                parse_lit!(tau),
                parse_lit!(w)
            )
        });
        // Scalar diagnostic: each energy retains its own sign, integer node
        // and fixed shift. Nonconstant OSEs must contribute their derivatives.
        let q = parse_lit!(sigma * (1 + series_t) ^ (1 / 2) + z * M + x);
        let r = parse_lit!(tau * (4 + 2 * series_t) ^ (1 / 2) + w * M + y);
        let q0 = parse_lit!(sigma + z * M + x);
        let r0 = parse_lit!(2 * tau + w * M + y);
        let coefficients = [
            &q0 * &r0,
            &q0 * parse_lit!(tau / 2) + parse_lit!(sigma / 2) * &r0,
            -&q0 * parse_lit!(tau / 16) + parse_lit!(sigma * tau / 4) - parse_lit!(sigma / 8) * &r0,
        ];
        let truncated_numerator =
            &keys[0] + &keys[1] * Atom::var(t) + &keys[2] * Atom::var(t).pow(2);
        let denominator = parse_lit!(series_t ^ 2 * (1 - series_t));
        let retained = (truncated_numerator / &denominator)
            .series_preserving_factors(t, Atom::Zero.as_view(), 0, &keys)
            .unwrap();
        assert_eq!(
            retained,
            &keys[0] * parse_lit!(1 / series_t ^ 2 + 1 / series_t + 1)
                + &keys[1] * parse_lit!(1 / series_t + 1)
                + &keys[2]
        );
        let rebuilt = keys
            .iter()
            .zip(&coefficients)
            .fold(retained, |atom, (key, body)| {
                atom.replace(key.to_pattern()).with(body.to_pattern())
            });
        let source = q * r / denominator;
        let native = source.series(t, Atom::Zero, 0).unwrap().to_atom();
        assert!((&rebuilt - native).expand().is_zero());
        // The first row makes the first energy vanish at the expansion point.
        // Other rows keep the two repeated-occurrence nodes independent.
        for values in [
            [1, -1, 1, 2, 1, 0, 0],
            [-1, 2, 1, -3, 2, 3, -1],
            [1, 0, -1, 2, -2, 1, 4],
        ] {
            let variables = ["sigma", "z", "tau", "w", "M", "x", "y"].map(|name| symbol!(name));
            let (specialized, rebuilt) = variables.into_iter().zip(values).fold(
                (source.clone(), rebuilt.clone()),
                |(source, result), (variable, value)| {
                    (
                        source.replace(variable).with(Atom::num(value)),
                        result.replace(variable).with(Atom::num(value)),
                    )
                },
            );
            assert!(
                (rebuilt - specialized.series(t, Atom::Zero, 0).unwrap().to_atom())
                    .expand()
                    .is_zero()
            );
        }
    }

    #[test]
    fn series_rejects_invalid_numerator_families() {
        let t = symbol!("series_t");
        let key = parse_lit!(n0);
        for keys in [
            vec![parse_lit!(n0 + n1)],
            vec![parse_lit!(n0(series_t))],
            vec![key.clone(), key.clone()],
            vec![key.clone(), parse_lit!(outer(n0))],
        ] {
            assert!(
                key.series_preserving_factors(t, Atom::Zero.as_view(), 0, &keys)
                    .is_err()
            );
        }
        for source in [
            parse_lit!(n0 ^ 2),
            parse_lit!((n0 + 1) ^ 1000),
            parse_lit!(1 / n0),
            parse_lit!(outer(n0)),
        ] {
            assert!(
                source
                    .series_preserving_factors(
                        t,
                        Atom::Zero.as_view(),
                        0,
                        std::slice::from_ref(&key),
                    )
                    .is_err()
            );
        }
    }

    #[test]
    fn series_with_families_preserves_native_errors() {
        let t = symbol!("series_t");
        let source = parse_lit!(n0 * singular_argument(1 / series_t));
        let native = source.series(t, Atom::Zero, 0).unwrap_err();
        let error = source
            .series_preserving_factors(t, Atom::Zero.as_view(), 0, &[parse_lit!(n0)])
            .unwrap_err();
        assert_eq!(
            error.downcast_ref::<symbolica::poly::series::SeriesError>(),
            Some(&native)
        );
    }

    #[test]
    fn dummy_parsing() {
        test_initialise().unwrap();

        let e_mass = parse_lit!(M_e);

        let m2 = &e_mass * &e_mass;

        let mink: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());

        let e = EdgeIndex(0);

        let sqrt = symbol!("sqrt_scalar", tag = SPENSO_TAG.broadcast);

        let a = function!(
            sqrt,
            (GS.emr_vec_index(e, mink.to_atom()) * GS.emr_vec_index(e, mink.to_atom()) + m2)
                .pow(Atom::num(2))
        );

        let net = a.parse_into_net().unwrap();

        println!("{}", net.dot_pretty())
    }

    #[test]
    fn canonize_color() {
        test_initialise().unwrap();
        let gls: Vec<Graph> = dot!(
            digraph{
            num = "1";

            ext0 [style=invis];
            2:0-> ext0 [id=0 dir=none is_cut=0  particle="a"];
            ext1 [style=invis];
            ext1-> 3:1 [id=1 dir=none is_cut=0  particle="a"];
            0:2-> 1:3 [id=2   particle="d"];
            0:4-> 1:5 [id=3 dir=none   particle="g"];
            3:6-> 0:7 [id=4   particle="d"];
            1:8-> 2:9 [id=5   particle="d"];
            2:10-> 3:11 [id=6   particle="d"];
        }

        digraph GL8{
            num = 1;
        0[int_id=V_74];
        1[int_id=V_74];
        2[int_id=V_71];
        3[int_id=V_71];
        ext0 [style=invis];
        2:0-> ext0 [id=0 dir=none is_cut=0  particle=a];
        ext1 [style=invis];
        ext1-> 3:1 [id=1 dir=none is_cut=0  particle=a];
        0:2-> 1:3 [id=2   particle=d];
        0:4-> 1:5 [id=3 dir=none   particle=g];
        0:6-> 3:7 [id=4 dir=back   particle="d~"];
        1:8-> 2:9 [id=5   particle=d];
        2:10-> 3:11 [id=6   particle=d];
        }

        )
        .unwrap();

        for g in gls {
            let mut numerator = g.numerator(&g.no_dummy(), &g.empty_subgraph());

            // TODO Check if we include overall factor in main
            numerator.state.expr *= &g.global_prefactor.num * &g.global_prefactor.projector; // * &gl5.overall_factor;
            // numerator.state.expr = numerator.state.expr.replace_multiple(&cpl_reps);

            let numerator_color_simplified = numerator
                .clone()
                .color_simplify()
                .get_single_atom()
                .unwrap()
                .canonize(Aind::Dummy)
                .expect("test expression should canonicalize");

            println!("numerator_color_simplified:{numerator_color_simplified}");
            println!("numerator:{}", numerator.state.expr);
        }
    }

    #[test]
    fn canonizations() {
        test_initialise().unwrap();

        let a = parse_lit!(
            ((-2 * spenso::projp(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(2))
            ) + spenso::projm(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(2))
            )) * -1𝑖
                / 6
                * UFO::sw
                ^ 2 + -1𝑖 / 2 * UFO::cw
                ^ 2 * spenso::projm(
                    spenso::bis(4, gammalooprs::edge(0)),
                    spenso::bis(4, gammalooprs::hedge(2))
                ))
                * ((-2
                    * spenso::projp(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    )
                    + spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    ))
                    * -1𝑖
                    / 6
                    * UFO::sw
                    ^ 2 + -1𝑖 / 2 * UFO::cw
                    ^ 2 * spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    ))
                * (-1 * UFO::MZ
                    ^ 2 * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    ) + gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(4)))
                        * gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(5))))
                * (-1 * gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                * 1
                / 3
                * UFO::MZ
                ^ (-2) * UFO::cw
                ^ (-2) * UFO::ee
                ^ 4 * UFO::sw
                ^ (-2)
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::e(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ebar(0, spenso::mink(4, gammalooprs::hedge(0)))
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::mink(4, gammalooprs::hedge(1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::mink(4, gammalooprs::edge(2, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::bis(4, gammalooprs::edge(0)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::mink(4, gammalooprs::edge(4, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::mink(4, gammalooprs::edge(1, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::hedge(0))
                    )
        );
        println!("a:{}", a);
        let b = parse_lit!(
            ((-2 * spenso::projp(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(6))
            ) + spenso::projm(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(6))
            )) * -1𝑖
                / 6
                * UFO::sw
                ^ 2 + -1𝑖 / 2 * UFO::cw
                ^ 2 * spenso::projm(
                    spenso::bis(4, gammalooprs::edge(0)),
                    spenso::bis(4, gammalooprs::hedge(6))
                ))
                * ((-2
                    * spenso::projp(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    )
                    + spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    ))
                    * -1𝑖
                    / 6
                    * UFO::sw
                    ^ 2 + -1𝑖 / 2 * UFO::cw
                    ^ 2 * spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    ))
                * (-1 * UFO::MZ
                    ^ 2 * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    ) + gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(4)))
                        * gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(5))))
                * (-1 * gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                * 1
                / 3
                * UFO::MZ
                ^ (-2) * UFO::cw
                ^ (-2) * UFO::ee
                ^ 4 * UFO::sw
                ^ (-2)
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::e(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ebar(0, spenso::mink(4, gammalooprs::hedge(0)))
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::mink(4, gammalooprs::hedge(0))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::bis(4, gammalooprs::edge(0)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::mink(4, gammalooprs::edge(2, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::mink(4, gammalooprs::edge(4, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::hedge(1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::mink(4, gammalooprs::edge(1, 1))
                    )
        );

        println!("b:{}", b);

        println!("ratio:{}", &a / &b);

        let ac = a
            .canonize(Aind::Dummy)
            .expect("test expression should canonicalize");
        let bc = b
            .canonize(Aind::Dummy)
            .expect("test expression should canonicalize");
        println!("ac:{}", ac);
        println!("bc:{}", bc);
        println!("ratio canonized:{}", ac / bc);
    }

    // #[test]
    // fn test_can() {
    //     let a = parse_lit!(T(a, b, c) * T(c, d, e) * T(d, b, f)(K(e) + P(e)) * (K(f) + P(f)));
    //     let b = parse_lit!(T(a, d) * T(d, c));

    //     let indices = vec![
    //         (parse_lit!(a), 1),
    //         (parse_lit!(a), 1),
    //         (parse_lit!(a), 1),
    //         (parse_lit!(b), 1),
    //         (parse_lit!(b), 1),
    //         (parse_lit!(c), 1),
    //         (parse_lit!(d), 1),
    //         (parse_lit!(e), 1),
    //         (parse_lit!(f), 1),
    //     ];

    //     let ac = a.canonize_tensors(&indices);
    //     println!("{}", ac.unwrap());
    //     let bc = b.canonize_tensors(&indices);
    //     println!("{}", bc.unwrap());
    // }
}
