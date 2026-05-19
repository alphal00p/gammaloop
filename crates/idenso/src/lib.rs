#![allow(uncommon_codepoints)]

use std::collections::BTreeMap;

use eyre::eyre;
use linnet::half_edge::subgraph::{BaseSubgraph, ModifySubSet, SuBitGraph, SubSetLike};
use shorthands::metric::{list_dangling_impl, wrap_dummies_impl, wrap_indices_impl};
use spenso::{
    network::{graph::NetworkEdge, library::function_lib::INBUILTS, parsing::ParseSettings},
    structure::{
        HasName, TensorStructure,
        representation::RepName,
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    symbol,
};
use thiserror::Error;

use crate::{
    dirac::{AGS, GammaSimplifier},
    rep_symbols::RS,
    representations::Bispinor,
    shorthands::metric::MetricSimplifier,
    tensor::{SymbolicNetExt, SymbolicNetParse},
};

pub mod tensor;

pub mod color;
pub mod cook;
pub mod dirac;
pub mod epsilon;
// pub mod parsing_ind;
#[cfg(feature = "python")]
pub mod python;
#[cfg(any(test, feature = "reference-cases"))]
pub mod reference_cases;
pub mod rep_symbols;
pub mod representations;
pub mod selective_expand;
pub mod shorthands;
#[cfg(test)]
pub(crate) mod test_support;

pub use cook::{
    CookMode, CookOutputTags, CookSettings, CookSourceFilter, CookTagFilter, Cookable, CookingError,
};

#[macro_export]
macro_rules!  symbol_set {
        // Identifier symbols with an explicit namespace.
        ($struct_name:ident, $static_name:ident, namespace = $namespace:literal; $($char:ident)*) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $char: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($char: ::symbolica::symbol!(concat!($namespace, "::", stringify!($char))),)*
            });
        };

        // Identifier symbols with explicit struct and static names
        ($struct_name:ident, $static_name:ident; $($char:ident)*) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $char: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($char: ::symbolica::symbol!(stringify!($char)),)*
            });
        };

        // String literals as individual statics with an explicit namespace.
        (statics, namespace = $namespace:literal; $($field:ident : $string:literal),* $(,)?) => {
            $(
                #[allow(non_upper_case_globals)]
                pub static $field: ::std::sync::LazyLock<::symbolica::atom::Symbol> =
                    ::std::sync::LazyLock::new(|| ::symbolica::symbol!(concat!($namespace, "::", $string)));
            )*
        };

        // String literals as individual statics
        (statics; $($field:ident : $string:literal),* $(,)?) => {
            $(
                #[allow(non_upper_case_globals)]
                pub static $field: ::std::sync::LazyLock<::symbolica::atom::Symbol> =
                    ::std::sync::LazyLock::new(|| ::symbolica::symbol!($string));
            )*
        };

        // String literals grouped in a struct with an explicit namespace.
        ($struct_name:ident, $static_name:ident, namespace = $namespace:literal; $($field:ident : $string:literal),* $(,)?) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $field: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($field: ::symbolica::symbol!(concat!($namespace, "::", $string)),)*
            });
        };

        // String literals grouped in a struct
        ($struct_name:ident, $static_name:ident; $($field:ident : $string:literal),* $(,)?) => {
            #[allow(dead_code, non_snake_case)]
            pub struct $struct_name {
                $(pub $field: ::symbolica::atom::Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($field: ::symbolica::symbol!($string),)*
            });
        };
    }

symbol_set!(Wildcards, W_;
    a_ b_ c_ d_ e_ f_ g_ h_ i_ j_ k_ l_ m_ n_ o_ p_ q_ r_ s_ t_ u_ v_ w_ x_ y_ z_
    a__ b__ c__ d__ e__ f__ g__ h__ i__ j__ k__ l__ m__ n__ o__ p__ q__ r__ s__ t__ u__ v__ w__ x__ y__ z__
    a___ b___ c___ d___ e___ f___ g___ h___ i___ j___ k___ l___ m___ n___ o___ p___ q___ r___ s___ t___ u___ v___ w___ x___ y___ z___
);

/// Builds a stripped or indexed bispinor representation atom.
#[macro_export]
macro_rules! bis {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::Bispinor {}; $($args)*)
    };
}

/// Builds a stripped or indexed spin-fundamental representation atom.
#[macro_export]
macro_rules! spf {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::SpinFundamental {}; $($args)*)
    };
}

/// Builds a stripped or indexed spin-antifundamental representation atom.
#[macro_export]
macro_rules! spaf {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::SpinAntiFundamental {}; $($args)*)
    };
}

/// Builds a stripped or indexed color-fundamental representation atom.
#[macro_export]
macro_rules! cof {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::ColorFundamental {}; $($args)*)
    };
}

/// Builds a stripped or indexed color-antifundamental representation atom.
#[macro_export]
macro_rules! coaf {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::ColorAntiFundamental {}; $($args)*)
    };
}

/// Builds a stripped or indexed color-sextet representation atom.
#[macro_export]
macro_rules! cos {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::ColorSextet {}; $($args)*)
    };
}

/// Builds a stripped or indexed color-antisextet representation atom.
#[macro_export]
macro_rules! coas {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::ColorAntiSextet {}; $($args)*)
    };
}

/// Builds a stripped or indexed color-adjoint representation atom.
#[macro_export]
macro_rules! coad {
    ($($args:tt)*) => {
        spenso::spenso_rep_atom!($crate::representations::ColorAdjoint {}; $($args)*)
    };
}

/// Defines operations related to manipulating abstract indices within symbolic expressions,
/// particularly relevant for physics calculations involving tensor structures and diagrams.
///
/// This trait provides methods for conjugating expressions, wrapping indices, and identifying
/// external ("dangling") indices.
pub trait IndexTooling {
    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        new_dummy: impl FnMut(usize) -> Aind,
    ) -> Atom;
    /// Wraps all abstract indices within the expression using a specified header symbol.
    ///
    /// This transforms indices like `mink(dim,idx)` into `mink(dim,header(idx))`. Useful for distinguishing
    /// between different copies of an expression, e.g., an amplitude and its complex conjugate.
    ///
    /// # Arguments
    /// * `header` - The [`Symbol`] to use as the wrapping function name.
    ///
    /// # Returns
    /// A new [`Atom`] with all indices wrapped.
    fn wrap_indices(&self, header: Symbol) -> Atom;

    fn spenso_conj(&self) -> Atom;

    /// Wraps only the dummy (contracted) abstract indices within the expression using a header symbol.
    ///
    /// Identifies indices that appear contracted (e.g., one covariant, one contravariant)
    /// and wraps only those, leaving external indices unchanged. Transforms `idx -> header(idx)`
    /// for dummy indices `idx`.
    ///
    /// # Arguments
    /// * `header` - The [`Symbol`] to use as the wrapping function name for dummy indices.
    ///
    /// # Returns
    /// A new [`Atom`] with only dummy indices wrapped.
    fn wrap_dummies<Aind: AbsInd + DummyAind + ParseableAind>(&self, header: Symbol) -> Atom;

    /// Computes the physics-aware conjugate of the expression.
    ///
    /// Applies conjugation rules specific to physics objects like spinors, gamma matrices,
    /// color representations, and the imaginary unit `i`. See implementation details for
    /// specific rules applied.
    ///
    /// # Returns
    /// A new [`Atom`] representing the conjugated expression.
    fn dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(&self) -> eyre::Result<Atom>;

    fn conjugate_transpose(&self, rep: impl RepName) -> Atom;

    /// Identifies and returns a list of dangling (external, uncontracted) indices.
    ///
    /// Analyzes the expression to find indices that are not summed over. Returns them
    /// as `Atom`s. Note that dual indices might be represented wrapped in a `dind` function.
    ///
    /// # Returns
    /// A `Vec<Atom>` where each `Atom` represents a dangling index.
    fn list_dangling<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Vec<Atom>;
}

impl IndexTooling for Atom {
    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        new_dummy: impl FnMut(usize) -> Aind,
    ) -> Atom {
        self.as_view().canonize(new_dummy)
    }

    fn spenso_conj(&self) -> Atom {
        self.as_view().spenso_conj()
    }
    fn wrap_indices(&self, header: Symbol) -> Atom {
        self.as_view().wrap_indices(header)
    }
    fn wrap_dummies<Aind: AbsInd + DummyAind + ParseableAind>(&self, header: Symbol) -> Atom {
        self.as_view().wrap_dummies::<Aind>(header)
    }
    fn dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(&self) -> eyre::Result<Atom> {
        self.as_view().dirac_adjoint::<Aind>()
    }

    fn conjugate_transpose(&self, rep: impl RepName) -> Atom {
        self.as_view().conjugate_transpose(rep)
    }
    fn list_dangling<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Vec<Atom> {
        self.as_view().list_dangling::<Aind>()
    }
}

#[cfg(test)]
mod syntax_macro_tests {
    #[allow(unused_imports)]
    use crate::{bis, coad, coaf, cof, color_d, color_d33, f, t};
    use spenso::shadowing::symbolica_utils::AtomCoreExt;

    #[test]
    fn representation_macros_build_surface_syntax() {
        crate::representations::initialize();

        insta::assert_snapshot!(bis!(4, i).to_bare_ordered_string(), @"bis(4,i)");
        insta::assert_snapshot!(cof!(Nc, i).to_bare_ordered_string(), @"cof(Nc,i)");
        insta::assert_snapshot!(coaf!(Nc, i).to_bare_ordered_string(), @"dind(cof(Nc,i))");
        insta::assert_snapshot!(coad!(Na, a).to_bare_ordered_string(), @"coad(Na,a)");
    }

    #[test]
    fn color_macros_build_surface_syntax() {
        crate::representations::initialize();

        let a = coad!(Na, a);
        let b = coad!(Na, b);
        let c = coad!(Na, c);

        insta::assert_snapshot!(t!(a.clone()).to_bare_ordered_string(), @"t(coad(Na,a),in,out)");
        insta::assert_snapshot!(f!(a.clone(), b.clone(), c.clone()).to_bare_ordered_string(), @"f(coad(Na,a),coad(Na,b),coad(Na,c))");
        insta::assert_snapshot!(color_d!(coad!(Na), a, b, c).to_bare_ordered_string(), @"d(coad(Na),coad(Na,a),coad(Na,b),coad(Na,c))");
        insta::assert_snapshot!(color_d33!(coad!(Na), coad!(Na)).to_bare_ordered_string(), @"d33(coad(Na),coad(Na))");
    }
}

#[derive(Error, Debug)]
pub enum AdjointError {
    #[error("Dummies already present:{0}")]
    DummiesAlready(Atom),
}

impl IndexTooling for AtomView<'_> {
    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        mut new_dummy: impl FnMut(usize) -> Aind,
    ) -> Atom {
        let mut net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
            .unwrap();

        // println!("{}", net.dot_pretty());

        let mut redual_reps = vec![];

        for t in net.store.tensors.iter_mut() {
            let mut reps = vec![];

            let mut pat = FunctionBuilder::new(t.name().unwrap());
            let mut rhs = pat.clone();
            for (i, s) in t.structure.external_structure_iter().enumerate() {
                if !s.rep_name().is_self_dual() && s.rep_name().is_dual() {
                    pat = pat.add_arg(
                        s.rep()
                            .dual()
                            .to_symbolic([Atom::var(symbol!(format!("a{i}_",)))]),
                    );

                    reps.push(Replacement::new(
                        s.rep().to_symbolic([Atom::var(RS.a_)]).to_pattern(),
                        s.rep().dual().to_symbolic([Atom::var(RS.a_)]),
                    ));
                } else {
                    pat = pat.add_arg(s.rep().to_symbolic([Atom::var(symbol!(format!("a{i}_",)))]));
                }
                rhs = rhs.add_arg(s.rep().to_symbolic([Atom::var(symbol!(format!("a{i}_",)))]));
            }
            if !reps.is_empty() {
                let rep = Replacement::new(pat.finish().to_pattern(), rhs.finish());
                // println!("{}", rep);
                redual_reps.push(rep);
                t.expression = t.expression.replace_multiple(&reps);
            }
        }

        let mut dummies = BTreeMap::new();

        for (p, _, d) in net.graph.graph.iter_edges() {
            if p.is_paired()
                && let NetworkEdge::Slot(s) = d.data
            {
                let (ind, group) = if s.rep_name().is_dual() {
                    (s.dual().to_atom(), s.rep().dual())
                } else {
                    (s.to_atom(), s.rep())
                };
                // println!("{}:{:?}", ind, group);
                dummies.insert(ind, group);
            }
        }

        let expr = net.simple_execute::<()>();

        let a = expr.canonize_tensors(dummies).unwrap();

        let mut reps = vec![];

        for (i, (d, r)) in a.dummy_indices.into_iter().enumerate() {
            reps.push(Replacement::new(
                d.to_pattern(),
                r.slot::<Aind, Aind>(new_dummy(i)).to_atom(),
            ));
        }

        a.canonical_form
            .replace_multiple(&reps)
            .replace_multiple(&redual_reps)
    }
    fn spenso_conj(&self) -> Atom {
        self.conj()
            .replace(Atom::var(RS.a__).conj())
            .with(INBUILTS.conj(RS.a__))
    }

    fn conjugate_transpose(&self, rep: impl RepName) -> Atom {
        let transpose_pat = function!(
            RS.a_,
            RS.a___,
            rep.to_symbolic([RS.d_, RS.i_]),
            rep.to_symbolic([RS.d_, RS.j_]),
            RS.b___
        )
        .to_pattern();

        let transpose_rhs = function!(
            RS.a_,
            RS.a___,
            rep.to_symbolic([RS.d_, RS.j_]),
            rep.to_symbolic([RS.d_, RS.i_]),
            RS.b___
        )
        .to_pattern();
        self.conj().replace(transpose_pat).with(transpose_rhs)
    }

    fn dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(&self) -> eyre::Result<Atom> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                take_first_term_from_sum: true,
                ..Default::default()
            })
            .unwrap();

        let bis_dangling: Vec<_> = net
            .graph
            .dangling_indices()
            .into_iter()
            .filter(|a| a.rep_name() == Bispinor {}.into())
            .collect();

        // println!("{}", net.dot_pretty());

        let mut a = self.spenso_conj();

        for i in bis_dangling {
            let mut dummy = i;
            dummy.aind = Aind::new_dummy();

            a = a.replace(i.to_atom()).with(dummy.to_atom())
                * function!(AGS.gamma0, i.to_atom(), dummy.to_atom());
        }
        let bis_graph = SuBitGraph::from_filter(&net.graph.graph, |e| {
            if let NetworkEdge::Slot(s) = e {
                s.rep_name() == Bispinor {}.into()
            } else {
                false
            }
        });

        let con = net.graph.graph.connected_components(&bis_graph);

        for c in con {
            let mut dangling: SuBitGraph = net.graph.graph.empty_subgraph();
            for i in c.included_iter() {
                if net.graph.graph.is_dangling(i) {
                    dangling.add(i);
                }
            }

            match dangling.n_included() {
                0 => {}
                1 => {}
                2 => {
                    let mut iter = dangling.included_iter();
                    let NetworkEdge::Slot(i) =
                        net.graph.graph[net.graph.graph[&iter.next().unwrap()]]
                    else {
                        break;
                    };

                    let NetworkEdge::Slot(j) =
                        net.graph.graph[net.graph.graph[&iter.next().unwrap()]]
                    else {
                        break;
                    };

                    a = a.replace_multiple(&[
                        Replacement::new(i.to_atom().to_pattern(), j.to_atom()),
                        Replacement::new(j.to_atom().to_pattern(), i.to_atom()),
                    ]);
                }
                _ => {
                    return Err(eyre!(
                        "Too many dangling bispinors for complex conjugation {}",
                        net.dot_pretty()
                    ));
                }
            }
        }
        Ok(a.simplify_gamma_conj::<Aind>()?
            .simplify_gamma0()
            .simplify_metrics())
    }

    fn wrap_indices(&self, header: Symbol) -> Atom {
        wrap_indices_impl(*self, header)
    }
    fn wrap_dummies<Aind: AbsInd + DummyAind + ParseableAind>(&self, header: Symbol) -> Atom {
        wrap_dummies_impl::<Aind>(*self, header)
    }
    fn list_dangling<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Vec<Atom> {
        list_dangling_impl::<Aind>(*self)
    }
}

#[cfg(test)]
pub mod test {
    use insta::assert_snapshot;
    use spenso::{
        p, shadowing::symbolica_utils::AtomCoreExt, slot, structure::abstract_index::AbstractIndex,
    };
    use symbolica::{atom::AtomCore, parse_lit, printer::CanonicalOrderingSettings};

    use crate::{Cookable, IndexTooling, gamma, test_support::test_initialize, u};

    #[test]
    fn gamma_conj() {
        let a = test_initialize();
        let mink4 = a.mink4;

        let bis4 = a.bis4;

        let p1 = p!(1, slot!(mink4, 3));

        let ubgu = u!(2, slot!(bis4, 2))
            * (gamma!(3, 1, 2) * p1 + bis4.g(1, 2))
            * (u!(1, slot!(bis4, 1))
                .dirac_adjoint::<AbstractIndex>()
                .unwrap());

        assert_snapshot!(
            ubgu.dirac_adjoint::<AbstractIndex>()
                .unwrap()
                .canonize(AbstractIndex::Dummy)
                .to_canonically_ordered_string(CanonicalOrderingSettings {
                    include_attributes: false,
                    include_namespace: false,hide_namespace:None,
                }),@"(conj(p(1,mink(4,d_0)))*gamma(bis(4,d_1),bis(4,d_2),mink(4,d_0))+g(bis(4,d_1),bis(4,d_2)))*conj(u(2,bis(4,d_3)))*gamma0(bis(4,d_1),bis(4,d_3))*u(1,bis(4,d_2))"
        );
    }

    #[test]
    fn canonize_color() {
        test_initialize();
        let expr = parse_lit!(
            f(coad(8, hedge(4)), coad(8, hedge(8)), coad(8, hedge(12)))
                * t(coad(8, hedge(4)), cof(3, hedge(16)), dind(cof(3, hedge(2))))
                * t(coad(8, hedge(8)), cof(3, hedge(2)), dind(cof(3, hedge(10))))
                * t(
                    coad(8, hedge(12)),
                    cof(3, hedge(10)),
                    dind(cof(3, hedge(16)))
                ),
            default_namespace = "spenso"
        );

        let can = expr.cook_indices();

        println!("{}", can);
        let can = can.canonize::<AbstractIndex>(AbstractIndex::Dummy);

        assert_snapshot!(
            can.to_bare_ordered_string(),@"f(coad(8,d_0),coad(8,d_1),coad(8,d_2))*t(coad(8,d_0),cof(3,d_3),dind(cof(3,d_4)))*t(coad(8,d_1),cof(3,d_4),dind(cof(3,d_5)))*t(coad(8,d_2),cof(3,d_5),dind(cof(3,d_3)))"
        );

        //  gives
        //
        // f(coad(8,dummy(0)),coad(8,dummy(1)),coad(8,dummy(2)))*t(coad(8,dummy(0)),cof(3,dummy(3)),cof(3,hedge(2)))*t(coad(8,dummy(1)),cof(3,hedge(2)),cof(3,hedge(10)))*t(coad(8,dummy(2)),cof(3,hedge(10)),cof(3,dummy(3)))
    }
}
