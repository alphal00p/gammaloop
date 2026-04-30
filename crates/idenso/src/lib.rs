#![allow(uncommon_codepoints)]

use std::collections::BTreeMap;

use base64::{Engine, engine::general_purpose::URL_SAFE};
use eyre::eyre;
use linnet::half_edge::subgraph::{BaseSubgraph, ModifySubSet, SuBitGraph, SubSetLike};
use metric::{
    CookingError, cook_function_view, cook_indices_impl, list_dangling_impl, wrap_dummies_impl,
    wrap_indices_impl,
};
use spenso::{
    network::{
        graph::NetworkEdge,
        library::function_lib::INBUILTS,
        parsing::{NetworkParse, ParseSettings},
    },
    structure::{
        HasName, TensorStructure,
        representation::RepName,
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, UserData},
    function,
    id::Replacement,
    symbol, tag,
};
use thiserror::Error;

use crate::{
    gamma::{AGS, GammaSimplifier},
    metric::MetricSimplifier,
    rep_symbols::RS,
    representations::Bispinor,
    tensor::{SymbolicNetExt, SymbolicNetParse},
};

pub mod tensor;

pub mod color;
pub mod gamma;
pub mod metric;
pub mod parsing_ind;
#[cfg(feature = "python")]
pub mod python;
pub mod rep_symbols;
pub mod representations;
pub mod schoonschip;

#[macro_export]
macro_rules!  symbol_set {
        // Identifier symbols with explicit struct and static names
        ($struct_name:ident, $static_name:ident; $($char:ident)*) => {
            #[allow(non_snake_case)]
            pub struct $struct_name {
                $(pub $char: Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($char: symbol!(stringify!($char)),)*
            });
        };

        // String literals as individual statics
        (statics; $($field:ident : $string:literal),* $(,)?) => {
            $(
                #[allow(non_upper_case_globals)]
                pub static $field: ::std::sync::LazyLock<Symbol> = ::std::sync::LazyLock::new(|| symbol!($string));
            )*
        };

        // String literals grouped in a struct
        ($struct_name:ident, $static_name:ident; $($field:ident : $string:literal),* $(,)?) => {
            #[allow(non_snake_case)]
            pub struct $struct_name {
                $(pub $field: Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($field: symbol!($string),)*
            });
        };
    }

symbol_set!(Wildcards, W_;
    a_ b_ c_ d_ e_ f_ g_ h_ i_ j_ k_ l_ m_ n_ o_ p_ q_ r_ s_ t_ u_ v_ w_ x_ y_ z_
    a__ b__ c__ d__ e__ f__ g__ h__ i__ j__ k__ l__ m__ n__ o__ p__ q__ r__ s__ t__ u__ v__ w__ x__ y__ z__
    a___ b___ c___ d___ e___ f___ g___ h___ i___ j___ k___ l___ m___ n___ o___ p___ q___ r___ s___ t___ u___ v___ w___ x___ y___ z___
);

/// Defines operations related to manipulating abstract indices within symbolic expressions,
/// particularly relevant for physics calculations involving tensor structures and diagrams.
///
/// This trait provides methods for conjugating expressions, wrapping indices (both all and
/// only dummy/contracted ones), simplifying indices ("cooking"), and identifying external
/// ("dangling") indices.
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
    fn wrap_dummies<Aind: AbsInd + ParseableAind>(&self, header: Symbol) -> Atom;

    /// Simplifies structured indices within function arguments into flattened variable symbols.
    ///
    /// Replaces indices like `mink(4, mu)` inside function arguments (not top-level) with
    /// unique variable symbols like `var_mink_4_mu`. This can aid pattern matching.
    ///
    /// # Returns
    /// A new [`Atom`] with "cooked" indices.
    fn cook_indices(&self) -> Atom;

    /// Converts a single function [`Atom`] into a flattened variable symbol based on its name and arguments.
    ///
    /// Expects the input `Atom` to be a function. Returns a variable `Atom` whose symbol name
    /// encodes the original function and its arguments (e.g., `f(a, b)` might become `var_f_a_b`).
    /// Fails if the input is not a function or if arguments are not convertible (e.g., polynomials).
    ///
    /// # Returns
    /// `Ok(Atom)` containing the new variable symbol on success.
    /// `Err(CookingError)` if the input cannot be cooked.
    fn cook_function(&self) -> Result<Atom, CookingError>;

    fn cook(&self) -> Atom;

    fn uncook(&self) -> Atom;

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
    fn list_dangling<Aind: AbsInd + ParseableAind>(&self) -> Vec<Atom>;
}

impl IndexTooling for Atom {
    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        new_dummy: impl FnMut(usize) -> Aind,
    ) -> Atom {
        self.as_view().canonize(new_dummy)
    }

    fn cook(&self) -> Atom {
        self.as_view().cook()
    }

    fn uncook(&self) -> Atom {
        self.as_view().uncook()
    }
    fn spenso_conj(&self) -> Atom {
        self.as_view().spenso_conj()
    }
    fn wrap_indices(&self, header: Symbol) -> Atom {
        self.as_view().wrap_indices(header)
    }
    fn wrap_dummies<Aind: AbsInd + ParseableAind>(&self, header: Symbol) -> Atom {
        self.as_view().wrap_dummies::<Aind>(header)
    }
    fn cook_indices(&self) -> Atom {
        self.as_view().cook_indices()
    }

    fn cook_function(&self) -> Result<Atom, CookingError> {
        self.as_view().cook_function()
    }

    fn dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(&self) -> eyre::Result<Atom> {
        self.as_view().dirac_adjoint::<Aind>()
    }

    fn conjugate_transpose(&self, rep: impl RepName) -> Atom {
        self.as_view().conjugate_transpose(rep)
    }
    fn list_dangling<Aind: AbsInd + ParseableAind>(&self) -> Vec<Atom> {
        self.as_view().list_dangling::<Aind>()
    }
}

#[derive(Error, Debug)]
pub enum AdjointError {
    #[error("Dummies already present:{0}")]
    DummiesAlready(Atom),
}

impl IndexTooling for AtomView<'_> {
    fn cook(&self) -> Atom {
        self.replace_map(|a, _, out| {
            if let AtomView::Fun(f) = a {
                let hash = blake3::hash(a.get_data());
                let bytes12 = &hash.as_bytes()[..12];
                let a = symbol!(
                    URL_SAFE.encode(bytes12) + f.get_symbol().get_name(),
                    data = UserData::Atom(a.to_owned()),
                    tag = tag!("cooked")
                );
                **out = Atom::var(a);
            }
        })
    }

    fn uncook(&self) -> Atom {
        self.replace_map(|a, _, out| {
            if let AtomView::Var(s) = a {
                let s = s.get_symbol();
                if s.has_tag("idenso::cooked")
                    && let UserData::Atom(a) = s.get_data()
                {
                    **out = a.clone();
                }
            }
        })
    }
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

    fn cook_function(&self) -> Result<Atom, CookingError> {
        cook_function_view(*self)
    }
    fn wrap_indices(&self, header: Symbol) -> Atom {
        wrap_indices_impl(*self, header)
    }
    fn cook_indices(&self) -> Atom {
        cook_indices_impl(*self)
    }
    fn wrap_dummies<Aind: AbsInd + ParseableAind>(&self, header: Symbol) -> Atom {
        wrap_dummies_impl::<Aind>(*self, header)
    }
    fn list_dangling<Aind: AbsInd + ParseableAind>(&self) -> Vec<Atom> {
        list_dangling_impl::<Aind>(*self)
    }
}

#[cfg(test)]
pub mod test {
    use std::sync::LazyLock;

    pub struct TestSymbols {
        pub p: Symbol,
        pub a: Symbol,
        pub u: Symbol,
    }

    pub fn test_initialize() {
        let _ = TS.p;
        initialize();
    }

    pub static TS: LazyLock<TestSymbols> = LazyLock::new(|| TestSymbols {
        p: symbol!("p";Real),
        a: symbol!("a"),
        u: symbol!("u"),
    });

    use insta::assert_snapshot;
    use spenso::structure::{
        IndexlessNamedStructure,
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    };
    use symbolica::{
        atom::{Atom, AtomCore, Symbol},
        function, parse_lit,
        printer::CanonicalOrderingSettings,
        symbol,
    };

    use crate::{
        IndexTooling,
        gamma::AGS,
        metric::PermuteWithMetric,
        representations::{Bispinor, initialize},
    };
    pub fn gamma(
        i: impl Into<AbstractIndex>,
        j: impl Into<AbstractIndex>,
        mu: impl Into<AbstractIndex>,
    ) -> Atom {
        let gamma_strct = IndexlessNamedStructure::<Symbol, ()>::from_iter(
            [
                Bispinor {}.new_rep(4).to_lib(),
                Bispinor {}.new_rep(4).cast(),
                Minkowski {}.new_rep(4).cast(),
            ],
            AGS.gamma,
            None,
        );
        gamma_strct
            .reindex([i.into(), j.into(), mu.into()])
            .unwrap()
            .permute_with_metric()
    }

    pub fn u(i: usize, m: impl Into<AbstractIndex>) -> Atom {
        let m_atom: AbstractIndex = m.into();
        let m_atom: Atom = m_atom.into();
        let mink = Bispinor {}.new_rep(4);
        function!(TS.u, i, mink.to_symbolic([m_atom]))
    }

    #[allow(dead_code)]
    pub fn p(i: usize, m: impl Into<AbstractIndex>) -> Atom {
        let m_atom: AbstractIndex = m.into();
        let m_atom: Atom = m_atom.into();
        let mink = Minkowski {}.new_rep(4);
        function!(TS.p, i, mink.to_symbolic([m_atom]))
    }

    #[allow(dead_code)]
    pub fn a(i: usize, m: impl Into<AbstractIndex>, n: impl Into<AbstractIndex>) -> Atom {
        let m_atom: AbstractIndex = m.into();
        let m_atom: Atom = m_atom.into();

        let n_atom: AbstractIndex = n.into();
        let n_atom: Atom = n_atom.into();
        let mink = Bispinor {}.new_rep(4);
        function!(
            TS.a,
            i,
            mink.to_symbolic([m_atom]),
            mink.to_symbolic([n_atom])
        )
    }

    #[test]
    fn gamma_conj() {
        test_initialize();
        // let expr = gamma(1, 2, 3).dirac_adjoint::<AbstractIndex>().unwrap();
        // println!("{expr}");
        //

        let ubgu = u(2, 2)
            * (gamma(1, 2, 3) * p(1, 3) + Bispinor {}.metric_atom([4, 1], [4, 2]))
            * (u(1, 1).dirac_adjoint::<AbstractIndex>().unwrap());
        // println!("Start:\n{}", ubgu);
        // let expr = ubgu.dirac_adjoint::<AbstractIndex>().unwrap();
        // println!("{expr}");

        // println!("conj trans\n {exp}");
        // exp = exp.simplify_gamma_conj::<AbstractIndex>().unwrap();
        // println!("conj gamma simplify \n{exp}");
        // exp = exp.simplify_gamma0::<AbstractIndex>();
        // println!("gamma0 simplify \n{exp}");
        // exp = exp.simplify_metrics();
        // println!("simplify metrics \n{exp}");

        assert_snapshot!(
            ubgu.dirac_adjoint::<AbstractIndex>()
                .unwrap()
                .canonize(AbstractIndex::Dummy)
                .to_canonically_ordered_string(CanonicalOrderingSettings {
                    include_attributes: false,
                    include_namespace: false,hide_namespace:None,
                }),@"(g(bis(4,d_1),bis(4,d_2))+gamma(bis(4,d_1),bis(4,d_2),mink(4,d_0))*p(1,mink(4,d_0)))*conj(u(2,bis(4,d_3)))*gamma0(bis(4,d_1),bis(4,d_3))*u(1,bis(4,d_2))"
        );

        // let ubggu = u(2, 3)
        //     * gamma(1, 2, 3)
        //     * gamma(2, 3, 1)
        //     * (u(1, 1).dirac_adjoint::<AbstractIndex>().unwrap());

        // println!("Ubggu:\n{}", ubggu);
        // println!(
        //     "Conjugate:\n{}",
        //     ubggu.dirac_adjoint::<AbstractIndex>().unwrap()
        // );

        // let ub = u(1, 1).dirac_adjoint::<AbstractIndex>().unwrap();
        // let u = u(2, 2);

        // let ubgu = &ub * gamma(1, 2, 3) * &u;
        // let cubgu = ubgu.dirac_adjoint::<AbstractIndex>().unwrap();

        // let ccubgu = ub.dirac_adjoint::<AbstractIndex>().unwrap()
        //     * gamma(1, 2, 3).dirac_adjoint::<AbstractIndex>().unwrap()
        //     * u.dirac_adjoint::<AbstractIndex>().unwrap();
        // println!("Start:\n{}", ubgu);
        // println!("Oneshot:\n{}", cubgu);
        // println!(
        //     "Separate:\n{}",
        //     ccubgu
        //         .simplify_gamma0::<AbstractIndex>()
        //         .replace(function!(symbol!("spenso::u"), RS.a__).spenso_conj())
        //         .with(function!(symbol!("spenso::uconj"), RS.a__))
        //         .canonize(AbstractIndex::Dummy)
        // );

        // let a = a(1, 1, 2) * a(2, 2, 4) * a(1, 4, 3) + a(1, 1, 4) * a(3, 4, 5) * a(1, 5, 3);

        // println!("{a}");
        // println!("{}", a.canonize(AbstractIndex::Dummy));
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

        let can = expr
            .cook_indices()
            .canonize::<AbstractIndex>(AbstractIndex::Dummy);

        assert_snapshot!(
            can.to_canonically_ordered_string(CanonicalOrderingSettings {
                include_namespace: false,
                include_attributes: false,
                hide_namespace: None
            }),@"f(coad(8,d_0),coad(8,d_1),coad(8,d_2))*t(coad(8,d_0),cof(3,d_3),dind(cof(3,d_4)))*t(coad(8,d_1),cof(3,d_4),dind(cof(3,d_5)))*t(coad(8,d_2),cof(3,d_5),dind(cof(3,d_3)))"
        );

        //  gives
        //
        // f(coad(8,dummy(0)),coad(8,dummy(1)),coad(8,dummy(2)))*t(coad(8,dummy(0)),cof(3,dummy(3)),cof(3,hedge(2)))*t(coad(8,dummy(1)),cof(3,hedge(2)),cof(3,hedge(10)))*t(coad(8,dummy(2)),cof(3,hedge(10)),cof(3,dummy(3)))
    }
}
