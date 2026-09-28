#![allow(uncommon_codepoints)]

use std::collections::{BTreeMap, BTreeSet};

use linnet::{
    half_edge::subgraph::{BaseSubgraph, ModifySubSet, SuBitGraph, SubSetLike},
    tree::child_vec::ChildVecStore,
};
use shorthands::metric::{list_dangling_impl, wrap_dummies_impl, wrap_indices_impl};
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode},
        library::function_lib::INBUILTS,
        library::symbolic::ETS,
        parsing::{AtomStructureExt, ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
        store::TensorScalarStore,
        tags::SPENSO_TAG,
    },
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        HasName, OrderedStructure, TensorStructure, ToSymbolic,
        representation::{Euclidean, Lorentz, Minkowski, RepName},
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
    symbol_set,
    symbolica_init::in_symbolica_initializer,
    tensor_symbol,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    id::{AliasedAtom, Replacement},
    initialize,
    printer::AtomPrinter,
    symbol,
};
use thiserror::Error;

use crate::{
    color::{CS, color_conj_impl},
    dirac::AGS,
    epsilon::EPSILON_SYMBOL,
    rep_symbols::RS,
    representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet, SpinFundamental},
    shorthands::{UndoShorthands, chain::Chain},
    tensor::{SymbolicNetExt, SymbolicNetParse},
};

initialize!(|| {
    in_symbolica_initializer(|| {
        let _ = Minkowski {}.to_symbolic([Atom::Zero]);
        let _ = Euclidean {}.to_symbolic([Atom::Zero]);
        let _ = Lorentz {}.to_symbolic([Atom::Zero]);
        let _ = SpinFundamental {}.to_symbolic([Atom::Zero]);
        let _ = Bispinor {}.to_symbolic([Atom::Zero]);
        let _ = ColorAdjoint {}.to_symbolic([Atom::Zero]);
        let _ = ColorFundamental {}.to_symbolic([Atom::Zero]);
        let _ = ColorSextet {}.to_symbolic([Atom::Zero]);
        let _ = RS.force_in_initializer().a_;
        let _ = AGS.force_in_initializer().gamma;
        let _ = *EPSILON_SYMBOL.force_in_initializer();
        let _ = CS.force_in_initializer().cf;
        let _ = ETS.force_in_initializer().delta;
        let _ = ETS.force_in_initializer().metric;
    });
});

pub mod tensor;

pub mod color;
pub mod cook;
pub mod dirac;
pub mod epsilon;
// pub mod parsing_ind;
#[cfg(any(test, feature = "reference-cases"))]
pub mod reference_cases;
pub mod rep_symbols;
pub mod representations;
pub mod shorthands;
#[cfg(test)]
pub(crate) mod test_support;

pub use cook::{
    CookMode, CookOutputTags, CookSettings, CookSourceFilter, CookTagFilter, Cookable, CookingError,
};

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

/// Errors produced while inspecting or rewriting tensor indices.
#[derive(Debug, Error)]
pub enum IndexToolingError {
    #[error("cannot list dangling indices: {reason}")]
    ListDangling { reason: String },
    #[error("cannot wrap dummy indices: {reason}")]
    WrapDummies { reason: String },
}

/// Errors produced while parsing or evaluating symbolic tensor networks.
#[derive(Debug, Error)]
pub enum NetworkToolingError {
    #[error("cannot parse tensor network: {reason}")]
    Parse { reason: String },
    #[error("cannot execute tensor network: {reason}")]
    Execute { reason: String },
    #[error("cannot extract tensor-network result: {reason}")]
    Result { reason: String },
}

/// Errors produced while canonicalizing a symbolic tensor expression.
#[derive(Debug, Error)]
pub enum CanonicalizationError {
    #[error(transparent)]
    Tensor(#[from] tensor::inference::TensorInferenceError),
    #[error(transparent)]
    Network(#[from] NetworkToolingError),
    #[error("cannot prepare tensor for canonicalization: {reason}")]
    Prepare { reason: String },
    #[error("cannot canonicalize tensor indices: {reason}")]
    Indices { reason: String },
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
    ) -> Result<Atom, CanonicalizationError>;
    /// Wraps all abstract indices within the expression using a specified header symbol.
    ///
    /// This transforms `mink(dim,idx)` into `mink(dim,spenso::index_scope(header,idx))`, retaining a valid abstract index. Useful for distinguishing
    /// between different copies of an expression, e.g., an amplitude and its complex conjugate.
    ///
    /// # Arguments
    /// * `header` - The [`Symbol`] identifying the independent index scope.
    ///
    /// # Returns
    /// A new [`Atom`] with all indices wrapped.
    fn wrap_indices(&self, header: Symbol) -> Atom;

    fn spenso_conj(&self) -> Atom;
    fn spenso_print<'a>(&'a self, settings: &SpensoPrintSettings) -> AtomPrinter<'a>;
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
    /// A new [`Atom`] with only dummy indices wrapped, or an error when the expression cannot be
    /// parsed as a tensor network.
    fn wrap_dummies<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        header: Symbol,
    ) -> Result<Atom, IndexToolingError>;

    /// Computes the physics-aware conjugate of the expression.
    ///
    /// Applies conjugation rules specific to physics objects like spinors, gamma matrices,
    /// color representations, and the imaginary unit `i`. See implementation details for
    /// specific rules applied.
    ///
    /// # Returns
    /// A new [`Atom`] representing the conjugated expression.
    /// When `preserve_indices` is true, external labels remain attached to the
    /// same amplitude legs. Gamma-zero boundary factors and conjugation still
    /// apply, but open-chain endpoints are not exchanged. This convention is
    /// linear across sums with different external fermion pairings.
    fn dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(
        &self,
        preserve_indices: bool,
    ) -> Result<Atom, AdjointError>;

    fn conjugate_transpose(&self, rep: impl RepName) -> Atom;

    /// Identifies and returns a list of dangling (external, uncontracted) indices.
    ///
    /// Analyzes the expression to find indices that are not summed over. Returns them
    /// as `Atom`s. Note that dual indices might be represented wrapped in a `dind` function.
    ///
    /// # Returns
    /// A `Vec<Atom>` where each `Atom` represents a dangling index, or an error when the
    /// expression cannot be parsed as a tensor network.
    fn list_dangling<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Vec<Atom>, IndexToolingError>;

    fn alias_subtensors(&self, tensor_name: &str) -> AliasedAtom;
}

impl IndexTooling for Atom {
    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        new_dummy: impl FnMut(usize) -> Aind,
    ) -> Result<Atom, CanonicalizationError> {
        self.as_view().canonize(new_dummy)
    }

    fn alias_subtensors(&self, tensor_name: &str) -> AliasedAtom {
        self.as_view().alias_subtensors(tensor_name)
    }

    fn spenso_print<'a>(&'a self, settings: &SpensoPrintSettings) -> AtomPrinter<'a> {
        self.printer(settings.nice_symbolica())
    }

    fn spenso_conj(&self) -> Atom {
        self.as_view().spenso_conj()
    }
    fn wrap_indices(&self, header: Symbol) -> Atom {
        self.as_view().wrap_indices(header)
    }
    fn wrap_dummies<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        header: Symbol,
    ) -> Result<Atom, IndexToolingError> {
        self.as_view().wrap_dummies::<Aind>(header)
    }
    fn dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(
        &self,
        preserve_indices: bool,
    ) -> Result<Atom, AdjointError> {
        self.as_view().dirac_adjoint::<Aind>(preserve_indices)
    }

    fn conjugate_transpose(&self, rep: impl RepName) -> Atom {
        self.as_view().conjugate_transpose(rep)
    }
    fn list_dangling<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Vec<Atom>, IndexToolingError> {
        self.as_view().list_dangling::<Aind>()
    }
}

#[cfg(test)]
mod syntax_macro_tests {
    #[allow(unused_imports)]
    use crate::{bis, coad, coaf, cof, color_d, f, t};
    use symbolica_utils::AtomPrintExt;

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
    }
}

#[derive(Error, Debug)]
pub enum AdjointError {
    #[error("Dummies already present:{0}")]
    DummiesAlready(Atom),
    #[error(transparent)]
    Network(#[from] NetworkToolingError),
    #[error(transparent)]
    Tensor(#[from] tensor::inference::TensorInferenceError),
    #[error("cannot encode Dirac-adjoint indices: {reason}")]
    IndexEncoding { reason: String },
    #[error(
        "cannot construct Dirac adjoint: bispinor component has {count} dangling indices (expected at most 2)"
    )]
    TooManyDanglingBispinors { count: usize },
}

impl IndexTooling for AtomView<'_> {
    fn spenso_print(&self, settings: &SpensoPrintSettings) -> AtomPrinter<'_> {
        self.printer(settings.nice_symbolica())
    }
    fn alias_subtensors(&self, tensor_name: &str) -> AliasedAtom {
        let tensor_symbol = tensor_symbol!(&tensor_name);
        self.alias_subexpressions(|a, _count, i| {
            if a.has_attributes_of(SPENSO_TAG.rep_) || a.has_attributes_of(SPENSO_TAG.tensor_) {
                None
            } else {
                match a.infer_structure::<OrderedStructure>() {
                    Ok(a) => Some(a.to_symbolic_with(
                        tensor_symbol,
                        &[Atom::num(i), Atom::num(_count)],
                        None,
                    )),
                    Err(_) => Some(function!(symbol!("se"), i, _count)),
                }
            }
        })
    }

    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        mut new_dummy: impl FnMut(usize) -> Aind,
    ) -> Result<Atom, CanonicalizationError> {
        // Admit the complete source before algebraic zero pruning. Callers
        // using structured index payloads must cook_indices at this boundary,
        // consistently for both vanishing and surviving expressions.
        let mut net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.remove_antisymmetric_zero_terms();

        let mut redual_reps = vec![];

        for t in net.store.tensors.iter_mut() {
            if t.expression.is_zero() {
                continue;
            }
            let mut reps = vec![];

            let name = t.name().ok_or_else(|| CanonicalizationError::Prepare {
                reason: format!("tensor `{}` has no function name", t.expression),
            })?;
            let mut pat = FunctionBuilder::new(name);
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
                t.invalidate_proofs();
                t.expression = t.expression.replace_multiple(&reps);
            }
        }

        let mut dummies = BTreeMap::new();
        for (p, _, d) in net.graph.graph.iter_edges() {
            if p.is_paired()
                && let NetworkEdge::Slot(s) = d.data
            {
                let slot = if s.rep_name().is_dual() { s.dual() } else { *s };
                // println!("{}:{:?}", ind, group);
                dummies.insert(slot.to_atom(), slot.rep());
            }
        }

        let external = net.graph.dangling_indices();
        // Reconstruct the admitted operation graph, without running tensor algebra
        // merely to hand its syntax to Symbolica's canonicalizer.
        let tree = net.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = net.graph.graph.node_id(net.graph.head());
        let expr = net
            .graph
            .to_expression_at(&tree, root, &mut |_, node| match node {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    Ok(Some(net.store.tensors[*index].expression.clone()))
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                    Ok(Some(net.store.get_scalar_ref(*index).clone()))
                }
                NetworkNode::Op(_) => Ok(None),
                _ => Err(spenso::network::TensorNetworkError::Other(eyre::eyre!(
                    "canonicalization requires stored symbolic leaves"
                ))),
            })
            .map_err(|error| CanonicalizationError::Prepare {
                reason: error.to_string(),
            })?;
        // Symbolica rejects tensor-containing negative powers even when their
        // base is a closed scalar network. Canonicalize those independent
        // contractions recursively, then keep them opaque during outer labeling.
        let mut aliases = BTreeMap::new();
        let mut used = expr.get_all_symbols(true);
        let mut serial = 0;
        let mut failure = None;
        let expr = expr.replace_map(|atom, _, out| {
            let AtomView::Pow(power) = atom else {
                return;
            };
            let (base, exponent) = power.get_base_exp();
            if i64::try_from(exponent).is_ok_and(|n| n > 0) {
                return;
            }
            let canonical =
                match base.canonize::<Aind>(&mut new_dummy as &mut dyn FnMut(usize) -> Aind) {
                    Ok(canonical) => canonical.pow(exponent),
                    Err(error) => {
                        failure = Some(error);
                        return;
                    }
                };
            let alias = aliases.entry(canonical).or_insert_with(|| {
                loop {
                    let name = symbol!(&format!("idenso::canonical_scalar_{serial}"));
                    serial += 1;
                    if used.insert(name) {
                        break Atom::var(name);
                    }
                }
            });
            **out = alias.clone();
        });
        if let Some(error) = failure {
            return Err(error);
        }
        let mut scalar_definitions = AliasedAtom::from(Atom::Zero);
        for (original, alias) in aliases {
            scalar_definitions.register_alias(alias, original);
        }
        // Spenso supplies validated incidence and external slots. Symbolica
        // canonizes the reconstructed expression, including dummy pairs inside
        // nested sums, without distributing any products of sums.
        let canonical =
            expr.canonize_tensors(dummies)
                .map_err(|error| CanonicalizationError::Indices {
                    reason: error.to_string(),
                })?;

        // Sum branches also have paired copies of genuinely external slots.
        // Preserve those labels after Symbolica has validated their incidence.
        let mut reserved = BTreeSet::new();
        for slot in external {
            let slot = if slot.rep_name().is_dual() {
                slot.dual()
            } else {
                slot
            };
            let slot = slot.to_atom();
            reserved.insert(slot);
        }

        let mut reps = vec![];
        let mut next_dummy = 0;

        // Only relabel names that Symbolica validated as contracted indices;
        // keep their replacements distinct from the genuine external slots.
        // Canceled terms can leave unused names in the validated index pool.
        for (d, r) in canonical
            .dummy_indices
            .into_iter()
            .filter(|(index, _)| canonical.canonical_form.contains(index))
        {
            let target = loop {
                let candidate = r.slot::<Aind, Aind>(new_dummy(next_dummy)).to_atom();
                next_dummy += 1;
                if reserved.insert(candidate.clone()) {
                    break candidate;
                }
            };
            reps.push(Replacement::new(d.to_pattern(), target));
        }

        Ok(scalar_definitions
            .map_root(|_| {
                canonical
                    .canonical_form
                    .replace_multiple(&reps)
                    .replace_multiple(&redual_reps)
            })
            .into_inner())
    }
    fn spenso_conj(&self) -> Atom {
        let conjugate = self
            .conj()
            .replace(Atom::var(RS.a__).conj())
            .with(INBUILTS.conj(RS.a__));
        color_conj_impl(conjugate.as_view())
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

    fn dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(
        &self,
        preserve_indices: bool,
    ) -> Result<Atom, AdjointError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: false,
                        expand_inside_chains: true,
                        expand_schoonship: true,
                    },
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        // Non-spinor words need only their existing conjugation. Retain compact
        // color chains and scalar spectators after the same checked parse.
        // Closed Dirac traces still carry internal bispinor edges.
        if !net.graph.graph.iter_edges().any(|(_, _, edge)| {
            matches!(edge.data, NetworkEdge::Slot(slot)
                if slot.rep_name() == Bispinor {}.into())
        }) {
            return Ok(self.spenso_conj());
        }

        let bis_dangling: Vec<_> = net
            .graph
            .dangling_indices()
            .into_iter()
            .filter(|a| a.rep_name() == Bispinor {}.into())
            .collect();

        // println!("{}", net.dot_pretty());

        let mut transpositions = Vec::new();
        // Only transposed output indices require a choice of matrix channels.
        // A physical adjoint keeps every external leg, including sums whose
        // alternatives pair those legs differently.
        if !preserve_indices {
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

                        transpositions.push((i, j));
                    }
                    _ => {
                        return Err(AdjointError::TooManyDanglingBispinors {
                            count: dangling.n_included(),
                        });
                    }
                }
            }
        }
        // Materialize matrix channels in this same parse before conjugation.
        // Scalar inner products stay compact, and every sum branch is retained.
        let explicit = net.simple_execute::<()>()?;
        let mut a = explicit.spenso_conj();
        for i in bis_dangling {
            let mut dummy = i;
            dummy.aind = Aind::new_dummy();
            a = a.replace(i.to_atom()).with(dummy.to_atom())
                * function!(AGS.gamma0, i.to_atom(), dummy.to_atom());
        }
        for (i, j) in transpositions {
            a = a.replace_multiple(&[
                Replacement::new(i.to_atom().to_pattern(), j.to_atom()),
                Replacement::new(j.to_atom().to_pattern(), i.to_atom()),
            ]);
        }
        a = dirac::DiracSimplifier::conjugate_matrices::<Aind>(a.as_view());
        // The generic index encoding is reversible; shared algebra sees the
        // same strict slot grammar as every other SymbolicTensor operation.
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let factored = dirac::DiracSimplifier::factor_gamma_zero(a.as_view());
        let cooked =
            cooking
                .try_cook(factored.as_view())
                .map_err(|error| AdjointError::IndexEncoding {
                    reason: format!("{error:?}"),
                })?;
        // Formal spin dimensions belong to chain channels. Conjugation and
        // gamma0 factoring have already acted on explicit matrices; normalize
        // their channels before the same strict typed admission used elsewhere.
        let tensor = tensor::SymbolicTensor::infer(cooked.chainify(Bispinor {}.into()))?;
        let resolved = if preserve_indices {
            tensor
                .simplify_gamma(dirac::GammaSimplifySettings {
                    gamma0: true,
                    evaluate_traces: false,
                    ..Default::default()
                })?
                .resolved()?
        } else {
            tensor.contract(Default::default())?.resolved()?
        };
        let expression = cooking.uncook(resolved.expression.as_view());
        if preserve_indices {
            Ok(expression.undo_chain::<Aind>()?)
        } else {
            // Restore one-word matrix heads at this raw expression boundary.
            // Their existing symmetry normalizers (notably gamma0) own physical
            // endpoint ordering; longer formal chains remain factored.
            Ok(expression.undo_single_length())
        }
    }

    fn wrap_indices(&self, header: Symbol) -> Atom {
        wrap_indices_impl(*self, header)
    }
    fn wrap_dummies<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        header: Symbol,
    ) -> Result<Atom, IndexToolingError> {
        wrap_dummies_impl::<Aind>(*self, header)
    }
    fn list_dangling<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Vec<Atom>, IndexToolingError> {
        list_dangling_impl::<Aind>(*self)
    }
}

#[cfg(test)]
pub mod test {
    use insta::assert_snapshot;
    use spenso::{network::tags::SPENSO_TAG, p, slot, structure::abstract_index::AbstractIndex};
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
        parse_lit,
        printer::CanonicalOrderingSettings,
        symbol,
    };
    use symbolica_utils::AtomPrintExt;

    use crate::{
        AdjointError, CanonicalizationError, Cookable, IndexTooling, NetworkToolingError, gamma,
        test_support::test_initialize, u,
    };

    #[test]
    fn dirac_adjoint_retains_nonspinor_conjugation_and_reduces_closed_spinor_words() {
        test_initialize();
        for coefficient in ["1", "z^2", "z^-1", "(z+w)^2"] {
            let word = Atom::parse(
                format!(
                    "chain(cof(3,i),dind(cof(3,j)),{coefficient}*t(coad(8,a),in,out),t(coad(8,b),in,out))"
                ),
                "spenso",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap();
            for preserve_indices in [false, true] {
                assert_eq!(
                    word.dirac_adjoint::<AbstractIndex>(preserve_indices)
                        .unwrap(),
                    word.spenso_conj(),
                );
            }
        }
        let closed = Atom::i()
            * gamma!([4, 11], [4, 12], spenso::mink!(4, 13))
            * gamma!([4, 12], [4, 11], spenso::mink!(4, 14));
        let expected =
            -Atom::num(4) * Atom::i() * spenso::g!(spenso::mink!(4, 13), spenso::mink!(4, 14));
        for preserve_indices in [false, true] {
            let adjoint = closed
                .dirac_adjoint::<AbstractIndex>(preserve_indices)
                .unwrap();
            let result = crate::tensor::SymbolicTensor::infer(adjoint)
                .unwrap()
                .simplify_gamma(Default::default())
                .unwrap()
                .resolved()
                .unwrap();
            assert_eq!(result.expression(), &expected);
        }
    }

    #[test]
    fn dirac_adjoint_preserves_physical_legs_across_different_pairings() {
        test_initialize();
        let direct = gamma!(1, 2, 5) * gamma!(3, 4, 5);
        let crossed = gamma!(1, 4, 6) * gamma!(3, 2, 6);
        let first = (Atom::one() + Atom::num(2) * Atom::i()) * &direct;
        let second = (Atom::num(3) - Atom::i()) * &crossed;
        let amplitude = &first + &second;
        let expected = (Atom::one() - Atom::num(2) * Atom::i()) * gamma!(2, 1, 5) * gamma!(4, 3, 5)
            + (Atom::num(3) + Atom::i()) * gamma!(4, 1, 6) * gamma!(2, 3, 6);
        // Domain collection gives internal contractions occurrence-local names.
        // Compare exact alpha normal forms while preserving every physical leg.
        let canonical = |value: &Atom| value.canonize(AbstractIndex::Dummy).unwrap();
        let adjoint = amplitude.dirac_adjoint::<AbstractIndex>(true).unwrap();
        assert_eq!(canonical(&adjoint), canonical(&expected));
        let separate = first.dirac_adjoint::<AbstractIndex>(true).unwrap()
            + second.dirac_adjoint::<AbstractIndex>(true).unwrap();
        assert_eq!(canonical(&adjoint), canonical(&separate));
        assert_eq!(
            canonical(&adjoint.dirac_adjoint::<AbstractIndex>(true).unwrap()),
            canonical(&amplitude)
        );
        assert_eq!(
            canonical(&direct.dirac_adjoint::<AbstractIndex>(false).unwrap()),
            canonical(&direct)
        );

        let malformed = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(Atom::var(symbol!("preserved_adjoint_operand")))
            .finish();
        assert!(matches!(
            malformed.dirac_adjoint::<AbstractIndex>(true),
            Err(AdjointError::Network(NetworkToolingError::Parse { .. }))
        ));
    }

    #[test]
    fn dirac_adjoint_lowers_brackets_and_compact_matrix_channels() {
        use crate::{representations::Bispinor, shorthands::chain::Chain};
        test_initialize();
        let p = SPENSO_TAG.rank_one_tensor_symbol("adjoint_compact_p");
        let explicit = gamma!(1, 2, 5)
            * gamma!(2, 3, 6)
            * FunctionBuilder::new(p)
                .add_arg(spenso::mink!(4, 6))
                .finish();
        let slashed = FunctionBuilder::new(crate::dirac::AGS.gamma)
            .add_arg(crate::bis!(4, 2))
            .add_arg(crate::bis!(4, 3))
            .add_arg(FunctionBuilder::new(p).add_arg(spenso::mink!(4)).finish())
            .finish();
        let compact = (gamma!(1, 2, 5) * slashed).chainify(Bispinor {}.into());
        let bracket = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&explicit)
            .finish();
        let canonical = |value: &Atom| value.canonize(AbstractIndex::Dummy).unwrap();
        let expected = explicit.dirac_adjoint::<AbstractIndex>(true).unwrap();
        for source in [bracket, compact] {
            let actual = source.dirac_adjoint::<AbstractIndex>(true).unwrap();
            assert_eq!(canonical(&actual), canonical(&expected));
            assert_eq!(
                canonical(&actual.dirac_adjoint::<AbstractIndex>(true).unwrap()),
                canonical(&explicit)
            );
        }
    }

    #[test]
    fn physical_adjoint_retains_scalar_and_foreign_factors() {
        test_initialize();
        let x = Atom::var(symbol!("adjoint_factored_x"; Real));
        let y = Atom::var(symbol!("adjoint_factored_y"; Real));
        let left = Atom::one() + x;
        let right = Atom::one() + y;
        let spectator = &left * &right;
        let foreign = spenso::tensor!(adjoint_foreign, crate::coad!(8, adjoint_color));
        let direct = gamma!(1, 2, 5) * gamma!(3, 4, 5);
        let crossed = gamma!(1, 4, 6) * gamma!(3, 2, 6);
        let amplitude = &spectator * &foreign * (direct + Atom::num(2) * crossed);
        let expected = spectator
            * foreign.spenso_conj()
            * (gamma!(2, 1, 5) * gamma!(4, 3, 5)
                + Atom::num(2) * gamma!(4, 1, 6) * gamma!(2, 3, 6));
        let actual = amplitude.dirac_adjoint::<AbstractIndex>(true).unwrap();
        // Neither spectator sum may be distributed or absorbed into a callback.
        assert!(actual.contains(&left) && actual.contains(&right));
        assert!(actual.contains(foreign.spenso_conj()));
        assert_eq!(
            actual.canonize(AbstractIndex::Dummy).unwrap(),
            expected.canonize(AbstractIndex::Dummy).unwrap()
        );
    }

    #[test]
    fn formal_adjoint_channels_specialize_to_four_dimensional_matrices() {
        use crate::shorthands::UndoShorthands;
        let reps = test_initialize();
        let expression = gamma!(
            reps.bis_d
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(99601)),
            reps.bis_d
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(99602)),
            reps.mink_d
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(99603))
        );
        // Standalone gamma factories remain restricted to their declared
        // matrix representation; the generic adjoint owns formal channels.
        assert!(crate::tensor::SymbolicTensor::infer(expression.clone()).is_err());
        let dimension = reps.bis_d.dim.to_symbolic();
        let specialize = |value: &Atom| {
            value
                .replace(dimension.to_pattern())
                .with(Atom::num(4).to_pattern())
        };
        let concrete = specialize(&expression);
        for preserve_indices in [false, true] {
            let result = expression
                .dirac_adjoint::<AbstractIndex>(preserve_indices)
                .unwrap();
            let expected = concrete
                .dirac_adjoint::<AbstractIndex>(preserve_indices)
                .unwrap();
            assert_eq!(
                specialize(&result)
                    .undo_chain::<AbstractIndex>()
                    .unwrap()
                    .canonize(AbstractIndex::Dummy)
                    .unwrap(),
                expected
                    .undo_chain::<AbstractIndex>()
                    .unwrap()
                    .canonize(AbstractIndex::Dummy)
                    .unwrap()
            );
        }
    }

    #[test]
    fn malformed_dirac_adjoint_returns_a_parse_error() {
        test_initialize();
        let malformed = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(Atom::var(symbol!("malformed_adjoint_operand")))
            .finish();

        let error = malformed
            .dirac_adjoint::<AbstractIndex>(false)
            .expect_err("one-argument dot notation should not have a Dirac adjoint");

        assert!(matches!(
            &error,
            AdjointError::Network(NetworkToolingError::Parse { reason })
                if reason.contains("Invalid dot function")
        ));
        assert!(error.to_string().contains("cannot parse tensor network"));
    }

    #[test]
    fn malformed_canonicalization_returns_a_parse_error() {
        test_initialize();
        let malformed = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(Atom::var(symbol!("malformed_canonicalization_operand")))
            .finish();

        let error = malformed
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .expect_err("one-argument dot notation should not canonicalize");

        assert!(matches!(
            &error,
            CanonicalizationError::Network(NetworkToolingError::Parse { reason })
                if reason.contains("Invalid dot function")
        ));
        assert!(error.to_string().contains("cannot parse tensor network"));
    }

    #[test]
    fn gamma_conj() {
        let a = test_initialize();
        let mink4 = a.mink4;

        let bis4 = a.bis4;

        let p1 = p!(1, slot!(mink4, 3));

        let ubgu = u!(2, slot!(bis4, 2))
            * (gamma!(1, 2, 3) * p1 + bis4.g(1, 2))
            * (u!(1, slot!(bis4, 1))
                .dirac_adjoint::<AbstractIndex>(false)
                .unwrap());

        assert_snapshot!(
            ubgu.dirac_adjoint::<AbstractIndex>(false)
                .unwrap()
                .canonize(AbstractIndex::Dummy)
                .expect("test expression should canonicalize")
                .to_canonically_ordered_string(
                    CanonicalOrderingSettings::new()
                        .include_attributes(false)
                        .include_namespace(false)
                ),@"(conj(p(1,mink(4,d_0)))*gamma(bis(4,d_1),bis(4,d_2),mink(4,d_0))+g(bis(4,d_1),bis(4,d_2)))*conj(u(2,bis(4,d_3)))*gamma0(bis(4,d_1),bis(4,d_3))*u(1,bis(4,d_2))"
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
        let can = can
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .expect("test expression should canonicalize");

        assert_snapshot!(
            can.to_bare_ordered_string(),@"-1*f(coad(8,d_0),coad(8,d_1),coad(8,d_2))*t(coad(8,d_0),cof(3,d_3),dind(cof(3,d_4)))*t(coad(8,d_1),cof(3,d_5),dind(cof(3,d_3)))*t(coad(8,d_2),cof(3,d_4),dind(cof(3,d_5)))"
        );

        //  gives
        //
        // f(coad(8,dummy(0)),coad(8,dummy(1)),coad(8,dummy(2)))*t(coad(8,dummy(0)),cof(3,dummy(3)),cof(3,hedge(2)))*t(coad(8,dummy(1)),cof(3,hedge(2)),cof(3,hedge(10)))*t(coad(8,dummy(2)),cof(3,hedge(10)),cof(3,dummy(3)))
    }
}
