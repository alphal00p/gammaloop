#![allow(uncommon_codepoints)]

use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::sync::atomic::{AtomicU64, Ordering};

use eyre::eyre;
use linnet::half_edge::subgraph::{BaseSubgraph, ModifySubSet, SuBitGraph, SubSetLike};
use linnet::tree::child_vec::ChildVecStore;
use shorthands::metric::{list_dangling_impl, wrap_dummies_impl, wrap_indices_impl};
use spenso::network::TensorNetworkError;
use spenso::network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp};
use spenso::network::library::function_lib::INBUILTS;
use spenso::network::library::symbolic::ETS;
use spenso::network::parsing::{
    ParseSettings, ParseState, SchoonschipExpansionMode, ShorthandParsing,
};
use spenso::network::store::TensorScalarStore;
use spenso::network::tags::SPENSO_TAG;
use spenso::network::tensor_symmetry::TensorSlotSymmetry;
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::representation::{Euclidean, LibraryRep, Lorentz, Minkowski, RepName};
use spenso::structure::slot::{
    AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind, SlotMatch, SlotMatcher,
};
use spenso::symbol_set;
use spenso::symbolica_init::in_symbolica_initializer;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol, SymbolBuilder};
use symbolica::id::{AliasedAtom, Pattern, Replacement};
use symbolica::printer::AtomPrinter;
use symbolica::{function, initialize, wrap_symbol};
use thiserror::Error;

use crate::color::{CS, color_conj_impl};
use crate::dirac::AGS;
use crate::epsilon::EPSILON_SYMBOL;
use crate::rep_symbols::RS;
use crate::representations::{
    Bispinor, ColorAdjoint, ColorFundamental, ColorSextet, SpinFundamental,
};
use crate::shorthands::UndoShorthands;
use crate::shorthands::chain::Chain;
use crate::tensor::inference::InterfaceInference;
use crate::tensor::{CollectionMode, SymbolicNetExt, SymbolicNetParse, SymbolicTensor};

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
        let _ = CS.force_in_initializer().cas;
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
    /// Unbound chain endpoints have no named index and produce an error.
    ///
    /// # Returns
    /// A `Vec<Atom>` where each `Atom` represents a dangling index, or an error when the
    /// expression cannot be parsed as a tensor network.
    fn list_dangling<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Vec<Atom>, IndexToolingError>;
}

impl IndexTooling for Atom {
    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        new_dummy: impl FnMut(usize) -> Aind,
    ) -> Result<Atom, CanonicalizationError> {
        self.as_view().canonize(new_dummy)
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
    #[error(
        "tensor words require an explicit adjoint policy; adjoint their operands and choose their order through the typed word constructor"
    )]
    ImplicitTensorWordAdjoint,
    #[error(transparent)]
    Composition(#[from] tensor::composition::TensorCompositionError),
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

impl tensor::SymbolicTensor<PartialStructure> {
    /// Prepare one admitted positive-power scope through the collection owner.
    /// The caller's native index type and operation-wide fresh-index reservation
    /// survive reversible encoding into the existing typed collection engine.
    fn prepare_canonical_power<Aind: DummyAind + ParseableAind + AbsInd>(
        value: AtomView<'_>,
        state: &ParseState<Aind>,
    ) -> Result<Atom, CanonicalizationError> {
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let result = (|| {
            let cooked = cooking.try_cook(value).map_err(|error| error.to_string())?;
            let mut slots = SlotMatcher::default();
            let mut written = HashSet::new();
            cooked.visitor(&mut |part| {
                if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(part) {
                    written.insert(slot.aind);
                    return false;
                }
                true
            });
            if written.is_empty() {
                return Ok(None);
            }
            let tensor = Self::infer(cooked.clone()).map_err(|error| error.to_string())?;
            let prepared = tensor
                .collect_with_map(
                    CollectionMode::Factored,
                    None,
                    |_| true,
                    |value, _, _| Ok(value),
                )
                .map_err(|error| error.to_string())?;
            if prepared.expression() == &cooked {
                return Ok(None);
            }
            let mut fresh = HashMap::new();
            let expression = prepared.expression().replace_map(|part, _, out| {
                if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(part)
                    && !written.contains(&slot.aind)
                {
                    let index = *fresh
                        .entry(slot.aind)
                        .or_insert_with(|| state.fresh_index());
                    **out = slot.rep().slot::<Aind, Aind>(index).to_atom();
                }
            });
            Ok::<_, String>(Some(cooking.uncook(expression.as_view())))
        })();
        result
            .map(|prepared| prepared.unwrap_or_else(|| value.to_owned()))
            .map_err(|reason| CanonicalizationError::Prepare { reason })
    }

    /// Dirac adjunction with the source's ordered logical axes. The external
    /// representations are dualized; zero retains those axes without parsing.
    /// As for axis permutation, unresolved axes acquire independent graph labels
    /// before the operation; their public logical order remains unchanged.
    pub fn dirac_adjoint(&self, preserve_indices: bool) -> Result<Self, AdjointError> {
        self.explicit_index_occurrences()?;
        let materialized;
        let source = if self.structure().open_positions().is_empty() {
            self
        } else {
            materialized = self.materialized()?;
            &materialized
        };
        let interface = PartialStructure::from_logical_slots(
            source
                .structure()
                .logical_slots()
                .into_iter()
                .map(|slot| slot.dual()),
        );
        if source.is_zero() {
            return Ok(Self::from_declared_parts(Atom::Zero, interface)?);
        }
        let (expression, has_spinor_edges) = Self::prepare_dirac_adjoint::<AbstractIndex>(
            source.expression().as_view(),
            preserve_indices,
        )?;
        let expression = if has_spinor_edges {
            expression.chainify(Bispinor {}.into())
        } else {
            expression
        };
        let value = Self::from_declared_parts(expression, interface)?;
        if has_spinor_edges {
            Ok(value.finish_dirac_adjoint(preserve_indices)?)
        } else {
            Ok(value)
        }
    }

    fn finish_dirac_adjoint(
        self,
        preserve_indices: bool,
    ) -> Result<Self, tensor::inference::TensorInferenceError> {
        if preserve_indices {
            self.simplify_algebra(&crate::tensor::AlgebraSettings {
                gamma: Some(dirac::GammaSimplifySettings {
                    gamma0: true,
                    evaluate_traces: false,
                    ..Default::default()
                }),
                epsilon: true,
                ..Default::default()
            })?
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
        } else {
            self.contract(Default::default())
        }
    }

    // The boolean records whether the prepared value needs Dirac cleanup.
    // Pure scalar/color conjugation keeps the raw entry point's existing form.
    fn prepare_dirac_adjoint<Aind: DummyAind + ParseableAind + AbsInd>(
        source: AtomView<'_>,
        preserve_indices: bool,
    ) -> Result<(Atom, bool), AdjointError> {
        let mut tensor_word = false;
        source.visitor(&mut |value| {
            if let AtomView::Fun(function) = value {
                if function.get_symbol().is_scalar() {
                    return false;
                }
                tensor_word |= function
                    .get_symbol()
                    .has_tag(&spenso::network::tags::SPENSO_TAG.tensor_word);
            }
            !tensor_word
        });
        if tensor_word {
            return Err(AdjointError::ImplicitTensorWordAdjoint);
        }
        let net = source
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
            return Ok((source.spenso_conj(), false));
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
        // Include the parser's newly materialized channel indices, then share
        // one collision-free allocator across boundary and matrix conjugation.
        let dummies = ParseState::<Aind>::default();
        dummies.reserve_indices(explicit.as_view());
        let mut a = explicit.spenso_conj();
        for i in bis_dangling {
            let mut dummy = i;
            dummy.aind = dummies.fresh_index();
            a = a
                .replace(Pattern::Literal(i.to_atom()))
                .with(Pattern::Literal(dummy.to_atom()))
                * function!(AGS.gamma0, i.to_atom(), dummy.to_atom());
        }
        for (i, j) in transpositions {
            a = a.replace_multiple(&[
                Replacement::new(Pattern::Literal(i.to_atom()), Pattern::Literal(j.to_atom())),
                Replacement::new(Pattern::Literal(j.to_atom()), Pattern::Literal(i.to_atom())),
            ]);
        }
        a = dirac::DiracSimplifier::conjugate_matrices::<Aind>(a.as_view(), &dummies);
        Ok((dirac::DiracSimplifier::factor_gamma_zero(a.as_view()), true))
    }
}

impl IndexTooling for AtomView<'_> {
    fn spenso_print(&self, settings: &SpensoPrintSettings) -> AtomPrinter<'_> {
        self.printer(settings.nice_symbolica())
    }
    fn canonize<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        mut new_dummy: impl FnMut(usize) -> Aind,
    ) -> Result<Atom, CanonicalizationError> {
        // Admit the complete source before algebraic zero pruning. The chosen
        // index parser owns accepted label syntax; canonicalization uses those
        // parsed identities without imposing another spelling restriction.
        // Retain each native scalar operand until its graph role has been
        // handled; a closed tensor word must not absorb its scalar coefficient.
        let settings = ParseSettings {
            precontract_scalars: false,
            ..ParseSettings::default()
        };
        let mut net = self
            .parse_to_symbolic_net::<Aind>(&settings)
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        // Native graph leaves own their private payloads. Alias scalar leaves
        // only at that boundary; identical syntax inside tensor metadata remains
        // untouched. Fresh registrations reject all existing global attributes.
        static NEXT_SCALAR_ALIAS: AtomicU64 = AtomicU64::new(0);
        let aliases = RefCell::new(BTreeMap::new());
        let scalar_alias = |original: Atom| {
            aliases
                .borrow_mut()
                .entry(original)
                .or_insert_with(|| {
                    loop {
                        let serial = NEXT_SCALAR_ALIAS.fetch_add(1, Ordering::Relaxed);
                        let builder = SymbolBuilder::new(wrap_symbol!(format!(
                            "idenso::canonical_scalar_{serial}"
                        )));
                        if let Ok(heads) = SymbolBuilder::build_group(vec![builder]) {
                            break Atom::var(heads[0]);
                        }
                    }
                })
                .clone()
        };
        let opaque_scalar = |value: AtomView<'_>| {
            !matches!(value, AtomView::Num(_) | AtomView::Var(_))
                && !InterfaceInference::has_structured_syntax(value)
        };
        let state = ParseState::<Aind>::default();
        state.reserve_indices(*self);
        let nested_dummy = RefCell::new(&mut new_dummy as &mut dyn FnMut(usize) -> Aind);
        let canonical_power = |base: AtomView<'_>, exponent: AtomView<'_>| {
            if !base
                .list_dangling::<Aind>()
                .map_err(|error| CanonicalizationError::Prepare {
                    reason: error.to_string(),
                })?
                .is_empty()
            {
                return Err(CanonicalizationError::Prepare {
                    reason: "a non-scalar tensor power requires a non-negative integer exponent"
                        .into(),
                });
            }
            let canonical = base.canonize::<Aind>(&mut **nested_dummy.borrow_mut())?;
            Ok::<_, CanonicalizationError>(
                SymbolicTensor::<PartialStructure>::scalar_power_expression(canonical, None)
                    .pow(exponent),
            )
        };
        let native_error = |error| TensorNetworkError::Other(eyre!("{error}"));
        let tree = net.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = net.graph.graph.node_id(net.graph.head());
        let prepared = net
            .graph
            .to_expression_at_with_operation(
                &tree,
                root,
                &mut |_, node| match node {
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                        Ok(Some(net.store.tensors[*index].expression.clone()))
                    }
                    NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                        let scalar = net.store.get_scalar_ref(*index);
                        let value = match scalar.as_view() {
                            AtomView::Num(_) | AtomView::Var(_) => scalar.clone(),
                            AtomView::Pow(power)
                                if !i64::try_from(power.get_exp()).is_ok_and(|n| n > 0) =>
                            {
                                let (base, exponent) = power.get_base_exp();
                                scalar_alias(canonical_power(base, exponent).map_err(native_error)?)
                            }
                            _ if InterfaceInference::has_structured_syntax(scalar.as_view()) => {
                                scalar.clone()
                            }
                            _ => scalar_alias(scalar.clone()),
                        };
                        Ok(Some(value))
                    }
                    NetworkNode::Op(_) => Ok(None),
                    _ => Err(TensorNetworkError::Other(eyre!(
                        "canonicalization requires stored symbolic leaves"
                    ))),
                },
                &mut |_, operation, children| {
                    let NetworkOp::Power(power) = operation else {
                        return Ok(None);
                    };
                    let [base] = children else {
                        unreachable!(
                            "native emitter checks unary arity before the operation callback"
                        )
                    };
                    if *power > 1 {
                        let expression =
                            SymbolicTensor::<PartialStructure>::scalar_power_expression(
                                base.clone(),
                                None,
                            )
                            .pow(*power);
                        return SymbolicTensor::prepare_canonical_power::<Aind>(
                            expression.as_view(),
                            &state,
                        )
                        .map(Some)
                        .map_err(native_error);
                    }
                    if *power < 0 {
                        let exponent = Atom::num(*power);
                        return canonical_power(base.as_view(), exponent.as_view())
                            .map(|power| Some(scalar_alias(power)))
                            .map_err(native_error);
                    }
                    Ok(None)
                },
            )
            .map_err(|error| CanonicalizationError::Prepare {
                reason: error.to_string(),
            })?;
        net = prepared
            .parse_to_symbolic_net::<Aind>(&settings)
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;
        net.remove_antisymmetric_zero_terms();

        // Keep duality wrappers in the expression graph. The canonicalizer
        // follows the same base index inside dind(...), while the wrapper
        // preserves which occurrence is covariant. Erasing that distinction
        // can collapse opposite color metrics into an invalid tensor square.
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
        // A variadic tensor word exposes only its free ports to the outer
        // network. Give Symbolica its internal slot inventory as well; its
        // canonicalizer owns deciding which of these slots are contracted.
        let mut slots = SlotMatcher::default();
        prepared.visitor(&mut |value| {
            let AtomView::Fun(function) = value else {
                return true;
            };
            if function.get_symbol().has_tag(&SPENSO_TAG.tensor_word) {
                value.visitor(&mut |part| {
                    if let SlotMatch::Explicit(literal) = slots.classify(part)
                        && let Ok(slot) = slots.parse::<LibraryRep, Aind>(part)
                    {
                        // Codecs can shorten equivalent index spellings, such as
                        // owner(1, 0) to owner(1). Symbolica's inventory matches
                        // literal slots, so retain the written payloads here.
                        let index = slot
                            .rep_name()
                            .base()
                            .to_symbolic([literal.dimension(), literal.index()]);
                        dummies.insert(index, slot.rep().base());
                        return false;
                    }
                    !matches!(part, AtomView::Fun(function) if function.get_symbol().is_scalar())
                });
                return false;
            }
            !function.get_symbol().is_scalar()
        });

        let external = net.graph.dangling_indices();
        // Reconstruct the admitted operation graph, without running tensor algebra
        // merely to hand its syntax to Symbolica's canonicalizer.
        let tree = net.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = net.graph.graph.node_id(net.graph.head());
        let expr = net
            .graph
            .to_expression_at(&tree, root, &mut |_, node| match node {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    let expression = &net.store.tensors[*index].expression;
                    let encoded = if let AtomView::Fun(function) = expression.as_view() {
                        TensorSlotSymmetry::canonical_view(function)
                            .map_err(|error| TensorNetworkError::Other(eyre!(error)))?
                    } else {
                        expression.clone()
                    };
                    Ok(Some(encoded))
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                    let scalar = net.store.get_scalar_ref(*index);
                    Ok(Some(if opaque_scalar(scalar.as_view()) {
                        scalar_alias(scalar.clone())
                    } else {
                        scalar.clone()
                    }))
                }
                NetworkNode::Op(_) => Ok(None),
                _ => Err(TensorNetworkError::Other(eyre!(
                    "canonicalization requires stored symbolic leaves"
                ))),
            })
            .map_err(|error| CanonicalizationError::Prepare {
                reason: error.to_string(),
            })?;
        let mut scalar_definitions = AliasedAtom::from(Atom::Zero);
        for (original, alias) in aliases.into_inner() {
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

        let mut reps = BTreeMap::new();
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
                let candidate = r.slot::<Aind, Aind>(new_dummy(next_dummy));
                next_dummy += 1;
                if reserved.insert(candidate.to_atom()) {
                    break candidate.aind();
                }
            };
            // The inventory retains literal spelling for Symbolica, while the
            // typed mapper uses parsed identities (owner(1, 0) and owner(1)
            // denote the same index). Key both sides through the slot parser.
            let source = slots
                .parse::<LibraryRep, Aind>(d.as_view())
                .map_err(|error| CanonicalizationError::Prepare {
                    reason: error.to_string(),
                })?;
            let source = if source.rep_name().is_dual() {
                source.dual()
            } else {
                source
            };
            reps.insert(source, target);
        }

        let relabeled = SlotMatcher::default().map_expression_indices::<Aind>(
            canonical.canonical_form.as_view(),
            |slot| {
                let slot = if slot.rep_name().is_dual() {
                    slot.dual()
                } else {
                    *slot
                };
                reps.get(&slot).copied()
            },
        );
        let canonical = TensorSlotSymmetry::restore_canonical_view(&relabeled)
            .map_err(|reason| CanonicalizationError::Prepare { reason })?;
        Ok(scalar_definitions.map_root(|_| canonical).into_inner())
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
        let (factored, has_spinor_edges) =
            tensor::SymbolicTensor::<PartialStructure>::prepare_dirac_adjoint::<Aind>(
                *self,
                preserve_indices,
            )?;
        if !has_spinor_edges {
            return Ok(factored);
        }
        // The generic index encoding is reversible; shared algebra sees the
        // same strict slot grammar as every other SymbolicTensor operation.
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
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
        let resolved = tensor.finish_dirac_adjoint(preserve_indices)?;
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
    use crate::tensor::SymbolicTensor;
    use crate::tensor::inference::InterfaceInference;
    use insta::assert_snapshot;
    use spenso::mink;
    use spenso::network::library::function_lib::INBUILTS;
    use spenso::structure::representation::{LibraryRep, RepName};
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::{
        network::tags::SPENSO_TAG,
        p, slot,
        structure::{abstract_index::AbstractIndex, slot::ParseableAind},
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
        function, parse_lit,
        printer::CanonicalOrderingSettings,
        symbol,
    };
    use symbolica_utils::AtomPrintExt;

    use crate::{
        AdjointError, CanonicalizationError, Cookable, IndexTooling, NetworkToolingError, gamma,
        test_support::test_initialize, u,
    };

    #[test]
    fn canonicalization_of_closed_scalar_power_matches_fresh_contracted_copies() {
        test_initialize();
        let a = SPENSO_TAG.tensor_symbol("canonical_power::a");
        let b = SPENSO_TAG.tensor_symbol("canonical_power::b");
        let index = mink!(4, 701);
        let norm = SymbolicTensor::infer(function!(a, &index) * function!(b, &index)).unwrap();
        for power in [2, 3] {
            let powered = norm
                .pow(&SymbolicTensor::infer(Atom::num(power)).unwrap())
                .unwrap();
            let mut explicit = SymbolicTensor::infer(Atom::one()).unwrap();
            for copy in 0..power {
                explicit = explicit
                    .multiply(
                        &norm
                            .wrap_indices(symbol!(&format!("canonical_power::copy_{copy}")), true)
                            .unwrap(),
                    )
                    .unwrap();
            }
            let expected = explicit.canonize().unwrap();
            let actual = powered.canonize().unwrap();
            assert_eq!(actual.expression(), expected.expression());
            assert_eq!(
                powered
                    .expanded(None, false)
                    .unwrap()
                    .canonize()
                    .unwrap()
                    .expression(),
                expected.expression()
            );
            assert_eq!(actual.canonize().unwrap().expression(), actual.expression());
            assert_ne!(actual.expression(), norm.canonize().unwrap().expression());
            assert_eq!(
                powered.expression().canonize(AbstractIndex::Dummy).unwrap(),
                *expected.expression()
            );
        }
    }

    #[test]
    fn canonicalization_of_complex_norm_power_in_scalar_density() {
        test_initialize();
        let h = SPENSO_TAG.tensor_symbol("canonical_complex_power::h");
        let rep = LibraryRep::new_dual("canonical_complex_power::weak")
            .unwrap()
            .new_rep(2);
        let scope = symbol!("canonical_complex_power::scope");
        let index = AbstractIndex::Normal(901).scoped(scope);
        let norm = SymbolicTensor::infer(
            INBUILTS.conj(function!(
                h,
                rep.dual().slot::<AbstractIndex, _>(index).to_atom()
            )) * function!(h, rep.slot::<AbstractIndex, _>(index).to_atom()),
        )
        .unwrap();
        let square = norm
            .pow(&SymbolicTensor::infer(Atom::num(2)).unwrap())
            .unwrap();
        let explicit = norm
            .wrap_indices(symbol!("canonical_complex_power::left"), true)
            .unwrap()
            .multiply(
                &norm
                    .wrap_indices(symbol!("canonical_complex_power::right"), true)
                    .unwrap(),
            )
            .unwrap();
        let coefficient = SymbolicTensor::infer(
            Atom::num((-1, 2)) * Atom::var(symbol!("canonical_complex_power::lambda"; Scalar)),
        )
        .unwrap();
        let density = square
            .multiply(&coefficient)
            .unwrap()
            .try_add(&norm)
            .unwrap();
        let expected = explicit
            .multiply(&coefficient)
            .unwrap()
            .try_add(&norm)
            .unwrap();
        let callback = symbol!("canonical_complex_power::callback", norm = |_, _| {});
        let wrapped_callback = INBUILTS.conj(function!(callback, &index.to_atom()));
        assert!(!InterfaceInference::normalization_is_intrinsic(
            wrapped_callback.as_view()
        ));
        assert_eq!(
            density.canonize().unwrap().expression(),
            expected.canonize().unwrap().expression()
        );
        assert_eq!(
            density
                .expanded(None, false)
                .unwrap()
                .canonize()
                .unwrap()
                .expression(),
            expected.canonize().unwrap().expression()
        );
    }

    #[test]
    fn canonicalization_of_word_internal_dummies_keeps_external_names_and_signs() {
        test_initialize();
        let head = symbol!("canonical_word::exterior"; Flat, Linear, Antisymmetric;
            tag = &SPENSO_TAG.tensor_word);
        let a = SPENSO_TAG.tensor_symbol("canonical_word::a");
        let b = SPENSO_TAG.tensor_symbol("canonical_word::b");
        let index = symbol!("canonical_word::owned_index", tag = &SPENSO_TAG.index);
        let density = |scope: usize, external: usize, reverse: bool| {
            let internal = function!(index, scope, 0);
            let left =
                SymbolicTensor::infer(function!(a, bis!(4, 71), cof!(3, &internal))).unwrap();
            let right =
                SymbolicTensor::infer(function!(b, bis!(4, 72), coaf!(3, &internal))).unwrap();
            let operands = if reverse {
                [right, left]
            } else {
                [left, right]
            };
            let word = SymbolicTensor::word(head, &operands).unwrap();
            word.multiply(
                &SymbolicTensor::infer(gamma!(bis!(4, 71), bis!(4, 72), mink!(4, external as u64)))
                    .unwrap(),
            )
            .unwrap()
        };
        let canonical = density(1, 73, false).canonize().unwrap();
        let renamed = density(2, 73, false).canonize().unwrap();
        assert_eq!(canonical.expression(), renamed.expression());
        assert_eq!(canonical.rank(), 1);
        assert_eq!(
            canonical.expression(),
            &-density(2, 73, true).canonize().unwrap().expression()
        );
        assert_ne!(
            canonical.expression(),
            density(2, 74, false).canonize().unwrap().expression()
        );
        assert_eq!(
            canonical.canonize().unwrap().expression(),
            canonical.expression()
        );
        assert!(
            !canonical
                .expression()
                .get_all_symbols(true)
                .contains(&index)
        );
    }

    #[test]
    fn canonicalization_of_powered_sums_keeps_factorization_and_outer_ports() {
        test_initialize();
        let index = mink!(4, 801);
        let pair = |left, right| {
            SymbolicTensor::infer(
                function!(SPENSO_TAG.tensor_symbol(left), &index)
                    * function!(SPENSO_TAG.tensor_symbol(right), &index),
            )
            .unwrap()
        };
        let sum = pair("canonical_sum_power::a", "canonical_sum_power::b")
            .try_add(&pair("canonical_sum_power::c", "canonical_sum_power::d"))
            .unwrap();
        let squared = sum
            .pow(&SymbolicTensor::infer(Atom::num(2)).unwrap())
            .unwrap();
        let explicit = sum
            .wrap_indices(symbol!("canonical_sum_power::left"), true)
            .unwrap()
            .multiply(
                &sum.wrap_indices(symbol!("canonical_sum_power::right"), true)
                    .unwrap(),
            )
            .unwrap();
        let free = mink!(4, AbstractIndex::Dummy(0).to_atom());
        let spectator = (Atom::var(symbol!("canonical_sum_power::x"; Scalar))
            + Atom::var(symbol!("canonical_sum_power::y"; Scalar)))
        .pow(7);
        let outer = SymbolicTensor::infer(
            function!(
                SPENSO_TAG.tensor_symbol("canonical_sum_power::outer"),
                &free
            ) * &spectator,
        )
        .unwrap();
        let actual = squared.multiply(&outer).unwrap().canonize().unwrap();
        let expected = explicit.multiply(&outer).unwrap().canonize().unwrap();
        assert_eq!(actual.expression(), expected.expression());
        assert_eq!(actual.rank(), 1);
        assert!(actual.expression().contains(&free));
        assert!(actual.expression().contains(&spectator));
        assert_ne!(actual.expression(), &actual.expression().expand());
        assert_eq!(actual.canonize().unwrap().expression(), actual.expression());
    }

    #[test]
    fn typed_adjoint_keeps_native_logical_axes_and_ranked_zero() {
        use crate::{CookMode, CookSettings, representations::Bispinor, tensor::SymbolicTensor};
        use spenso::{
            iterators::IteratableTensor,
            network::{
                ExecutionResult, Network, Sequential, SmallestDegree,
                library::{
                    panicing::ErroringLibrary,
                    symbolic::{ExplicitKey, TensorLibrary},
                },
                parsing::{ParseSettings, ShadowedStructure},
                store::NetworkStore,
            },
            structure::{
                Canonicalized, TensorStructure,
                partial::{PartialStructure, PartialStructureExt},
                representation::{LibraryRep, RepName},
                slot::DualSlotTo,
            },
            tensors::parametric::ParamTensor,
        };
        use std::collections::BTreeMap;
        type Tensor = ParamTensor<ShadowedStructure<AbstractIndex>>;
        type Lib = TensorLibrary<ParamTensor<ExplicitKey<AbstractIndex>>, AbstractIndex>;
        type Net = Network<
            NetworkStore<Tensor, Atom>,
            ExplicitKey<AbstractIndex>,
            symbolica::atom::Symbol,
        >;

        test_initialize();
        let matrix = spenso::tensor_symbol!("typed_adjoint::matrix"; Real);
        let label = symbol!("typed_adjoint::index", tags = [SPENSO_TAG.index.clone()]);
        let bis = LibraryRep::from(Bispinor {}).new_rep(4);
        let raw = function!(
            matrix,
            bis.to_symbolic([function!(label, 3)]),
            bis.to_symbolic([function!(label, 3, 1)])
        );
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let native = cooking.try_cook(raw.as_view()).unwrap();
        assert_eq!(native, raw, "valid native index labels are preserved");
        let inferred = SymbolicTensor::infer(native).unwrap();
        let mut slots = inferred.structure().logical_slots();
        slots.reverse();
        let source = SymbolicTensor::from_declared_parts(
            inferred.expression().clone(),
            PartialStructure::from_logical_slots(slots),
        )
        .unwrap();
        let interface = PartialStructure::from_logical_slots(
            source
                .structure()
                .logical_slots()
                .into_iter()
                .map(|slot| slot.dual()),
        );
        let mut library = Lib::new();
        library
            .insert_explicit_sparse(
                ExplicitKey::from_iter([bis, bis], matrix, None),
                (0usize..4).flat_map(|row| {
                    (0usize..4)
                        .map(move |column| (vec![row, column], Atom::num(1 + row * 4 + column)))
                }),
                Atom::Zero,
            )
            .unwrap();
        library
            .insert_explicit_sparse(
                ExplicitKey::from_iter([bis, bis], crate::dirac::AGS.gamma0, None),
                [(0, 2), (1, 3), (2, 0), (3, 1)]
                    .into_iter()
                    .map(|(row, column)| (vec![row, column], Atom::one())),
                Atom::Zero,
            )
            .unwrap();
        let components = |value: SymbolicTensor<PartialStructure>| {
            let lowered: Canonicalized<Net> = value
                .into_network(
                    &library,
                    &ErroringLibrary::new(),
                    &ParseSettings::default(),
                    None,
                )
                .unwrap();
            let (mut network, layout) = lowered.into_parts();
            network
                .execute::<Sequential, SmallestDegree, _, _, _>(&library, &ErroringLibrary::new())
                .unwrap();
            let ExecutionResult::Val(tensor) = network.result_tensor(&library).unwrap() else {
                panic!("an open adjoint, including zero, retains two axes");
            };
            assert_eq!(
                tensor.external_reps_iter().collect::<Vec<_>>(),
                vec![bis, bis]
            );
            let mut entries = (0..4)
                .flat_map(|row| (0..4).map(move |column| (vec![row, column], Atom::Zero)))
                .collect::<BTreeMap<_, _>>();
            for (indices, value) in tensor.iter_expanded() {
                entries.insert(layout.canonical_to_logical(&indices), value.to_owned());
            }
            entries
        };
        let mut nonzero = Vec::new();
        for preserve_indices in [false, true] {
            let expected = raw
                .dirac_adjoint::<AbstractIndex>(preserve_indices)
                .unwrap();
            let expected = SymbolicTensor::from_declared_parts(
                cooking.try_cook(expected.as_view()).unwrap(),
                interface.clone(),
            )
            .unwrap();
            let actual = source.dirac_adjoint(preserve_indices).unwrap();
            assert_eq!(actual.structure(), &interface);
            let entries = components(actual);
            assert_eq!(entries, components(expected));
            nonzero.push(entries);
            let zero = source
                .with_checked_expression(Atom::Zero)
                .unwrap()
                .dirac_adjoint(preserve_indices)
                .unwrap();
            assert_eq!(zero.structure(), &interface);
            assert!(components(zero).values().all(|value| value.is_zero()));
        }
        assert_ne!(
            nonzero[0], nonzero[1],
            "the matrix control must distinguish endpoint transposition"
        );

        let anonymous =
            SymbolicTensor::infer(function!(matrix, bis.to_symbolic([]), bis.to_symbolic([])))
                .unwrap();
        assert_eq!(anonymous.structure().open_positions().len(), 2);
        let named = anonymous.materialized().unwrap();
        let named_dual = PartialStructure::from_logical_slots(
            named
                .structure()
                .logical_slots()
                .into_iter()
                .map(|slot| slot.dual()),
        );
        for preserve_indices in [false, true] {
            let expected = named
                .expression()
                .dirac_adjoint::<AbstractIndex>(preserve_indices)
                .unwrap();
            let expected =
                SymbolicTensor::from_declared_parts(expected, named_dual.clone()).unwrap();
            let actual = anonymous.dirac_adjoint(preserve_indices).unwrap();
            assert_eq!(actual.rank(), 2);
            assert_eq!(components(actual), components(expected));
            let zero = anonymous
                .with_checked_expression(Atom::Zero)
                .unwrap()
                .dirac_adjoint(preserve_indices)
                .unwrap();
            assert_eq!(zero.rank(), 2);
            assert!(components(zero).values().all(|value| value.is_zero()));
        }
        assert_eq!(anonymous.structure().open_positions().len(), 2);

        for raw in [
            Atom::i() * Atom::num(2),
            parse_lit!(spenso::g(
                spenso::cof(3, 1),
                spenso::dind(spenso::cof(3, 2))
            )),
        ] {
            let source = SymbolicTensor::infer(raw.clone()).unwrap();
            for preserve_indices in [false, true] {
                let actual = source.dirac_adjoint(preserve_indices).unwrap();
                assert_eq!(
                    actual.expression(),
                    &raw.dirac_adjoint::<AbstractIndex>(preserve_indices)
                        .unwrap()
                );
                assert_eq!(
                    actual.structure().logical_slots(),
                    source
                        .structure()
                        .logical_slots()
                        .into_iter()
                        .map(|slot| slot.dual())
                        .collect::<Vec<_>>()
                );
            }
        }
    }

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
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    gamma: Some(Default::default()),
                    epsilon: true,
                    ..Default::default()
                })
                .unwrap()
                .contract(crate::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap();
            assert_eq!(result.expression(), &expected);
        }
    }

    #[test]
    fn canonicalization_keeps_requested_dummy_names_literal() {
        test_initialize();
        let left = spenso::tensor_symbol!("canonical_literal_left");
        let right = spenso::tensor_symbol!("canonical_literal_right");
        let source = spenso::mink!(4, Atom::var(symbol!("canonical_input_")));
        let target_index = AbstractIndex::Symbol(symbol!("canonical_output_").into());
        let target = spenso::mink!(4, target_index.to_atom());
        let expression = function!(left, &source) * function!(right, &source);
        assert_eq!(
            expression.canonize(|_| target_index).unwrap(),
            function!(left, &target) * function!(right, &target)
        );
    }

    #[test]
    fn dirac_adjoint_treats_admitted_wildcard_spelling_as_literal_labels() {
        test_initialize();
        let original = [
            crate::bis!(4, 97301),
            crate::bis!(4, 97302),
            spenso::mink!(4, 97303),
        ];
        let literal = [
            crate::bis!(4, Atom::var(symbol!("adjoint_literal_i_"))),
            crate::bis!(4, Atom::var(symbol!("adjoint_literal_j_"))),
            spenso::mink!(4, Atom::var(symbol!("adjoint_literal_mu_"))),
        ];
        let relabel = |expression: &Atom| {
            expression.replace_map(|node, _, out| {
                if let Some(position) = original.iter().position(|slot| slot.as_view() == node) {
                    **out = literal[position].clone();
                }
            })
        };
        let source = gamma!(&original[0], &original[1], &original[2]);
        let labeled = relabel(&source);
        crate::tensor::SymbolicTensor::infer(labeled.clone()).unwrap();
        for preserve in [true, false] {
            let expected = relabel(&source.dirac_adjoint::<AbstractIndex>(preserve).unwrap());
            let actual = labeled.dirac_adjoint::<AbstractIndex>(preserve).unwrap();
            assert_eq!(
                actual.canonize(AbstractIndex::Dummy).unwrap(),
                expected.canonize(AbstractIndex::Dummy).unwrap()
            );
            let actual = crate::tensor::SymbolicTensor::infer(actual).unwrap();
            let expected = crate::tensor::SymbolicTensor::infer(expected).unwrap();
            assert_eq!(actual.structure(), expected.structure());
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
            can.to_bare_ordered_string(),@"-1*f(coad(8,d_0),coad(8,d_1),coad(8,d_2))*t(coad(8,d_0),cof(3,d_4),dind(cof(3,d_3)))*t(coad(8,d_1),cof(3,d_5),dind(cof(3,d_4)))*t(coad(8,d_2),cof(3,d_3),dind(cof(3,d_5)))"
        );

        //  gives
        //
        // f(coad(8,dummy(0)),coad(8,dummy(1)),coad(8,dummy(2)))*t(coad(8,dummy(0)),cof(3,dummy(3)),cof(3,hedge(2)))*t(coad(8,dummy(1)),cof(3,hedge(2)),cof(3,hedge(10)))*t(coad(8,dummy(2)),cof(3,hedge(10)),cof(3,dummy(3)))
    }
}
