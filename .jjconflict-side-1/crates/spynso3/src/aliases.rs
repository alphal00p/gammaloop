//! Python conversion and result wrapping for Idenso's typed alias registry.

use std::{collections::HashMap, sync::Arc};

use idenso::tensor::{SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError};
use pyo3::prelude::*;
use spenso::structure::partial::PartialStructure;
use symbolica::{
    api::python::{PythonExpression, PythonExpressionEvaluator},
    atom::{AliasedAtom, Atom, Symbol},
};

use crate::{ModuleInit, expression::TensorExpression};

type Descriptor = (Option<Symbol>, Vec<Atom>);

/// A symbolic tensor whose repeated subexpressions have named definitions.
///
/// Tensor simplification uses this class to retain shared work instead of
/// expanding large intermediate expressions. ``root`` uses those names and
/// ``aliases`` supplies their definitions. Scalar results use the same class.
/// Call ``to_expression()`` to resolve definitions into a TensorExpression,
/// or ``expand()`` to resolve them and distribute the scalar algebra.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import AliasedTensorExpression
/// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
/// >>> aliases.to_expression().rank
/// 2
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "AliasedTensorExpression",
    module = "symbolica.community.tensor"
)]
pub struct AliasedTensorExpression {
    pub(crate) value: Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>,
    descriptor: Descriptor,
    descriptors: Arc<HashMap<Atom, (Descriptor, Descriptor)>>,
}

impl ModuleInit for AliasedTensorExpression {}

impl AliasedTensorExpression {
    pub(crate) fn from_parts(
        value: Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>,
        name: Option<Symbol>,
        arguments: Vec<Atom>,
    ) -> Self {
        Self {
            value,
            descriptor: (name, arguments),
            descriptors: Arc::new(HashMap::new()),
        }
    }

    fn with_value(&self, value: Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>) -> Self {
        Self {
            value,
            descriptor: self.descriptor.clone(),
            descriptors: Arc::clone(&self.descriptors),
        }
    }

    fn wrap(
        py: Python<'_>,
        tensor: SymbolicTensor<PartialStructure>,
        descriptor: &Descriptor,
    ) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_shared(py, Arc::new(tensor), descriptor.0, descriptor.1.clone())
    }

    fn descriptor(value: &PyRef<'_, TensorExpression>) -> Descriptor {
        (value.name, value.name_args.clone())
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl AliasedTensorExpression {
    /// Build a tensor expression with explicitly supplied shared definitions.
    ///
    /// Parameters
    /// ----------
    /// root : TensorExpression
    ///     Expression containing the alias handles.
    /// aliases : sequence of (TensorExpression, TensorExpression), default []
    ///     Opaque handle and definition pairs, as returned by from_expression().
    ///     A scalar alias may use a scalar variable. Each pair must have compatible tensor
    ///     interfaces and the definitions must form a consistent dependency graph.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     The root and shared definitions, without eager substitution.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> original = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> value = AliasedTensorExpression(original.root, original.aliases)
    /// >>> value.to_expression().rank
    /// 2
    #[new]
    #[pyo3(signature = (root, aliases=Vec::new()))]
    fn new(
        py: Python<'_>,
        root: PyRef<'_, TensorExpression>,
        aliases: Vec<(Py<TensorExpression>, Py<TensorExpression>)>,
    ) -> PyResult<Self> {
        let mut descriptors = HashMap::new();
        let aliases = aliases
            .into_iter()
            .map(|(handle, body)| {
                let handle = handle.borrow(py);
                let body = body.borrow(py);
                descriptors.insert(
                    handle.atom().clone(),
                    (Self::descriptor(&handle), Self::descriptor(&body)),
                );
                (
                    TensorExpression::structured(&handle).clone(),
                    TensorExpression::structured(&body).clone(),
                )
            })
            .collect::<Vec<_>>();
        Ok(Self {
            value: Arc::new(
                TensorExpression::structured(&root)
                    .clone()
                    .with_aliases(aliases)
                    .map_err(TensorExpression::inference_error)?,
            ),
            descriptor: Self::descriptor(&root),
            descriptors: Arc::new(descriptors),
        })
    }

    /// Name an entire tensor expression with one fresh alias.
    ///
    /// Parameters
    /// ----------
    /// expression : TensorExpression
    ///     Tensor to retain as a shared definition.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     A fresh handle as the root and one definition containing expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> len(aliases.aliases)
    /// 1
    #[staticmethod]
    fn from_expression(expression: PyRef<'_, TensorExpression>) -> PyResult<Self> {
        let body = TensorExpression::structured(&expression).clone();
        let handle = body
            .alias_handle()
            .map_err(TensorExpression::inference_error)?;
        let descriptor = Self::descriptor(&expression);
        let descriptors = HashMap::from([(
            handle.expression().clone(),
            ((None, Vec::new()), descriptor.clone()),
        )]);
        Ok(Self {
            value: Arc::new(
                handle
                    .clone()
                    .with_aliases([(handle, body)])
                    .map_err(TensorExpression::inference_error)?,
            ),
            descriptor,
            descriptors: Arc::new(descriptors),
        })
    }

    /// The tensor expression before substituting its shared definitions.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Root containing alias handles.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> aliases.root.rank
    /// 2
    #[getter]
    fn root(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(py, self.value.root(), &self.descriptor)
    }

    /// Whether symbolic metric/vector contraction is certified complete.
    ///
    /// Returns
    /// -------
    /// bool
    ///     False also covers a value not yet contracted, or a contraction stopped
    ///     before completion. It does not indicate whether aliases are resolved.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> complete = aliases.contraction_complete
    #[getter]
    fn contraction_complete(&self) -> bool {
        self.value.contraction_complete()
    }

    /// The tensor handles and their shared definitions.
    ///
    /// Returns
    /// -------
    /// list of (TensorExpression, TensorExpression)
    ///     Each pair maps a handle to its definition.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> len(aliases.aliases)
    /// 1
    #[getter]
    fn aliases(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<(Py<TensorExpression>, Py<TensorExpression>)>> {
        self.value
            .aliases()
            .map_err(TensorExpression::inference_error)?
            .into_iter()
            .map(|(handle, body)| {
                let empty = ((None, Vec::new()), (None, Vec::new()));
                let (handle_descriptor, body_descriptor) =
                    self.descriptors.get(handle.expression()).unwrap_or(&empty);
                Ok((
                    Self::wrap(py, handle, handle_descriptor)?,
                    Self::wrap(py, body, body_descriptor)?,
                ))
            })
            .collect()
    }

    /// Substitute the shared definitions into the root.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Materialized tensor algebra; this is not a raw Symbolica Expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> aliases.to_expression().rank
    /// 2
    fn to_expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(
            py,
            self.value
                .resolved()
                .map_err(TensorExpression::inference_error)?,
            &self.descriptor,
        )
    }

    /// Resolve aliases and expand the resulting scalar algebra.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Expanded tensor algebra. This can be much larger than the shared form.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> aliases.expand().rank
    /// 2
    fn expand(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        Self::wrap(
            py,
            self.value
                .expanded()
                .map_err(TensorExpression::inference_error)?,
            &self.descriptor,
        )
    }

    /// Measure the stored symbolic payload of the root and definitions.
    ///
    /// Returns
    /// -------
    /// int
    ///     Byte count of the symbolic payload; not the full Python allocation
    ///     or total process memory.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> aliases.get_byte_size() > 0
    /// True
    fn get_byte_size(&self) -> usize {
        self.value.expression().get_byte_size()
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> value = sp.AliasedTensorExpression.from_expression(sp.TensorExpression((x + 1)**2))
    /// >>> text = repr(value)
    fn __repr__(&self) -> String {
        format!(
            "AliasedTensorExpression(root={}, definitions={}, bytes={})",
            crate::display::format_structured(&self.value.root(), false),
            self.value.expression().get_aliases().len(),
            self.get_byte_size(),
        )
    }

    /// Transform each shared definition while retaining its tensor interface.
    ///
    /// Parameters
    /// ----------
    /// callback : callable
    ///     Receives a definition as TensorExpression and returns a
    ///     TensorExpression with a compatible interface.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     Updated definitions under the same handles; the original is unchanged.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> mapped = aliases.map_aliases(lambda body: body.expand())
    /// >>> mapped.to_expression().rank
    /// 2
    fn map_aliases(&self, py: Python<'_>, callback: Py<PyAny>) -> PyResult<Self> {
        let mut callback_error = None;
        let mut descriptors = (*self.descriptors).clone();
        let value = self.value.map_aliases(|handle, body| {
            let result = (|| {
                let empty = ((None, Vec::new()), (None, Vec::new()));
                let descriptor = self.descriptors.get(handle.expression()).unwrap_or(&empty);
                let argument = Self::wrap(py, body, &descriptor.1)?;
                let output = callback.call1(py, (argument,))?;
                let tensor = output.extract::<PyRef<'_, TensorExpression>>(py)?;
                descriptors.insert(
                    handle.expression().clone(),
                    (descriptor.0.clone(), Self::descriptor(&tensor)),
                );
                Ok::<_, PyErr>(TensorExpression::structured(&tensor).clone())
            })();
            result.map_err(|error| {
                callback_error = Some(error);
                TensorInferenceError::Invalid("alias callback failed".into())
            })
        });
        if let Some(error) = callback_error {
            return Err(error);
        }
        Ok(Self {
            value: Arc::new(value.map_err(TensorExpression::inference_error)?),
            descriptor: self.descriptor.clone(),
            descriptors: Arc::new(descriptors),
        })
    }

    /// Apply one rule or an ordered list of whole-tensor replacement rules.
    ///
    /// Parameters
    /// ----------
    /// rules : TensorRule or list of TensorRule
    ///     Rules in priority order. The first match at a source tensor is used;
    ///     newly produced replacements are not matched again during this call.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     Rewritten value with a checked external interface and shared definitions.
    ///
    /// Notes
    /// -----
    /// The traversal follows tensor sums and products, leaving scalar key
    /// arguments opaque. Matched tensors need explicit axes; unsupported tensor
    /// powers and unresolved interfaces are rejected. Normalization/evaluation
    /// callbacks remaining in a replacement are not supported. For scalar
    /// component substitutions in a network, use TensorNetwork.replace instead.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import TensorRule
    /// >>> indexed = A("i", "j")
    /// >>> rule = TensorRule(indexed, 2 * indexed)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(indexed)
    /// >>> replaced = aliases.replace(rule)
    #[pyo3(signature=(rules))]
    fn replace(&self, rules: &Bound<'_, PyAny>) -> PyResult<Self> {
        let rules = crate::tensor_rule::PyTensorRule::extract_rules(rules)?;
        let value = self
            .value
            .replace_rules(&rules)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Report completion under the last reduction settings and budget.
    ///
    /// Returns
    /// -------
    /// ReductionStatus
    ///     Complete, Deferred, or Capped. Unfinished results remain exact and can be
    ///     passed to contract() or simplify_algebra() again.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, ReductionStatus
    /// >>> TensorExpression(3).contract().reduction_status == ReductionStatus.Complete
    /// True
    #[getter]
    fn reduction_status(&self) -> crate::simplification::PyReductionStatus {
        self.value.reduction_status().into()
    }

    /// Contract compatible indices and collect compact matrix products.
    ///
    /// Parameters
    /// ----------
    /// representations : sequence of Representation or RepresentationName, optional
    ///     Filter connecting axes by representation name: None permits all, [] none.
    ///     Compatible dimensions and orientations are still required. Intrinsic
    ///     constructor normalization, including symmetric dots, is unaffected.
    /// metrics : bool, default True
    ///     Contract compatible metric tensors.
    /// rank_one : bool, default True
    ///     Contract compatible rank-one tensors into scalar products.
    /// collect_chains : bool, default True
    ///     Collect compatible matrix products into ordered chains.
    /// collect_traces : bool, default True
    ///     Collect closed chains without evaluating traces. False preserves indexed closure.
    /// order : sequence of int, optional
    ///     Requested factor contraction order; None uses the planner's default.
    /// max_steps_per_domain : int, default 16
    ///     Maximum successful transformations per root or alias definition; zero is allowed.
    ///     Inspect reduction_status for exact unfinished work.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     An exact result retaining shared definitions. Inspect reduction_status for
    ///     unfinished work; to_expression() resolves definitions without distributing sums.
    ///
    /// Notes
    /// -----
    /// This performs symbolic structural contraction. Use simplify_algebra() for Dirac,
    /// color, or Levi-Civita identities, or to_tensor() to evaluate stored components.
    /// Polynomial expansion remains an explicit operation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(3)
    /// >>> p = TensorName.vector("p")(space)
    /// >>> expression = TensorExpression.g(space)("i", "j") * p("j")
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> expression = AliasedTensorExpression.from_expression(expression)
    /// >>> result = expression.contract()
    /// >>> result.to_expression().rank
    /// 1
    #[pyo3(signature=(*, representations=None, metrics=true, rank_one=true, collect_chains=true, collect_traces=true, order=None, max_steps_per_domain=16))]
    #[allow(clippy::too_many_arguments)]
    fn contract(
        &self,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Sequence[Representation | RepresentationName]]", imports=("typing")))]
        representations: Option<Vec<Bound<'_, PyAny>>>,
        metrics: bool,
        rank_one: bool,
        collect_chains: bool,
        collect_traces: bool,
        order: Option<Vec<usize>>,
        max_steps_per_domain: usize,
    ) -> PyResult<Self> {
        let representations = crate::simplification::contraction_representations(representations)?;
        let settings = idenso::tensor::ContractSettings {
            representations: representations.as_deref(),
            metrics,
            rank_one,
            collect_chains,
            collect_traces,
            order: order.as_deref(),
            max_passes: max_steps_per_domain,
        };
        let value = self
            .value
            .contract(settings)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Write compact metric products using dot notation.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     An equivalent value in the requested notation; aliases remain shared.
    ///
    /// Notes
    /// -----
    /// This changes notation; it does not itself contract explicit vector indices.
    /// Use contract() first when those contractions are needed.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> converted = aliases.to_dots()
    fn to_dots(&self) -> PyResult<Self> {
        self.value
            .to_dots()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Write dot products as explicit indexed contractions.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     An equivalent value in the requested notation; aliases remain shared.
    ///
    /// Notes
    /// -----
    /// Other compact forms and factored sums are retained. This does not
    /// evaluate component data or expand a polynomial.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> aliases = AliasedTensorExpression.from_expression(A("i", "j"))
    /// >>> converted = aliases.undo_dots()
    fn undo_dots(&self) -> PyResult<Self> {
        self.value
            .undo_dots()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Expose collected chain factors with fresh compatible dummy indices.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     Equivalent algebra in the requested notation, retaining shared definitions
    ///     when called on an AliasedTensorExpression.
    ///
    /// Notes
    /// -----
    /// This changes notation without distributing polynomial sums or requesting
    /// Dirac/color identities. Use undo_trace().undo_chain() to expose indexed trace factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, chain, trace
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> expression = chain(space("i"), space("j"), A, A)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> expression = AliasedTensorExpression.from_expression(expression)
    /// >>> opened = expression.undo_chain()
    fn undo_chain(&self) -> PyResult<Self> {
        self.value
            .undo_chain()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Open compact traces as closed chains without evaluating them.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     Equivalent algebra in the requested notation, retaining shared definitions
    ///     when called on an AliasedTensorExpression.
    ///
    /// Notes
    /// -----
    /// This changes notation without distributing polynomial sums or requesting
    /// Dirac/color identities. Use undo_trace().undo_chain() to expose indexed trace factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, chain, trace
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> expression = trace(space, A, A)
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> expression = AliasedTensorExpression.from_expression(expression)
    /// >>> opened = expression.undo_trace()
    fn undo_trace(&self) -> PyResult<Self> {
        self.value
            .undo_trace()
            .map(|value| self.with_value(value))
            .map_err(TensorExpression::inference_error)
    }

    /// Apply selected tensor identities while retaining shared subexpressions.
    ///
    /// Parameters
    /// ----------
    /// gamma, color, epsilon : bool, default False
    ///     Enable only the selected identity families and their prerequisite contractions.
    /// gamma_output : {"reduced", "chains"}, optional
    ///     Apply Clifford identities (default), or only collect Dirac chains.
    /// gamma_ordering : {"repeated_pairs", "canonical"}, optional
    ///     Open-chain ordering; defaults to repeated-pair reduction.
    /// gamma_evaluate_traces : bool, optional
    ///     Evaluate closed Dirac traces; defaults to True for reduced output.
    /// gamma0 : bool, optional
    ///     Reduce four-dimensional gamma-zero sandwiches; defaults to False.
    /// gamma_conjugate : bool, optional
    ///     Rewrite complex-conjugated Dirac matrices; defaults to False.
    /// gamma_expand_three_gamma_epsilon : bool, optional
    ///     Expand three four-dimensional gammas in a gamma5/epsilon basis; defaults
    ///     to False. This does not enable epsilon identities.
    /// color_evaluate_traces : bool, optional
    ///     Evaluate supported closed color traces; defaults to True.
    /// color_expand_fierz : bool, optional
    ///     Apply cross-chain Fierz identities, possibly introducing sums; defaults to True.
    /// color_substitute_cof_dimension_invariants : bool, optional
    ///     Substitute fundamental SU(N) invariants with T_R=1/2; defaults to False.
    /// max_steps_per_domain : int, default 16
    ///     Maximum successful transformations per root or alias definition; zero is allowed.
    ///     Capped results retain exact unfinished work. Family-specific options require
    ///     the corresponding family to be enabled. None selects its native default.
    ///
    /// Returns
    /// -------
    /// AliasedTensorExpression
    ///     Exact reduced algebra with shared definitions and a reduction_status.
    ///     to_expression() resolves definitions; expand() additionally distributes sums.
    ///
    /// Notes
    /// -----
    /// Each enabled family includes its necessary contractions. Brackets are transparent
    /// and global polynomial expansion is not required. Use contract() for structural
    /// contractions independent of an algebra family. A capped result retains exact
    /// unfinished work and can be reduced again.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> gamma = TensorExpression.dirac_gamma(4)
    /// >>> expression = gamma("a", "b", "mu") * gamma("b", "a", "nu")
    /// >>> from symbolica.community.tensor import AliasedTensorExpression
    /// >>> expression = AliasedTensorExpression.from_expression(expression)
    /// >>> result = expression.simplify_algebra(gamma=True)
    /// >>> result.to_expression().rank
    /// 2
    #[pyo3(signature=(*, gamma=false, color=false, epsilon=false, gamma_output=None, gamma_ordering=None, gamma_evaluate_traces=None, gamma0=None, gamma_conjugate=None, gamma_expand_three_gamma_epsilon=None, color_evaluate_traces=None, color_expand_fierz=None, color_substitute_cof_dimension_invariants=None, max_steps_per_domain=16))]
    #[allow(clippy::too_many_arguments)]
    fn simplify_algebra(
        &self,
        gamma: bool,
        color: bool,
        epsilon: bool,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Literal['reduced', 'chains']]", imports=("typing")))]
        gamma_output: Option<&str>,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Literal['repeated_pairs', 'canonical']]", imports=("typing")))]
        gamma_ordering: Option<&str>,
        gamma_evaluate_traces: Option<bool>,
        gamma0: Option<bool>,
        gamma_conjugate: Option<bool>,
        gamma_expand_three_gamma_epsilon: Option<bool>,
        color_evaluate_traces: Option<bool>,
        color_expand_fierz: Option<bool>,
        color_substitute_cof_dimension_invariants: Option<bool>,
        max_steps_per_domain: usize,
    ) -> PyResult<Self> {
        let settings = crate::simplification::algebra_settings(
            gamma,
            color,
            epsilon,
            gamma_output,
            gamma_ordering,
            gamma_evaluate_traces,
            gamma0,
            gamma_conjugate,
            gamma_expand_three_gamma_epsilon,
            color_evaluate_traces,
            color_expand_fierz,
            color_substitute_cof_dimension_invariants,
            max_steps_per_domain,
        )?;
        let value = self
            .value
            .simplify_algebra(&settings)
            .map_err(TensorExpression::inference_error)?;
        Ok(self.with_value(value))
    }

    /// Prepare numerical evaluation of the symbolic algebra without component expansion.
    ///
    /// Parameters
    /// ----------
    /// params : sequence of Expression
    ///     Symbolic leaves to treat as numerical inputs, in input-row order.
    /// iterations : int, default 1
    ///     Horner-optimization iterations.
    /// n_cores : int, default 1
    ///     Number of optimization cores.
    ///
    /// Returns
    /// -------
    /// symbolica.Evaluator
    ///     Evaluator for the symbolic expression and its shared definitions.
    ///
    /// Notes
    /// -----
    /// This does not resolve tensor components from a TensorLibrary. To evaluate
    /// component formulas, first use to_tensor() and then Tensor.evaluator().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorExpression, AliasedTensorExpression
    /// >>> x = S("x")
    /// >>> expression = AliasedTensorExpression.from_expression(TensorExpression(x**2 + 1))
    /// >>> evaluator = expression.evaluator([x])
    #[pyo3(signature = (params, *, iterations=1, n_cores=1))]
    fn evaluator(
        &self,
        py: Python<'_>,
        params: Vec<PythonExpression>,
        iterations: usize,
        n_cores: usize,
    ) -> PyResult<PythonExpressionEvaluator> {
        let params = params
            .into_iter()
            .map(|value| value.expr)
            .collect::<Vec<_>>();
        let builder = self
            .value
            .evaluator(&params)
            .map_err(TensorExpression::inference_error)?;
        TensorExpression::build_evaluator(py, builder, iterations, n_cores)
    }
}
