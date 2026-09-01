use crate::color::ColorSimplifier;
use crate::dirac::GammaSimplifier;
use crate::selective_expand::SelectiveExpand;
use crate::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};

use crate::{Cookable, IndexTooling};
use pyo3::{
    Bound, PyResult, Python,
    exceptions::PyValueError,
    pyfunction,
    types::{PyAnyMethods, PyDictMethods, PyListMethods, PyModule, PyModuleMethods},
    wrap_pyfunction,
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pyfunction;
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::atom::Symbol;

use symbolica::api::python::PythonExpression;

mod algebra;
mod expansion;
mod tooling;

pub use algebra::{
    GammaConjugationError, PyColorCasimirSettings, PyColorSimplifySettings, PyGammaChainOrdering,
    PyGammaSimplifySettings, collect_color, collect_color_constants, collect_gamma_chains,
    simplify_epsilon, simplify_gamma_conjugate, simplify_gamma0, to_cof_dimension_invariants,
    to_color_casimir, wrap_color,
};
pub use expansion::{PythonTerm, expand_in_patterns};
pub use tooling::{
    CanonicalizationError, CookingError, DiracAdjointError, DotExpansionError, NetworkToolingError,
    PyCookMode, PyCookSettings, PyCookSourceFilter, PyCookTagFilter, PySchoonschipContractionOrder,
    PySchoonschipMode, PySchoonschipSettings, PySchoonschipTraversal, RegisteredRepresentation,
    alias_subtensors, canonize, chainify, collect_chains, conjugate_transpose, cook, expand_dots,
    metric_shorthand_to_dot, normalize_chains, normalize_dots, schoonschip, schoonschip_net,
    spenso_conjugate, uncook, undo_all, undo_chain, undo_dots, undo_schoonschip,
    undo_single_length, undo_trace,
};

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Construct the physics-aware Dirac adjoint of a tensor expression.
///
/// Idenso takes the symbolic complex conjugate, reverses compatible open bispinor chains, and
/// inserts the registered `gamma0` factors required at dangling bispinor slots.
/// The input must use the representation-aware Spenso forms registered on import.
/// Raises `DiracAdjointError` when the tensor network does not define a consistent adjoint.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import dirac_adjoint, list_dangling
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> # Built-in representations are registered automatically on import.
/// >>> # Re-importing the module does not require explicit re-registration.
/// >>> bispinor = Representation.bis(4)
/// >>> spinor = TensorName("u")(bispinor("alpha")).to_expression()
/// >>> adjoint = dirac_adjoint(spinor)
/// >>> len(list_dangling(adjoint)) == 1
/// True
/// >>> "gamma0" in str(adjoint)
/// True
/// ```
///
/// # Arguments
/// - `expression`: a Spenso-compatible tensor expression.
///
/// # Returns
/// The representation-aware Dirac adjoint.
pub fn dirac_adjoint(expression: &PythonExpression) -> PyResult<PythonExpression> {
    expression
        .expr
        .dirac_adjoint::<AbstractIndex>()
        .map(Into::into)
        .map_err(|error| tooling::DiracAdjointError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Expand products around factors carrying registered Minkowski indices.
///
/// This is a selective symbolic expansion: Minkowski-bearing factors become polynomial
/// variables while unrelated sectors remain coefficients. It does not substitute explicit
/// four-vector components or choose a metric signature.
///
/// # Arguments
/// - `expression`: a factorized Spenso-compatible expression.
///
/// # Returns
/// `(structure, coefficient)` pairs distributed around Minkowski-bearing factors.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_mink
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> # Built-in representations are registered automatically on import.
/// >>> minkowski = Representation.mink(4)
/// >>> mu, nu = minkowski("mu"), minkowski("nu")
/// >>> p, q, r = TensorName("p"), TensorName("q"), TensorName("r")
/// >>> p_mu = p(mu).to_expression()
/// >>> q_nu, r_nu = q(nu).to_expression(), r(nu).to_expression()
/// >>> factorized = p_mu * (q_nu + r_nu)
/// >>> terms = expand_mink(factorized)
/// >>> sum(structure * coefficient for structure, coefficient in terms) == p_mu * q_nu + p_mu * r_nu
/// True
/// ```
pub fn expand_mink(expression: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(expression.expr.expand_mink())
}

/// Expand products around factors carrying registered bispinor indices.
///
/// # Arguments
/// - `expression`: a factorized Spenso-compatible expression.
///
/// # Returns
/// `(structure, coefficient)` pairs distributed around bispinor-bearing factors. No explicit
/// spinor components are substituted.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_bis
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> # Built-in representations are registered automatically on import.
/// >>> bispinor = Representation.bis(4)
/// >>> alpha, beta = bispinor("alpha"), bispinor("beta")
/// >>> u, v, w = TensorName("u"), TensorName("v"), TensorName("w")
/// >>> u_alpha = u(alpha).to_expression()
/// >>> v_beta, w_beta = v(beta).to_expression(), w(beta).to_expression()
/// >>> factorized = u_alpha * (v_beta + w_beta)
/// >>> terms = expand_bis(factorized)
/// >>> sum(structure * coefficient for structure, coefficient in terms) == u_alpha * v_beta + u_alpha * w_beta
/// True
/// ```
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_bis(expression: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(expression.expr.expand_bis())
}

/// Expand products around factors carrying Minkowski or bispinor indices.
///
/// This combines the selection patterns of `expand_mink()` and `expand_bis()` in one
/// coefficient pass. Other representation families remain in the coefficient sector.
///
/// # Arguments
/// - `expression`: a factorized Spenso-compatible expression.
///
/// # Returns
/// `(structure, coefficient)` pairs distributed around both selected representation families.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_mink_bis
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> # Built-in representations are registered automatically on import.
/// >>> minkowski, bispinor = Representation.mink(4), Representation.bis(4)
/// >>> p_mu = TensorName("p")(minkowski("mu")).to_expression()
/// >>> q_mu = TensorName("q")(minkowski("mu")).to_expression()
/// >>> u_a = TensorName("u")(bispinor("a")).to_expression()
/// >>> v_a = TensorName("v")(bispinor("a")).to_expression()
/// >>> factorized = (p_mu + q_mu) * (u_a + v_a)
/// >>> expected = p_mu * u_a + p_mu * v_a + q_mu * u_a + q_mu * v_a
/// >>> terms = expand_mink_bis(factorized)
/// >>> sum(structure * coefficient for structure, coefficient in terms) == expected
/// True
/// ```
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_mink_bis(expression: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(expression.expr.expand_mink_bis())
}

/// Expand products around registered color factors.
///
/// Fundamental, antifundamental, adjoint, color-chain, color-trace, and supported invariant
/// factors form the selected sector. This only distributes the symbolic expression; use
/// `simplify_color()` separately to apply SU(N) identities.
///
/// # Arguments
/// - `expression`: a factorized Spenso-compatible expression.
///
/// # Returns
/// `(structure, coefficient)` pairs distributed around color-bearing factors.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_color
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> # Built-in representations are registered automatically on import.
/// >>> adjoint, fundamental = Representation.coad(8), Representation.cof(3)
/// >>> antifundamental = fundamental.dual()
/// >>> generator = TensorName.t()
/// >>> t_a = generator(
/// ...     adjoint("a"), fundamental("i"), antifundamental("j")
/// ... ).to_expression()
/// >>> t_b = generator(
/// ...     adjoint("b"), fundamental("k"), antifundamental("l")
/// ... ).to_expression()
/// >>> t_c = generator(
/// ...     adjoint("c"), fundamental("m"), antifundamental("n")
/// ... ).to_expression()
/// >>> factorized = t_a * (t_b + t_c)
/// >>> terms = expand_color(factorized)
/// >>> sum(structure * coefficient for structure, coefficient in terms) == t_a * t_b + t_a * t_c
/// True
/// ```
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_color(expression: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(ColorSimplifier::expand_color(&expression.expr))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Expand products around registered metric tensors.
///
/// This is a structural expansion only. It neither contracts the metrics nor substitutes a
/// dimension or signature; call `simplify_metrics()` separately for supported contractions.
///
/// # Arguments
/// - `expression`: a factorized Spenso-compatible expression.
///
/// # Returns
/// `(structure, coefficient)` pairs distributed around metric factors.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_metrics
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> # Built-in representations are registered automatically on import.
/// >>> minkowski = Representation.mink(4)
/// >>> metric = TensorName.g()
/// >>> g_mn = metric(minkowski("mu"), minkowski("nu")).to_expression()
/// >>> g_rs = metric(minkowski("rho"), minkowski("sigma")).to_expression()
/// >>> g_ab = metric(minkowski("alpha"), minkowski("beta")).to_expression()
/// >>> factorized = g_mn * (g_rs + g_ab)
/// >>> terms = expand_metrics(factorized)
/// >>> sum(structure * coefficient for structure, coefficient in terms) == g_mn * g_rs + g_mn * g_ab
/// True
/// ```
pub fn expand_metrics(expression: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(expression.expr.expand_metrics())
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Wrap all abstract indices with a header symbol
///
/// # Arguments
/// - `expression`: input expression containing tensor indices
/// - `header`: symbol to use as the wrapper function for all indices
///
/// # Returns
/// Expression with all indices wrapped by the header symbol.
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorName, Slot, Representation
/// import symbolica as sp
/// from symbolica.community.idenso import wrap_indices
///
/// T = TensorName("T")
/// rep = Representation.euc(3)
/// # With slots (creates TensorExpression)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(x, mu, nu)  # T(x; mu, nu)
/// print(tensor_with_args)
/// print(wrap_indices(tensor_with_args.to_expression(), sp.S("wrap")))
///
/// ```
pub fn wrap_indices(expression: &PythonExpression, header: Symbol) -> PythonExpression {
    expression.expr.wrap_indices(header).into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
/// Convert complex nested index structures into flattened symbolic names.
///
/// Transforms hierarchical index expressions within tensor function arguments
/// into simplified, flat symbolic representations. This "cooking" process is
/// essential for pattern matching, simplification, and computational efficiency
/// when dealing with complex tensor expressions.
///
/// **Index Cooking Transformation:**
/// - Nested structure: `mink(4, f(g(h(μ))))` → `mink(4, f_g_h_mu)`
/// - Function chains: `lorentz(up(mu))` → `lorentz(up_mu)`
/// - Complex arguments: `tensor(rep(dim,type(idx)))` → `tensor(rep(dim,type_idx))`
///
/// **Scope:**
/// - Only affects indices appearing as function arguments
/// - Preserves top-level function structure
/// # Arguments
/// - `expression`: expression containing complex nested index structures
///
/// # Returns
/// Expression with flattened, simplified index names.
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorName, Slot, Representation
/// import symbolica as sp
/// from symbolica.community.idenso import wrap_indices, cook_indices
///
/// T = TensorName("T")
/// rep = Representation.euc(3)
/// # With slots (creates TensorExpression)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(x, mu, nu)  # T(x; mu, nu)
/// print(tensor_with_args)
/// print(
///     cook_indices(wrap_indices(tensor_with_args.to_expression(), sp.S("wrap")))
/// )
/// ```
pub fn cook_indices(
    expression: &PythonExpression,
    settings: Option<&tooling::PyCookSettings>,
) -> PyResult<PythonExpression> {
    let settings = tooling::PyCookSettings::indices_or(settings);
    settings
        .try_cook_indices(expression.expr.as_view())
        .map(Into::into)
        .map_err(|error| tooling::CookingError::new_err(format!("cannot cook indices: {error:?}")))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
/// Convert a single function call into a flattened variable symbol.
///
/// Transforms a function expression with arguments into a single symbolic variable
/// whose name encodes both the function name and its arguments. This is the
/// atomic version of `cook_indices()`, operating on individual function calls
/// rather than complete expressions.
///
/// **Function Cooking Transform:**
/// - Simple function: `f(a, b)` → `f_a_b`
/// - Nested arguments: `tensor(rep(mu))` → `tensor_rep_mu`
/// - Multiple arguments: `gamma(alpha, beta, mu)` → `gamma_alpha_beta_mu`
/// - Complex names: `my_function(x, y)` → `my_function_x_y`
///
///
/// **Constraints:**
/// - Input must be a single function call (not sum, product, etc.)
/// - Arguments must be cookable (symbols, numbers, simple functions)
/// - Cannot cook expressions containing polynomials or complex structures
///
/// # Arguments
/// - `expression`: expression representing a single function call to cook
///
/// # Returns
/// Expression containing the flattened variable symbol.
///
/// # Raises
/// `TypeError` if input is not a cookable function or contains invalid argument types.
///
/// # Examples:
/// ```python
/// import symbolica as sp
/// from symbolica.community.idenso import cook_function
///
/// # Simple function cooking
/// f = sp.S('f')
/// a, b = sp.S('a','b')
///
/// cooked = cook_function(f(a, b))
/// print(cooked)  # Outputs: f_a_b
/// ```
pub fn cook_function(
    expression: &PythonExpression,
    settings: Option<&tooling::PyCookSettings>,
) -> PyResult<PythonExpression> {
    let settings = tooling::PyCookSettings::flattened_or(settings);
    expression
        .expr
        .cook_function_with_settings(&settings)
        .map_err(|error| tooling::CookingError::new_err(format!("cannot cook: {error:?}")))
        .map(|a| a.into())
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Wraps only the dummy (contracted) indices within the expression using a header symbol.
///
/// Similar to `wrap_indices`, but selectively identifies and wraps only contracted
/// indices (those appearing once upstairs and once downstairs, or twice in a
/// self-dual representation), leaving external (dangling) indices untouched.
/// This is crucial for proper index management in tensor calculations.
///
/// Contracted indices are those that:
/// - Appear in both upper and lower positions (for dualizable reps)
/// - Appear twice in the same position (for self-dual reps)
/// - Are summed over (Einstein summation convention)
///
/// # Arguments
/// - `expression`: input expression containing both dummy and free indices
/// - `header`: symbol to use as wrapper function name for dummy indices only
///
/// # Returns
/// A new expression with only contracted indices wrapped.
///
/// # Raises
/// `ValueError` when the expression cannot be parsed as a tensor network.
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorName, Slot, Representation
/// import symbolica as sp
/// from symbolica.community.idenso import simplify_metrics, wrap_dummies
///
/// T = TensorName("T")
/// rep = Representation.euc(3)
/// # With slots (creates TensorExpression)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(x, mu, nu, nu)  # T(x; mu, nu, nu)
/// # print(tensor_with_args)
/// print(wrap_dummies(tensor_with_args.to_expression(), sp.S("wrap")))
///
/// ```
pub fn wrap_dummies(expression: &PythonExpression, header: Symbol) -> PyResult<PythonExpression> {
    expression
        .expr
        .wrap_dummies::<AbstractIndex>(header)
        .map(Into::into)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Lists the dangling (external, uncontracted) indices present in the expression.
///
/// Identifies and returns all indices that are not summed over (i.e., not dummy
/// indices). These are the "free" indices that appear in the final result and
/// determine the tensor rank of the expression. For dualizable representations,
/// downstairs indices are represented wrapped in `dind(...)`.
///
/// This is essential for:
/// - Verifying index conservation in tensor equations
/// - Determining the rank and structure of tensor expressions
/// - Debugging index contractions
///
/// # Arguments
/// - `expression`: tensor expression to analyze
///
/// # Returns
/// A list of expressions, each representing a free (dangling) index.
///
/// # Raises
/// `ValueError` when the expression cannot be parsed as a tensor network.
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorName, Slot, Representation
/// import symbolica as sp
/// from symbolica.community.idenso import (
///     list_dangling,
/// )
///
/// T = TensorName("T")
/// rep = Representation.euc(3)
/// # With slots (creates TensorExpression)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(x, mu, nu, nu)  # T(x; mu, nu, nu)
/// # print(tensor_with_args)
/// print(list_dangling(tensor_with_args.to_expression()))
/// ```
pub fn list_dangling(expression: &PythonExpression) -> PyResult<Vec<PythonExpression>> {
    expression
        .expr
        .list_dangling::<AbstractIndex>()
        .map(|indices| indices.into_iter().map(Into::into).collect())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
/// Simplify registered Spenso gamma chains and traces with Idenso's default rules.
///
/// The dimension-generic part applies compatible Clifford anticommutation, adjacent
/// contractions, and ordinary odd/even trace recursion. Chisholm identities, gamma-five
/// anticommutation and traces, gamma-zero conjugation, and chiral-projector rules are applied
/// only when the expression carries explicit four-dimensional Minkowski and bispinor
/// representations.
///
/// This function does not select or implement a dimensional-regularization gamma-five scheme.
/// Its gamma-five rules are strictly four-dimensional, and the default Python entry point does
/// not enable the optional three-gamma epsilon expansion available through `GammaSimplifySettings`.
/// Gamma factors must use the Spenso representation-aware forms registered on import;
/// unrecognized plain Symbolica functions are left unchanged.
///
/// # Examples
/// ```python
/// >>> from symbolica import E
/// >>> from symbolica.community.idenso import simplify_gamma
/// >>> # Built-in representations are registered automatically on import.
/// >>> trace = E('''
/// ...     gamma(bis(4,a),bis(4,b),mink(4,mu))
/// ...     * gamma(bis(4,b),bis(4,a),mink(4,nu))
/// ... ''', default_namespace="spenso")
/// >>> simplified = simplify_gamma(trace)
/// >>> "gamma(" not in str(simplified) and "g(" in str(simplified)
/// True
/// ```
///
/// The native gamma argument order is `bis(dim,alpha), bis(dim,beta), mink(dim,mu)`:
/// `alpha` and `beta` are spinor indices, followed by the Lorentz index `mu`.
/// These forms can also be constructed through the HEP tensor library.
///
/// # Arguments
/// - `expression`: expression containing gamma matrix products and traces
///
/// # Returns
/// The simplified expression with gamma algebra applied.
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorLibrary, TensorName
/// from symbolica.community.idenso import simplify_gamma
/// from symbolica import S, Expression
/// # Get HEP library with standard tensors
/// hep_lib = TensorLibrary.hep_lib()
/// # Access standard tensors like gamma matrices
/// gamma_structure = hep_lib[S("spenso::gamma")]
/// print(gamma_structure)
/// print(simplify_gamma(gamma_structure(3, 4, 7) * gamma_structure(7, 4, 3)))
/// ```
pub fn simplify_gamma(
    expression: &PythonExpression,
    settings: Option<&algebra::PyGammaSimplifySettings>,
) -> PythonExpression {
    match settings {
        Some(settings) => expression.expr.simplify_gamma_with(settings.rust()),
        None => expression.expr.simplify_gamma(),
    }
    .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Converts contracted Lorentz/Minkowski indices into dot product notation.
///
/// Automatically identifies and converts patterns like `p(mink(D, mu)) * q(mink(D, mu))`
/// into the compact, representation-carrying `dot(p(mink(D)), q(mink(D)))` notation. This
/// simplification is essential for physics calculations involving four-vectors.
///
/// The function recognizes:
/// - Contracted vector indices: `pᵘqᵤ → p·q`
/// - Multiple contractions: `pᵘqᵤrᵛsᵥ → (p·q)(r·s)`
/// - Self-contractions: `pᵘpᵤ → p²`
///
/// # Arguments
/// - `expression`: expression containing contracted Minkowski vector indices
///
/// # Returns
/// The expression with vector contractions converted to dot products.
///
/// # Examples:
/// ```python
/// from symbolica.community.idenso import to_dots
/// from symbolica.community.spenso import Representation, TensorName
/// p = TensorName("p")
/// q = TensorName("q")
/// rep = Representation.euc(3)
/// # With slots (creates TensorExpression)
/// mu = rep("mu")
/// nu = rep("nu")
///
/// print(to_dots( p(mu)*q(mu)))
/// ```
pub fn to_dots(expression: &PythonExpression) -> PythonExpression {
    expression.expr.to_dots().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Simplifies contractions involving metric tensors and identity tensors.
///
/// Applies fundamental tensor algebra rules for metric and identity tensors:
///
/// **Metric tensor rules:**
/// - `gᵘᵛ pᵥ → pᵘ` (index raising/lowering)
/// - `gᵘᵛ gᵥρ → gᵘρ` or `δᵘρ` (metric composition)
/// - `gᵘᵤ → D` (dimension of spacetime)
/// - `ηᵘᵛ pᵥ → pᵘ` (flat metric contractions)
///
/// **Identity tensor rules:**
/// - `δᵘᵛ pᵥ → pᵘ` (Kronecker delta contraction)
/// - `δᵘᵤ → D` (trace of identity)
///
/// The function recognizes metrics as `spenso::g(...)`
///
/// # Arguments
/// - `expression`: expression containing metric/identity tensor contractions
///
/// # Returns
/// The simplified expression with metric rules applied.
///
/// # Examples:
/// ```python
/// from symbolica.community.idenso import simplify_metrics, to_dots
/// from symbolica.community.spenso import Representation, TensorExpression, TensorName
/// q = TensorName("q")
/// rep = Representation.euc(3)
/// g = TensorExpression.g(rep)
/// # With slots (creates TensorExpression)
/// mu = rep("mu")
/// nu = rep("nu")
/// print(simplify_metrics(g("mu", "nu") * q(mu)))
/// ```
pub fn simplify_metrics(expression: &PythonExpression) -> PythonExpression {
    expression.expr.simplify_metrics().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
/// Simplify registered Spenso color chains, traces, generators, and structure constants.
///
/// With the default Python settings, the simplifier evaluates supported closed traces and
/// expands contractions between generators on separate fundamental chains. Its normalization
/// conventions are
///
/// - `Tr(T^a T^b) = TR δ^{ab}`;
/// - `Σ_a (T^a)_i^j (T^a)_k^l = TR (δ_i^l δ_k^j - δ_i^j δ_k^l/Nc)`;
/// - `Σ_a (T^a)_i^j (T^a)_j^k = CF δ_i^k`;
/// - `Σ_{c,d} f^{acd} f^{bcd} = CA δ^{ab}`.
/// Antisymmetry and Jacobi identities apply to the registered structure constants.
///
/// `CA = Nc`, `CF = (Nc² - 1)/(2Nc)`, and `TR = 1/2` are the conventional fundamental
/// SU(Nc) specialization, not identities imposed on every input. The default simplifier keeps
/// representation invariants symbolic where possible; explicit dimension substitution is a
/// separate `ColorSimplifySettings` option.
///
/// # Examples
/// ```python
/// >>> from symbolica import E
/// >>> from symbolica.community.idenso import simplify_color
/// >>> # Built-in representations are registered automatically on import.
/// >>> generators = E('''
/// ...     t(coad(Nc^2-1,a),cof(Nc,i),dind(cof(Nc,j)))
/// ...     * t(coad(Nc^2-1,a),cof(Nc,k),dind(cof(Nc,l)))
/// ... ''', default_namespace="spenso")
/// >>> simplified = simplify_color(generators)
/// >>> "t(" not in str(simplified) and "g(" in str(simplified)
/// True
/// ```
///
/// **Representation invariants:**
/// Use `Representation.dimension`, `.casimir()`, `.dynkin_index()`, and `.gram(...)`
/// to construct the scalar invariants associated with explicitly typed color structures.
///
/// # Arguments
/// - `expression`: expression containing SU(N) color structures
///
/// # Returns
/// The simplified expression, reduced to representation-owned scalar invariants when possible.
/// Unsupported or open indexed structures may remain explicitly in the result; their presence
/// is not an error.
///
/// # Notes
/// Only representation-aware Spenso color forms are recognized. Plain Symbolica functions with
/// similar names are left unchanged.
///
pub fn simplify_color(
    expression: &PythonExpression,
    settings: Option<&algebra::PyColorSimplifySettings>,
) -> PythonExpression {
    match settings {
        Some(settings) => expression.expr.simplify_color_with(settings.rust()),
        None => expression.expr.simplify_color(),
    }
    .into()
}

pub struct IdensoModule;

macro_rules! define_idenso_python_surface {
    ($($function:ident),+ $(,)?; $($registered:literal),+ $(,)?) => {
        pub(crate) fn initialize_alg_simp(m: &Bound<'_, PyModule>) -> PyResult<()> {
            algebra::register(m)?;
            tooling::register(m)?;
            expansion::register(m)?;
            $(m.add_function(wrap_pyfunction!($function, m)?)?;)+
            let exports = m
                .dict()
                .keys()
                .iter()
                .filter_map(|key| key.extract::<String>().ok())
                .filter(|name| name != "initialize" && !name.starts_with('_'))
                .collect::<Vec<_>>();
            m.add("__all__", exports)?;
            Ok(())
        }

        /// The functions, classes, and exceptions registered on `symbolica.community.idenso`.
        #[cfg(feature = "python_stubgen")]
        pub const PYTHON_STUB_SURFACE: &[&str] = &[$(stringify!($function),)+ $($registered,)+];
    };
}

impl symbolica::api::python::SymbolicaCommunityModule for IdensoModule {
    fn get_name() -> String {
        "idenso".into()
    }

    fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
        initialize_alg_simp(m)
    }

    fn initialize(_py: Python) -> PyResult<()> {
        crate::representations::initialize();
        Ok(())
    }
}

define_idenso_python_surface! {
    simplify_gamma,
    to_dots,
    simplify_metrics,
    simplify_color,
    wrap_indices,
    cook_indices,
    cook_function,
    wrap_dummies,
    list_dangling,
    dirac_adjoint,
    expand_bis,
    expand_mink_bis,
    expand_mink,
    expand_metrics,
    expand_color,
    ;
    "CanonicalizationError",
    "ColorCasimirSettings",
    "ColorSimplifySettings",
    "CookMode",
    "CookSettings",
    "CookSourceFilter",
    "CookTagFilter",
    "CookingError",
    "DiracAdjointError",
    "DotExpansionError",
    "GammaChainOrdering",
    "GammaConjugationError",
    "GammaSimplifySettings",
    "NetworkToolingError",
    "SchoonschipContractionOrder",
    "SchoonschipMode",
    "SchoonschipSettings",
    "SchoonschipTraversal",
    "alias_subtensors",
    "canonize",
    "chainify",
    "collect_chains",
    "collect_color",
    "collect_color_constants",
    "collect_gamma_chains",
    "conjugate_transpose",
    "cook",
    "expand_dots",
    "expand_in_patterns",
    "metric_shorthand_to_dot",
    "normalize_chains",
    "normalize_dots",
    "schoonschip",
    "schoonschip_net",
    "simplify_epsilon",
    "simplify_gamma0",
    "simplify_gamma_conjugate",
    "spenso_conjugate",
    "to_cof_dimension_invariants",
    "to_color_casimir",
    "uncook",
    "undo_all",
    "undo_chain",
    "undo_dots",
    "undo_schoonschip",
    "undo_single_length",
    "undo_trace",
    "wrap_color",
}

#[cfg(test)]
mod tests {
    use pyo3::IntoPyObject;
    use pyo3::types::{PyAnyMethods, PyList};
    use spenso::network::tags::SPENSO_TAG;
    use symbolica::{
        atom::{Atom, FunctionBuilder},
        symbol,
    };

    use super::*;

    const PUBLIC_API: &[&str] = &[
        "CanonicalizationError",
        "ColorCasimirSettings",
        "ColorSimplifySettings",
        "CookMode",
        "CookSettings",
        "CookSourceFilter",
        "CookTagFilter",
        "CookingError",
        "DiracAdjointError",
        "DotExpansionError",
        "GammaChainOrdering",
        "GammaConjugationError",
        "GammaSimplifySettings",
        "NetworkToolingError",
        "SchoonschipContractionOrder",
        "SchoonschipMode",
        "SchoonschipSettings",
        "SchoonschipTraversal",
        "alias_subtensors",
        "canonize",
        "chainify",
        "collect_chains",
        "collect_color",
        "collect_color_constants",
        "collect_gamma_chains",
        "conjugate_transpose",
        "cook",
        "cook_function",
        "cook_indices",
        "dirac_adjoint",
        "expand_bis",
        "expand_color",
        "expand_dots",
        "expand_in_patterns",
        "expand_metrics",
        "expand_mink",
        "expand_mink_bis",
        "list_dangling",
        "metric_shorthand_to_dot",
        "normalize_chains",
        "normalize_dots",
        "schoonschip",
        "schoonschip_net",
        "simplify_color",
        "simplify_epsilon",
        "simplify_gamma",
        "simplify_gamma0",
        "simplify_gamma_conjugate",
        "simplify_metrics",
        "spenso_conjugate",
        "to_cof_dimension_invariants",
        "to_color_casimir",
        "to_dots",
        "uncook",
        "undo_all",
        "undo_chain",
        "undo_dots",
        "undo_schoonschip",
        "undo_single_length",
        "undo_trace",
        "wrap_color",
        "wrap_dummies",
        "wrap_indices",
    ];

    #[test]
    fn registers_exact_public_python_surface() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "idenso").unwrap();
            initialize_alg_simp(&module).unwrap();

            let mut actual = module
                .getattr("__all__")
                .unwrap()
                .cast_into::<PyList>()
                .unwrap()
                .extract::<Vec<String>>()
                .unwrap();
            let mut expected = PUBLIC_API
                .iter()
                .map(|name| (*name).to_string())
                .collect::<Vec<_>>();
            actual.sort();
            expected.sort();
            assert_eq!(actual, expected);
            assert!(module.getattr("initialize").is_err());
            assert!(module.getattr("initialize_module").is_err());
        });
    }

    #[test]
    fn network_tooling_failures_are_python_value_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "idenso")?;
            initialize_alg_simp(&module)?;
            let malformed = PythonExpression {
                expr: FunctionBuilder::new(SPENSO_TAG.dot)
                    .add_arg(Atom::var(symbol!("malformed_dot_operand")))
                    .finish(),
            }
            .into_pyobject(py)?;

            for name in ["undo_dots", "schoonschip_net"] {
                let error = module
                    .getattr(name)?
                    .call1((&malformed,))
                    .expect_err("malformed dot notation should return an error");
                assert!(error.is_instance_of::<NetworkToolingError>(py));
                assert!(error.is_instance_of::<PyValueError>(py));
                assert!(error.to_string().contains("cannot parse tensor network"));
                assert!(error.to_string().contains("Invalid dot function"));
            }

            let error = module
                .getattr("canonize")?
                .call1((&malformed,))
                .expect_err("malformed dot notation should not canonicalize");
            assert!(error.is_instance_of::<CanonicalizationError>(py));
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("cannot parse tensor network"));
            assert!(error.to_string().contains("Invalid dot function"));

            let error = module
                .getattr("dirac_adjoint")?
                .call1((&malformed,))
                .expect_err("malformed dot notation should not have a Dirac adjoint");
            assert!(error.is_instance_of::<DiracAdjointError>(py));
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("cannot parse tensor network"));
            assert!(error.to_string().contains("Invalid dot function"));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn public_python_surface_has_runtime_docstrings() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "idenso").unwrap();
            initialize_alg_simp(&module).unwrap();

            for name in PUBLIC_API {
                let documentation = module
                    .getattr(*name)
                    .unwrap()
                    .getattr("__doc__")
                    .unwrap()
                    .extract::<Option<String>>()
                    .unwrap();
                assert!(
                    documentation.is_some_and(|documentation| !documentation.trim().is_empty()),
                    "{name} is missing its Python docstring"
                );
            }
        });
    }

    #[test]
    fn python_signatures_use_public_names_and_concrete_defaults() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "idenso").unwrap();
            initialize_alg_simp(&module).unwrap();
            let inspect_signature = PyModule::import(py, "inspect")
                .unwrap()
                .getattr("signature")
                .unwrap();

            for (name, expected) in [
                (
                    "GammaSimplifySettings",
                    "(*, chain_ordering=None, evaluate_traces=True, expand_three_gamma_epsilon=False)",
                ),
                (
                    "CookSettings",
                    "(*, mode=None, source=None, output_tags=None, preserve_tags=False)",
                ),
                (
                    "SchoonschipSettings",
                    "(*, depth_limit=1, mode=None, traversal=None, expand_contracted_sums=False, simplify_chain_like_functions=False, schoonschip_rank1_tensors=True, contraction_order=None)",
                ),
            ] {
                let class = module.getattr(name).unwrap();
                let signature = class
                    .getattr("__text_signature__")
                    .unwrap()
                    .extract::<String>()
                    .unwrap();
                assert_eq!(signature, expected, "unexpected signature for {name}");
                let inspected = inspect_signature
                    .call1((&class,))
                    .unwrap()
                    .str()
                    .unwrap()
                    .to_string();
                assert!(
                    !inspected.contains("..."),
                    "inspect.signature retained an ellipsis for {name}: {inspected}"
                );
            }

            for name in [
                "cook_function",
                "cook_indices",
                "dirac_adjoint",
                "expand_bis",
                "expand_color",
                "expand_metrics",
                "expand_mink",
                "expand_mink_bis",
                "list_dangling",
                "simplify_color",
                "simplify_gamma",
                "simplify_metrics",
                "to_dots",
                "wrap_dummies",
                "wrap_indices",
            ] {
                let signature = module
                    .getattr(name)
                    .unwrap()
                    .getattr("__text_signature__")
                    .unwrap()
                    .extract::<String>()
                    .unwrap();
                assert!(
                    signature.contains("expression") && !signature.contains("self_"),
                    "unexpected signature for {name}: {signature}"
                );
            }
        });
    }
}
