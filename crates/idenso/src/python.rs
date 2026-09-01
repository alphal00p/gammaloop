use crate::color::ColorSimplifier;
use crate::dirac::GammaSimplifier;
use crate::selective_expand::SelectiveExpand;
use crate::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};

use crate::{Cookable, IndexTooling};
use pyo3::{
    Bound, PyResult, Python, pyfunction,
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
    CanonicalizationError, CookingError, DiracAdjointError, DotExpansionError, PyCookMode,
    PyCookSettings, PyCookSourceFilter, PyCookTagFilter, PySchoonschipContractionOrder,
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
pub fn dirac_adjoint(self_: &PythonExpression) -> PyResult<PythonExpression> {
    self_
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
/// Expands factorized terms containing Minkowski spacetime indices.
///
/// Finds and expands factorized expressions involving Minkowski tensors and vectors,
/// unfolding multiplicative structures into expanded sums for subsequent simplification.
/// Does not expand into explicit components but rather expands nested factorizations.
///
/// **Factorization Expansion:**
/// - `A(μ) * (B(ν) + C(ν)) → A(μ)*B(ν) + A(μ)*C(ν)`
/// - Nested products with Minkowski indices get distributed
/// - Parenthesized expressions are expanded algebraically
/// - Prepares expressions for metric simplification and contractions
///
/// **Applications:**
/// - Expanding factorized tensor expressions before simplification
/// - Preparing for index contraction algorithms
/// - Unfolding nested products in field theory calculations
/// - Algebraic manipulation of relativistic expressions
///
/// # Arguments
/// - `self_`: expression containing factorized terms with Minkowski indices
///
/// # Returns
/// `(structure, coefficient)` pairs with factorizations unfolded.
///
/// # Examples:
/// ```python
/// import symbolica as sp
/// from symbolica.community.idenso import expand_mink
///
/// # Expand factorized vector expression
/// p = sp.S('p')
/// q = sp.S('q')
/// r = sp.S('r')
/// mu = sp.S('mu')
/// factorized = p(mu) * (q(mu) + r(mu))
/// terms = expand_mink(factorized)
/// for structure, coefficient in terms:
///     print(structure, coefficient)
///
/// # Complex factorization
/// A = sp.S('A')
/// expr = A * (p(mu) * q(mu) + r(mu))
/// terms = expand_mink(expr)
/// ```
pub fn expand_mink(self_: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(self_.expr.expand_mink())
}

/// Expands factorized terms containing Dirac bispinor indices.
///
/// Finds and expands factorized expressions involving bispinor tensors and spinors,
/// unfolding multiplicative structures into expanded sums for subsequent simplification.
/// Does not expand into explicit components but rather expands nested factorizations.
///
/// **Factorization Expansion:**
/// - `ψ(α) * (γ(μ,α,β) + σ(μ,α,β)) → ψ(α)*γ(μ,α,β) + ψ(α)*σ(μ,α,β)`
/// - Nested products with bispinor indices get distributed
/// - Parenthesized expressions are expanded algebraically
/// - Prepares expressions for gamma matrix simplification
///
///
/// # Arguments
/// - `self_`: expression containing factorized terms with bispinor indices
///
/// # Returns
/// `(structure, coefficient)` pairs with factorizations unfolded.
///
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_bis(self_: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(self_.expr.expand_bis())
}

/// Expands factorized terms containing both Minkowski and bispinor indices.
///
/// Combines the functionality of `expand_mink()` and `expand_bis()` to perform
/// simultaneous expansion of factorized expressions involving both spacetime
/// and spinor indices.
///
/// # Arguments
/// - `self_`: expression containing factorized terms with both index types
///
/// # Returns
/// `(structure, coefficient)` pairs with all factorizations unfolded.
///
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_mink_bis(self_: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(self_.expr.expand_mink_bis())
}

/// Expands factorized terms containing SU(N) color indices.
///
/// Finds and expands factorized expressions involving color tensors and fields,
/// unfolding multiplicative structures into expanded sums for subsequent simplification.
/// Does not expand into explicit components but rather expands nested factorizations.
///
/// **Factorization Expansion:**
/// - `q(a) * (T(b,a,c) + S(b,a,c)) → q(a)*T(b,a,c) + q(a)*S(b,a,c)`
/// - Nested products with color indices get distributed
/// - Parenthesized expressions are expanded algebraically
/// - Prepares expressions for color algebra simplification
///
/// **Applications:**
/// - Expanding factorized QCD expressions before simplification
/// - Preparing for SU(N) algebra algorithms
/// - Unfolding nested products in gauge theory calculations
/// - Algebraic manipulation of color structures
///
/// # Arguments
/// - `self_`: expression containing factorized terms with color indices
///
/// # Returns
/// `(structure, coefficient)` pairs with factorizations unfolded.
///
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_color(self_: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(ColorSimplifier::expand_color(&self_.expr))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Expands factorized terms containing metric tensors.
///
/// Finds and expands factorized expressions involving metric tensors and related
/// geometric objects, unfolding multiplicative structures for subsequent simplification.
///
/// # Arguments
/// - `self_`: expression containing factorized metric terms
///
/// # Returns
/// `(structure, coefficient)` pairs with metric factorizations unfolded.
pub fn expand_metrics(self_: &PythonExpression) -> Vec<PythonTerm> {
    expansion::python_terms(self_.expr.expand_metrics())
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Wrap all abstract indices with a header symbol
///
/// # Arguments
/// - `self_`: input expression containing tensor indices
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
pub fn wrap_indices(self_: &PythonExpression, header: Symbol) -> PythonExpression {
    self_.expr.wrap_indices(header).into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (self_, settings = None))]
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
/// - `self_`: expression containing complex nested index structures
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
    self_: &PythonExpression,
    settings: Option<&tooling::PyCookSettings>,
) -> PyResult<PythonExpression> {
    let settings = tooling::PyCookSettings::indices_or(settings);
    settings
        .try_cook_indices(self_.expr.as_view())
        .map(Into::into)
        .map_err(|error| tooling::CookingError::new_err(format!("cannot cook indices: {error:?}")))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (self_, settings = None))]
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
/// - Multiple arguments: `gamma(mu, alpha, beta)` → `gamma_mu_alpha_beta`
/// - Complex names: `my_function(x, y)` → `my_function_x_y`
///
///
/// **Constraints:**
/// - Input must be a single function call (not sum, product, etc.)
/// - Arguments must be cookable (symbols, numbers, simple functions)
/// - Cannot cook expressions containing polynomials or complex structures
///
/// # Arguments
/// - `self_`: expression representing a single function call to cook
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
    self_: &PythonExpression,
    settings: Option<&tooling::PyCookSettings>,
) -> PyResult<PythonExpression> {
    let settings = tooling::PyCookSettings::flattened_or(settings);
    self_
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
/// - `self_`: input expression containing both dummy and free indices
/// - `header`: symbol to use as wrapper function name for dummy indices only
///
/// # Returns
/// A new expression with only contracted indices wrapped.
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
pub fn wrap_dummies(self_: &PythonExpression, header: Symbol) -> PythonExpression {
    self_.expr.wrap_dummies::<AbstractIndex>(header).into()
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
/// - `self_`: tensor expression to analyze
///
/// # Returns
/// A list of expressions, each representing a free (dangling) index.
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
pub fn list_dangling(self_: &PythonExpression) -> Vec<PythonExpression> {
    self_
        .expr
        .list_dangling::<AbstractIndex>()
        .into_iter()
        .map(|a| a.into())
        .collect()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (self_, settings = None))]
/// Applies Clifford algebra rules and trace identities to simplify gamma matrix expressions.
///
/// Performs comprehensive simplifications of Dirac gamma matrix algebra including:
/// - **Anticommutation relations**: `{γᵘ, γᵛ} = 2gᵘᵛ`
/// - **Trace identities**: `Tr(γᵘ) = 0`, `Tr(γᵘγᵛ) = 4gᵘᵛ`, etc.
/// - **Gamma5 properties**: `{γ₅, γᵘ} = 0`, `(γ₅)² = 1`
/// - **Chain simplifications**: Reduces products of gamma matrices
/// - **Contraction rules**: Simplifies contracted gamma matrix products
///
/// The function recognizes gamma matrices represented as `spenso::gamma(spenso::mink(dim,mu), spenso::bis(dim,alpha), spenso::bis(dim,beta))`
/// where `mu` is the Lorentz index and `alpha`, `beta` are spinor indices.
/// These can be easily created using the hep_lib.
///
/// # Arguments
/// - `self_`: expression containing gamma matrix products and traces
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
/// print(simplify_gamma(gamma_structure(7, 3, 4) * gamma_structure(3, 7, 4)))
/// ```
pub fn simplify_gamma(
    self_: &PythonExpression,
    settings: Option<&algebra::PyGammaSimplifySettings>,
) -> PythonExpression {
    match settings {
        Some(settings) => self_.expr.simplify_gamma_with(settings.rust()),
        None => self_.expr.simplify_gamma(),
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
/// - `self_`: expression containing contracted Minkowski vector indices
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
pub fn to_dots(self_: &PythonExpression) -> PythonExpression {
    self_.expr.to_dots().into()
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
/// - `self_`: expression containing metric/identity tensor contractions
///
/// # Returns
/// The simplified expression with metric rules applied.
///
/// # Examples:
/// ```python
/// from symbolica.community.idenso import simplify_metrics, to_dots
/// from symbolica.community.spenso import Representation, TensorName
/// q = TensorName("q")
/// g = TensorName.g()
/// rep = Representation.euc(3)
/// # With slots (creates TensorExpression)
/// mu = rep("mu")
/// nu = rep("nu")
/// print(simplify_metrics(g(mu, nu) * q(mu)))
/// ```
pub fn simplify_metrics(self_: &PythonExpression) -> PythonExpression {
    self_.expr.simplify_metrics().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (self_, settings = None))]
/// Applies SU(N) color algebra rules to simplify color group structures.
///
/// Performs comprehensive simplifications of SU(N) color algebra including:
///
/// **Structure constants:**
/// - Contractions expressed through the adjoint representation's Casimir
/// - Antisymmetry and Jacobi identities
/// - Antisymmetry: `f^{abc} = -f^{bac}`
///
/// **Generators and traces:**
/// - Generator traces expressed through the representation's Dynkin index
/// - Fierz identities using the supplied fundamental and adjoint dimensions
/// - Generator squares expressed through the representation's Casimir
///
/// **Representation invariants:**
/// Use `Representation.dimension`, `.casimir()`, `.dynkin_index()`, and `.gram(...)`
/// to construct the scalar invariants associated with explicitly typed color structures.
///
/// # Arguments
/// - `self_`: expression containing SU(N) color structures
///
/// # Returns
/// Simplified expression with color algebra reduced to representation-owned scalar
/// invariants when possible.
///
/// # Notes
/// If explicit color indices remain after simplification, it indicates the expression
/// could not be fully reduced to color-scalar form.
///
pub fn simplify_color(
    self_: &PythonExpression,
    settings: Option<&algebra::PyColorSimplifySettings>,
) -> PythonExpression {
    match settings {
        Some(settings) => self_.expr.simplify_color_with(settings.rust()),
        None => self_.expr.simplify_color(),
    }
    .into()
}

pub struct IdensoModule;

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

pub(crate) fn initialize_alg_simp(m: &Bound<'_, PyModule>) -> PyResult<()> {
    algebra::register(m)?;
    tooling::register(m)?;
    expansion::register(m)?;
    m.add_function(wrap_pyfunction!(simplify_gamma, m)?)?;
    m.add_function(wrap_pyfunction!(to_dots, m)?)?;
    m.add_function(wrap_pyfunction!(simplify_metrics, m)?)?;
    m.add_function(wrap_pyfunction!(simplify_color, m)?)?;
    m.add_function(wrap_pyfunction!(wrap_indices, m)?)?;
    m.add_function(wrap_pyfunction!(cook_indices, m)?)?;
    m.add_function(wrap_pyfunction!(cook_function, m)?)?;
    m.add_function(wrap_pyfunction!(wrap_dummies, m)?)?;
    m.add_function(wrap_pyfunction!(list_dangling, m)?)?;
    m.add_function(wrap_pyfunction!(dirac_adjoint, m)?)?;
    m.add_function(wrap_pyfunction!(expand_bis, m)?)?;
    m.add_function(wrap_pyfunction!(expand_mink_bis, m)?)?;
    m.add_function(wrap_pyfunction!(expand_mink, m)?)?;
    m.add_function(wrap_pyfunction!(expand_metrics, m)?)?;
    m.add_function(wrap_pyfunction!(expand_color, m)?)?;

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

#[cfg(test)]
mod tests {
    use pyo3::types::{PyAnyMethods, PyList};

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
}
