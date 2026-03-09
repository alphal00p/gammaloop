use crate::color::{ColorSimplifier, SelectiveExpand};
use crate::gamma::GammaSimplifier;
use crate::metric::MetricSimplifier;
use crate::representations::initialize;

use pyo3::{
    Bound, PyResult, Python,
    exceptions::PyTypeError,
    pyfunction,
    types::{PyModule, PyModuleMethods},
    wrap_pyfunction,
};
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::atom::Atom;

use crate::IndexTooling;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pyfunction;
use symbolica::atom::Symbol;

use symbolica::api::python::PythonExpression;

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn dirac_adjoint(self_: &PythonExpression) -> PythonExpression {
    self_.expr.dirac_adjoint::<AbstractIndex>().unwrap().into()
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
/// # Args:
///     self_: Expression containing factorized terms with Minkowski indices
///
/// # Returns:
///     Expanded expression with factorizations unfolded
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
/// expanded = expand_mink(factorized)  # p(mu)*q(mu) + p(mu)*r(mu)
///
/// # Complex factorization
/// A = sp.S('A')
/// expr = A * (p(mu) * q(mu) + r(mu))
/// result = expand_mink(expr)  # A*p(mu)*q(mu) + A*r(mu)
/// ```
pub fn expand_mink(self_: &PythonExpression) -> PythonExpression {
    self_
        .expr
        .expand_mink()
        .iter()
        .fold(Atom::Zero, |a, (c, s)| a + c * s)
        .into()
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
/// # Args:
///     self_: Expression containing factorized terms with bispinor indices
///
/// # Returns:
///     Expanded expression with factorizations unfolded
///
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_bis(self_: &PythonExpression) -> PythonExpression {
    self_
        .expr
        .expand_bis()
        .iter()
        .fold(Atom::Zero, |a, (c, s)| a + c * s)
        .into()
}

/// Expands factorized terms containing both Minkowski and bispinor indices.
///
/// Combines the functionality of `expand_mink()` and `expand_bis()` to perform
/// simultaneous expansion of factorized expressions involving both spacetime
/// and spinor indices.
///
/// # Args:
///     self_: Expression containing factorized terms with both index types
///
/// # Returns:
///     Expanded expression with all factorizations unfolded
///
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_mink_bis(self_: &PythonExpression) -> PythonExpression {
    self_
        .expr
        .expand_mink_bis()
        .iter()
        .fold(Atom::Zero, |a, (c, s)| a + c * s)
        .into()
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
/// # Args:
///     self_: Expression containing factorized terms with color indices
///
/// # Returns:
///     Expanded expression with factorizations unfolded
///
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_color(self_: &PythonExpression) -> PythonExpression {
    self_
        .expr
        .expand_color()
        .iter()
        .fold(Atom::Zero, |a, (c, s)| a + c * s)
        .into()
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
/// # Args:
///     self_ (Expression): The expression containing factorized metric terms
///
/// # Returns:
///     Expression: The expanded expression with metric factorizations unfolded
pub fn expand_metrics(self_: &PythonExpression) -> PythonExpression {
    self_
        .expr
        .expand_metrics()
        .iter()
        .fold(Atom::Zero, |a, (c, s)| a + c * s)
        .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Wrap all abstract indices with a header symbol
///
/// # Args:
///     self_: The input expression containing tensor indices
///     header: Symbol to use as the wrapper function for all indices
///
/// # Returns:
///     Expression with all indices wrapped by the header symbol
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorName, Slot, Representation
/// import symbolica as sp
/// from symbolica.community.idenso import wrap_indices
///
/// T = TensorName("T")
/// rep = Representation.euc(3)
/// # With slots (creates TensorIndices)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(mu, nu, x)  # T(mu, nu; x)
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
#[pyfunction]
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
/// # Args:
///     self_: Expression containing complex nested index structures
///
/// # Returns:
///     Expression with flattened, simplified index names
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorName, Slot, Representation
/// import symbolica as sp
/// from symbolica.community.idenso import wrap_indices, cook_indices
///
/// T = TensorName("T")
/// rep = Representation.euc(3)
/// # With slots (creates TensorIndices)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(mu, nu, x)  # T(mu, nu; x)
/// print(tensor_with_args)
/// print(
///     cook_indices(wrap_indices(tensor_with_args.to_expression(), sp.S("wrap")))
/// )
/// ```
pub fn cook_indices(self_: &PythonExpression) -> PythonExpression {
    self_.expr.cook_indices().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
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
/// # Args:
///     self_: Expression representing a single function call to cook
///
/// # Returns:
///     Expression containing the flattened variable symbol
///
/// # Raises:
///     TypeError: If input is not a cookable function or contains invalid argument types
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
pub fn cook_function(self_: &PythonExpression) -> PyResult<PythonExpression> {
    self_
        .expr
        .cook_function()
        .map_err(|a| PyTypeError::new_err(format!("cannot cook: {a:?}")))
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
/// # Args:
///     self_ (Expression): The input expression containing both dummy and free indices
///     header (Symbol): The symbol to use as wrapper function name for dummy indices only
///
/// # Returns:
///     Expression: A new expression with only contracted indices wrapped
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorName, Slot, Representation
/// import symbolica as sp
/// from symbolica.community.idenso import simplify_metrics, wrap_dummies
///
/// T = TensorName("T")
/// rep = Representation.euc(3)
/// # With slots (creates TensorIndices)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(mu, nu, nu, x)  # T(mu, nu; x)
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
/// # Args:
///     self_ (Expression): The tensor expression to analyze
///
/// # Returns:
///     list[Expression]: A list of expressions, each representing a free (dangling) index
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
/// # With slots (creates TensorIndices)
/// mu = rep("mu")
/// nu = rep("nu")
/// x = sp.S("x")
/// tensor_with_args = T(mu, nu, nu, x)  # T(mu, nu; x)
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
#[pyfunction]
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
/// # Args:
///     self_ (Expression): The expression containing gamma matrix products and traces
///
/// # Returns:
///     Expression: The simplified expression with gamma algebra applied
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
pub fn simplify_gamma(self_: &PythonExpression) -> PythonExpression {
    self_.expr.simplify_gamma().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Converts contracted Lorentz/Minkowski indices into dot product notation.
///
/// Automatically identifies and converts patterns like `p(mink(D, mu)) * q(mink(D, mu))`
/// into the more compact and physically meaningful `dot(p, q)` notation. This
/// simplification is essential for physics calculations involving four-vectors.
///
/// The function recognizes:
/// - Contracted vector indices: `pᵘqᵤ → p·q`
/// - Multiple contractions: `pᵘqᵤrᵛsᵥ → (p·q)(r·s)`
/// - Self-contractions: `pᵘpᵤ → p²`
///
/// # Args:
///     self_ (Expression): The expression containing contracted Minkowski vector indices
///
/// # Returns:
///     Expression: The expression with vector contractions converted to dot products
///
/// # Examples:
/// ```python
/// from symbolica.community.idenso import to_dots
/// from symbolica.community.spenso import Representation, TensorName
/// p = TensorName("p")
/// q = TensorName("q")
/// rep = Representation.euc(3)
/// # With slots (creates TensorIndices)
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
/// # Args:
///     self_ (Expression): The expression containing metric/identity tensor contractions
///
/// # Returns:
///     Expression: The simplified expression with metric rules applied
///
/// # Examples:
/// ```python
/// from symbolica.community.idenso import simplify_metrics, to_dots
/// from symbolica.community.spenso import Representation, TensorName
/// q = TensorName("q")
/// g = TensorName.g
/// rep = Representation.euc(3)
/// # With slots (creates TensorIndices)
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
#[pyfunction]
/// Applies SU(N) color algebra rules to simplify color group structures.
///
/// Performs comprehensive simplifications of SU(N) color algebra including:
///
/// **Structure constants:**
/// - `f^{abc}f^{ade} = CA δ^{bc}` (Casimir relations)
/// - `f^{abc}f^{bcd} = CA f^{acd}` (Jacobi identities)
/// - Antisymmetry: `f^{abc} = -f^{bac}`
///
/// **Generators and traces:**
/// - `Tr(T^a T^b) = TR δ^{ab}` (orthogonality)
/// - `T^a_{ij} T^a_{kl} = δ_{il}δ_{jk}/Nc - δ_{ij}δ_{kl}/Nc²` (Fierz identity)
/// - `(T^a)² = CF` (fundamental Casimir)
///
/// **Color factors:**
/// - `Nc`: Number of colors
/// - `CA = Nc`: Adjoint Casimir
/// - `CF = (Nc² - 1)/(2Nc)`: Fundamental Casimir
/// - `TR = 1/2`: Normalization factor
///
/// # Args:
///     self_ (Expression): The expression containing SU(N) color structures
///
/// # Returns:
///     Expression: Simplified expression with color algebra reduced to scalar factors
///                 (Nc, CA, CF, TR) when possible
///
/// # Notes:
///     If explicit color indices remain after simplification, it indicates the
///     expression could not be fully reduced to color-scalar form.
///
pub fn simplify_color(self_: &PythonExpression) -> PythonExpression {
    self_.expr.simplify_color().into()
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
    m.add_function(wrap_pyfunction!(initialize, m)?)?;
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

    Ok(())
}
