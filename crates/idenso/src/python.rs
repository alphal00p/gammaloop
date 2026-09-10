use crate::color::ColorSimplifier;
use crate::dirac::GammaSimplifier;
use crate::representations::initialize;
use crate::selective_expand::SelectiveExpand;
use crate::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};

use pyo3::{
    Bound, PyResult, Python,
    exceptions::PyTypeError,
    pyfunction,
    types::{PyModule, PyModuleMethods},
    wrap_pyfunction,
};
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::atom::Atom;

use crate::{Cookable, IndexTooling};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pyfunction;
use symbolica::atom::Symbol;

use symbolica::api::python::PythonExpression;

/// Return the Dirac adjoint of a Symbolica tensor expression.
///
/// Idenso takes the symbolic complex conjugate, reverses compatible open bispinor chains, and
/// inserts the registered `gamma0` factors required at dangling bispinor slots.
/// The input must use the representation-aware Spenso forms registered by `initialize()`.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import dirac_adjoint, initialize, list_dangling
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> initialize()
/// >>> initialize() is None  # Registration is idempotent.
/// True
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
/// - `self_`: a Spenso-compatible tensor expression.
///
/// # Returns
/// The representation-aware Dirac adjoint.
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
/// Expand products around factors carrying registered Minkowski indices.
///
/// This is a selective symbolic expansion: Minkowski-bearing factors become polynomial
/// variables while unrelated sectors remain coefficients. It does not substitute explicit
/// four-vector components or choose a metric signature.
///
/// # Arguments
/// - `self_`: a factorized Spenso-compatible expression.
///
/// # Returns
/// The expression distributed around its Minkowski-bearing factors.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_mink, initialize
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> initialize()
/// >>> minkowski = Representation.mink(4)
/// >>> mu, nu = minkowski("mu"), minkowski("nu")
/// >>> p, q, r = TensorName("p"), TensorName("q"), TensorName("r")
/// >>> p_mu = p(mu).to_expression()
/// >>> q_nu, r_nu = q(nu).to_expression(), r(nu).to_expression()
/// >>> factorized = p_mu * (q_nu + r_nu)
/// >>> expand_mink(factorized) == p_mu * q_nu + p_mu * r_nu
/// True
/// ```
pub fn expand_mink(self_: &PythonExpression) -> PythonExpression {
    self_
        .expr
        .expand_mink()
        .iter()
        .fold(Atom::Zero, |a, (c, s)| a + c * s)
        .into()
}

/// Expand products around factors carrying registered bispinor indices.
///
/// # Arguments
/// - `self_`: a factorized Spenso-compatible expression.
///
/// # Returns
/// The expression distributed around its bispinor-bearing factors. No explicit spinor
/// components are substituted.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_bis, initialize
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> initialize()
/// >>> bispinor = Representation.bis(4)
/// >>> alpha, beta = bispinor("alpha"), bispinor("beta")
/// >>> u, v, w = TensorName("u"), TensorName("v"), TensorName("w")
/// >>> u_alpha = u(alpha).to_expression()
/// >>> v_beta, w_beta = v(beta).to_expression(), w(beta).to_expression()
/// >>> factorized = u_alpha * (v_beta + w_beta)
/// >>> expand_bis(factorized) == u_alpha * v_beta + u_alpha * w_beta
/// True
/// ```
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

/// Expand products around factors carrying Minkowski or bispinor indices.
///
/// This combines the selection patterns of `expand_mink()` and `expand_bis()` in one
/// coefficient pass. Other representation families remain in the coefficient sector.
///
/// # Arguments
/// - `self_`: a factorized Spenso-compatible expression.
///
/// # Returns
/// The expression distributed around both selected representation families.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_mink_bis, initialize
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> initialize()
/// >>> minkowski, bispinor = Representation.mink(4), Representation.bis(4)
/// >>> p_mu = TensorName("p")(minkowski("mu")).to_expression()
/// >>> q_mu = TensorName("q")(minkowski("mu")).to_expression()
/// >>> u_a = TensorName("u")(bispinor("a")).to_expression()
/// >>> v_a = TensorName("v")(bispinor("a")).to_expression()
/// >>> factorized = (p_mu + q_mu) * (u_a + v_a)
/// >>> expected = p_mu * u_a + p_mu * v_a + q_mu * u_a + q_mu * v_a
/// >>> expand_mink_bis(factorized) == expected
/// True
/// ```
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

/// Expand products around registered color factors.
///
/// Fundamental, antifundamental, adjoint, color-chain, color-trace, and supported invariant
/// factors form the selected sector. This only distributes the symbolic expression; use
/// `simplify_color()` separately to apply SU(N) identities.
///
/// # Arguments
/// - `self_`: a factorized Spenso-compatible expression.
///
/// # Returns
/// The expression distributed around its color-bearing factors.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_color, initialize
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> initialize()
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
/// >>> expand_color(factorized) == t_a * t_b + t_a * t_c
/// True
/// ```
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
/// Expand products around registered metric tensors.
///
/// This is a structural expansion only. It neither contracts the metrics nor substitutes a
/// dimension or signature; call `simplify_metrics()` separately for supported contractions.
///
/// # Arguments
/// - `self_`: a factorized Spenso-compatible expression.
///
/// # Returns
/// The expression distributed around its metric factors.
///
/// # Examples
/// ```python
/// >>> from symbolica.community.idenso import expand_metrics, initialize
/// >>> from symbolica.community.spenso import Representation, TensorName
/// >>> initialize()
/// >>> minkowski = Representation.mink(4)
/// >>> metric = TensorName.g()
/// >>> g_mn = metric(minkowski("mu"), minkowski("nu")).to_expression()
/// >>> g_rs = metric(minkowski("rho"), minkowski("sigma")).to_expression()
/// >>> g_ab = metric(minkowski("alpha"), minkowski("beta")).to_expression()
/// >>> factorized = g_mn * (g_rs + g_ab)
/// >>> expand_metrics(factorized) == g_mn * g_rs + g_mn * g_ab
/// True
/// ```
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
/// not enable the optional three-gamma epsilon expansion available through Rust settings.
/// Gamma factors must use the Spenso representation-aware forms registered by `initialize()`;
/// unrecognized plain Symbolica functions are left unchanged.
///
/// # Examples
/// ```python
/// >>> from symbolica import E
/// >>> from symbolica.community.idenso import initialize, simplify_gamma
/// >>> initialize()
/// >>> trace = E('''
/// ...     gamma(bis(4,a),bis(4,b),mink(4,mu))
/// ...     * gamma(bis(4,b),bis(4,a),mink(4,nu))
/// ... ''', default_namespace="spenso")
/// >>> simplified = simplify_gamma(trace)
/// >>> "gamma(" not in str(simplified) and "g(" in str(simplified)
/// True
/// ```
///
/// # Arguments
/// - `self_`: expression containing gamma matrix products and traces
///
/// # Returns
/// The simplified expression with gamma algebra applied.
///
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
///
/// `CA = Nc`, `CF = (Nc² - 1)/(2Nc)`, and `TR = 1/2` are the conventional fundamental
/// SU(Nc) specialization, not identities imposed on every input. The default simplifier keeps
/// representation invariants symbolic where possible; explicit dimension substitution is a
/// separate Rust setting.
///
/// # Examples
/// ```python
/// >>> from symbolica import E
/// >>> from symbolica.community.idenso import initialize, simplify_color
/// >>> initialize()
/// >>> generators = E('''
/// ...     t(coad(Nc^2-1,a),cof(Nc,i),dind(cof(Nc,j)))
/// ...     * t(coad(Nc^2-1,a),cof(Nc,k),dind(cof(Nc,l)))
/// ... ''', default_namespace="spenso")
/// >>> simplified = simplify_color(generators)
/// >>> "t(" not in str(simplified) and "g(" in str(simplified)
/// True
/// ```
///
/// # Arguments
/// - `self_`: expression containing SU(N) color structures
///
/// # Returns
/// The simplified expression. Unsupported or open indexed structures may remain explicitly in
/// the result; their presence is not an error.
///
/// # Notes
/// Only representation-aware Spenso color forms are recognized. Plain Symbolica functions with
/// similar names are left unchanged.
///
pub fn simplify_color(self_: &PythonExpression) -> PythonExpression {
    self_.expr.simplify_color().into()
}

pub struct IdensoModule;

macro_rules! define_idenso_python_surface {
    ($($function:ident),+ $(,)?) => {
        pub(crate) fn initialize_alg_simp(m: &Bound<'_, PyModule>) -> PyResult<()> {
            $(m.add_function(wrap_pyfunction!($function, m)?)?;)+
            Ok(())
        }

        /// The functions registered on `symbolica.community.idenso`.
        #[cfg(feature = "python_stubgen")]
        pub const PYTHON_STUB_SURFACE: &[&str] = &[$(stringify!($function),)+];
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
    initialize,
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
}
