use pyo3::prelude::*;
use symbolica::{api::python::PythonExpression, atom::Atom, symbol};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Canonical expression heads for diagram and tensor algebra.
///
/// Use these references for pattern matching and momentum construction without
/// depending on internal namespaces. Model parameters and couplings belong to
/// their model; obtain those through ``model.parameter(name).symbol`` or
/// ``model.coupling(name).symbol`` instead.
///
/// Examples
/// --------
/// >>> from symbolica.community import hep
/// >>> P = hep.Symbols.external_momentum
/// >>> incoming = [P(0), P(1)]
/// >>> mass = hep.Model.standard_model().particle("e-").mass_expression
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Symbols", module = "symbolica.community.feynkit", frozen)]
pub struct PySymbols;

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySymbols {
    /// External momentum family indexed by physical leg.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.external_momentum
    #[classattr]
    fn external_momentum() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::external_momentum()).into()
    }

    /// Loop momentum family indexed by loop basis position.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.loop_momentum
    #[classattr]
    fn loop_momentum() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::loop_momentum()).into()
    }

    /// Momentum family indexed by graph edge.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.edge_momentum
    #[classattr]
    fn edge_momentum() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::momentum()).into()
    }

    /// Tagged propagator denominator head; its fourth argument is the inverse denominator.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.denominator
    #[classattr]
    fn denominator() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::denominator()).into()
    }

    /// Lorentz dimension used by generated diagram expressions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.dimension
    #[classattr]
    fn dimension() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::dimension()).into()
    }

    /// Half-edge index family used to match external tensor slots.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.half_edge
    #[classattr]
    fn half_edge() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::hedge_index()).into()
    }

    /// Vector polarization wavefunction head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.polarization
    #[classattr]
    fn polarization() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::epsilon()).into()
    }

    /// Conjugated vector polarization wavefunction head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.polarization_conjugate
    #[classattr]
    fn polarization_conjugate() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::epsilonbar()).into()
    }

    /// Metric tensor head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.metric
    #[classattr]
    fn metric() -> PythonExpression {
        Atom::var(symbol!("spenso::g")).into()
    }

    /// Dirac gamma tensor head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.gamma
    #[classattr]
    fn gamma() -> PythonExpression {
        Atom::var(symbol!("spenso::gamma")).into()
    }

    /// Antisymmetric color tensor head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.color_f
    #[classattr]
    fn color_f() -> PythonExpression {
        Atom::var(symbol!("spenso::f")).into()
    }

    /// Lorentz representation head, accepting dimension and optionally index.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.lorentz
    #[classattr]
    fn lorentz() -> PythonExpression {
        Atom::var(symbol!("spenso::mink")).into()
    }

    /// Bispinor representation head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.spinor
    #[classattr]
    fn spinor() -> PythonExpression {
        Atom::var(symbol!("spenso::bis")).into()
    }

    /// Fundamental color representation head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.color_fundamental
    #[classattr]
    fn color_fundamental() -> PythonExpression {
        Atom::var(symbol!("spenso::cof")).into()
    }

    /// Adjoint color representation head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.color_adjoint
    #[classattr]
    fn color_adjoint() -> PythonExpression {
        Atom::var(symbol!("spenso::coad")).into()
    }

    /// Number of colors used by symbolic color simplification.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.color_number
    #[classattr]
    fn color_number() -> PythonExpression {
        Atom::var(symbol!("spenso::Nc")).into()
    }

    /// Conjugation wrapper used in tensor expressions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.conjugate
    #[classattr]
    fn conjugate() -> PythonExpression {
        Atom::var(symbol!("spenso::conj")).into()
    }

    /// Scalar product head used by contracted tensor expressions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.dot
    #[classattr]
    fn dot() -> PythonExpression {
        Atom::var(symbol!("spenso::dot")).into()
    }

    /// Color Casimir invariant head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.casimir
    #[classattr]
    fn casimir() -> PythonExpression {
        Atom::var(symbol!("spenso::cas")).into()
    }

    /// Color representation index invariant head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.color_index
    #[classattr]
    fn color_index() -> PythonExpression {
        Atom::var(symbol!("spenso::idx")).into()
    }

    /// Metric head in model propagator definitions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.ufo_metric
    #[classattr]
    fn ufo_metric() -> PythonExpression {
        Atom::var(symbol!("UFO::Metric")).into()
    }

    /// Index placeholder head in model propagator definitions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.ufo_index
    #[classattr]
    fn ufo_index() -> PythonExpression {
        Atom::var(symbol!("UFO::idx")).into()
    }

    /// Momentum head in model propagator definitions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.ufo_momentum
    #[classattr]
    fn ufo_momentum() -> PythonExpression {
        Atom::var(symbol!("UFO::P")).into()
    }
    /// Incoming placeholder for an open gamma chain.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.chain_in
    #[classattr]
    fn chain_in() -> PythonExpression {
        Atom::var(symbol!("spenso::in")).into()
    }
    /// Outgoing placeholder for an open gamma chain.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.chain_out
    #[classattr]
    fn chain_out() -> PythonExpression {
        Atom::var(symbol!("spenso::out")).into()
    }
    /// Compact gamma-chain wrapper.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.gamma_chain
    #[classattr]
    fn gamma_chain() -> PythonExpression {
        Atom::var(symbol!("spenso::chain")).into()
    }
    /// Antisymmetric Levi-Civita tensor head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.levi_civita
    #[classattr]
    fn levi_civita() -> PythonExpression {
        Atom::var(symbol!("spenso::epsilon")).into()
    }
    /// Dirac charge-conjugation matrix head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.charge_conjugation
    #[classattr]
    fn charge_conjugation() -> PythonExpression {
        Atom::var(symbol!("spenso::charge_conjugation")).into()
    }
    /// Cyclic symmetry wrapper.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.cyclic
    #[classattr]
    fn cyclic() -> PythonExpression {
        Atom::var(symbol!("spenso::cyclic")).into()
    }
    /// Symmetric tensor wrapper.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.symmetric
    #[classattr]
    fn symmetric() -> PythonExpression {
        Atom::var(symbol!("spenso::sym")).into()
    }
    /// Antisymmetric tensor wrapper.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.antisymmetric
    #[classattr]
    fn antisymmetric() -> PythonExpression {
        Atom::var(symbol!("spenso::antisym")).into()
    }
    /// Time-component gamma matrix head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.gamma_zero
    #[classattr]
    fn gamma_zero() -> PythonExpression {
        Atom::var(symbol!("spenso::gamma0")).into()
    }
    /// Compact tensor trace wrapper.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> reference = hep.Symbols.trace
    #[classattr]
    fn trace() -> PythonExpression {
        Atom::var(symbol!("spenso::trace")).into()
    }
    /// Complex-conjugation helper used in imported model formulas.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> model = hep.Model.standard_model()
    /// >>> conjugate = hep.Symbols.model_conjugate(model.parameter("CKM1x1").symbol)
    #[classattr]
    fn model_conjugate() -> PythonExpression {
        Atom::var(symbol!("UFO::complexconjugate")).into()
    }
}
