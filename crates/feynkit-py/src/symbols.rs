use pyo3::prelude::*;
use symbolica::{api::python::PythonExpression, atom::Atom, symbol};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Canonical expression heads owned by diagrams and imported models.
///
/// Tensor and representation vocabulary belongs to `symbolica.community.tensor`.
/// Use these references for diagram patterns and momentum construction without
/// depending on internal namespaces. External and loop momenta belong to
/// ``Kinematics``. Model parameters and couplings belong to
/// their model; obtain those through ``model.parameter(name).symbol`` or
/// ``model.coupling(name).symbol`` instead.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> Q = hep.Symbols.edge_momentum
/// >>> edge_momentum = Q(0)
/// >>> mass = hep.Model.standard_model().particle("e-").mass
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Symbols", module = "symbolica.community.hepkit", frozen)]
pub struct PySymbols;

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PySymbols {
    /// Momentum family indexed by graph edge.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.edge_momentum
    #[classattr]
    fn edge_momentum() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::momentum()).into()
    }

    /// Tagged propagator denominator head; its fourth argument is the inverse denominator.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.denominator
    #[classattr]
    fn denominator() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::denominator()).into()
    }

    /// Lorentz dimension used by generated diagram expressions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.dimension
    #[classattr]
    fn dimension() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::dimension()).into()
    }

    /// Half-edge index family used to match external tensor slots.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.half_edge
    #[classattr]
    fn half_edge() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::hedge_index()).into()
    }

    /// Vector polarization wavefunction head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.polarization
    #[classattr]
    fn polarization() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::epsilon()).into()
    }

    /// Conjugated vector polarization wavefunction head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.polarization_conjugate
    #[classattr]
    fn polarization_conjugate() -> PythonExpression {
        Atom::var(feynkit_graph::symbols::epsilonbar()).into()
    }

    /// Metric head in model propagator definitions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.ufo_metric
    #[classattr]
    fn ufo_metric() -> PythonExpression {
        Atom::var(symbol!("UFO::Metric")).into()
    }

    /// Index placeholder head in model propagator definitions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.ufo_index
    #[classattr]
    fn ufo_index() -> PythonExpression {
        Atom::var(symbol!("UFO::idx")).into()
    }

    /// Momentum head in model propagator definitions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> reference = hep.Symbols.ufo_momentum
    #[classattr]
    fn ufo_momentum() -> PythonExpression {
        Atom::var(symbol!("UFO::P")).into()
    }
    /// Complex-conjugation helper used in imported model formulas.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hepkit as hep
    /// >>> model = hep.Model.standard_model()
    /// >>> conjugate = hep.Symbols.model_conjugate(model.parameter("CKM1x1").symbol)
    #[classattr]
    fn model_conjugate() -> PythonExpression {
        Atom::var(symbol!("UFO::complexconjugate")).into()
    }
}
