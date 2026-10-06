use feynkit_kinematics::Wavefunction;
use pyo3::{prelude::*, types::PyComplex};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// A fixed numerical external state in GammaLoop's MadGraph phase convention.
///
/// Obtain states from ``FourMomentum.wavefunction(kind, helicity)``. Vector
/// components use ``(E,x,y,z)`` and signature ``+---``; spinors use the chiral
/// gamma-matrix basis. These external states have four components independently
/// of the dimension used for internal symbolic Lorentz/Dirac algebra. A scalar
/// has one component. No helicity sum, spin average or coupling is included.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> momentum = hep.FourMomentum(150.0, 0.0, 0.0, 150.0)
/// >>> state = momentum.wavefunction("epsilon", hep.Helicity.PLUS)
/// >>> assert state.kind == "epsilon" and len(state) == 4
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Wavefunction",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyWavefunction {
    pub(crate) inner: Wavefunction<f64>,
}

impl PyWavefunction {
    /// Borrow the shared native state for other HEPKit backends.
    pub fn as_wavefunction(&self) -> &Wavefunction<f64> {
        &self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyWavefunction {
    /// One of ``scalar``, ``epsilon``, ``epsilon_bar``, ``u``, ``u_bar``, ``v``, ``v_bar``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wavefunction`` class example:
    ///
    /// >>> assert state.kind == "epsilon"
    #[getter]
    fn kind(&self) -> String {
        self.inner.kind().to_string()
    }

    /// Return a copy of the numerical components as native Python complex values.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wavefunction`` class example:
    ///
    /// >>> values = state.components
    /// >>> assert len(values) == 4 and values[0] == 0j
    /// >>> assert abs(values[1] + 2**-0.5) < 1e-14
    #[getter]
    fn components<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyComplex>> {
        self.inner
            .components()
            .iter()
            .map(|c| PyComplex::from_doubles(py, c.re, c.im))
            .collect()
    }

    /// Conjugate a scalar/vector or take the chiral Dirac adjoint of a spinor.
    ///
    /// Calling this twice restores both the original components and state kind.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wavefunction`` class example:
    ///
    /// >>> assert state.bar().kind == "epsilon_bar"
    /// >>> assert state.bar().bar() == state
    fn bar(&self) -> Self {
        Self {
            inner: self.inner.bar(),
        }
    }

    /// Return one for scalar states and four for vector or spinor states.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wavefunction`` class example:
    ///
    /// >>> assert len(state) == 4
    fn __len__(&self) -> usize {
        self.inner.components().len()
    }
    /// Describe the numerical state kind and its components.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wavefunction`` class example:
    ///
    /// >>> assert "Wavefunction" in repr(state)
    fn __repr__(&self) -> String {
        format!(
            "Wavefunction(kind={:?}, components={:?})",
            self.inner.kind().to_string(),
            self.inner.components()
        )
    }
    /// Compare both the external-state kind and its numerical components.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wavefunction`` class example:
    ///
    /// >>> assert state == state.bar().bar()
    ///
    /// Parameters
    /// ----------
    /// other : Wavefunction
    ///     State to compare with this one.
    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}
