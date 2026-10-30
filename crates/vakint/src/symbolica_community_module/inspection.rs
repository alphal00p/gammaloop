//! Read-only Python views of the shared numerical and inspection records.

use std::collections::BTreeMap;

use feynkit_py::PyFeynmanDiagram;
use pyo3::{exceptions::PyTypeError, prelude::*};
use symbolica::{
    api::python::PythonExpression,
    atom::Atom,
    domains::float::{PythonComplexFloat, PythonFloat},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;

use crate::{
    AvailablePrecision, CoefficientPrecision, ComponentPrecision, ConstantSymbol, DataProvenance,
    ErrorQuality, MasterCoefficient, MasterSymbol, MasterSymbolicExpression, NumericalCoefficient,
    PrecisionSource, PrecisionWarning, VakintConstant, VakintMaster,
    presentation::{Content, Present, Presentation},
};

use super::{presentation::render, vakint_to_python_error};

macro_rules! record {
    // Keep the Python name as a token tree so both PyO3 and its stub generator
    // see the literal, rather than an opaque macro-fragment group.
    ($wrapper:ident, $python:tt, $native:ty, {$($methods:tt)*}) => {
        /// Read-only record from Vakint's precision-aware data snapshot.
        #[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
        #[pyclass(name = $python, module = "symbolica.community.hepkit.vakint", frozen, from_py_object)]
        #[derive(Clone)]
        pub struct $wrapper {
            pub(super) value: $native,
        }

        impl From<$native> for $wrapper {
            fn from(value: $native) -> Self { Self { value } }
        }

        impl super::ModuleInit for $wrapper {}

        #[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
        #[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
        #[pymethods]
        impl $wrapper {
            /// Structured terminal representation; respects terminal color settings.
            fn __repr__(&self) -> String { self.value.presentation().to_string() }

            fn __str__(&self) -> String { self.__repr__() }

            /// Render the stored record without numerical reevaluation.
            fn to_html(&self, py: Python<'_>) -> PyResult<String> {
                render(py, &self.value.presentation(), None)
            }

            fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> { self.to_html(py) }

            $($methods)*
        }
    };
}

impl Present for ErrorQuality {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Numerical error quality");
        card.field("Quality", format!("{self:?}"));
        card
    }
}

record!(ErrorQualityWrapper, "ErrorQuality", ErrorQuality, {
    fn __eq__(&self, other: &Self) -> bool {
        self.name() == other.name()
    }
    fn __hash__(&self) -> u64 {
        u64::from(self.name() == "Estimated")
    }
    #[getter]
    fn name(&self) -> &'static str {
        match self.value {
            ErrorQuality::Bounded => "Bounded",
            ErrorQuality::Estimated => "Estimated",
        }
    }
    #[classattr]
    #[pyo3(name = "Bounded")]
    fn bounded() -> ErrorQualityWrapper {
        ErrorQuality::Bounded.into()
    }
    #[classattr]
    #[pyo3(name = "Estimated")]
    fn estimated() -> ErrorQualityWrapper {
        ErrorQuality::Estimated.into()
    }
});

record!(
    AvailablePrecisionWrapper,
    "AvailablePrecision",
    AvailablePrecision,
    {
        fn __eq__(&self, other: &Self) -> bool {
            self.kind() == other.kind()
                && self.max_decimal_precision() == other.max_decimal_precision()
        }
        fn __hash__(&self) -> u64 {
            self.max_decimal_precision()
                .map_or(0, |digits| u64::from(digits) + 1)
        }
        #[getter]
        fn kind(&self) -> &'static str {
            match self.value {
                AvailablePrecision::ArbitraryPrecision => "ArbitraryPrecision",
                AvailablePrecision::FixedPrecision(_) => "FixedPrecision",
            }
        }
        #[getter]
        fn max_decimal_precision(&self) -> Option<u32> {
            match self.value {
                AvailablePrecision::ArbitraryPrecision => None,
                AvailablePrecision::FixedPrecision(digits) => Some(digits),
            }
        }
        #[classattr]
        #[pyo3(name = "ArbitraryPrecision")]
        fn arbitrary_precision() -> AvailablePrecisionWrapper {
            AvailablePrecision::ArbitraryPrecision.into()
        }
        #[staticmethod]
        #[pyo3(name = "FixedPrecision")]
        fn fixed_precision(max_decimal_precision: u32) -> AvailablePrecisionWrapper {
            AvailablePrecision::FixedPrecision(max_decimal_precision).into()
        }
    }
);

record!(
    ComponentPrecisionWrapper,
    "ComponentPrecision",
    ComponentPrecision,
    {
        /// Estimated significant decimal digits, or None when no relative accuracy is available.
        #[getter]
        fn significant_digits(&self) -> Option<u32> {
            self.value.significant_digits()
        }
        #[getter]
        fn significant_bits(&self) -> Option<u32> {
            self.value.significant_bits()
        }
        /// Absolute error radius as a native arbitrary-precision Symbolica Float.
        #[getter]
        fn absolute_error(&self) -> PythonFloat {
            PythonFloat(self.value.absolute_error().clone())
        }
        #[getter]
        fn exact_zero(&self) -> bool {
            self.value.exact_zero()
        }
        #[getter]
        fn quality(&self) -> ErrorQualityWrapper {
            self.value.quality().into()
        }
    }
);

record!(
    CoefficientPrecisionWrapper,
    "CoefficientPrecision",
    CoefficientPrecision,
    {
        #[getter]
        fn real(&self) -> ComponentPrecisionWrapper {
            self.value.real().clone().into()
        }
        #[getter]
        fn imaginary(&self) -> ComponentPrecisionWrapper {
            self.value.imaginary().clone().into()
        }
    }
);

record!(
    PrecisionSourceWrapper,
    "PrecisionSource",
    PrecisionSource,
    {
        #[getter]
        fn identity(&self) -> &str {
            &self.value.identity
        }
        #[getter]
        fn master_id(&self) -> Option<&str> {
            self.value.master_id.as_deref()
        }
        #[getter]
        fn epsilon_power(&self) -> Option<i64> {
            self.value.epsilon_power
        }
        #[getter]
        fn available_decimal_precision(&self) -> Option<u32> {
            self.value.available_decimal_precision
        }
        #[getter]
        fn source_sha256(&self) -> &str {
            &self.value.source_sha256
        }
    }
);

record!(
    PrecisionWarningWrapper,
    "PrecisionWarning",
    PrecisionWarning,
    {
        #[getter]
        fn epsilon_power(&self) -> Option<i64> {
            self.value.epsilon_power
        }
        #[getter]
        fn requested_decimal_precision(&self) -> u32 {
            self.value.requested_decimal_precision
        }
        #[getter]
        fn achieved_decimal_precision(&self) -> Option<u32> {
            self.value.achieved_decimal_precision
        }
        #[getter]
        fn sources(&self) -> Vec<PrecisionSourceWrapper> {
            self.value.sources.iter().cloned().map(Into::into).collect()
        }
        #[getter]
        fn data_version(&self) -> &str {
            &self.value.data_version
        }
        #[getter]
        fn message(&self) -> &str {
            &self.value.message
        }
    }
);

record!(DataProvenanceWrapper, "DataProvenance", DataProvenance, {
    #[getter]
    fn version(&self) -> &str {
        self.value.version()
    }
    #[getter]
    fn source_sha256(&self) -> &str {
        self.value.source_sha256()
    }
    #[getter]
    fn master_id(&self) -> Option<&str> {
        self.value.master_id()
    }
    #[getter]
    fn epsilon_power(&self) -> Option<i64> {
        self.value.epsilon_power()
    }
});

record!(
    NumericalCoefficientWrapper,
    "NumericalCoefficient",
    NumericalCoefficient,
    {
        /// Raw central value. Inspect precision before interpreting its digits.
        #[getter]
        fn numerical_value(&self) -> PythonComplexFloat {
            PythonComplexFloat(self.value.value().clone())
        }
        #[getter]
        fn precision(&self) -> CoefficientPrecisionWrapper {
            self.value.precision().clone().into()
        }
        /// Encode independently achieved component precision; rejects uncertain zeros.
        fn to_expression(&self) -> PyResult<PythonExpression> {
            self.value
                .to_atom()
                .map(Into::into)
                .map_err(vakint_to_python_error)
        }
    }
);

fn symbolic_identity(value: &Bound<'_, PyAny>) -> PyResult<Atom> {
    if let Ok(expression) = value.extract::<PythonExpression>() {
        return Ok(expression.expr);
    }
    if let Ok(name) = value.extract::<String>() {
        return crate::vakint_parse!(&name).map_err(PyTypeError::new_err);
    }
    Err(PyTypeError::new_err(
        "symbol must be a Symbolica expression or a name string",
    ))
}

macro_rules! identifier {
    ($wrapper:ident, $python:tt, $native:ty) => {
        impl Present for $native {
            fn presentation(&self) -> Presentation {
                let mut card = Presentation::new($python);
                card.field("Identity", Content::Expression(self.as_atom().clone()));
                card
            }
        }
        record!($wrapper, $python, $native, {
            fn __eq__(&self, other: &Self) -> bool { self.value == other.value }
            fn __hash__(&self) -> u64 {
                use std::hash::{Hash, Hasher};
                let mut hasher = std::hash::DefaultHasher::new();
                self.value.hash(&mut hasher);
                hasher.finish()
            }
            #[new]
            fn new(
                #[gen_stub(override_type(type_repr="symbolica.Expression | str", imports=("symbolica")))]
                symbol: &Bound<'_, PyAny>,
            ) -> PyResult<Self> {
                <$native>::try_from(symbolic_identity(symbol)?).map(Into::into).map_err(vakint_to_python_error)
            }
            #[getter]
            fn expression(&self) -> PythonExpression { self.value.as_atom().clone().into() }
        });
    };
}

identifier!(ConstantSymbolWrapper, "ConstantSymbol", ConstantSymbol);
identifier!(MasterSymbolWrapper, "MasterSymbol", MasterSymbol);

record!(
    MasterSymbolicExpressionWrapper,
    "MasterSymbolicExpression",
    MasterSymbolicExpression,
    {
        #[getter]
        fn kind(&self) -> &'static str {
            match &self.value {
                MasterSymbolicExpression::Expression(_) => "Expression",
                MasterSymbolicExpression::NumericalOnly(_) => "NumericalOnly",
            }
        }
        #[getter]
        fn expression(&self) -> Option<PythonExpression> {
            match &self.value {
                MasterSymbolicExpression::Expression(value) => Some(value.clone().into()),
                MasterSymbolicExpression::NumericalOnly(_) => None,
            }
        }
        #[getter]
        fn max_decimal_precision(&self) -> Option<u32> {
            match &self.value {
                MasterSymbolicExpression::Expression(_) => None,
                MasterSymbolicExpression::NumericalOnly(digits) => Some(*digits),
            }
        }
    }
);

record!(VakintConstantWrapper, "VakintConstant", VakintConstant, {
    /// Encode achieved precision in a native expression; rejects uncertain zeros.
    fn to_expression(&self) -> PyResult<PythonExpression> {
        self.value
            .to_atom()
            .map(Into::into)
            .map_err(vakint_to_python_error)
    }
    #[getter]
    fn requested_decimal_precision(&self) -> u32 {
        self.value.requested_decimal_precision()
    }
    #[getter]
    fn symbol(&self) -> ConstantSymbolWrapper {
        self.value.symbol().clone().into()
    }
    /// The defining (multiple) integral in Vakint's inert `Integral`,
    /// `LoopIntegral` and `SeriesCoefficient` heads; never evaluated.
    #[getter]
    fn definition(&self) -> PythonExpression {
        self.value.definition().clone().into()
    }
    #[getter]
    fn available_precision(&self) -> AvailablePrecisionWrapper {
        self.value.available_precision().into()
    }
    #[getter]
    fn numerical_value(&self) -> PythonFloat {
        PythonFloat(self.value.numerical_value().clone())
    }
    #[getter]
    fn precision(&self) -> CoefficientPrecisionWrapper {
        self.value.precision().clone().into()
    }
    #[getter]
    fn warnings(&self) -> Vec<PrecisionWarningWrapper> {
        self.value
            .warnings()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }
    #[getter]
    fn provenance(&self) -> DataProvenanceWrapper {
        self.value.provenance().clone().into()
    }
});

record!(
    MasterCoefficientWrapper,
    "MasterCoefficient",
    MasterCoefficient,
    {
        /// Encode achieved precision in a native expression; rejects uncertain zeros.
        fn to_expression(&self) -> PyResult<PythonExpression> {
            self.value
                .to_atom()
                .map(Into::into)
                .map_err(vakint_to_python_error)
        }
        #[getter]
        fn requested_decimal_precision(&self) -> u32 {
            self.value.requested_decimal_precision()
        }
        #[getter]
        fn symbol(&self) -> PythonExpression {
            self.value.symbol().clone().into()
        }
        #[getter]
        fn epsilon_power(&self) -> i64 {
            self.value.epsilon_power()
        }
        #[getter]
        fn available_precision(&self) -> AvailablePrecisionWrapper {
            self.value.available_precision().into()
        }
        #[getter]
        fn numerical_value(&self) -> PythonFloat {
            PythonFloat(self.value.numerical_value().clone())
        }
        #[getter]
        fn precision(&self) -> CoefficientPrecisionWrapper {
            self.value.precision().clone().into()
        }
        #[getter]
        fn warnings(&self) -> Vec<PrecisionWarningWrapper> {
            self.value
                .warnings()
                .iter()
                .cloned()
                .map(Into::into)
                .collect()
        }
        #[getter]
        fn provenance(&self) -> DataProvenanceWrapper {
            self.value.provenance().clone().into()
        }
        #[getter]
        fn symbolic_expression(&self) -> MasterSymbolicExpressionWrapper {
            self.value.symbolic_expression().clone().into()
        }
    }
);

/// All available coefficients and the authenticated graph of one declared master.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "VakintMaster",
    module = "symbolica.community.hepkit.vakint",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct VakintMasterWrapper {
    pub(super) value: VakintMaster,
}

impl From<VakintMaster> for VakintMasterWrapper {
    fn from(value: VakintMaster) -> Self {
        Self { value }
    }
}
impl super::ModuleInit for VakintMasterWrapper {}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl VakintMasterWrapper {
    #[getter]
    fn requested_decimal_precision(&self) -> u32 {
        self.value.requested_decimal_precision()
    }
    #[getter]
    fn symbol(&self) -> MasterSymbolWrapper {
        self.value.symbol().clone().into()
    }
    #[getter]
    fn master_id(&self) -> &str {
        self.value.master_id()
    }
    #[getter]
    fn family_sha256(&self) -> &str {
        self.value.family_sha256()
    }
    #[getter]
    fn loop_count(&self) -> usize {
        self.value.loop_count()
    }
    #[getter]
    fn normalization(&self) -> &str {
        self.value.normalization()
    }
    #[getter]
    fn provenance(&self) -> DataProvenanceWrapper {
        self.value.provenance().clone().into()
    }
    #[getter]
    fn graph(&self) -> PyFeynmanDiagram {
        self.value.graph().clone().into()
    }
    #[getter]
    fn coefficients(&self) -> BTreeMap<i64, MasterCoefficientWrapper> {
        self.value
            .coefficients()
            .iter()
            .map(|(power, value)| (*power, value.clone().into()))
            .collect()
    }
    #[getter]
    fn is_computed(&self) -> bool {
        self.value.is_computed()
    }
    fn __repr__(&self) -> String {
        self.value.to_string()
    }
    fn __str__(&self) -> String {
        self.__repr__()
    }
    /// Render coefficients and the normal FeynKit graph without reevaluating them.
    fn to_html(&self, py: Python<'_>) -> PyResult<String> {
        render(py, &self.value.presentation(), Some(self.value.graph()))
    }
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.to_html(py)
    }
}
