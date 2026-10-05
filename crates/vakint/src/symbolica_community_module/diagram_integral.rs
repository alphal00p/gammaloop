use feynkit_py::{PyFeynmanDiagram, PyIntegralFamily};
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyInt, PyMapping, PyMappingMethods},
};
use symbolica::{api::python::PythonExpression, domains::integer::Integer};

use super::{VakintExpressionWrapper, vakint_to_python_error};
use crate::{VakintExpression, diagram_integral::DiagramIntegralOptions};

/// Build a native VakintExpression from a routed, uncut vacuum graph.
///
/// The family must retain the ordered physical denominator prefix in its stored
/// loop basis; auxiliary powers must be nonpositive. Numerators use native
/// Kinematics.scalar_product notation. Spectator order defines numerical vector
/// IDs. Parameter substitutions are simultaneous and must already be reflected
/// in the family's physical denominators. No graph factor is inserted.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
#[pyo3(signature = (diagram, family, numerator, *, powers=None, parameter_substitutions=None, external_momenta=Vec::new()), text_signature = "(diagram, family, numerator, *, powers=None, parameter_substitutions=None, external_momenta=())")]
pub(crate) fn integral_from_diagram(
    diagram: &PyFeynmanDiagram,
    family: &PyIntegralFamily,
    numerator: PythonExpression,
    #[gen_stub(override_type(type_repr="typing.Sequence[int] | None", imports=("typing")))]
    powers: Option<&Bound<'_, PyAny>>,
    #[gen_stub(override_type(type_repr="typing.Mapping[symbolica.Expression, symbolica.Expression] | None", imports=("typing", "symbolica")))]
    parameter_substitutions: Option<&Bound<'_, PyMapping>>,
    #[gen_stub(override_type(type_repr="typing.Sequence[symbolica.Expression]", imports=("typing", "symbolica")))]
    external_momenta: Vec<PythonExpression>,
) -> PyResult<VakintExpressionWrapper> {
    super::record_usage();
    let powers = powers
        .map(|values| {
            values
                .try_iter()?
                .map(|value| {
                    let value = value?;
                    if !value.is_exact_instance_of::<PyInt>() {
                        return Err(PyValueError::new_err(
                            "powers must be one integer (not bool) per family denominator",
                        ));
                    }
                    value.extract::<Integer>()
                })
                .collect::<PyResult<Vec<_>>>()
        })
        .transpose()?;
    let parameter_substitutions = parameter_substitutions
        .map(|values| {
            values
                .items()?
                .try_iter()?
                .map(|item| {
                    let (source, target) =
                        item?.extract::<(PythonExpression, PythonExpression)>()?;
                    Ok((source.expr, target.expr))
                })
                .collect::<PyResult<Vec<_>>>()
        })
        .transpose()?
        .unwrap_or_default();
    let options = DiagramIntegralOptions {
        powers,
        parameter_substitutions,
        external_momenta: external_momenta
            .into_iter()
            .map(|value| value.expr)
            .collect(),
    };
    Ok(VakintExpressionWrapper {
        value: VakintExpression::from_diagram(
            diagram.as_diagram()?,
            family.as_family(),
            &numerator.expr,
            &options,
        )
        .map_err(vakint_to_python_error)?,
    })
}
