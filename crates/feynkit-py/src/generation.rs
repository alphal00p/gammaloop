use std::{collections::BTreeMap, sync::Arc};

use feynkit_generator::{
    CancellationToken, DiagramGroup, GenerationFilter, GenerationOptions, GenerationReport,
    GenerationResult, GenerationType, Generator, GraphGroupingOptions, GroupMember,
    NumeratorGrouping, ParticleSelector, Process, SelfEnergyFilterOptions, SewnFilterOptions,
    SnailFilterOptions, TadpoleFilterOptions,
};
use feynkit_model::Model;
use pyo3::{
    FromPyObject,
    prelude::*,
    types::{PyAny, PyModule},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods},
};

use crate::{
    error,
    graph::{PyFeynmanDiagram, parse_symbolic_annotation},
    model::PyModel,
};
use symbolica::api::python::PythonExpression;

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "GenerationType",
    module = "symbolica.community.feynkit",
    rename_all = "SCREAMING_SNAKE_CASE",
    frozen,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyGenerationType {
    Amplitude,
    CrossSection,
}

impl From<GenerationType> for PyGenerationType {
    fn from(value: GenerationType) -> Self {
        match value {
            GenerationType::Amplitude => Self::Amplitude,
            GenerationType::CrossSection => Self::CrossSection,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationType {
    #[getter]
    fn value(&self) -> &'static str {
        match self {
            Self::Amplitude => "amplitude",
            Self::CrossSection => "cross_section",
        }
    }

    fn __str__(&self) -> &'static str {
        self.value()
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other) = other.cast::<Self>() {
            self == other.get()
        } else {
            other
                .extract::<String>()
                .is_ok_and(|other| self.value() == other)
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ParticleSelector",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PyParticleSelector {
    inner: ParticleSelector,
}

impl From<ParticleSelector> for PyParticleSelector {
    fn from(inner: ParticleSelector) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParticleSelector {
    #[staticmethod]
    fn by_name(name: String) -> Self {
        Self {
            inner: ParticleSelector::Name(name),
        }
    }

    #[staticmethod]
    fn by_pdg(pdg: i64) -> Self {
        Self {
            inner: ParticleSelector::Pdg(pdg),
        }
    }

    #[getter]
    fn name(&self) -> Option<&str> {
        match &self.inner {
            ParticleSelector::Name(name) => Some(name),
            ParticleSelector::Pdg(_) => None,
        }
    }

    #[getter]
    fn pdg(&self) -> Option<i64> {
        match &self.inner {
            ParticleSelector::Name(_) => None,
            ParticleSelector::Pdg(pdg) => Some(*pdg),
        }
    }

    #[getter]
    fn is_name(&self) -> bool {
        matches!(&self.inner, ParticleSelector::Name(_))
    }

    #[getter]
    fn is_pdg(&self) -> bool {
        matches!(&self.inner, ParticleSelector::Pdg(_))
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            ParticleSelector::Name(name) => format!("ParticleSelector.by_name('{name}')"),
            ParticleSelector::Pdg(pdg) => format!("ParticleSelector.by_pdg({pdg})"),
        }
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other) = other.cast::<Self>() {
            self == other.get()
        } else if let Ok(name) = other.extract::<String>() {
            self.inner == ParticleSelector::Name(name)
        } else {
            other
                .extract::<i64>()
                .is_ok_and(|pdg| self.inner == ParticleSelector::Pdg(pdg))
        }
    }
}

#[derive(FromPyObject)]
enum SelectorInput {
    Selector(PyParticleSelector),
    Name(String),
    Pdg(i64),
}

impl From<SelectorInput> for ParticleSelector {
    fn from(value: SelectorInput) -> Self {
        match value {
            SelectorInput::Selector(selector) => selector.inner,
            SelectorInput::Name(name) => Self::Name(name),
            SelectorInput::Pdg(pdg) => Self::Pdg(pdg),
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for SelectorInput {
    fn type_input() -> TypeInfo {
        PyParticleSelector::type_input() | String::type_input() | i64::type_input()
    }

    fn type_output() -> TypeInfo {
        PyParticleSelector::type_output() | String::type_output() | i64::type_output()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Process",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyProcess {
    pub(crate) inner: Process,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyProcess {
    #[staticmethod]
    fn amplitude(incoming: Vec<SelectorInput>, outgoing: Vec<SelectorInput>) -> Self {
        Self {
            inner: Process::amplitude(
                incoming.into_iter().map(ParticleSelector::from),
                outgoing.into_iter().map(ParticleSelector::from),
            ),
        }
    }

    #[staticmethod]
    fn cross_section(incoming: Vec<SelectorInput>, outgoing: Vec<SelectorInput>) -> Self {
        Self {
            inner: Process::cross_section(
                incoming.into_iter().map(ParticleSelector::from),
                outgoing.into_iter().map(ParticleSelector::from),
            ),
        }
    }

    fn with_loop_count(&self, minimum: usize, maximum: usize) -> PyResult<Self> {
        self.inner
            .clone()
            .with_loop_count(minimum, maximum)
            .map(|inner| Self { inner })
            .map_err(error::process)
    }

    fn with_final_state_alternatives(
        &self,
        alternatives: Vec<Vec<SelectorInput>>,
    ) -> PyResult<Self> {
        self.inner
            .clone()
            .with_final_state_alternatives(
                alternatives
                    .into_iter()
                    .map(|state| state.into_iter().map(ParticleSelector::from)),
            )
            .map(|inner| Self { inner })
            .map_err(error::process)
    }

    #[pyo3(signature = (*, initial=false, final_state=false, left_right=false, external_fermions=false))]
    fn with_symmetrization(
        &self,
        initial: bool,
        final_state: bool,
        left_right: bool,
        external_fermions: bool,
    ) -> Self {
        Self {
            inner: self
                .inner
                .clone()
                .symmetrize_initial(initial)
                .symmetrize_final(final_state)
                .symmetrize_left_right(left_right)
                .symmetrize_external_fermions(external_fermions),
        }
    }

    #[getter]
    fn generation_type(&self) -> PyGenerationType {
        self.inner.generation_type().into()
    }

    #[getter]
    fn incoming(&self) -> Vec<PyParticleSelector> {
        self.inner
            .incoming()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn outgoing_alternatives(&self) -> Vec<Vec<PyParticleSelector>> {
        self.inner
            .outgoing_alternatives()
            .iter()
            .map(|state| state.iter().cloned().map(Into::into).collect())
            .collect()
    }

    #[getter]
    fn loop_count(&self) -> (usize, usize) {
        let range = self.inner.loop_count();
        (*range.start(), *range.end())
    }

    #[getter]
    fn symmetrizes_initial(&self) -> bool {
        self.inner.symmetrizes_initial()
    }

    #[getter]
    fn symmetrizes_final(&self) -> bool {
        self.inner.symmetrizes_final()
    }

    #[getter]
    fn symmetrizes_left_right(&self) -> bool {
        self.inner.symmetrizes_left_right()
    }

    #[getter]
    fn symmetrizes_external_fermions(&self) -> bool {
        self.inner.symmetrizes_external_fermions()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CancellationToken",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyCancellationToken {
    inner: CancellationToken,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCancellationToken {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn cancel(&self) {
        self.inner.cancel();
    }

    #[getter]
    fn is_cancelled(&self) -> bool {
        self.inner.is_cancelled()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationOptions",
    module = "symbolica.community.feynkit",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyGenerationOptions {
    pub(crate) inner: GenerationOptions,
}

impl PyGenerationOptions {
    fn add_graph_filter(&mut self, filter: GenerationFilter) {
        self.inner = self.inner.clone().with_graph_filter(filter);
    }

    fn add_cut_amplitude_filter(&mut self, filter: GenerationFilter) {
        self.inner = self.inner.clone().with_cut_amplitude_filter(filter);
    }

    fn grouping_options(
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) -> GraphGroupingOptions {
        GraphGroupingOptions {
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationOptions {
    #[new]
    #[pyo3(signature = (*, threads=None, max_vertices=None, allow_self_loops=false, allow_zero_flow_edges=false, graph_prefix=None))]
    fn new(
        threads: Option<usize>,
        max_vertices: Option<usize>,
        allow_self_loops: bool,
        allow_zero_flow_edges: bool,
        graph_prefix: Option<String>,
    ) -> Self {
        let mut inner = GenerationOptions::default()
            .allow_self_loops(allow_self_loops)
            .allow_zero_flow_edges(allow_zero_flow_edges);
        if let Some(threads) = threads {
            inner = inner.threads(threads);
        }
        if let Some(max_vertices) = max_vertices {
            inner = inner.max_vertices(max_vertices);
        }
        if let Some(prefix) = graph_prefix {
            inner = inner.graph_prefix(prefix);
        }
        Self { inner }
    }

    fn set_cancellation_token(&mut self, token: &PyCancellationToken) {
        self.inner = self.inner.clone().cancellation_token(token.inner.clone());
    }

    fn add_particle_veto(&mut self, pdg_codes: Vec<i64>) {
        self.add_graph_filter(GenerationFilter::ParticleVeto(pdg_codes));
    }

    fn add_vertex_allow(&mut self, vertices: Vec<String>) {
        self.add_graph_filter(GenerationFilter::VertexAllow(vertices));
    }

    fn add_vertex_veto(&mut self, vertices: Vec<String>) {
        self.add_graph_filter(GenerationFilter::VertexVeto(vertices));
    }

    fn set_maximum_bridges(&mut self, maximum: usize) {
        self.add_graph_filter(GenerationFilter::MaxNumberOfBridges(maximum));
    }

    #[pyo3(signature = (*, veto_massive=true, veto_massless=true, only_scaleless=false))]
    fn set_self_energy_filter(
        &mut self,
        veto_massive: bool,
        veto_massless: bool,
        only_scaleless: bool,
    ) {
        self.add_graph_filter(GenerationFilter::SelfEnergy(SelfEnergyFilterOptions {
            veto_massive,
            veto_massless,
            only_scaleless,
        }));
    }

    #[pyo3(signature = (*, veto_attached_to_massive=true, veto_attached_to_massless=true, only_scaleless=false))]
    fn set_tadpole_filter(
        &mut self,
        veto_attached_to_massive: bool,
        veto_attached_to_massless: bool,
        only_scaleless: bool,
    ) {
        self.add_graph_filter(GenerationFilter::Tadpoles(TadpoleFilterOptions {
            veto_attached_to_massive,
            veto_attached_to_massless,
            only_scaleless,
        }));
    }

    #[pyo3(signature = (*, veto_attached_to_massive=false, veto_attached_to_massless=true, only_scaleless=false))]
    fn set_zero_snail_filter(
        &mut self,
        veto_attached_to_massive: bool,
        veto_attached_to_massless: bool,
        only_scaleless: bool,
    ) {
        self.add_graph_filter(GenerationFilter::ZeroSnails(SnailFilterOptions {
            veto_attached_to_massive,
            veto_attached_to_massless,
            only_scaleless,
        }));
    }

    fn set_coupling_orders(&mut self, orders: BTreeMap<String, (usize, Option<usize>)>) {
        self.add_graph_filter(GenerationFilter::CouplingOrders(orders));
    }

    fn set_loop_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::LoopCountRange((minimum, maximum)));
    }

    fn set_fermion_loop_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::FermionLoopCountRange((minimum, maximum)));
    }

    fn set_factorized_loop_topologies_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::FactorizedLoopTopologiesCountRange((
            minimum, maximum,
        )));
    }

    fn set_blob_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::BlobRange(minimum..=maximum));
    }

    fn set_spectator_range(&mut self, minimum: usize, maximum: usize) {
        self.add_graph_filter(GenerationFilter::SpectatorRange(minimum..=maximum));
    }

    fn set_perturbative_orders(&mut self, orders: BTreeMap<String, usize>) {
        self.add_graph_filter(GenerationFilter::PerturbativeOrders(orders));
    }

    #[pyo3(signature = (*, filter_tadpoles=true))]
    fn set_sewn_filter(&mut self, filter_tadpoles: bool) {
        self.add_graph_filter(GenerationFilter::Sewn(SewnFilterOptions {
            filter_tadpoles,
        }));
    }

    fn set_cut_amplitude_coupling_orders(
        &mut self,
        orders: BTreeMap<String, (usize, Option<usize>)>,
    ) {
        self.add_cut_amplitude_filter(GenerationFilter::CouplingOrders(orders));
    }

    fn set_cut_amplitude_loop_count_range(&mut self, minimum: usize, maximum: usize) {
        self.add_cut_amplitude_filter(GenerationFilter::LoopCountRange((minimum, maximum)));
    }

    fn disable_numerator_grouping(&mut self) {
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::None);
    }

    fn detect_zero_numerators(&mut self) {
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::OnlyDetectZeroes);
    }

    #[pyo3(signature = (*, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false))]
    fn group_identical_numerators(
        &mut self,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) {
        let options = Self::grouping_options(
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        );
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::Identical(options));
    }

    #[pyo3(signature = (*, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false))]
    fn group_numerators_up_to_sign(
        &mut self,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) {
        let options = Self::grouping_options(
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        );
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::UpToSign(options));
    }

    #[pyo3(signature = (*, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false))]
    fn group_numerators_up_to_scalar(
        &mut self,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
    ) {
        let options = Self::grouping_options(
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
        );
        self.inner = self
            .inner
            .clone()
            .numerator_grouping(NumeratorGrouping::UpToScalar(options));
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationReport",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerationReport {
    inner: GenerationReport,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationReport {
    #[getter]
    fn topology_count(&self) -> usize {
        self.inner.topology_count
    }
    #[getter]
    fn interaction_assignment_count(&self) -> usize {
        self.inner.interaction_assignment_count
    }
    #[getter]
    fn retained_count(&self) -> usize {
        self.inner.retained_count
    }
    #[getter]
    fn zero_numerator_count(&self) -> usize {
        self.inner.zero_numerator_count
    }
    #[getter]
    fn completed(&self) -> bool {
        self.inner.completed
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GroupMember",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGroupMember {
    inner: GroupMember,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGroupMember {
    #[getter]
    fn source_diagram(&self) -> usize {
        self.inner.source_diagram
    }
    #[getter]
    fn diagram(&self) -> usize {
        self.inner.diagram
    }
    #[getter]
    fn ratio(&self) -> &str {
        &self.inner.ratio
    }
    fn ratio_expression(&self) -> PyResult<PythonExpression> {
        parse_symbolic_annotation(&self.inner.ratio).map_err(error::GenerationError::new_err)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramGroup",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramGroup {
    inner: DiagramGroup,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramGroup {
    #[getter]
    fn master(&self) -> usize {
        self.inner.master
    }
    #[getter]
    fn members(&self) -> Vec<PyGroupMember> {
        self.inner
            .members
            .iter()
            .cloned()
            .map(|inner| PyGroupMember { inner })
            .collect()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationResult",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerationResult {
    inner: GenerationResult,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerationResult {
    #[getter]
    fn diagrams(&self) -> Vec<PyFeynmanDiagram> {
        self.inner
            .diagrams
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn groups(&self) -> Vec<PyDiagramGroup> {
        self.inner
            .groups
            .iter()
            .cloned()
            .map(|inner| PyDiagramGroup { inner })
            .collect()
    }

    #[getter]
    fn report(&self) -> PyGenerationReport {
        PyGenerationReport {
            inner: self.inner.report.clone(),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.diagrams.len()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Generator",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerator {
    inner: Generator,
    model: Arc<Model>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyGenerator {
    #[new]
    fn new(model: &PyModel) -> Self {
        Self {
            inner: Generator::new(model.inner.clone()),
            model: model.inner.clone(),
        }
    }

    #[getter]
    fn model(&self) -> PyModel {
        self.model.clone().into()
    }

    #[pyo3(signature = (process, options=None))]
    fn generate(
        &self,
        py: Python<'_>,
        process: &PyProcess,
        options: Option<&PyGenerationOptions>,
    ) -> PyResult<PyGenerationResult> {
        let generator = self.inner.clone();
        let process = process.inner.clone();
        let options = options.cloned().unwrap_or_default().inner;
        py.detach(move || generator.generate(&process, &options))
            .map(|inner| PyGenerationResult { inner })
            .map_err(error::generation)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGenerationType>()?;
    module.add_class::<PyParticleSelector>()?;
    module.add_class::<PyProcess>()?;
    module.add_class::<PyCancellationToken>()?;
    module.add_class::<PyGenerationOptions>()?;
    module.add_class::<PyGenerationReport>()?;
    module.add_class::<PyGroupMember>()?;
    module.add_class::<PyDiagramGroup>()?;
    module.add_class::<PyGenerationResult>()?;
    module.add_class::<PyGenerator>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;

    use super::*;

    #[test]
    fn process_metadata_and_selectors_are_typed_and_round_trip() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
            crate::initialize_feynkit(&module).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            let code = CString::new(
                r#"
by_name = fk.ParticleSelector.by_name("1")
by_pdg = fk.ParticleSelector.by_pdg(1)

assert by_name.name == "1"
assert by_name.pdg is None
assert by_name.is_name and not by_name.is_pdg
assert by_pdg.name is None
assert by_pdg.pdg == 1
assert by_pdg.is_pdg and not by_pdg.is_name
assert by_name == "1"
assert by_pdg == 1
assert by_name != by_pdg
try:
    by_pdg.pdg = 2
except AttributeError:
    pass
else:
    raise AssertionError("particle selectors must be immutable")

process = fk.Process.amplitude(
    [by_name, by_pdg],
    [fk.ParticleSelector.by_name("out"), fk.ParticleSelector.by_pdg(-1)],
).with_loop_count(1, 2)
assert process.generation_type == fk.GenerationType.AMPLITUDE
assert process.generation_type == "amplitude"
assert isinstance(process.generation_type, fk.GenerationType)
assert all(isinstance(selector, fk.ParticleSelector) for selector in process.incoming)
assert process.incoming == [by_name, by_pdg]
assert process.outgoing_alternatives[0][0].name == "out"
assert process.outgoing_alternatives[0][1].pdg == -1

round_tripped = fk.Process.amplitude(
    process.incoming,
    process.outgoing_alternatives[0],
).with_loop_count(*process.loop_count)
assert round_tripped.generation_type == process.generation_type
assert round_tripped.incoming == process.incoming
assert round_tripped.outgoing_alternatives == process.outgoing_alternatives
assert round_tripped.loop_count == process.loop_count

cross_section = fk.Process.cross_section([by_pdg], [by_name])
cross_section = cross_section.with_final_state_alternatives([[by_name], [by_pdg]])
cross_section = cross_section.with_symmetrization(
    initial=True,
    final_state=True,
    left_right=True,
    external_fermions=True,
)
assert cross_section.generation_type == fk.GenerationType.CROSS_SECTION
assert cross_section.generation_type.value == "cross_section"
assert cross_section.outgoing_alternatives == [[by_name], [by_pdg]]
assert cross_section.symmetrizes_initial
assert cross_section.symmetrizes_final
assert cross_section.symmetrizes_left_right
assert cross_section.symmetrizes_external_fermions
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
