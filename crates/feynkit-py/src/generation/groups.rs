use super::*;

/// One diagram and its numerator ratio inside a grouped result.
///
/// A group member points into ``BuiltDiagramSet.diagrams`` and expresses its
/// numerator relative to the group's master diagram.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> result = process.amplitude(loops=hep.Loops.total(1)).generate().build()
/// >>> diagram = result.diagrams[0]
/// >>> group = result.groups[0]
/// >>> member = group.members[0]
/// >>> diagram = result.diagrams[member.diagram]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GroupMember",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGroupMember {
    // Box is required for alignment: the source ID contains a u128 requiring
    // 16-byte alignment, but CPython's wasm32 object allocator guarantees only 8.
    // Rust allocates the boxed payload with the required alignment; do not inline it.
    pub(super) inner: Box<GroupMember>,
}

#[cfg(target_arch = "wasm32")]
const _: () = {
    assert!(std::mem::align_of::<PyGroupMember>() <= 8);
    assert!(
        std::mem::align_of::<<PyGroupMember as pyo3::impl_::pyclass::PyClassImpl>::Layout>() <= 8
    );
};

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyGroupMember {
    /// Return the generated-order index from before zero-numerator removal.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``GroupMember`` class example:
    ///
    /// >>> source_index = member.source_diagram
    /// >>> source_identity = (source_index, member.source_id, member.source_name)
    #[getter]
    fn source_diagram(&self) -> usize {
        self.inner.source_diagram
    }
    /// Return the requested source ID, or None for minimal generation output.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``GroupMember`` class example:
    ///
    /// >>> source_ids = [entry.source_id for entry in group.members]
    #[getter]
    fn source_id(&self) -> Option<String> {
        self.inner
            .source_identity
            .as_ref()
            .map(|identity| identity.id.to_string())
    }
    /// Return the requested source name, or None for minimal generation output.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``GroupMember`` class example:
    ///
    /// >>> source_names = [entry.source_name for entry in group.members]
    #[getter]
    fn source_name(&self) -> Option<&str> {
        self.inner
            .source_identity
            .as_ref()
            .map(|identity| identity.name.as_str())
    }
    /// Return the collapsed master index in ``BuiltDiagramSet.diagrams``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``GroupMember`` class example:
    ///
    /// >>> result.diagrams[result.groups[0].members[0].diagram]
    #[getter]
    fn diagram(&self) -> usize {
        self.inner.diagram
    }
    /// Return the member numerator divided by the group master numerator.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``GroupMember`` class example:
    ///
    /// >>> relative_weights = [entry.ratio for entry in group.members]
    #[getter]
    fn ratio(&self) -> String {
        self.inner.ratio.to_plain_string()
    }
    /// Parse the numerator ratio as a native Symbolica expression.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``GroupMember`` class example:
    ///
    /// >>> ratio = result.groups[0].members[0].ratio_expression()
    fn ratio_expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.ratio.clone(),
        }
    }
    /// Parse the source diagram's numerator-independent factor.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``GroupMember`` class example:
    ///
    /// >>> group = result.groups[0]
    /// >>> member = group.members[0]
    /// >>> master = result.diagrams[group.master]
    /// >>> reconstructed = (member.overall_factor_expression()
    /// ...                  * member.ratio_expression()
    /// ...                  * master.numerator_expression())
    /// >>> reconstructed
    fn overall_factor_expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.overall_factor.clone(),
        }
    }
}

/// Diagrams whose numerators are related by known scalar factors.
///
/// Grouping lets amplitude calculations evaluate one master numerator and
/// reconstruct related diagrams from each member's symbolic ratio.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> result = process.amplitude(loops=hep.Loops.total(1)).generate().build()
/// >>> diagram = result.diagrams[0]
/// >>> group = result.groups[0]
/// >>> master_diagram = result.diagrams[group.master]
/// >>> ratios = [member.ratio for member in group.members]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramGroup",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramGroup {
    pub(super) inner: DiagramGroup,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyDiagramGroup {
    /// Return the retained-diagram index used as the numerator reference.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramGroup`` class example:
    ///
    /// >>> result.diagrams[result.groups[0].master]
    #[getter]
    fn master(&self) -> usize {
        self.inner.master
    }
    /// Return the deterministically ordered group members, including the master.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``DiagramGroup`` class example:
    ///
    /// >>> weighted_diagrams = [(result.diagrams[m.diagram], m.ratio) for m in group.members]
    #[getter]
    fn members(&self) -> Vec<PyGroupMember> {
        self.inner
            .members
            .iter()
            .cloned()
            .map(|inner| PyGroupMember {
                inner: Box::new(inner),
            })
            .collect()
    }
}
