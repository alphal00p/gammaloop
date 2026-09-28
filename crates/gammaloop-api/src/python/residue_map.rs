//! Detached Python views of the persisted native three-dimensional expression.

use std::collections::HashMap;

use gammalooprs::cff::{esurface::Esurface, hsurface::Hsurface};
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use pyo3::{
    prelude::*,
    types::{PyDict, PyList, PyTuple},
};
use symbolica::{atom::AtomCore, domains::rational::Rational};
use three_dimensional_reps::{
    expression::{CFFVariant, OrientationData, OrientationExpression},
    generation::CffEnergyFactorOwnership,
    surface::{HybridSurfaceID, LinearSurfaceKind, SurfaceOrigin},
    GeneratedThreeDExpression, LinearEnergyExpr,
};

use super::py_builtin_from_json_value;

fn fraction<'py>(py: Python<'py>, value: &Rational) -> PyResult<Bound<'py, PyAny>> {
    // Fraction parses arbitrary-size integer ratios exactly, without a float boundary.
    py.import("fractions")?
        .getattr("Fraction")?
        .call1((value.to_string(),))
}

fn energy_terms<'py>(
    py: Python<'py>,
    terms: &[(EdgeIndex, Rational)],
) -> PyResult<Bound<'py, PyTuple>> {
    let terms = terms
        .iter()
        .map(|(edge, coefficient)| Ok((edge.0, fraction(py, coefficient)?)))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, terms)
}

fn surface_reference(surface: HybridSurfaceID) -> (&'static str, Option<usize>) {
    match surface {
        HybridSurfaceID::Esurface(id) => ("esurface", Some(id.0)),
        HybridSurfaceID::Hsurface(id) => ("hsurface", Some(id.0)),
        HybridSurfaceID::Linear(id) => ("linear", Some(id.0)),
        HybridSurfaceID::Unit => ("unit", None),
        HybridSurfaceID::Infinite => ("infinite", None),
    }
}

fn energy_factor_ownership(ownership: CffEnergyFactorOwnership) -> &'static str {
    match ownership {
        CffEnergyFactorOwnership::GlobalSourceProduct => "global_source_product",
        CffEnergyFactorOwnership::VariantLocal => "variant_local",
    }
}

/// Immutable exact affine energy expression, retaining native term order.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "LinearEnergyExpression", frozen, eq, hash, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PyLinearEnergyExpression {
    inner: LinearEnergyExpr,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyLinearEnergyExpression {
    /// Ordered ``(internal_edge_id, Fraction)`` terms multiplying on-shell energies.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[tuple[int, fractions.Fraction], ...]", imports = ("fractions")))]
    fn internal_terms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        energy_terms(py, &self.inner.internal_terms)
    }

    /// Ordered ``(external_edge_id, Fraction)`` terms multiplying external energies.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[tuple[int, fractions.Fraction], ...]", imports = ("fractions")))]
    fn external_terms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        energy_terms(py, &self.inner.external_terms)
    }

    /// Exact coefficient of the independent numerator sampling scale M.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "fractions.Fraction", imports = ("fractions")))]
    fn uniform_scale_coeff<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        fraction(py, &self.inner.uniform_scale_coeff)
    }

    /// Exact constant energy shift.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "fractions.Fraction", imports = ("fractions")))]
    fn constant<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        fraction(py, &self.inner.constant)
    }

    /// Native symbolic form for diagnostics; equality uses the exact stored terms.
    #[getter]
    fn canonical_string(&self) -> String {
        self.inner.to_atom(&[]).to_canonical_string()
    }
}

/// Immutable native residue identity, independent of IDs, labels and scalar variants.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "ResidueMapKey", frozen, eq, hash, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PyResidueMapKey {
    directions: EdgeVec<Orientation>,
    loop_energy_map: Vec<LinearEnergyExpr>,
    edge_energy_map: Vec<LinearEnergyExpr>,
}

impl From<&OrientationExpression> for PyResidueMapKey {
    fn from(expression: &OrientationExpression) -> Self {
        Self {
            directions: expression.data.orientation.clone(),
            loop_energy_map: expression.loop_energy_map.clone(),
            edge_energy_map: expression.edge_energy_map.clone(),
        }
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyResidueMapKey {
    /// Complete native edge-direction vector: default 1, reversed -1, undirected 0.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int, ...]"))]
    fn directions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            self.directions
                .iter()
                .map(|(_, direction)| match direction {
                    Orientation::Default => 1,
                    Orientation::Reversed => -1,
                    Orientation::Undirected => 0,
                }),
        )
    }

    /// Exact affine energies in native loop-coordinate order.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[LinearEnergyExpression, ...]"))]
    fn loop_energy_map<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let maps = self
            .loop_energy_map
            .iter()
            .map(|inner| {
                Py::new(
                    py,
                    PyLinearEnergyExpression {
                        inner: inner.clone(),
                    },
                )
            })
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, maps)
    }

    /// Exact affine energies in native internal-edge order.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[LinearEnergyExpression, ...]"))]
    fn edge_energy_map<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let maps = self
            .edge_energy_map
            .iter()
            .map(|inner| {
                Py::new(
                    py,
                    PyLinearEnergyExpression {
                        inner: inner.clone(),
                    },
                )
            })
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, maps)
    }

    /// The native ``residue_map_key()`` diagnostic string.
    #[getter]
    fn canonical_string(&self) -> String {
        OrientationExpression {
            data: OrientationData::new(self.directions.clone()),
            loop_energy_map: self.loop_energy_map.clone(),
            edge_energy_map: self.edge_energy_map.clone(),
            variants: Vec::new(),
        }
        .residue_map_key()
    }
}

/// One scalar variant of a native residue, with unexpanded denominator topology.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "ResidueVariant", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct PyResidueVariant {
    inner: CFFVariant,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyResidueVariant {
    #[getter]
    fn origin(&self) -> Option<String> {
        self.inner.origin.clone()
    }

    /// Exact scalar prefactor, including the separately recorded absorbed signs.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "fractions.Fraction", imports = ("fractions")))]
    fn prefactor<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        fraction(py, &self.inner.prefactor)
    }

    /// Ordered edge occurrences supplying factors ``1/(2 E_edge)``; repetitions remain.
    #[getter]
    fn half_edges(&self) -> Vec<usize> {
        self.inner.half_edges.iter().map(|edge| edge.0).collect()
    }

    /// Ordered original denominator-edge occurrences, including repeated edges.
    #[getter]
    fn denominator_edges(&self) -> Vec<usize> {
        self.inner
            .denominator_edges
            .iter()
            .map(|edge| edge.0)
            .collect()
    }

    /// Absorbed orientation signs keyed by tagged surface references.
    #[getter]
    fn denominator_surface_signs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        for (surface, sign) in &self.inner.denominator_surface_signs {
            result.set_item(surface_reference(*surface), *sign)?;
        }
        Ok(result)
    }

    /// Absorbed routing signs keyed by ordered tuples of denominator-edge support.
    #[getter]
    fn denominator_edge_support_signs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        for (edges, sign) in &self.inner.denominator_edge_support_signs {
            result.set_item(PyTuple::new(py, edges.iter().map(|edge| edge.0))?, *sign)?;
        }
        Ok(result)
    }

    /// Power of ``1/M`` multiplying this variant.
    #[getter]
    fn uniform_scale_power(&self) -> usize {
        self.inner.uniform_scale_power
    }

    /// Ordered numerator-surface references, preserving multiplicities.
    #[getter]
    fn numerator_surfaces(&self) -> Vec<(&'static str, Option<usize>)> {
        self.inner
            .numerator_surfaces
            .iter()
            .copied()
            .map(surface_reference)
            .collect()
    }

    /// Native tree: root ID (or None) and nodes with surface, children and parent IDs.
    #[getter]
    fn denominator<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let tree = PyDict::new(py);
        let nodes = PyList::empty(py);
        for node in self.inner.denominator.iter_nodes() {
            let record = PyDict::new(py);
            record.set_item("node_id", node.node_id.0)?;
            record.set_item("surface", surface_reference(node.data))?;
            record.set_item(
                "children",
                node.children.iter().map(|id| id.0).collect::<Vec<_>>(),
            )?;
            record.set_item("parent", node.parent.map(|id| id.0))?;
            nodes.append(record)?;
        }
        tree.set_item(
            "root",
            self.inner
                .denominator
                .iter_nodes()
                .find(|node| node.parent.is_none())
                .map(|node| node.node_id.0),
        )?;
        tree.set_item("nodes", nodes)?;
        Ok(tree)
    }
}

/// One expression-local residue entry; several native IDs can share the same key.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Residue", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct PyResidue {
    native_id: usize,
    inner: OrientationExpression,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyResidue {
    /// Native expression-local ID, independent of runtime execution slots.
    #[getter]
    fn native_id(&self) -> usize {
        self.native_id
    }
    #[getter]
    fn label(&self) -> Option<String> {
        self.inner.data.label.clone()
    }
    #[getter]
    fn numerator_map_index(&self) -> Option<usize> {
        self.inner.data.numerator_map_index
    }
    #[getter]
    fn variants(&self) -> Vec<PyResidueVariant> {
        self.inner
            .variants
            .iter()
            .map(|inner| PyResidueVariant {
                inner: inner.clone(),
            })
            .collect()
    }
}

/// Detached snapshot of one graph's persisted native 3D representation.
///
/// This is the generated graph-level residue map, not the evaluated UV/threshold
/// catalogue or a partition into independently integrable LTD contributions.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "ResidueMap", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct PyResidueMap {
    graph_name: String,
    generated: GeneratedThreeDExpression<Esurface, Hsurface>,
}

impl PyResidueMap {
    pub(crate) fn new(
        graph_name: String,
        generated: &GeneratedThreeDExpression<Esurface, Hsurface>,
    ) -> Self {
        Self {
            graph_name,
            generated: generated.clone(),
        }
    }

    fn grouped_entries(&self) -> Vec<(PyResidueMapKey, Vec<PyResidue>)> {
        let mut indices = HashMap::new();
        let mut entries: Vec<(PyResidueMapKey, Vec<PyResidue>)> = Vec::new();
        for (id, expression) in self.generated.expression.orientations.iter_enumerated() {
            let key = PyResidueMapKey::from(expression);
            let index = *indices.entry(key.clone()).or_insert_with(|| {
                entries.push((key, Vec::new()));
                entries.len() - 1
            });
            entries[index].1.push(PyResidue {
                native_id: id.0,
                inner: expression.clone(),
            });
        }
        entries
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyResidueMap {
    #[getter]
    fn graph_name(&self) -> String {
        self.graph_name.clone()
    }
    #[getter]
    fn representation(&self) -> String {
        self.generated.representation.to_string()
    }

    /// New dict keyed by immutable native residue identities; every native ID is retained.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "dict[ResidueMapKey, list[Residue]]"))]
    fn entries<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        for (key, residues) in self.grouped_entries() {
            result.set_item(Py::new(py, key)?, residues)?;
        }
        Ok(result)
    }

    /// Shared native surface records, keyed by tagged ``(kind, id)`` references.
    /// Unit and infinite surfaces use literal references with ID None.
    #[getter]
    fn surfaces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        let cache = &self.generated.expression.surfaces;
        for (id, surface) in cache.esurface_cache.iter_enumerated() {
            let record = PyDict::new(py);
            record.set_item("kind", "esurface")?;
            record.set_item(
                "energies",
                surface
                    .energies
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>(),
            )?;
            let shift = surface
                .external_shift
                .iter()
                .map(|(edge, sign)| (*edge, Rational::from(*sign)))
                .collect::<Vec<_>>();
            record.set_item("external_shift", energy_terms(py, &shift)?)?;
            let vertices = serde_json::to_value(surface.vertex_set)
                .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
            record.set_item("vertex_set", py_builtin_from_json_value(py, &vertices)?)?;
            result.set_item(surface_reference(HybridSurfaceID::Esurface(id)), record)?;
        }
        for (id, surface) in cache.hsurface_cache.iter_enumerated() {
            let record = PyDict::new(py);
            record.set_item("kind", "hsurface")?;
            record.set_item(
                "positive_energies",
                surface
                    .positive_energies
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>(),
            )?;
            record.set_item(
                "negative_energies",
                surface
                    .negative_energies
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>(),
            )?;
            let shift = surface
                .external_shift
                .iter()
                .map(|(edge, sign)| (*edge, Rational::from(*sign)))
                .collect::<Vec<_>>();
            record.set_item("external_shift", energy_terms(py, &shift)?)?;
            let vertices = serde_json::to_value(surface.vertex_set)
                .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
            record.set_item("vertex_set", py_builtin_from_json_value(py, &vertices)?)?;
            result.set_item(surface_reference(HybridSurfaceID::Hsurface(id)), record)?;
        }
        for (id, surface) in cache.linear_surface_cache.iter_enumerated() {
            let record = PyDict::new(py);
            record.set_item(
                "kind",
                match surface.kind {
                    LinearSurfaceKind::Esurface => "esurface",
                    LinearSurfaceKind::Hsurface => "hsurface",
                },
            )?;
            record.set_item(
                "origin",
                match surface.origin {
                    SurfaceOrigin::Physical => "physical",
                    SurfaceOrigin::Helper => "helper",
                },
            )?;
            record.set_item("numerator_only", surface.numerator_only)?;
            record.set_item(
                "expression",
                PyLinearEnergyExpression {
                    inner: surface.expression.clone(),
                },
            )?;
            result.set_item(surface_reference(HybridSurfaceID::Linear(id)), record)?;
        }
        for kind in ["unit", "infinite"] {
            let record = PyDict::new(py);
            record.set_item("kind", kind)?;
            result.set_item((kind, None::<usize>), record)?;
        }
        Ok(result)
    }

    /// Unintegrated four-dimensional denominator records, preserving edge ID and power.
    #[getter]
    fn residual_denominators<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let result = PyList::empty(py);
        for denominator in &self.generated.expression.residual_denominators {
            let record = PyDict::new(py);
            record.set_item("edge_id", denominator.edge_id.0)?;
            record.set_item("power", denominator.power)?;
            record.set_item("origin", denominator.origin.clone())?;
            result.append(record)?;
        }
        Ok(result)
    }

    #[getter]
    fn energy_factor_ownership(&self) -> &'static str {
        energy_factor_ownership(self.generated.energy_factor_ownership)
    }

    /// Source components with their edge IDs, energy ownership and both prefactor signs.
    #[getter]
    fn energy_factor_components<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let result = PyList::empty(py);
        for component in &self.generated.energy_factor_components {
            let record = PyDict::new(py);
            record.set_item("internal_edge_ids", component.internal_edge_ids.clone())?;
            record.set_item("ownership", energy_factor_ownership(component.ownership))?;
            record.set_item(
                "denominator_only_global_prefactor_sign",
                component.denominator_only_global_prefactor_sign.factor(),
            )?;
            record.set_item(
                "core_global_prefactor_sign",
                component.core_global_prefactor_sign.factor(),
            )?;
            result.append(record)?;
        }
        Ok(result)
    }

    #[getter]
    fn denominator_only_global_prefactor_sign(&self) -> i64 {
        self.generated
            .denominator_only_global_prefactor_sign
            .factor()
    }
    #[getter]
    fn core_global_prefactor_sign(&self) -> i64 {
        self.generated.core_global_prefactor_sign.factor()
    }
}

#[cfg(test)]
mod tests {
    use std::{
        collections::BTreeMap,
        hash::{DefaultHasher, Hash, Hasher},
    };

    use three_dimensional_reps::{
        expression::ThreeDExpression,
        generation::CffGlobalPrefactorSign,
        surface::{LinearSurfaceID, SurfaceCache},
        tree::{NodeId, Tree},
        RepresentationMode,
    };

    use super::*;

    fn native_residue() -> OrientationExpression {
        let large = Rational::from(i64::MAX) * Rational::from(i64::MAX) / Rational::from(7);
        let energy = LinearEnergyExpr {
            internal_terms: vec![(EdgeIndex(1), large)],
            external_terms: vec![(EdgeIndex(0), Rational::from((2, 3)))],
            uniform_scale_coeff: Rational::from((5, 7)),
            constant: Rational::from((-11, 13)),
        };
        let surface = HybridSurfaceID::Linear(LinearSurfaceID(0));
        let mut denominator = Tree::from_root(HybridSurfaceID::Unit);
        denominator.insert_node(NodeId(0), surface);
        denominator.insert_node(NodeId(1), surface);
        OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([
                Orientation::Undirected,
                Orientation::Default,
                Orientation::Reversed,
            ])),
            loop_energy_map: vec![energy.clone()],
            edge_energy_map: vec![LinearEnergyExpr::zero(), energy.clone(), energy],
            variants: vec![CFFVariant {
                origin: Some("native snapshot test".to_owned()),
                prefactor: Rational::from((-2, 3)),
                half_edges: vec![EdgeIndex(1), EdgeIndex(1)],
                denominator_edges: vec![EdgeIndex(1), EdgeIndex(1)],
                denominator_surface_signs: BTreeMap::from([(surface, -1)]),
                denominator_edge_support_signs: BTreeMap::from([(
                    vec![EdgeIndex(1), EdgeIndex(1)],
                    -1,
                )]),
                uniform_scale_power: 2,
                numerator_surfaces: vec![surface, surface],
                denominator,
            }],
        }
    }

    #[test]
    fn native_residue_key_identity_keeps_exact_affine_coefficients() {
        let native = native_residue();
        let key = PyResidueMapKey::from(&native);
        assert_eq!(key.canonical_string(), native.residue_map_key());
        assert_eq!(
            key.loop_energy_map[0].internal_terms[0].1.clone() * Rational::from(7),
            Rational::from(i64::MAX) * Rational::from(i64::MAX)
        );
        let mut relabeled = native.clone();
        relabeled.data.label = Some("another native ID's diagnostic label".to_owned());
        relabeled.data.numerator_map_index = Some(97);
        relabeled.variants[0].prefactor = Rational::from(99);
        let same_key = PyResidueMapKey::from(&relabeled);
        assert_eq!(key, same_key);
        let hash = |value: &PyResidueMapKey| {
            let mut hasher = DefaultHasher::new();
            value.hash(&mut hasher);
            hasher.finish()
        };
        assert_eq!(hash(&key), hash(&same_key));
        for part in 0..4 {
            let mut different = native.clone();
            match part {
                0 => different.loop_energy_map[0].uniform_scale_coeff += Rational::from(1),
                1 => different.edge_energy_map[1].constant += Rational::from(1),
                2 => different.loop_energy_map[0].external_terms[0].1 += Rational::from((1, 17)),
                _ => different.data.orientation[EdgeIndex(1)] = Orientation::Reversed,
            }
            assert_ne!(key, PyResidueMapKey::from(&different));
            assert_ne!(native.residue_map_key(), different.residue_map_key());
        }
    }

    #[test]
    fn native_residue_snapshot_keeps_duplicate_ids_and_repeated_denominators() {
        let first = native_residue();
        let mut second = first.clone();
        second.data.label = Some("second native ID".to_owned());
        let mut third = first.clone();
        third.edge_energy_map[1].uniform_scale_coeff += Rational::from(1);
        let mut generated = GeneratedThreeDExpression {
            representation: RepresentationMode::Ltd,
            expression: ThreeDExpression {
                orientations: vec![first, second, third].into(),
                surfaces: SurfaceCache::<Esurface, Hsurface>::new(),
                residual_denominators: Vec::new(),
            },
            energy_factor_ownership: CffEnergyFactorOwnership::VariantLocal,
            energy_factor_components: Vec::new(),
            source_energy_degree_bounds: Vec::new(),
            denominator_only_global_prefactor_sign: CffGlobalPrefactorSign::default(),
            core_global_prefactor_sign: CffGlobalPrefactorSign::default(),
        };
        let snapshot = PyResidueMap::new("test_graph".to_owned(), &generated);
        generated.expression.orientations.clear();
        let entries = snapshot.grouped_entries();
        assert_eq!(entries.len(), 2);
        assert_eq!(
            entries[0]
                .1
                .iter()
                .map(|residue| residue.native_id)
                .collect::<Vec<_>>(),
            [0, 1]
        );
        assert_eq!(entries[1].1[0].native_id, 2);
        let variant = &entries[0].1[0].inner.variants[0];
        assert_eq!(variant.half_edges, [EdgeIndex(1), EdgeIndex(1)]);
        assert_eq!(variant.denominator_edges, [EdgeIndex(1), EdgeIndex(1)]);
        let nodes = variant.denominator.iter_nodes().collect::<Vec<_>>();
        assert_eq!(nodes.len(), 3);
        assert_eq!(nodes[1].data, nodes[2].data);
        assert_eq!(nodes[1].children, [NodeId(2)]);
        assert_eq!(nodes[2].parent, Some(NodeId(1)));
    }
}
