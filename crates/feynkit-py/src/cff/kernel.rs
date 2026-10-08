//! One inspection model for topology-only and bounded-numerator CFF generation.
use super::*;
use crate::energy::integration::EnergyGraph;
use feynkit_cff::generalized::{
    CffEnergyFactorOwnership, CffGlobalPrefactorSign, Generate3DExpressionOptions, HybridSurfaceID,
    RepresentationMode, generate_3d_expression,
};
use linnet::half_edge::involution::{EdgeIndex, Orientation};

#[derive(Clone)]
pub(super) struct CffTerm {
    pub path: Vec<usize>,
    pub coefficient: Atom,
    pub numerator: Atom,
    pub energies: Vec<usize>,
    pub numerator_factors: Vec<usize>,
    pub energy_map: BTreeMap<usize, Atom>,
    pub origin: Option<String>,
}
impl CffTerm {
    pub fn prefactor(&self) -> Atom {
        self.energies
            .iter()
            .fold(self.coefficient.clone(), |value, edge| {
                value / (Atom::num(2) * feynkit_cff::symbols::on_shell_atom(EdgeIndex(*edge)))
            })
    }
    pub fn numerator_atom(&self, surfaces: &[PyEnergySurface]) -> Atom {
        self.numerator_factors
            .iter()
            .fold(self.numerator.clone(), |value, id| {
                value * surfaces[*id].atom(false)
            })
    }
    pub fn atom(&self, surfaces: &[PyEnergySurface]) -> Atom {
        self.path.iter().fold(
            self.prefactor() * self.numerator_atom(surfaces),
            |value, id| value / surfaces[*id].atom(false),
        )
    }
}
#[derive(Clone)]
pub(super) struct CffOrientationData {
    pub id: usize,
    pub directions: BTreeMap<usize, EdgeOrientation>,
    pub terms: Vec<CffTerm>,
    pub expression: Atom,
}
impl CffOrientationData {
    pub fn atom(&self) -> Atom {
        self.expression.clone()
    }
}
#[derive(Clone)]
pub(super) struct CffKernel {
    pub orientations: Vec<CffOrientationData>,
    pub surfaces: Vec<PyEnergySurface>,
    pub report: CffReport,
}
impl CffKernel {
    pub fn to_atom(&self) -> Atom {
        self.orientations
            .iter()
            .fold(Atom::Zero, |sum, o| sum + o.atom())
    }
    pub fn term_count(&self) -> usize {
        self.orientations.iter().map(|o| o.terms.len()).sum()
    }

    pub fn from_topology(result: CffResult) -> Self {
        let surfaces: Vec<_> = (0..result.surfaces.energy_surfaces().len())
            .map(|i| SurfaceId::Energy(feynkit_cff::EnergySurfaceId(i)))
            .chain(
                (0..result.surfaces.h_surfaces().len())
                    .map(|i| SurfaceId::H(feynkit_cff::HSurfaceId(i))),
            )
            .map(|id| PyEnergySurface::from_cff(id, &result.surfaces))
            .collect();
        let keys: BTreeMap<_, _> = surfaces
            .iter()
            .enumerate()
            .map(|(i, s)| (s.cff_id.unwrap(), i))
            .collect();
        let orientations = result
            .expression
            .orientations()
            .iter()
            .map(|o| CffOrientationData {
                id: o.id.index(),
                expression: o.expression.to_atom_inverse(),
                directions: o
                    .data
                    .orientation
                    .iter()
                    .map(|(e, d)| (e.index(), *d))
                    .collect(),
                terms: o
                    .denominator_products()
                    .into_iter()
                    .filter(|path| !path.contains(&SurfaceId::Infinite))
                    .map(|path| CffTerm {
                        path: path
                            .into_iter()
                            .filter(|s| *s != SurfaceId::Unit)
                            .map(|s| keys[&s])
                            .collect(),
                        coefficient: Atom::num(1),
                        numerator: Atom::num(1),
                        energies: Vec::new(),
                        numerator_factors: Vec::new(),
                        energy_map: BTreeMap::new(),
                        origin: None,
                    })
                    .collect(),
            })
            .collect();
        Self {
            orientations,
            surfaces,
            report: result.report,
        }
    }

    pub fn from_numerator(
        diagram: &feynkit_graph::FeynmanDiagram,
        numerator: &Atom,
    ) -> PyResult<(Self, Atom, BTreeMap<usize, usize>)> {
        let graph = EnergyGraph::new(diagram)?;
        let bounds = graph.numerator_bounds(numerator)?;
        let generated = generate_3d_expression(
            &graph.parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Some(bounds.clone()),
                ..Default::default()
            },
        )
        .map_err(|e| error::CffError::new_err(e.to_string()))?;
        if !generated.expression.residual_denominators.is_empty() {
            return Err(error::CffError::new_err(
                "residual four-dimensional denominators are not supported by the CFF Python constructor",
            ));
        }
        if generated
            .expression
            .surfaces
            .linear_surface_cache
            .iter()
            .any(|s| s.expression.uses_uniform_scale())
        {
            return Err(error::CffError::new_err(
                "uniform numerator sampling scales require an explicit scale",
            ));
        }
        // Native metadata owns each independent rational contour frame. Do not
        // guess a global parity from the final number of displayed factors.
        let normalization = Atom::num(generated.energy_factor_components.iter().fold(
            1,
            |sign, component| {
                let frame = match component.ownership {
                    CffEnergyFactorOwnership::GlobalSourceProduct => {
                        component.core_global_prefactor_sign
                    }
                    CffEnergyFactorOwnership::VariantLocal => {
                        component.denominator_only_global_prefactor_sign
                    }
                };
                sign * CffGlobalPrefactorSign::from_exponent(component.internal_edge_ids.len())
                    .product(frame)
                    .factor()
            },
        ));
        let surfaces: Vec<_> = generated
            .expression
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, s)| {
                let definition = graph.physical_energy(&s.expression);
                let boundary: BTreeSet<_> = definition
                    .internal_terms
                    .iter()
                    .filter(|(_, c)| !c.is_zero())
                    .map(|(e, _)| e.0)
                    .collect();
                let mut components = Vec::<BTreeSet<usize>>::new();
                let mut remaining: BTreeSet<_> = diagram.vertices().map(|(v, _)| v.0).collect();
                while let Some(&start) = remaining.first() {
                    let mut members = BTreeSet::from([start]);
                    let mut queue = vec![start];
                    remaining.remove(&start);
                    while let Some(v) = queue.pop() {
                        for (_, ends, _) in
                            diagram.edges().filter(|(e, _, _)| !boundary.contains(&e.0))
                        {
                            if let (Some(a), Some(b)) = (ends.source, ends.target) {
                                let next = if a.0 == v {
                                    Some(b.0)
                                } else if b.0 == v {
                                    Some(a.0)
                                } else {
                                    None
                                };
                                if let Some(next) = next
                                    && remaining.remove(&next)
                                {
                                    members.insert(next);
                                    queue.push(next);
                                }
                            }
                        }
                    }
                    components.push(members);
                }
                let mut vertices: Vec<usize> = components
                    .into_iter()
                    .find(|members| {
                        let crossing: BTreeSet<_> = diagram
                            .edges()
                            .filter_map(|(e, ends, _)| match (ends.source, ends.target) {
                                (Some(a), Some(b))
                                    if members.contains(&a.0) != members.contains(&b.0) =>
                                {
                                    Some(e.0)
                                }
                                _ => None,
                            })
                            .collect();
                        !boundary.is_empty() && crossing == boundary
                    })
                    .map(|v| v.into_iter().collect())
                    .unwrap_or_default();
                if s.origin != feynkit_cff::generalized::surface::SurfaceOrigin::Physical
                    || s.numerator_only
                {
                    vertices.clear();
                }
                PyEnergySurface::from_linear(id, s.kind, &definition)
                    .with_vertices(vertices)
                    .with_provenance(s.origin, s.numerator_only)
            })
            .collect();
        let mut orientations = Vec::<CffOrientationData>::new();
        for o in &generated.expression.orientations {
            let directions: BTreeMap<_, _> = graph
                .indices
                .internal
                .iter()
                .map(|(local, physical)| {
                    (
                        *physical,
                        match o.data.orientation[EdgeIndex(*local)] {
                            Orientation::Default => EdgeOrientation::Default,
                            Orientation::Reversed => EdgeOrientation::Reversed,
                            Orientation::Undirected => EdgeOrientation::Undirected,
                        },
                    )
                })
                .collect();
            let position = orientations
                .iter()
                .position(|other| other.directions == directions)
                .unwrap_or_else(|| {
                    let id = orientations.len();
                    orientations.push(CffOrientationData {
                        id,
                        directions,
                        terms: Vec::new(),
                        expression: Atom::Zero,
                    });
                    id
                });
            let value = graph.numerator_at(numerator, &o.edge_energy_map);
            let energy_map = graph
                .indices
                .internal
                .iter()
                .map(|(local, physical)| {
                    (
                        *physical,
                        graph
                            .physical_energy(&o.edge_energy_map[*local])
                            .to_atom(&[]),
                    )
                })
                .collect::<BTreeMap<_, _>>();
            for variant in &o.variants {
                if variant.uniform_scale_power != 0 {
                    return Err(error::CffError::new_err(
                        "uniform numerator sampling scales require an explicit scale",
                    ));
                }
                let ids = |list: Vec<HybridSurfaceID>| -> PyResult<Vec<usize>> {
                    list.into_iter()
                        .filter(|s| *s != HybridSurfaceID::Unit)
                        .map(|s| match s {
                            HybridSurfaceID::Linear(id) => Ok(id.0),
                            _ => Err(error::CffError::new_err(
                                "expected affine generated CFF surface",
                            )),
                        })
                        .collect()
                };
                let weight = CffTerm {
                    path: Vec::new(),
                    coefficient: Atom::num(variant.prefactor.clone()),
                    numerator: value.clone(),
                    energies: variant
                        .half_edges
                        .iter()
                        .map(|e| graph.indices.internal[&e.0])
                        .collect(),
                    numerator_factors: ids(variant.numerator_surfaces.clone())?,
                    energy_map: BTreeMap::new(),
                    origin: None,
                };
                orientations[position].expression +=
                    weight.atom(&surfaces) * variant.denominator.to_atom_inv();
                for leaf in variant.denominator.get_bottom_layer() {
                    let mut node = variant.denominator.get_node(leaf);
                    let mut path = Vec::new();
                    loop {
                        path.push(node.data);
                        let Some(parent) = node.parent else { break };
                        node = variant.denominator.get_node(parent);
                    }
                    if path.contains(&HybridSurfaceID::Infinite) {
                        continue;
                    }
                    path.reverse();
                    orientations[position].terms.push(CffTerm {
                        path: ids(path)?,
                        coefficient: Atom::num(variant.prefactor.clone()),
                        numerator: value.clone(),
                        energies: variant
                            .half_edges
                            .iter()
                            .map(|e| graph.indices.internal[&e.0])
                            .collect(),
                        numerator_factors: ids(variant.numerator_surfaces.clone())?,
                        energy_map: energy_map.clone(),
                        origin: variant.origin.clone(),
                    });
                }
            }
        }
        let report = CffReport {
            candidate_orientations: orientations.len(),
            acyclic_orientations: orientations.len(),
            unfolded_terms: orientations.iter().map(|o| o.terms.len()).sum(),
            interned_surfaces: surfaces.len(),
        };
        Ok((
            Self {
                orientations,
                surfaces,
                report,
            },
            normalization,
            bounds
                .into_iter()
                .map(|(i, d)| (graph.indices.internal[&i], d))
                .collect(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_graph::FeynmanDiagram;
    use feynkit_model::Model;
    use symbolica::symbol;

    #[test]
    fn polynomial_cff_matches_exact_one_loop_contours() {
        Python::initialize();
        let model = Arc::new(
            Model::from_json(include_str!("../../tests/fixtures/scalars_2p_3p.json")).unwrap(),
        );
        let diagram = FeynmanDiagram::from_dot(model, "digraph triangle { edge [particle=scalar_0]; ext [style=invis]; ext -> a; b -> ext; c -> ext; a -> b; b -> c; c -> a [lmb_id=0]; }").unwrap();
        let graph = EnergyGraph::new(&diagram).unwrap();
        let q = feynkit_cff::symbols::external_energy_atom(EdgeIndex(graph.indices.internal[&0]));
        let t = Atom::var(symbol!("cff_test::t"));
        for degree in 0..=4 {
            let numerator = (q.clone() + Atom::num((2, 3))).pow(degree);
            let (kernel, normalization, _) =
                CffKernel::from_numerator(&diagram, &numerator).unwrap();
            let actual = (normalization * kernel.to_atom()).replace_multiple(
                kernel
                    .surfaces
                    .iter()
                    .map(|s| symbolica::id::Replacement::new(s.atom(false), s.atom(true))),
            );
            let momenta: Vec<_> = graph
                .parsed
                .internal_edges
                .iter()
                .map(|edge| {
                    let sign = edge.signature.loop_signature[0];
                    let shift = edge.signature.external_signature.iter().enumerate().fold(
                        Atom::Zero,
                        |sum, (i, c)| {
                            sum + Atom::num(*c)
                                * feynkit_cff::symbols::external_energy_atom(EdgeIndex(
                                    graph.indices.external[&i],
                                ))
                        },
                    );
                    (sign, shift)
                })
                .collect();
            let routed_numerator = numerator.replace_multiple(graph.indices.internal.iter().map(
                |(local, physical)| {
                    symbolica::id::Replacement::new(
                        feynkit_cff::symbols::external_energy_atom(EdgeIndex(*physical)),
                        Atom::num(momenta[*local].0) * &t + &momenta[*local].1,
                    )
                },
            ));
            let mut expected = Atom::Zero;
            for (local, (sign, shift)) in momenta.iter().enumerate() {
                let energy =
                    feynkit_cff::symbols::on_shell_atom(EdgeIndex(graph.indices.internal[&local]));
                let root = &energy - Atom::num(*sign) * shift;
                let mut term = -routed_numerator
                    .replace_multiple([symbolica::id::Replacement::new(t.clone(), root.clone())])
                    / (Atom::num(2) * energy);
                for (other, (s, p)) in momenta.iter().enumerate() {
                    if other != local {
                        term /= (Atom::num(*s) * &root + p).pow(2)
                            - feynkit_cff::symbols::on_shell_atom(EdgeIndex(
                                graph.indices.internal[&other],
                            ))
                            .pow(2);
                    }
                }
                expected += term;
            }
            assert!((actual - expected).together().is_zero(), "degree {degree}");
        }
    }
}
