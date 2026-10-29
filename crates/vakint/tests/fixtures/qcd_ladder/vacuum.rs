//! Mechanical vacuum incidence plus a native verified simultaneous routing.

use super::*;
use feynkit_graph::IntegralMapping;

pub(super) struct Line {
    pub momentum: Atom,
    pub original_edges: Vec<EdgeId>,
}

pub(super) struct Vacuum {
    pub diagram: FeynmanDiagram,
    pub family: IntegralFamily,
    pub mapping: IntegralMapping,
    pub lines: Vec<Line>,
    pub mass: Atom,
}

impl Vacuum {
    pub fn new(source: &FeynmanDiagram, dimension: &Atom) -> Result<Self> {
        let loop_count = source.loop_count();
        let rows = loop_rows(source);
        let attachments = source
            .edges()
            .filter_map(|(_, ends, edge)| {
                edge.external
                    .as_ref()
                    .and_then(|_| ends.source.or(ends.target))
            })
            .collect::<BTreeSet<_>>();
        if attachments.len() != 2 {
            return Err("the source needs two distinct external attachment vertices".into());
        }
        let mut removed = BTreeSet::new();
        let mut records = Vec::new();
        // With no rung, both current insertions belong to the same two-line
        // bubble. Collapse them together to one dotted self-loop; independently
        // collapsing both insertions would remove every vacuum vertex.
        let retained_attachment = if loop_count == 1 {
            let internal = source
                .edges()
                .filter(|(_, ends, _)| !is_external_edge(source, *ends))
                .collect::<Vec<_>>();
            let [left, right] = internal.as_slice() else {
                return Err("the zero-rung source is not a two-propagator bubble".into());
            };
            let p = &rows[&left.0];
            let q = &rows[&right.0];
            if !(p - q).expand().is_zero() && !(p + q).expand().is_zero() {
                return Err("the zero-rung propagators differ at zero external momentum".into());
            }
            let vertex = *attachments.first().expect("two attachments checked");
            records.push((
                vertex,
                vertex,
                Line {
                    momentum: p.clone(),
                    original_edges: vec![left.0, right.0],
                },
            ));
            removed.extend([left.0, right.0]);
            Some(vertex)
        } else {
            None
        };
        for vertex in attachments.iter().filter(|_| loop_count > 1) {
            let adjacent = source
                .edges()
                .filter(|(_, ends, _)| {
                    !is_external_edge(source, *ends)
                        && (ends.source == Some(*vertex) || ends.target == Some(*vertex))
                })
                .collect::<Vec<_>>();
            let [left, right] = adjacent.as_slice() else {
                return Err("amputation must leave a degree-two insertion vertex".into());
            };
            if adjacent.iter().any(|(_, _, edge)| {
                source
                    .model()
                    .particle_by_id(edge.particle)
                    .unwrap()
                    .pdg_code
                    .abs()
                    != 1
            }) {
                return Err("only the two fermion propagators adjacent to each external gamma insertion may be merged".into());
            }
            let p = &rows[&left.0];
            let q = &rows[&right.0];
            if !(p - q).expand().is_zero() && !(p + q).expand().is_zero() {
                return Err("the two merged massless propagators are not equal at q=0".into());
            }
            let other = |ends: &feynkit_graph::EdgeEndpoints| {
                if ends.source == Some(*vertex) {
                    ends.target.expect("internal edge")
                } else {
                    ends.source.expect("internal edge")
                }
            };
            records.push((
                other(&left.1),
                other(&right.1),
                Line {
                    momentum: if left.1.target == Some(*vertex) {
                        p.clone()
                    } else {
                        -p
                    },
                    original_edges: vec![left.0, right.0],
                },
            ));
            removed.extend([left.0, right.0]);
        }
        for (id, ends, _) in source.edges() {
            if !is_external_edge(source, ends) && !removed.contains(&id) {
                records.push((
                    ends.source.expect("internal edge"),
                    ends.target.expect("internal edge"),
                    Line {
                        momentum: rows[&id].clone(),
                        original_edges: vec![id],
                    },
                ));
            }
        }
        let expected_incidence_edges = if loop_count == 1 {
            1
        } else {
            3 * loop_count - 3
        };
        let expected_insertion_dots = if loop_count == 1 { 1 } else { 2 };
        if records.len() != expected_incidence_edges
            || records
                .iter()
                .filter(|r| r.2.original_edges.len() == 2)
                .count()
                != expected_insertion_dots
            || records
                .iter()
                .map(|r| r.2.original_edges.len())
                .sum::<usize>()
                != 3 * loop_count - 1
        {
            return Err(
                "vacuum topology failed its loop-dependent incidence and total-power census".into(),
            );
        }
        let model = Model::from_json(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/phi3.json"
        )))
        .expect("valid phi3 fixture");
        let particle = model.particle_id("phi")?;
        let mass = model.particle_by_id(particle)?.symbolic_mass(&model);
        let mut builder = FeynmanDiagram::builder(model, "projected_qcd_ladder_vacuum");
        let vertices = source
            .vertices()
            .filter(|(id, _)| !attachments.contains(id) || Some(*id) == retained_attachment)
            .map(|(id, _)| {
                (
                    id,
                    builder.add_vertex(DiagramVertex {
                        name: format!("v{}", id.0),
                        interaction: None,
                        numerator: SymbolicTensor::infer(Atom::one()).expect("one is scalar"),
                    }),
                )
            })
            .collect::<BTreeMap<_, _>>();
        for (index, (left, right, _)) in records.iter().enumerate() {
            let id = builder.add_edge(
                vertices[left],
                vertices[right],
                DiagramEdge::new(particle, false),
            )?;
            if id != EdgeId(index) {
                return Err("native builder changed fixture edge insertion identities".into());
            }
        }
        let diagram = builder.build()?;
        if diagram.loop_count() != loop_count {
            return Err("amputation/collapse changed physical loop count".into());
        }
        let loops = (0..loop_count)
            .map(|axis| symbols::loop_momentum().call(axis))
            .collect::<Vec<_>>();
        let kin = Kinematics::in_dimension(dimension)?.with_momenta(loops.clone())?;
        let source_family = IntegralFamily::new(
            loops.clone(),
            vec![],
            records
                .iter()
                .map(|(_, _, line)| {
                    kin.scalar_product(&line.momentum, &line.momentum)
                        .map(|square| square - mass.clone().pow(2))
                })
                .collect::<std::result::Result<Vec<_>, _>>()?,
            &kin,
        )?;
        let family = diagram.propagator_family(&kin)?;
        let target_rows = loop_rows(&diagram);
        // Each original native chord survives, possibly as its equal partner
        // at an insertion. Reading those unit rows determines one witness;
        // the native family API verifies every denominator and its Jacobian.
        let images = loops
            .iter()
            .map(|axis| {
                records
                    .iter()
                    .enumerate()
                    .find_map(|(slot, (_, _, line))| {
                        if (&line.momentum - axis).expand().is_zero() {
                            Some(target_rows[&EdgeId(slot)].clone())
                        } else if (&line.momentum + axis).expand().is_zero() {
                            Some(-&target_rows[&EdgeId(slot)])
                        } else {
                            None
                        }
                    })
                    .ok_or_else(|| {
                        "an original chord disappeared during exact equal-line collapse".into()
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        let mapping = source_family
            .mapping_to(&family, &images)?
            .ok_or("native all-row/unit-measure vacuum transport verification failed")?;
        if mapping.denominator_map() != (0..records.len()).collect::<Vec<_>>() {
            return Err(
                "the verified map did not preserve the declared physical line identities".into(),
            );
        }
        Ok(Self {
            diagram,
            family,
            mapping,
            lines: records.into_iter().map(|r| r.2).collect(),
            mass,
        })
    }

    pub fn original_slot(&self, edge: EdgeId) -> Result<usize> {
        self.lines
            .iter()
            .position(|line| line.original_edges.contains(&edge))
            .ok_or_else(|| "original internal edge has no vacuum line".into())
    }

    /// Count exact quadratic forms independently of graph-edge incidence.
    /// Native ingress must combine every coincident fermion pair.
    pub fn validate_initial_ingress(&self, evidence: &Path) -> Result<()> {
        let initial = self
            .lines
            .iter()
            .map(|line| line.original_edges.len() as i32)
            .collect::<Vec<_>>();
        let mut forms = BTreeMap::<Atom, i32>::new();
        for (form, power) in self.family.denominators().iter().zip(&initial) {
            *forms.entry(form.expand()).or_default() += power;
        }
        let mut expected = forms.values().copied().collect::<Vec<_>>();
        expected.sort();
        let loop_count = self.diagram.loop_count();
        let expected_powers = std::iter::repeat_n(1, loop_count - 1)
            .chain(std::iter::repeat_n(2, loop_count))
            .collect::<Vec<_>>();
        if forms.len() != 2 * loop_count - 1 || expected != expected_powers {
            return Err(
                "native quadratic-form census does not match gluon rungs and doubled fermion lines"
                    .into(),
            );
        }
        let marker = Atom::var(symbol!("qcd_ladder_control::retained_numerator"));
        let input = self.input(&marker, &initial)?;
        let expression = VakintExpression::try_from(input.as_view())?;
        if expression.0.len() != 1 || expression.0[0].numerator != marker {
            return Err(
                "native repeated-denominator ingress changed the retained scalar numerator".into(),
            );
        }
        let [id, left, right, momentum, mass, power] =
            ["id_", "left_", "right_", "momentum_", "mass_", "power_"]
                .map(|name| symbol!(&format!("qcd_ladder_control::{name}")));
        let s = &S;
        let pattern = s.prop.call((
            Atom::var(id),
            s.edge.call((Atom::var(left), Atom::var(right))),
            Atom::var(momentum),
            Atom::var(mass),
            Atom::var(power),
        ));
        let records = expression.0[0]
            .integral
            .pattern_match(&pattern.to_pattern(), None, None)
            .collect::<Vec<_>>();
        let mut graph = Graph::<(), ()>::new();
        let mut nodes = BTreeMap::new();
        let mut powers = Vec::new();
        for record in &records {
            let source = *nodes
                .entry(record[&left].clone())
                .or_insert_with(|| graph.add_node(()));
            let target = *nodes
                .entry(record[&right].clone())
                .or_insert_with(|| graph.add_node(()));
            graph.add_edge(source, target, false, ())?;
            powers.push(i32::try_from(i64::try_from(record[&power].as_view())?)?);
        }
        powers.sort();
        if records.len() != 2 * loop_count - 1
            || powers != expected_powers
            || !graph.is_connected()
            || graph.num_loops() != loop_count
        {
            return Err(
                "native repeated-denominator coalescing changed the physical scalar graph".into(),
            );
        }
        // Even an exactly cancellable numerator factor stays in the numerator
        // at this boundary. Only equal denominator forms are coalesced.
        let raw_momentum = &self.lines[0].momentum;
        let raw_square = self
            .family
            .kinematics()
            .scalar_product(raw_momentum, raw_momentum)?;
        let cancellable = &marker * (raw_square - self.mass.clone().pow(2));
        let retained = VakintExpression::try_from(self.input(&cancellable, &initial)?.as_view())?;
        let first = records
            .iter()
            .find(|record| record[&id].is_zero())
            .ok_or("initial ingress lost its first physical denominator")?;
        let expected =
            &marker * (S.dot.call((&first[&momentum], &first[&momentum])) - &first[&mass]);
        if retained.0.len() != 1
            || retained.0[0].integral != expression.0[0].integral
            || !(&retained.0[0].numerator - expected)
                .expand()
                .together()
                .is_zero()
        {
            return Err(
                "benchmark ingress canceled a physical denominator against its numerator".into(),
            );
        }
        std::fs::write(
            evidence.join("native-vacuum-counts.txt"),
            format!(
                "Original: {} internal edges, {} vertices, {loop_count} loops.\nAfter external insertion collapse: {} incidence edges, {} vertices, A={}.\nExact native quadratic forms: {} distinct ({} fermion forms of power2, {} gluon forms of power1).\nNative ingress: {} propagators, {loop_count} loops, numerator retained exactly.\n",
                3 * loop_count - 1,
                2 * loop_count,
                self.lines.len(),
                self.diagram.vertices().count(),
                3 * loop_count - 1,
                forms.len(),
                loop_count,
                loop_count - 1,
                records.len(),
            ),
        )?;
        Ok(())
    }

    pub fn input(&self, numerator: &Atom, powers: &[i32]) -> Result<Atom> {
        let numerator = self.mapping.apply(numerator);
        let muv = Atom::var(symbol!("vakint::mUV"));
        // Native ingress authenticates the already-substituted family against
        // the original graph plus these simultaneous parameter replacements.
        let family = IntegralFamily::new(
            self.family.loop_momenta().to_vec(),
            vec![],
            self.family
                .denominators()
                .iter()
                .map(|denominator| denominator.replace(self.mass.clone()).with(muv.clone()))
                .collect(),
            self.family.kinematics(),
        )?;
        let powers = self
            .mapping
            .map_powers(powers)?
            .into_iter()
            .map(Integer::from)
            .collect();
        Ok(VakintExpression::from_diagram(
            &self.diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions {
                powers: Some(powers),
                parameter_substitutions: vec![(self.mass.clone(), muv)],
                ..Default::default()
            },
        )?
        .into())
    }
}

#[test]
fn ingress_retains_cancelable_numerator_and_physical_powers() -> Result<()> {
    let generated = dressed_diagram()?;
    let vacuum = Vacuum::new(&generated.diagrams[0], &dimension()?)?;
    let evidence = tempfile::tempdir()?;
    vacuum.validate_initial_ingress(evidence.path())
}
