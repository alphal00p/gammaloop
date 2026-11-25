//! Reconstruct the renormalized massless self-energy in a basis of scalar
//! bubble chains. Each coefficient stays inside the pole map of its graph.
use super::*;
use linnet::half_edge::subgraph::SubSetLike;

static BUBBLE_POWER: LazyLock<Atom> =
    LazyLock::new(|| Atom::var(symbol!("renormalized_insertion_bubble_power")));

static RADIAL_X: LazyLock<Symbol> = LazyLock::new(|| symbol!("factored_insertion_radial_x"));
static RADIAL_Q2: LazyLock<Symbol> = LazyLock::new(|| symbol!("factored_insertion_radial_q2"));

type ReplacementTerm = (
    FeynmanDiagram,
    Vec<DiagramRelabeling>,
    Vec<(VertexId, usize)>,
);

impl PropagatorInsertionEvaluator {
    fn insertion_coefficients(
        &self,
        insertion: &FeynmanDiagram,
        vakint: &Vakint,
        depth: usize,
    ) -> Result<Vec<Atom>> {
        let stored_depth = if crate::uv::delta_store::installed().is_some() {
            depth.max(6)
        } else {
            depth
        };
        let mut values = self.insertion_coefficients_at(insertion, vakint, stored_depth)?;
        if stored_depth > depth {
            for (j, value) in values.iter_mut().enumerate() {
                *value = Self::truncate_insertion(value, depth as i64 - insertion.loop_count() as i64 + j as i64 - 1)?;
            }
        }
        Ok(values)
    }

    pub(super) fn renormalized_key(&self, depth: usize) -> String {
        format!("renormalized-sigma:attached-v1:depth={depth}:closure={:?}:{}", self.closure, self.settings_key)
    }

    pub(super) fn insertion_coefficients_at(
        &self, insertion: &FeynmanDiagram, vakint: &Vakint, depth: usize,
    ) -> Result<Vec<Atom>> {
        let n = insertion.loop_count();
        let compute = || -> Result<Vec<Atom>> {
            if n == 1 {
                return Ok(vec![
                    Atom::Zero,
                    self.normalized_insertion(insertion, vakint, depth, insertion.name())?,
                ]);
            }
            self.renormalized_insertion_uncached(insertion, vakint, depth)
        };
        let Some(cache) = &self.counterterm_cache else {
            return compute();
        };
        let graph = insertion.to_gamma_loop_graph(None, true)?;
        let boundary = graph.external_filter::<SuBitGraph>();
        let support = graph.no_dummy().subtract(&boundary);
        // One tensor-valued polynomial is one atomic cache row. Its degree
        // labels the chain length and the key fixes the exact epsilon depth.
        let packed = cache.counterterm_with_kind(
            &graph,
            &support,
            insertion.model(),
            self.renormalized_key(depth),
            crate::uv::counterterm_cache::KeyKind::RenormalizedInsertion,
            || {
                Ok(compute()?
                    .into_iter()
                    .enumerate()
                    .map(|(j, c)| c * BUBBLE_POWER.pow(j))
                    .sum())
            },
        )?;
        let (_, routing) = RStarIntegrand::physical_external_routing(&graph, &boundary)?;
        let mut values = vec![Atom::Zero; n + 1];
        for (monomial, coefficient) in packed.coefficient_list::<i8>(&[&*BUBBLE_POWER]) {
            let j = (0..=n)
                .find(|j| monomial == BUBBLE_POWER.pow(*j))
                .ok_or_else(|| eyre!("invalid renormalized insertion degree"))?;
            values[j] =
                crate::uv::rstar::linearize::linearize(coefficient.replace_multiple(&routing));
        }
        Ok(values)
    }

    /// Sum projected counterterms of the genuine replacement graphs. The
    /// caller supplies the graph-specific external projection. No physical
    /// preparation of the original parent is evaluated.
    pub fn renormalized_counterterm(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
        closing: bool,
        project: impl FnMut(&Graph, &Model, &RStarPlan) -> Result<Atom>,
    ) -> Result<Option<Atom>> {
        self.project_insertions(diagram, vakint, closing, false, project)
    }

    /// Evaluate the physical skeleton with a single retained R(σ) tensor at each vertex.
    /// Experimental: the projector can reject non-polynomial dressed vertices.
    pub fn factored_renormalized_counterterm(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
        closing: bool,
        project: impl FnMut(&Graph, &Model, &RStarPlan) -> Result<Atom>,
    ) -> Result<Option<Atom>> {
        self.project_insertions(diagram, vakint, closing, true, project)
    }

    fn project_insertions(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
        closing: bool,
        factored: bool,
        mut project: impl FnMut(&Graph, &Model, &RStarPlan) -> Result<Atom>,
    ) -> Result<Option<Atom>> {
        ensure!(
            diagram.cuts().is_empty(),
            "insertion extraction requires an uncut amplitude"
        );
        let insertions = self.selected_r_insertions(diagram)?;
        if insertions.is_empty() {
            return Ok(None);
        }
        let terms = self.replacement_terms(diagram, vakint, &insertions, factored)?;
        if crate::uv::rstar::plan::uv_children_log() || insertion_identity_log() {
            eprintln!(
                "RENORMALIZED_INSERTIONS graph={} insertions={} insertion_loops={} replacement_graphs={}",
                diagram.name(),
                insertions.len(),
                insertions.iter().map(|s| s.loops).sum::<usize>(),
                terms.len()
            );
        }
        let mut sum = Atom::Zero;
        for (term, maps, holes) in terms {
            // A zero R(sigma) annihilates this entire skeleton, including its
            // forest terms. Its zero tensor has no external slots to project.
            if holes.iter().any(|(id, _)| {
                term.vertex(*id)
                    .is_some_and(|v| v.numerator.expression().is_zero())
            }) {
                continue;
            }
            let graph = term.to_gamma_loop_graph(None, true)?;
            let mut plan = RStarPlan::with_automatic_edge(&graph, term.model(), &self.settings)?;
            if let Some(cache) = &self.counterterm_cache {
                plan = plan.with_counterterm_cache(cache.clone());
            }
            let provenance = crate::uv::rstar::uncontracted::k5::provenance(
                diagram,
                holes
                    .iter()
                    .zip(&insertions)
                    .map(|((vertex, length), sigma)| (vertex.0, *length, &sigma.support)),
            )?;
            let mut tagged = crate::uv::rstar::uncontracted::k5::scoped(
                provenance.as_ref(),
                &graph.name,
                closing && plan.closes_root(&graph),
                &self.settings,
                || {
                    crate::uv::rstar::expression::in_dressed_scope(|| {
                        project(&graph, term.model(), &plan)
                    })
                },
            )?;
            for map in maps.iter().rev() {
                tagged = map.restore_atom(&tagged)?;
            }
            sum += tagged;
        }
        Ok(Some(sum))
    }
    /// The physical skeleton with one unexpanded R(σ) tensor per hole.
    /// Its integration requires momentum-dependent analytic propagator powers.
    pub fn factored_renormalized_skeleton(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
    ) -> Result<Option<FeynmanDiagram>> {
        ensure!(
            diagram.cuts().is_empty(),
            "factored insertions require an uncut graph"
        );
        let holes = self.selected_r_insertions(diagram)?;
        if holes.is_empty() {
            return Ok(None);
        }
        let mut terms = self.replacement_terms(diagram, vakint, &holes, true)?;
        ensure!(terms.len() == 1, "expected one factored skeleton");
        Ok(terms.pop().map(|(diagram, _, _)| diagram))
    }

    /// Direct UV pole of a massless logarithmic scalar bubble dressed by
    /// two-point insertions. The physical external momentum keeps its IR finite.
    /// Each nontrivial scalar vertex must cancel one boundary propagator.
    /// Return None when the one-loop radial kernel does not apply.
    pub fn factored_scalar_bubble_counterterm(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
    ) -> Result<Option<Atom>> {
        use spenso::{
            network::tags::SPENSO_TAG,
            structure::representation::{Minkowski, RepName},
        };
        let holes = self.selected_r_insertions(diagram)?;
        let degree: usize = holes.iter().map(|s| s.loops).sum();
        if degree == 0 || diagram.loop_count() != degree + 1 {
            return Ok(None);
        }
        let Some(skeleton) = self.factored_renormalized_skeleton(diagram, vakint)? else {
            return Ok(None);
        };
        let support = skeleton.internal_subgraph();
        let external = skeleton
            .edges()
            .filter(|(_, _, e)| e.external.is_some())
            .filter_map(|(_, ends, _)| ends.source.or(ends.target))
            .collect::<Vec<_>>();
        if !matches!(external.as_slice(), [a,b] if a != b)
            || !skeleton.underlying().is_connected(&support)
            || !skeleton.underlying().bridges_of(&support).is_empty()
        {
            return Ok(None);
        }
        let x = Atom::var(*RADIAL_X);
        let t = Atom::var(*RADIAL_Q2);
        let eps = Atom::var(GS.dim_epsilon);
        let mut momenta = Vec::new();
        for (id, _, edge) in skeleton.edges().filter(|(_, _, e)| e.external.is_none()) {
            if edge.is_dummy
                || edge.propagator_power != 1
                || !skeleton
                    .model()
                    .particle_by_id(edge.particle)?
                    .is_massless(skeleton.model())
                || !NumeratorInterface::export(&edge.numerator).0.is_one()
            {
                return Ok(None);
            }
            let q = function!(
                GS.emr_mom,
                id.0,
                Minkowski {}.new_rep(GS.dim).to_symbolic([])
            );
            let q2 = function!(SPENSO_TAG.dot, &q, &q);
            momenta.push(((-&q2 / Atom::var(GS.mu_r_sq)).pow(-&eps), q2));
        }
        let mut powers = 0;
        let mut factors = Vec::new();
        for (_, v) in skeleton.vertices() {
            let mut local = NumeratorInterface::export(&v.numerator).0;
            for (power, _) in &momenta {
                local = local.replace(power.to_pattern()).with(x.to_pattern());
            }
            for (_, q2) in &momenta {
                local = local.replace(q2.to_pattern()).with(t.to_pattern());
            }
            if local.contains_symbol(GS.loop_mom)
                || local.contains_symbol(GS.external_mom)
                || local.contains_symbol(GS.emr_mom)
                || local.contains_symbol(GS.mu_r_sq)
                || super::super::regularity::has_tensor_structure(&local)
            {
                return Ok(None);
            }
            if local.contains_symbol(*RADIAL_Q2) {
                let at_one = local.replace(t.to_pattern()).with(Atom::one());
                if !(&local - &t * &at_one).expand().is_zero() {
                    return Ok(None);
                }
                local = at_one;
                powers += 1;
            }
            factors.push(
                crate::uv::approx::integrated::Integrated::regulated_for_expansion(
                    &local,
                    "factored_radial_vertex",
                )?
                .into_owned(),
            );
        }
        if powers + 2 != momenta.len() {
            return Ok(None);
        }
        // UV radial integration is (1/ε) ∫_0^1 dx. Integration of the
        // product at rational interpolation nodes is exact through 'degree'.
        let mut mean = Atom::Zero;
        for i in 0..=degree {
            let node = Atom::num((i as i64, degree as i64));
            let mut basis = vec![Atom::one()];
            for j in (0..=degree).filter(|j| *j != i) {
                let inv = Atom::one() / Atom::num(i as i64 - j as i64);
                let mut next = vec![Atom::Zero; basis.len() + 1];
                for (k, c) in basis.iter().enumerate() {
                    next[k] -= c * &inv * Atom::num(j as i64);
                    next[k + 1] += c * &inv * Atom::num(degree as i64);
                }
                basis = next;
            }
            let weight: Atom = basis
                .into_iter()
                .enumerate()
                .map(|(k, c)| c / Atom::num(k as i64 + 1))
                .sum();
            let value: Atom = factors
                .iter()
                .map(|f| f.replace(x.to_pattern()).with(node.to_pattern()))
                .product();
            mean += weight * value;
        }
        // Surface area relative to the unit G(1,1) integral. The first
        // missing term is O(ε^(degree+1)) and cannot contribute a pole.
        let unit = self.unit_scalar_bubble(diagram, vakint, degree + 2)?;
        let log_ratio: Atom = (2..=degree)
            .map(|n| {
                let coeff =
                    Atom::num(2).pow(n) - Atom::num(3) - Atom::num(if n % 2 == 0 { 1 } else { -1 });
                coeff * function!(symbolica::transcendental::zeta(), n) * eps.pow(n)
                    / Atom::num(n as i64)
            })
            .sum();
        let ratio = function!(Symbol::EXP, log_ratio)
            .series(GS.dim_epsilon, 0, degree as i64)?
            .to_atom()
            * (Atom::one() - Atom::num(2) * &eps)
            / (Atom::one() - &eps);
        Ok(Some(
            -(mean * unit * ratio)
                .series(GS.dim_epsilon, 0, -1)?
                .to_atom()
                .expand(),
        ))
    }

    fn replacement_terms(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
        insertions: &[PropagatorInsertion],
        factored: bool,
    ) -> Result<Vec<ReplacementTerm>> {
        let mut terms = vec![(
            diagram.clone(),
            Vec::<DiagramRelabeling>::new(),
            Vec::<(VertexId, usize)>::new(),
        )];
        for insertion in insertions {
            let excised = diagram.excise_with_mapping(&insertion.support, &[])?;
            let coefficients =
                self.insertion_coefficients(&excised.diagram, vakint, diagram.loop_count())?;
            let coefficients = if factored {
                vec![self.factored_insertion(
                    &excised.diagram,
                    vakint,
                    diagram.loop_count(),
                    &coefficients,
                )?]
            } else {
                coefficients
            };
            let mut next_terms = Vec::new();
            for (base, maps, holes) in terms {
                let mut current = insertion.clone();
                for map in &maps {
                    current.support = map.translate_subgraph(&current.support)?;
                    current.boundary = [
                        map.translate_half_edge(current.boundary[0])?,
                        map.translate_half_edge(current.boundary[1])?,
                    ];
                }
                let [near, far] = current.boundary;
                let particle = base
                    .edges()
                    .find(|(id, _, _)| *id == near.edge)
                    .ok_or_else(|| eyre!("missing insertion boundary"))?
                    .2
                    .particle;
                for (length, c) in coefficients
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| factored || !c.is_zero())
                {
                    let mut c = excised.relabeling.restore_atom(c)?;
                    for map in &maps {
                        c = map.translate_atom(&c, None)?;
                    }
                    let sign = if near.endpoint == far.endpoint { -1 } else { 1 };
                    c = c
                        .replace(function!(GS.emr_mom, far.edge.0, W_.x___))
                        .with(Atom::num(sign) * function!(GS.emr_mom, near.edge.0, W_.x___));
                    let c = open_slashes(&freshen_insertion_indices(
                        &c,
                        &[
                            NumeratorInterface::export(base.numerator()).0,
                            NumeratorInterface::export(base.projector()).0,
                            NumeratorInterface::export(base.numerator_prefactor()).0,
                        ],
                    )?);
                    let replacement =
                        base.scalar_bubble_surrogate_with_length(&current, particle, length)?;
                    let graph = replacement.with_coefficient(&SymbolicTensor::infer(c)?)?;
                    let mut term_maps = maps.clone();
                    let mut term_holes = holes
                        .iter()
                        .map(|(v, n)| {
                            Ok((
                                replacement
                                    .relabeling
                                    .vertex(*v)
                                    .ok_or_else(|| eyre!("lost an insertion vertex"))?,
                                *n,
                            ))
                        })
                        .collect::<Result<Vec<_>>>()?;
                    term_holes.push((replacement.coefficient_vertex, length));
                    term_maps.push(replacement.relabeling);
                    next_terms.push((graph, term_maps, term_holes));
                }
            }
            terms = next_terms;
            if terms.is_empty() {
                break;
            }
        }
        Ok(terms)
    }

    fn unit_scalar_bubble(
        &self,
        insertion: &FeynmanDiagram,
        vakint: &Vakint,
        depth: usize,
    ) -> Result<Atom> {
        use spenso::{
            network::tags::SPENSO_TAG,
            structure::representation::{Minkowski, RepName},
        };
        let q0 = function!(GS.emr_mom, 0, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let unit = self
            .integrate_cached(&scalar_bubble(insertion)?, vakint, depth)?
            .replace(function!(SPENSO_TAG.dot, &q0, &q0))
            .with(Atom::num(-1))
            .replace(Atom::var(GS.mu_r_sq))
            .with(Atom::one());
        ensure!(
            !unit.contains_symbol(GS.emr_mom),
            "unit-scale scalar bubble retains its external momentum: {unit}"
        );
        Ok(unit)
    }

    fn factored_insertion(
        &self,
        insertion: &FeynmanDiagram,
        vakint: &Vakint,
        depth: usize,
        coefficients: &[Atom],
    ) -> Result<Atom> {
        use spenso::{
            network::tags::SPENSO_TAG,
            structure::representation::{Minkowski, RepName},
        };
        if coefficients.iter().all(|value| value.is_zero()) {
            return Ok(Atom::Zero);
        }
        let external = insertion
            .edges()
            .find_map(|(id, _, e)| e.external.as_ref().map(|_| id))
            .ok_or_else(|| eyre!("missing insertion momentum"))?;
        let q = |edge: usize| {
            function!(
                GS.emr_mom,
                edge,
                Minkowski {}.new_rep(GS.dim).to_symbolic([])
            )
        };
        let qh = q(external.0);
        let unit = self.unit_scalar_bubble(insertion, vakint, depth)?;
        let eps = GS.dim_epsilon;
        let pole = |a: &Atom| -> Result<Atom> { Ok(a.series(eps, 0, -1)?.to_atom().expand()) };
        let kb = pole(&unit)?;
        let x = (-function!(SPENSO_TAG.dot, &qh, &qh) / Atom::var(GS.mu_r_sq)).pow(-Atom::var(eps));
        let y = unit * x - &kb;
        // An attached coefficient e contributes e y^j + (e K(B)-K(e B)) y^(j-1).
        // e is regular and B-K(B) is finite. Horner evaluation keeps the line factor intact.
        let mut factor = Atom::Zero;
        let mut correction = Atom::Zero;
        for (j, e) in coefficients.iter().enumerate().rev() {
            factor = factor * &y + e + &correction;
            if j != 0 {
                let product = e * &kb;
                correction = &product - pole(&product)?;
            }
        }
        Ok(factor)
    }

    fn renormalized_insertion_uncached(
        &self,
        insertion: &FeynmanDiagram,
        vakint: &Vakint,
        parent_loops: usize,
    ) -> Result<Vec<Atom>> {
        let loops = insertion.loop_count();
        let started = Instant::now();
        let bare = self.integrate_cached(insertion, vakint, parent_loops)?;
        let bubble = scalar_bubble(insertion)?;
        let bubble_integral = self.integrate_cached(&bubble, vakint, parent_loops)?;
        let external = insertion
            .edges()
            .find_map(|(id, _, data)| data.external.as_ref().map(|_| id))
            .ok_or_else(|| eyre!("excised insertion has no external momentum"))?;
        let bubble_integral = bubble_integral
            .replace(function!(GS.emr_mom, 0, W_.x___))
            .with(function!(GS.emr_mom, external.0, W_.x___));
        let pole = |value: &Atom| -> Result<Atom> {
            Ok(value.series(GS.dim_epsilon, 0, -1)?.to_atom().expand())
        };
        // Scalar dimensions (traces, Dirac algebra) must be continued before
        // any Laurent operation; representation labels stay metadata.
        let continued = |value: &Atom| -> Result<Atom> {
            Ok(
                crate::uv::approx::integrated::Integrated::regulated_for_expansion(
                    value,
                    "renormalized_sigma",
                )?
                .into_owned(),
            )
        };
        let bare = continued(&bare)?;
        let bubble_integral = continued(&bubble_integral)?;
        let graph = insertion.to_gamma_loop_graph(None, true)?;
        let mut plan = RStarPlan::with_automatic_edge(&graph, insertion.model(), &self.settings)?;
        if let Some(cache) = &self.counterterm_cache {
            plan = plan.with_counterterm_cache(cache.clone());
        }
        let groups = plan.physical_uv_terms_by_loops(
            &graph,
            insertion.model(),
            vakint,
            &self.settings,
            parent_loops,
        )?;
        let groups = groups
            .into_iter()
            .map(|(k, value)| Ok((k, continued(&value)?)))
            .collect::<Result<std::collections::BTreeMap<_, _>>>()?;
        let mut total = bare.clone();
        for value in groups.values() {
            total += value;
        }
        let delta = -pole(&total)?;
        // ctilde_k, normalized by B^k.
        let mut normalized = vec![Atom::Zero; loops + 1];
        normalized[loops] = &bare / bubble_integral.pow(loops as i64);
        for (&k, value) in &groups {
            ensure!(
                k < loops,
                "a proper forest remainder has as many loops as the insertion"
            );
            normalized[k] = &normalized[k] + value / bubble_integral.pow(k as i64);
        }
        normalized[0] = &normalized[0] + &delta;
        let kb = pole(&bubble_integral)?;
        let binomial = |n: usize, k: usize| -> i64 {
            (0..k).fold(1i64, |acc, i| acc * (n - i) as i64 / (i + 1) as i64)
        };
        let mut result: Vec<Atom> = Vec::with_capacity(loops + 1);
        for j in 0..=loops {
            let mut d = Atom::Zero;
            for k in j..=loops {
                d = d + Atom::num(binomial(k, j)) * kb.pow((k - j) as i64) * &normalized[k];
            }
            // Exact through eps^(L-N+j-1) (DESIGN.md 3.4); higher orders of d_j
            // would need deeper insertion data and are never used by H_j.
            let d = d
                .series(
                    GS.dim_epsilon,
                    0,
                    parent_loops as i64 - loops as i64 + j as i64 - 1,
                )?
                .to_atom()
                .expand()
                .together();
            ensure!(
                d.series(GS.dim_epsilon, 0, -1)?
                    .to_atom()
                    .expand()
                    .together()
                    .is_zero(),
                "summed insertion coefficient d_{j} has epsilon poles: {d}"
            );
            ensure!(
                !d.contains_symbol(Symbol::LOG) && !d.contains_symbol(GS.mu_r_sq),
                "summed insertion coefficient d_{j} retains a scale: {d}"
            );
            result.push(d);
        }
        {
            // Coefficient on the first necklace vertex (the existing surrogate
            // topology): the necklace then represents e y^j + delta(e) y^(j-1),
            // delta(e) = e K(B) - K(e B), regular and local. Back-substitute so
            // that sum_j e_j [y^j + ...] = sum_j d_j y^j (DESIGN.md 3.1).
            let truncate = |value: Atom, j: usize| -> Result<Atom> {
                Ok(value
                    .series(
                        GS.dim_epsilon,
                        0,
                        parent_loops as i64 - loops as i64 + j as i64 - 1,
                    )?
                    .to_atom()
                    .expand()
                    .together())
            };
            let mut e = vec![Atom::Zero; loops + 1];
            e[loops] = result[loops].clone();
            for j in (0..loops).rev() {
                let next = &e[j + 1];
                let delta = (next * &kb - pole(&(next * &bubble_integral))?).expand();
                e[j] = truncate(&result[j] - delta, j)?;
                ensure!(
                    e[j].series(GS.dim_epsilon, 0, -1)?
                        .to_atom()
                        .expand()
                        .together()
                        .is_zero(),
                    "attached necklace coefficient e_{j} has epsilon poles"
                );
            }
            result = e;
        }
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            loops = loops,
            parent_loops = parent_loops,
            groups = groups.len(),
            seconds = started.elapsed().as_secs_f64(),
            "Computed a renormalized insertion in the necklace basis"
        );
        Ok(result)
    }
}
