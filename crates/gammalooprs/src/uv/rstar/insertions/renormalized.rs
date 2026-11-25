//! Reconstruct the renormalized massless self-energy in a basis of scalar
//! bubble chains. Each coefficient stays inside the pole map of its graph.
use super::*;

static BUBBLE_POWER: LazyLock<Atom> =
    LazyLock::new(|| Atom::var(symbol!("renormalized_insertion_bubble_power")));

impl PropagatorInsertionEvaluator {
    fn insertion_coefficients(
        &self,
        insertion: &FeynmanDiagram,
        vakint: &Vakint,
        depth: usize,
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
            format!(
                "renormalized-sigma:attached-v1:depth={depth}:closure={:?}:{}",
                self.closure, self.settings_key
            ),
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
        mut project: impl FnMut(&Graph, &Model, &RStarPlan) -> Result<Atom>,
    ) -> Result<Option<Atom>> {
        ensure!(
            diagram.cuts().is_empty(),
            "insertion extraction requires an uncut amplitude"
        );
        let insertions = selected_insertions(diagram)?;
        if insertions.is_empty() {
            return Ok(None);
        }
        let mut terms = vec![(
            diagram.clone(),
            Vec::<DiagramRelabeling>::new(),
            Vec::<(VertexId, usize)>::new(),
        )];
        for insertion in &insertions {
            let excised = diagram.excise_with_mapping(&insertion.support, &[])?;
            let coefficients =
                self.insertion_coefficients(&excised.diagram, vakint, diagram.loop_count())?;
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
                    .filter(|(_, c)| !c.is_zero())
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
