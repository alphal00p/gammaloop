//! Direct recognition of five-loop master inputs.
//!
//! A scalar term that is a package master in any momentum routing, line
//! labelling, numerator spelling or with its non-positive lines given either
//! as zero/negative-power propagators or contracted (the form a graph
//! generator emits), is rewritten as that master's own common-basis
//! representation before witness selection, so it is reduced as the master
//! itself instead of by an IBP walk from a symmetric image.
//!
//! The witness enumeration only sees the automorphisms of the graph that the
//! input is written on. A contracted input is matched to the census embedding
//! of its graph class and a parent-form input to its own parent; when the
//! master lives in another embedding of the same graph (BASIS-v2 P387: census
//! `I5L_32745_pinch_4_6` versus the master's slots), no enumerated witness is
//! the master and the terminal-first selection cannot help.
//!
//! Recognition: contract every non-positive line (negative powers become
//! numerator factors), compare the power-colored canonical incidence graph
//! with those of the package masters, and for every verified power-preserving
//! unimodular momentum map onto a candidate master lower the transported
//! numerator over the master's denominators. The term is rewritten only when
//! that lowering lies on remaining terminals only and the ordinary witness
//! selection of the rewritten term confirms it. Every other input keeps its
//! original path; the rewritten term is the same integral (FeynKit verifies
//! every denominator and the Jacobian of the map).

use super::*;
use crate::topologies::five_loop::data as census;
use crate::utils::vakint_macros::vk_parse;

/// One package master in its own contracted long form.
pub(crate) struct MasterForm {
    key: IntegralKey,
    parent: usize,
    canonical: symbolica::graph::Graph<u8, symbolica::graph::Empty>,
    contracted: Atom,
}

fn first_parent(key: &IntegralKey) -> Option<usize> {
    census::PARENTS.iter().position(|parent| {
        key.powers()
            .iter()
            .enumerate()
            .all(|(slot, power)| *power <= 0 || parent.physical_slots.contains(&slot))
    })
}

fn momentum(slot: usize) -> Atom {
    common_momentum_rows()[slot]
        .iter()
        .enumerate()
        .map(|(axis, coefficient)| Atom::num(*coefficient) * S.k.call(axis + 1))
        .sum()
}

fn incidence(parent: usize, slot: usize) -> Option<(usize, usize)> {
    let vertices = &census::PARENTS[parent].vertices;
    let signed = (slot + 1) as i64;
    Some((
        vertices.iter().position(|vertex| vertex.contains(&signed))?,
        vertices.iter().position(|vertex| vertex.contains(&-signed))?,
    ))
}

fn find<K: Ord + Copy>(root: &mut BTreeMap<K, K>, vertex: K) -> K {
    let mut current = vertex;
    while let Some(&next) = root.get(&current) {
        if next == current {
            break;
        }
        current = next;
    }
    root.insert(vertex, current);
    current
}

fn union<K: Ord + Copy>(root: &mut BTreeMap<K, K>, left: K, right: K) {
    let (left, right) = (find(root, left), find(root, right));
    if left != right {
        root.insert(left, right);
    }
}

/// The master's positive lines in its first defining parent with every other
/// physical line contracted, at unit mass.
fn contracted_form(key: &IntegralKey, parent: usize) -> Option<Atom> {
    let physical = &census::PARENTS[parent].physical_slots;
    let mut root = BTreeMap::new();
    for slot in physical {
        if key.powers()[*slot] <= 0 {
            let (left, right) = incidence(parent, *slot)?;
            union(&mut root, left, right);
        }
    }
    let mut labels = BTreeMap::new();
    let mut integral = Atom::one();
    let mut id = 0usize;
    for slot in physical {
        let power = key.powers()[*slot];
        if power <= 0 {
            continue;
        }
        let (left, right) = incidence(parent, *slot)?;
        let (left, right) = (find(&mut root, left), find(&mut root, right));
        let next = labels.len() + 1;
        let left = *labels.entry(left).or_insert(next);
        let next = labels.len() + 1;
        let right = *labels.entry(right).or_insert(next);
        id += 1;
        integral *= S.prop.call((
            id,
            S.edge.call((left, right)),
            momentum(*slot),
            Atom::one(),
            Atom::num(power),
        ));
    }
    Some(S.topo.call(integral))
}

/// The master's defining-parent long form (every physical line, non-positive
/// ones at power zero) at the input's mass, and its positive base key.
fn parent_form(key: &IntegralKey, parent: usize, mass: &Atom) -> Option<(Atom, IntegralKey)> {
    let mut integral = Atom::one();
    for (index, slot) in census::PARENTS[parent].physical_slots.iter().enumerate() {
        let (left, right) = incidence(parent, *slot)?;
        integral *= S.prop.call((
            index + 1,
            S.edge.call((left + 1, right + 1)),
            momentum(*slot),
            mass,
            Atom::num(key.powers()[*slot].max(0)),
        ));
    }
    let base = IntegralKey::try_new(
        key.powers()
            .iter()
            .map(|power| (*power).max(0))
            .collect::<Vec<_>>(),
    )
    .ok()?;
    Some((S.topo.call(integral), base))
}

/// The input with every non-positive line contracted (vertices merged) and
/// negative powers moved into the scalar-product numerator, with its common
/// mass; `None` if it is not a unit-shaped long form with integer powers and
/// one nonzero mass, or if its numerator keeps open loop-momentum indices.
fn contract(
    settings: &VakintSettings,
    integral: AtomView,
    numerator: AtomView,
) -> Result<Option<(Atom, Atom, Atom)>, VakintError> {
    let pattern = vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,pow_)")
        .unwrap()
        .to_pattern();
    let mut lines = Vec::new();
    let mut original = Atom::one();
    for matched in integral.pattern_match(&pattern, None, None) {
        let (Some(id), Some(left), Some(right), Some(power)) = (
            crate::get_integer_from_atom(matched[&vk_symbol!("id_")].as_view()),
            crate::get_integer_from_atom(matched[&vk_symbol!("nl_")].as_view()),
            crate::get_integer_from_atom(matched[&vk_symbol!("nr_")].as_view()),
            crate::get_integer_from_atom(matched[&vk_symbol!("pow_")].as_view()),
        ) else {
            return Ok(None);
        };
        let momentum = matched[&vk_symbol!("q_")].clone();
        let mass = matched[&vk_symbol!("m_")].clone();
        original *= S.prop.call((
            Atom::num(id),
            S.edge.call((Atom::num(left), Atom::num(right))),
            &momentum,
            &mass,
            Atom::num(power),
        ));
        lines.push((left, right, momentum, mass, power));
    }
    if lines.is_empty() || S.topo.call(original).as_view() != integral {
        return Ok(None);
    }
    let mass = lines[0].3.clone();
    if mass.is_zero() || lines.iter().any(|line| line.3 != mass) {
        return Ok(None);
    }
    let mut numerator = VakintTerm::scalar_product_form(settings, numerator)?;
    let mut open_index = false;
    numerator.as_view().visitor(&mut |value| {
        if let AtomView::Fun(function) = value
            && function.get_symbol() == S.k
            && function.get_nargs() != 1
        {
            open_index = true;
        }
        true
    });
    if open_index {
        return Ok(None);
    }
    let mut root = BTreeMap::new();
    for (left, right, momentum, _, power) in &lines {
        if *power <= 0 {
            union(&mut root, *left, *right);
        }
        if *power < 0 {
            let exponent = u32::try_from(power.unsigned_abs()).map_err(failure)?;
            numerator =
                numerator * (S.dot.call((momentum, momentum)) - &mass).pow(exponent);
        }
    }
    let mut labels = BTreeMap::new();
    let mut contracted = Atom::one();
    let mut id = 0usize;
    for (left, right, momentum, _, power) in &lines {
        if *power <= 0 {
            continue;
        }
        let (left, right) = (find(&mut root, *left), find(&mut root, *right));
        let next = labels.len() + 1;
        let left = *labels.entry(left).or_insert(next);
        let next = labels.len() + 1;
        let right = *labels.entry(right).or_insert(next);
        id += 1;
        contracted *= S.prop.call((
            id,
            S.edge.call((left, right)),
            momentum,
            Atom::one(),
            Atom::num(*power),
        ));
    }
    Ok(Some((S.topo.call(contracted), numerator, mass)))
}

/// Resource exhaustion (RAM/time budget, explicit diagnostic counts) is the
/// caller's error; every other failure only disables the recognition.
fn resource_error(error: &VakintError) -> bool {
    matches!(
        error,
        VakintError::ResourceLimit(_)
            | VakintError::MomentumMatchResourceLimit { .. }
            | VakintError::RustRedEvaluation(RustRedEvaluationError::ResourceLimit(_))
    )
}

/// Substitute the loop momenta k(1..) of `numerator` simultaneously.
fn transport_numerator(numerator: &Atom, images: &[Atom]) -> Atom {
    let placeholder = |axis: usize| Atom::var(vk_symbol!(format!("vakint_master_transport_{axis}")));
    let mut out = numerator.clone();
    for axis in 1..=images.len() {
        out = out
            .replace(S.k.call(axis).to_pattern())
            .with(placeholder(axis).to_pattern());
    }
    for (index, image) in images.iter().enumerate() {
        out = out
            .replace(placeholder(index + 1).to_pattern())
            .with(image.to_pattern());
    }
    out
}

impl SavedApplication {
    /// Remaining terminals of the package's own family.
    fn own_terminals(&self) -> Result<BTreeSet<IntegralKey>, RustRedEvaluationError> {
        let prepared = self.prepared.lock().map_err(failure)?;
        Ok(prepared
            .remaining_terminals()
            .iter()
            .filter(|terminal| terminal.family_fingerprint() == prepared.family().fingerprint())
            .map(|terminal| terminal.integral().clone())
            .collect())
    }

    fn master_forms(&self) -> &[MasterForm] {
        self.master_forms.get_or_init(|| {
            let started = std::time::Instant::now();
            let terminals = self.own_terminals().unwrap_or_default();
            let mut forms = Vec::new();
            let mut skipped = 0;
            for key in terminals {
                let form = first_parent(&key).and_then(|parent| {
                    let contracted = contracted_form(&key, parent)?;
                    let canonical =
                        crate::topologies::power_colored_canonical_graph(contracted.as_view(), 5)
                            .ok()
                            .flatten()?;
                    Some(MasterForm {
                        key: key.clone(),
                        parent,
                        canonical,
                        contracted,
                    })
                });
                match form {
                    Some(form) => forms.push(form),
                    None => skipped += 1,
                }
            }
            log::debug!(target: "vakint::saved_application", "phase=master_forms masters={} skipped={skipped} elapsed_seconds={:.6}", forms.len(), started.elapsed().as_secs_f64());
            forms
        })
    }
}

impl Vakint {
    /// The exact lowered keys of a scalar numerator under one witness, without
    /// the saved entry-union admission check of [`Self::matched_input_key_with`]
    /// (a remaining terminal reduces to itself whether or not its own key is
    /// inside the recorded entry predicates).
    pub(super) fn witness_keys_unchecked(
        &self,
        settings: &VakintSettings,
        application: &SavedApplication,
        numerator: AtomView,
        matched: &ReplacementRules,
    ) -> Result<Vec<IntegralKey>, VakintError> {
        let mut matched = matched.clone();
        matched.apply_replacement_rules()?;
        let mut term = crate::VakintTerm {
            integral: Atom::one(),
            numerator: numerator.to_owned(),
            vectors: Vec::new(),
        };
        term.apply_scalar_numerator_replacement_rules(&matched, settings)?;
        let descriptor = select(&matched.canonical_topology)?;
        let input = MatchedScalarInput::prepare(
            term.numerator.as_view(),
            &matched,
            &descriptor.momenta,
            descriptor.physical_slots,
            15,
            settings,
            ScalarNumeratorForm::PreserveProducts,
        )?;
        let prepared = application.prepared.lock().map_err(failure)?;
        let scalar = FamilyScalarNumeratorService::try_new(
            prepared.family(),
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(native_failure)?;
        let lowered = scalar
            .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
            .map_err(native_failure)?;
        Ok(lowered
            .terms()
            .iter()
            .map(|term| term.integral().clone())
            .collect())
    }

    /// Rewrite a five-loop scalar term that equals a package master as that
    /// master's own common-basis representation (see the module
    /// documentation) and return the witness to use for the rewritten term:
    /// the ordinary first (incidence) witness of the rewritten term, which is
    /// the identity map onto its defining parent's census topology, accepted
    /// only if the numerator lowers under it onto remaining terminals only. A
    /// terminal-only key of one integral is unique, so this is the key the
    /// terminal-first selection returns, without enumerating automorphisms.
    /// `None` leaves the term untouched.
    pub(crate) fn transport_rustred_master_term(
        &self,
        settings: &VakintSettings,
        integral: &mut Atom,
        numerator: &mut Atom,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        // Recognition is an optimization: a structural failure on an input
        // the ordinary path accepts leaves the term untouched; only resource
        // limits propagate.
        match self.transport_rustred_master_term_checked(settings, integral, numerator) {
            Err(error) if !resource_error(&error) => {
                log::debug!(target: "vakint::saved_application", "phase=master_transport outcome=skipped error={error}");
                Ok(None)
            }
            result => result,
        }
    }

    fn transport_rustred_master_term_checked(
        &self,
        settings: &VakintSettings,
        integral: &mut Atom,
        numerator: &mut Atom,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        if settings.rustred_artifact_path.is_none() && !self.has_rustred_loop(5) {
            return Ok(None);
        }
        let started = std::time::Instant::now();
        let outcome = |outcome: &str, candidates: usize, maps: usize| {
            log::debug!(target: "vakint::saved_application", "phase=master_transport outcome={outcome} candidates={candidates} maps={maps} elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        };
        let Some((contracted, contracted_numerator, mass)) =
            contract(settings, integral.as_view(), numerator.as_view())?
        else {
            outcome("not_a_long_form", 0, 0);
            return Ok(None);
        };
        let canonical = match crate::topologies::power_colored_canonical_graph(contracted.as_view(), 5) {
            Ok(Some(canonical)) => canonical,
            Ok(None) | Err(VakintError::InvalidIntegralFormat(_)) => {
                outcome("no_canonical_form", 0, 0);
                return Ok(None);
            }
            Err(error) => return Err(error),
        };
        let application = self.saved_application(settings)?;
        let candidates = application
            .master_forms()
            .iter()
            .filter(|form| form.canonical == canonical)
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            outcome("no_candidate", 0, 0);
            return Ok(None);
        }
        let terminals = application.own_terminals()?;
        let mut tried = 0;
        for form in &candidates {
            let Some((target, _)) = parent_form(&form.key, form.parent, &mass) else {
                continue;
            };
            let Some(witness) = self
                .topologies
                .match_topologies_to_user_input(target.as_view(), settings.allow_unknown_integrals)?
                .filter(|witness| matches!(witness.canonical_topology, Topology::FiveLoop(_)))
            else {
                continue;
            };
            // The isomorphisms are offered lazily: the first map whose
            // transported numerator lowers onto remaining terminals only ends
            // the search (highly symmetric product graphs have many).
            let mut accepted = None;
            let found = crate::topologies::find_power_preserving_momentum_map(
                contracted.as_view(),
                form.contracted.as_view(),
                5,
                |images| {
                    tried += 1;
                    let transported = transport_numerator(&contracted_numerator, images);
                    let keys = match self.witness_keys_unchecked(
                        settings,
                        &application,
                        transported.as_view(),
                        &witness,
                    ) {
                        Ok(keys) => keys,
                        Err(VakintError::RustRedEvaluation(RustRedEvaluationError::ResourceLimit(
                            limit,
                        ))) => return Err(RustRedEvaluationError::ResourceLimit(limit).into()),
                        Err(_) => return Ok(false),
                    };
                    if !keys.is_empty() && keys.iter().all(|key| terminals.contains(key)) {
                        accepted = Some((transported, keys.len()));
                        return Ok(true);
                    }
                    Ok(false)
                },
            );
            match found {
                Ok(Some(true)) => {}
                Ok(_) | Err(VakintError::InvalidIntegralFormat(_)) => continue,
                Err(error) => return Err(error),
            }
            let Some((transported, keys)) = accepted else {
                continue;
            };
            log::debug!(target: "vakint::five_loop_match", "phase=master_transport_selected master={:?} keys={keys}", form.key.powers());
            *integral = target;
            *numerator = transported;
            outcome("master", candidates.len(), tried);
            return Ok(Some(witness));
        }
        outcome("not_a_master", candidates.len(), tried);
        Ok(None)
    }
}

/// Contracted long form (unit mass) of a four-loop program's integral with
/// `powers` on its physical slots, every non-positive line contracted.
fn four_loop_contracted(edges: &[(i64, i64)], momenta: &[Atom], powers: &[i64]) -> Atom {
    let mut root = BTreeMap::new();
    for ((left, right), power) in edges.iter().zip(powers) {
        if *power <= 0 {
            union(&mut root, *left, *right);
        }
    }
    let mut labels = BTreeMap::new();
    let mut integral = Atom::one();
    let mut id = 0usize;
    for (((left, right), momentum), power) in edges.iter().zip(momenta).zip(powers) {
        if *power <= 0 {
            continue;
        }
        let (left, right) = (find(&mut root, *left), find(&mut root, *right));
        let next = labels.len() + 1;
        let left = *labels.entry(left).or_insert(next);
        let next = labels.len() + 1;
        let right = *labels.entry(right).or_insert(next);
        id += 1;
        integral *= S.prop.call((id, S.edge.call((left, right)), momentum, Atom::one(), Atom::num(*power)));
    }
    S.topo.call(integral)
}

/// The program's own long form of an integral (every physical slot, in slot
/// order, at mass `mass`), as its shipped test inputs write it.
fn four_loop_long_form(edges: &[(i64, i64)], momenta: &[Atom], powers: &[i64], mass: &Atom) -> Atom {
    let mut integral = Atom::one();
    for (slot, (((left, right), momentum), power)) in edges.iter().zip(momenta).zip(powers).enumerate() {
        integral *= S.prop.call((
            slot + 1,
            S.edge.call((Atom::num(*left), Atom::num(*right))),
            momentum,
            mass,
            Atom::num(*power),
        ));
    }
    S.topo.call(integral)
}

impl Vakint {
    /// Rewrite a four-loop scalar integral that equals a shipped four-loop
    /// master (positive powers only, no numerator) under a verified
    /// power-preserving momentum map as that master's own long form, so the
    /// ordinary matcher sees the program's slot order (input propagator ids in
    /// another order made matching take ~20-30 s and could select a symmetric
    /// image that the program reduces by a walk). Returns whether it did.
    pub(crate) fn transport_four_loop_master(
        &self,
        settings: &VakintSettings,
        integral: &mut Atom,
        numerator: AtomView,
    ) -> Result<bool, VakintError> {
        match self.transport_four_loop_master_checked(settings, integral, numerator) {
            Err(error) if !resource_error(&error) => {
                log::debug!(target: "vakint::saved_application", "phase=four_loop_master_transport outcome=skipped error={error}");
                Ok(false)
            }
            result => result,
        }
    }

    fn transport_four_loop_master_checked(
        &self,
        settings: &VakintSettings,
        integral: &mut Atom,
        numerator: AtomView,
    ) -> Result<bool, VakintError> {
        let Some(data) = self.rustred_data.as_ref() else {
            return Ok(false);
        };
        let (Some(programs), Some(pack)) = (data.four_loop.as_ref(), data.loops.get(&4)) else {
            return Ok(false);
        };
        let started = std::time::Instant::now();
        // Already written in a program's own slot order: nothing to do.
        if let AtomView::Fun(topo) = integral.as_view()
            && topo.get_symbol() == S.topo
            && topo.get_nargs() == 1
        {
            let pattern = vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,pow_)").unwrap().to_pattern();
            let mut lines = integral
                .as_view()
                .pattern_match(&pattern, None, None)
                .map(|matched| {
                    (
                        crate::get_integer_from_atom(matched[&vk_symbol!("id_")].as_view()),
                        crate::get_integer_from_atom(matched[&vk_symbol!("nl_")].as_view()),
                        crate::get_integer_from_atom(matched[&vk_symbol!("nr_")].as_view()),
                        matched[&vk_symbol!("q_")].clone(),
                    )
                })
                .collect::<Vec<_>>();
            lines.sort_by_key(|line| line.0);
            let mut index = 0;
            while let (Some((_, edges)), Some(momenta)) =
                (programs.physical_edges(index), programs.physical_momenta(index))
            {
                if lines.len() == momenta.len()
                    && lines.iter().zip(edges.iter().zip(momenta)).enumerate().all(
                        |(slot, ((id, left, right, q), ((l, r), momentum)))| {
                            *id == Some(slot as i64 + 1) && *left == Some(*l) && *right == Some(*r) && q == momentum
                        },
                    )
                {
                    return Ok(false);
                }
                index += 1;
            }
        }
        let Some((contracted, contracted_numerator, mass)) =
            contract(settings, integral.as_view(), numerator)?
        else {
            return Ok(false);
        };
        if contracted_numerator.contains_symbol(S.k) {
            return Ok(false);
        }
        let canonical = match crate::topologies::power_colored_canonical_graph(contracted.as_view(), 4) {
            Ok(Some(canonical)) => canonical,
            Ok(None) | Err(VakintError::InvalidIntegralFormat(_)) => return Ok(false),
            Err(error) => return Err(error),
        };
        let forms = programs.master_forms.get_or_init(|| {
            let mut forms = Vec::new();
            for master in &pack.manifest.masters {
                let (Some(alias), Some(powers)) = (
                    master.get("alias").and_then(|alias| alias.as_str()),
                    master
                        .get("powers")
                        .and_then(|powers| serde_json::from_value::<Vec<i64>>(powers.clone()).ok()),
                ) else {
                    continue;
                };
                let label = format!("four-loop-{}", alias.trim_start_matches("4L_"));
                let mut index = 0;
                while let Some((program_label, edges)) = programs.physical_edges(index) {
                    if program_label == label {
                        let momenta = programs.physical_momenta(index).unwrap_or_default();
                        if powers.len() >= momenta.len()
                            && powers[momenta.len()..].iter().all(|power| *power == 0)
                        {
                            let contracted = four_loop_contracted(edges, momenta, &powers[..momenta.len()]);
                            if let Ok(Some(canonical)) =
                                crate::topologies::power_colored_canonical_graph(contracted.as_view(), 4)
                            {
                                forms.push(super::super::four_loop::FourLoopMasterForm {
                                    program: index,
                                    powers: powers[..momenta.len()].to_vec(),
                                    canonical,
                                    contracted,
                                });
                            }
                        }
                        break;
                    }
                    index += 1;
                }
            }
            log::debug!(target: "vakint::saved_application", "phase=four_loop_master_forms masters={}", forms.len());
            forms
        });
        for form in forms.iter().filter(|form| form.canonical == canonical) {
            let found = crate::topologies::find_power_preserving_momentum_map(
                contracted.as_view(),
                form.contracted.as_view(),
                4,
                |_| Ok(true),
            );
            match found {
                Ok(Some(true)) => {}
                Ok(_) | Err(VakintError::InvalidIntegralFormat(_)) => continue,
                Err(error) => return Err(error),
            }
            let (Some((_, edges)), Some(momenta)) = (
                programs.physical_edges(form.program),
                programs.physical_momenta(form.program),
            ) else {
                continue;
            };
            *integral = four_loop_long_form(edges, momenta, &form.powers, &mass);
            log::debug!(target: "vakint::saved_application", "phase=four_loop_master_transport program={} powers={:?} elapsed_seconds={:.6}", form.program, form.powers, started.elapsed().as_secs_f64());
            return Ok(true);
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A common-basis integrand in its first defining parent, as the B1-K
    /// probe and `admission_census_tests::input_for_key` build it.
    fn parent_input(key: &IntegralKey) -> Option<(Atom, Atom)> {
        let (integral, _) = parent_form(key, first_parent(key)?, &Atom::one())?;
        let mut numerator = Atom::one();
        for (slot, power) in key.powers().iter().enumerate() {
            if *power < 0 {
                numerator *= (S.dot.call((momentum(slot), momentum(slot))) - Atom::one())
                    .pow(power.unsigned_abs() as u32);
            }
        }
        Some((integral, numerator))
    }

    /// A unimodular loop-momentum relabelling: k1->k2, k2->k1+k3, k3->-k3,
    /// k4->k5, k5->k4-k2.
    fn relabel(atom: &Atom) -> Atom {
        let k = |axis: usize| S.k.call(axis);
        transport_numerator(atom, &[k(2), k(1) + k(3), -k(3), k(5), k(4) - k(2)])
    }

    fn canonical(input: &Atom) -> symbolica::graph::Graph<u8, symbolica::graph::Empty> {
        crate::topologies::power_colored_canonical_graph(input.as_view(), 5)
            .unwrap()
            .unwrap()
    }

    /// Package-free: contraction, power-colored canonical forms and the
    /// verified power-preserving maps of the BASIS-v2 master P387.
    #[test]
    fn power_preserving_maps_find_relabelled_copies_and_respect_powers() {
        Vakint::initialize_vakint_symbols();
        let key = IntegralKey::try_new(vec![1, 1, 0, 1, 2, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1]).unwrap();
        let parent = first_parent(&key).unwrap();
        let contracted = contracted_form(&key, parent).unwrap();
        // Contracting the parent form with the negative power written on its
        // line reproduces the master's contracted graph and moves that power
        // into the numerator, as the B1-K probe writes it.
        let (_, numerator) = parent_input(&key).unwrap();
        let mut signed = Atom::one();
        for (index, slot) in census::PARENTS[parent].physical_slots.iter().enumerate() {
            let (left, right) = incidence(parent, *slot).unwrap();
            signed *= S.prop.call((
                index + 1,
                S.edge.call((left + 1, right + 1)),
                momentum(*slot),
                Atom::one(),
                Atom::num(key.powers()[*slot]),
            ));
        }
        let signed = S.topo.call(signed);
        let one = Atom::one();
        let (from_parent, moved, mass) =
            contract(&VakintSettings::default(), signed.as_view(), one.as_view())
                .unwrap()
                .unwrap();
        assert_eq!(mass, Atom::one());
        assert_eq!(canonical(&from_parent), canonical(&contracted));
        assert_eq!(moved, numerator);
        // A relabelled copy maps back, by maps that are verified and unimodular.
        let relabelled = relabel(&contracted);
        assert_eq!(canonical(&relabelled), canonical(&contracted));
        let maps = crate::topologies::power_preserving_momentum_maps(
            relabelled.as_view(),
            contracted.as_view(),
            5,
        )
        .unwrap()
        .unwrap();
        assert!(!maps.is_empty());
        // Transporting the relabelled numerator by any map gives back a
        // numerator over the master's lines (the inverse relabelling).
        assert!(maps.iter().all(|images| images.len() == 5));
        // Another power pattern on the same graph has a different colored
        // canonical form and no power-preserving map.
        let other = IntegralKey::try_new(vec![1, 1, 0, 1, 3, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1]).unwrap();
        let other = contracted_form(&other, first_parent(&other).unwrap()).unwrap();
        assert_ne!(canonical(&other), canonical(&contracted));
        assert!(
            crate::topologies::power_preserving_momentum_maps(other.as_view(), contracted.as_view(), 5)
                .unwrap()
                .unwrap()
                .is_empty()
        );
    }

    /// Every package master in its parent form, contracted, and relabelled in
    /// both forms is rewritten as the master itself (its identity witness
    /// lowers to exactly the master key); the P387-sector non-masters V1-V3
    /// (numerator squared, extra dots) keep their ordinary path; any rewrite
    /// of a dotted or numerator variant of a master lands on a single master.
    #[test]
    #[ignore = "loads the official five-loop saved application (minutes)"]
    fn every_master_input_form_is_transported_to_the_master_and_variants_are_not() {
        let directory = super::super::super::test_data::external_path_through(5).join("5L/ibps");
        let mut settings = VakintSettings::default();
        let engine = super::super::super::test_data::external_engine()
            .with_rustred_application(&mut settings, &directory)
            .unwrap();
        let application = engine.saved_application(&settings).unwrap();
        let terminals = application.own_terminals().unwrap();
        let transported_keys = |integral: &Atom, numerator: &Atom| {
            let (mut integral, mut numerator) = (integral.clone(), numerator.clone());
            engine
                .transport_rustred_master_term(&settings, &mut integral, &mut numerator)
                .unwrap()
                .map(|witness| {
                    engine
                        .witness_keys_unchecked(&settings, &application, numerator.as_view(), &witness)
                        .unwrap()
                })
        };
        let mut checked = 0;
        for key in &terminals {
            let (integral, numerator) = parent_input(key).unwrap();
            let (contracted, contracted_numerator, _) =
                contract(&settings, integral.as_view(), numerator.as_view())
                    .unwrap()
                    .unwrap();
            for (form, (integral, numerator)) in [
                (integral.clone(), numerator.clone()),
                (contracted.clone(), contracted_numerator.clone()),
                (relabel(&contracted), relabel(&contracted_numerator)),
                (relabel(&integral), relabel(&numerator)),
            ]
            .iter()
            .enumerate()
            {
                assert_eq!(
                    transported_keys(integral, numerator),
                    Some(vec![key.clone()]),
                    "master {:?}, form {form}",
                    key.powers()
                );
            }
            let first_positive = key.powers().iter().position(|power| *power > 0).unwrap();
            let first_free = key.powers().iter().position(|power| *power <= 0).unwrap();
            for (slot, delta) in [(first_positive, 1), (first_free, -1)] {
                let mut powers = key.powers().to_vec();
                powers[slot] += delta;
                let variant = IntegralKey::try_new(powers).unwrap();
                let Some((integral, numerator)) = parent_input(&variant) else {
                    continue;
                };
                if let Some(keys) = transported_keys(&integral, &numerator) {
                    assert!(
                        keys.len() == 1 && terminals.contains(&keys[0]),
                        "variant {:?} rewritten onto {keys:?}",
                        variant.powers()
                    );
                }
            }
            checked += 1;
        }
        assert_eq!(checked, terminals.len());
        for powers in [
            [1, 1, 0, 1, 2, -2, 1, 1, 1, 1, 0, 1, 0, 0, 1],
            [1, 1, 0, 1, 3, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1],
            [2, 1, 0, 1, 2, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1],
        ] {
            let variant = IntegralKey::try_new(powers.to_vec()).unwrap();
            let (integral, numerator) = parent_input(&variant).unwrap();
            let (contracted, contracted_numerator, _) =
                contract(&settings, integral.as_view(), numerator.as_view())
                    .unwrap()
                    .unwrap();
            assert_eq!(transported_keys(&integral, &numerator), None, "{powers:?}");
            assert_eq!(
                transported_keys(&contracted, &contracted_numerator),
                None,
                "{powers:?} contracted"
            );
        }
    }

    /// Every shipped four-loop master with reversed propagator ids or under a
    /// unimodular loop relabelling is rewritten into the program form of a
    /// shipped master (the data ships some integrals in several families), the
    /// program form itself is left untouched, and a dotted variant is
    /// rewritten only onto a shipped master's form.
    #[test]
    #[ignore = "loads the official four-loop data"]
    fn every_four_loop_master_form_is_written_in_its_program_form() {
        let engine = super::super::super::test_data::external_engine();
        let settings = VakintSettings::default();
        let data = engine.rustred_data.as_ref().unwrap();
        let programs = data.four_loop.as_ref().unwrap();
        let k = |axis: usize| S.k.call(axis);
        let relabel4 = |atom: &Atom| transport_numerator(atom, &[k(2), k(1) + k(3), -k(3), k(4) - k(1)]);
        let mut checked = 0;
        for master in &data.loops[&4].manifest.masters {
            let alias = master["alias"].as_str().unwrap();
            let powers: Vec<i64> = serde_json::from_value(master["powers"].clone()).unwrap();
            let label = format!("four-loop-{}", alias.trim_start_matches("4L_"));
            let index = (0..4).find(|index| programs.physical_edges(*index).unwrap().0 == label).unwrap();
            let (_, edges) = programs.physical_edges(index).unwrap();
            let momenta = programs.physical_momenta(index).unwrap();
            let own = four_loop_long_form(edges, momenta, &powers[..momenta.len()], &Atom::one());
            let reversed = {
                let n = momenta.len();
                let mut integral = Atom::one();
                for (slot, (((left, right), momentum), power)) in
                    edges.iter().zip(momenta).zip(&powers).enumerate()
                {
                    integral *= S.prop.call((
                        n - slot,
                        S.edge.call((Atom::num(*left), Atom::num(*right))),
                        momentum,
                        Atom::one(),
                        Atom::num(*power),
                    ));
                }
                S.topo.call(integral)
            };
            // Program forms of every shipped master (the data ships some integrals
            // in several families: a rewrite may land on any of them).
            let shipped = |integral: &Atom| {
                programs.master_forms.get().unwrap().iter().any(|form| {
                    let (_, edges) = programs.physical_edges(form.program).unwrap();
                    let momenta = programs.physical_momenta(form.program).unwrap();
                    four_loop_long_form(edges, momenta, &form.powers, &Atom::one()) == *integral
                })
            };
            for (form, input) in [own.clone(), reversed, relabel4(&own)].into_iter().enumerate() {
                let mut integral = input;
                let one = Atom::one();
                let moved = engine.transport_four_loop_master(&settings, &mut integral, one.as_view()).unwrap();
                if form == 0 {
                    // The program form itself is left untouched.
                    assert!(!moved && integral == own, "{alias} {powers:?} form 0");
                } else {
                    assert!(moved && shipped(&integral), "{alias} {powers:?} form {form}: {integral}");
                }
            }
            // A dotted variant is rewritten only onto the program form of a
            // shipped master (it can equal one only under a symmetry).
            let mut dotted = powers[..momenta.len()].to_vec();
            let first = dotted.iter().position(|power| *power > 0).unwrap();
            dotted[first] += 1;
            let mut integral = four_loop_long_form(edges, momenta, &dotted, &Atom::one());
            let one = Atom::one();
            if engine.transport_four_loop_master(&settings, &mut integral, one.as_view()).unwrap() {
                assert!(shipped(&integral), "{alias} dotted {dotted:?} rewritten onto a non-master");
            }
            checked += 1;
        }
        assert_eq!(checked, data.loops[&4].manifest.masters.len());
    }
}
