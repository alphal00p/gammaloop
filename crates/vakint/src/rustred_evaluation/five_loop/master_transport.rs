//! Direct recognition of five-loop and four-loop master inputs.
//!
//! A scalar term that is a package master in any momentum routing, line
//! labelling, numerator spelling, graph embedding or with its non-positive
//! lines given either as zero/negative-power propagators or contracted (the
//! form a graph generator emits), is rewritten as that master's own
//! common-basis representation before witness selection, so it is reduced as
//! the master itself instead of by an IBP walk from a symmetric image.
//!
//! The witness enumeration only sees the automorphisms of the graph that the
//! input is written on. A contracted input is matched to the census embedding
//! of its graph class and a parent-form input to its own parent; when the
//! master lives in another embedding of the same graph (BASIS-v2 P387: census
//! `I5L_32745_pinch_4_6` versus the master's slots), or in a 2-isomorphic one
//! (a Whitney twist at a 2-vertex cut: same line momenta, non-isomorphic
//! incidence graph), no enumerated witness is the master and the
//! terminal-first selection cannot help.
//!
//! Recognition: contract every non-positive line (negative powers become
//! numerator factors) and compare the power-colored canonical incidence graph
//! with those of the package masters; masters on another graph with the same
//! power-colored line-matroid profile are compared at the matroid level. For
//! every verified power-preserving unimodular momentum map onto a candidate
//! (FeynKit verifies every denominator and the Jacobian), the transported
//! numerator is lowered over the master's denominators under the identity
//! witness of the master's defining parent. The term is rewritten as that
//! parent form, with that witness, when the lowering lies on remaining
//! terminals only; no automorphism enumeration and no entry admission are
//! involved, and a request made of remaining terminals only is then built
//! without the native walk (see `SavedApplication::evaluate_ordered_batch`).
//! The rewritten term is the same integral. Every other input keeps its
//! original path; a candidate whose first lowering has another numerator
//! degree than the master is abandoned (every map preserves the degree), and
//! non-master outcomes are memoized per contracted input.

use super::*;
use crate::topologies::five_loop::data as census;
use crate::utils::vakint_macros::vk_parse;

/// One package master in its own contracted long form.
pub(crate) struct MasterForm {
    key: IntegralKey,
    parent: usize,
    canonical: symbolica::graph::Graph<u8, symbolica::graph::Empty>,
    /// Power-colored profile of its line matroid (matroid-level candidates).
    profile: Option<crate::topologies::MatroidProfile>,
    /// Numerator degree in scalar products: the sum of its negative powers.
    degree: usize,
    contracted: Atom,
}

/// Whether a request for `loops`-loop integrals is evaluated by RustRed under
/// `settings` (the first backend in the evaluation order that can own such an
/// integral is RustRed). The input rewrites below are value-neutral but only
/// serve RustRed's matcher and reduction, so other backends keep the input.
pub(crate) fn rustred_evaluates_loop(settings: &VakintSettings, loops: usize) -> bool {
    for method in &settings.evaluation_order.0 {
        match method {
            crate::EvaluationMethod::RustRed(_) => return true,
            crate::EvaluationMethod::AlphaLoop(_) | crate::EvaluationMethod::MATAD(_)
                if loops <= 3 =>
            {
                return false;
            }
            crate::EvaluationMethod::FMFT(_) if loops == 4 => return false,
            crate::EvaluationMethod::PySecDec(_) => return false,
            _ => {}
        }
    }
    false
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
        vertices
            .iter()
            .position(|vertex| vertex.contains(&signed))?,
        vertices
            .iter()
            .position(|vertex| vertex.contains(&-signed))?,
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
            numerator *= (S.dot.call((momentum, momentum)) - &mass).pow(exponent);
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
    let placeholder =
        |axis: usize| Atom::var(vk_symbol!(format!("vakint_master_transport_{axis}")));
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

/// Largest numerator degree (in scalar products) among lowered keys, measured
/// from the positive `base` of the target form: a monomial `D^a` lowers to
/// `base - a`; `None` if a key lies above the base.
fn lowered_degree(base: &IntegralKey, keys: &[IntegralKey]) -> Option<usize> {
    keys.iter()
        .map(|key| {
            base.powers()
                .iter()
                .zip(key.powers())
                .try_fold(0usize, |sum, (base, power)| {
                    usize::try_from(base - power).ok().map(|step| sum + step)
                })
        })
        .try_fold(0usize, |degree, key| key.map(|key| degree.max(key)))
}

/// Memo key of a contracted scalar input (integral, numerator, mass).
fn memo_key(contracted: &Atom, numerator: &Atom, mass: &Atom) -> [u8; 32] {
    let mut hasher = Sha256::new();
    for atom in [contracted, numerator, mass] {
        hasher.update(atom.to_canonical_string().as_bytes());
        hasher.update([0u8]);
    }
    hasher.finalize().into()
}

impl SavedApplication {
    /// Remaining terminals of the package's own family.
    pub(super) fn own_terminals(&self) -> Result<BTreeSet<IntegralKey>, RustRedEvaluationError> {
        let prepared = self.prepared.lock().map_err(failure)?;
        Ok(prepared
            .remaining_terminals()
            .iter()
            .filter(|terminal| terminal.family_fingerprint() == prepared.family().fingerprint())
            .map(|terminal| terminal.integral().clone())
            .collect())
    }

    /// The package's own (common-basis) family.
    fn own_family(&self) -> Result<Arc<rustred::family::IntegralFamily>, RustRedEvaluationError> {
        let prepared = self.prepared.lock().map_err(failure)?;
        prepared
            .output_family_shared(prepared.family().fingerprint())
            .ok_or_else(|| failure("saved application has no output family"))
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
                    let profile =
                        crate::topologies::power_colored_matroid_profile(contracted.as_view(), 5)
                            .ok()
                            .flatten();
                    let degree = key
                        .powers()
                        .iter()
                        .filter(|power| **power < 0)
                        .map(|power| power.unsigned_abs() as usize)
                        .sum();
                    Some(MasterForm {
                        key: key.clone(),
                        parent,
                        canonical,
                        profile,
                        degree,
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
    /// inside the recorded entry predicates). Used by the tests; evaluation
    /// shares one lowering service per call ([`Self::witness_keys_with`]).
    #[cfg(test)]
    pub(super) fn witness_keys_unchecked(
        &self,
        settings: &VakintSettings,
        application: &SavedApplication,
        numerator: AtomView,
        matched: &ReplacementRules,
    ) -> Result<Vec<IntegralKey>, VakintError> {
        let family = application.own_family()?;
        let scalar = FamilyScalarNumeratorService::try_new(
            &family,
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(native_failure)?;
        self.witness_keys_with(settings, &scalar, numerator, matched)
    }

    /// The exact lowered keys of a scalar numerator under one witness, without
    /// the saved entry-union admission check, with a prepared lowering service.
    fn witness_keys_with(
        &self,
        settings: &VakintSettings,
        scalar: &FamilyScalarNumeratorService<'_>,
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
        let canonical =
            match crate::topologies::power_colored_canonical_graph(contracted.as_view(), 5) {
                Ok(Some(canonical)) => canonical,
                Ok(None) | Err(VakintError::InvalidIntegralFormat(_)) => {
                    outcome("no_canonical_form", 0, 0);
                    return Ok(None);
                }
                Err(error) => return Err(error),
            };
        let application = self.saved_application(settings)?;
        let forms = application.master_forms();
        // Masters on the same colored graph first; then masters on another
        // graph with the same colored line-matroid profile (2-isomorphic
        // embeddings), compared by matroid-level maps.
        let profile =
            match crate::topologies::power_colored_matroid_profile(contracted.as_view(), 5) {
                Ok(profile) => profile,
                Err(VakintError::InvalidIntegralFormat(_)) => None,
                Err(error) => return Err(error),
            };
        let candidates = forms
            .iter()
            .filter(|form| form.canonical == canonical)
            .map(|form| (form, false))
            .chain(
                forms
                    .iter()
                    .filter(|form| {
                        form.canonical != canonical && profile.is_some() && form.profile == profile
                    })
                    .map(|form| (form, true)),
            )
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            outcome("no_candidate", 0, 0);
            return Ok(None);
        }
        let memo = memo_key(&contracted, &contracted_numerator, &mass);
        if application
            .master_transport_misses
            .lock()
            .map_err(failure)?
            .contains(&memo)
        {
            outcome("not_a_master_memoized", candidates.len(), 0);
            return Ok(None);
        }
        let terminals = application.own_terminals()?;
        // One lowering service for every map of this call.
        let family = application.own_family()?;
        let scalar = FamilyScalarNumeratorService::try_new(
            &family,
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(native_failure)?;
        let mut tried = 0;
        for (form, matroid) in &candidates {
            let Some((target, base)) = parent_form(&form.key, form.parent, &mass) else {
                continue;
            };
            let Some(witness) = self
                .topologies
                .match_topologies_to_user_input(target.as_view(), settings.allow_unknown_integrals)?
                .filter(|witness| matches!(witness.canonical_topology, Topology::FiveLoop(_)))
            else {
                continue;
            };
            // The maps are offered lazily: the first map whose transported
            // numerator lowers onto remaining terminals only ends the search
            // (highly symmetric product graphs have many). Every map
            // preserves the numerator degree, so a first lowering of another
            // degree than the master's ends this candidate.
            let mut accepted = None;
            let mut first = true;
            let mut accept = |images: &[Atom]| -> Result<bool, VakintError> {
                tried += 1;
                let transported = transport_numerator(&contracted_numerator, images);
                let keys = match self.witness_keys_with(
                    settings,
                    &scalar,
                    transported.as_view(),
                    &witness,
                ) {
                    Ok(keys) => keys,
                    Err(error) if resource_error(&error) => return Err(error),
                    Err(_) => return Ok(false),
                };
                if !keys.is_empty() && keys.iter().all(|key| terminals.contains(key)) {
                    accepted = Some((transported, keys.len()));
                    return Ok(true);
                }
                Ok(std::mem::take(&mut first) && lowered_degree(&base, &keys) != Some(form.degree))
            };
            let found = if *matroid {
                crate::topologies::find_power_preserving_matroid_momentum_map(
                    contracted.as_view(),
                    form.contracted.as_view(),
                    5,
                    &mut accept,
                )
            } else {
                crate::topologies::find_power_preserving_momentum_map(
                    contracted.as_view(),
                    form.contracted.as_view(),
                    5,
                    &mut accept,
                )
            };
            match found {
                Ok(_) | Err(VakintError::InvalidIntegralFormat(_)) => {}
                Err(error) => return Err(error),
            }
            let Some((transported, keys)) = accepted else {
                continue;
            };
            log::debug!(target: "vakint::five_loop_match", "phase=master_transport_selected master={:?} keys={keys} matroid={matroid}", form.key.powers());
            *integral = target;
            *numerator = transported;
            outcome("master", candidates.len(), tried);
            return Ok(Some(witness));
        }
        application
            .master_transport_misses
            .lock()
            .map_err(failure)?
            .insert(memo);
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
        integral *= S.prop.call((
            id,
            S.edge.call((left, right)),
            momentum,
            Atom::one(),
            Atom::num(*power),
        ));
    }
    S.topo.call(integral)
}

/// The program's own long form of an integral (every physical slot, in slot
/// order, at mass `mass`), as its shipped test inputs write it.
fn four_loop_long_form(
    edges: &[(i64, i64)],
    momenta: &[Atom],
    powers: &[i64],
    mass: &Atom,
) -> Atom {
    let mut integral = Atom::one();
    for (slot, (((left, right), momentum), power)) in
        edges.iter().zip(momenta).zip(powers).enumerate()
    {
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

/// The denominators of a long form with positive powers: sorted (line
/// momentum up to sign as integer coefficients over `loops` loop momenta,
/// power), or `None` for another shape.
fn denominator_signature(integral: AtomView, loops: usize) -> Option<Vec<(Vec<i128>, i64)>> {
    let pattern = vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,pow_)")
        .unwrap()
        .to_pattern();
    let variables = (1..=loops).map(|axis| S.k.call(axis)).collect::<Vec<_>>();
    let mut lines = integral
        .pattern_match(&pattern, None, None)
        .map(|matched| {
            let power = crate::get_integer_from_atom(matched[&vk_symbol!("pow_")].as_view())
                .filter(|power| *power > 0)?;
            let mut vector =
                crate::topologies::integer_line_vector(&matched[&vk_symbol!("q_")], &variables)?;
            if vector
                .iter()
                .find(|value| **value != 0)
                .is_some_and(|first| *first < 0)
            {
                vector.iter_mut().for_each(|value| *value = -*value);
            }
            Some((vector, power))
        })
        .collect::<Option<Vec<_>>>()?;
    lines.sort_unstable();
    Some(lines)
}

impl Vakint {
    /// Rewrite a four-loop scalar integral that equals a shipped four-loop
    /// master (positive powers only, no numerator) as a master's own long
    /// form, so the ordinary matcher sees the program's slot order (input
    /// propagator ids in another order made matching take ~20-30 s and could
    /// select a symmetric image that the program reduces by a walk). Returns
    /// whether it did.
    ///
    /// The data ships some integrals as masters of several families. The
    /// rewrite is deterministic: (1) an input already in a program's slot
    /// order is untouched; (2) an input with exactly the denominators of a
    /// shipped master's program form (up to propagator ids, vertex labels,
    /// orientations and contracted zero lines) is written as that master: the
    /// identity map, the input's own family; (3) otherwise the input is
    /// written as the first shipped master in manifest order that equals it
    /// under a verified power-preserving momentum map of the line matroid:
    /// one fixed representative per equivalence class. Only when RustRed
    /// evaluates four-loop integrals under `settings`.
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
        if !rustred_evaluates_loop(settings, 4) {
            return Ok(false);
        }
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
            let pattern = vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,pow_)")
                .unwrap()
                .to_pattern();
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
            while let (Some((_, edges)), Some(momenta)) = (
                programs.physical_edges(index),
                programs.physical_momenta(index),
            ) {
                if lines.len() == momenta.len()
                    && lines.iter().zip(edges.iter().zip(momenta)).enumerate().all(
                        |(slot, ((id, left, right, q), ((l, r), momentum)))| {
                            *id == Some(slot as i64 + 1)
                                && *left == Some(*l)
                                && *right == Some(*r)
                                && q == momentum
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
                        // Scalar masters only: a master with a numerator
                        // (negative power) is not a scalar input's image.
                        if powers.len() >= momenta.len()
                            && powers[momenta.len()..].iter().all(|power| *power == 0)
                            && powers[..momenta.len()].iter().all(|power| *power >= 0)
                        {
                            let powers = powers[..momenta.len()].to_vec();
                            let contracted = four_loop_contracted(edges, momenta, &powers);
                            if let (Ok(Some(profile)), Some(lines)) = (
                                crate::topologies::power_colored_matroid_profile(
                                    contracted.as_view(),
                                    4,
                                ),
                                denominator_signature(contracted.as_view(), 4),
                            ) {
                                forms.push(super::super::four_loop::FourLoopMasterForm {
                                    program: index,
                                    powers,
                                    profile,
                                    lines,
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
        let rewrite = |integral: &mut Atom,
                       form: &super::super::four_loop::FourLoopMasterForm,
                       rule: &str| {
            let (Some((_, edges)), Some(momenta)) = (
                programs.physical_edges(form.program),
                programs.physical_momenta(form.program),
            ) else {
                return false;
            };
            *integral = four_loop_long_form(edges, momenta, &form.powers, &mass);
            log::debug!(target: "vakint::saved_application", "phase=four_loop_master_transport rule={rule} program={} powers={:?} elapsed_seconds={:.6}", form.program, form.powers, started.elapsed().as_secs_f64());
            true
        };
        // (2) The input's own family: a master with exactly its denominators.
        if let Some(lines) = denominator_signature(contracted.as_view(), 4)
            && let Some(form) = forms.iter().find(|form| form.lines == lines)
        {
            return Ok(rewrite(integral, form, "own_denominators"));
        }
        // (3) The class representative: the first equal master in manifest
        // order under a verified map of the line matroid.
        let profile =
            match crate::topologies::power_colored_matroid_profile(contracted.as_view(), 4) {
                Ok(Some(profile)) => profile,
                Ok(None) | Err(VakintError::InvalidIntegralFormat(_)) => return Ok(false),
                Err(error) => return Err(error),
            };
        for form in forms.iter().filter(|form| form.profile == profile) {
            match crate::topologies::find_power_preserving_matroid_momentum_map(
                contracted.as_view(),
                form.contracted.as_view(),
                4,
                |_| Ok(true),
            ) {
                Ok(Some(true)) => return Ok(rewrite(integral, form, "class_representative")),
                Ok(_) | Err(VakintError::InvalidIntegralFormat(_)) => continue,
                Err(error) => return Err(error),
            }
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

    /// The lines `(id, left, right, momentum, mass, power)` of a long form.
    fn lines(integral: &Atom) -> Vec<(i64, i64, i64, Atom, Atom, Atom)> {
        let pattern = vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,pow_)")
            .unwrap()
            .to_pattern();
        integral
            .pattern_match(&pattern, None, None)
            .map(|matched| {
                let integer = |name: &str| {
                    crate::get_integer_from_atom(matched[&vk_symbol!(name)].as_view()).unwrap()
                };
                (
                    integer("id_"),
                    integer("nl_"),
                    integer("nr_"),
                    matched[&vk_symbol!("q_")].clone(),
                    matched[&vk_symbol!("m_")].clone(),
                    matched[&vk_symbol!("pow_")].clone(),
                )
            })
            .collect()
    }

    fn long_form(lines: &[(i64, i64, i64, Atom, Atom, Atom)]) -> Atom {
        S.topo.call(
            lines
                .iter()
                .map(|(id, left, right, q, m, power)| {
                    S.prop
                        .call((*id, S.edge.call((*left, *right)), q, m, power))
                })
                .product::<Atom>(),
        )
    }

    /// A deterministic pseudo-random graph-generator spelling of a long form
    /// and its numerator: vertices relabelled, propagator ids permuted, line
    /// orientations flipped (momentum negated) and a unimodular loop
    /// relabelling applied to both.
    fn generator_form(integral: &Atom, numerator: &Atom, loops: usize, seed: u64) -> (Atom, Atom) {
        let mut state = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
        let mut next = move |modulus: usize| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 33) as usize) % modulus.max(1)
        };
        let mut images = (1..=loops).map(|axis| S.k.call(axis)).collect::<Vec<_>>();
        for _ in 0..8 {
            let (i, j) = (next(loops), next(loops));
            match next(3) {
                0 if i != j => {
                    let sign = if next(2) == 0 { 1 } else { -1 };
                    images[i] = &images[i] + Atom::num(sign) * &images[j];
                }
                1 => images.swap(i, j),
                _ => images[i] = -&images[i],
            }
        }
        let mut lines = lines(integral);
        let vertices = lines
            .iter()
            .flat_map(|line| [line.1, line.2])
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let mut labels = vertices.clone();
        for index in (1..labels.len()).rev() {
            labels.swap(index, next(index + 1));
        }
        let relabel = |vertex: i64| labels[vertices.iter().position(|v| *v == vertex).unwrap()];
        let mut ids = (1..=lines.len() as i64).collect::<Vec<_>>();
        for index in (1..ids.len()).rev() {
            ids.swap(index, next(index + 1));
        }
        for (line, id) in lines.iter_mut().zip(ids) {
            line.0 = id;
            let (left, right) = (relabel(line.1), relabel(line.2));
            if next(2) == 0 {
                (line.1, line.2) = (left, right);
            } else {
                (line.1, line.2) = (right, left);
                line.3 = -&line.3;
            }
        }
        (
            transport_numerator(&long_form(&lines), &images),
            transport_numerator(numerator, &images),
        )
    }

    /// The Whitney twists of a contracted long form (a 2-vertex cut {u, v} and
    /// a component C of the rest: the C-side lines exchange u and v and are
    /// reversed, momenta unchanged; the same denominators, a valid routing)
    /// whose power-colored incidence graph is not isomorphic to the input's.
    fn whitney_twists(integral: &Atom) -> Vec<Atom> {
        let lines = lines(integral);
        let vertices = lines
            .iter()
            .flat_map(|line| [line.1, line.2])
            .collect::<BTreeSet<_>>();
        let original = canonical(integral);
        let mut twists = Vec::new();
        for &u in &vertices {
            for &v in vertices.range(u + 1..) {
                let mut root = BTreeMap::new();
                for line in &lines {
                    if ![u, v].contains(&line.1) && ![u, v].contains(&line.2) {
                        union(&mut root, line.1, line.2);
                    }
                }
                let components = vertices
                    .iter()
                    .filter(|vertex| ![u, v].contains(vertex))
                    .map(|vertex| find(&mut root, *vertex))
                    .collect::<BTreeSet<_>>();
                if components.len() < 2 {
                    continue;
                }
                for component in components {
                    let mut twisted = lines.clone();
                    for line in &mut twisted {
                        let inside = |vertex: i64, root: &mut BTreeMap<i64, i64>| {
                            ![u, v].contains(&vertex) && find(root, vertex) == component
                        };
                        if inside(line.1, &mut root) || inside(line.2, &mut root) {
                            let swap = |vertex: i64| match vertex {
                                x if x == u => v,
                                x if x == v => u,
                                x => x,
                            };
                            (line.1, line.2) = (swap(line.2), swap(line.1));
                        }
                    }
                    let twisted = long_form(&twisted);
                    if canonical(&twisted) != original
                        && !twists
                            .iter()
                            .any(|other| canonical(other) == canonical(&twisted))
                    {
                        twists.push(twisted);
                    }
                }
            }
        }
        twists
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
        let other =
            IntegralKey::try_new(vec![1, 1, 0, 1, 3, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1]).unwrap();
        let other = contracted_form(&other, first_parent(&other).unwrap()).unwrap();
        assert_ne!(canonical(&other), canonical(&contracted));
        assert!(
            crate::topologies::power_preserving_momentum_maps(
                other.as_view(),
                contracted.as_view(),
                5
            )
            .unwrap()
            .unwrap()
            .is_empty()
        );
        // The pseudo-random graph-generator spelling is the same colored graph.
        let (generated, _) = generator_form(&contracted, &Atom::one(), 5, 7);
        assert_eq!(canonical(&generated), canonical(&contracted));
    }

    /// Package-free: the matroid-level maps recognise Whitney-twisted
    /// embeddings of masters (no graph isomorphism exists), every map they
    /// offer is verified, and they respect powers.
    #[test]
    fn matroid_maps_recognise_whitney_twists_and_respect_powers() {
        Vakint::initialize_vakint_symbols();
        for powers in [
            // BASIS-v2 P472, P463, P536, P426 (VSPEED review, out/C.json).
            [1, 1, 2, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1],
            [1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1],
            [2, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1],
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0],
        ] {
            let key = IntegralKey::try_new(powers.to_vec()).unwrap();
            let contracted = contracted_form(&key, first_parent(&key).unwrap()).unwrap();
            let twists = whitney_twists(&contracted);
            assert!(!twists.is_empty(), "{powers:?} has a non-isomorphic twist");
            let profile = |input: &Atom| {
                crate::topologies::power_colored_matroid_profile(input.as_view(), 5)
                    .unwrap()
                    .unwrap()
            };
            for twist in twists.iter().chain([&contracted]) {
                let (generated, _) = generator_form(twist, &Atom::one(), 5, 11);
                for input in [twist.clone(), generated] {
                    assert_eq!(profile(&input), profile(&contracted));
                    let mut maps = 0;
                    let found = crate::topologies::find_power_preserving_matroid_momentum_map(
                        input.as_view(),
                        contracted.as_view(),
                        5,
                        |images| {
                            maps += 1;
                            assert_eq!(images.len(), 5);
                            Ok(false)
                        },
                    )
                    .unwrap();
                    assert_eq!(found, Some(false));
                    assert!(maps > 0, "{powers:?}: no matroid map from {input}");
                    // The colored-graph search cannot see a twist.
                    if canonical(&input) != canonical(&contracted) {
                        assert_eq!(
                            crate::topologies::find_power_preserving_momentum_map(
                                input.as_view(),
                                contracted.as_view(),
                                5,
                                |_| Ok(true),
                            )
                            .unwrap(),
                            Some(false)
                        );
                    }
                }
            }
            // Another power on one line: no power-preserving map.
            let mut dotted = powers;
            let slot = dotted.iter().position(|power| *power == 1).unwrap();
            dotted[slot] = 2;
            let dotted = IntegralKey::try_new(dotted.to_vec()).unwrap();
            let dotted = contracted_form(&dotted, first_parent(&dotted).unwrap()).unwrap();
            assert_eq!(
                crate::topologies::find_power_preserving_matroid_momentum_map(
                    twists[0].as_view(),
                    dotted.as_view(),
                    5,
                    |_| Ok(true),
                )
                .unwrap(),
                Some(false)
            );
        }
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
                        .witness_keys_unchecked(
                            &settings,
                            &application,
                            numerator.as_view(),
                            &witness,
                        )
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
            let mut forms = vec![
                (integral.clone(), numerator.clone()),
                (contracted.clone(), contracted_numerator.clone()),
                (relabel(&contracted), relabel(&contracted_numerator)),
                (relabel(&integral), relabel(&numerator)),
            ];
            forms.extend(
                whitney_twists(&contracted)
                    .into_iter()
                    .map(|twist| (twist, contracted_numerator.clone())),
            );
            for (form, (integral, numerator)) in forms.iter().enumerate() {
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

    /// End to end, as an evaluation runs: every package master given in its
    /// parent form, contracted (a graph generator's form), in a pseudo-random
    /// graph-generator spelling (ids, vertex labels, orientations, loop
    /// momenta) and as a non-isomorphic Whitney twist of its contraction is
    /// reduced as the master itself, through `evaluate_rustred_symbolic` (Retain)
    /// and through the public `evaluate` (Numerical), with no native walk:
    /// either no native request (a request of remaining terminals only) or at
    /// most two processed keys. Every form gives the parent form's value; a
    /// non-product master's Retain value is that master alone.
    #[test]
    #[ignore = "loads the official five-loop package and evaluates every master (minutes)"]
    fn every_master_input_form_is_reduced_as_the_master_without_a_walk() {
        let root = super::super::super::test_data::external_path_through(5);
        let mut engine = Vakint::new().unwrap();
        engine.load_external_data(&root, 5).unwrap();
        let settings = |mode: crate::RustRedMasterMode| VakintSettings {
            form_exe_path: "/master-transport-test-must-not-use-form".into(),
            evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
                crate::RustRedEvaluationOptions { master_mode: mode },
            )]),
            number_of_terms_in_epsilon_expansion: 6,
            ..Default::default()
        };
        let retain = settings(crate::RustRedMasterMode::Retain);
        let numerical = settings(crate::RustRedMasterMode::Numerical);
        let application = engine.saved_application(&retain).unwrap();
        let (fingerprint, products) = {
            let prepared = application.prepared.lock().unwrap();
            let fingerprint = prepared.family().fingerprint().to_owned();
            let products = prepared
                .remaining_terminals()
                .iter()
                .filter(|terminal| {
                    terminal.family_fingerprint() == fingerprint
                        && (prepared.terminal_product(terminal).is_some()
                            || prepared.terminal_product_chart(terminal).is_some())
                })
                .map(|terminal| terminal.integral().clone())
                .collect::<BTreeSet<_>>();
            (fingerprint, products)
        };
        let terminals = application.own_terminals().unwrap();
        let masters = vk_parse!("RustRedMaster(args__)").unwrap().to_pattern();
        // One evaluation with a fresh observer: (value, native processed keys).
        // A public numerical result that needs master orders the package does
        // not ship fails closed naming the master: that error is the value.
        let evaluate = |public: bool, input: &Atom| -> (Result<Atom, String>, Option<usize>) {
            let mut callback = |_: &rustred::solver::CandidateRoutedApplicationEnvelope<'_>| {
                Ok(rustred::solver::CandidateRoutedObserverControl::Continue)
            };
            let mut observer = CandidateRoutedApplicationObserver::new(
                rustred::solver::CandidateRoutedObservationOptions::default(),
                &mut callback,
            );
            let value = if public {
                match engine.evaluate_observed(
                    &numerical,
                    input.as_view(),
                    Some(&mut observer),
                    None,
                ) {
                    Ok(value) => Ok(value),
                    Err(VakintError::MissingMasterOrders(missing)) => Err(format!("{missing:?}")),
                    Err(error) => panic!("public evaluation of {input}: {error}"),
                }
            } else {
                Ok(engine
                    .evaluate_rustred_symbolic_observed(
                        &retain,
                        input.as_view(),
                        Some(&mut observer),
                        None,
                    )
                    .unwrap()
                    .value)
            };
            (
                value,
                observer
                    .last_ordered_progress()
                    .map(|progress| progress.processed_keys),
            )
        };
        let mut checked = 0;
        for (index, key) in terminals.iter().enumerate() {
            let (integral, numerator) = parent_input(key).unwrap();
            let (contracted, contracted_numerator, _) =
                contract(&retain, integral.as_view(), numerator.as_view())
                    .unwrap()
                    .unwrap();
            let mut forms = vec![
                ("parent", integral.clone() * &numerator),
                ("contracted", contracted.clone() * &contracted_numerator),
            ];
            let (generated, generated_numerator) =
                generator_form(&contracted, &contracted_numerator, 5, index as u64 + 1);
            forms.push(("generator", generated * generated_numerator));
            if let Some(twist) = whitney_twists(&contracted).into_iter().next() {
                forms.push(("twist", twist * &contracted_numerator));
            }
            let mut expected: [Option<Result<Atom, String>>; 2] = [None, None];
            for (name, input) in &forms {
                for public in [false, true] {
                    let _ = engine.clear_rustred_application_witness_cache();
                    let (value, processed) = evaluate(public, input);
                    assert!(
                        processed.is_none_or(|processed| processed <= 2),
                        "master {:?} form {name} public={public}: {processed:?} processed keys",
                        key.powers()
                    );
                    let slot = &mut expected[usize::from(public)];
                    match slot {
                        None => {
                            if let Ok(value) = &value
                                && !public
                                && !products.contains(key)
                            {
                                let symbol = SavedApplication::symbolic_integral(
                                    &fingerprint,
                                    key,
                                    vk_symbol!("RustRedMaster"),
                                    &mut BTreeMap::new(),
                                );
                                assert!(
                                    value.contains(symbol.as_view())
                                        && value.pattern_match(&masters, None, None).count() == 1,
                                    "master {:?}: Retain value {value}",
                                    key.powers()
                                );
                            }
                            *slot = Some(value);
                        }
                        Some(expected) => assert_eq!(
                            &value,
                            expected,
                            "master {:?} form {name} public={public}",
                            key.powers()
                        ),
                    }
                }
            }
            checked += 1;
        }
        assert_eq!(checked, terminals.len());
    }

    /// Every shipped scalar four-loop master: its program form is left
    /// untouched; with reversed propagator ids, or contracted with relabelled
    /// vertices and flipped orientations (its own denominators) it is written
    /// as exactly its own program form; under loop relabellings it is written
    /// as one fixed shipped master (the class representative, the same for
    /// every relabelling) that equals it under a verified map; a dotted variant
    /// is rewritten only onto a shipped master equal to it under such a map.
    #[test]
    #[ignore = "loads the official four-loop data"]
    fn every_four_loop_master_form_is_written_in_its_program_form() {
        let engine = super::super::super::test_data::external_engine();
        let settings = VakintSettings {
            evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
                crate::RustRedEvaluationOptions::default(),
            )]),
            ..Default::default()
        };
        let data = engine.rustred_data.as_ref().unwrap();
        let programs = data.four_loop.as_ref().unwrap();
        let k = |axis: usize| S.k.call(axis);
        let relabel_a =
            |atom: &Atom| transport_numerator(atom, &[k(2), k(1) + k(3), -k(3), k(4) - k(1)]);
        let relabel_b =
            |atom: &Atom| transport_numerator(atom, &[-k(4), k(3), k(1) - k(2), k(2) + k(4)]);
        let one = Atom::one();
        let rewritten = |input: &Atom| {
            let mut integral = input.clone();
            let moved = engine
                .transport_four_loop_master(&settings, &mut integral, one.as_view())
                .unwrap();
            (moved, integral)
        };
        let shipped = |integral: &Atom| {
            programs.master_forms.get().unwrap().iter().any(|form| {
                let (_, edges) = programs.physical_edges(form.program).unwrap();
                let momenta = programs.physical_momenta(form.program).unwrap();
                four_loop_long_form(edges, momenta, &form.powers, &Atom::one()) == *integral
            })
        };
        // A verified power-preserving matroid-level map between two inputs.
        let equal = |left: &Atom, right: &Atom| {
            let contracted = |input: &Atom| {
                contract(&settings, input.as_view(), one.as_view())
                    .unwrap()
                    .unwrap()
                    .0
            };
            crate::topologies::find_power_preserving_matroid_momentum_map(
                contracted(left).as_view(),
                contracted(right).as_view(),
                4,
                |_| Ok(true),
            )
            .unwrap()
                == Some(true)
        };
        let mut checked = 0;
        for master in &data.loops[&4].manifest.masters {
            let alias = master["alias"].as_str().unwrap();
            let powers: Vec<i64> = serde_json::from_value(master["powers"].clone()).unwrap();
            let label = format!("four-loop-{}", alias.trim_start_matches("4L_"));
            let index = (0..4)
                .find(|index| programs.physical_edges(*index).unwrap().0 == label)
                .unwrap();
            let (_, edges) = programs.physical_edges(index).unwrap();
            let momenta = programs.physical_momenta(index).unwrap();
            let powers = &powers[..momenta.len()];
            let own = four_loop_long_form(edges, momenta, powers, &Atom::one());
            // The program form itself is left untouched.
            assert_eq!(
                rewritten(&own),
                (false, own.clone()),
                "{alias} {powers:?} program form"
            );
            if powers.iter().any(|power| *power < 0) {
                // Masters with numerators are not scalar images.
                continue;
            }
            let reversed = {
                let mut lines = lines(&own);
                let n = lines.len() as i64;
                for line in &mut lines {
                    line.0 = n + 1 - line.0;
                }
                long_form(&lines)
            };
            let (contracted_own, _, _) = contract(&settings, own.as_view(), one.as_view())
                .unwrap()
                .unwrap();
            let (generated, _) = {
                // Vertex labels, ids and orientations only: the own denominators.
                let mut lines = lines(&contracted_own);
                let vertices = lines
                    .iter()
                    .flat_map(|line| [line.1, line.2])
                    .max()
                    .unwrap();
                for (position, line) in lines.iter_mut().enumerate() {
                    (line.1, line.2) = (vertices + 1 - line.2, vertices + 1 - line.1);
                    if position % 2 == 0 {
                        line.3 = -&line.3;
                        (line.1, line.2) = (line.2, line.1);
                    }
                }
                (long_form(&lines), ())
            };
            for (form, input) in [("reversed", reversed), ("generator", generated)] {
                assert_eq!(
                    rewritten(&input),
                    (true, own.clone()),
                    "{alias} {powers:?} {form}"
                );
            }
            let (moved_a, representative) = rewritten(&relabel_a(&own));
            let (moved_b, other) = rewritten(&relabel_b(&own));
            assert!(moved_a && moved_b, "{alias} {powers:?} relabelled");
            assert_eq!(
                representative, other,
                "{alias} {powers:?}: one representative"
            );
            assert!(
                shipped(&representative),
                "{alias} {powers:?}: {representative}"
            );
            assert!(
                equal(&own, &representative),
                "{alias} {powers:?}: {representative}"
            );
            // A dotted variant is rewritten only onto a shipped master that
            // equals it (it can equal one only under a symmetry).
            let mut dotted = powers.to_vec();
            let first = dotted.iter().position(|power| *power > 0).unwrap();
            dotted[first] += 1;
            let input = four_loop_long_form(edges, momenta, &dotted, &Atom::one());
            let (moved, integral) = rewritten(&relabel_a(&input));
            if moved {
                assert!(
                    shipped(&integral) && equal(&input, &integral),
                    "{alias} dotted {dotted:?} rewritten onto {integral}"
                );
            }
            checked += 1;
        }
        assert!(checked > 0);
    }
}
