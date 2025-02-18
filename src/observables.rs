use crate::momentum::FourMomentum;
use crate::utils::{FloatLike, F};
use crate::Settings;
use itertools::Itertools;
#[allow(unused_imports)]
use libc::{c_double, c_int, c_void};
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use lorentz_vector::LorentzVector;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use spenso::complex::Complex;
use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum PhaseSpaceSelectorSettings {
    #[serde(rename = "jet")]
    Jet(JetSliceSettings),
    #[serde(rename = "ranged")]
    RangeFilter(RangeFilterSettings),
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum ObservableSettings {
    #[serde(rename = "jet1pt")]
    Jet1PT(Jet1PTSettings),
    #[serde(rename = "one_particle_obs")]
    SingleParticleObservable(SingleParticleObservableSettings),
    AFB(AFBSettings),
    #[serde(rename = "cross_section")]
    CrossSection,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct RangeFilterSettings {
    pub pdgs: Vec<isize>,
    pub filter: FilterQuantity,
    pub min_value: f64,
    pub max_value: f64,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct JetSliceSettings {
    pub min_jets: usize,
    pub max_jets: usize,
    pub min_j1pt: f64,
    pub max_j1pt: f64,
    pub dR: f64,
    pub min_jpt: f64,
    pub use_fastjet: bool,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct Jet1PTSettings {
    pub x_min: f64,
    pub x_max: f64,
    pub dR: f64,
    pub min_jpt: f64,
    pub n_bins: usize,
    pub write_to_file: bool,
    pub filename: String,
    pub use_fastjet: bool,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct AFBSettings {
    pub x_min: f64,
    pub x_max: f64,
    pub n_bins: usize,
    pub write_to_file: bool,
    pub filename: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct SingleParticleObservableSettings {
    pub x_min: f64,
    pub x_max: f64,
    pub n_bins: usize,
    #[serde(default = "default_true")]
    pub write_to_file: bool,
    #[serde(default = "default_false")]
    pub log_obs: bool,
    #[serde(default = "default_false")]
    pub log_x_axis: bool,
    #[serde(default = "default_true")]
    pub log_y_axis: bool,
    pub filename: String,
    pub pdgs: Vec<isize>,
    pub quantity: FilterQuantity,
}

fn default_true() -> bool {
    true
}

fn default_false() -> bool {
    false
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum FilterQuantity {
    #[serde(rename = "E")]
    Energy,
    #[serde(rename = "CosTheta")]
    CosThetaP,
    #[serde(rename = "PT")]
    PT,
}

impl fmt::Display for FilterQuantity {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FilterQuantity::Energy => {
                write!(f, "Energy")
            }
            FilterQuantity::CosThetaP => {
                write!(f, "CosTheta")
            }
            FilterQuantity::PT => {
                write!(f, "Pt")
            }
        }
    }
}

#[cfg(feature = "fjcore")]
mod fjcore {
    use libc::{c_double, c_int, c_void};

    #[link(name = "fjcore", kind = "static")]
    extern "C" {
        pub fn fastjet_workspace() -> *mut c_void;
        //pub fn fastjet_free(workspace: *mut c_void);
        pub fn fastjetppgenkt_(
            workspace: *mut c_void,
            p: *const c_double,
            npart: *const c_int,
            R: *const c_double,
            ptjet_min: *const c_double,
            palg: *const c_double,
            jets: *mut c_double,
            njets: *mut c_int,
            whichjets: *mut c_int,
        );
    }
}

#[derive(Debug, Default, Clone)]
pub struct EventInfo {
    pub accepted_event_counter: usize,
    pub rejected_event_counter: usize,
    pub no_phase_space_counter: usize,
    pub zero_eval_counter: usize,
}

#[derive(Debug, Default, Clone)]
pub struct AverageAndErrorAccumulator {
    sum: f64,
    sum_sq: f64,
    weight_sum: f64,
    avg_sum: f64,
    pub avg: f64,
    pub err: f64,
    guess: f64,
    chi_sq: f64,
    chi_sum: f64,
    chi_sq_sum: f64,
    num_samples: usize,
    cur_iter: usize,
}

impl AverageAndErrorAccumulator {
    pub fn add_sample(&mut self, sample: f64) {
        self.sum += sample;
        self.sum_sq += sample * sample;
        self.num_samples += 1;
    }

    pub fn merge_samples(&mut self, other: &mut AverageAndErrorAccumulator) {
        self.sum += other.sum;
        self.sum_sq += other.sum_sq;
        self.num_samples += other.num_samples;

        // reset the other
        other.sum = 0.;
        other.sum_sq = 0.;
        other.num_samples = 0;
    }

    pub fn update_iter(&mut self) {
        // TODO: we could be throwing away events that are very rare
        if self.num_samples < 2 {
            self.cur_iter += 1;
            return;
        }

        let n = self.num_samples as f64;
        self.sum /= n;
        self.sum_sq /= n * n;
        let mut w = (self.sum_sq * n).sqrt();

        w = ((w + self.sum) * (w - self.sum)) / (n - 1.);
        if w == 0. {
            w = f64::EPSILON;
        }
        w = 1. / w;

        self.weight_sum += w;
        self.avg_sum += w * self.sum;
        let sigsq = 1. / self.weight_sum;
        self.avg = sigsq * self.avg_sum;
        self.err = sigsq.sqrt();
        if self.cur_iter == 0 {
            self.guess = self.sum;
        }
        w *= self.sum - self.guess;
        self.chi_sum += w;
        self.chi_sq_sum += w * self.sum;
        self.chi_sq = self.chi_sq_sum - self.avg * self.chi_sum;

        // reset
        self.sum = 0.;
        self.sum_sq = 0.;
        self.num_samples = 0;
        self.cur_iter += 1;
    }
}

#[derive(Default)]
pub struct EventManager {
    pub event_selector: Vec<Selectors>,
    pub observables: Vec<Observables>,
    pub event_buffer: Vec<Event>,
    pub track_events: bool,
    pub accepted_event_counter: usize,
    pub rejected_event_counter: usize,
    pub no_phase_space_counter: usize,
    pub zero_eval_counter: usize,
    pub time_integrand_evaluation: bool,
    pub integrand_evaluation_timing: u128,
    pub integrand_evaluation_timing_start: Option<std::time::Instant>,
    pub event_group_counter: usize,
}

#[allow(unused)]
impl EventManager {
    pub fn new(track_events: bool, settings: Settings) -> EventManager {
        let mut observables = vec![];
        for o in &settings.observables {
            match o {
                ObservableSettings::Jet1PT(settings) => {
                    if track_events {
                        observables.push(Observables::Jet1PT(Jet1PTObservable::new(
                            settings.x_min,
                            settings.x_max,
                            settings.n_bins,
                            settings.dR,
                            settings.min_jpt,
                            settings.write_to_file,
                            settings.filename.clone(),
                            settings.use_fastjet,
                        )));
                    }
                }
                ObservableSettings::CrossSection => {
                    observables.push(Observables::CrossSection(CrossSectionObservable::default()));
                }
                ObservableSettings::AFB(settings) => {
                    if track_events {
                        observables.push(Observables::AFB(AFBObservable::new(
                            settings.x_min,
                            settings.x_max,
                            settings.n_bins,
                            settings.write_to_file,
                            settings.filename.clone(),
                        )));
                    }
                }
                ObservableSettings::SingleParticleObservable(settings) => {
                    if track_events {
                        observables.push(Observables::SingleParticle(
                            SingleParticleObservable::new(
                                settings.x_min,
                                settings.x_max,
                                settings.n_bins,
                                settings.write_to_file,
                                settings.filename.clone(),
                                settings.quantity,
                                settings.pdgs.clone(),
                                settings.log_obs,
                                settings.log_x_axis,
                                settings.log_y_axis,
                            ),
                        ));
                    }
                }
            }
        }

        let mut selectors = vec![];
        for s in &settings.selectors {
            match s {
                PhaseSpaceSelectorSettings::RangeFilter(settings) => {
                    selectors.push(Selectors::RangeFilter(RangedSelector::new(settings)));
                }
                PhaseSpaceSelectorSettings::Jet(settings) => {
                    selectors.push(Selectors::Jet(JetSelector::new(settings)));
                }
            }
        }

        EventManager {
            event_selector: selectors,
            observables,
            event_buffer: Vec::with_capacity(100),
            track_events,
            accepted_event_counter: 0,
            rejected_event_counter: 0,
            no_phase_space_counter: 0,
            zero_eval_counter: 0,
            time_integrand_evaluation: false,
            integrand_evaluation_timing: 0,
            integrand_evaluation_timing_start: None,
            event_group_counter: 0,
        }
    }

    pub fn create_event<T: FloatLike>(
        &self,
        orig_incoming_momenta: Vec<FourMomentum<F<T>>>,
        cut_momenta: Vec<FourMomentum<F<T>>>,
    ) -> Event {
        let mut incoming_momenta = SmallVec::new();
        for p in orig_incoming_momenta.iter() {
            incoming_momenta.push(p.to_f64())
        }
        let mut outgoing_momenta: SmallVec<[FourMomentum<F<f64>>; 4]> = SmallVec::new();
        for p in cut_momenta.iter() {
            outgoing_momenta.push(p.to_f64())
        }
        let final_state_particle_ids = SmallVec::from(vec![0; outgoing_momenta.len()]);
        Event {
            kinematic_configuration: (incoming_momenta, outgoing_momenta),
            final_state_particle_ids,
            integrand: Complex::new(F(1.), F(0.)),
            weights: SmallVec::from_slice(&[F(0.)]),
        }
    }

    pub fn pass_selection(&mut self, event: &mut Event) -> bool {
        for f in &mut self.event_selector {
            if !f.process_event(event) {
                return false;
            }
        }
        true
    }

    pub fn add_event(&mut self, mut event: Event) -> bool {
        if !self.track_events && self.event_selector.is_empty() {
            return true;
        }

        // run the event through the selectors
        for f in &mut self.event_selector {
            if !f.process_event(&mut event) {
                if self.track_events {
                    self.rejected_event_counter += 1;
                }
                return false;
            }
        }

        if self.track_events {
            self.event_buffer.push(event);
            self.accepted_event_counter += 1;
        }
        true
    }

    pub fn merge_samples(&mut self, other: &mut EventManager) {
        for (o1, o2) in self
            .observables
            .iter_mut()
            .zip(other.observables.iter_mut())
        {
            o1.merge_samples(o2);
        }

        self.accepted_event_counter += other.accepted_event_counter;
        self.rejected_event_counter += other.rejected_event_counter;
        self.no_phase_space_counter += other.no_phase_space_counter;
        self.zero_eval_counter += other.zero_eval_counter;
        self.event_group_counter += other.event_group_counter;
        self.integrand_evaluation_timing = other.integrand_evaluation_timing;
        other.accepted_event_counter = 0;
        other.rejected_event_counter = 0;
        other.event_group_counter = 0;
        other.no_phase_space_counter = 0;
        other.zero_eval_counter = 0;
        other.integrand_evaluation_timing = 0;
    }

    pub fn update_result(&mut self, iter: usize) {
        for o in &mut self.observables {
            o.update_result(self.event_group_counter, iter);
        }
    }

    pub fn update_live_result(&mut self) {
        for o in &mut self.observables {
            // for now, only live update the cross section
            #[allow(clippy::single_match)]
            match o {
                Observables::CrossSection(c) => c.update_live_result(),
                _ => {}
            }
        }
    }

    pub fn clear(&mut self, count_as_rejected: bool) {
        self.accepted_event_counter -= self.event_buffer.len();
        if count_as_rejected {
            self.rejected_event_counter += self.event_buffer.len();
        }
        self.event_buffer.clear();
    }

    pub fn process_events(&mut self, integrand_result: Complex<f64>, integrand_weight: f64) {
        if self.track_events {
            // give the events to an observable function
            for o in &mut self.observables {
                o.process_event_group(&self.event_buffer, integrand_weight);
            }
            self.event_buffer.clear();
        } else {
            for o in &mut self.observables {
                if let Observables::CrossSection(c) = o {
                    c.add_sample(integrand_result, integrand_weight);
                    break;
                }
            }
        }

        self.event_group_counter += 1;
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Event {
    #[allow(clippy::type_complexity)]
    pub kinematic_configuration: (
        SmallVec<[FourMomentum<F<f64>>; 2]>,
        SmallVec<[FourMomentum<F<f64>>; 4]>,
    ),
    pub final_state_particle_ids: SmallVec<[isize; 5]>,
    pub integrand: Complex<F<f64>>,
    pub weights: SmallVec<[F<f64>; 1]>,
}

#[allow(unused)]
pub enum Selectors {
    All(NoEventSelector),
    Jet(JetSelector),
    RangeFilter(RangedSelector),
}

#[allow(unused)]
impl Selectors {
    #[inline]
    fn process_event(&mut self, event: &mut Event) -> bool {
        match self {
            Selectors::All(f) => f.process_event(event),
            Selectors::Jet(f) => f.process_event(event),
            Selectors::RangeFilter(f) => f.process_event(event),
        }
    }
}

pub trait EventSelector {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event(&mut self, event: &mut Event) -> bool;
}

#[derive(Default)]
pub struct NoEventSelector {}

impl EventSelector for NoEventSelector {
    #[inline]
    fn process_event(&mut self, _event: &mut Event) -> bool {
        true
    }
}

#[derive(Debug, Clone)]
pub struct RangedSelector {
    pdgs: Vec<isize>,
    filter: FilterQuantity,
    min_value: f64,
    max_value: f64,
}

impl EventSelector for RangedSelector {
    #[inline]
    fn process_event(&mut self, event: &mut Event) -> bool {
        for (pdg, mom) in event
            .final_state_particle_ids
            .iter()
            .zip(&event.kinematic_configuration.1)
        {
            if self.pdgs.contains(pdg) {
                let value = match self.filter {
                    FilterQuantity::Energy => mom.temporal.value,
                    FilterQuantity::CosThetaP => mom.spatial.pz / mom.spatial.norm(),
                    FilterQuantity::PT => mom.pt(),
                };

                if value < F(self.min_value) || value > F(self.max_value) {
                    return false;
                }
            }
        }
        true
    }
}

#[allow(dead_code)]
impl RangedSelector {
    pub fn new(s: &RangeFilterSettings) -> RangedSelector {
        RangedSelector {
            pdgs: s.pdgs.clone(),
            filter: s.filter,
            min_value: s.min_value,
            max_value: s.max_value,
        }
    }
}

#[derive(Debug, Clone)]
pub struct JetClustering {
    //ordered_jets: Vec<LorentzVector<f64>>,
    ordered_pt: Vec<f64>,
    //jet_structure: Vec<usize>, // map from original momenta to jets
    d_r: f64,
    min_jpt: f64,
    use_fastjet: bool,
    fastjet_jets_in: Vec<f64>,
    fastjet_jets_out: Vec<f64>,
    fastjet_jets_map: Vec<c_int>,
    #[cfg(feature = "fjcore")]
    fastjet_workspace: *mut c_void,
}

unsafe impl std::marker::Send for JetClustering {}

#[allow(dead_code)]
impl JetClustering {
    pub fn new(use_fastjet: bool, d_r: f64, min_jpt: f64) -> JetClustering {
        #[cfg(feature = "fjcore")]
        let fastjet_workspace = unsafe { fjcore::fastjet_workspace() }; // TODO: free

        JetClustering {
            //ordered_jets: vec![],
            //jet_structure: vec![],
            ordered_pt: vec![],
            fastjet_jets_in: vec![],
            fastjet_jets_out: vec![],
            fastjet_jets_map: vec![],
            d_r,
            min_jpt,
            #[cfg(feature = "fjcore")]
            fastjet_workspace,
            use_fastjet,
        }
    }

    pub fn cluster_fastjet(&mut self, event: &Event) {
        self.fastjet_jets_in.clear();

        #[allow(unused_variables)]
        let mut len: c_int = 0;
        for (e, id) in event
            .kinematic_configuration
            .1
            .iter()
            .zip_eq(&event.final_state_particle_ids)
        {
            // filter for jet particles: u, d, c, s, d, g, QCD_ghost
            //TODO make it a hyperparam of the observable!
            if id.abs() < 6
                || *id == 21
                || id.abs() == 82
                || (id.abs() >= 3370 && id.abs() < 3380)
                || id.abs() == 0
            {
                self.fastjet_jets_in.extend(e.into_iter().map(|x| x.0));
                len += 1;
            }
        }

        self.fastjet_jets_out.clear();
        self.fastjet_jets_out.resize(self.fastjet_jets_in.len(), 0.);

        self.fastjet_jets_map.clear();
        self.fastjet_jets_map.resize(self.fastjet_jets_in.len(), 0);

        #[allow(unused_mut)]
        let mut actual_len: c_int = 0;
        #[allow(unused_variables)]
        let palg = -1.0;
        #[allow(unused_variables)]
        let clustering_ptjet_min = 0.;
        #[cfg(feature = "fjcore")]
        {
            if len > 0 {
                unsafe {
                    fjcore::fastjetppgenkt_(
                        self.fastjet_workspace,
                        &self.fastjet_jets_in[0] as *const c_double,
                        &len as *const c_int,
                        &self.d_r as *const c_double,
                        &clustering_ptjet_min as *const c_double,
                        &palg as *const c_double,
                        &mut self.fastjet_jets_out[0] as *mut c_double,
                        &mut actual_len as *mut c_int,
                        &mut self.fastjet_jets_map[0] as *mut c_int,
                    );
                }
            }
        }

        self.ordered_pt.clear();
        for i in 0..actual_len as usize {
            let jet = LorentzVector::from_slice(&self.fastjet_jets_out[i * 4..(i + 1) * 4]);
            self.ordered_pt.push(jet.pt());
        }

        // Filter order_pt jets by removing those below the necessary min_jpt
        let min_pt = self.min_jpt;
        self.ordered_pt.retain(|&pt| pt >= min_pt);
    }

    pub fn cluster(&mut self, event: &Event) {
        if self.use_fastjet {
            self.cluster_fastjet(event)
        } else {
            self.cluster_custom(event)
        }
    }

    pub fn cluster_custom(&mut self, event: &Event) {
        //self.ordered_jets.clear();
        self.ordered_pt.clear();
        //self.jet_structure.clear();

        // group jets based on their dR
        let ext = &event.kinematic_configuration.1;

        let index_offset = 2;
        if ext.len() == 5 {
            let pt1 = ext[index_offset].pt();
            let pt2 = ext[index_offset + 1].pt();
            let pt3 = ext[index_offset + 2].pt();
            let dr12 = ext[index_offset].delta_r(&ext[index_offset + 1]);
            let dr13 = ext[index_offset].delta_r(&ext[index_offset + 2]);
            let dr23 = ext[index_offset + 1].delta_r(&ext[index_offset + 2]);

            // we have at least 2 jets
            if dr12 < F(self.d_r) {
                let m12 = ext[index_offset] + ext[index_offset + 1];
                let pt12 = m12.pt();
                self.ordered_pt.push(pt12.0);
                self.ordered_pt.push(pt3.0);
            } else if dr13 < F(self.d_r) {
                let m13 = ext[index_offset] + ext[index_offset + 2];
                let pt13 = m13.pt();
                self.ordered_pt.push(pt13.0);
                self.ordered_pt.push(pt2.0);
            } else if dr23 < F(self.d_r) {
                let m23 = ext[index_offset + 1] + ext[index_offset + 2];
                let pt23 = m23.pt();
                self.ordered_pt.push(pt23.0);
                self.ordered_pt.push(pt1.0);
            } else {
                self.ordered_pt.push(pt1.0);
                self.ordered_pt.push(pt2.0);
                self.ordered_pt.push(pt3.0);
            }

            // sort from large pt to small
            self.ordered_pt
                .sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
        } else {
            // treat every external momentum as its own jet for now
            for m in ext {
                self.ordered_pt.push(m.pt().0);
            }
            self.ordered_pt
                .sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
        }

        // Filter order_pt jets by removing those below the necessary min_jpt
        let min_pt = self.min_jpt;
        self.ordered_pt.retain(|&pt| pt >= min_pt);
    }
}

pub struct JetSelector {
    jet_selector_settings: JetSliceSettings,
    clustering: JetClustering,
}

impl JetSelector {
    pub fn new(settings: &JetSliceSettings) -> JetSelector {
        JetSelector {
            jet_selector_settings: settings.clone(),
            clustering: JetClustering::new(settings.use_fastjet, settings.dR, settings.min_jpt),
        }
    }
}

impl EventSelector for JetSelector {
    #[inline]
    fn process_event(&mut self, event: &mut Event) -> bool {
        self.clustering.cluster(event);
        //info!("Event: {:#?}",event);
        //info!("clustering: {:#?}",self.clustering);
        // if self.clustering.ordered_pt[0] > 30. {
        //     info!("Event: {:#?}", event);
        //     info!("clustering: {:#?}", self.clustering);
        //     panic!("Ordered pt too large!: {:?}", self.clustering.ordered_pt);
        // }
        self.clustering.ordered_pt.len() >= self.jet_selector_settings.min_jets
            && self.clustering.ordered_pt.len() <= self.jet_selector_settings.max_jets
            && (self.clustering.ordered_pt.is_empty()
                || (self.clustering.ordered_pt[0] >= self.jet_selector_settings.min_j1pt
                    && (self.jet_selector_settings.max_j1pt < 0.0
                        || (self.clustering.ordered_pt[0] <= self.jet_selector_settings.max_j1pt))))
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum Observables {
    CrossSection(CrossSectionObservable),
    Jet1PT(Jet1PTObservable),
    AFB(AFBObservable),
    SingleParticle(SingleParticleObservable),
}

#[allow(dead_code)]
impl Observables {
    #[inline]
    pub fn process_event_group(&mut self, event: &[Event], integrator_weight: f64) {
        use self::Observables::*;
        match self {
            CrossSection(o) => o.process_event_group(event, integrator_weight),
            Jet1PT(o) => o.process_event_group(event, integrator_weight),
            AFB(o) => o.process_event_group(event, integrator_weight),
            SingleParticle(o) => o.process_event_group(event, integrator_weight),
        }
    }

    #[inline]
    pub fn merge_samples(&mut self, other: &mut Observables) {
        use self::Observables::*;
        match (self, other) {
            (CrossSection(o1), CrossSection(o2)) => o1.merge_samples(o2),
            (Jet1PT(o1), Jet1PT(o2)) => o1.merge_samples(o2),
            (AFB(o1), AFB(o2)) => o1.merge_samples(o2),
            (SingleParticle(o1), SingleParticle(o2)) => o1.merge_samples(o2),
            (o1, o2) => panic!(
                "Cannot merge observables of different types: {:?} vs {:?}",
                o1, o2
            ),
        }
    }

    /// Produce the result (histogram, etc.) of the observable from all processed event groups.
    #[inline]
    pub fn update_result(&mut self, total_events: usize, iter: usize) {
        use self::Observables::*;
        match self {
            CrossSection(o) => o.update_result(total_events, iter),
            Jet1PT(o) => o.update_result(total_events, iter),
            AFB(o) => o.update_result(total_events, iter),
            SingleParticle(o) => o.update_result(total_events, iter),
        }
    }
}

pub trait Observable {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event_group(&mut self, event: &[Event], integrator_weight: f64);

    fn merge_samples(&mut self, other: &mut Self)
    where
        Self: Sized;

    /// Produce the result (histogram, etc.) of the observable from all processed event groups.
    fn update_result(&mut self, total_events: usize, iter: usize);
}

#[derive(Debug, Clone, Default)]
pub struct CrossSectionObservable {
    pub re: AverageAndErrorAccumulator,
    pub im: AverageAndErrorAccumulator,
}

impl CrossSectionObservable {
    pub fn add_sample(&mut self, integrand: Complex<f64>, integrator_weight: f64) {
        self.re.add_sample(integrand.re * integrator_weight);
        self.im.add_sample(integrand.im * integrator_weight);
    }

    /// Give a live update on a copy of the statistics
    pub fn update_live_result(&self) {
        let mut re = self.re.clone();
        re.update_iter();
        let mut im = self.im.clone();
        im.update_iter();
    }
}

impl Observable for CrossSectionObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        let mut integrand = Complex::<F<f64>>::new_zero();
        for e in events {
            integrand += e.integrand;
        }

        self.re.add_sample(integrand.re.0 * integrator_weight);
        self.im.add_sample(integrand.im.0 * integrator_weight);
    }

    fn merge_samples(&mut self, other: &mut CrossSectionObservable) {
        self.re.merge_samples(&mut other.re);
        self.im.merge_samples(&mut other.im);
    }

    fn update_result(&mut self, _total_events: usize, _iter: usize) {
        self.re.update_iter();
        self.im.update_iter();
    }
}

#[derive(Debug, Clone)]
pub struct Jet1PTObservable {
    x_min: f64,
    x_max: f64,
    jet_clustering: JetClustering,
    index_event_accumulator: Vec<(isize, f64)>,
    bins: Vec<AverageAndErrorAccumulator>,
    write_to_file: bool,
    filename: String,
    num_events: usize,
}

impl Jet1PTObservable {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        x_min: f64,
        x_max: f64,
        num_bins: usize,
        d_r: f64,
        min_jpt: f64,
        write_to_file: bool,
        filename: String,
        use_fastjet: bool,
    ) -> Jet1PTObservable {
        Jet1PTObservable {
            x_min,
            x_max,
            index_event_accumulator: Vec::with_capacity(20),
            bins: vec![AverageAndErrorAccumulator::default(); num_bins],
            jet_clustering: JetClustering::new(use_fastjet, d_r, min_jpt),
            write_to_file,
            filename,
            num_events: 0,
        }
    }
}

impl Observable for Jet1PTObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        // add events in a correlated manner such that cancellations are realized in the error
        self.index_event_accumulator.clear();
        for e in events {
            self.jet_clustering.cluster_fastjet(e);
            let max_pt = self.jet_clustering.ordered_pt[0];

            let index = ((max_pt - self.x_min) / (self.x_max - self.x_min) * self.bins.len() as f64)
                as isize;
            if index >= 0 && index < self.bins.len() as isize {
                let mut new = true;
                for i in self.index_event_accumulator.iter_mut() {
                    if i.0 == index {
                        *i = (index, i.1 + e.integrand.re.0 * integrator_weight);
                        new = false;
                        break;
                    }
                }
                if new {
                    self.index_event_accumulator
                        .push((index, e.integrand.re.0 * integrator_weight));
                }
            }
        }

        for i in &self.index_event_accumulator {
            self.bins[i.0 as usize].add_sample(i.1);
        }
    }

    fn merge_samples(&mut self, other: &mut Jet1PTObservable) {
        // TODO: check if bins are compatible?
        for (b, bo) in self.bins.iter_mut().zip_eq(&mut other.bins) {
            b.merge_samples(bo);
        }
    }

    fn update_result(&mut self, total_events: usize, iter: usize) {
        let diff = total_events - self.num_events;
        self.num_events = total_events;

        // rescale the entries in each bin by the number of samples
        // per bin over the complete number of events (including rejected) ones
        // this is not the same as recaling the average after recombining
        // with previous iterations
        for b in &mut self.bins {
            let scaling = b.num_samples as f64 / diff as f64;
            b.sum *= scaling;
            b.sum_sq *= scaling * scaling;
            b.update_iter();
        }

        if self.write_to_file {
            let file_name = {
                let split = self.filename.split('.').collect::<Vec<_>>();
                if split.len() == 1 {
                    format!("{}_iter_{}", self.filename, iter)
                } else {
                    format!(
                        "{}_iter_{}.{}",
                        split.as_slice()[..split.len() - 1].join("."),
                        iter,
                        split.last().unwrap()
                    )
                }
            };
            let mut f = BufWriter::new(File::create(file_name).expect("Could not create log file"));

            writeln!(f, "##& xmin & xmax & central value & dy &\n").unwrap();
            writeln!(
                f,
                "<histogram> {} \"j1 pT |X_AXIS@LIN |Y_AXIS@LOG |TYPE@AL\"",
                self.bins.len()
            )
            .unwrap();

            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64
                    + self.x_min;

                writeln!(
                    f,
                    "  {:.8e}   {:.8e}   {:.8e}   {:.8e}",
                    c1, c2, b.avg, b.err,
                )
                .unwrap();
            }

            writeln!(f, "<\\histogram>").unwrap();
        } else {
            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64;
                info!("{}={}: {} +/ {}", c1, c2, b.avg, b.err,);
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct AFBObservable {
    x_min: f64,
    x_max: f64,
    event_accumulator: Vec<(isize, f64)>,
    bins: Vec<AverageAndErrorAccumulator>,
    write_to_file: bool,
    filename: String,
    num_events: usize,
}

impl AFBObservable {
    pub fn new(
        x_min: f64,
        x_max: f64,
        num_bins: usize,
        write_to_file: bool,
        filename: String,
    ) -> AFBObservable {
        AFBObservable {
            x_min,
            x_max,
            event_accumulator: Vec::with_capacity(20),
            bins: vec![AverageAndErrorAccumulator::default(); num_bins],
            write_to_file,
            filename,
            num_events: 0,
        }
    }
}

impl Observable for AFBObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        // add events in a correlated manner such that cancellations are realized in the error
        self.event_accumulator.clear();
        for e in events {
            let pin = e.kinematic_configuration.0[1];

            let pout = e
                .final_state_particle_ids
                .iter()
                .zip(&e.kinematic_configuration.1)
                .find(|(id, _mom)| 0 < **id && **id <= 6)
                .map(|(_, mom)| *mom)
                .unwrap();

            let costheta =
                pin.spatial * (pout.spatial) / (pin.spatial.norm() * pout.spatial.norm());

            let index = ((costheta.0 - self.x_min) / (self.x_max - self.x_min)
                * self.bins.len() as f64) as isize;

            if index >= 0 && index < self.bins.len() as isize {
                let mut new = true;
                for i in self.event_accumulator.iter_mut() {
                    if i.0 == index {
                        *i = (index, i.1 + e.integrand.re.0 * integrator_weight);
                        new = false;
                        break;
                    }
                }
                if new {
                    self.event_accumulator
                        .push((index, e.integrand.re.0 * integrator_weight));
                }
            }
        }

        for i in &self.event_accumulator {
            self.bins[i.0 as usize].add_sample(i.1);
        }
    }

    fn merge_samples(&mut self, other: &mut AFBObservable) {
        // TODO: check if bins are compatible?
        for (b, bo) in self.bins.iter_mut().zip_eq(&mut other.bins) {
            b.merge_samples(bo);
        }
    }

    fn update_result(&mut self, total_events: usize, iter: usize) {
        let diff = total_events - self.num_events;
        self.num_events = total_events;

        // rescale the entries in each bin by the number of samples
        // per bin over the complete number of events (including rejected) ones
        // this is not the same as recaling the average after recombining
        // with previous iterations
        for b in &mut self.bins {
            let scaling = b.num_samples as f64 / diff as f64;
            b.sum *= scaling;
            b.sum_sq *= scaling * scaling;
            b.update_iter();
        }

        if self.write_to_file {
            let file_name = {
                let split = self.filename.split('.').collect::<Vec<_>>();
                if split.len() == 1 {
                    format!("{}_iter_{}", self.filename, iter)
                } else {
                    format!(
                        "{}_iter_{}.{}",
                        split.as_slice()[..split.len() - 1].join("."),
                        iter,
                        split.last().unwrap()
                    )
                }
            };
            let mut f =
                BufWriter::new(File::create(file_name).expect("Could not create .HwU file"));
            writeln!(f, "##& xmin & xmax & central value & dy &\n").unwrap();
            writeln!(
                f,
                "<histogram> {} \"Angular distribution |X_AXIS@LIN |Y_AXIS@LOG |TYPE@AL\"",
                self.bins.len()
            )
            .unwrap();

            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64
                    + self.x_min;

                writeln!(
                    f,
                    "  {:.8e}   {:.8e}   {:.8e}   {:.8e}",
                    c1, c2, b.avg, b.err,
                )
                .unwrap();
            }

            writeln!(f, "<\\histogram>").unwrap();
        } else {
            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64;
                info!("{}={}: {} +/ {}", c1, c2, b.avg, b.err,);
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct SingleParticleObservable {
    x_min: f64,
    x_max: f64,
    event_accumulator: Vec<(isize, f64)>,
    bins: Vec<AverageAndErrorAccumulator>,
    write_to_file: bool,
    filename: String,
    num_events: usize,
    quantity: FilterQuantity,
    pdgs: Vec<isize>,
    log_obs: bool,
    log_x_axis: bool,
    log_y_axis: bool,
}

impl SingleParticleObservable {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        x_min: f64,
        x_max: f64,
        num_bins: usize,
        write_to_file: bool,
        filename: String,
        quantity: FilterQuantity,
        pdgs: Vec<isize>,
        log_obs: bool,
        log_x_axis: bool,
        log_y_axis: bool,
    ) -> SingleParticleObservable {
        SingleParticleObservable {
            x_min,
            x_max,
            event_accumulator: Vec::with_capacity(20),
            bins: vec![AverageAndErrorAccumulator::default(); num_bins],
            write_to_file,
            filename,
            num_events: 0,
            quantity,
            pdgs,
            log_obs,
            log_x_axis,
            log_y_axis,
        }
    }
}

impl Observable for SingleParticleObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        // add events in a correlated manner such that cancellations are realized in the error
        self.event_accumulator.clear();
        for e in events {
            let pin = e.kinematic_configuration.0[1];

            let pout_iter = e
                .final_state_particle_ids
                .iter()
                .zip(&e.kinematic_configuration.1)
                .filter(|(id, _mom)| self.pdgs.contains(*id))
                .map(|(_, mom)| *mom);

            for pout in pout_iter {
                let mut obs_evaluated = match self.quantity {
                    FilterQuantity::CosThetaP => {
                        pin.spatial * pout.spatial / (pin.spatial.norm() * pout.spatial.norm())
                    }
                    FilterQuantity::Energy => pout.temporal.value,
                    FilterQuantity::PT => pout.pt(),
                };
                if self.log_obs {
                    obs_evaluated = obs_evaluated.log10();
                }
                let index = ((obs_evaluated.0 - self.x_min) / (self.x_max - self.x_min)
                    * self.bins.len() as f64) as isize;

                if index >= 0 && index < self.bins.len() as isize {
                    let mut new = true;
                    for i in self.event_accumulator.iter_mut() {
                        if i.0 == index {
                            *i = (index, i.1 + e.integrand.re.0 * integrator_weight);
                            new = false;
                            break;
                        }
                    }
                    if new {
                        self.event_accumulator
                            .push((index, e.integrand.re.0 * integrator_weight));
                    }
                }
            }
        }

        for i in &self.event_accumulator {
            self.bins[i.0 as usize].add_sample(i.1);
        }
    }

    fn merge_samples(&mut self, other: &mut SingleParticleObservable) {
        // TODO: check if bins are compatible?
        for (b, bo) in self.bins.iter_mut().zip_eq(&mut other.bins) {
            b.merge_samples(bo);
        }
    }

    fn update_result(&mut self, total_events: usize, iter: usize) {
        let diff = total_events - self.num_events;
        self.num_events = total_events;

        // rescale the entries in each bin by the number of samples
        // per bin over the complete number of events (including rejected) ones
        // this is not the same as recaling the average after recombining
        // with previous iterations
        for b in &mut self.bins {
            let scaling = b.num_samples as f64 / diff as f64;
            b.sum *= scaling;
            b.sum_sq *= scaling * scaling;
            b.update_iter();
        }

        //let title = format!("{}({})", self.quantity, self.pdgs.iter().map(|&pdg| pdg.to_string()).join(",") );
        let title = Path::new(&self.filename)
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();
        let x_axis_mode = if self.log_x_axis { "LOG" } else { "LIN" };
        let y_axis_mode = if self.log_y_axis { "LOG" } else { "LIN" };
        if self.write_to_file {
            let file_name = {
                let split = self.filename.split('.').collect::<Vec<_>>();
                if split.len() == 1 {
                    format!("{}_iter_{}", self.filename, iter)
                } else {
                    format!(
                        "{}_iter_{}.{}",
                        split.as_slice()[..split.len() - 1].join("."),
                        iter,
                        split.last().unwrap()
                    )
                }
            };
            let mut f = BufWriter::new(File::create(file_name).expect("Could not create HwU file"));
            writeln!(f, "##& xmin & xmax & central value & dy &\n").unwrap();
            writeln!(
                f,
                "<histogram> {} \"{} |X_AXIS@{} |Y_AXIS@{} |TYPE@AL\"",
                self.bins.len(),
                title,
                x_axis_mode,
                y_axis_mode
            )
            .unwrap();

            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64
                    + self.x_min;

                writeln!(
                    f,
                    "  {:.8e}   {:.8e}   {:.8e}   {:.8e}",
                    c1, c2, b.avg, b.err,
                )
                .unwrap();
            }

            writeln!(f, "<\\histogram>").unwrap();
        } else {
            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64;
                info!("{}={}: {} +/ {}", c1, c2, b.avg, b.err,);
            }
        }
    }
}
