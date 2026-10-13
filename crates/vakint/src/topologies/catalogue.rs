//! Immutable optional matching data, prepared at the saved-artifact boundary.

use std::{sync::Arc, time::Instant};
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::{Topologies, Topology, five_loop};
use crate::{
    Integral, VakintError, get_integer_from_atom,
    symbols::S,
    utils::vakint_macros::{vk_parse, vk_symbol},
};

/// One-time catalogue preparation, separate from artifact loading and reduction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiveLoopCatalogueStatistics {
    pub topology_count: usize,
    pub preparation_seconds: f64,
}

#[derive(Debug)]
pub(crate) struct FiveLoopCatalogue {
    topologies: Vec<Topology>,
    statistics: FiveLoopCatalogueStatistics,
}

impl FiveLoopCatalogue {
    pub(crate) fn prepare() -> Result<Arc<Self>, VakintError> {
        let started = Instant::now();
        log::debug!(target: "vakint::saved_application", "phase=five_loop_catalogue_begin");
        let topologies = five_loop::generate()?;
        let statistics = FiveLoopCatalogueStatistics {
            topology_count: topologies.len(),
            preparation_seconds: started.elapsed().as_secs_f64(),
        };
        log::debug!(target: "vakint::saved_application", "phase=five_loop_catalogue_end topology_count={} elapsed_seconds={:.6}", statistics.topology_count, statistics.preparation_seconds);
        Ok(Arc::new(Self {
            topologies,
            statistics,
        }))
    }
}

impl Topologies {
    /// This only decides whether artifact loading is relevant, never entry
    /// admission. Reuse the ordinary Unknown matcher's validated loop census
    /// for long forms, including five-loop forms outside the saved catalogue.
    pub(crate) fn needs_five_loop_catalogue(input: AtomView) -> Result<bool, VakintError> {
        if five_loop::is_known_short_form(input) {
            return Ok(true);
        }
        // The legacy Unknown matcher assumes at least one loop momentum.
        // Avoid invoking it for constant rows or unrelated short-form heads.
        let momentum_pattern = vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,power_)")
            .unwrap()
            .to_pattern();
        if !input
            .pattern_match(&momentum_pattern, None, None)
            .any(|matched| matched[&vk_symbol!("q_")].contains_symbol(S.k))
        {
            return Ok(false);
        }
        let Some(matched) =
            Topology::Unknown(Integral::default()).match_topology_to_user_input(input)?
        else {
            return Ok(false);
        };
        let loops = matched
            .canonical_expression_substitutions
            .get(&vk_parse!("n_loops").unwrap())
            .and_then(|value| get_integer_from_atom(value.as_view()))
            .ok_or_else(|| {
                VakintError::InvalidIntegralFormat(
                    "native unknown-topology match omitted its loop count".into(),
                )
            })?;
        Ok(loops == 5)
    }

    /// Observe preparation without creating a catalogue or opening an artifact.
    pub fn five_loop_catalogue_statistics(&self) -> Option<FiveLoopCatalogueStatistics> {
        self.1.get().map(|catalogue| catalogue.statistics)
    }

    pub(crate) fn five_loop_catalogue(&self) -> Option<&Arc<FiveLoopCatalogue>> {
        self.1.get()
    }

    pub(crate) fn install_five_loop_catalogue(&self, catalogue: Arc<FiveLoopCatalogue>) {
        // The compiled catalogue is independent of application/cache options.
        // Publication is immutable and matching never holds the loader mutex.
        let _ = self.1.set(catalogue);
    }

    pub(super) fn iter(&self) -> impl Iterator<Item = &Topology> {
        // Preserve custom/base order and keep the unknown fallback last.
        self.0
            .iter()
            .filter(|t| !matches!(t, Topology::Unknown(_)))
            .chain(self.1.get().into_iter().flat_map(|c| c.topologies.iter()))
            .chain(self.0.iter().filter(|t| matches!(t, Topology::Unknown(_))))
    }

    /// Independent native census tests do not load reduction publications.
    #[cfg(test)]
    pub(crate) fn prepare_five_loop_catalogue_for_tests(&self) -> Result<(), VakintError> {
        if self.1.get().is_none() {
            self.install_five_loop_catalogue(FiveLoopCatalogue::prepare()?);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
