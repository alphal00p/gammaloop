//! Persistent, portable insertion requests for a given parent pole depth.
use super::*;
use crate::uv::delta_store;

/// An excised occurrence in its original labelled boundary frame and the exact
/// persistent-store keys shared by isomorphic occurrences at this depth.
pub struct RTableRequest {
    pub bare_key: [u8; 32],
    pub renormalized_key: [u8; 32],
    pub graph_sha: [u8; 32],
    pub loops: usize,
    pub depth: usize,
    diagram: FeynmanDiagram,
}

impl RTableRequest {
    /// Read the exact accepted rows, using the store's normal tier and integrity checks.
    pub fn persistent_value_hashes(&self) -> Result<[[u8; 32]; 2]> {
        let (_, store) = delta_store::installed()
            .ok_or_else(|| eyre!("insertion table requires a persistent store"))?;
        let read = |key: &[u8; 32]| match store.get(key) {
            delta_store::StoreGet::Hit(row) => Ok(row.value_sha),
            delta_store::StoreGet::Miss => Err(eyre!(
                "no accepted persistent row for insertion key {key:02x?}"
            )),
            delta_store::StoreGet::Error(error) => {
                Err(eyre!("cannot verify persisted insertion: {error}"))
            }
        };
        Ok([read(&self.bare_key)?, read(&self.renormalized_key)?])
    }
}

impl PropagatorInsertionEvaluator {
    pub(super) fn selected_r_insertions(
        &self,
        diagram: &FeynmanDiagram,
    ) -> Result<Vec<PropagatorInsertion>> {
        let rustred = matches!(
            self.settings
                .vakint
                .true_settings()
                .evaluation_order
                .0
                .first(),
            Some(vakint::EvaluationMethod::RustRed(_))
        );
        select_insertions_through(diagram, if rustred { 5 } else { 4 })
    }

    pub(super) fn truncate_insertion(value: &Atom, top: i64) -> Result<Atom> {
        Ok(value.series(GS.dim_epsilon, 0, top)?.to_atom().expand())
    }

    /// Inventory the same disjoint typed insertions as the R-subtraction route.
    /// Requires an installed store context, including the actual master data.
    pub fn r_table_requests(
        &self,
        parent: &FeynmanDiagram,
        depth: usize,
    ) -> Result<Vec<RTableRequest>> {
        ensure!(
            depth >= parent.loop_count() && depth <= 6,
            "table depth must cover the parent and be at most six"
        );
        self.selected_r_insertions(parent)?
            .into_iter()
            .map(|sigma| {
                let diagram = parent.excise_with_mapping(&sigma.support, &[])?.diagram;
                let graph = diagram.to_gamma_loop_graph(None, true)?;
                let support = graph
                    .no_dummy()
                    .subtract(&graph.external_filter::<SuBitGraph>());
                let certificate =
                    CountertermIdentity::for_insertion(&graph, &support, diagram.model())?;
                let (model, bytes) = certificate
                    .key
                    .portable_bytes()
                    .ok_or_else(|| eyre!("non-portable insertion certificate"))?;
                let context = delta_store::context(&self.renormalized_key(depth))?;
                let (renormalized_key, graph_sha, _) =
                    delta_store::key_sha("renormalized-sigma", &model, &context.sha, &bytes);
                let context = delta_store::context(&self.bare_key(depth, false))?;
                let (bare_key, _, _) =
                    delta_store::key_sha("bare-sigma", &model, &context.sha, &bytes);
                Ok(RTableRequest {
                    bare_key,
                    renormalized_key,
                    graph_sha,
                    loops: sigma.loops,
                    depth,
                    diagram,
                })
            })
            .collect()
    }

    /// Compute/read both rows with the production cache and boundary transport.
    /// An optional deeper run checks every retained coefficient independently.
    pub fn fill_r_table_request(
        &self,
        req: &RTableRequest,
        vakint: &Vakint,
        verify_depth: bool,
    ) -> Result<()> {
        ensure!(
            self.counterterm_cache.is_some(),
            "insertion table preparation requires a counterterm cache"
        );
        let check = |value: &Atom, top: i64| -> Result<()> {
            // The public guard expects its caller to project to the requested
            // order. UNKNOWN terms solely in the unused tail are admissible.
            let retained = Self::truncate_insertion(value, top)?;
            let unknown = vakint.unknown_master_coefficients(
                &self.settings.vakint.true_settings(),
                retained.as_view(),
                Some(top),
            )?;
            ensure!(
                unknown.is_empty(),
                "unknown master coefficients in the retained insertion through epsilon^{top}: {unknown:?}"
            );
            Ok(())
        };
        let bare = self.integrate_cached(&req.diagram, vakint, req.depth)?;
        let coefficients = self.insertion_coefficients_at(&req.diagram, vakint, req.depth)?;
        let top = req.depth as i64 - req.loops as i64 - 1;
        check(&bare, top)?;
        for (j, e) in coefficients.iter().enumerate() {
            check(e, top + j as i64)?;
        }
        if verify_depth {
            let deeper = self.integrate_cached(&req.diagram, vakint, req.depth + 2)?;
            check(&deeper, top)?;
            ensure!(
                insertion_coefficients_agree(
                    &Self::truncate_insertion(&bare, top)?,
                    &Self::truncate_insertion(&deeper, top)?
                )?
                .1,
                "bare insertion changes at depth N+2"
            );
            let deeper = self.insertion_coefficients_at(&req.diagram, vakint, req.depth + 2)?;
            ensure!(
                deeper.len() == coefficients.len(),
                "renormalized basis changes at depth N+2"
            );
            for (j, (a, b)) in coefficients.iter().zip(&deeper).enumerate() {
                let top = top + j as i64;
                check(b, top)?;
                ensure!(
                    insertion_coefficients_agree(
                        &Self::truncate_insertion(a, top)?,
                        &Self::truncate_insertion(b, top)?
                    )?
                    .1,
                    "renormalized coefficient {j} changes at depth N+2"
                );
            }
        }
        Ok(())
    }
}
