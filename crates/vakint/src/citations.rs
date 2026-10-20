use std::sync::atomic::{AtomicU8, Ordering};

#[cfg(feature = "symbolica_community_module")]
use symbolica::api::python::Citation;

pub(crate) static USED_CITATIONS: CitationUsage = CitationUsage(AtomicU8::new(0));

#[derive(Clone, Copy)]
#[repr(u8)]
pub(crate) enum CitationSource {
    Vakint = 16,
    Form = 1,
    Matad = 2,
    Fmft = 4,
    PySecDec = 8,
    RustRed = 32,
}

/// Monotonic, process-wide usage; reporting never clears another caller's flags.
#[derive(Default)]
pub(crate) struct CitationUsage(AtomicU8);

impl CitationUsage {
    #[inline]
    pub(crate) fn record(&self, source: CitationSource) {
        let bit = source as u8;
        // After first use this is only a relaxed read. No data is published by
        // these flags, so no acquire/release synchronization is required.
        if self.0.load(Ordering::Relaxed) & bit == 0 {
            self.0.fetch_or(bit, Ordering::Relaxed);
        }
    }

    #[cfg(feature = "symbolica_community_module")]
    pub(crate) fn get_citations(&self) -> Vec<Citation> {
        let used = self.0.load(Ordering::Relaxed);
        if used == 0 {
            return Vec::new();
        }
        // Authorship: docs/products/registry.toml. Backend papers:
        // docs/products/vakint/content/evaluation.typ, "Methods and software to cite".
        let mut citations = vec![Citation {
            id: "https://github.com/alphal00p/vakint#vakint".into(),
            url: "https://github.com/alphal00p/vakint".into(),
            reference: "Lucien Huber, Valentin Hirschi. Vakint (2026).".into(),
            bibtex: r#"@software{vakint,
  author = {Lucien Huber and Valentin Hirschi},
  title = {Vakint},
  year = {2026},
  url = {https://github.com/alphal00p/vakint}
}"#
            .into(),
            reasons: vec!["Provides the Vakint functionality in this community module.".into()],
            description: "".into(),
            relevance: None,
        }];
        for source in [
            CitationSource::Form,
            CitationSource::Matad,
            CitationSource::Fmft,
            CitationSource::PySecDec,
            CitationSource::RustRed,
        ] {
            if used & source as u8 != 0 {
                citations.push(source.citation());
            }
        }
        citations
    }
}

#[cfg(feature = "symbolica_community_module")]
impl CitationSource {
    fn citation(self) -> Citation {
        match self {
            Self::Vakint => unreachable!("the package credit is emitted separately"),
            Self::Form => Citation {
                id: "arXiv:1203.6543".into(),
                url: "https://arxiv.org/abs/1203.6543".into(),
                reference: "J. Kuipers and T. Ueda and J. A. M. Vermaseren and J. Vollinga. FORM version 4.0 (2012).".into(),
                bibtex: r#"@article{vakint_form,
  author = {J. Kuipers and T. Ueda and J. A. M. Vermaseren and J. Vollinga},
  title = {{FORM version 4.0}},
  year = {2012},
  eprint = {1203.6543},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Vakint used FORM for symbolic evaluation or tensor reduction.".into()],
                description: "".into(),
                relevance: None,
            },
            Self::Matad => Citation {
                id: "arXiv:hep-ph/0009029".into(),
                url: "https://arxiv.org/abs/hep-ph/0009029".into(),
                reference: "M. Steinhauser. MATAD: a program package for the computation of MAssive TADpoles (2000).".into(),
                bibtex: r#"@article{vakint_matad,
  author = {M. Steinhauser},
  title = {{MATAD: a program package for the computation of MAssive TADpoles}},
  year = {2000},
  eprint = {hep-ph/0009029},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Vakint used the MATAD massive-tadpole backend.".into()],
                description: "".into(),
                relevance: None,
            },
            Self::Fmft => Citation {
                id: "arXiv:1707.01710".into(),
                url: "https://arxiv.org/abs/1707.01710".into(),
                reference: "Andrey Pikelner. FMFT: Fully Massive Four-loop Tadpoles (2017).".into(),
                bibtex: r#"@article{vakint_fmft,
  author = {Andrey Pikelner},
  title = {{FMFT: Fully Massive Four-loop Tadpoles}},
  year = {2017},
  eprint = {1707.01710},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Vakint used the FMFT four-loop backend or its master-basis evaluations.".into()],
                description: "".into(),
                relevance: None,
            },
            Self::RustRed => Citation {
                id: "https://github.com/alphal00p/rustred".into(),
                url: "https://github.com/alphal00p/rustred".into(),
                reference: "RustRed: native integration-by-parts reduction.".into(),
                bibtex: "@software{rustred, title={RustRed}, url={https://github.com/alphal00p/rustred}}".into(),
                reasons: vec!["Vakint used native RustRed reduction with its shipped, checked artifacts.".into()],
                description: "Four-loop candidate evaluation is checked pointwise, not a full-family closure certificate.".into(),
                relevance: None,
            },
            Self::PySecDec => Citation {
                id: "arXiv:1703.09692".into(),
                url: "https://arxiv.org/abs/1703.09692".into(),
                reference: "S. Borowka and G. Heinrich and S. Jahn and S. P. Jones and M. Kerner and J. Schlenk and T. Zirke. pySecDec: a toolbox for the numerical evaluation of multi-scale integrals (2017).".into(),
                bibtex: r#"@article{vakint_pysecdec,
  author = {S. Borowka and G. Heinrich and S. Jahn and S. P. Jones and M. Kerner and J. Schlenk and T. Zirke},
  title = {{pySecDec: a toolbox for the numerical evaluation of multi-scale integrals}},
  year = {2017},
  eprint = {1703.09692},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Vakint used pySecDec numerical sector decomposition, including reused output.".into()],
                description: "".into(),
                relevance: None,
            },
        }
    }
}

#[cfg(all(test, feature = "symbolica_community_module"))]
mod tests {
    use super::*;

    #[test]
    fn native_catalog_usage_does_not_claim_form_execution() {
        let usage = CitationUsage::default();
        usage.record(CitationSource::RustRed);
        usage.record(CitationSource::Fmft);
        let ids: Vec<_> = usage.get_citations().into_iter().map(|c| c.id).collect();
        assert!(
            ids.iter()
                .any(|id| id == "https://github.com/alphal00p/rustred")
        );
        assert!(ids.iter().any(|id| id == "arXiv:1707.01710"));
        assert!(!ids.iter().any(|id| id == "arXiv:1203.6543"));
    }

    #[test]
    fn backend_citations_are_lazy_cumulative_and_deduplicated() {
        let usage = CitationUsage::default();
        assert!(usage.get_citations().is_empty());
        usage.record(CitationSource::PySecDec);
        usage.record(CitationSource::PySecDec);
        let citations = usage.get_citations();
        assert_eq!(citations.len(), 2);
        assert_eq!(citations[1].id, "arXiv:1703.09692");
        assert!(!citations[1].bibtex.is_empty());
        usage.record(CitationSource::Form);
        usage.record(CitationSource::Matad);
        usage.record(CitationSource::Fmft);
        assert_eq!(usage.get_citations().len(), 5);
        assert_eq!(usage.get_citations().len(), 5);
    }

    #[test]
    fn concurrent_backend_use_preserves_every_citation() {
        let usage = CitationUsage::default();
        std::thread::scope(|scope| {
            for source in [
                CitationSource::Form,
                CitationSource::Matad,
                CitationSource::Fmft,
                CitationSource::PySecDec,
            ] {
                let usage = &usage;
                scope.spawn(move || {
                    for _ in 0..1000 {
                        usage.record(source);
                    }
                });
            }
        });
        let ids: Vec<_> = usage.get_citations().into_iter().map(|c| c.id).collect();
        assert_eq!(
            &ids[1..],
            &[
                "arXiv:1203.6543",
                "arXiv:hep-ph/0009029",
                "arXiv:1707.01710",
                "arXiv:1703.09692"
            ]
        );
    }
}
