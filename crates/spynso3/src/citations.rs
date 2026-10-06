//! Cumulative attribution for successful operations in the community binding.

use std::sync::atomic::{AtomicU16, Ordering};

use idenso::tensor::AlgebraSettings;
use symbolica::api::python::Citation;

pub(crate) static USED: CitationUsage = CitationUsage(AtomicU16::new(0));

#[derive(Clone, Copy)]
#[repr(u16)]
pub(crate) enum Usage {
    Tensor = 1,
    Network = 2,
    Evaluation = 4,
    Contraction = 8,
    Algebra = 16,
    Gamma = 32,
    Color = 64,
    Epsilon = 128,
    Canonicalization = 256,
    DiracAdjoint = 512,
    Notation = 1024,
}

impl Usage {
    pub(crate) fn record(self) {
        USED.record(self);
    }

    fn reason(self) -> &'static str {
        match self {
            Self::Tensor => "Symbolic tensor expressions and representations.",
            Self::Network => "Tensor-network contractions.",
            Self::Evaluation => "Numerical evaluation of tensor components.",
            Self::Contraction => "Symbolic tensor contractions.",
            Self::Algebra => "Tensor algebra simplification.",
            Self::Gamma => "Dirac gamma algebra.",
            Self::Color => "Color algebra.",
            Self::Epsilon => "Levi-Civita identities.",
            Self::Canonicalization => "Canonical tensor expressions and dummy indices.",
            Self::DiracAdjoint => "Dirac adjoints.",
            Self::Notation => "Tensor index and contraction notation.",
        }
    }
}

/// Flags are monotonic: reporting neither resets usage nor loses concurrent calls.
#[derive(Default)]
pub(crate) struct CitationUsage(AtomicU16);

impl CitationUsage {
    fn record(&self, usage: Usage) {
        let bit = usage as u16;
        if self.0.load(Ordering::Relaxed) & bit == 0 {
            self.0.fetch_or(bit, Ordering::Relaxed);
        }
    }

    pub(crate) fn record_algebra(&self, settings: &AlgebraSettings) {
        self.record(Usage::Algebra);
        if settings.gamma.is_some() {
            self.record(Usage::Gamma);
        }
        if settings.color.is_some() {
            self.record(Usage::Color);
        }
        if settings.epsilon {
            self.record(Usage::Epsilon);
        }
    }

    pub(crate) fn get_citations(&self) -> Vec<Citation> {
        let used = self.0.load(Ordering::Relaxed);
        if used == 0 {
            return Vec::new();
        }
        let reasons = |operations: &[Usage]| {
            operations
                .iter()
                .copied()
                .filter(|usage| used & *usage as u16 != 0)
                .map(|usage| usage.reason().to_owned())
                .collect::<Vec<_>>()
        };
        // Authorship and software DOIs are maintained in docs/products/registry.toml.
        let mut spenso_reasons = reasons(&[Usage::Tensor, Usage::Network, Usage::Evaluation]);
        let idenso_reasons = reasons(&[
            Usage::Contraction,
            Usage::Algebra,
            Usage::Gamma,
            Usage::Color,
            Usage::Epsilon,
            Usage::Canonicalization,
            Usage::DiracAdjoint,
            Usage::Notation,
        ]);
        if !idenso_reasons.is_empty() {
            spenso_reasons.push("Tensor-expression structure and display for Idenso.".into());
        }
        let mut citations = vec![Citation {
            id: "10.5281/zenodo.18248388".into(),
            reference: "Lucien Huber. Spenso (2026).".into(),
            bibtex: r#"@software{spenso,
  author = {Lucien Huber},
  title = {Spenso},
  year = {2026},
  url = {https://github.com/alphal00p/spenso},
  doi = {10.5281/zenodo.18248388}
}"#
            .into(),
            reasons: spenso_reasons,
            description:
                "Symbolic tensors, tensor networks, and mathematical display in Symbolica.".into(),
            relevance: None,
        }];
        if !idenso_reasons.is_empty() {
            citations.push(Citation {
                id: "10.5281/zenodo.18248409".into(),
                reference: "Lucien Huber, Ben Ruijl. Idenso (2026).".into(),
                bibtex: r#"@software{idenso,
  author = {Lucien Huber and Ben Ruijl},
  title = {Idenso},
  year = {2026},
  url = {https://github.com/alphal00p/spenso},
  doi = {10.5281/zenodo.18248409}
}"#.into(),
                reasons: idenso_reasons,
                description: "Symbolic tensor algebra, including Dirac matrices, color factors, and Levi-Civita identities.".into(),
                relevance: None,
            });
        }
        if used & Usage::Gamma as u16 != 0 {
            citations.extend([
                Citation {
                    id: "arXiv:2601.19982".into(),
                    reference: "J. Davies, T. Kaneko, C. Marinissen, T. Ueda, J. A. M. Vermaseren. FORM Version 5.0 (2026).".into(),
                    bibtex: r#"@article{form_5,
  author = {J. Davies and T. Kaneko and C. Marinissen and T. Ueda and J. A. M. Vermaseren},
  title = {{FORM Version 5.0}},
  year = {2026},
  doi = {10.48550/arXiv.2601.19982},
  eprint = {2601.19982},
  archivePrefix = {arXiv},
  primaryClass = {hep-ph}
}"#.into(),
                    reasons: reasons(&[Usage::Gamma]),
                    description: "Basis for the Dirac gamma algebra rules in Idenso.".into(),
                    relevance: None,
                },
                Citation {
                    id: "arXiv:1707.06453".into(),
                    reference: "B. Ruijl, T. Ueda, J. A. M. Vermaseren. FORM version 4.2 (2017).".into(),
                    bibtex: r#"@article{form_4_2,
  author = {Ben Ruijl and Takahiro Ueda and Jos Vermaseren},
  title = {{FORM version 4.2}},
  year = {2017},
  doi = {10.48550/arXiv.1707.06453},
  eprint = {1707.06453},
  archivePrefix = {arXiv},
  primaryClass = {hep-ph}
}"#.into(),
                    reasons: reasons(&[Usage::Gamma]),
                    description: "Basis for the Dirac gamma algebra rules in Idenso.".into(),
                    relevance: None,
                },
            ]);
        }
        if used & Usage::Color as u16 != 0 {
            citations.push(Citation {
                id: "arXiv:hep-ph/9802376".into(),
                reference: "T. van Ritbergen, A. N. Schellekens, J. A. M. Vermaseren. Group theory factors for Feynman diagrams (1999).".into(),
                bibtex: r#"@article{vanRitbergen:1998pn,
  author = {T. van Ritbergen and A. N. Schellekens and J. A. M. Vermaseren},
  title = {Group theory factors for {Feynman} diagrams},
  journal = {International Journal of Modern Physics A},
  volume = {14},
  pages = {41--96},
  year = {1999},
  doi = {10.1142/S0217751X99000038},
  eprint = {hep-ph/9802376},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: reasons(&[Usage::Color]),
                description: "Basis for the color algebra rules in Idenso, through color.h.".into(),
                relevance: None,
            });
        }
        citations
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn citations_accumulate_reasons_without_conflating_packages() {
        let usage = CitationUsage::default();
        assert!(usage.get_citations().is_empty());
        usage.record(Usage::Tensor);
        usage.record(Usage::Tensor);
        let citations = usage.get_citations();
        assert_eq!(citations.len(), 1);
        assert_eq!(citations[0].reasons.len(), 1);

        usage.record_algebra(&AlgebraSettings {
            gamma: None,
            color: Some(Default::default()),
            epsilon: false,
            ..Default::default()
        });
        let citations = usage.get_citations();
        assert_eq!(citations.len(), 3);
        assert_eq!(citations[2].id, "arXiv:hep-ph/9802376");
        assert_eq!(
            citations[1].reasons,
            [Usage::Algebra.reason(), Usage::Color.reason()]
        );
        for citation in &citations {
            let identifier = citation.id.strip_prefix("arXiv:").unwrap_or(&citation.id);
            assert!(citation.bibtex.contains(identifier));
            assert!(!citation.description.is_empty());
            assert_eq!(citation.to_bibtex(), citation.bibtex);
            assert!(citation.to_markdown(false).contains(&citation.reasons[0]));
        }
        assert_eq!(usage.get_citations()[1].reasons, citations[1].reasons);

        usage.record(Usage::Gamma);
        usage.record(Usage::Color);
        usage.record(Usage::Gamma);
        let citations = usage.get_citations();
        assert_eq!(citations.len(), 5);
        assert_eq!(citations[2].id, "arXiv:2601.19982");
        assert_eq!(citations[2].reasons, [Usage::Gamma.reason()]);
        assert_eq!(citations[3].id, "arXiv:1707.06453");
        assert_eq!(citations[3].reasons, [Usage::Gamma.reason()]);
        assert_eq!(citations[4].id, "arXiv:hep-ph/9802376");
        assert_eq!(citations[4].reasons, [Usage::Color.reason()]);
        assert!(
            citations[2]
                .to_bibtex()
                .contains("10.48550/arXiv.2601.19982")
        );
        assert!(
            citations[3]
                .to_bibtex()
                .contains("10.48550/arXiv.1707.06453")
        );
        assert!(
            citations[4]
                .to_bibtex()
                .contains("10.1142/S0217751X99000038")
        );
    }

    #[test]
    fn concurrent_operations_preserve_all_reasons() {
        let usage = CitationUsage::default();
        std::thread::scope(|scope| {
            for operation in [
                Usage::Tensor,
                Usage::Network,
                Usage::Evaluation,
                Usage::Contraction,
                Usage::Canonicalization,
                Usage::DiracAdjoint,
                Usage::Notation,
            ] {
                let usage = &usage;
                scope.spawn(move || usage.record(operation));
            }
        });
        let citations = usage.get_citations();
        assert_eq!(citations[0].reasons.len(), 4);
        assert_eq!(citations[1].reasons.len(), 4);
    }

    #[test]
    fn idenso_alone_always_credits_spenso_structure_and_display() {
        let usage = CitationUsage::default();
        usage.record(Usage::DiracAdjoint);
        let citations = usage.get_citations();
        assert_eq!(citations.len(), 2);
        assert_eq!(citations[0].reference, "Lucien Huber. Spenso (2026).");
        assert!(citations[0].reasons[0].contains("Tensor-expression structure and display"));
        assert_eq!(
            citations[1].reference,
            "Lucien Huber, Ben Ruijl. Idenso (2026)."
        );
    }
}
