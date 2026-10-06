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
            Self::Tensor => {
                "Constructed or manipulated typed tensors, tensor expressions, or tensor structures."
            }
            Self::Network => "Executed a tensor network using component data.",
            Self::Evaluation => "Evaluated tensor components numerically through Symbolica.",
            Self::Contraction => {
                "Ran symbolic tensor contraction, including metric and chain/trace reduction."
            }
            Self::Algebra => "Ran the tensor-algebra simplification scheduler.",
            Self::Gamma => "Enabled Dirac gamma identities in tensor-algebra simplification.",
            Self::Color => "Enabled color identities in tensor-algebra simplification.",
            Self::Epsilon => "Enabled epsilon identities in tensor-algebra simplification.",
            Self::Canonicalization => {
                "Canonicalized symbolic tensor contractions and dummy indices."
            }
            Self::DiracAdjoint => "Constructed a Dirac adjoint of a tensor expression.",
            Self::Notation => "Rewrote symbolic tensor indices or dot, chain, and trace notation.",
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
            spenso_reasons.push(
                "Provided the Symbolica tensor-expression structure and display used by Idenso."
                    .into(),
            );
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
}"#.into(),
            reasons: spenso_reasons,
            description: "Spenso defines tensor-expression structures and their display in Symbolica, representation-aware tensors, dense and sparse component storage, and tensor-network execution. Reasons summarize successful operations across the current process.".into(),
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
                description: "Idenso implements symbolic tensor contraction, canonicalization, Dirac adjoints, and Dirac, color, and epsilon identities. Enabled algebra families describe the requested configuration; they do not certify that a particular identity changed the result.".into(),
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
        assert_eq!(citations.len(), 2);
        assert_eq!(
            citations[1].reasons,
            [Usage::Algebra.reason(), Usage::Color.reason()]
        );
        for citation in &citations {
            assert!(citation.bibtex.contains(&citation.id));
            assert!(!citation.description.is_empty());
            assert_eq!(citation.to_bibtex(), citation.bibtex);
            assert!(citation.to_markdown(false).contains(&citation.reasons[0]));
        }
        assert_eq!(usage.get_citations()[1].reasons, citations[1].reasons);
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
        assert!(citations[0].reasons[0].contains("tensor-expression structure and display"));
        assert_eq!(
            citations[1].reference,
            "Lucien Huber, Ben Ruijl. Idenso (2026)."
        );
    }
}
