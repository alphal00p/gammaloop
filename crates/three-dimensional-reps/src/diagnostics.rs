use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProfileWarningRow {
    pub label: String,
    pub family: String,
    pub has_repeated_propagators: bool,
}

pub fn profile_warnings(rows: &[ProfileWarningRow]) -> Vec<String> {
    let mut labels = rows
        .iter()
        .filter(|row| row.family == "pure_ltd" && row.has_repeated_propagators)
        .map(|row| row.label.clone())
        .collect::<Vec<_>>();
    labels.sort();
    labels.dedup();

    if labels.is_empty() {
        return Vec::new();
    }

    vec![format!(
        "Reminder: pure_ltd profiling rows for repeated-propagator graphs ({}) do not include numerator derivatives and are not the valid repeated-propagator LTD formula. Use the public ltd representation for derivative-free repeated propagators; pure_ltd rows are only mass-shift diagnostics.",
        labels.join(", ")
    )]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_warning_mentions_repeated_pure_ltd_rows_only() {
        let warnings = profile_warnings(&[
            ProfileWarningRow {
                label: "pure-ltd-ref".to_string(),
                family: "pure_ltd".to_string(),
                has_repeated_propagators: true,
            },
            ProfileWarningRow {
                label: "ltd".to_string(),
                family: "ltd".to_string(),
                has_repeated_propagators: true,
            },
            ProfileWarningRow {
                label: "pure-ltd-plain".to_string(),
                family: "pure_ltd".to_string(),
                has_repeated_propagators: false,
            },
        ]);

        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("pure-ltd-ref"));
        assert!(warnings[0].contains("public ltd representation"));
        assert!(warnings[0].contains("numerator derivatives"));
        assert!(!warnings[0].contains("hybrid"));
    }
}
