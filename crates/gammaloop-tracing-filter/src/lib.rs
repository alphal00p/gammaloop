use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use tracing::{Event, Metadata, Subscriber, field::Visit, level_filters::LevelFilter};
use tracing_subscriber::layer::{Context, Filter};
use tracing_subscriber::registry::LookupSpan;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GammaLogFilter {
    directives: Vec<Directive>,
    max_level_hint: LevelFilter,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Directive {
    target: Option<String>,
    tags: Vec<TagRequirement>,
    level: LevelFilter,
    order: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum TagRequirement {
    Present(String),
    Absent(String),
    Value { name: String, expected: String },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Specificity {
    target_len: usize,
    tag_count: usize,
    order: usize,
}

#[derive(Debug, Default)]
struct BoolFieldVisitor {
    present_fields: BTreeSet<String>,
    field_values: BTreeMap<String, String>,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ParseError {
    #[error("invalid filter directive")]
    InvalidDirective,
    #[error("invalid filter directive: {0}")]
    InvalidDirectiveMessage(String),
    #[error("unknown log level '{0}'")]
    UnknownLevel(String),
}

impl GammaLogFilter {
    pub fn parse(user_spec: &str) -> Result<Self, ParseError> {
        Self::parse_with_default(user_spec, None)
    }

    pub fn parse_with_default(
        user_spec: &str,
        default_level: Option<LevelFilter>,
    ) -> Result<Self, ParseError> {
        let mut directives = Vec::new();
        for (order, raw_directive) in split_directives(user_spec)?.into_iter().enumerate() {
            if raw_directive.trim().is_empty() {
                continue;
            }
            directives.push(parse_directive(raw_directive, order)?);
        }

        if directives.is_empty()
            && let Some(default_level) = default_level
        {
            directives.push(Directive {
                target: None,
                tags: Vec::new(),
                level: default_level,
                order: 0,
            });
        }

        let max_level_hint = directives
            .iter()
            .map(|directive| directive.level)
            .max_by_key(|level| level_rank(*level))
            .unwrap_or(LevelFilter::OFF);

        Ok(Self {
            directives,
            max_level_hint,
        })
    }

    pub fn max_level_hint(&self) -> Option<LevelFilter> {
        Some(self.max_level_hint)
    }

    pub fn is_effectively_off(user_spec: &str) -> bool {
        let Ok(parts) = split_directives(user_spec) else {
            return false;
        };
        let mut saw_any = false;
        for part in parts {
            let trimmed = part.trim();
            if trimmed.is_empty() {
                continue;
            }
            saw_any = true;
            let Ok(directive) = parse_directive(trimmed, 0) else {
                return false;
            };
            if directive.level != LevelFilter::OFF {
                return false;
            }
        }
        saw_any
    }

    #[cfg(test)]
    fn enabled_for_test(&self, target: &str, level: LevelFilter, tags: &[(&str, bool)]) -> bool {
        let present_fields = tags
            .iter()
            .map(|(name, _)| (*name).to_string())
            .collect::<BTreeSet<_>>();
        let tag_map = tags
            .iter()
            .map(|(name, value)| ((*name).to_string(), value.to_string()))
            .collect::<BTreeMap<_, _>>();
        self.selected_level(target, Some(&present_fields), Some(&tag_map))
            .is_some_and(|selected| level_rank(selected) >= level_rank(level))
    }

    fn selected_level(
        &self,
        target: &str,
        present_fields: Option<&BTreeSet<String>>,
        tags: Option<&BTreeMap<String, String>>,
    ) -> Option<LevelFilter> {
        self.directives
            .iter()
            .filter(|directive| directive.matches_target(target))
            .filter(|directive| match (present_fields, tags) {
                (Some(present_fields), Some(tags)) => directive.matches_tags(present_fields, tags),
                _ => directive.tags.is_empty(),
            })
            .max_by_key(|directive| directive.specificity())
            .map(|directive| directive.level)
    }

    fn candidate_level(&self, target: &str) -> Option<LevelFilter> {
        self.directives
            .iter()
            .filter(|directive| directive.matches_target(target))
            .map(|directive| directive.level)
            .max_by_key(|level| level_rank(*level))
    }

    fn callsite_interest(
        &self,
        metadata: &'static Metadata<'static>,
    ) -> tracing::subscriber::Interest {
        let has_dynamic = self
            .directives
            .iter()
            .filter(|directive| directive.matches_target(metadata.target()))
            .any(|directive| directive.could_match_dynamically(metadata));

        if has_dynamic {
            return tracing::subscriber::Interest::sometimes();
        }

        let selected = self
            .directives
            .iter()
            .filter(|directive| directive.matches_target(metadata.target()))
            .filter(|directive| directive.matches_metadata_fields(metadata))
            .max_by_key(|directive| directive.specificity())
            .map(|directive| directive.level);

        match selected {
            Some(level) if event_level_enabled(metadata, level) => {
                tracing::subscriber::Interest::always()
            }
            _ => tracing::subscriber::Interest::never(),
        }
    }
}

impl Directive {
    fn matches_target(&self, event_target: &str) -> bool {
        let Some(target) = self.target.as_deref() else {
            return true;
        };
        event_target == target
            || event_target
                .strip_prefix(target)
                .is_some_and(|rest| rest.starts_with("::"))
    }

    fn matches_tags(
        &self,
        present_fields: &BTreeSet<String>,
        tags: &BTreeMap<String, String>,
    ) -> bool {
        self.tags.iter().all(|requirement| match requirement {
            TagRequirement::Present(name) => present_fields.contains(name),
            TagRequirement::Absent(name) => !present_fields.contains(name),
            TagRequirement::Value { name, expected } => {
                tags.get(name).is_some_and(|actual| actual == expected)
            }
        })
    }

    fn matches_metadata_fields(&self, metadata: &'static Metadata<'static>) -> bool {
        self.tags.iter().all(|requirement| match requirement {
            TagRequirement::Present(name) => metadata.fields().field(name).is_some(),
            TagRequirement::Absent(name) => metadata.fields().field(name).is_none(),
            TagRequirement::Value { .. } => false,
        })
    }

    fn could_match_dynamically(&self, metadata: &'static Metadata<'static>) -> bool {
        self.tags.iter().all(|requirement| match requirement {
            TagRequirement::Present(name) => metadata.fields().field(name).is_some(),
            TagRequirement::Absent(name) => metadata.fields().field(name).is_none(),
            TagRequirement::Value { name, .. } => metadata.fields().field(name).is_some(),
        }) && self.tags.iter().any(TagRequirement::is_value_match)
    }

    fn specificity(&self) -> Specificity {
        Specificity {
            target_len: self.target.as_ref().map_or(0, String::len),
            tag_count: self.tags.len(),
            order: self.order,
        }
    }
}

impl Visit for BoolFieldVisitor {
    fn record_bool(&mut self, field: &tracing::field::Field, value: bool) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), value.to_string());
    }

    fn record_f64(&mut self, field: &tracing::field::Field, value: f64) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), value.to_string());
    }

    fn record_i64(&mut self, field: &tracing::field::Field, value: i64) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), value.to_string());
    }

    fn record_u64(&mut self, field: &tracing::field::Field, value: u64) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), value.to_string());
    }

    fn record_i128(&mut self, field: &tracing::field::Field, value: i128) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), value.to_string());
    }

    fn record_u128(&mut self, field: &tracing::field::Field, value: u128) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), value.to_string());
    }

    fn record_str(&mut self, field: &tracing::field::Field, value: &str) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), value.to_string());
    }

    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        self.present_fields.insert(field.name().to_string());
        self.field_values
            .insert(field.name().to_string(), format!("{value:?}"));
    }
}

impl<S> Filter<S> for GammaLogFilter
where
    S: Subscriber + for<'lookup> LookupSpan<'lookup>,
{
    fn enabled(&self, metadata: &Metadata<'_>, _cx: &Context<'_, S>) -> bool {
        self.candidate_level(metadata.target())
            .is_some_and(|selected| event_level_enabled(metadata, selected))
    }

    fn event_enabled(&self, event: &Event<'_>, _cx: &Context<'_, S>) -> bool {
        let mut visitor = BoolFieldVisitor::default();
        event.record(&mut visitor);
        self.selected_level(
            event.metadata().target(),
            Some(&visitor.present_fields),
            Some(&visitor.field_values),
        )
        .is_some_and(|selected| event_level_enabled(event.metadata(), selected))
    }

    fn callsite_enabled(
        &self,
        metadata: &'static Metadata<'static>,
    ) -> tracing::subscriber::Interest {
        self.callsite_interest(metadata)
    }

    fn max_level_hint(&self) -> Option<LevelFilter> {
        self.max_level_hint()
    }
}

fn event_level_enabled(metadata: &Metadata<'_>, selected: LevelFilter) -> bool {
    selected != LevelFilter::OFF
        && level_rank(selected) >= level_rank(LevelFilter::from_level(*metadata.level()))
}

fn level_rank(level: LevelFilter) -> u8 {
    match level {
        LevelFilter::OFF => 0,
        LevelFilter::ERROR => 1,
        LevelFilter::WARN => 2,
        LevelFilter::INFO => 3,
        LevelFilter::DEBUG => 4,
        LevelFilter::TRACE => 5,
    }
}

fn parse_directive(raw_directive: &str, order: usize) -> Result<Directive, ParseError> {
    let trimmed = raw_directive.trim();
    if trimmed.is_empty() {
        return Err(ParseError::InvalidDirective);
    }

    let (head, level) = if let Some(index) = find_top_level_equals(trimmed)? {
        let (head, tail) = trimmed.split_at(index);
        let level = parse_level(tail[1..].trim())?;
        (head.trim(), level)
    } else if !trimmed.contains('[') && parse_level(trimmed).is_ok() {
        ("", parse_level(trimmed)?)
    } else {
        (trimmed, LevelFilter::TRACE)
    };

    let (target, tags) = parse_head(head)?;
    Ok(Directive {
        target,
        tags,
        level,
        order,
    })
}

fn parse_head(head: &str) -> Result<(Option<String>, Vec<TagRequirement>), ParseError> {
    let trimmed = head.trim();
    if trimmed.is_empty() {
        return Ok((None, Vec::new()));
    }

    let Some(open_bracket) = trimmed.find('[') else {
        return Ok((Some(trimmed.to_string()), Vec::new()));
    };

    let close_bracket = trimmed
        .rfind(']')
        .ok_or_else(|| ParseError::InvalidDirectiveMessage(trimmed.to_string()))?;
    if close_bracket != trimmed.len() - 1 {
        return Err(ParseError::InvalidDirectiveMessage(trimmed.to_string()));
    }

    let target = trimmed[..open_bracket].trim();
    let bracket_body = trimmed[open_bracket + 1..close_bracket].trim();
    let tags = parse_tag_group(bracket_body)?;
    let target = (!target.is_empty()).then(|| target.to_string());
    Ok((target, tags))
}

fn parse_tag_group(raw: &str) -> Result<Vec<TagRequirement>, ParseError> {
    let inner = raw
        .strip_prefix('{')
        .and_then(|rest| rest.strip_suffix('}'))
        .ok_or_else(|| ParseError::InvalidDirectiveMessage(raw.to_string()))?
        .trim();
    if inner.is_empty() {
        return Err(ParseError::InvalidDirectiveMessage(raw.to_string()));
    }

    inner
        .split(',')
        .map(str::trim)
        .map(parse_tag_requirement)
        .collect()
}

fn parse_tag_requirement(raw: &str) -> Result<TagRequirement, ParseError> {
    if let Some(name) = raw.strip_prefix('!') {
        let name = name.trim();
        if !is_valid_ident(name) {
            return Err(ParseError::InvalidDirectiveMessage(raw.to_string()));
        }
        return Ok(TagRequirement::Absent(name.to_string()));
    }

    if let Some((name, expected)) = raw.split_once('=') {
        let name = name.trim();
        let expected = expected.trim();
        if !is_valid_ident(name) {
            return Err(ParseError::InvalidDirectiveMessage(raw.to_string()));
        }
        return Ok(TagRequirement::Value {
            name: name.to_string(),
            expected: parse_expected_value(expected),
        });
    }

    if !is_valid_ident(raw) {
        return Err(ParseError::InvalidDirectiveMessage(raw.to_string()));
    }

    Ok(TagRequirement::Present(raw.to_string()))
}

impl TagRequirement {
    fn is_value_match(&self) -> bool {
        matches!(self, TagRequirement::Value { .. })
    }
}

fn parse_expected_value(raw: &str) -> String {
    if raw.len() >= 2
        && ((raw.starts_with('"') && raw.ends_with('"'))
            || (raw.starts_with('\'') && raw.ends_with('\'')))
    {
        return raw[1..raw.len() - 1].to_string();
    }

    raw.to_string()
}

fn parse_level(raw: &str) -> Result<LevelFilter, ParseError> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "off" => Ok(LevelFilter::OFF),
        "error" => Ok(LevelFilter::ERROR),
        "warn" | "warning" => Ok(LevelFilter::WARN),
        "info" => Ok(LevelFilter::INFO),
        "debug" => Ok(LevelFilter::DEBUG),
        "trace" => Ok(LevelFilter::TRACE),
        _ => Err(ParseError::UnknownLevel(raw.trim().to_string())),
    }
}

fn split_directives(spec: &str) -> Result<Vec<&str>, ParseError> {
    let mut directives = Vec::new();
    let mut start = 0usize;
    let mut bracket_depth = 0usize;
    let mut brace_depth = 0usize;

    for (index, ch) in spec.char_indices() {
        match ch {
            '[' => bracket_depth += 1,
            ']' => {
                bracket_depth = bracket_depth
                    .checked_sub(1)
                    .ok_or(ParseError::InvalidDirective)?;
            }
            '{' => brace_depth += 1,
            '}' => {
                brace_depth = brace_depth
                    .checked_sub(1)
                    .ok_or(ParseError::InvalidDirective)?;
            }
            ',' if bracket_depth == 0 && brace_depth == 0 => {
                directives.push(&spec[start..index]);
                start = index + 1;
            }
            _ => {}
        }
    }

    if bracket_depth != 0 || brace_depth != 0 {
        return Err(ParseError::InvalidDirective);
    }

    directives.push(&spec[start..]);
    Ok(directives)
}

fn find_top_level_equals(spec: &str) -> Result<Option<usize>, ParseError> {
    let mut bracket_depth = 0usize;
    let mut brace_depth = 0usize;
    let mut found = None;

    for (index, ch) in spec.char_indices() {
        match ch {
            '[' => bracket_depth += 1,
            ']' => {
                bracket_depth = bracket_depth
                    .checked_sub(1)
                    .ok_or(ParseError::InvalidDirective)?;
            }
            '{' => brace_depth += 1,
            '}' => {
                brace_depth = brace_depth
                    .checked_sub(1)
                    .ok_or(ParseError::InvalidDirective)?;
            }
            '=' if bracket_depth == 0 && brace_depth == 0 => found = Some(index),
            _ => {}
        }
    }

    if bracket_depth != 0 || brace_depth != 0 {
        return Err(ParseError::InvalidDirective);
    }

    Ok(found)
}

fn is_valid_ident(name: &str) -> bool {
    let mut chars = name.chars();
    let Some(first) = chars.next() else {
        return false;
    };
    if !(first == '_' || first.is_ascii_alphabetic()) {
        return false;
    }
    chars.all(|ch| ch == '_' || ch.is_ascii_alphanumeric())
}

#[cfg(test)]
mod tests {
    use std::{
        collections::{BTreeMap, BTreeSet},
        fs,
        path::{Path, PathBuf},
    };

    use super::GammaLogFilter;
    use tracing::level_filters::LevelFilter;

    #[test]
    fn parses_documented_examples() {
        for spec in [
            "gammalooprs::uv::forest[{generation,uv}]=debug",
            "gammalooprs::uv::forest[{generation,uv,orientation,dump}]=debug",
            "gammalooprs::uv::forest[{generation,uv,term}]=debug",
            "gammalooprs::integrands::process::amplitude[{integration,subtraction}]=debug",
            "gammalooprs::integrands::process::cross_section[{integration,sample,inspect}]=debug",
            "gammalooprs::integrands::process::cross_section[{integration,cut,solver}]=debug",
            "gammalooprs::integrands::process::cross_section[{integration,!inspect}]=debug",
            "gammalooprs::integrands::process::cross_section[{integration,inspect=true}]=debug",
        ] {
            assert!(
                GammaLogFilter::parse(spec).is_ok(),
                "filter example should parse: {spec}"
            );
        }
    }

    #[test]
    fn all_run_card_log_directives_parse() {
        let run_cards_root =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/resources/run_cards");
        let mut files = Vec::new();
        collect_toml_files(&run_cards_root, &mut files);

        for path in files {
            let content = fs::read_to_string(&path).unwrap();
            for (key, spec) in extract_log_specs(&content) {
                assert!(
                    GammaLogFilter::parse(&spec).is_ok(),
                    "{key} in {} should parse: {spec}",
                    path.display()
                );
            }
        }
    }

    #[test]
    fn supports_full_negation_for_tags() {
        let filter = GammaLogFilter::parse(
            "gammalooprs::integrands::process::cross_section=info,gammalooprs::integrands::process::cross_section[{integration,!inspect}]=debug",
        )
        .unwrap();

        assert!(filter.enabled_for_test(
            "gammalooprs::integrands::process::cross_section",
            LevelFilter::DEBUG,
            &[("integration", true)],
        ));
        assert!(!filter.enabled_for_test(
            "gammalooprs::integrands::process::cross_section",
            LevelFilter::DEBUG,
            &[("integration", true), ("inspect", true)],
        ));
        assert!(filter.enabled_for_test(
            "gammalooprs::integrands::process::cross_section",
            LevelFilter::INFO,
            &[("integration", true), ("inspect", false)],
        ));
    }

    #[test]
    fn supports_explicit_boolean_value_matching() {
        let filter = GammaLogFilter::parse(
            "gammalooprs=info,gammalooprs::integrands::process::cross_section[{integration,inspect=true}]=debug",
        )
        .unwrap();

        assert!(filter.enabled_for_test(
            "gammalooprs::integrands::process::cross_section",
            LevelFilter::DEBUG,
            &[("integration", true), ("inspect", true)],
        ));
        assert!(!filter.enabled_for_test(
            "gammalooprs::integrands::process::cross_section",
            LevelFilter::DEBUG,
            &[("integration", true), ("inspect", false)],
        ));
        assert!(filter.enabled_for_test(
            "gammalooprs::integrands::process::cross_section",
            LevelFilter::INFO,
            &[("integration", true), ("inspect", false)],
        ));
    }

    #[test]
    fn supports_arbitrary_value_matching() {
        let filter = GammaLogFilter::parse(
            "gammalooprs=info,gammalooprs::integrands::process::cross_section[{integration,mode=summary}]=debug",
        )
        .unwrap();

        let present_fields = ["integration", "mode"]
            .into_iter()
            .map(str::to_string)
            .collect::<BTreeSet<_>>();
        let tags = [
            ("integration".to_string(), "true".to_string()),
            ("mode".to_string(), "summary".to_string()),
        ]
        .into_iter()
        .collect::<BTreeMap<_, _>>();
        assert!(
            filter
                .selected_level(
                    "gammalooprs::integrands::process::cross_section",
                    Some(&present_fields),
                    Some(&tags),
                )
                .is_some_and(|level| level == LevelFilter::DEBUG)
        );
    }

    #[test]
    fn base_target_filter_applies_when_specific_tags_do_not_match() {
        let filter = GammaLogFilter::parse(
            "gammalooprs=info,gammalooprs::uv::forest[{generation,uv,orientation,dump}]=debug",
        )
        .unwrap();

        assert!(filter.enabled_for_test(
            "gammalooprs::uv::forest",
            LevelFilter::DEBUG,
            &[
                ("generation", true),
                ("uv", true),
                ("orientation", true),
                ("dump", true),
            ],
        ));
        assert!(!filter.enabled_for_test(
            "gammalooprs::uv::forest",
            LevelFilter::DEBUG,
            &[("generation", true), ("uv", true)],
        ));
        assert!(filter.enabled_for_test(
            "gammalooprs::uv::forest",
            LevelFilter::INFO,
            &[("generation", true), ("uv", true)],
        ));
    }

    #[test]
    fn empty_specs_use_supplied_default() {
        let filter = GammaLogFilter::parse_with_default("", Some(LevelFilter::WARN)).unwrap();
        assert_eq!(filter.max_level_hint(), Some(LevelFilter::WARN));
    }

    #[test]
    fn effectively_off_accounts_for_tag_group_commas() {
        assert!(GammaLogFilter::is_effectively_off("off"));
        assert!(GammaLogFilter::is_effectively_off(
            "gammalooprs=off,symbolica=off"
        ));
        assert!(!GammaLogFilter::is_effectively_off(
            "gammalooprs=info,gammalooprs::uv::forest[{generation,uv}]=debug"
        ));
    }

    fn collect_toml_files(dir: &Path, out: &mut Vec<PathBuf>) {
        for entry in fs::read_dir(dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_dir() {
                collect_toml_files(&path, out);
            } else if path.extension().is_some_and(|ext| ext == "toml") {
                out.push(path);
            }
        }
    }

    fn extract_log_specs(content: &str) -> Vec<(String, String)> {
        content
            .lines()
            .filter_map(|line| {
                let trimmed = line.trim();
                ["display_directive", "logfile_directive"]
                    .into_iter()
                    .find_map(|key| {
                        let prefix = format!("{key} = ");
                        let rest = trimmed.strip_prefix(&prefix)?;
                        let quoted = rest.strip_prefix('"')?;
                        let value = quoted.split('"').next()?;
                        Some((key.to_string(), value.to_string()))
                    })
            })
            .collect()
    }
}
