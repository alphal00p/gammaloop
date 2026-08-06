use symbolica::prelude::*;

use color_eyre::Result;
use eyre::eyre;

fn escaped_dot_string(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for c in value.chars() {
        match c {
            '\\' => escaped.push_str("\\\\"),
            '"' => escaped.push_str("\\\""),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            _ => escaped.push(c),
        }
    }
    escaped
}

/// Escapes raw text as a complete quoted DOT attribute value, including outer quotes.
pub(crate) fn dot_attr_value(value: &str) -> String {
    format!("\"{}\"", escaped_dot_string(value))
}

/// Escapes raw text for a DOT statement value whose writer supplies the outer quotes.
pub(crate) fn dot_statement_value(value: &str) -> String {
    escaped_dot_string(value)
}

/// Decodes the quoted-string escapes retained by the DOT parser.
///
/// Unknown Graphviz escapes are preserved literally. This matters for embedded languages such as
/// TOML, where a backslash may itself be meaningful after the DOT layer has been removed.
pub(crate) fn decode_dot_string(value: &str) -> Result<String> {
    let value = value
        .strip_prefix('"')
        .and_then(|value| value.strip_suffix('"'))
        .unwrap_or(value);
    let mut decoded = String::with_capacity(value.len());
    let mut chars = value.chars();
    while let Some(character) = chars.next() {
        if character != '\\' {
            decoded.push(character);
            continue;
        }

        let escaped = chars
            .next()
            .ok_or_else(|| eyre!("DOT string ends with an incomplete escape"))?;
        match escaped {
            '\\' => decoded.push('\\'),
            '"' => decoded.push('"'),
            'n' => decoded.push('\n'),
            'r' => decoded.push('\r'),
            't' => decoded.push('\t'),
            other => {
                decoded.push('\\');
                decoded.push(other);
            }
        }
    }
    Ok(decoded)
}

pub trait ToQuoted {
    fn to_quoted(&self) -> String;
}

// impl ToQuoted for SmartString<LazyCompact> {
//     fn to_quoted(&self) -> String {
//         format!("\"{}\"", self)
//     }
// }

pub trait ToOrderedSimple {
    fn to_ordered_simple(&self) -> String;
}

impl<A> ToOrderedSimple for A
where
    A: AtomCore,
{
    fn to_ordered_simple(&self) -> String {
        self.to_canonically_ordered_string(
            CanonicalOrderingSettings::new()
                .include_namespace(false)
                .include_attributes(false)
                .hide_namespace(None),
        )
    }
}
impl<A> ToQuoted for A
where
    A: AtomCore,
{
    fn to_quoted(&self) -> String {
        // let mut opts = PrintOptions::file();
        // opts.hide_namespace = Some("gammalooprs");
        self.to_canonically_ordered_string(
            CanonicalOrderingSettings::new()
                .include_namespace(true)
                .include_attributes(false)
                .hide_namespace(Some("gammalooprs")),
        )
        .to_string()
        //.printer(opts))
    }
}

pub trait FromStripedStr: Sized {
    fn strip_from(string: &str) -> Result<Self>;
}

pub trait StripParse {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T>;
}

impl StripParse for String {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T> {
        T::strip_from(self.as_str())
    }
}
impl StripParse for &String {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T> {
        T::strip_from(self.as_str())
    }
}
impl StripParse for &str {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T> {
        T::strip_from(self)
    }
}
impl FromStripedStr for String {
    fn strip_from(string: &str) -> Result<Self> {
        Ok(string
            .strip_prefix('"')
            .unwrap_or(string)
            .strip_suffix('"')
            .unwrap_or(string)
            .into())
    }
}

impl FromStripedStr for Atom {
    fn strip_from(string: &str) -> Result<Self> {
        let a = string
            .strip_prefix('"')
            .unwrap_or(string)
            .strip_suffix('"')
            .unwrap_or(string);
        try_parse!(a).map_err(|e| eyre!("Symbolica parsing error: {e}"))
    }
}

impl FromStripedStr for bool {
    fn strip_from(s: &str) -> Result<Self> {
        Ok(s.strip_prefix('"')
            .unwrap_or(s)
            .strip_suffix('"')
            .unwrap_or(s)
            .parse()?)
    }
}
impl FromStripedStr for i32 {
    fn strip_from(s: &str) -> Result<Self> {
        Ok(s.strip_prefix('"')
            .unwrap_or(s)
            .strip_suffix('"')
            .unwrap_or(s)
            .parse()?)
    }
}
impl FromStripedStr for usize {
    fn strip_from(s: &str) -> Result<Self> {
        Ok(s.strip_prefix('"')
            .unwrap_or(s)
            .strip_suffix('"')
            .unwrap_or(s)
            .parse()?)
    }
}

#[cfg(test)]
mod tests {
    use linnet::parser::DotGraph;

    use super::{decode_dot_string, dot_attr_value, dot_statement_value};

    const RAW_DOT_VALUE: &str = "quote:\" slash:\\ lf:\n cr:\r tab:\t";
    const ESCAPED_DOT_VALUE: &str = r#"quote:\" slash:\\ lf:\n cr:\r tab:\t"#;

    #[test]
    fn dot_values_escape_quoted_string_content() {
        assert_eq!(dot_statement_value(RAW_DOT_VALUE), ESCAPED_DOT_VALUE);
        assert_eq!(
            dot_attr_value(RAW_DOT_VALUE),
            format!(r#""{ESCAPED_DOT_VALUE}""#)
        );
    }

    #[test]
    fn dot_statement_value_survives_dot_round_trip() {
        let mut graph: DotGraph = DotGraph::from_string("digraph G {}").unwrap();
        graph
            .global_data
            .statements
            .insert("escaped".to_string(), dot_statement_value(RAW_DOT_VALUE));

        let reparsed: DotGraph = DotGraph::from_string(graph.debug_dot()).unwrap();

        assert_eq!(
            reparsed.global_data.statements["escaped"],
            ESCAPED_DOT_VALUE
        );
        assert_eq!(
            decode_dot_string(&reparsed.global_data.statements["escaped"]).unwrap(),
            RAW_DOT_VALUE
        );
    }

    #[test]
    fn dot_codec_round_trips_unicode_and_unknown_escapes() {
        let raw = "eta:η theta:θ regex:\\d quote:\" newline:\n";
        let encoded = dot_statement_value(raw);

        assert_eq!(decode_dot_string(&encoded).unwrap(), raw);
    }
}
