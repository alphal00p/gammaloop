use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

use eyre::{Result, bail};

/// A dictionary-valued Typst configuration assembled from closed native values.
///
/// This type deliberately has no raw-expression constructor. Executable Typst
/// can enter a configuration only through a statically imported module symbol.
#[derive(Clone, Debug, PartialEq)]
pub struct TypstConfig {
    fields: BTreeMap<String, TypstValue>,
}

impl TypstConfig {
    pub fn new(fields: BTreeMap<String, TypstValue>) -> Result<Self> {
        let config = Self { fields };
        config.validate()?;
        Ok(config)
    }

    pub(crate) fn source(&self) -> Result<String> {
        self.validate()?;
        dictionary_source(&self.fields)
    }

    pub(crate) fn module_aliases(&self) -> BTreeSet<&str> {
        let mut aliases = BTreeSet::new();
        for value in self.fields.values() {
            value.collect_module_aliases(&mut aliases);
        }
        aliases
    }

    fn validate(&self) -> Result<()> {
        for value in self.fields.values() {
            value.validate()?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TypstLengthUnit {
    Pt,
    Mm,
    Cm,
    In,
    Em,
}

impl TypstLengthUnit {
    fn name(self) -> &'static str {
        match self {
            Self::Pt => "pt",
            Self::Mm => "mm",
            Self::Cm => "cm",
            Self::In => "in",
            Self::Em => "em",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TypstAngleUnit {
    Deg,
    Rad,
}

impl TypstAngleUnit {
    fn name(self) -> &'static str {
        match self {
            Self::Deg => "deg",
            Self::Rad => "rad",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TypstMathScript {
    Symbol(String),
    Integer(i64),
}

impl TypstMathScript {
    fn validate(&self) -> Result<()> {
        if let Self::Symbol(value) = self {
            validate_math_token(value, "math script")?;
        }
        Ok(())
    }

    fn source(&self) -> String {
        match self {
            Self::Symbol(value) => value.clone(),
            Self::Integer(value) => value.to_string(),
        }
    }
}

/// Closed native Typst values accepted by the renderer.
#[derive(Clone, Debug, PartialEq)]
pub enum TypstValue {
    None,
    Auto,
    Bool(bool),
    Integer(i64),
    Float(f64),
    String(String),
    Array(Vec<Self>),
    Dictionary(BTreeMap<String, Self>),
    Length(f64, TypstLengthUnit),
    Ratio(f64),
    RelativeLength {
        ratio: Option<f64>,
        length: Option<(f64, TypstLengthUnit)>,
    },
    Angle(f64, TypstAngleUnit),
    Fraction(f64),
    NamedColor(String),
    HexColor(String),
    Rgba(u8, u8, u8, Option<u8>),
    Text {
        text: String,
        options: BTreeMap<String, Self>,
    },
    Math {
        symbol: String,
        subscript: Option<TypstMathScript>,
        superscript: Option<TypstMathScript>,
    },
    ModuleSymbol {
        alias: String,
        path: Vec<String>,
    },
    Call {
        callee: Box<Self>,
        positional: Vec<Self>,
        named: BTreeMap<String, Self>,
    },
    Bind {
        function: Box<Self>,
        positional: Vec<Self>,
        named: BTreeMap<String, Self>,
    },
}

impl TypstValue {
    fn validate(&self) -> Result<()> {
        match self {
            Self::Float(value)
            | Self::Length(value, _)
            | Self::Ratio(value)
            | Self::Angle(value, _)
            | Self::Fraction(value) => validate_finite(*value)?,
            Self::RelativeLength { ratio, length } => {
                if ratio.is_none() && length.is_none() {
                    bail!("a relative length needs a ratio or length component");
                }
                if let Some(value) = ratio {
                    validate_finite(*value)?;
                }
                if let Some((value, _)) = length {
                    validate_finite(*value)?;
                }
            }
            Self::Array(values) => {
                for value in values {
                    value.validate()?;
                }
            }
            Self::Dictionary(values) => {
                for value in values.values() {
                    value.validate()?;
                }
            }
            Self::NamedColor(value) => validate_named_color(value)?,
            Self::HexColor(value) => validate_hex_color(value)?,
            Self::Text { options, .. } => {
                for (name, value) in options {
                    validate_identifier(name, "text option")?;
                    value.validate()?;
                }
            }
            Self::Math {
                symbol,
                subscript,
                superscript,
            } => {
                validate_math_token(symbol, "math symbol")?;
                if let Some(script) = subscript {
                    script.validate()?;
                }
                if let Some(script) = superscript {
                    script.validate()?;
                }
            }
            Self::ModuleSymbol { alias, path } => {
                validate_module_alias(alias)?;
                if path.is_empty() {
                    bail!("a Typst module symbol needs an export path");
                }
                for part in path {
                    validate_identifier(part, "module export")?;
                }
            }
            Self::Call {
                callee,
                positional,
                named,
            } => {
                if !matches!(
                    callee.as_ref(),
                    Self::ModuleSymbol { .. } | Self::Bind { .. }
                ) {
                    bail!("a Typst call must target an imported module function or binding");
                }
                callee.validate()?;
                validate_arguments(positional, named)?;
            }
            Self::Bind {
                function,
                positional,
                named,
            } => {
                if !matches!(function.as_ref(), Self::ModuleSymbol { .. }) {
                    bail!("a Typst binding must target an imported module function");
                }
                function.validate()?;
                validate_arguments(positional, named)?;
            }
            Self::None
            | Self::Auto
            | Self::Bool(_)
            | Self::Integer(_)
            | Self::String(_)
            | Self::Rgba(..) => {}
        }
        Ok(())
    }

    fn source(&self) -> Result<String> {
        self.validate()?;
        Ok(match self {
            Self::None => "none".to_owned(),
            Self::Auto => "auto".to_owned(),
            Self::Bool(value) => value.to_string(),
            Self::Integer(value) => value.to_string(),
            Self::Float(value) => typst_float(*value),
            Self::String(value) => typst_string(value),
            Self::Array(values) => tuple_source(values)?,
            Self::Dictionary(values) => dictionary_source(values)?,
            Self::Length(value, unit) => {
                format!("{}{unit}", typst_float(*value), unit = unit.name())
            }
            Self::Ratio(value) => format!("{}%", typst_float(*value)),
            Self::RelativeLength { ratio, length } => {
                let mut terms = Vec::with_capacity(2);
                if let Some(value) = ratio {
                    terms.push(format!("{}%", typst_float(*value)));
                }
                if let Some((value, unit)) = length {
                    terms.push(format!("{}{unit}", typst_float(*value), unit = unit.name()));
                }
                format!("({})", terms.join(" + "))
            }
            Self::Angle(value, unit) => {
                format!("{}{unit}", typst_float(*value), unit = unit.name())
            }
            Self::Fraction(value) => format!("{}fr", typst_float(*value)),
            Self::NamedColor(value) => value.clone(),
            Self::HexColor(value) => format!("rgb({})", typst_string(value)),
            Self::Rgba(red, green, blue, alpha) => {
                let mut values = vec![red.to_string(), green.to_string(), blue.to_string()];
                if let Some(alpha) = alpha {
                    values.push(alpha.to_string());
                }
                format!("rgb({})", values.join(", "))
            }
            Self::Text { text, options } => {
                let mut arguments = vec![typst_string(text)];
                for (name, value) in options {
                    arguments.push(format!("{name}: {}", value.source()?));
                }
                format!("text({})", arguments.join(", "))
            }
            Self::Math {
                symbol,
                subscript,
                superscript,
            } => {
                let mut math = symbol.clone();
                if let Some(script) = subscript {
                    let _ = write!(math, "_({})", script.source());
                }
                if let Some(script) = superscript {
                    let _ = write!(math, "^({})", script.source());
                }
                format!("${math}$")
            }
            Self::ModuleSymbol { alias, path } => format!("{alias}.{}", path.join(".")),
            Self::Call {
                callee,
                positional,
                named,
            } => format!(
                "{}({})",
                callee.source()?,
                arguments_source(positional, named)?
            ),
            Self::Bind {
                function,
                positional,
                named,
            } => format!(
                "{}.with({})",
                function.source()?,
                arguments_source(positional, named)?
            ),
        })
    }

    fn collect_module_aliases<'a>(&'a self, aliases: &mut BTreeSet<&'a str>) {
        match self {
            Self::Array(values) => {
                for value in values {
                    value.collect_module_aliases(aliases);
                }
            }
            Self::Dictionary(values) => {
                for value in values.values() {
                    value.collect_module_aliases(aliases);
                }
            }
            Self::Text { options, .. } => {
                for value in options.values() {
                    value.collect_module_aliases(aliases);
                }
            }
            Self::ModuleSymbol { alias, .. } => {
                aliases.insert(alias);
            }
            Self::Call {
                callee,
                positional,
                named,
            }
            | Self::Bind {
                function: callee,
                positional,
                named,
            } => {
                callee.collect_module_aliases(aliases);
                for value in positional {
                    value.collect_module_aliases(aliases);
                }
                for value in named.values() {
                    value.collect_module_aliases(aliases);
                }
            }
            _ => {}
        }
    }
}

fn validate_arguments(
    positional: &[TypstValue],
    named: &BTreeMap<String, TypstValue>,
) -> Result<()> {
    for value in positional {
        value.validate()?;
    }
    for (name, value) in named {
        validate_identifier(name, "named argument")?;
        value.validate()?;
    }
    Ok(())
}

fn arguments_source(
    positional: &[TypstValue],
    named: &BTreeMap<String, TypstValue>,
) -> Result<String> {
    let mut values = positional
        .iter()
        .map(TypstValue::source)
        .collect::<Result<Vec<_>>>()?;
    for (name, value) in named {
        values.push(format!("{name}: {}", value.source()?));
    }
    Ok(values.join(", "))
}

fn tuple_source(values: &[TypstValue]) -> Result<String> {
    let values = values
        .iter()
        .map(TypstValue::source)
        .collect::<Result<Vec<_>>>()?;
    Ok(if values.len() == 1 {
        format!("({},)", values[0])
    } else {
        format!("({})", values.join(", "))
    })
}

fn dictionary_source(values: &BTreeMap<String, TypstValue>) -> Result<String> {
    let entries = values
        .iter()
        .map(|(name, value)| Ok(format!("({}): {}", typst_string(name), value.source()?)))
        .collect::<Result<Vec<_>>>()?;
    Ok(if entries.is_empty() {
        "(:)".to_owned()
    } else {
        format!("({})", entries.join(", "))
    })
}

fn typst_float(value: f64) -> String {
    let value = value.to_string();
    if value.contains(['.', 'e', 'E']) {
        value
    } else {
        format!("{value}.0")
    }
}

fn typst_string(value: &str) -> String {
    let mut output = String::with_capacity(value.len() + 2);
    output.push('"');
    for character in value.chars() {
        match character {
            '\\' => output.push_str("\\\\"),
            '"' => output.push_str("\\\""),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            character if character.is_control() => {
                let _ = write!(output, "\\u{{{:x}}}", character as u32);
            }
            character => output.push(character),
        }
    }
    output.push('"');
    output
}

fn validate_finite(value: f64) -> Result<()> {
    if !value.is_finite() {
        bail!("Typst numbers must be finite");
    }
    Ok(())
}

fn validate_identifier(value: &str, what: &str) -> Result<()> {
    let mut characters = value.chars();
    let valid = characters
        .next()
        .is_some_and(|first| first == '_' || first.is_alphabetic())
        && characters
            .all(|character| character == '_' || character == '-' || character.is_alphanumeric());
    if !valid {
        bail!("invalid Typst {what} {value:?}");
    }
    Ok(())
}

fn validate_module_alias(value: &str) -> Result<()> {
    let mut characters = value.chars();
    let valid = characters
        .next()
        .is_some_and(|first| first == '_' || first.is_ascii_alphabetic())
        && characters.all(|character| character == '_' || character.is_ascii_alphanumeric());
    if !valid {
        bail!("invalid Typst module alias {value:?}");
    }
    Ok(())
}

fn validate_math_token(value: &str, what: &str) -> Result<()> {
    if value.is_empty()
        || !value
            .chars()
            .all(|character| character == '_' || character == '-' || character.is_alphanumeric())
    {
        bail!("invalid Typst {what} {value:?}");
    }
    Ok(())
}

fn validate_named_color(value: &str) -> Result<()> {
    const COLORS: &[&str] = &[
        "black", "gray", "silver", "white", "navy", "blue", "aqua", "teal", "eastern", "purple",
        "fuchsia", "maroon", "red", "orange", "yellow", "olive", "green", "lime",
    ];
    if !COLORS.contains(&value) {
        bail!("unknown Typst named color {value:?}");
    }
    Ok(())
}

fn validate_hex_color(value: &str) -> Result<()> {
    let Some(digits) = value.strip_prefix('#') else {
        bail!("invalid Typst hexadecimal color {value:?}");
    };
    if !matches!(digits.len(), 3 | 4 | 6 | 8)
        || !digits
            .chars()
            .all(|character| character.is_ascii_hexdigit())
    {
        bail!("invalid Typst hexadecimal color {value:?}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serializes_closed_values_and_module_calls() {
        let function = TypstValue::ModuleSymbol {
            alias: "styles".to_owned(),
            path: vec!["edge_style".to_owned()],
        };
        let config = TypstConfig::new(BTreeMap::from([
            (
                "title".to_owned(),
                TypstValue::String("a \"title\"".to_owned()),
            ),
            (
                "style".to_owned(),
                TypstValue::Call {
                    callee: Box::new(function),
                    positional: vec![TypstValue::Length(1.0, TypstLengthUnit::Pt)],
                    named: BTreeMap::new(),
                },
            ),
        ]))
        .unwrap();
        assert_eq!(config.module_aliases(), BTreeSet::from(["styles"]));
        let source = config.source().unwrap();
        assert!(source.contains("styles.edge_style(1.0pt)"));
        assert!(source.contains("a \\\"title\\\""));
    }

    #[test]
    fn rejects_non_finite_values_and_injected_symbols() {
        assert!(
            TypstConfig::new(BTreeMap::from([(
                "value".to_owned(),
                TypstValue::Float(f64::NAN),
            )]))
            .is_err()
        );
        assert!(
            TypstConfig::new(BTreeMap::from([(
                "value".to_owned(),
                TypstValue::ModuleSymbol {
                    alias: "styles; panic()".to_owned(),
                    path: vec!["value".to_owned()],
                },
            )]))
            .is_err()
        );
    }

    #[test]
    fn escapes_strings_without_creating_typst_syntax() {
        let value = TypstValue::String("slash\\quote\"line\nreturn\rtab\tbell\u{7}".to_owned());
        assert_eq!(
            value.source().unwrap(),
            "\"slash\\\\quote\\\"line\\nreturn\\rtab\\tbell\\u{7}\""
        );
    }
}
