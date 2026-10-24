use crate::{DisplayText, MathMl};
use std::borrow::Cow;
use symbolica::atom::{Symbol, UserData};
use thiserror::Error;

/// Invalid structured display input. See [`IndexDisplay::symbol`] for an example.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("invalid index display: {0}")]
pub struct IndexDisplayError(pub String);

/// A safe, printer-neutral description of one mathematical index label.
///
/// The tree is deliberately small: plugins may deserialize it from Symbolica
/// symbol metadata and turn it into Typst source without evaluating arbitrary
/// source supplied by an Atom payload.
/// See [`Self::to_mathml`] for an example rendering a structured label.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum IndexDisplay {
    Symbol(String),
    /// Upright text embedded in a mathematical label.
    Text(String),
    Number(i64),
    Sequence(Vec<IndexDisplay>),
    /// A nested argument/row list. Unlike `Sequence`, entries retain their
    /// comma/row boundary when used by a safe mathematical display node.
    List(Vec<IndexDisplay>),
    /// A validated mathematical display node.
    ///
    /// `head` is never emitted directly. Rendering dispatches through a fixed
    /// allowlist, so imported symbol metadata cannot inject Typst source.
    Math {
        head: String,
        arguments: Vec<IndexDisplay>,
    },
    Attach {
        base: Box<IndexDisplay>,
        top: Option<Box<IndexDisplay>>,
        bottom: Option<Box<IndexDisplay>>,
    },
}

impl IndexDisplay {
    /// Construct a nonempty label of at most 128 non-control characters.
    ///
    /// # Errors
    /// Returns an error for empty, oversized or control-containing names.
    /// See [`Self::to_mathml`] for a runnable example.
    pub fn symbol(name: impl Into<String>) -> Result<Self, IndexDisplayError> {
        let name = name.into();
        if name.is_empty() {
            return Err(IndexDisplayError(
                "an index label cannot be empty".to_owned(),
            ));
        }
        if name.chars().count() > 128 || name.chars().any(char::is_control) {
            return Err(IndexDisplayError(
                "an index label must contain at most 128 non-control characters".to_owned(),
            ));
        }
        Ok(Self::Symbol(name))
    }

    /// Attach a lower script, retaining an existing upper script. See [`Self::to_mathml`].
    #[must_use]
    pub fn with_bottom(self, bottom: IndexDisplay) -> Self {
        match self {
            Self::Attach {
                base,
                top,
                bottom: None,
            } => Self::Attach {
                base,
                top,
                bottom: Some(Box::new(bottom)),
            },
            display => Self::Attach {
                base: Box::new(display),
                top: None,
                bottom: Some(Box::new(bottom)),
            },
        }
    }

    /// Attach an upper script, retaining an existing lower script.
    ///
    /// ```
    /// use symbolica_notation::IndexDisplay;
    /// assert_eq!(IndexDisplay::Number(2).with_top(IndexDisplay::Number(3)).to_typst_source(), "attach(2,t:3)");
    /// ```
    #[must_use]
    pub fn with_top(self, top: IndexDisplay) -> Self {
        match self {
            Self::Attach {
                base,
                top: None,
                bottom,
            } => Self::Attach {
                base,
                top: Some(Box::new(top)),
                bottom,
            },
            display => Self::Attach {
                base: Box::new(display),
                top: Some(Box::new(top)),
                bottom: None,
            },
        }
    }

    /// Construct upright text of at most 1024 non-control characters.
    ///
    /// # Errors
    /// Returns an error for oversized or control-containing text.
    ///
    /// ```
    /// use symbolica_notation::{IndexDisplay, IndexDisplayError};
    /// assert_eq!(IndexDisplay::text("rate")?.to_native_string(), "rate");
    /// # Ok::<(), IndexDisplayError>(())
    /// ```
    pub fn text(value: impl Into<String>) -> Result<Self, IndexDisplayError> {
        let value = value.into();
        if value.chars().count() > 1024 || value.chars().any(char::is_control) {
            return Err(IndexDisplayError(
                "display text must contain at most 1024 non-control characters".to_owned(),
            ));
        }
        Ok(Self::Text(value))
    }

    /// Construct an allowlisted mathematical operation with at most 64 arguments.
    ///
    /// # Errors
    /// Returns an error for unsupported heads or too many arguments. Arity is
    /// checked by the target renderer.
    ///
    /// ```
    /// use symbolica_notation::{IndexDisplay, IndexDisplayError};
    /// let fraction = IndexDisplay::math("frac", vec![IndexDisplay::Number(1), IndexDisplay::Number(2)])?;
    /// assert_eq!(fraction.to_typst_source(), "frac(1,2)");
    /// # Ok::<(), IndexDisplayError>(())
    /// ```
    pub fn math(
        head: impl Into<String>,
        arguments: Vec<IndexDisplay>,
    ) -> Result<Self, IndexDisplayError> {
        let head = head.into();
        if !matches!(
            head.as_str(),
            "arg"
                | "add"
                | "sub"
                | "plus"
                | "neg"
                | "times"
                | "dot"
                | "factorial"
                | "mul"
                | "()"
                | "group"
                | "lr"
                | "pow"
                | "frac"
                | "root"
                | "call"
                | "op-call"
                | "op"
                | "accent"
                | "underline"
                | "overline"
                | "underbrace"
                | "overbrace"
                | "underbracket"
                | "overbracket"
                | "underparen"
                | "overparen"
                | "undershell"
                | "overshell"
                | "cancel"
                | "vec"
                | "mat"
                | "cases"
                | "class"
                | "mid"
                | "scripts"
                | "limits"
                | "stretch"
        ) {
            return Err(IndexDisplayError(format!(
                "unsupported mathematical display node {head:?}"
            )));
        }
        if arguments.len() > 64 {
            return Err(IndexDisplayError(
                "a mathematical display node cannot contain more than 64 arguments".to_owned(),
            ));
        }
        Ok(Self::Math { head, arguments })
    }

    fn escaped_typst_string(value: &str) -> String {
        value
            .replace('\\', "\\\\")
            .replace('"', "\\\"")
            .replace('\n', "\\n")
            .replace('\r', "\\r")
    }

    fn literal_text(&self) -> Option<&str> {
        match self {
            Self::Symbol(value) | Self::Text(value) => Some(value),
            _ => None,
        }
    }

    fn typst_math_source(head: &str, arguments: &[Self]) -> String {
        let sources = arguments
            .iter()
            .map(Self::to_typst_source)
            .collect::<Vec<_>>();
        let joined = |separator: &str| sources.join(separator);
        let unary = |prefix: &str, suffix: &str| {
            sources
                .first()
                .map(|body| format!("{prefix}{body}{suffix}"))
                .unwrap_or_else(|| r#"upright("?")"#.to_owned())
        };
        let call = |name: &str| format!("{name}({})", joined(","));

        match head {
            "arg" => joined(","),
            "add" => joined(" + "),
            "sub" => joined(" - "),
            "plus" => unary("+", ""),
            "neg" => unary("-", ""),
            "times" => joined(" times "),
            "dot" => joined(" dot "),
            "factorial" => unary("", "!"),
            "mul" => joined(" "),
            "()" | "group" | "lr" => unary("lr((", "))"),
            "pow" if sources.len() == 2 => {
                format!("attach({},t:{})", sources[0], sources[1])
            }
            "frac" if sources.len() == 2 => format!("frac({},{})", sources[0], sources[1]),
            "root" if sources.len() == 1 => format!("root({})", sources[0]),
            "root" if sources.len() == 2 => {
                format!("root({},{})", sources[0], sources[1])
            }
            "call" | "op-call" if !sources.is_empty() => {
                let body = sources[1..].join(",");
                format!("{}({body})", sources[0])
            }
            "op" if arguments.len() == 1 => {
                let text = arguments[0].literal_text().unwrap_or("?");
                format!(r#"op("{}")"#, Self::escaped_typst_string(text))
            }
            "accent" if arguments.len() == 2 => {
                let accent = arguments[1].literal_text().unwrap_or("?");
                format!(
                    r#"accent({},"{}")"#,
                    sources[0],
                    Self::escaped_typst_string(accent)
                )
            }
            "underline" | "overline" | "cancel" | "mid" | "scripts" | "limits" | "stretch"
                if sources.len() == 1 =>
            {
                call(head)
            }
            "underbrace" | "overbrace" | "underbracket" | "overbracket" | "underparen"
            | "overparen" | "undershell" | "overshell"
                if (1..=2).contains(&sources.len()) =>
            {
                call(head)
            }
            "vec" => call("vec"),
            "mat" => {
                let rows = arguments
                    .iter()
                    .map(|row| match row {
                        Self::List(entries) => entries
                            .iter()
                            .map(Self::to_typst_source)
                            .collect::<Vec<_>>()
                            .join(","),
                        value => value.to_typst_source(),
                    })
                    .collect::<Vec<_>>()
                    .join(";");
                format!("mat({rows})")
            }
            "cases" => call("cases"),
            "class" if arguments.len() == 2 => {
                let class = arguments[0].literal_text().unwrap_or("normal");
                format!(
                    r#"class("{}",{})"#,
                    Self::escaped_typst_string(class),
                    sources[1]
                )
            }
            // The constructor rejects unknown nodes. This fixed fallback also
            // keeps forged/imported user data from becoming executable source.
            _ => r#"upright("?")"#.to_owned(),
        }
    }

    fn is_typst_math_name(name: &str) -> bool {
        matches!(
            name,
            "alpha"
                | "beta"
                | "gamma"
                | "delta"
                | "epsilon"
                | "zeta"
                | "eta"
                | "theta"
                | "iota"
                | "kappa"
                | "lambda"
                | "mu"
                | "nu"
                | "xi"
                | "omicron"
                | "pi"
                | "rho"
                | "sigma"
                | "tau"
                | "upsilon"
                | "phi"
                | "chi"
                | "psi"
                | "omega"
        )
    }

    /// Render escaped Typst math source. See [`Self::math`] for an example.
    #[must_use]
    pub fn to_typst_source(&self) -> String {
        match self {
            Self::Symbol(name)
                if name.chars().count() == 1
                    && name
                        .chars()
                        .next()
                        .is_some_and(|character| character.is_alphanumeric()) =>
            {
                name.clone()
            }
            Self::Symbol(name) if Self::is_typst_math_name(name) => name.clone(),
            Self::Symbol(name) => {
                format!(r#"italic("{}")"#, Self::escaped_typst_string(name))
            }
            Self::Text(value) => {
                format!(r#"upright("{}")"#, Self::escaped_typst_string(value))
            }
            Self::Number(number) => number.to_string(),
            Self::Sequence(items) => items
                .iter()
                .map(Self::to_typst_source)
                .collect::<Vec<_>>()
                .join(" "),
            Self::List(items) => items
                .iter()
                .map(Self::to_typst_source)
                .collect::<Vec<_>>()
                .join(","),
            Self::Math { head, arguments } => Self::typst_math_source(head, arguments),
            Self::Attach { base, top, bottom } => {
                let mut source = format!("attach({}", base.to_typst_source());
                if let Some(top) = top {
                    source.push_str(",t:");
                    source.push_str(&top.to_typst_source());
                }
                if let Some(bottom) = bottom {
                    source.push_str(",b:");
                    source.push_str(&bottom.to_typst_source());
                }
                source.push(')');
                source
            }
        }
    }

    /// Render the native structured label with LaTeX metacharacters escaped.
    ///
    /// ```
    /// use symbolica_notation::IndexDisplay;
    /// let label = IndexDisplay::text("A&B")?;
    /// assert_eq!(label.to_latex_source(), "A\\&B");
    /// # Ok::<(), symbolica_notation::IndexDisplayError>(())
    /// ```
    #[must_use]
    pub fn to_latex_source(&self) -> String {
        DisplayText::escape_latex(&self.to_native_string())
    }

    /// Render the plain label with explicit script delimiters. See [`Self::text`].
    #[must_use]
    pub fn to_native_string(&self) -> String {
        match self {
            Self::Symbol(name) => name.clone(),
            Self::Text(value) => value.clone(),
            Self::Number(number) => number.to_string(),
            Self::Sequence(items) => items
                .iter()
                .map(Self::to_native_string)
                .collect::<Vec<_>>()
                .join(" "),
            Self::List(items) => items
                .iter()
                .map(Self::to_native_string)
                .collect::<Vec<_>>()
                .join(","),
            Self::Math { .. } => self.to_typst_source(),
            Self::Attach { base, top, bottom } => {
                let mut output = base.to_native_string();
                if let Some(top) = top {
                    output.push_str("^(");
                    output.push_str(&top.to_native_string());
                    output.push(')');
                }
                if let Some(bottom) = bottom {
                    output.push_str("_(");
                    output.push_str(&bottom.to_native_string());
                    output.push(')');
                }
                output
            }
        }
    }

    /// Render a complete MathML fragment, without a `math` root or source paths.
    ///
    /// This is a bounded notation offer: `None` declines unsupported or malformed
    /// display nodes, invalid XML characters, or trees exceeding depth 16, 256
    /// nodes, 64 children per sequence, or 4096 UTF-8 source/output bytes. Callers
    /// may then retain the index's ordinary symbolic label. Public enum variants
    /// can bypass constructors, so the borrowed tree is checked before allocating
    /// recursively rendered children.
    ///
    /// # Examples
    ///
    /// ```
    /// use symbolica_notation::{IndexDisplay, IndexDisplayError};
    ///
    /// let label = IndexDisplay::symbol("mu")?.with_bottom(IndexDisplay::Number(2));
    /// assert_eq!(label.to_mathml().as_deref(), Some("<msub><mi>μ</mi><mn>2</mn></msub>"));
    /// # Ok::<(), IndexDisplayError>(())
    /// ```
    pub fn to_mathml(&self) -> Option<String> {
        self.mathml_fits(0, &mut 256, &mut 4096)
            .then(|| self.mathml_source())?
    }

    fn mathml_fits(&self, depth: usize, nodes: &mut usize, bytes: &mut usize) -> bool {
        if depth > 16 || *nodes == 0 {
            return false;
        }
        *nodes -= 1;
        match self {
            Self::Symbol(value) | Self::Text(value) => {
                if value.len() > *bytes
                    || value
                        .chars()
                        .any(|c| c.is_control() || matches!(c, '\u{fffe}' | '\u{ffff}'))
                {
                    return false;
                }
                *bytes -= value.len();
                true
            }
            Self::Number(_) => true,
            Self::Sequence(items) | Self::List(items) => {
                items.len() <= 64
                    && items
                        .iter()
                        .all(|item| item.mathml_fits(depth + 1, nodes, bytes))
            }
            Self::Math { head, arguments } => {
                head.len() <= 16
                    && arguments.len() <= 64
                    && arguments
                        .iter()
                        .all(|item| item.mathml_fits(depth + 1, nodes, bytes))
            }
            Self::Attach { base, top, bottom } => {
                base.mathml_fits(depth + 1, nodes, bytes)
                    && top
                        .as_deref()
                        .is_none_or(|item| item.mathml_fits(depth + 1, nodes, bytes))
                    && bottom
                        .as_deref()
                        .is_none_or(|item| item.mathml_fits(depth + 1, nodes, bytes))
            }
        }
    }

    fn mathml_symbol(name: &str) -> &str {
        match name {
            "alpha" => "α",
            "beta" => "β",
            "gamma" => "γ",
            "delta" => "δ",
            "epsilon" => "ε",
            "zeta" => "ζ",
            "eta" => "η",
            "theta" => "θ",
            "iota" => "ι",
            "kappa" => "κ",
            "lambda" => "λ",
            "mu" => "μ",
            "nu" => "ν",
            "xi" => "ξ",
            "omicron" => "ο",
            "pi" => "π",
            "rho" => "ρ",
            "sigma" => "σ",
            "tau" => "τ",
            "upsilon" => "υ",
            "phi" => "φ",
            "chi" => "χ",
            "psi" => "ψ",
            "omega" => "ω",
            _ => name,
        }
    }

    fn mathml_group(source: &str) -> String {
        format!("<mrow><mo>(</mo>{source}<mo>)</mo></mrow>")
    }

    fn mathml_operand<'a>(&self, source: &'a str) -> Cow<'a, str> {
        if matches!(self, Self::Sequence(items) | Self::List(items) if items.len() > 1)
            || matches!(self, Self::Math { head, .. }
                if matches!(head.as_str(), "arg" | "add" | "sub" | "plus" | "neg" | "times" | "dot" | "mul"))
        {
            Cow::Owned(Self::mathml_group(source))
        } else {
            Cow::Borrowed(source)
        }
    }

    fn mathml_script_base<'a>(&self, source: &'a str) -> Cow<'a, str> {
        if matches!(self, Self::Attach { .. } | Self::Number(..=-1))
            || matches!(self, Self::Math { head, .. } if head == "pow")
        {
            Cow::Owned(Self::mathml_group(source))
        } else {
            self.mathml_operand(source)
        }
    }

    fn mathml_attachment(base: &Self, top: Option<&Self>, bottom: Option<&Self>) -> Option<String> {
        let source = base.mathml_source()?;
        let (tag, scripts) = match (top, bottom) {
            (Some(top), Some(bottom)) => (
                "msubsup",
                format!("{}{}", bottom.mathml_source()?, top.mathml_source()?),
            ),
            (Some(top), None) => ("msup", top.mathml_source()?),
            (None, Some(bottom)) => ("msub", bottom.mathml_source()?),
            (None, None) => return Some(source),
        };
        let base = base.mathml_script_base(&source);
        Some(format!("<{tag}>{base}{scripts}</{tag}>"))
    }

    fn mathml_table_row(&self) -> Option<String> {
        let cells = match self {
            Self::List(items) => items
                .iter()
                .map(Self::mathml_source)
                .collect::<Option<Vec<_>>>()?,
            value => vec![value.mathml_source()?],
        };
        Some(format!(
            "<mtr><mtd>{}</mtd></mtr>",
            cells.join("</mtd><mtd>")
        ))
    }

    fn mathml_source(&self) -> Option<String> {
        let source = match self {
            Self::Symbol(name) => format!("<mi>{}</mi>", MathMl::escape(Self::mathml_symbol(name))),
            Self::Text(value) => format!("<mtext>{}</mtext>", MathMl::escape(value)),
            Self::Number(number) => format!("<mn>{number}</mn>"),
            Self::Sequence(items) | Self::List(items) => {
                let sources = items
                    .iter()
                    .map(Self::mathml_source)
                    .collect::<Option<Vec<_>>>()?;
                let separator = if matches!(self, Self::List(_)) {
                    "<mo>,</mo>"
                } else {
                    ""
                };
                format!("<mrow>{}</mrow>", sources.join(separator))
            }
            Self::Attach { base, top, bottom } => {
                Self::mathml_attachment(base, top.as_deref(), bottom.as_deref())?
            }
            Self::Math { head, arguments } => Self::mathml_math_source(head, arguments)?,
        };
        (source.len() <= 4096).then_some(source)
    }

    fn mathml_math_source(head: &str, arguments: &[Self]) -> Option<String> {
        if head == "pow" {
            let [base, top] = arguments else { return None };
            return Self::mathml_attachment(base, Some(top), None);
        }
        if matches!(head, "vec" | "mat" | "cases") {
            let rows = if head == "vec" {
                arguments
                    .iter()
                    .map(|argument| {
                        Some(format!(
                            "<mtr><mtd>{}</mtd></mtr>",
                            argument.mathml_source()?
                        ))
                    })
                    .collect::<Option<Vec<_>>>()?
                    .join("")
            } else {
                arguments
                    .iter()
                    .map(Self::mathml_table_row)
                    .collect::<Option<Vec<_>>>()?
                    .join("")
            };
            return Some(if head == "cases" {
                format!("<mrow><mo>{{</mo><mtable columnalign=\"left\">{rows}</mtable></mrow>")
            } else {
                format!("<mrow><mo>(</mo><mtable>{rows}</mtable><mo>)</mo></mrow>")
            });
        }
        let sources = arguments
            .iter()
            .map(Self::mathml_source)
            .collect::<Option<Vec<_>>>()?;
        let joined = |separator: &str| format!("<mrow>{}</mrow>", sources.join(separator));
        let unary = || (arguments.len() == 1).then(|| arguments[0].mathml_operand(&sources[0]));
        Some(match head {
            "arg" => joined("<mo>,</mo>"),
            "add" => joined("<mo>+</mo>"),
            "sub" if !sources.is_empty() => {
                let mut terms = vec![Cow::Borrowed(sources[0].as_str())];
                terms.extend(
                    arguments[1..]
                        .iter()
                        .zip(&sources[1..])
                        .map(|(arg, source)| arg.mathml_operand(source)),
                );
                format!("<mrow>{}</mrow>", terms.join("<mo>−</mo>"))
            }
            "times" | "dot" | "mul" => {
                let terms = arguments
                    .iter()
                    .zip(&sources)
                    .map(|(arg, source)| arg.mathml_operand(source))
                    .collect::<Vec<_>>();
                let separator = match head {
                    "times" => "<mo>×</mo>",
                    "dot" => "<mo>⋅</mo>",
                    _ => "<mo>\u{2062}</mo>",
                };
                format!("<mrow>{}</mrow>", terms.join(separator))
            }
            "plus" => format!("<mrow><mo>+</mo>{}</mrow>", unary()?),
            "neg" => format!("<mrow><mo>−</mo>{}</mrow>", unary()?),
            "factorial" => format!("<mrow>{}<mo>!</mo></mrow>", unary()?),
            "()" | "group" | "lr" if sources.len() == 1 => Self::mathml_group(&sources[0]),
            "frac" if sources.len() == 2 => format!("<mfrac>{}{}</mfrac>", sources[0], sources[1]),
            "root" if sources.len() == 1 => format!("<msqrt>{}</msqrt>", sources[0]),
            "root" if sources.len() == 2 => format!("<mroot>{}{}</mroot>", sources[1], sources[0]),
            "call" | "op-call" if !sources.is_empty() => format!(
                "<mrow>{}<mo>\u{2061}</mo>{}</mrow>",
                sources[0],
                Self::mathml_group(&sources[1..].join("<mo>,</mo>"))
            ),
            "op" if sources.len() == 1 => format!(
                "<mi mathvariant=\"normal\">{}</mi>",
                MathMl::escape(arguments[0].literal_text()?)
            ),
            "accent" if sources.len() == 2 => {
                let accent = match arguments[1].literal_text()? {
                    "macron" | "bar" => "¯",
                    "hat" | "circumflex" => "^",
                    "tilde" => "~",
                    "arrow" | "arrow.r" | "vec" => "→",
                    "dot" => "˙",
                    "ddot" => "¨",
                    value => value,
                };
                format!(
                    "<mover accent=\"true\">{}<mo>{}</mo></mover>",
                    sources[0],
                    MathMl::escape(accent)
                )
            }
            "underline" | "overline" if sources.len() == 1 => {
                let (tag, attribute, mark) = if head == "underline" {
                    ("munder", "accentunder", "_")
                } else {
                    ("mover", "accent", "¯")
                };
                format!(
                    "<{tag} {attribute}=\"true\">{}<mo>{mark}</mo></{tag}>",
                    sources[0]
                )
            }
            "underbrace" | "overbrace" | "underbracket" | "overbracket" | "underparen"
            | "overparen" | "undershell" | "overshell"
                if (1..=2).contains(&sources.len()) =>
            {
                let (tag, attribute) = if head.starts_with("under") {
                    ("munder", "accentunder")
                } else {
                    ("mover", "accent")
                };
                let mark = match head {
                    "underbrace" => "⏟",
                    "overbrace" => "⏞",
                    "underbracket" => "⎵",
                    "overbracket" => "⎴",
                    "underparen" => "⏝",
                    "overparen" => "⏜",
                    "undershell" => "⏡",
                    _ => "⏠",
                };
                let decorated = format!(
                    "<{tag} {attribute}=\"true\">{}<mo>{mark}</mo></{tag}>",
                    sources[0]
                );
                if sources.len() == 2 {
                    format!("<{tag}>{decorated}{}</{tag}>", sources[1])
                } else {
                    decorated
                }
            }
            "cancel" if sources.len() == 1 => format!(
                "<menclose notation=\"updiagonalstrike\">{}</menclose>",
                sources[0]
            ),
            "class" if sources.len() == 2 && arguments[0].literal_text().is_some() => {
                format!("<mrow>{}</mrow>", sources[1])
            }
            "scripts" | "limits" if sources.len() == 1 => {
                if head == "limits"
                    && let Some((top, bottom)) = match &arguments[0] {
                        Self::Attach { top, bottom, .. } => Some((top.is_some(), bottom.is_some())),
                        Self::Math { head, arguments } if head == "pow" && arguments.len() == 2 => {
                            Some((true, false))
                        }
                        _ => None,
                    }
                {
                    let tags = match (top, bottom) {
                        (true, true) => Some(("msubsup", "munderover")),
                        (true, false) => Some(("msup", "mover")),
                        (false, true) => Some(("msub", "munder")),
                        (false, false) => None,
                    };
                    if let Some((old, new)) = tags {
                        let body = sources[0]
                            .strip_prefix(&format!("<{old}>"))?
                            .strip_suffix(&format!("</{old}>"))?;
                        format!("<{new}>{body}</{new}>")
                    } else {
                        sources.into_iter().next()?
                    }
                } else {
                    format!("<mrow>{}</mrow>", sources[0])
                }
            }
            "mid" | "stretch" if sources.len() == 1 => {
                if let Some(value) = arguments[0].literal_text() {
                    format!(
                        "<mo stretchy=\"true\">{}</mo>",
                        MathMl::escape(Self::mathml_symbol(value))
                    )
                } else {
                    format!("<mrow>{}</mrow>", sources[0])
                }
            }
            _ => return None,
        })
    }

    /// Encode one display node for embedding in existing portable metadata.
    /// See [`Self::from_node_user_data`] for a round-trip example.
    pub fn node_user_data(&self) -> UserData {
        match self {
            Self::Symbol(name) => UserData::List(vec![
                UserData::String("symbol".to_owned()),
                UserData::String(name.clone()),
            ]),
            Self::Text(value) => UserData::List(vec![
                UserData::String("text".to_owned()),
                UserData::String(value.clone()),
            ]),
            Self::Number(number) => UserData::List(vec![
                UserData::String("number".to_owned()),
                UserData::Integer(*number),
            ]),
            Self::Sequence(items) => UserData::List(vec![
                UserData::String("sequence".to_owned()),
                UserData::List(items.iter().map(Self::node_user_data).collect()),
            ]),
            Self::List(items) => UserData::List(vec![
                UserData::String("list".to_owned()),
                UserData::List(items.iter().map(Self::node_user_data).collect()),
            ]),
            Self::Math { head, arguments } => UserData::List(vec![
                UserData::String("math".to_owned()),
                UserData::String(head.clone()),
                UserData::List(arguments.iter().map(Self::node_user_data).collect()),
            ]),
            Self::Attach { base, top, bottom } => UserData::List(vec![
                UserData::String("attach".to_owned()),
                base.node_user_data(),
                top.as_deref()
                    .map(Self::node_user_data)
                    .unwrap_or(UserData::None),
                bottom
                    .as_deref()
                    .map(Self::node_user_data)
                    .unwrap_or(UserData::None),
            ]),
        }
    }

    /// Decode a bounded metadata node at the supplied nesting depth.
    /// Returns `None` for malformed nodes or exceeded bounds.
    ///
    /// ```
    /// use symbolica_notation::{IndexDisplay, IndexDisplayError};
    /// let label = IndexDisplay::symbol("mu")?;
    /// assert_eq!(IndexDisplay::from_node_user_data(&label.node_user_data(), 0), Some(label));
    /// # Ok::<(), IndexDisplayError>(())
    /// ```
    pub fn from_node_user_data(data: &UserData, depth: usize) -> Option<Self> {
        Self::parse_node_user_data(data, depth, &mut 256)
    }

    fn parse_node_user_data(data: &UserData, depth: usize, remaining: &mut usize) -> Option<Self> {
        if depth > 16 || *remaining == 0 {
            return None;
        }
        *remaining -= 1;
        let UserData::List(fields) = data else {
            return None;
        };
        match fields.as_slice() {
            [UserData::String(kind), UserData::String(name)] if kind == "symbol" => {
                if name.len() > 512 {
                    return None;
                }
                Self::symbol(name.clone()).ok()
            }
            [UserData::String(kind), UserData::String(value)] if kind == "text" => {
                if value.len() > 4096 {
                    return None;
                }
                Self::text(value.clone()).ok()
            }
            [UserData::String(kind), UserData::Integer(number)] if kind == "number" => {
                Some(Self::Number(*number))
            }
            [UserData::String(kind), UserData::List(items)] if kind == "sequence" => {
                if items.len() > 64 {
                    return None;
                }
                Some(Self::Sequence(
                    items
                        .iter()
                        .map(|item| Self::parse_node_user_data(item, depth + 1, remaining))
                        .collect::<Option<Vec<_>>>()?,
                ))
            }
            [UserData::String(kind), UserData::List(items)] if kind == "list" => {
                if items.len() > 64 {
                    return None;
                }
                Some(Self::List(
                    items
                        .iter()
                        .map(|item| Self::parse_node_user_data(item, depth + 1, remaining))
                        .collect::<Option<Vec<_>>>()?,
                ))
            }
            [
                UserData::String(kind),
                UserData::String(head),
                UserData::List(arguments),
            ] if kind == "math" => {
                if head.len() > 16 || arguments.len() > 64 {
                    return None;
                }
                let arguments = arguments
                    .iter()
                    .map(|argument| Self::parse_node_user_data(argument, depth + 1, remaining))
                    .collect::<Option<Vec<_>>>()?;
                Self::math(head.clone(), arguments).ok()
            }
            [UserData::String(kind), base, top, bottom] if kind == "attach" => Some(Self::Attach {
                base: Box::new(Self::parse_node_user_data(base, depth + 1, remaining)?),
                top: match top {
                    UserData::None => None,
                    value => Some(Box::new(Self::parse_node_user_data(
                        value,
                        depth + 1,
                        remaining,
                    )?)),
                },
                bottom: match bottom {
                    UserData::None => None,
                    value => Some(Box::new(Self::parse_node_user_data(
                        value,
                        depth + 1,
                        remaining,
                    )?)),
                },
            }),
            _ => None,
        }
    }

    /// Encode the unchanged `spenso::index-display-v1` symbol attachment.
    /// See [`Self::from_symbol`] for a round-trip example.
    pub fn symbol_user_data(&self) -> UserData {
        UserData::List(vec![
            UserData::String("spenso::index-display-v1".to_owned()),
            self.node_user_data(),
        ])
    }

    /// Decode an existing symbol attachment, declining malformed metadata.
    ///
    /// ```
    /// use std::error::Error;
    /// use symbolica::atom::SymbolBuilder;
    /// use symbolica::wrap_symbol;
    /// use symbolica_notation::IndexDisplay;
    /// let label = IndexDisplay::symbol("mu")?;
    /// let symbol = SymbolBuilder::new(wrap_symbol!("notation_doc::label"))
    ///     .with_user_data(label.symbol_user_data()).build().map_err(|error| error.to_string())?;
    /// assert_eq!(IndexDisplay::from_symbol(symbol), Some(label));
    /// # Ok::<(), Box<dyn Error>>(())
    /// ```
    pub fn from_symbol(symbol: Symbol) -> Option<Self> {
        match symbol.get_data() {
            UserData::List(fields) => match fields.as_slice() {
                [UserData::String(kind), display] if kind == "spenso::index-display-v1" => {
                    Self::from_node_user_data(display, 0)
                }
                _ => None,
            },
            _ => None,
        }
    }
}

#[cfg(test)]
mod index_mathml_tests {
    use super::IndexDisplay as D;

    #[test]
    fn mathml_escapes_literal_text_and_converts_greek_symbols() {
        assert_eq!(
            D::Symbol("<&>\"'".into()).to_mathml().unwrap(),
            "<mi>&lt;&amp;&gt;&quot;&#39;</mi>"
        );
        assert_eq!(
            D::Text("mu <&>\"'".into()).to_mathml().unwrap(),
            "<mtext>mu &lt;&amp;&gt;&quot;&#39;</mtext>"
        );
        assert_eq!(D::Symbol("mu".into()).to_mathml().unwrap(), "<mi>μ</mi>");
        assert_eq!(
            D::Symbol("Gamma".into()).to_mathml().unwrap(),
            "<mi>Gamma</mi>"
        );
        assert_eq!(D::Symbol("Mu".into()).to_mathml().unwrap(), "<mi>Mu</mi>");
        assert_eq!(D::Symbol("Γ".into()).to_mathml().unwrap(), "<mi>Γ</mi>");
        assert_eq!(D::Symbol("Μ".into()).to_mathml().unwrap(), "<mi>Μ</mi>");
        assert_eq!(D::Number(-12).to_mathml().unwrap(), "<mn>-12</mn>");
    }

    #[test]
    fn mathml_scripts_and_fractions_preserve_grouping() {
        let sum = D::math("add", vec![D::Symbol("x".into()), D::Symbol("y".into())]).unwrap();
        let display = sum
            .with_top(D::Number(2))
            .with_bottom(D::Symbol("mu".into()));
        assert_eq!(
            display.to_mathml().unwrap(),
            "<msubsup><mrow><mo>(</mo><mrow><mi>x</mi><mo>+</mo><mi>y</mi></mrow><mo>)</mo></mrow><mi>μ</mi><mn>2</mn></msubsup>"
        );
        let fraction = D::math("frac", vec![display, D::Number(3)]).unwrap();
        assert!(
            fraction
                .to_mathml()
                .unwrap()
                .starts_with("<mfrac><msubsup>")
        );
        assert!(
            fraction
                .to_mathml()
                .unwrap()
                .ends_with("</msubsup><mn>3</mn></mfrac>")
        );
        assert_eq!(
            D::math("root", vec![D::Number(3), D::Symbol("x".into())])
                .unwrap()
                .to_mathml()
                .unwrap(),
            "<mroot><mi>x</mi><mn>3</mn></mroot>"
        );
        assert_eq!(
            D::math("pow", vec![D::Number(-2), D::Number(2)])
                .unwrap()
                .to_mathml()
                .unwrap(),
            "<msup><mrow><mo>(</mo><mn>-2</mn><mo>)</mo></mrow><mn>2</mn></msup>"
        );
        let nested = D::math(
            "pow",
            vec![D::Symbol("x".into()).with_top(D::Number(2)), D::Number(3)],
        )
        .unwrap();
        assert_eq!(
            nested.to_mathml().unwrap(),
            "<msup><mrow><mo>(</mo><msup><mi>x</mi><mn>2</mn></msup><mo>)</mo></mrow><mn>3</mn></msup>"
        );
    }

    #[test]
    fn mathml_matrix_vector_and_cases_retain_cell_boundaries() {
        let rows = vec![
            D::List(vec![D::Symbol("alpha".into()), D::Number(1)]),
            D::List(vec![D::Number(0), D::Text("x<y".into())]),
        ];
        assert_eq!(
            D::math("mat", rows.clone()).unwrap().to_mathml().unwrap(),
            "<mrow><mo>(</mo><mtable><mtr><mtd><mi>α</mi></mtd><mtd><mn>1</mn></mtd></mtr><mtr><mtd><mn>0</mn></mtd><mtd><mtext>x&lt;y</mtext></mtd></mtr></mtable><mo>)</mo></mrow>"
        );
        let cases = D::math("cases", rows).unwrap().to_mathml().unwrap();
        assert!(cases.starts_with("<mrow><mo>{</mo><mtable columnalign=\"left\">"));
        assert!(!cases.contains("<mo>,</mo>"));
        assert_eq!(
            D::math("vec", vec![D::Number(1), D::Number(2)])
                .unwrap()
                .to_mathml()
                .unwrap(),
            "<mrow><mo>(</mo><mtable><mtr><mtd><mn>1</mn></mtd></mtr><mtr><mtd><mn>2</mn></mtd></mtr></mtable><mo>)</mo></mrow>"
        );
    }

    #[test]
    fn mathml_safe_display_heads_have_direct_renderings() {
        for head in [
            "arg",
            "add",
            "sub",
            "plus",
            "neg",
            "times",
            "dot",
            "factorial",
            "mul",
            "()",
            "group",
            "lr",
            "root",
            "call",
            "op-call",
            "op",
            "underline",
            "overline",
            "underbrace",
            "overbrace",
            "underbracket",
            "overbracket",
            "underparen",
            "overparen",
            "undershell",
            "overshell",
            "cancel",
            "mid",
            "scripts",
            "limits",
            "stretch",
        ] {
            assert!(
                D::math(head, vec![D::Symbol("x".into())])
                    .unwrap()
                    .to_mathml()
                    .is_some(),
                "{head}"
            );
        }
        for head in ["pow", "frac", "accent", "class"] {
            assert!(
                D::math(head, vec![D::Symbol("x".into()), D::Symbol("y".into())])
                    .unwrap()
                    .to_mathml()
                    .is_some(),
                "{head}"
            );
        }
        assert_eq!(
            D::math(
                "accent",
                vec![D::Symbol("p".into()), D::Symbol("vec".into())]
            )
            .unwrap()
            .to_mathml()
            .unwrap(),
            "<mover accent=\"true\"><mi>p</mi><mo>→</mo></mover>"
        );
        assert_eq!(
            D::math(
                "limits",
                vec![
                    D::Symbol("x".into())
                        .with_top(D::Number(2))
                        .with_bottom(D::Number(1))
                ]
            )
            .unwrap()
            .to_mathml()
            .unwrap(),
            "<munderover><mi>x</mi><mn>1</mn><mn>2</mn></munderover>"
        );
    }

    #[test]
    fn mathml_limits_treat_power_and_attachment_spellings_identically() {
        let power = D::math("pow", vec![D::Symbol("x".into()), D::Number(2)]).unwrap();
        let attachment = D::Symbol("x".into()).with_top(D::Number(2));
        assert_eq!(power.to_mathml(), attachment.to_mathml());
        let power_limits = D::math("limits", vec![power]).unwrap();
        let attachment_limits = D::math("limits", vec![attachment]).unwrap();
        assert_eq!(
            power_limits.to_typst_source(),
            attachment_limits.to_typst_source()
        );
        assert_eq!(power_limits.to_mathml(), attachment_limits.to_mathml());
        assert_eq!(
            power_limits.to_mathml().unwrap(),
            "<mover><mi>x</mi><mn>2</mn></mover>"
        );
    }

    #[test]
    fn mathml_rejects_forged_unsupported_and_unbounded_trees() {
        assert!(
            D::Math {
                head: "raw-code".into(),
                arguments: vec![]
            }
            .to_mathml()
            .is_none()
        );
        assert!(
            D::Math {
                head: "frac".into(),
                arguments: vec![D::Number(1)]
            }
            .to_mathml()
            .is_none()
        );
        assert!(D::Sequence(vec![D::Number(1); 65]).to_mathml().is_none());
        let mut deep = D::Number(1);
        for _ in 0..17 {
            deep = D::Sequence(vec![deep]);
        }
        assert!(deep.to_mathml().is_none());
        assert!(
            D::Sequence(vec![D::Sequence(vec![D::Number(1); 64]); 5])
                .to_mathml()
                .is_none()
        );
        assert!(D::Text("é".repeat(2040)).to_mathml().is_some());
        assert!(D::Text("é".repeat(2041)).to_mathml().is_none());
        assert_eq!(D::Text("x".repeat(4081)).to_mathml().unwrap().len(), 4096);
        assert!(D::Text("x".repeat(4082)).to_mathml().is_none());
        assert!(D::Text("&".repeat(817)).to_mathml().is_none());
        assert!(D::Text("\u{0}".into()).to_mathml().is_none());
        assert!(D::Text("\u{ffff}".into()).to_mathml().is_none());
    }
}

#[cfg(test)]
mod metadata_tests {
    use super::IndexDisplay;
    use symbolica::atom::{NamespacedSymbol, SymbolBuilder};
    #[test]
    fn mathematical_display_nodes_render_and_round_trip_as_safe_user_data() {
        let display = IndexDisplay::math(
            "add",
            vec![
                IndexDisplay::math(
                    "accent",
                    vec![
                        IndexDisplay::symbol("p").unwrap(),
                        IndexDisplay::symbol("⃗").unwrap(),
                    ],
                )
                .unwrap(),
                IndexDisplay::math(
                    "frac",
                    vec![IndexDisplay::symbol("q").unwrap(), IndexDisplay::Number(2)],
                )
                .unwrap(),
                IndexDisplay::text("soft").unwrap(),
            ],
        )
        .unwrap();

        assert_eq!(
            display.to_typst_source(),
            r#"accent(p,"⃗") + frac(q,2) + upright("soft")"#
        );
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(
            "spenso_rep_metadata_tests::MathDisplayRoundTrip",
        ))
        .with_user_data(display.symbol_user_data())
        .build()
        .unwrap();
        assert_eq!(IndexDisplay::from_symbol(symbol), Some(display));
    }

    #[test]
    fn unknown_math_display_heads_cannot_be_constructed() {
        assert!(IndexDisplay::math("raw-code", vec![]).is_err());
    }

    #[test]
    fn mathematical_display_metadata_parsing_has_a_total_node_budget() {
        let wide = IndexDisplay::Sequence(vec![
            IndexDisplay::Sequence(vec![
                IndexDisplay::Number(1);
                64
            ]);
            5
        ]);
        assert!(IndexDisplay::from_node_user_data(&wide.node_user_data(), 0).is_none());
        let too_many_arguments = IndexDisplay::Math {
            head: "add".into(),
            arguments: vec![IndexDisplay::Number(1); 65],
        };
        assert!(
            IndexDisplay::from_node_user_data(&too_many_arguments.node_user_data(), 0).is_none()
        );
    }
}
