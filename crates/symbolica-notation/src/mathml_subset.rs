use crate::MathMl;
use roxmltree::{Document, Node};
use thiserror::Error;

/// A bounded, validated presentation MathML fragment with no executable markup.
///
/// Unlike [`crate::IndexDisplay`], this preserves presentation attributes and
/// arbitrary combinations of the allowed MathML layout elements.
///
/// ```
/// use symbolica_notation::MathMlSubset;
/// let head = MathMlSubset::parse("<msub><mi>A</mi><mn>2</mn></msub>")?;
/// assert!(head.to_mathml().contains("<msub>"));
/// assert!(MathMlSubset::parse("<img src='x' onerror='alert(1)'/>").is_err());
/// # Ok::<(), symbolica_notation::MathMlSubsetError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MathMlSubset(Vec<PresentationNode>);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum PresentationNode {
    Text(String),
    Element {
        name: String,
        attributes: Vec<(String, String)>,
        children: Vec<Self>,
    },
}

/// A rejected MathML element, attribute, value, syntax, or resource bound.
/// See [`MathMlSubset`] for examples.
#[derive(Debug, Error)]
pub enum MathMlSubsetError {
    /// Disallowed presentation markup or a resource bound. See [`MathMlSubset`].
    #[error("invalid MathML print mapping: {0}")]
    Invalid(String),
    /// Malformed XML, including unknown entities. See [`MathMlSubset`].
    #[error("invalid MathML print mapping: {0}")]
    Syntax(#[from] roxmltree::Error),
}

impl MathMlSubset {
    /// Parse at most 4096 bytes, 256 nodes, and 16 levels of presentation MathML.
    ///
    /// Namespaces, comments, processing instructions, CDATA, declarations, and
    /// non-XML entities are rejected. Output is reconstructed from validated
    /// nodes and escaped text, never copied from the input markup.
    ///
    /// # Errors
    /// Returns [`MathMlSubsetError`] naming the rejected construct or bound.
    /// See [`Self`] for examples.
    pub fn parse(source: &str) -> Result<Self, MathMlSubsetError> {
        if source.len() > 4096 {
            return Err(MathMlSubsetError::Invalid(
                "source exceeds 4096 bytes".to_owned(),
            ));
        }
        if source.contains("<!") || source.contains("<?") {
            return Err(MathMlSubsetError::Invalid(
                "declarations, comments, processing instructions and CDATA are forbidden"
                    .to_owned(),
            ));
        }
        // Identify unsupported element names before XML syntax validation, so
        // malformed HTML such as an unquoted img attribute still names img in
        // the Python error. The XML parser remains authoritative for syntax.
        for opening in source.split('<').skip(1) {
            let opening = opening.strip_prefix('/').unwrap_or(opening);
            let name = opening
                .split(|character: char| {
                    character.is_whitespace() || matches!(character, '>' | '/')
                })
                .next()
                .unwrap_or_default();
            if !PresentationNode::allowed_element(name) {
                return Err(MathMlSubsetError::Invalid(format!(
                    "forbidden element {name}; namespace prefixes are not allowed"
                )));
            }
        }
        let wrapped = format!("<mrow>{source}</mrow>");
        let document = Document::parse(&wrapped)?;
        let mut remaining = 256;
        let nodes = document
            .root_element()
            .children()
            .map(|node| PresentationNode::parse(node, 0, &mut remaining))
            .collect::<Result<Vec<_>, _>>()?;
        let value = Self(nodes);
        if value.to_mathml().len() > 16384 {
            return Err(MathMlSubsetError::Invalid(
                "escaped output exceeds 16384 bytes".to_owned(),
            ));
        }
        Ok(value)
    }

    /// Whether raising this head to a power requires explicit parentheses.
    /// Transparent rows are inspected structurally, including heads
    /// carrying both an upper and a lower script.
    ///
    /// ```
    /// use symbolica_notation::MathMlSubset;
    /// let head = MathMlSubset::parse("<msubsup><mi>T</mi><mi>i</mi><mi>j</mi></msubsup>")?;
    /// assert!(head.needs_power_fence());
    /// let accent = MathMlSubset::parse("<mover><mi>ϵ</mi><mo>¯</mo></mover>")?;
    /// assert!(!accent.needs_power_fence());
    /// # Ok::<(), symbolica_notation::MathMlSubsetError>(())
    /// ```
    #[must_use]
    pub fn needs_power_fence(&self) -> bool {
        PresentationNode::needs_power_fence(&self.0)
    }

    /// Reconstruct escaped MathML, grouped into one script base.
    /// See [`Self`] for an example.
    #[must_use]
    pub fn to_mathml(&self) -> String {
        let mut output = String::from("<mrow>");
        for node in &self.0 {
            node.write(&mut output);
        }
        output.push_str("</mrow>");
        output
    }
}

impl PresentationNode {
    fn needs_power_fence(nodes: &[Self]) -> bool {
        let mut nodes = nodes
            .iter()
            .filter(|node| !matches!(node, Self::Text(text) if text.trim().is_empty()));
        let Some(node) = nodes.next() else {
            return false;
        };
        if nodes.next().is_some() {
            return true;
        }
        match node {
            Self::Element { name, children, .. } if name == "mrow" => {
                Self::needs_power_fence(children)
            }
            Self::Element { name, .. } => {
                matches!(name.as_str(), "msup" | "msubsup" | "munderover")
            }
            Self::Text(_) => false,
        }
    }

    fn allowed_element(name: &str) -> bool {
        matches!(
            name,
            "mi" | "mn"
                | "mo"
                | "mtext"
                | "mrow"
                | "msub"
                | "msup"
                | "msubsup"
                | "mover"
                | "munder"
                | "munderover"
                | "mfrac"
                | "msqrt"
                | "mroot"
                | "mspace"
                | "menclose"
        )
    }

    fn parse(
        node: Node<'_, '_>,
        depth: usize,
        remaining: &mut usize,
    ) -> Result<Self, MathMlSubsetError> {
        if depth > 16 || *remaining == 0 {
            return Err(MathMlSubsetError::Invalid(
                "depth or node count exceeds 16/256".to_owned(),
            ));
        }
        *remaining -= 1;
        if node.is_text() {
            let text = node.text().unwrap_or_default();
            // XML parsing has validated characters, including XML whitespace.
            return Ok(Self::Text(text.to_owned()));
        }
        let name = node.tag_name().name();
        if !node.is_element()
            || node.tag_name().namespace().is_some()
            || node.namespaces().len() != 0
        {
            return Err(MathMlSubsetError::Invalid(format!(
                "forbidden namespace or node on element {name}"
            )));
        }
        if !Self::allowed_element(name) {
            return Err(MathMlSubsetError::Invalid(format!(
                "forbidden element {name}"
            )));
        }
        let mut attributes = Vec::new();
        for attribute in node.attributes() {
            let key = attribute.name();
            let value = attribute.value();
            let allowed = attribute.namespace().is_none()
                && match key {
                    "mathvariant" => matches!(
                        value,
                        "normal"
                            | "bold"
                            | "italic"
                            | "bold-italic"
                            | "double-struck"
                            | "bold-fraktur"
                            | "script"
                            | "bold-script"
                            | "fraktur"
                            | "sans-serif"
                            | "bold-sans-serif"
                            | "sans-serif-italic"
                            | "sans-serif-bold-italic"
                            | "monospace"
                            | "initial"
                            | "tailed"
                            | "looped"
                            | "stretched"
                    ),
                    "accent" => {
                        matches!(name, "mover" | "munderover") && matches!(value, "true" | "false")
                    }
                    "accentunder" => {
                        matches!(name, "munder" | "munderover") && matches!(value, "true" | "false")
                    }
                    "stretchy" => name == "mo" && matches!(value, "true" | "false"),
                    "width" => name == "mspace" && Self::valid_width(value),
                    "notation" => {
                        name == "menclose"
                            && !value.is_empty()
                            && value.split_whitespace().all(|word| {
                                matches!(
                                    word,
                                    "longdiv"
                                        | "actuarial"
                                        | "radical"
                                        | "box"
                                        | "roundedbox"
                                        | "circle"
                                        | "left"
                                        | "right"
                                        | "top"
                                        | "bottom"
                                        | "updiagonalstrike"
                                        | "downdiagonalstrike"
                                        | "verticalstrike"
                                        | "horizontalstrike"
                                        | "madruwb"
                                        | "phasorangle"
                                        | "updiagonalarrow"
                                )
                            })
                    }
                    "mathcolor" => {
                        !value.is_empty()
                            && value.len() <= 32
                            && (value.bytes().all(|byte| byte.is_ascii_alphabetic())
                                || value.strip_prefix('#').is_some_and(|hex| {
                                    matches!(hex.len(), 3 | 6)
                                        && hex.bytes().all(|byte| byte.is_ascii_hexdigit())
                                }))
                    }
                    _ => false,
                };
            if !allowed {
                return Err(MathMlSubsetError::Invalid(format!(
                    "forbidden attribute {key}={value:?} on element {name}"
                )));
            }
            attributes.push((key.to_owned(), value.to_owned()));
        }
        attributes.sort();
        let children = node
            .children()
            .map(|child| Self::parse(child, depth + 1, remaining))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self::Element {
            name: name.to_owned(),
            attributes,
            children,
        })
    }

    fn valid_width(value: &str) -> bool {
        ["em", "ex", "px", "pt", "pc", "in", "cm", "mm", "%"]
            .into_iter()
            .any(|unit| {
                value.strip_suffix(unit).is_some_and(|number| {
                    !number.is_empty()
                        && number
                            .bytes()
                            .all(|byte| byte.is_ascii_digit() || byte == b'.')
                        && number.parse::<f64>().is_ok_and(|number| {
                            number.is_finite() && (0.0..=10000.0).contains(&number)
                        })
                })
            })
    }

    fn write(&self, output: &mut String) {
        match self {
            Self::Text(text) => output.push_str(&MathMl::escape(text)),
            Self::Element {
                name,
                attributes,
                children,
            } => {
                output.push('<');
                output.push_str(name);
                for (key, value) in attributes {
                    output.push(' ');
                    output.push_str(key);
                    output.push_str("=\"");
                    output.push_str(&MathMl::escape(value));
                    output.push('"');
                }
                output.push('>');
                for child in children {
                    child.write(output);
                }
                output.push_str("</");
                output.push_str(name);
                output.push('>');
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allows_layout_and_reconstructs_escaped_text() {
        assert!(MathMlSubset::parse("<mrow>\n  <mi>A</mi>\n</mrow>").is_ok());
        let parsed = MathMlSubset::parse(
            "<mover accent='true'><mi mathvariant='bold'>A&amp;B</mi><mo>¯</mo></mover>",
        )
        .unwrap();
        assert_eq!(
            parsed.to_mathml(),
            "<mrow><mover accent=\"true\"><mi mathvariant=\"bold\">A&amp;B</mi><mo>¯</mo></mover></mrow>"
        );
        assert_eq!(
            MathMlSubset::parse(&parsed.to_mathml())
                .unwrap()
                .to_mathml(),
            format!("<mrow>{}</mrow>", parsed.to_mathml())
        );
    }

    #[test]
    fn rejects_executable_and_opaque_xml() {
        for (source, offending) in [
            ("<img src='x' onerror='alert(1)'/>", "img"),
            ("<img src=x onerror=boom>", "img"),
            ("<mi onclick='x'>A</mi>", "onclick"),
            ("<mi href='x'>A</mi>", "href"),
            ("<mi style='color:red'>A</mi>", "style"),
            ("<mi class='x'>A</mi>", "class"),
            ("<mi id='x'>A</mi>", "id"),
            ("<mi src='x'>A</mi>", "src"),
            ("<mi mathvariant='invalid'>A</mi>", "mathvariant"),
            ("<mspace width='url(x)'/>", "width"),
            ("<x:mi xmlns:x='urn:x'>A</x:mi>", "namespace"),
            ("<!--x--><mi>A</mi>", "comments"),
            ("<?x y?><mi>A</mi>", "processing"),
            ("<![CDATA[<mi>A</mi>]]>", "CDATA"),
            ("<mi>&unknown;</mi>", "unknown"),
        ] {
            let error = MathMlSubset::parse(source).unwrap_err().to_string();
            assert!(error.contains(offending), "{source}: {error}");
        }
        assert!(MathMlSubset::parse(&"x".repeat(4097)).is_err());
        assert!(
            MathMlSubset::parse(&format!("{}x{}", "<mrow>".repeat(18), "</mrow>".repeat(18)))
                .is_err()
        );
        assert!(MathMlSubset::parse(&"<mi>x</mi>".repeat(130)).is_err());
    }
}
