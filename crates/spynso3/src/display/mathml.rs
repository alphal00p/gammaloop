//! Let CSS wrap whole additive terms, without changing their internal math layout.

/// Typst emits XML-compatible MathML inside an HTML document. Preserve the HTML
/// and exact MathML source, inserting only presentation groups at safe breaks.
pub(super) fn wrap_sums(fragment: &str) -> String {
    let mut output = String::with_capacity(fragment.len());
    let mut rest = fragment;
    while let Some(start) = rest.find("<math") {
        let Some(end) = rest[start..].find("</math>") else {
            break;
        };
        let end = start + end + "</math>".len();
        output.push_str(&rest[..start]);
        let source = &rest[start..end];
        if let Some(wrapped) = wrap_sum(source) {
            output.push_str(&wrapped);
        } else {
            output.push_str(source);
        }
        rest = &rest[end..];
    }
    output.push_str(rest);
    output
}

fn wrap_sum(source: &str) -> Option<String> {
    // Unknown/non-XML fragments keep their original presentation. This should
    // never prevent a custom notation from being displayed.
    let document = roxmltree::Document::parse(source).ok()?;
    let math = document.root_element();
    if math.tag_name().name() != "math"
        || math.attribute("display") != Some("block")
        || math.has_attribute("data-spenso-wrap")
    {
        return None;
    }
    let children: Vec<_> = math.children().filter(|node| node.is_element()).collect();
    let mut breaks = Vec::new();
    let mut depth = 0usize;
    for (index, node) in children.iter().enumerate() {
        if node.tag_name().name() != "mo" {
            continue;
        }
        let text = node.text().unwrap_or_default().trim();
        match text {
            "(" | "[" | "{" | "⟨" | "⌈" | "⌊" => depth += 1,
            ")" | "]" | "}" | "⟩" | "⌉" | "⌋" => depth = depth.checked_sub(1)?,
            "|" | "‖" => return None, // Ambiguous fences: keep them together.
            "+" | "−" | "-" if depth == 0 && index > 0 && index + 1 < children.len() => {
                let previous = children[index - 1];
                if node.attribute("form") != Some("prefix")
                    && (previous.tag_name().name() != "mo"
                        || matches!(previous.text(), Some(")" | "]" | "}" | "⟩" | "⌉" | "⌋")))
                {
                    breaks.push(*node);
                }
            }
            _ if node.attribute("fence") == Some("true") => return None,
            _ => {}
        }
    }
    if breaks.is_empty() || depth != 0 {
        return None;
    }

    let mut output = String::with_capacity(source.len() + 80 * breaks.len());
    output.push_str("<math data-spenso-wrap=\"\"");
    let mut position = "<math".len();
    let first = children.first()?.range().start;
    output.push_str(&source[position..first]);
    output.push_str("<mrow data-spenso-term=\"\">");
    position = first;
    for operator in breaks {
        let start = operator.range().start;
        output.push_str(&source[position..start]);
        output.push_str("</mrow><mrow data-spenso-term=\"\">");
        position = start;
        // Grouping moves this binary sign to the beginning of an mrow. Keep
        // its original infix spacing instead of turning it into a unary sign.
        if !operator.has_attribute("form") {
            output.push_str("<mo form=\"infix\"");
            position += "<mo".len();
        }
    }
    let end = children.last()?.range().end;
    output.push_str(&source[position..end]);
    output.push_str("</mrow>");
    output.push_str(&source[end..]);
    Some(output)
}

#[cfg(test)]
mod tests {
    use super::wrap_sums;

    #[test]
    fn groups_only_top_level_binary_signs_and_preserves_nested_math() {
        let fraction = "<mfrac><mi>a</mi><mrow><mi>b</mi><mo>+</mo><mi>c</mi></mrow></mfrac>";
        let html = format!(
            "<style>.x{{}}</style><math display=\"block\"><mo>−</mo>{fraction}<mo>+</mo><msub><mi>T</mi><mi>i</mi></msub><mo>−</mo><mi>z</mi></math>"
        );
        let wrapped = wrap_sums(&html);
        assert!(wrapped.starts_with("<style>.x{}</style>"));
        assert!(wrapped.contains(fraction));
        assert_eq!(wrapped.matches("data-spenso-term").count(), 3);
        assert_eq!(wrapped.matches("form=\"infix\"").count(), 2);
        roxmltree::Document::parse(&wrapped[wrapped.find("<math").unwrap()..]).unwrap();
        assert_eq!(wrap_sums(&wrapped), wrapped);
    }

    #[test]
    fn leaves_fenced_sums_products_and_inline_math_intact() {
        for html in [
            "<math display=\"block\"><mo>(</mo><mi>a</mi><mo>+</mo><mi>b</mi><mo>)</mo><mi>c</mi></math>",
            "<math display=\"block\"><mi>a</mi><mi>b</mi></math>",
            "<math><mi>a</mi><mo>+</mo><mi>b</mi></math>",
        ] {
            assert_eq!(wrap_sums(html), html);
        }
    }
}
