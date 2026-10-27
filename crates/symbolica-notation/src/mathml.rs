/// Shared escaping and bounded embedding of MathML fragments.
///
/// ```
/// use symbolica_notation::MathMl;
/// assert_eq!(MathMl::escape("<&"), "&lt;&amp;");
/// assert_eq!(MathMl::fragment("<mi>x</mi>".into(), 4096).as_deref(), Some("<mi>x</mi>"));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MathMl;

impl MathMl {
    /// Escape XML text and attributes. See the [`MathMl`] example.
    #[must_use]
    pub fn escape(value: &str) -> String {
        let mut output = String::new();
        for character in value.chars() {
            match character {
                '&' => output.push_str("&amp;"),
                '<' => output.push_str("&lt;"),
                '>' => output.push_str("&gt;"),
                '"' => output.push_str("&quot;"),
                '\'' => output.push_str("&#39;"),
                character => output.push(character),
            }
        }
        output
    }

    /// Strip Symbolica viewer attributes from a bounded, complete fragment.
    ///
    /// Returns `None` for oversized, omitted, root-wrapped or malformed input.
    /// This is metadata removal, not sanitization of untrusted HTML.
    /// See the [`MathMl`] example.
    pub fn fragment(mut value: String, max_bytes: usize) -> Option<String> {
        if value.len() > max_bytes || value.contains("<mtext>⋯</mtext>") || value.contains("<math")
        {
            return None;
        }
        let bytes = value.as_bytes();
        let mut cursor = 0;
        let mut attributes = Vec::new();
        while let Some(start) = value[cursor..].find('<').map(|offset| cursor + offset) {
            let tail = &value[start..];
            if tail.starts_with("<!--") {
                cursor = start + tail.find("-->")? + 3;
                continue;
            }
            if tail.starts_with("<![CDATA[") {
                cursor = start + tail.find("]]>")? + 3;
                continue;
            }
            if tail.starts_with("<?") {
                cursor = start + tail.find("?>")? + 2;
                continue;
            }
            cursor = start + 1;
            if matches!(bytes.get(cursor), Some(b'/' | b'!')) {
                cursor = start + tail.find('>')? + 1;
                continue;
            }
            // Token text is opaque. Only names at opening-tag attribute boundaries
            // can own viewer metadata; quoted values may themselves contain `>`.
            while let Some(byte) = bytes.get(cursor)
                && !byte.is_ascii_whitespace()
                && !matches!(byte, b'/' | b'>')
            {
                cursor += 1;
            }
            loop {
                let whitespace = cursor;
                while bytes.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                    cursor += 1;
                }
                if bytes.get(cursor) == Some(&b'>') {
                    cursor += 1;
                    break;
                }
                if bytes.get(cursor..cursor + 2) == Some(b"/>".as_slice()) {
                    cursor += 2;
                    break;
                }
                let name = cursor;
                while let Some(byte) = bytes.get(cursor)
                    && !byte.is_ascii_whitespace()
                    && !matches!(byte, b'=' | b'/' | b'>')
                {
                    cursor += 1;
                }
                if cursor == name {
                    return None;
                }
                let symbolic = value[name..cursor].starts_with("data-symbolica-");
                while bytes.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                    cursor += 1;
                }
                if bytes.get(cursor) != Some(&b'=') {
                    return None;
                }
                cursor += 1;
                while bytes.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                    cursor += 1;
                }
                let quote = *bytes.get(cursor)?;
                if !matches!(quote, b'\'' | b'"') {
                    return None;
                }
                cursor += 1;
                cursor += value[cursor..].find(char::from(quote))? + 1;
                if symbolic {
                    attributes.push(whitespace..cursor);
                }
            }
        }
        for attribute in attributes.into_iter().rev() {
            value.replace_range(attribute, "");
        }
        value.trim_start().starts_with('<').then_some(value)
    }
}
