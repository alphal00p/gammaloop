pub(crate) fn escape_html(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            _ => escaped.push(character),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::escape_html;

    #[test]
    fn escapes_html_metadata() {
        assert_eq!(
            escape_html("<script data-x='a&b'>\"x\"</script>"),
            "&lt;script data-x=&#39;a&amp;b&#39;&gt;&quot;x&quot;&lt;/script&gt;"
        );
    }
}
