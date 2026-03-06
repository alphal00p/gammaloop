#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CommandLineParseError {
    IncompleteShellSyntax,
}

pub fn split_command_line(input: &str) -> Result<Vec<String>, CommandLineParseError> {
    let normalized = escape_process_ref_hash_ids(input);
    shlex::split(&normalized).ok_or(CommandLineParseError::IncompleteShellSyntax)
}

fn escape_process_ref_hash_ids(input: &str) -> String {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Mode {
        Unquoted,
        SingleQuoted,
        DoubleQuoted,
    }

    if !input.contains('#') {
        return input.to_string();
    }

    let mut escaped = String::with_capacity(input.len());
    let mut chars = input.chars().peekable();
    let mut mode = Mode::Unquoted;
    let mut at_word_start = true;
    let mut unquoted_escape = false;
    let mut double_quoted_escape = false;

    while let Some(ch) = chars.next() {
        match mode {
            Mode::Unquoted => {
                if unquoted_escape {
                    escaped.push(ch);
                    unquoted_escape = false;
                    at_word_start = false;
                    continue;
                }

                match ch {
                    '\\' => {
                        escaped.push(ch);
                        unquoted_escape = true;
                        at_word_start = false;
                    }
                    '\'' => {
                        escaped.push(ch);
                        mode = Mode::SingleQuoted;
                        at_word_start = false;
                    }
                    '"' => {
                        escaped.push(ch);
                        mode = Mode::DoubleQuoted;
                        at_word_start = false;
                    }
                    c if c.is_whitespace() => {
                        escaped.push(c);
                        at_word_start = true;
                    }
                    '#' if at_word_start
                        && chars.peek().is_some_and(|next| next.is_ascii_digit()) =>
                    {
                        escaped.push('\\');
                        escaped.push('#');
                        at_word_start = false;
                    }
                    _ => {
                        escaped.push(ch);
                        at_word_start = false;
                    }
                }
            }
            Mode::SingleQuoted => {
                escaped.push(ch);
                if ch == '\'' {
                    mode = Mode::Unquoted;
                }
            }
            Mode::DoubleQuoted => {
                escaped.push(ch);
                if double_quoted_escape {
                    double_quoted_escape = false;
                    continue;
                }
                match ch {
                    '\\' => double_quoted_escape = true,
                    '"' => mode = Mode::Unquoted,
                    _ => {}
                }
            }
        }
    }

    escaped
}

#[cfg(test)]
mod test {
    use super::split_command_line;

    #[test]
    fn split_supports_multiline_quoted_values() {
        let line = "set process -p epem_a_tth -i LO string '[integrator]\nn_start = 1000\n'";
        let parts = split_command_line(line).unwrap();
        assert_eq!(
            parts,
            vec![
                "set",
                "process",
                "-p",
                "epem_a_tth",
                "-i",
                "LO",
                "string",
                "[integrator]\nn_start = 1000\n",
            ]
        );
    }

    #[test]
    fn split_preserves_hash_process_references() {
        let line = "display integrands -p #12";
        let parts = split_command_line(line).unwrap();
        assert_eq!(parts, vec!["display", "integrands", "-p", "#12"]);
    }

    #[test]
    fn split_keeps_shell_style_comments() {
        let line = "display processes # trailing comment";
        let parts = split_command_line(line).unwrap();
        assert_eq!(parts, vec!["display", "processes"]);
    }
}
