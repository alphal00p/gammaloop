#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CommandLineParseError {
    IncompleteShellSyntax,
}

pub fn split_command_line(input: &str) -> Result<Vec<String>, CommandLineParseError> {
    let normalized = escape_process_ref_hash_ids(input);
    shlex::split(&normalized).ok_or(CommandLineParseError::IncompleteShellSyntax)
}

pub fn normalize_clap_args(args: Vec<String>) -> Vec<String> {
    let mut normalized = Vec::with_capacity(args.len());
    let mut index = 0usize;
    while index < args.len() {
        let current = &args[index];
        if matches!(current.as_str(), "--target" | "-t")
            && index + 2 < args.len()
            && !args[index + 1].contains('=')
            && is_shared_target_component(&args[index + 1])
            && is_shared_target_component(&args[index + 2])
        {
            normalized.push(format!("{current}={},{}", args[index + 1], args[index + 2]));
            index += 3;
            continue;
        }

        normalized.push(current.clone());
        index += 1;
    }
    normalized
}

pub fn split_command_list(input: &str) -> Result<Vec<String>, CommandLineParseError> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Mode {
        Unquoted,
        SingleQuoted,
        DoubleQuoted,
    }

    let mut parts = Vec::new();
    let mut current = String::new();
    let chars = input.chars().peekable();
    let mut mode = Mode::Unquoted;
    let mut unquoted_escape = false;
    let mut double_quoted_escape = false;

    for ch in chars {
        match mode {
            Mode::Unquoted => {
                if unquoted_escape {
                    current.push(ch);
                    unquoted_escape = false;
                    continue;
                }

                match ch {
                    '\\' => {
                        current.push(ch);
                        unquoted_escape = true;
                    }
                    '\'' => {
                        current.push(ch);
                        mode = Mode::SingleQuoted;
                    }
                    '"' => {
                        current.push(ch);
                        mode = Mode::DoubleQuoted;
                    }
                    ';' => {
                        let trimmed = current.trim();
                        if !trimmed.is_empty() {
                            parts.push(trimmed.to_string());
                        }
                        current.clear();
                    }
                    _ => current.push(ch),
                }
            }
            Mode::SingleQuoted => {
                current.push(ch);
                if ch == '\'' {
                    mode = Mode::Unquoted;
                }
            }
            Mode::DoubleQuoted => {
                current.push(ch);
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

    if unquoted_escape || double_quoted_escape || mode != Mode::Unquoted {
        return Err(CommandLineParseError::IncompleteShellSyntax);
    }

    let trimmed = current.trim();
    if !trimmed.is_empty() {
        parts.push(trimmed.to_string());
    }

    Ok(parts)
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

fn is_shared_target_component(token: &str) -> bool {
    token.parse::<f64>().is_ok()
}

#[cfg(test)]
mod test {
    use super::{normalize_clap_args, split_command_line, split_command_list};

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
        let line = "display integrand -p #12";
        let parts = split_command_line(line).unwrap();
        assert_eq!(parts, vec!["display", "integrand", "-p", "#12"]);
    }

    #[test]
    fn split_keeps_shell_style_comments() {
        let line = "display processes # trailing comment";
        let parts = split_command_line(line).unwrap();
        assert_eq!(parts, vec!["display", "processes"]);
    }

    #[test]
    fn split_preserves_generate_process_tokens() {
        let parts = split_command_line("generate g g > h").unwrap();
        assert_eq!(parts, vec!["generate", "g", "g", ">", "h"]);
    }

    #[test]
    fn split_preserves_generate_process_tokens_before_flags() {
        let parts = split_command_line("generate --only-diagrams g g > h").unwrap();
        assert_eq!(
            parts,
            vec!["generate", "--only-diagrams", "g", "g", ">", "h"]
        );
    }

    #[test]
    fn split_command_list_handles_top_level_semicolons() {
        let parts = split_command_list("display processes; display settings global ;display model")
            .unwrap();
        assert_eq!(
            parts,
            vec![
                "display processes",
                "display settings global",
                "display model"
            ]
        );
    }

    #[test]
    fn split_command_list_preserves_quoted_semicolons() {
        let parts = split_command_list(
            "set process -p foo string 'a; b'; run block_a -c \"display settings global; display model -a\"",
        )
        .unwrap();
        assert_eq!(
            parts,
            vec![
                "set process -p foo string 'a; b'",
                "run block_a -c \"display settings global; display model -a\"",
            ]
        );
    }

    #[test]
    fn split_command_list_rejects_incomplete_shell_syntax() {
        let err =
            split_command_list("display processes; set process string 'unterminated").unwrap_err();
        assert_eq!(err, super::CommandLineParseError::IncompleteShellSyntax);
    }

    #[test]
    fn normalize_clap_args_rewrites_shared_target_components() {
        let args = vec![
            "integrate".to_string(),
            "--target".to_string(),
            "-1.0214510394091818e-6".to_string(),
            "0.0".to_string(),
            "--restart".to_string(),
        ];

        assert_eq!(
            normalize_clap_args(args),
            vec![
                "integrate".to_string(),
                "--target=-1.0214510394091818e-6,0.0".to_string(),
                "--restart".to_string(),
            ]
        );
    }

    #[test]
    fn normalize_clap_args_preserves_repeated_keyed_targets() {
        let args = vec![
            "integrate".to_string(),
            "--target".to_string(),
            "aa_aa@1L=-1.0,0.0".to_string(),
            "--target".to_string(),
            "aa_aa@1L_up=-1.0,0.0".to_string(),
        ];

        assert_eq!(normalize_clap_args(args.clone()), args);
    }
}
