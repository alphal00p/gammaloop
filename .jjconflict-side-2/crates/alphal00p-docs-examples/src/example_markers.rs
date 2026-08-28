use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

use eyre::{ContextCompat, Result, ensure};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ExampleMarker {
    pub(crate) mode: String,
    pub(crate) id: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct FencedBlock {
    pub(crate) language: String,
    pub(crate) marker: Option<ExampleMarker>,
    pub(crate) line: usize,
    pub(crate) code: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct NamedBlock {
    pub(crate) source: PathBuf,
    pub(crate) mode: String,
    pub(crate) line: usize,
}

pub(crate) fn fenced_blocks(source: &str, path: &Path) -> Result<Vec<FencedBlock>> {
    let lines = source.lines().collect::<Vec<_>>();
    let mut blocks = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        let Some(language) = lines[index].strip_prefix("```") else {
            index += 1;
            continue;
        };
        let language = language.trim().to_owned();
        ensure!(
            !language.is_empty(),
            "untyped fenced block at {}:{}",
            path.display(),
            index + 1
        );
        let marker = index
            .checked_sub(1)
            .and_then(|line| lines[line].trim().strip_prefix("// docs-example:"))
            .map(|marker| parse_marker(marker.trim_start(), path, index))
            .transpose()?;
        let start = index + 1;
        index += 1;
        let mut code = Vec::new();
        while index < lines.len() && lines[index] != "```" {
            code.push(lines[index]);
            index += 1;
        }
        ensure!(
            index < lines.len(),
            "unterminated {language} block at {}:{}",
            path.display(),
            start
        );
        blocks.push(FencedBlock {
            language,
            marker,
            line: start,
            code: code.join("\n"),
        });
        index += 1;
    }
    Ok(blocks)
}

pub(crate) fn named_blocks(sources: &[(PathBuf, String)]) -> Result<BTreeMap<String, NamedBlock>> {
    let mut named = BTreeMap::new();
    for (source_path, source) in sources {
        for block in fenced_blocks(source, source_path)? {
            let Some(ExampleMarker { mode, id: Some(id) }) = block.marker else {
                continue;
            };
            let candidate = NamedBlock {
                source: source_path.clone(),
                mode,
                line: block.line,
            };
            if let Some(previous) = named.insert(id.clone(), candidate.clone()) {
                ensure!(
                    false,
                    "duplicate documentation example marker {id:?} at {}:{} and {}:{}",
                    previous.source.display(),
                    previous.line,
                    candidate.source.display(),
                    candidate.line
                );
            }
        }
    }
    Ok(named)
}

pub(crate) fn require_named_block(
    named: &BTreeMap<String, NamedBlock>,
    id: &str,
    source: &Path,
    mode: &str,
) -> Result<()> {
    let block = named
        .get(id)
        .wrap_err_with(|| format!("example {id} has no named documentation block"))?;
    ensure!(
        block.source == source,
        "example {id} marker is in {}, expected {}",
        block.source.display(),
        source.display()
    );
    ensure!(
        block.mode == mode,
        "example {id} marker uses verification mode {}, expected {mode}",
        block.mode
    );
    Ok(())
}

fn parse_marker(marker: &str, path: &Path, line: usize) -> Result<ExampleMarker> {
    let mut fields = marker.split_whitespace();
    let mode = fields.next().unwrap_or_default();
    ensure!(
        !mode.is_empty(),
        "empty documentation example mode at {}:{line}",
        path.display()
    );
    let id = fields.next().map(str::to_owned);
    ensure!(
        fields.next().is_none(),
        "documentation example marker at {}:{line} must be MODE or MODE ID",
        path.display()
    );
    Ok(ExampleMarker {
        mode: mode.to_owned(),
        id,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_legacy_and_named_markers() {
        let blocks = fenced_blocks(
            "// docs-example: run\n```rust\nfn main() {}\n```\n\
             // docs-example: compile canonical-id\n```python\npass\n```\n",
            Path::new("manual.typ"),
        )
        .unwrap();

        assert_eq!(
            blocks[0].marker,
            Some(ExampleMarker {
                mode: "run".to_owned(),
                id: None,
            })
        );
        assert_eq!(
            blocks[1].marker,
            Some(ExampleMarker {
                mode: "compile".to_owned(),
                id: Some("canonical-id".to_owned()),
            })
        );
    }

    #[test]
    fn rejects_ambiguous_or_duplicate_named_markers() {
        let empty = fenced_blocks(
            "// docs-example:\n```rust\nfn main() {}\n```\n",
            Path::new("manual.typ"),
        )
        .unwrap_err();
        assert!(
            empty
                .to_string()
                .contains("empty documentation example mode")
        );

        let malformed = fenced_blocks(
            "// docs-example: run one extra\n```rust\nfn main() {}\n```\n",
            Path::new("manual.typ"),
        )
        .unwrap_err();
        assert!(malformed.to_string().contains("must be MODE or MODE ID"));

        let sources = vec![
            (
                PathBuf::from("first.typ"),
                "// docs-example: run duplicate\n```rust\nfn main() {}\n```\n".to_owned(),
            ),
            (
                PathBuf::from("second.typ"),
                "// docs-example: compile duplicate\n```python\npass\n```\n".to_owned(),
            ),
        ];
        let duplicate = named_blocks(&sources).unwrap_err();
        assert!(
            duplicate
                .to_string()
                .contains("duplicate documentation example marker")
        );
    }

    #[test]
    fn requires_the_registered_source_and_mode() {
        let sources = vec![(
            PathBuf::from("actual.typ"),
            "// docs-example: compile canonical-id\n```python\npass\n```\n".to_owned(),
        )];
        let named = named_blocks(&sources).unwrap();

        require_named_block(&named, "canonical-id", Path::new("actual.typ"), "compile").unwrap();
        assert!(
            require_named_block(
                &named,
                "canonical-id",
                Path::new("elsewhere.typ"),
                "compile"
            )
            .unwrap_err()
            .to_string()
            .contains("expected elsewhere.typ")
        );
        assert!(
            require_named_block(&named, "canonical-id", Path::new("actual.typ"), "run")
                .unwrap_err()
                .to_string()
                .contains("expected run")
        );
    }
}
