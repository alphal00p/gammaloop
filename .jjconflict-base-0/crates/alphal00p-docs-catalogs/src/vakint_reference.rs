//! Vakint topology and external-tool metadata exported from its runtime generator.

use std::collections::BTreeSet;

use eyre::{ContextCompat, Result, ensure};
use regex::Regex;
use vakint::Vakint;

use crate::generated::{
    DependencyReference, GENERATED_REFERENCE_SCHEMA, GeneratedSource, TopologyReference,
    VakintReference,
};

const VAKINT_SOURCE: &str = include_str!("../../vakint/src/lib.rs");

pub fn export() -> Result<VakintReference> {
    let rendered = Vakint::new()?.topologies.to_string();
    let ansi = Regex::new(r"\x1b\[[0-9;]*m")?;
    let rendered = ansi.replace_all(&rendered, "");
    let topology = Regex::new(
        r"name='(?P<name>[^']+)', n_loops=(?P<loops>[0-9]+), n_props_top_topo=(?P<props>[0-9]+)",
    )?;
    let mut topologies = Vec::new();
    for captures in topology.captures_iter(&rendered) {
        let name = captures["name"].to_owned();
        if name.eq_ignore_ascii_case("unknown") {
            continue;
        }
        topologies.push(TopologyReference {
            name,
            loops: captures["loops"].parse()?,
            propagator_slots: captures["props"].parse()?,
        });
    }
    ensure!(
        !topologies.is_empty(),
        "Vakint's runtime topology generator exported no supported topologies"
    );
    let mut names = BTreeSet::new();
    for topology in &topologies {
        ensure!(
            names.insert(topology.name.as_str()),
            "Vakint generated duplicate topology name {}",
            topology.name
        );
    }

    let dependencies = vec![
        dependency("FORM", "MINIMAL_FORM_VERSION")?,
        dependency("pySecDec", "MINIMAL_PYSECDEC_VERSION")?,
    ];
    Ok(VakintReference {
        schema_version: GENERATED_REFERENCE_SCHEMA,
        sources: vec![
            GeneratedSource {
                kind: "runtime-topology-generator".to_owned(),
                path: "crates/vakint/src/topologies.rs".to_owned(),
                symbol: "vakint::topologies::Topologies::generate_topologies".to_owned(),
            },
            GeneratedSource {
                kind: "source-version-constant".to_owned(),
                path: "crates/vakint/src/lib.rs".to_owned(),
                symbol: "MINIMAL_FORM_VERSION, MINIMAL_PYSECDEC_VERSION".to_owned(),
            },
        ],
        dependencies,
        topologies,
    })
}

fn dependency(name: &str, symbol: &str) -> Result<DependencyReference> {
    let expression = format!(r#"static {symbol}: &str = "(?P<version>[^"]+)";"#);
    let pattern = Regex::new(&expression)?;
    let version = pattern
        .captures(VAKINT_SOURCE)
        .and_then(|captures| captures.name("version"))
        .map(|capture| capture.as_str().to_owned())
        .with_context(|| format!("could not find Vakint dependency constant {symbol}"))?;
    Ok(DependencyReference {
        name: name.to_owned(),
        minimum_version: version,
        source_symbol: symbol.to_owned(),
    })
}

#[cfg(test)]
mod tests {
    use super::export;

    #[test]
    fn runtime_topologies_and_versions_are_exported() {
        let reference = export().unwrap();
        assert!(
            reference
                .topologies
                .iter()
                .any(|topology| topology.name == "I1L")
        );
        assert!(
            reference
                .topologies
                .iter()
                .any(|topology| topology.loops == 4)
        );
        assert_eq!(reference.dependencies[0].minimum_version, "4.2.1");
        assert_eq!(reference.dependencies[1].minimum_version, "1.6.4");
    }
}
