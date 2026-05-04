use std::time::Instant;

mod common;

fn selected(filter: &str, name: &str) -> bool {
    name.contains(filter)
}

fn report(name: &str, start: Instant, result: &symbolica::atom::Atom) {
    let elapsed = start.elapsed();
    println!(
        "{name}: elapsed={elapsed:.3?} terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
}

fn main() {
    let Some(filter) = std::env::args().nth(1) else {
        eprintln!(
            "usage: cargo bench --package idenso --bench vertex_algebra_once --profile dev-optim -- <filter>"
        );
        eprintln!(
            "available: bare_vertex_substitution_5, bare_vertex_substitution_8, network_vertex_substitution_8, network_full_algebra_N where N is 5..8"
        );
        return;
    };

    if selected(&filter, "bare_vertex_substitution_5") {
        let start = Instant::now();
        let result = common::bare_vertex_substitution(common::bare_vertex_fixture_5());
        report("bare_vertex_substitution_5", start, &result);
    }

    if selected(&filter, "bare_vertex_substitution_8") {
        let start = Instant::now();
        let result = common::bare_vertex_substitution(common::bare_vertex_fixture_8());
        report("bare_vertex_substitution_8", start, &result);
    }

    if selected(&filter, "network_vertex_substitution_8") {
        let start = Instant::now();
        let result = common::network_vertex_substitution(common::network_vertex_fixture_8());
        report("network_vertex_substitution_8", start, &result);
    }

    for vertex_count in 5..=8 {
        let name = format!("network_full_algebra_{vertex_count}");
        if selected(&filter, &name) {
            let start = Instant::now();
            let result =
                common::network_schoonschip_substituted(common::network_vertex_substitution(
                    common::network_vertex_fixture_with_count(vertex_count),
                ));
            report(&name, start, &result);
        }
    }
}
