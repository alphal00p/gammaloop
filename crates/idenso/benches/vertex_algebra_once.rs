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
        eprintln!("available: network_vertex_substitution_8, network_full_algebra_8");
        return;
    };

    if selected(&filter, "network_vertex_substitution_8") {
        let start = Instant::now();
        let result = common::network_vertex_substitution(common::network_vertex_fixture_8());
        report("network_vertex_substitution_8", start, &result);
    }

    if selected(&filter, "network_full_algebra_8") {
        let start = Instant::now();
        let result = common::network_full_algebra_8(common::network_vertex_fixture_8());
        common::assert_no_network_internal_indices(&result);
        report("network_full_algebra_8", start, &result);
    }
}
