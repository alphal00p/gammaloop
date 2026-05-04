use std::time::Instant;

use idenso::schoonschip::SchoonschipContractionOrder;

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

fn order_from_name(name: &str) -> Option<SchoonschipContractionOrder> {
    match name {
        "smallest_degree" => Some(SchoonschipContractionOrder::SmallestDegree),
        "largest_degree" => Some(SchoonschipContractionOrder::LargestDegree),
        "min_largest_operand_bytes" => Some(SchoonschipContractionOrder::MinLargestOperandBytes),
        "min_product_terms" => Some(SchoonschipContractionOrder::MinProductTerms),
        "min_product_bytes" => Some(SchoonschipContractionOrder::MinProductBytes),
        "smallest_degree_min_largest_operand_bytes" => {
            Some(SchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes)
        }
        "smallest_degree_min_product_terms" => {
            Some(SchoonschipContractionOrder::SmallestDegreeMinProductTerms)
        }
        "smallest_degree_min_product_bytes" => {
            Some(SchoonschipContractionOrder::SmallestDegreeMinProductBytes)
        }
        _ => None,
    }
}

fn order_name(order: SchoonschipContractionOrder) -> &'static str {
    match order {
        SchoonschipContractionOrder::SmallestDegree => "smallest_degree",
        SchoonschipContractionOrder::LargestDegree => "largest_degree",
        SchoonschipContractionOrder::MinLargestOperandBytes => "min_largest_operand_bytes",
        SchoonschipContractionOrder::MinProductTerms => "min_product_terms",
        SchoonschipContractionOrder::MinProductBytes => "min_product_bytes",
        SchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes => {
            "smallest_degree_min_largest_operand_bytes"
        }
        SchoonschipContractionOrder::SmallestDegreeMinProductTerms => {
            "smallest_degree_min_product_terms"
        }
        SchoonschipContractionOrder::SmallestDegreeMinProductBytes => {
            "smallest_degree_min_product_bytes"
        }
    }
}

fn main() {
    let Some(filter) = std::env::args().nth(1) else {
        eprintln!(
            "usage: cargo bench --package idenso --bench vertex_algebra_once --profile dev-optim -- <filter>"
        );
        eprintln!(
            "available: bare_vertex_substitution_5, bare_vertex_substitution_8, network_vertex_substitution_8, network_full_algebra_N where N is 5..8"
        );
        eprintln!(
            "optional second arg for network_full_algebra_N: smallest_degree, largest_degree, min_largest_operand_bytes, min_product_terms, min_product_bytes, smallest_degree_min_largest_operand_bytes, smallest_degree_min_product_terms, smallest_degree_min_product_bytes"
        );
        return;
    };
    let order = std::env::args()
        .nth(2)
        .as_deref()
        .and_then(order_from_name)
        .unwrap_or_default();
    let assert_full_simplification = std::env::args().nth(3).as_deref() != Some("no_assert");

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
            let result = common::network_schoonschip_substituted_with_order(
                common::network_vertex_substitution(common::network_vertex_fixture_with_count(
                    vertex_count,
                )),
                order,
            );
            if vertex_count == 8 {
                let residual = common::residual_network_internal_indices(&result, vertex_count);
                if residual.is_empty() {
                    println!("network_full_algebra_8 residual_internal_indices=[]");
                } else {
                    println!("network_full_algebra_8 residual_internal_indices={residual:?}");
                    if assert_full_simplification {
                        panic!("internal indices remained: {residual:?}");
                    }
                }
            }
            report(&format!("{name}_{}", order_name(order)), start, &result);
        }
    }
}
