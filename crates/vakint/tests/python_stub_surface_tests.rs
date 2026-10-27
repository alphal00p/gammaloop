#![cfg(all(feature = "python_stubgen", feature = "symbolica_community_module"))]

use pyo3::prelude::*;

#[test]
fn graph_adapter_and_rich_output_are_present_in_valid_python_stubs() {
    Python::initialize();
    let info = vakint::stub_info().expect("build Python stub inventory");
    let module = &info.modules["symbolica.community.hepkit.vakint"];
    let graph_adapter = &module.function["integral_from_diagram"];
    assert_eq!(graph_adapter.len(), 1);
    assert_eq!(graph_adapter[0].r#return.name, "VakintExpression");
    assert_eq!(
        graph_adapter[0]
            .parameters
            .positional_or_keyword
            .iter()
            .map(|parameter| parameter.name)
            .collect::<Vec<_>>(),
        ["diagram", "family", "numerator"]
    );
    assert_eq!(
        graph_adapter[0]
            .parameters
            .keyword_only
            .iter()
            .map(|parameter| parameter.name)
            .collect::<Vec<_>>(),
        [
            "powers",
            "physical_slots",
            "parameter_substitutions",
            "external_momenta"
        ]
    );

    for name in ["VakintExpression", "VakintNumericalResult"] {
        let class = module
            .class
            .values()
            .find(|class| class.name == name)
            .unwrap_or_else(|| panic!("missing {name} stub"));
        let methods = &class.methods["formatted"];
        assert_eq!(methods.len(), 1);
        assert_eq!(methods[0].r#return.name, "symbolica.FormattedOutput");
        assert!(methods[0].r#return.import.contains(&"symbolica".into()));
    }

    // Parse the emitted source too: Rust defaults on overridden parameter types
    // must still produce valid Python syntax in the stub.
    Python::attach(|py| {
        py.import("ast")
            .expect("import Python's ast module")
            .call_method1("parse", (module.to_string(),))
            .expect("the Vakint stub must be valid Python");
    });
}
