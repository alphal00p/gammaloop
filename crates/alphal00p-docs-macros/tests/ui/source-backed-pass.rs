use alphal00p_docs::{func, macro_item, trait_item, ty};
use alphal00p_docs_schema::{DocFormat, DocItemKind, DocMemberKind};

const SOURCE: &str = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt";

#[func(
    id = "compute",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::compute"
)]
fn compute_marker() {}

#[func(
    id = "ExternalRecord::increment",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::ExternalRecord::increment"
)]
fn method_marker() {}

#[ty(
    id = "ExternalRecord",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::ExternalRecord"
)]
fn type_marker() {}

#[trait_item(
    id = "ExternalOperation",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::ExternalOperation"
)]
fn trait_marker() {}

#[macro_item(
    id = "external_record",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::external_record"
)]
fn macro_marker() {}

#[macro_item(
    id = "external_bang",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::external_bang"
)]
fn bang_marker() {}

#[macro_item(
    id = "external_attribute",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::external_attribute"
)]
fn attribute_marker() {}

#[macro_item(
    id = "ExternalDerive",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::ExternalDerive"
)]
fn derive_marker() {}

#[ty(
    id = "Undocumented",
    summary = "Fallback prose comes from the required adapter summary.",
    format = "rust-markdown",
    source = "crates/alphal00p-docs-macros/tests/ui/source-backed-source.txt",
    source_id = "fixture::Undocumented"
)]
fn undocumented_marker() {}

fn main() {
    let function = __alphal00p_docs_func_compute_marker();
    assert_eq!(function.id, "compute");
    assert_eq!(function.kind, DocItemKind::Function);
    assert_eq!(function.params[0].name, "value");
    assert!(function.signature.as_deref().unwrap().contains("value : u32"));
    assert!(function.returns.is_some());
    assert_eq!(function.docs.as_ref().unwrap().format, DocFormat::RustMarkdown);
    assert!(function.docs.unwrap().body.contains("second paragraph"));
    let source = function.source.unwrap();
    assert_eq!(source.identifier, "fixture::compute");
    assert_eq!(source.file, SOURCE);
    assert_eq!(source.line, 4);

    let method = __alphal00p_docs_func_method_marker();
    assert_eq!(method.kind, DocItemKind::Method);
    assert_eq!(method.params[0].name, "self");
    assert_eq!(method.params[1].name, "amount");

    let ty = __alphal00p_docs_ty_type_marker();
    assert_eq!(ty.members.len(), 1);
    assert_eq!(ty.members[0].kind, DocMemberKind::Field);
    assert_eq!(ty.members[0].name, "value");

    let documented_trait = __alphal00p_docs_trait_trait_marker();
    assert_eq!(documented_trait.members[0].kind, DocMemberKind::Method);
    let signature = documented_trait.members[0].signature.as_deref().unwrap();
    assert!(signature.ends_with("{ … }"));
    assert!(!signature.contains("{ value }"));

    let documented_macro = __alphal00p_docs_macro_macro_marker();
    assert_eq!(documented_macro.kind, DocItemKind::ExportedMacro);

    let function_like = __alphal00p_docs_macro_bang_marker();
    assert!(
        function_like
            .signature
            .as_deref()
            .unwrap()
            .starts_with("external_bang!(...)")
    );
    assert_eq!(function_like.source.unwrap().line, 46);

    let attribute = __alphal00p_docs_macro_attribute_marker();
    assert!(
        attribute
            .signature
            .as_deref()
            .unwrap()
            .starts_with("#[external_attribute]")
    );
    assert_eq!(attribute.source.unwrap().line, 52);

    let derive = __alphal00p_docs_macro_derive_marker();
    assert!(
        derive
            .signature
            .as_deref()
            .unwrap()
            .starts_with("#[derive(ExternalDerive)]")
    );
    assert_eq!(derive.source.unwrap().line, 58);

    let fallback = __alphal00p_docs_ty_undocumented_marker();
    assert_eq!(
        fallback.docs.unwrap().body,
        "Fallback prose comes from the required adapter summary."
    );
}
