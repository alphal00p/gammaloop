use alphal00p_docs::{func, macro_item, scope, trait_item, ty};
use alphal00p_docs_schema::{DocFormat, DocItemKind, DocMemberKind, DocScope};

#[func(id = "add", keyword = "arithmetic")]
/// Adds two values. #emph[The full body remains Typst markup.]
///
/// ```rust
/// assert_eq!(add(2, 3), 5);
/// ```
pub fn add(left: u32, right: u32) -> u32 {
    left + right
}

#[ty]
/// A pair of values. The descriptor includes its fields.
pub struct Pair {
    /// Left value. #emph[Stored as Typst markup.]
    pub left: u32,
    /// Right value.
    pub right: u32,
    private: u32,
    pub(crate) crate_visible: u32,
}

#[ty]
/// Selection mode. Variants remain ordered.
pub enum Mode {
    /// Selects the #strong[fast] path.
    Fast,
    /// Retries failed work.
    Careful {
        /// Maximum #emph[number] of attempts.
        attempts: usize,
    },
}

#[trait_item]
/// A named object. Implementations supply the name.
pub trait Named {
    /// Prefixes every #emph[name].
    const PREFIX: &'static str;
    /// Reports naming #emph[failures].
    type Error;
    /// Returns this object's #emph[name].
    fn name(&self) -> Result<String, Self::Error>;
    /// Returns the default #emph[name].
    fn default_name() -> String {
        "fallback".to_owned()
    }
}

pub struct Counter(u32);

impl Counter {
    #[func]
    /// Multiplies the counter. #emph[Method details remain Typst markup.]
    pub fn scaled(&self, factor: u32) -> u32 {
        self.0 * factor
    }

    #[func(id = "origin", owner = "Counter")]
    /// Creates a zero counter. #emph[Associated function details remain Typst markup.]
    pub fn zero() -> Self {
        Self(0)
    }
}

pub struct Gauge;

impl Gauge {
    #[func]
    /// Scales a gauge. The source identifier must not collide with `Counter::scaled`.
    pub fn scaled(&self, factor: u32) -> u32 {
        factor
    }
}

pub struct Wrapper<T>(T);

impl<T> Wrapper<T> {
    #[func]
    /// Returns the wrapped value. Generic arguments do not belong in Rustdoc paths.
    pub fn into_inner(self) -> T {
        self.0
    }

    #[func(owner = "Wrapper<T>")]
    /// Creates a wrapper. The declared generic owner is compile-time checked.
    pub fn new(value: T) -> Self {
        Self(value)
    }
}

#[macro_item]
/// Doubles one expression. The invocation contract is cataloged.
macro_rules! twice {
    ($value:expr) => {
        $value * 2
    };
}

#[scope(title = "Public API")]
/// Ordered public surface. Adapters insert descriptors explicitly.
pub mod api {}

fn main() {
    let function = __alphal00p_docs_func_add();
    assert_eq!(function.id, "add");
    assert_eq!(function.params[0].name, "left");
    assert_eq!(function.examples[0].language, "rust");
    assert!(function.returns.is_some());
    assert!(function.source.unwrap().identifier.ends_with("::add"));
    assert!(function.docs.unwrap().body.contains("full body"));

    let method = Counter::__alphal00p_docs_func_scaled();
    assert_eq!(method.id, "Counter::scaled");
    assert_eq!(method.kind, DocItemKind::Method);
    assert_eq!(method.params[0].name, "self");
    assert_eq!(method.docs.unwrap().format, DocFormat::TypstMarkup);
    let method_source = method.source.unwrap().identifier;
    assert!(
        method_source.ends_with("::Counter::scaled"),
        "{method_source}"
    );

    let associated = Counter::__alphal00p_docs_func_zero();
    assert_eq!(associated.id, "Counter::origin");
    assert_eq!(associated.kind, DocItemKind::Method);
    let associated_source = associated.source.unwrap().identifier;
    assert!(
        associated_source.ends_with("::Counter::zero"),
        "{associated_source}"
    );

    let gauge = Gauge::__alphal00p_docs_func_scaled();
    assert_eq!(gauge.id, "Gauge::scaled");
    let gauge_source = gauge.source.unwrap().identifier;
    assert!(gauge_source.ends_with("::Gauge::scaled"), "{gauge_source}");
    assert_ne!(method_source, gauge_source);

    let generic = Wrapper::<u8>::__alphal00p_docs_func_into_inner();
    assert_eq!(generic.id, "Wrapper::into_inner");
    let generic_source = generic.source.unwrap().identifier;
    assert!(
        generic_source.ends_with("::Wrapper::into_inner"),
        "{generic_source}"
    );
    let generic_associated = Wrapper::<u8>::__alphal00p_docs_func_new();
    assert_eq!(generic_associated.id, "Wrapper::new");
    let generic_associated_source = generic_associated.source.unwrap().identifier;
    assert!(
        generic_associated_source.ends_with("::Wrapper::new"),
        "{generic_associated_source}"
    );

    let mut methods = DocScope::new("methods", "Methods");
    methods
        .define_item(Counter::__alphal00p_docs_func_scaled())
        .unwrap();
    methods
        .define_item(Gauge::__alphal00p_docs_func_scaled())
        .unwrap();
    methods
        .define_item(Counter::__alphal00p_docs_func_zero())
        .unwrap();
    assert_eq!(
        methods.items.keys().map(String::as_str).collect::<Vec<_>>(),
        vec!["Counter::scaled", "Gauge::scaled", "Counter::origin"]
    );

    let pair = __alphal00p_docs_ty_Pair();
    assert_eq!(
        pair.members
            .iter()
            .map(|member| member.name.as_str())
            .collect::<Vec<_>>(),
        vec!["left", "right"]
    );
    assert_eq!(pair.members[0].kind, DocMemberKind::Field);
    assert_eq!(
        pair.members[0].docs.as_ref().unwrap().format,
        DocFormat::TypstMarkup
    );
    assert!(
        pair.members[0]
            .docs
            .as_ref()
            .unwrap()
            .body
            .contains("#emph")
    );

    let mode = __alphal00p_docs_ty_Mode();
    assert_eq!(mode.members[0].kind, DocMemberKind::Variant);
    assert_eq!(
        mode.members[0].docs.as_ref().unwrap().format,
        DocFormat::TypstMarkup
    );
    assert_eq!(
        mode.members[1].members[0].docs.as_ref().unwrap().format,
        DocFormat::TypstMarkup
    );

    let named = __alphal00p_docs_trait_Named();
    assert_eq!(named.members[0].kind, DocMemberKind::AssociatedConst);
    assert_eq!(named.members[1].kind, DocMemberKind::AssociatedType);
    assert_eq!(named.members[2].kind, DocMemberKind::Method);
    assert_eq!(named.members[3].kind, DocMemberKind::AssociatedFunction);
    assert!(
        named.members[2]
            .signature
            .as_deref()
            .unwrap()
            .ends_with(';')
    );
    let default_signature = named.members[3].signature.as_deref().unwrap();
    assert!(default_signature.ends_with("{ … }"));
    assert!(!default_signature.contains("fallback"));
    assert!(
        named
            .members
            .iter()
            .all(|member| { member.docs.as_ref().unwrap().format == DocFormat::TypstMarkup })
    );

    let documented_macro = __alphal00p_docs_macro_twice();
    let mut scope = __alphal00p_docs_scope_api();
    scope.define_item(pair).unwrap();
    scope.define_item(mode).unwrap();
    scope.define_item(named).unwrap();
    scope.define_item(documented_macro).unwrap();
    scope.validate().unwrap();
    assert_eq!(
        scope.items.keys().map(String::as_str).collect::<Vec<_>>(),
        vec!["Pair", "Mode", "Named", "twice"]
    );

    assert_eq!(add(2, 3), 5);
    assert_eq!(twice!(4), 8);
}
