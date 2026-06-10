use proc_macro::TokenStream;
use quote::quote;

/// Shorthand for debug-level tracing spans that carry GammaLoop span context.
///
/// The arguments are appended to the generated `fields(...)` list:
///
/// ```ignore
/// #[debug_instrument(graph = %graph.log_display(), current = %current.log_display())]
/// ```
///
/// expands to:
///
/// ```ignore
/// #[tracing::instrument(
///     skip_all,
///     level = "debug",
///     fields(
///         span_context = true,
///         graph = %graph.log_display(),
///         current = %current.log_display(),
///     )
/// )]
/// ```
#[proc_macro_attribute]
pub fn debug_instrument(attr: TokenStream, item: TokenStream) -> TokenStream {
    let fields = proc_macro2::TokenStream::from(attr);
    let item = proc_macro2::TokenStream::from(item);

    let extra_fields = if fields.is_empty() {
        quote! {}
    } else {
        quote! {, #fields}
    };

    quote! {
        #[tracing::instrument(skip_all, level = "debug", fields(span_context = true #extra_fields))]
        #item
    }
    .into()
}
