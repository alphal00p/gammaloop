use proc_macro::TokenStream;
use proc_macro2::TokenStream as TokenStream2;
use quote::{format_ident, quote};
use syn::{Attribute, Expr, ExprLit, Lit, Meta, Token, parse::Parser, punctuated::Punctuated};

#[derive(Default)]
pub(crate) struct CommonArgs {
    pub id: Option<String>,
    pub name: Option<String>,
    pub title: Option<String>,
    pub summary: Option<String>,
    pub since: Option<String>,
    pub format: Option<String>,
    pub owner: Option<String>,
    /// Workspace-relative source file read by an internal adapter descriptor.
    pub source: Option<String>,
    /// Stable public identifier selected from `source`.
    pub source_id: Option<String>,
    /// Optional item-kind override for source-backed adapters.
    pub kind: Option<String>,
    pub keywords: Vec<String>,
}

impl CommonArgs {
    pub fn parse(tokens: TokenStream) -> syn::Result<Self> {
        let metas = Punctuated::<Meta, Token![,]>::parse_terminated.parse(tokens)?;
        let mut args = Self::default();
        for meta in metas {
            let Meta::NameValue(value) = meta else {
                return Err(syn::Error::new_spanned(meta, "expected name = string"));
            };
            let Some(key) = value.path.get_ident().map(ToString::to_string) else {
                return Err(syn::Error::new_spanned(
                    value.path,
                    "expected a simple name",
                ));
            };
            let Expr::Lit(ExprLit {
                lit: Lit::Str(value),
                ..
            }) = value.value
            else {
                return Err(syn::Error::new_spanned(
                    value.value,
                    "expected a string literal",
                ));
            };
            let value = value.value();
            match key.as_str() {
                "id" => set_once(&mut args.id, value, &key)?,
                "name" => set_once(&mut args.name, value, &key)?,
                "title" => set_once(&mut args.title, value, &key)?,
                "summary" => set_once(&mut args.summary, value, &key)?,
                "since" => set_once(&mut args.since, value, &key)?,
                "format" => set_once(&mut args.format, value, &key)?,
                "owner" => set_once(&mut args.owner, value, &key)?,
                "source" => set_once(&mut args.source, value, &key)?,
                "source_id" => set_once(&mut args.source_id, value, &key)?,
                "kind" => set_once(&mut args.kind, value, &key)?,
                "keyword" => args.keywords.push(value),
                _ => {
                    return Err(syn::Error::new_spanned(
                        value,
                        format!("unknown documentation argument {key}"),
                    ));
                }
            }
        }
        Ok(args)
    }

    pub fn format_tokens(&self) -> syn::Result<TokenStream2> {
        match self.format.as_deref().unwrap_or("typst") {
            "typst" | "typst-markup" => Ok(quote!(::alphal00p_docs_schema::DocFormat::TypstMarkup)),
            "rust-markdown" | "markdown" => {
                Ok(quote!(::alphal00p_docs_schema::DocFormat::RustMarkdown))
            }
            "python-docstring" => Ok(quote!(::alphal00p_docs_schema::DocFormat::PythonDocstring)),
            "plain-text" => Ok(quote!(::alphal00p_docs_schema::DocFormat::PlainText)),
            value => Err(syn::Error::new(
                proc_macro2::Span::call_site(),
                format!("unknown documentation format {value}"),
            )),
        }
    }

    pub fn is_typst_markup(&self) -> bool {
        matches!(
            self.format.as_deref().unwrap_or("typst"),
            "typst" | "typst-markup"
        )
    }

    pub fn owner_type(&self) -> syn::Result<Option<syn::Type>> {
        self.owner
            .as_deref()
            .map(|owner| {
                syn::parse_str(owner).map_err(|error| {
                    syn::Error::new(
                        proc_macro2::Span::call_site(),
                        format!("invalid owner type `{owner}`: {error}"),
                    )
                })
            })
            .transpose()
    }

    pub fn source_item(&self) -> syn::Result<Option<(&str, &str)>> {
        match (self.source.as_deref(), self.source_id.as_deref()) {
            (None, None) => Ok(None),
            (Some(source), Some(identifier)) => Ok(Some((source, identifier))),
            _ => Err(syn::Error::new(
                proc_macro2::Span::call_site(),
                "`source` and `source_id` must be supplied together",
            )),
        }
    }

    pub fn reject_owner(&self, attribute: &str) -> syn::Result<()> {
        if self.owner.is_some() {
            Err(syn::Error::new(
                proc_macro2::Span::call_site(),
                format!("`owner` is only supported by #[alphal00p_docs::func], not {attribute}"),
            ))
        } else {
            Ok(())
        }
    }

    pub fn reject_source_item(&self, attribute: &str) -> syn::Result<()> {
        if self.source.is_some() || self.source_id.is_some() || self.kind.is_some() {
            Err(syn::Error::new(
                proc_macro2::Span::call_site(),
                format!("source-backed arguments are not supported by {attribute}"),
            ))
        } else {
            Ok(())
        }
    }
}

fn set_once(slot: &mut Option<String>, value: String, key: &str) -> syn::Result<()> {
    if slot.replace(value).is_some() {
        Err(syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("duplicate documentation argument {key}"),
        ))
    } else {
        Ok(())
    }
}

pub(crate) fn docs(attrs: &[Attribute]) -> String {
    attrs
        .iter()
        .filter_map(|attribute| {
            if !attribute.path().is_ident("doc") {
                return None;
            }
            let Meta::NameValue(value) = &attribute.meta else {
                return None;
            };
            let Expr::Lit(ExprLit {
                lit: Lit::Str(value),
                ..
            }) = &value.value
            else {
                return None;
            };
            Some(value.value().trim_start().to_owned())
        })
        .collect::<Vec<_>>()
        .join("\n")
        .trim()
        .to_owned()
}

pub(crate) fn take_docs(attrs: &mut Vec<Attribute>) -> String {
    let value = docs(attrs);
    attrs.retain(|attribute| !attribute.path().is_ident("doc"));
    value
}

pub(crate) fn replace_docs(attrs: &mut Vec<Attribute>, summary: &str) {
    if !summary.is_empty() {
        let summary = syn::LitStr::new(summary, proc_macro2::Span::call_site());
        attrs.push(syn::parse_quote!(#[doc = #summary]));
    }
}

pub(crate) fn first_sentence(docs: &str) -> String {
    let prose = docs
        .lines()
        .map(str::trim)
        .find(|line| {
            !line.is_empty() && !line.starts_with(['#', '=', '-', '*']) && !line.starts_with("~~~")
        })
        .unwrap_or_default();
    for (index, character) in prose.char_indices() {
        if character == '.'
            && prose[index + 1..]
                .chars()
                .next()
                .is_none_or(char::is_whitespace)
        {
            return prose[..=index].to_owned();
        }
    }
    prose.to_owned()
}

pub(crate) fn fenced_examples(docs: &str) -> Vec<(String, String)> {
    let mut examples = vec![];
    let mut lines = docs.lines();
    while let Some(line) = lines.next() {
        let trimmed = line.trim();
        let Some(language) = trimmed.strip_prefix("```") else {
            continue;
        };
        let mut code = vec![];
        for line in lines.by_ref() {
            if line.trim() == "```" {
                break;
            }
            code.push(line);
        }
        let code = code.join("\n").trim().to_owned();
        if !code.is_empty() {
            examples.push((language.trim().to_owned(), code));
        }
    }
    examples
}

pub(crate) fn docs_assignment(target: &str, docs: &str, format: &TokenStream2) -> TokenStream2 {
    if docs.is_empty() {
        quote!()
    } else {
        let target = format_ident!("{target}");
        quote! {
            #target.docs = Some(::alphal00p_docs_schema::DocText::new(#format, #docs));
        }
    }
}

pub(crate) fn string_assignment(target: &str, field: &str, value: Option<&str>) -> TokenStream2 {
    let Some(value) = value else {
        return quote!();
    };
    let target = format_ident!("{target}");
    let field = format_ident!("{field}");
    quote! {
        #target.#field = Some(#value.to_owned());
    }
}

#[cfg(test)]
mod tests {
    use super::{fenced_examples, first_sentence};

    #[test]
    fn extracts_first_prose_sentence() {
        assert_eq!(
            first_sentence("A complete sentence. More #emph[markup]."),
            "A complete sentence."
        );
        assert_eq!(
            first_sentence("= Heading\n\nThe useful sentence. Another."),
            "The useful sentence."
        );
    }

    #[test]
    fn extracts_fenced_examples_in_order() {
        assert_eq!(
            fenced_examples("Example.\n\n```rust\nlet answer = 42;\n```"),
            vec![("rust".to_owned(), "let answer = 42;".to_owned())]
        );
    }
}
