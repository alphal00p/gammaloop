//! Opt-in metadata attributes for the alphal00p documentation catalogs.

mod args;

use std::{env, fs, path::PathBuf};

use args::{
    CommonArgs, docs, docs_assignment, fenced_examples, first_sentence, replace_docs,
    string_assignment, take_docs,
};
use proc_macro::TokenStream;
use proc_macro2::TokenStream as TokenStream2;
use quote::{ToTokens, format_ident, quote};
use syn::{
    Attribute, Field, Fields, FnArg, ImplItem, ImplItemFn, Item, ItemEnum, ItemFn, ItemMacro,
    ItemMod, ItemStruct, ItemTrait, ItemType, Pat, ReturnType, TraitItem, Type, parse_macro_input,
};

#[derive(Clone, Copy)]
enum SourceContext<'a> {
    Module,
    Inherent {
        declared_owner: Option<&'a syn::Type>,
    },
    External {
        identifier: &'a str,
        file: &'a str,
        line: u32,
        column: u32,
    },
}

#[derive(Clone, Copy)]
enum ExternalKind {
    Function,
    Type,
    Trait,
    Macro,
}

/// Expand an adapter marker from the definition in another workspace source
/// file. The marker itself is consumed: only the static descriptor and an
/// `include_bytes!` dependency edge are emitted. This lets independently
/// publishable product crates remain free of internal documentation
/// dependencies while making their real comments and signatures authoritative.
fn expand_external_item(
    arguments: &CommonArgs,
    tokens: TokenStream2,
    expected: ExternalKind,
) -> syn::Result<Option<TokenStream2>> {
    let Some((source_file, source_id)) = arguments.source_item()? else {
        if arguments.kind.is_some() {
            return Err(syn::Error::new(
                proc_macro2::Span::call_site(),
                "`kind` requires `source` and `source_id`",
            ));
        }
        return Ok(None);
    };
    if arguments.owner.is_some() {
        return Err(syn::Error::new(
            proc_macro2::Span::call_site(),
            "`owner` cannot be combined with a source-backed adapter",
        ));
    }
    let marker = syn::parse2::<ItemFn>(tokens).map_err(|_| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            "a source-backed documentation attribute requires an empty marker function",
        )
    })?;
    if !marker.sig.inputs.is_empty()
        || !marker.sig.generics.params.is_empty()
        || marker.sig.generics.where_clause.is_some()
        || !matches!(marker.sig.output, ReturnType::Default)
        || marker.sig.constness.is_some()
        || marker.sig.asyncness.is_some()
        || marker.sig.unsafety.is_some()
        || marker.sig.abi.is_some()
        || marker.sig.variadic.is_some()
        || !marker.block.stmts.is_empty()
    {
        return Err(syn::Error::new_spanned(
            marker.sig,
            "a source-backed adapter must be an empty `fn marker() {}`",
        ));
    }
    let marker_ident = marker.sig.ident;
    let absolute = resolve_external_source(source_file)?;
    let source = fs::read_to_string(&absolute).map_err(|error| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!(
                "failed to read documentation source {}: {error}",
                absolute.display()
            ),
        )
    })?;
    let parsed = syn::parse_file(&source).map_err(|error| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("failed to parse documentation source {source_file}: {error}"),
        )
    })?;
    let source_name = source_id
        .rsplit("::")
        .next()
        .unwrap_or(source_id)
        .trim_end_matches('!');
    let dependency = format_ident!("__ALPHAL00P_DOCS_SOURCE_{}", marker_ident);
    let absolute = absolute.to_string_lossy().into_owned();

    let descriptor = match expected {
        ExternalKind::Function => external_function_descriptor(
            arguments,
            &parsed.items,
            &source,
            source_file,
            source_id,
            source_name,
            &marker_ident,
        )?,
        ExternalKind::Type => external_type_descriptor(
            arguments,
            &parsed.items,
            &source,
            source_file,
            source_id,
            source_name,
            &marker_ident,
        )?,
        ExternalKind::Trait => external_trait_descriptor(
            arguments,
            &parsed.items,
            &source,
            source_file,
            source_id,
            source_name,
            &marker_ident,
        )?,
        ExternalKind::Macro => external_macro_descriptor(
            arguments,
            &parsed.items,
            &source,
            source_file,
            source_id,
            source_name,
            &marker_ident,
        )?,
    };
    Ok(Some(quote! {
        #[doc(hidden)]
        #[allow(non_upper_case_globals)]
        const #dependency: usize = include_bytes!(#absolute).len();
        #descriptor
    }))
}

fn resolve_external_source(source: &str) -> syn::Result<PathBuf> {
    let relative = std::path::Path::new(source);
    if relative.is_absolute()
        || relative
            .components()
            .any(|component| component == std::path::Component::ParentDir)
    {
        return Err(syn::Error::new(
            proc_macro2::Span::call_site(),
            "source-backed documentation paths must be workspace-relative and cannot contain `..`",
        ));
    }
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").map_err(|error| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("CARGO_MANIFEST_DIR is unavailable: {error}"),
        )
    })?);
    for root in manifest.ancestors() {
        let candidate = root.join(relative);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }
    Err(syn::Error::new(
        proc_macro2::Span::call_site(),
        format!("workspace documentation source {source} does not exist"),
    ))
}

fn external_function_descriptor(
    arguments: &CommonArgs,
    items: &[Item],
    source: &str,
    source_file: &str,
    source_id: &str,
    name: &str,
    marker: &syn::Ident,
) -> syn::Result<TokenStream2> {
    if let Some(function) = items.iter().find_map(|item| match item {
        Item::Fn(function) if function.sig.ident == name => Some(function.clone()),
        _ => None,
    }) {
        let mut attributes = function.attrs;
        let docs = prepare_external_docs(arguments, &mut attributes);
        let signature = function.sig;
        let descriptor = format_ident!("__alphal00p_docs_func_{}", marker);
        let parameters = signature.inputs.iter().map(parameter);
        let returns = return_descriptor(&signature.output);
        let line = external_source_line(source, name, None, ExternalKind::Function)?;
        let kind = external_kind_tokens(arguments, "function")?;
        return Ok(item_descriptor(
            arguments,
            name.to_owned(),
            descriptor,
            kind,
            &docs,
            Some(signature.to_token_stream().to_string()),
            quote! {
                descriptor.params = vec![#(#parameters),*];
                #returns
            },
            SourceContext::External {
                identifier: source_id,
                file: source_file,
                line,
                column: 1,
            },
        ));
    }

    let owner = source_id.rsplit("::").nth(1).ok_or_else(|| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("source function {source_id} was not found"),
        )
    })?;
    let method = items.iter().find_map(|item| {
        let Item::Impl(item) = item else { return None };
        if external_type_name(&item.self_ty).as_deref() != Some(owner) {
            return None;
        }
        item.items.iter().find_map(|item| match item {
            ImplItem::Fn(function) if function.sig.ident == name => Some(function.clone()),
            _ => None,
        })
    });
    let mut method = method.ok_or_else(|| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("source method {source_id} was not found"),
        )
    })?;
    let docs = prepare_external_docs(arguments, &mut method.attrs);
    let signature = method.sig;
    let descriptor = format_ident!("__alphal00p_docs_func_{}", marker);
    let parameters = signature.inputs.iter().map(parameter);
    let returns = return_descriptor(&signature.output);
    let line = external_source_line(source, name, Some(owner), ExternalKind::Function)?;
    let kind = external_kind_tokens(arguments, "method")?;
    Ok(item_descriptor(
        arguments,
        name.to_owned(),
        descriptor,
        kind,
        &docs,
        Some(signature.to_token_stream().to_string()),
        quote! {
            descriptor.params = vec![#(#parameters),*];
            #returns
        },
        SourceContext::External {
            identifier: source_id,
            file: source_file,
            line,
            column: 1,
        },
    ))
}

fn return_descriptor(output: &ReturnType) -> TokenStream2 {
    match output {
        ReturnType::Default => quote!(),
        ReturnType::Type(_, ty) => {
            let signature = ty.to_token_stream().to_string();
            quote! {
                descriptor.returns = Some(::alphal00p_docs_schema::DocText::new(
                    ::alphal00p_docs_schema::DocFormat::PlainText,
                    #signature,
                ));
            }
        }
    }
}

fn external_type_name(ty: &Type) -> Option<String> {
    match ty {
        Type::Path(path) => path
            .path
            .segments
            .last()
            .map(|segment| segment.ident.to_string()),
        Type::Reference(reference) => external_type_name(&reference.elem),
        _ => None,
    }
}

fn external_type_descriptor(
    arguments: &CommonArgs,
    items: &[Item],
    source: &str,
    source_file: &str,
    source_id: &str,
    name: &str,
    marker: &syn::Ident,
) -> syn::Result<TokenStream2> {
    let descriptor = format_ident!("__alphal00p_docs_ty_{}", marker);
    let kind = external_kind_tokens(arguments, "type")?;
    let line = external_source_line(source, name, None, ExternalKind::Type)?;
    let context = SourceContext::External {
        identifier: source_id,
        file: source_file,
        line,
        column: 1,
    };
    for item in items {
        match item {
            Item::Struct(item) if item.ident == name => {
                let mut item = item.clone();
                let docs = prepare_external_docs(arguments, &mut item.attrs);
                let format = arguments.format_tokens()?;
                let members = fields(&mut item.fields, arguments, &format);
                let signature = item.to_token_stream().to_string();
                return Ok(item_descriptor(
                    arguments,
                    name.to_owned(),
                    descriptor,
                    kind,
                    &docs,
                    Some(signature),
                    quote!(descriptor.members = vec![#(#members),*];),
                    context,
                ));
            }
            Item::Enum(item) if item.ident == name => {
                let mut item = item.clone();
                let docs = prepare_external_docs(arguments, &mut item.attrs);
                let format = arguments.format_tokens()?;
                let variants = item
                    .variants
                    .iter_mut()
                    .map(|variant| {
                        let name = variant.ident.to_string();
                        let docs = prepare_member_docs(&mut variant.attrs, arguments, &format);
                        let fields = fields(&mut variant.fields, arguments, &format);
                        quote! {
                            {
                                let mut member = ::alphal00p_docs_schema::DocMember::new(
                                    #name,
                                    ::alphal00p_docs_schema::DocMemberKind::Variant,
                                );
                                #docs
                                member.members = vec![#(#fields),*];
                                member
                            }
                        }
                    })
                    .collect::<Vec<_>>();
                let signature = item.to_token_stream().to_string();
                return Ok(item_descriptor(
                    arguments,
                    name.to_owned(),
                    descriptor,
                    kind,
                    &docs,
                    Some(signature),
                    quote!(descriptor.members = vec![#(#variants),*];),
                    context,
                ));
            }
            Item::Type(item) if item.ident == name => {
                let mut item: ItemType = item.clone();
                let docs = prepare_external_docs(arguments, &mut item.attrs);
                let signature = item.to_token_stream().to_string();
                return Ok(item_descriptor(
                    arguments,
                    name.to_owned(),
                    descriptor,
                    kind,
                    &docs,
                    Some(signature),
                    quote!(),
                    context,
                ));
            }
            _ => {}
        }
    }
    Err(syn::Error::new(
        proc_macro2::Span::call_site(),
        format!("source type {source_id} was not found"),
    ))
}

fn external_trait_descriptor(
    arguments: &CommonArgs,
    items: &[Item],
    source: &str,
    source_file: &str,
    source_id: &str,
    name: &str,
    marker: &syn::Ident,
) -> syn::Result<TokenStream2> {
    let mut item = items
        .iter()
        .find_map(|item| match item {
            Item::Trait(item) if item.ident == name => Some(item.clone()),
            _ => None,
        })
        .ok_or_else(|| {
            syn::Error::new(
                proc_macro2::Span::call_site(),
                format!("source trait {source_id} was not found"),
            )
        })?;
    let docs = prepare_external_docs(arguments, &mut item.attrs);
    let format = arguments.format_tokens()?;
    let members = item
        .items
        .iter_mut()
        .filter_map(|item| match item {
            TraitItem::Fn(item) => {
                let kind = if item.sig.receiver().is_some() {
                    quote!(::alphal00p_docs_schema::DocMemberKind::Method)
                } else {
                    quote!(::alphal00p_docs_schema::DocMemberKind::AssociatedFunction)
                };
                let docs = prepare_member_docs(&mut item.attrs, arguments, &format);
                Some(trait_member(
                    item.sig.ident.to_string(),
                    kind,
                    item.to_token_stream().to_string(),
                    docs,
                ))
            }
            TraitItem::Type(item) => {
                let docs = prepare_member_docs(&mut item.attrs, arguments, &format);
                Some(trait_member(
                    item.ident.to_string(),
                    quote!(::alphal00p_docs_schema::DocMemberKind::AssociatedType),
                    item.to_token_stream().to_string(),
                    docs,
                ))
            }
            TraitItem::Const(item) => {
                let docs = prepare_member_docs(&mut item.attrs, arguments, &format);
                Some(trait_member(
                    item.ident.to_string(),
                    quote!(::alphal00p_docs_schema::DocMemberKind::AssociatedConst),
                    item.to_token_stream().to_string(),
                    docs,
                ))
            }
            _ => None,
        })
        .collect::<Vec<_>>();
    let signature = item.to_token_stream().to_string();
    let descriptor = format_ident!("__alphal00p_docs_trait_{}", marker);
    let kind = external_kind_tokens(arguments, "trait")?;
    let line = external_source_line(source, name, None, ExternalKind::Trait)?;
    Ok(item_descriptor(
        arguments,
        name.to_owned(),
        descriptor,
        kind,
        &docs,
        Some(signature),
        quote!(descriptor.members = vec![#(#members),*];),
        SourceContext::External {
            identifier: source_id,
            file: source_file,
            line,
            column: 1,
        },
    ))
}

fn external_macro_descriptor(
    arguments: &CommonArgs,
    items: &[Item],
    source: &str,
    source_file: &str,
    source_id: &str,
    name: &str,
    marker: &syn::Ident,
) -> syn::Result<TokenStream2> {
    let descriptor = format_ident!("__alphal00p_docs_macro_{}", marker);
    let kind = external_kind_tokens(arguments, "exported-macro")?;
    if let Some(mut item) = items.iter().find_map(|item| match item {
        Item::Macro(item) if item.ident.as_ref().is_some_and(|ident| ident == name) => {
            Some(item.clone())
        }
        _ => None,
    }) {
        let docs = prepare_external_docs(arguments, &mut item.attrs);
        let line = external_source_line(source, name, None, ExternalKind::Macro)?;
        return Ok(item_descriptor(
            arguments,
            name.to_owned(),
            descriptor,
            kind,
            &docs,
            Some(item.to_token_stream().to_string()),
            quote!(),
            SourceContext::External {
                identifier: source_id,
                file: source_file,
                line,
                column: 1,
            },
        ));
    }
    if let Some(mut item) = items.iter().find_map(|item| match item {
        Item::Fn(item) if exported_proc_macro_name(item).as_deref() == Some(name) => {
            Some(item.clone())
        }
        _ => None,
    }) {
        let docs = prepare_external_docs(arguments, &mut item.attrs);
        let line = external_proc_macro_source_line(source, &item, name)?;
        let signature = exported_proc_macro_signature(&item, name);
        return Ok(item_descriptor(
            arguments,
            name.to_owned(),
            descriptor,
            kind,
            &docs,
            Some(signature),
            quote!(),
            SourceContext::External {
                identifier: source_id,
                file: source_file,
                line,
                column: 1,
            },
        ));
    }
    Err(syn::Error::new(
        proc_macro2::Span::call_site(),
        format!("source macro {source_id} was not found"),
    ))
}

fn exported_proc_macro_name(item: &ItemFn) -> Option<String> {
    for attribute in &item.attrs {
        if attribute.path().is_ident("proc_macro_derive") {
            let tokens = attribute.meta.to_token_stream().to_string();
            let arguments = tokens.split_once('(')?.1.rsplit_once(')')?.0;
            return arguments
                .split(',')
                .next()
                .map(str::trim)
                .filter(|name| !name.is_empty())
                .map(ToOwned::to_owned);
        }
        if attribute.path().is_ident("proc_macro")
            || attribute.path().is_ident("proc_macro_attribute")
        {
            return Some(item.sig.ident.to_string());
        }
    }
    None
}

fn exported_proc_macro_signature(item: &ItemFn, exported_name: &str) -> String {
    let invocation = if item
        .attrs
        .iter()
        .any(|attribute| attribute.path().is_ident("proc_macro_derive"))
    {
        format!("#[derive({exported_name})]")
    } else if item
        .attrs
        .iter()
        .any(|attribute| attribute.path().is_ident("proc_macro_attribute"))
    {
        format!("#[{exported_name}]")
    } else {
        format!("{exported_name}!(...)")
    };
    format!("{invocation}\n{}", item.sig.to_token_stream())
}

fn external_kind_tokens(arguments: &CommonArgs, default: &str) -> syn::Result<TokenStream2> {
    let kind = arguments.kind.as_deref().unwrap_or(default);
    let variant = match kind {
        "function" => quote!(::alphal00p_docs_schema::DocItemKind::Function),
        "type" => quote!(::alphal00p_docs_schema::DocItemKind::Type),
        "trait" => quote!(::alphal00p_docs_schema::DocItemKind::Trait),
        "exported-macro" | "macro" => {
            quote!(::alphal00p_docs_schema::DocItemKind::ExportedMacro)
        }
        "method" => quote!(::alphal00p_docs_schema::DocItemKind::Method),
        "setting" => quote!(::alphal00p_docs_schema::DocItemKind::Setting),
        "command" => quote!(::alphal00p_docs_schema::DocItemKind::Command),
        _ => {
            return Err(syn::Error::new(
                proc_macro2::Span::call_site(),
                format!("unknown source-backed item kind {kind}"),
            ));
        }
    };
    Ok(variant)
}

fn external_source_line(
    source: &str,
    name: &str,
    owner: Option<&str>,
    kind: ExternalKind,
) -> syn::Result<u32> {
    let markers = match kind {
        ExternalKind::Function => vec![format!("fn {name}")],
        ExternalKind::Type => vec![
            format!("struct {name}"),
            format!("enum {name}"),
            format!("type {name}"),
        ],
        ExternalKind::Trait => vec![format!("trait {name}")],
        ExternalKind::Macro => vec![format!("macro_rules! {name}")],
    };
    let lines = external_code_lines(source);
    let mut candidates = lines
        .iter()
        .enumerate()
        .filter_map(|(index, line)| {
            markers
                .iter()
                .any(|marker| external_line_has_marker(line, marker))
                .then_some(index)
        })
        .collect::<Vec<_>>();
    if let Some(owner) = owner {
        candidates
            .retain(|candidate| external_candidate_belongs_to_impl(&lines, *candidate, owner));
    }
    if candidates.len() != 1 {
        return Err(syn::Error::new(
            proc_macro2::Span::call_site(),
            format!(
                "source declaration {name} resolved to {} candidate lines",
                candidates.len()
            ),
        ));
    }
    u32::try_from(candidates[0] + 1).map_err(|_| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("source declaration {name} has a line number larger than u32"),
        )
    })
}

fn external_proc_macro_source_line(
    source: &str,
    item: &ItemFn,
    exported_name: &str,
) -> syn::Result<u32> {
    let derive = item
        .attrs
        .iter()
        .any(|attribute| attribute.path().is_ident("proc_macro_derive"));
    let marker = if derive {
        format!("proc_macro_derive({exported_name}")
    } else {
        format!("fn {}", item.sig.ident)
    };
    let lines = external_code_lines(source);
    let mut candidates = lines
        .iter()
        .enumerate()
        .filter_map(|(index, line)| external_line_has_marker(line, &marker).then_some(index))
        .collect::<Vec<_>>();
    if !derive {
        for candidate in &mut candidates {
            *candidate = external_proc_macro_attribute_line(&lines, *candidate);
        }
        candidates.sort_unstable();
        candidates.dedup();
    }
    if candidates.len() != 1 {
        return Err(syn::Error::new(
            proc_macro2::Span::call_site(),
            format!(
                "source declaration {exported_name} resolved to {} candidate lines",
                candidates.len()
            ),
        ));
    }
    u32::try_from(candidates[0] + 1).map_err(|_| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("source declaration {exported_name} has a line number larger than u32"),
        )
    })
}

fn external_line_has_marker(line: &str, marker: &str) -> bool {
    line.match_indices(marker).any(|(start, _)| {
        let before = line[..start].chars().next_back();
        let after = line[start + marker.len()..].chars().next();
        before.is_none_or(|character| !(character.is_ascii_alphanumeric() || character == '_'))
            && after
                .is_none_or(|character| !(character.is_ascii_alphanumeric() || character == '_'))
    })
}

fn external_proc_macro_attribute_line(lines: &[String], function_line: usize) -> usize {
    for index in (0..=function_line).rev() {
        let line = lines[index].trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with("#[proc_macro") {
            return index;
        }
        if index != function_line && !line.starts_with("#[") {
            break;
        }
    }
    function_line
}

fn external_candidate_belongs_to_impl(lines: &[String], candidate: usize, owner: &str) -> bool {
    let Some(start) = (0..candidate).rev().find(|index| {
        let line = lines[*index].trim_start();
        line.starts_with("impl ") || line.starts_with("impl<")
    }) else {
        return false;
    };
    let mut header = String::new();
    for line in lines.iter().skip(start).take(8) {
        header.push_str(line.trim());
        header.push(' ');
        if line.contains('{') {
            break;
        }
    }
    header
        .split(|character: char| !(character.is_ascii_alphanumeric() || character == '_'))
        .any(|segment| segment == owner)
}

fn external_code_lines(source: &str) -> Vec<String> {
    let source = source.as_bytes();
    let mut code = source.to_vec();
    let mut index = 0;
    let mut line_comment = false;
    let mut block_comment_depth = 0_u32;
    let mut quoted_string = false;
    let mut raw_string_hashes = None;

    while index < source.len() {
        if line_comment {
            if source[index] == b'\n' {
                line_comment = false;
            } else {
                code[index] = b' ';
            }
            index += 1;
            continue;
        }
        if block_comment_depth > 0 {
            if source[index..].starts_with(b"/*") {
                code[index] = b' ';
                code[index + 1] = b' ';
                block_comment_depth += 1;
                index += 2;
            } else if source[index..].starts_with(b"*/") {
                code[index] = b' ';
                code[index + 1] = b' ';
                block_comment_depth -= 1;
                index += 2;
            } else {
                if source[index] != b'\n' && source[index] != b'\r' {
                    code[index] = b' ';
                }
                index += 1;
            }
            continue;
        }
        if let Some(hashes) = raw_string_hashes {
            if source[index] == b'"'
                && source
                    .get(index + 1..index + 1 + hashes)
                    .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
            {
                for byte in &mut code[index..index + 1 + hashes] {
                    *byte = b' ';
                }
                raw_string_hashes = None;
                index += 1 + hashes;
            } else {
                if source[index] != b'\n' && source[index] != b'\r' {
                    code[index] = b' ';
                }
                index += 1;
            }
            continue;
        }
        if quoted_string {
            if source[index] == b'\\' {
                code[index] = b' ';
                if let Some(next) = code.get_mut(index + 1)
                    && *next != b'\n'
                    && *next != b'\r'
                {
                    *next = b' ';
                }
                index += 2;
            } else {
                if source[index] == b'"' {
                    quoted_string = false;
                }
                if source[index] != b'\n' && source[index] != b'\r' {
                    code[index] = b' ';
                }
                index += 1;
            }
            continue;
        }

        if source[index..].starts_with(b"//") {
            code[index] = b' ';
            code[index + 1] = b' ';
            line_comment = true;
            index += 2;
        } else if source[index..].starts_with(b"/*") {
            code[index] = b' ';
            code[index + 1] = b' ';
            block_comment_depth = 1;
            index += 2;
        } else if let Some((delimiter_length, hashes)) = external_raw_string_start(source, index) {
            for byte in &mut code[index..index + delimiter_length] {
                *byte = b' ';
            }
            raw_string_hashes = Some(hashes);
            index += delimiter_length;
        } else if source[index] == b'"' {
            code[index] = b' ';
            quoted_string = true;
            index += 1;
        } else {
            index += 1;
        }
    }

    String::from_utf8(code)
        .expect("masking UTF-8 source with ASCII spaces preserves UTF-8")
        .lines()
        .map(ToOwned::to_owned)
        .collect()
}

fn external_raw_string_start(source: &[u8], index: usize) -> Option<(usize, usize)> {
    let mut cursor = match source.get(index..index + 2) {
        Some(b"br" | b"rb") => index + 2,
        _ if source.get(index) == Some(&b'r') => index + 1,
        _ => return None,
    };
    let hash_start = cursor;
    while source.get(cursor) == Some(&b'#') {
        cursor += 1;
    }
    (source.get(cursor) == Some(&b'"')).then_some((cursor - index + 1, cursor - hash_start))
}

/// Capture a function or inherent method in a documentation descriptor.
///
/// Methods with a `self` receiver infer their owner. A receiverless associated
/// function must declare its containing type, such as `owner = "Counter"` (or
/// `owner = "Wrapper<T>"` in a generic impl), because attribute macros cannot
/// otherwise distinguish it from a free function.
#[proc_macro_attribute]
pub fn func(arguments: TokenStream, input: TokenStream) -> TokenStream {
    let arguments = match CommonArgs::parse(arguments) {
        Ok(arguments) => arguments,
        Err(error) => return error.into_compile_error().into(),
    };
    let tokens = TokenStream2::from(input);
    match expand_external_item(&arguments, tokens.clone(), ExternalKind::Function) {
        Ok(Some(expanded)) => return expanded.into(),
        Ok(None) => {}
        Err(error) => return error.into_compile_error().into(),
    }
    let declared_owner = match arguments.owner_type() {
        Ok(owner) => owner,
        Err(error) => return error.into_compile_error().into(),
    };
    if let Ok(mut item) = syn::parse2::<ItemFn>(tokens.clone()) {
        let signature = item.sig.clone();
        let full_docs = prepare_item_docs(&arguments, &mut item.attrs);
        let source_context = if signature.receiver().is_some() || declared_owner.is_some() {
            SourceContext::Inherent {
                declared_owner: declared_owner.as_ref(),
            }
        } else {
            SourceContext::Module
        };
        let kind = match source_context {
            SourceContext::Inherent { .. } => {
                quote!(::alphal00p_docs_schema::DocItemKind::Method)
            }
            SourceContext::Module => quote!(::alphal00p_docs_schema::DocItemKind::Function),
            SourceContext::External { .. } => unreachable!("handled before normal expansion"),
        };
        let descriptor =
            function_descriptor(&arguments, &signature, &full_docs, kind, source_context);
        return quote!(#item #descriptor).into();
    }
    if let Ok(mut item) = syn::parse2::<ImplItemFn>(tokens) {
        let signature = item.sig.clone();
        let full_docs = prepare_item_docs(&arguments, &mut item.attrs);
        let descriptor = function_descriptor(
            &arguments,
            &signature,
            &full_docs,
            quote!(::alphal00p_docs_schema::DocItemKind::Method),
            SourceContext::Inherent {
                declared_owner: declared_owner.as_ref(),
            },
        );
        return quote!(#item #descriptor).into();
    }

    syn::Error::new(
        proc_macro2::Span::call_site(),
        "#[alphal00p_docs::func] requires a function or inherent method",
    )
    .into_compile_error()
    .into()
}

fn function_descriptor(
    arguments: &CommonArgs,
    signature: &syn::Signature,
    full_docs: &str,
    kind: TokenStream2,
    source_context: SourceContext<'_>,
) -> TokenStream2 {
    let ident = &signature.ident;
    let descriptor = format_ident!("__alphal00p_docs_func_{}", ident);
    let signature_text = signature.to_token_stream().to_string();
    let parameters = signature.inputs.iter().map(parameter);
    let returns = match &signature.output {
        ReturnType::Default => quote!(),
        ReturnType::Type(_, ty) => {
            let signature = ty.to_token_stream().to_string();
            quote! {
                descriptor.returns = Some(::alphal00p_docs_schema::DocText::new(
                    ::alphal00p_docs_schema::DocFormat::PlainText,
                    #signature,
                ));
            }
        }
    };
    item_descriptor(
        arguments,
        ident.to_string(),
        descriptor,
        kind,
        full_docs,
        Some(signature_text),
        quote! {
            descriptor.params = vec![#(#parameters),*];
            #returns
        },
        source_context,
    )
}

fn parameter(argument: &FnArg) -> TokenStream2 {
    match argument {
        FnArg::Receiver(receiver) => {
            let signature = receiver.to_token_stream().to_string();
            quote! {
                {
                    let mut parameter = ::alphal00p_docs_schema::DocParam::new("self");
                    parameter.ty = Some(#signature.to_owned());
                    parameter
                }
            }
        }
        FnArg::Typed(argument) => {
            let name = match argument.pat.as_ref() {
                Pat::Ident(ident) => ident.ident.to_string(),
                pattern => pattern.to_token_stream().to_string(),
            };
            let signature = argument.ty.to_token_stream().to_string();
            quote! {
                {
                    let mut parameter = ::alphal00p_docs_schema::DocParam::new(#name);
                    parameter.ty = Some(#signature.to_owned());
                    parameter
                }
            }
        }
    }
}

#[proc_macro_attribute]
pub fn ty(arguments: TokenStream, input: TokenStream) -> TokenStream {
    let arguments = match CommonArgs::parse(arguments) {
        Ok(arguments) => arguments,
        Err(error) => return error.into_compile_error().into(),
    };
    let tokens = TokenStream2::from(input);
    match expand_external_item(&arguments, tokens.clone(), ExternalKind::Type) {
        Ok(Some(expanded)) => return expanded.into(),
        Ok(None) => {}
        Err(error) => return error.into_compile_error().into(),
    }
    if let Err(error) = arguments.reject_owner("#[alphal00p_docs::ty]") {
        return error.into_compile_error().into();
    }

    if let Ok(mut item) = syn::parse2::<ItemStruct>(tokens.clone()) {
        let full_docs = prepare_item_docs(&arguments, &mut item.attrs);
        let format = match arguments.format_tokens() {
            Ok(format) => format,
            Err(error) => return error.into_compile_error().into(),
        };
        let ident = item.ident.clone();
        let descriptor = format_ident!("__alphal00p_docs_ty_{}", ident);
        let members = fields(&mut item.fields, &arguments, &format);
        let signature = item.to_token_stream().to_string();
        let generated = item_descriptor(
            &arguments,
            ident.to_string(),
            descriptor,
            quote!(::alphal00p_docs_schema::DocItemKind::Type),
            &full_docs,
            Some(signature),
            quote!(descriptor.members = vec![#(#members),*];),
            SourceContext::Module,
        );
        return quote!(#item #generated).into();
    }

    if let Ok(mut item) = syn::parse2::<ItemEnum>(tokens) {
        let full_docs = prepare_item_docs(&arguments, &mut item.attrs);
        let format = match arguments.format_tokens() {
            Ok(format) => format,
            Err(error) => return error.into_compile_error().into(),
        };
        let ident = item.ident.clone();
        let descriptor = format_ident!("__alphal00p_docs_ty_{}", ident);
        let variants = item
            .variants
            .iter_mut()
            .map(|variant| {
                let name = variant.ident.to_string();
                let docs = prepare_member_docs(&mut variant.attrs, &arguments, &format);
                let fields = fields(&mut variant.fields, &arguments, &format);
                quote! {
                    {
                        let mut member = ::alphal00p_docs_schema::DocMember::new(
                            #name,
                            ::alphal00p_docs_schema::DocMemberKind::Variant,
                        );
                        #docs
                        member.members = vec![#(#fields),*];
                        member
                    }
                }
            })
            .collect::<Vec<_>>();
        let signature = item.to_token_stream().to_string();
        let generated = item_descriptor(
            &arguments,
            ident.to_string(),
            descriptor,
            quote!(::alphal00p_docs_schema::DocItemKind::Type),
            &full_docs,
            Some(signature),
            quote!(descriptor.members = vec![#(#variants),*];),
            SourceContext::Module,
        );
        return quote!(#item #generated).into();
    }

    syn::Error::new(
        proc_macro2::Span::call_site(),
        "#[alphal00p_docs::ty] requires a struct or enum",
    )
    .into_compile_error()
    .into()
}

fn fields(fields: &mut Fields, arguments: &CommonArgs, format: &TokenStream2) -> Vec<TokenStream2> {
    fields
        .iter_mut()
        .enumerate()
        .map(|(index, field)| field_descriptor(index, field, arguments, format))
        .collect()
}

fn field_descriptor(
    index: usize,
    field: &mut Field,
    arguments: &CommonArgs,
    format: &TokenStream2,
) -> TokenStream2 {
    let name = field
        .ident
        .as_ref()
        .map(ToString::to_string)
        .unwrap_or_else(|| index.to_string());
    let signature = field.ty.to_token_stream().to_string();
    let docs = prepare_member_docs(&mut field.attrs, arguments, format);
    quote! {
        {
            let mut member = ::alphal00p_docs_schema::DocMember::new(
                #name,
                ::alphal00p_docs_schema::DocMemberKind::Field,
            );
            member.signature = Some(#signature.to_owned());
            #docs
            member
        }
    }
}

fn prepare_member_docs(
    attributes: &mut Vec<Attribute>,
    arguments: &CommonArgs,
    format: &TokenStream2,
) -> TokenStream2 {
    let full_docs = if arguments.is_typst_markup() {
        let full_docs = take_docs(attributes);
        replace_docs(attributes, &first_sentence(&full_docs));
        full_docs
    } else {
        docs(attributes)
    };
    if full_docs.is_empty() {
        quote!()
    } else {
        quote! {
            member.docs = Some(::alphal00p_docs_schema::DocText::new(
                #format,
                #full_docs,
            ));
        }
    }
}

#[proc_macro_attribute]
pub fn trait_item(arguments: TokenStream, input: TokenStream) -> TokenStream {
    let arguments = match CommonArgs::parse(arguments) {
        Ok(arguments) => arguments,
        Err(error) => return error.into_compile_error().into(),
    };
    let tokens = TokenStream2::from(input.clone());
    match expand_external_item(&arguments, tokens, ExternalKind::Trait) {
        Ok(Some(expanded)) => return expanded.into(),
        Ok(None) => {}
        Err(error) => return error.into_compile_error().into(),
    }
    if let Err(error) = arguments.reject_owner("#[alphal00p_docs::trait_item]") {
        return error.into_compile_error().into();
    }
    let mut item = parse_macro_input!(input as ItemTrait);
    let full_docs = prepare_item_docs(&arguments, &mut item.attrs);
    let format = match arguments.format_tokens() {
        Ok(format) => format,
        Err(error) => return error.into_compile_error().into(),
    };
    let ident = item.ident.clone();
    let descriptor = format_ident!("__alphal00p_docs_trait_{}", ident);
    let members = item
        .items
        .iter_mut()
        .filter_map(|item| match item {
            TraitItem::Fn(item) => {
                let kind = if item.sig.receiver().is_some() {
                    quote!(::alphal00p_docs_schema::DocMemberKind::Method)
                } else {
                    quote!(::alphal00p_docs_schema::DocMemberKind::AssociatedFunction)
                };
                let docs = prepare_member_docs(&mut item.attrs, &arguments, &format);
                Some(trait_member(
                    item.sig.ident.to_string(),
                    kind,
                    item.to_token_stream().to_string(),
                    docs,
                ))
            }
            TraitItem::Type(item) => {
                let docs = prepare_member_docs(&mut item.attrs, &arguments, &format);
                Some(trait_member(
                    item.ident.to_string(),
                    quote!(::alphal00p_docs_schema::DocMemberKind::AssociatedType),
                    item.to_token_stream().to_string(),
                    docs,
                ))
            }
            TraitItem::Const(item) => {
                let docs = prepare_member_docs(&mut item.attrs, &arguments, &format);
                Some(trait_member(
                    item.ident.to_string(),
                    quote!(::alphal00p_docs_schema::DocMemberKind::AssociatedConst),
                    item.to_token_stream().to_string(),
                    docs,
                ))
            }
            _ => None,
        })
        .collect::<Vec<_>>();
    let signature = item.to_token_stream().to_string();
    let generated = item_descriptor(
        &arguments,
        ident.to_string(),
        descriptor,
        quote!(::alphal00p_docs_schema::DocItemKind::Trait),
        &full_docs,
        Some(signature),
        quote!(descriptor.members = vec![#(#members),*];),
        SourceContext::Module,
    );
    quote!(#item #generated).into()
}

fn trait_member(
    name: String,
    kind: TokenStream2,
    signature: String,
    docs: TokenStream2,
) -> TokenStream2 {
    quote! {
        {
            let mut member = ::alphal00p_docs_schema::DocMember::new(#name, #kind);
            member.signature = Some(#signature.to_owned());
            #docs
            member
        }
    }
}

#[proc_macro_attribute]
pub fn macro_item(arguments: TokenStream, input: TokenStream) -> TokenStream {
    let arguments = match CommonArgs::parse(arguments) {
        Ok(arguments) => arguments,
        Err(error) => return error.into_compile_error().into(),
    };
    let tokens = TokenStream2::from(input);
    match expand_external_item(&arguments, tokens.clone(), ExternalKind::Macro) {
        Ok(Some(expanded)) => return expanded.into(),
        Ok(None) => {}
        Err(error) => return error.into_compile_error().into(),
    }
    if let Err(error) = arguments.reject_owner("#[alphal00p_docs::macro_item]") {
        return error.into_compile_error().into();
    }

    if let Ok(mut item) = syn::parse2::<ItemMacro>(tokens.clone()) {
        let Some(ident) = item.ident.clone() else {
            return macro_item_target_error(item).into_compile_error().into();
        };
        let full_docs = prepare_item_docs(&arguments, &mut item.attrs);
        let descriptor = format_ident!("__alphal00p_docs_macro_{}", ident);
        let signature = item.to_token_stream().to_string();
        let generated = item_descriptor(
            &arguments,
            ident.to_string(),
            descriptor,
            quote!(::alphal00p_docs_schema::DocItemKind::ExportedMacro),
            &full_docs,
            Some(signature),
            quote!(),
            SourceContext::Module,
        );
        return quote!(#item #generated).into();
    }

    if let Ok(item) = syn::parse2::<ItemFn>(tokens) {
        return match expand_proc_macro_function(&arguments, item) {
            Ok(expanded) => expanded.into(),
            Err(error) => error.into_compile_error().into(),
        };
    }

    macro_item_target_error(quote!(#[macro_item]))
        .into_compile_error()
        .into()
}

fn expand_proc_macro_function(
    arguments: &CommonArgs,
    mut item: ItemFn,
) -> syn::Result<TokenStream2> {
    if !is_exported_proc_macro(&item.attrs) {
        return Err(macro_item_target_error(&item.sig.ident));
    }

    let full_docs = prepare_item_docs(arguments, &mut item.attrs);
    let ident = item.sig.ident.clone();
    let descriptor = format_ident!("__alphal00p_docs_macro_{}", ident);
    let signature = item.sig.to_token_stream().to_string();
    let generated = item_descriptor_with_visibility(
        arguments,
        ident.to_string(),
        descriptor,
        quote!(::alphal00p_docs_schema::DocItemKind::ExportedMacro),
        &full_docs,
        Some(signature),
        quote!(),
        quote!(pub(crate)),
        SourceContext::Module,
    );
    Ok(quote!(#item #generated))
}

fn is_exported_proc_macro(attributes: &[Attribute]) -> bool {
    attributes.iter().any(|attribute| {
        attribute.path().is_ident("proc_macro")
            || attribute.path().is_ident("proc_macro_attribute")
            || attribute.path().is_ident("proc_macro_derive")
    })
}

fn macro_item_target_error(tokens: impl ToTokens) -> syn::Error {
    syn::Error::new_spanned(
        tokens,
        "#[alphal00p_docs::macro_item] requires macro_rules! or a function marked \
         #[proc_macro], #[proc_macro_attribute], or #[proc_macro_derive]",
    )
}

#[proc_macro_attribute]
pub fn scope(arguments: TokenStream, input: TokenStream) -> TokenStream {
    let arguments = match CommonArgs::parse(arguments) {
        Ok(arguments) => arguments,
        Err(error) => return error.into_compile_error().into(),
    };
    if let Err(error) = arguments.reject_owner("#[alphal00p_docs::scope]") {
        return error.into_compile_error().into();
    }
    if let Err(error) = arguments.reject_source_item("#[alphal00p_docs::scope]") {
        return error.into_compile_error().into();
    }
    let mut item = parse_macro_input!(input as ItemMod);
    let full_docs = prepare_item_docs(&arguments, &mut item.attrs);
    let ident = item.ident.clone();
    let descriptor = format_ident!("__alphal00p_docs_scope_{}", ident);
    let default_id = ident.to_string();
    let id = arguments.id.as_deref().unwrap_or(&default_id);
    let title = arguments
        .title
        .as_deref()
        .or(arguments.name.as_deref())
        .unwrap_or(id);
    let summary = arguments
        .summary
        .clone()
        .unwrap_or_else(|| first_sentence(&full_docs));
    let format = match arguments.format_tokens() {
        Ok(format) => format,
        Err(error) => return error.into_compile_error().into(),
    };
    let docs = docs_assignment("descriptor", &full_docs, &format);
    let summary = string_assignment(
        "descriptor",
        "summary",
        (!summary.is_empty()).then_some(summary.as_str()),
    );
    let since = string_assignment("descriptor", "since", arguments.since.as_deref());
    let keywords = &arguments.keywords;

    quote! {
        #item

        #[doc(hidden)]
        #[allow(non_snake_case)]
        pub fn #descriptor() -> ::alphal00p_docs_schema::DocScope {
            let mut descriptor = ::alphal00p_docs_schema::DocScope::new(#id, #title);
            #docs
            #summary
            #since
            descriptor.keywords = vec![#(#keywords.to_owned()),*];
            descriptor.source = Some(::alphal00p_docs_schema::SourceLocation::new(
                format!("{}::{}", module_path!(), #id),
                file!(),
                line!(),
                column!(),
            ));
            descriptor
        }
    }
    .into()
}

fn prepare_item_docs(arguments: &CommonArgs, attributes: &mut Vec<Attribute>) -> String {
    if arguments.is_typst_markup() {
        let full_docs = take_docs(attributes);
        let summary = arguments
            .summary
            .clone()
            .unwrap_or_else(|| first_sentence(&full_docs));
        replace_docs(attributes, &summary);
        full_docs
    } else {
        docs(attributes)
    }
}

/// Existing product Rustdoc is authoritative for source-backed descriptors.
/// Some supported items predate comprehensive comments, so the adapter's
/// required summary also serves as their complete prose until the source gains
/// its own documentation. This keeps strict catalogs useful without inventing
/// signatures or source locations in a parallel table.
fn prepare_external_docs(arguments: &CommonArgs, attributes: &mut Vec<Attribute>) -> String {
    let docs = prepare_item_docs(arguments, attributes);
    if docs.trim().is_empty() {
        arguments.summary.clone().unwrap_or_default()
    } else {
        docs
    }
}

#[allow(clippy::too_many_arguments)]
fn item_descriptor(
    arguments: &CommonArgs,
    default_name: String,
    descriptor: syn::Ident,
    kind: TokenStream2,
    full_docs: &str,
    signature: Option<String>,
    members: TokenStream2,
    source_context: SourceContext<'_>,
) -> TokenStream2 {
    item_descriptor_with_visibility(
        arguments,
        default_name,
        descriptor,
        kind,
        full_docs,
        signature,
        members,
        quote!(pub),
        source_context,
    )
}

#[allow(clippy::too_many_arguments)]
fn item_descriptor_with_visibility(
    arguments: &CommonArgs,
    default_name: String,
    descriptor: syn::Ident,
    kind: TokenStream2,
    full_docs: &str,
    signature: Option<String>,
    members: TokenStream2,
    visibility: TokenStream2,
    source_context: SourceContext<'_>,
) -> TokenStream2 {
    let id = arguments.id.as_deref().unwrap_or(&default_name);
    let name = arguments.name.as_deref().unwrap_or(&default_name);
    let title = arguments.title.as_deref().unwrap_or(name);
    let summary = arguments
        .summary
        .clone()
        .unwrap_or_else(|| first_sentence(full_docs));
    let format = match arguments.format_tokens() {
        Ok(format) => format,
        Err(error) => return error.into_compile_error(),
    };
    let docs = docs_assignment("descriptor", full_docs, &format);
    let summary = string_assignment(
        "descriptor",
        "summary",
        (!summary.is_empty()).then_some(summary.as_str()),
    );
    let signature = string_assignment("descriptor", "signature", signature.as_deref());
    let since = string_assignment("descriptor", "since", arguments.since.as_deref());
    let keywords = &arguments.keywords;
    let examples = fenced_examples(full_docs);
    let examples = examples
        .iter()
        .enumerate()
        .map(|(index, (language, code))| {
            let title = format!("Example {}", index + 1);
            quote! {
                ::alphal00p_docs_schema::DocExample::new(#title, #language, #code)
            }
        });
    let source_name = default_name.clone();
    let identity = match source_context {
        SourceContext::Module => quote! {
            (
                #id.to_owned(),
                format!("{}::{}", module_path!(), #source_name),
                file!().to_owned(),
                line!(),
                column!(),
            )
        },
        SourceContext::Inherent { declared_owner } => {
            let owner_check = declared_owner.map_or_else(TokenStream2::new, |declared_owner| {
                quote! {
                    trait __Alphal00pDocsOwnerMatches<Declared: ?Sized> {}
                    impl<Owner: ?Sized> __Alphal00pDocsOwnerMatches<Owner> for Owner {}
                    fn __alphal00p_docs_assert_owner<Actual: ?Sized, Declared: ?Sized>()
                    where
                        Actual: __Alphal00pDocsOwnerMatches<Declared>,
                    {}
                    __alphal00p_docs_assert_owner::<Self, #declared_owner>();
                }
            });
            quote! {
                {
                    #owner_check
                    let owner = ::core::any::type_name::<Self>();
                    let owner = owner.split_once('<').map_or(owner, |(owner, _)| owner);
                    let short_owner = owner.rsplit("::").next().unwrap_or(owner);
                    (
                        format!("{}::{}", short_owner, #id),
                        format!("{}::{}", owner, #source_name),
                        file!().to_owned(),
                        line!(),
                        column!(),
                    )
                }
            }
        }
        SourceContext::External {
            identifier,
            file,
            line,
            column,
        } => quote! {
            (
                #id.to_owned(),
                #identifier.to_owned(),
                #file.to_owned(),
                #line,
                #column,
            )
        },
    };

    quote! {
        #[doc(hidden)]
        #[allow(non_snake_case)]
        #visibility fn #descriptor() -> ::alphal00p_docs_schema::DocItem {
            let (id, source_identifier, source_file, source_line, source_column) = #identity;
            let mut descriptor =
                ::alphal00p_docs_schema::DocItem::new(id, #name, #title, #kind);
            #docs
            #summary
            #signature
            #since
            descriptor.keywords = vec![#(#keywords.to_owned()),*];
            descriptor.examples = vec![#(#examples),*];
            #members
            descriptor.source = Some(::alphal00p_docs_schema::SourceLocation::new(
                source_identifier,
                source_file,
                source_line,
                source_column,
            ));
            descriptor
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        CommonArgs, ExternalKind, expand_proc_macro_function, external_proc_macro_source_line,
        external_source_line, is_exported_proc_macro, prepare_item_docs, prepare_member_docs,
    };
    use crate::args::docs;
    use quote::ToTokens;
    use syn::{Item, ItemFn, ItemStruct, Visibility, parse_quote};

    #[test]
    fn only_typst_item_docs_are_shortened_for_rustdoc() {
        let full_docs =
            "First sentence. Second sentence.\n\n# Details\nThe original body remains intact.";

        for format in ["rust-markdown", "python-docstring", "plain-text"] {
            let arguments = CommonArgs {
                format: Some(format.to_owned()),
                summary: Some("Catalog summary.".to_owned()),
                ..CommonArgs::default()
            };
            let mut item: ItemFn = parse_quote! {
                #[doc = "First sentence. Second sentence."]
                #[doc = ""]
                #[doc = "# Details"]
                #[doc = "The original body remains intact."]
                fn documented() {}
            };

            assert_eq!(prepare_item_docs(&arguments, &mut item.attrs), full_docs);
            assert_eq!(docs(&item.attrs), full_docs);
        }

        let mut item: ItemFn = parse_quote! {
            #[doc = "First sentence. Second sentence."]
            #[doc = ""]
            #[doc = "# Details"]
            #[doc = "The original body remains intact."]
            fn documented() {}
        };
        assert_eq!(
            prepare_item_docs(&CommonArgs::default(), &mut item.attrs),
            full_docs
        );
        assert_eq!(docs(&item.attrs), "First sentence.");
    }

    #[test]
    fn only_typst_member_docs_are_shortened_for_rustdoc() {
        let full_docs = "Field first sentence. Field second sentence.\n\nField details.";

        for (arguments, expected) in [
            (CommonArgs::default(), "Field first sentence."),
            (
                CommonArgs {
                    format: Some("rust-markdown".to_owned()),
                    ..CommonArgs::default()
                },
                full_docs,
            ),
            (
                CommonArgs {
                    format: Some("python-docstring".to_owned()),
                    ..CommonArgs::default()
                },
                full_docs,
            ),
        ] {
            let mut item: ItemStruct = parse_quote! {
                struct Documented {
                    #[doc = "Field first sentence. Field second sentence."]
                    #[doc = ""]
                    #[doc = "Field details."]
                    field: u8,
                }
            };
            let format = arguments.format_tokens().unwrap();
            let _ = prepare_member_docs(
                &mut item.fields.iter_mut().next().unwrap().attrs,
                &arguments,
                &format,
            );

            assert_eq!(docs(&item.fields.iter().next().unwrap().attrs), expected);
        }
    }

    #[test]
    fn recognizes_all_exported_proc_macro_function_attributes() {
        for item in [
            parse_quote!(
                #[proc_macro]
                pub fn bang(input: Tokens) -> Tokens {
                    input
                }
            ),
            parse_quote!(
                #[proc_macro_attribute]
                pub fn attr(_: Tokens, item: Tokens) -> Tokens {
                    item
                }
            ),
            parse_quote!(
                #[proc_macro_derive(DeriveMe)]
                pub fn derive(input: Tokens) -> Tokens {
                    input
                }
            ),
        ] {
            let item: ItemFn = item;
            assert!(is_exported_proc_macro(&item.attrs));
        }

        let plain: ItemFn = parse_quote!(
            pub fn helper() {}
        );
        assert!(!is_exported_proc_macro(&plain.attrs));
    }

    #[test]
    fn proc_macro_function_gets_crate_visible_exported_macro_descriptor() {
        let item: ItemFn = parse_quote! {
            #[proc_macro_derive(DeriveMe)]
            /// Derives a value. #emph[Typst details.]
            pub fn derive(input: Tokens) -> Tokens { input }
        };
        let expanded = expand_proc_macro_function(&CommonArgs::default(), item).unwrap();
        let file: syn::File = syn::parse2(expanded).unwrap();
        let Item::Fn(descriptor) = file.items.last().unwrap() else {
            panic!("expected a descriptor function");
        };

        assert_eq!(descriptor.sig.ident, "__alphal00p_docs_macro_derive");
        assert!(matches!(descriptor.vis, Visibility::Restricted(_)));
        assert!(
            descriptor
                .to_token_stream()
                .to_string()
                .contains("DocItemKind :: ExportedMacro")
        );
    }

    #[test]
    fn external_source_line_resolves_macro_kinds_and_ignores_comments() {
        let source = r#"
// macro_rules! documented {
//     () => {};
// }
/*
#[proc_macro]
pub fn bang(input: Tokens) -> Tokens { input }

#[proc_macro_attribute]
pub fn traced(_: Tokens, item: Tokens) -> Tokens { item }

#[proc_macro_derive(DeriveMe)]
pub fn derive_me(input: Tokens) -> Tokens { input }
*/
const _: &str = "macro_rules! documented { fake }";

#[macro_export]
macro_rules! documented {
    () => {};
}
pub fn documented() {}

#[proc_macro]
pub fn bang(input: Tokens) -> Tokens { input }

#[proc_macro_attribute]
pub fn traced(_: Tokens, item: Tokens) -> Tokens { item }

#[proc_macro_derive(DeriveMe, attributes(skip_me))]
pub fn derive_me(input: Tokens) -> Tokens { input }
"#;
        let file = syn::parse_file(source).unwrap();
        let function = |name: &str| {
            file.items
                .iter()
                .find_map(|item| match item {
                    Item::Fn(item) if item.sig.ident == name => Some(item),
                    _ => None,
                })
                .unwrap()
        };
        let expected_line = |declaration: &str| {
            source
                .lines()
                .enumerate()
                .filter_map(|(index, line)| (line == declaration).then_some(index))
                .last()
                .map(|index| u32::try_from(index + 1).unwrap())
                .unwrap()
        };

        assert_eq!(
            external_source_line(source, "documented", None, ExternalKind::Macro).unwrap(),
            expected_line("macro_rules! documented {")
        );
        assert_eq!(
            external_proc_macro_source_line(source, function("bang"), "bang").unwrap(),
            expected_line("#[proc_macro]")
        );
        assert_eq!(
            external_proc_macro_source_line(source, function("traced"), "traced").unwrap(),
            expected_line("#[proc_macro_attribute]")
        );
        assert_eq!(
            external_proc_macro_source_line(source, function("derive_me"), "DeriveMe").unwrap(),
            expected_line("#[proc_macro_derive(DeriveMe, attributes(skip_me))]")
        );
    }
}
