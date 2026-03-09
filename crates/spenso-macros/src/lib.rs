use proc_macro::TokenStream;
use proc_macro2::Span;
use quote::{format_ident, quote, ToTokens}; // Ensure ToTokens is imported
use syn::{
    parse::Parser,
    parse_macro_input,
    punctuated::Punctuated,
    spanned::Spanned,
    Attribute,
    Data,
    DeriveInput,
    Expr,
    ExprLit,
    ExprPath,
    Ident,
    Lit,
    LitStr,
    Meta,
    Path, // Use Path and Meta
    Token,
};

// --- RepresentationAttrs struct and parse_representation_attributes function (keep as is) ---
#[derive(Debug)]
struct RepresentationAttrs {
    name: LitStr,
    is_self_dual: bool,
    custom_dual_name: Option<Ident>,
}

fn parse_representation_attributes(attrs: &[Attribute]) -> Result<RepresentationAttrs, syn::Error> {
    // ... (implementation is correct) ...
    let mut rep_name: Option<LitStr> = None;
    let mut is_self_dual = false;
    let mut custom_dual_name: Option<Ident> = None;

    let rep_attr = attrs
        .iter()
        .find(|attr| attr.path().is_ident("representation"))
        .ok_or_else(|| {
            syn::Error::new(
                Span::call_site(), // Consider spanning the struct ident if possible
                "Missing #[representation(...)] attribute",
            )
        })?;

    let meta = &rep_attr.meta;
    let list = match meta {
        Meta::List(list) => list,
        _ => {
            return Err(syn::Error::new_spanned(
                meta,
                "Expected #[representation(...)] format",
            ))
        }
    };

    let parser = Punctuated::<Meta, Token![,]>::parse_terminated;
    let nested_metas = parser.parse2(list.tokens.clone()).map_err(|e| {
        syn::Error::new(
            e.span(),
            format!("Failed to parse attribute arguments: {}", e),
        )
    })?;

    for meta_item in nested_metas.iter() {
        match meta_item {
            Meta::NameValue(nv) if nv.path.is_ident("name") => {
                if rep_name.is_some() {
                    return Err(syn::Error::new_spanned(nv, "Duplicate `name` specified"));
                }
                if let Expr::Lit(ExprLit {
                    lit: Lit::Str(lit_str),
                    ..
                }) = &nv.value
                {
                    rep_name = Some(lit_str.clone());
                } else {
                    return Err(syn::Error::new_spanned(
                        &nv.value,
                        "Expected string literal for `name`",
                    ));
                }
            }
            Meta::Path(path) if path.is_ident("self_dual") => {
                if is_self_dual {
                    return Err(syn::Error::new_spanned(
                        path,
                        "Duplicate `self_dual` specified",
                    ));
                }
                is_self_dual = true;
            }
            Meta::NameValue(nv) if nv.path.is_ident("dual_name") => {
                if custom_dual_name.is_some() {
                    return Err(syn::Error::new_spanned(
                        nv,
                        "Duplicate `dual_name` specified",
                    ));
                }
                match &nv.value {
                    Expr::Lit(ExprLit {
                        lit: Lit::Str(lit_str),
                        ..
                    }) => {
                        custom_dual_name = Some(Ident::new(&lit_str.value(), lit_str.span()));
                    }
                    Expr::Path(ExprPath { path, .. }) => {
                        if let Some(ident) = path.get_ident() {
                            custom_dual_name = Some(ident.clone());
                        } else {
                            return Err(syn::Error::new_spanned(
                                &nv.value,
                                "Expected simple identifier for `dual_name` (e.g., MyDualName)",
                            ));
                        }
                    }
                    _ => {
                        return Err(syn::Error::new_spanned(
                            &nv.value,
                            "Expected string literal or identifier for `dual_name`",
                        ));
                    }
                }
            }
            _ => {
                return Err(syn::Error::new_spanned(
                    meta_item,
                    "Unsupported item in #[representation(...)] attribute",
                ));
            }
        }
    }

    let name = rep_name.ok_or_else(|| {
        syn::Error::new_spanned(
            list.tokens.clone(),
            "Missing required `name = \"...\"` in #[representation(...)]",
        )
    })?;

    if is_self_dual && custom_dual_name.is_some() {
        let error_span = nested_metas
            .iter()
            .find(|m| matches!(m, Meta::NameValue(nv) if nv.path.is_ident("dual_name")))
            .map_or_else(|| list.tokens.span(), |m| m.span());

        return Err(syn::Error::new(
            error_span,
            "`dual_name` cannot be specified for a `self_dual` representation",
        ));
    }

    Ok(RepresentationAttrs {
        name,
        is_self_dual,
        custom_dual_name,
    })
}

// --- Revised get_filtered_derive_paths using Meta parsing ---
fn get_filtered_derive_paths(attrs: &[Attribute]) -> Result<Vec<Path>, syn::Error> {
    let mut derived_traits = Vec::new();

    for attr in attrs {
        // Ensure it's the derive attribute: #[derive(...)]
        if attr.path().is_ident("derive") {
            // Use the Meta parser for the arguments inside derive(...)
            match attr.parse_args_with(Punctuated::<Meta, Token![,]>::parse_terminated) {
                Ok(nested_metas) => {
                    // eprintln!("    Parsed derive args as Metas successfully: {:#?}", nested_metas.iter().map(|m| m.to_token_stream().to_string()).collect::<Vec<_>>());
                    for meta in nested_metas {
                        // Standard derives like Debug, Clone, Serialize should appear as Meta::Path
                        if let Meta::Path(path) = meta {
                            // Check if the last segment of the path is "SimpleRepresentation"
                            let is_target_derive = path
                                .segments
                                .last()
                                .is_some_and(|segment| segment.ident == "SimpleRepresentation");

                            if !is_target_derive {
                                derived_traits.push(path); // Keep the original Path struct
                            }
                        } else {
                            // If we find something else (like Meta::List or Meta::NameValue)
                            // in a derive list, it's unexpected for standard traits.
                            // We could ignore it or error. Let's error for clarity.
                            return Err(syn::Error::new_spanned(
                                 meta, // Span the problematic meta item
                                 "Expected simple trait paths (e.g., Debug, Clone) in derive attribute, found other meta item.",
                             ));
                        }
                    }
                }
                Err(e) => {
                    // If parsing as Metas fails, the derive attribute format is likely malformed.
                    // Return the error, pointing at the problematic attribute
                    return Err(syn::Error::new_spanned(
                         attr.to_token_stream(), // Span the whole attribute
                         format!("Failed to parse derive arguments: {}. Check syntax inside #[derive(...)].", e),
                     ));
                }
            }
        }
    }

    Ok(derived_traits)
}

#[proc_macro_derive(SimpleRepresentation, attributes(representation))]
pub fn derive_simple_representation(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);

    // --- Extract Struct Info ---
    let fields = match &input.data {
        Data::Struct(s) => s.fields.clone(),
        _ => {
            return syn::Error::new_spanned(
                &input.ident,
                "SimpleRepresentation can only be derived for structs",
            )
            .to_compile_error()
            .into();
        }
    };

    // --- Parse Attributes ---
    let vis = &input.vis;
    let repr_attrs = match parse_representation_attributes(&input.attrs) {
        Ok(attrs) => attrs,
        Err(e) => return e.to_compile_error().into(),
    };
    // Use the revised helper function
    let derived_traits = match get_filtered_derive_paths(&input.attrs) {
        Ok(traits) => traits,
        Err(e) => return e.to_compile_error().into(),
    };

    let base_type_ident = &input.ident;
    let name_lit = &repr_attrs.name;
    let is_self_dual = repr_attrs.is_self_dual;

    // --- Bounds needed for RepName impls ---
    let base_bounds = quote! { Default + Copy };
    let dual_bounds = quote! { Default + Copy };

    // --- Common RepName Implementation Parts ---
    let base_repname_common_impl = quote! {
        #[inline]
        fn from_library_rep(rep: ::spenso::structure::representation::LibraryRep) -> ::std::result::Result<Self, ::spenso::structure::representation::RepresentationError>{
            rep.try_into()
        }
        #[inline] fn base(&self) -> Self::Base where Self::Base: Default { Self::Base::default() }
        #[inline] fn is_base(&self) -> bool { ::std::any::TypeId::of::<Self>() == ::std::any::TypeId::of::<Self::Base>() }
    };

    // --- Display Impls ---
    let base_display_impl = quote! {
        impl ::std::fmt::Display for #base_type_ident where #base_type_ident: Copy + Into<::spenso::structure::representation::LibraryRep> {
            fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result { write!(f, "{}", ::spenso::structure::representation::LibraryRep::from(*self)) }
        }
    };

    // --- Branch based on duality ---
    let expanded = if is_self_dual {
        // ===== Self-Dual Case ===== (Remains Correct)
        let rep_new_call =
            quote! { ::spenso::structure::representation::LibraryRep::new_self_dual(#name_lit) };
        let base_from_impl = quote! {
            impl From<#base_type_ident> for ::spenso::structure::representation::LibraryRep
                where #base_type_ident: Copy
                {
                    fn from(_value: #base_type_ident) -> Self {
                        #rep_new_call.expect(concat!("Failed to create self-dual Rep for ", #name_lit))
                    }
                }
        };

        let base_try_from_impl = quote! {
            impl TryFrom<::spenso::structure::representation::LibraryRep> for #base_type_ident where #base_type_ident: Default {
                type Error = ::spenso::structure::representation::RepresentationError;

                fn try_from(rep: ::spenso::structure::representation::LibraryRep) -> ::std::result::Result<Self, Self::Error>   {
                    let expected_rep = #rep_new_call.expect(concat!("Failed to create self-dual Rep for ", #name_lit));
                    if rep == expected_rep {
                        ::std::result::Result::Ok(#base_type_ident::default())
                    } else {
                        ::std::result::Result::Err(::spenso::structure::representation::RepresentationError::WrongRepresentationError(#name_lit.to_owned(), rep.to_string()))
                    }
                }
            }
        };

        let base_repname_impl = quote! {
            impl ::spenso::structure::representation::RepName for #base_type_ident where #base_type_ident: #base_bounds {
                type Base = #base_type_ident;
                type Dual = #base_type_ident;

                #[inline]
                fn orientation(self) -> ::linnet::half_edge::involution::Orientation {
                    ::linnet::half_edge::involution::Orientation::Undirected
                }

                #base_repname_common_impl
                #[inline]
                fn is_dual(self) -> bool { true }
                #[inline] fn matches(&self, _other: &Self::Dual) -> bool { true }
                #[inline] fn dual(self) -> Self::Dual { self }
            }
        };
        quote! {
            impl #base_type_ident {
                pub const NAME: &'static str = #name_lit;
            }
            #base_from_impl
            #base_try_from_impl
            #base_repname_impl
            #base_display_impl
        }
    } else {
        // ===== Dualizable Case =====

        // 1. Determine the identifier for the Dual struct
        let dual_type_ident = match &repr_attrs.custom_dual_name {
            Some(custom_name) => custom_name.clone(),
            None => format_ident!("Dual{}", base_type_ident, span = base_type_ident.span()),
        };

        // 2. Generate the Dual struct definition (using the filtered derives)
        let derive_attr = if !derived_traits.is_empty() {
            // Pass the collected Vec<Path> directly to quote
            quote! { #[derive( #(#derived_traits),* )] }
        } else {
            quote! {}
        };
        let dual_struct_def = quote! {
            #derive_attr
            #vis struct #dual_type_ident #fields
        };

        // 3. Define Rep creation calls
        let rep_new_base_call =
            quote! { ::spenso::structure::representation::LibraryRep::new_dual(#name_lit) };
        let rep_new_dual_call = quote! { #rep_new_base_call.expect(concat!("Failed to create dual Rep for ", #name_lit)).dual() };

        // 4. Implementations for the BASE type (Remains Correct)
        let base_from_impl = quote! {
            impl From<#base_type_ident> for ::spenso::structure::representation::LibraryRep where #base_type_ident: Copy {
                fn from(_value: #base_type_ident) -> Self {
                    #rep_new_base_call.expect(concat!("Failed to create Rep for ", #name_lit))
                }
            }
        };
        let base_try_from_impl = quote! {
            impl TryFrom<::spenso::structure::representation::LibraryRep> for #base_type_ident where #base_type_ident: Default {
                type Error = ::spenso::structure::representation::RepresentationError;

                fn try_from(rep: ::spenso::structure::representation::LibraryRep) -> ::std::result::Result<Self, Self::Error> {
                    let expected_rep = #rep_new_base_call.expect(concat!("Failed to create Rep for ", #name_lit));
                    if rep == expected_rep {
                        ::std::result::Result::Ok(#base_type_ident::default())
                    } else {
                        ::std::result::Result::Err(::spenso::structure::representation::RepresentationError::WrongRepresentationError(#name_lit.to_owned(), rep.to_string()))
                    }
                }
            }
        };

        let base_repname_impl = quote! {
            impl ::spenso::structure::representation::RepName for #base_type_ident where #base_type_ident: #base_bounds, #dual_type_ident: #dual_bounds {
                type Base = #base_type_ident;
                type Dual = #dual_type_ident;


                #[inline]
                fn orientation(self) -> ::linnet::half_edge::involution::Orientation {
                    ::linnet::half_edge::involution::Orientation::Default
                }

                #base_repname_common_impl
                #[inline]
                fn is_dual(self) -> bool { false }
                #[inline]
                fn matches(&self, _other: &Self::Dual) -> bool { true }
                #[inline]
                fn dual(self) -> Self::Dual where Self::Dual: Default {
                    #dual_type_ident::default()
                }
            }
        };
        let base_impls = quote! {
            impl #base_type_ident {
                pub const NAME: &'static str = #name_lit;
            }
            #base_from_impl
            #base_try_from_impl
            #base_repname_impl
            #base_display_impl
        };

        // 5. Implementations for the generated DUAL type (Remains Correct)
        let dual_display_impl = quote! {
            impl ::std::fmt::Display for #dual_type_ident where #dual_type_ident: Copy + Into<::spenso::structure::representation::LibraryRep> {
                fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
                    write!(f, "{}", ::spenso::structure::representation::LibraryRep::from(*self))
                }
            }
        };
        let dual_from_impl = quote! {
            impl From<#dual_type_ident> for ::spenso::structure::representation::LibraryRep where #dual_type_ident: Copy {
                fn from(_value: #dual_type_ident) -> Self { #rep_new_dual_call }
            }
        };
        let dual_try_from_impl = quote! {
            impl TryFrom<::spenso::structure::representation::LibraryRep> for #dual_type_ident where #dual_type_ident: Default {
                type Error = ::spenso::structure::representation::RepresentationError; fn try_from(rep: ::spenso::structure::representation::LibraryRep) -> ::std::result::Result<Self, Self::Error> {
                    let base_rep = #rep_new_base_call.expect(concat!("Failed to create dual Rep for ", #name_lit));
                    let expected_rep = base_rep.dual();
                    if rep == expected_rep {
                        ::std::result::Result::Ok(#dual_type_ident::default())
                    } else {
                        ::std::result::Result::Err(::spenso::structure::representation::RepresentationError::WrongRepresentationError(expected_rep.to_string(), rep.to_string()))
                    }
                }
            }
        };
        let dual_repname_impl = quote! {
            impl ::spenso::structure::representation::RepName for #dual_type_ident where #dual_type_ident: #dual_bounds, #base_type_ident: #base_bounds {
                type Base = #base_type_ident;
                type Dual = #base_type_ident;

                #[inline]
                fn orientation(self) -> ::linnet::half_edge::involution::Orientation {
                    ::linnet::half_edge::involution::Orientation::Reversed
                }
                #base_repname_common_impl
                #[inline]
                fn dual(self) -> Self::Dual where Self::Dual: Default { #base_type_ident::default() }
                #[inline]
                fn is_dual(self) -> bool { true }
                #[inline]
                fn matches(&self, _other: &Self::Dual) -> bool { true }
                #[inline]
                fn is_neg(self, i: usize) -> bool where Self: Copy, Self::Dual: Copy + ::spenso::structure::representation::RepName {
                    self.dual().is_neg(i)
                }
            }
        };
        let dual_impls = quote! {
            #dual_from_impl
            #dual_try_from_impl
            #dual_repname_impl
            #dual_display_impl
        };

        // 6. Combine generated code
        quote! {
            #dual_struct_def
            #base_impls
            #dual_impls
        }
    };

    TokenStream::from(expanded)
}
