use symbolica::{
    atom::{Atom, AtomCore},
    printer::CanonicalOrderingSettings,
    try_parse,
};

use color_eyre::Result;
use eyre::eyre;

pub trait ToQuoted {
    fn to_quoted(&self) -> String;
}

// impl ToQuoted for SmartString<LazyCompact> {
//     fn to_quoted(&self) -> String {
//         format!("\"{}\"", self)
//     }
// }

pub trait ToOrderedSimple {
    fn to_ordered_simple(&self) -> String;
}

impl<A> ToOrderedSimple for A
where
    A: AtomCore,
{
    fn to_ordered_simple(&self) -> String {
        self.to_canonically_ordered_string(CanonicalOrderingSettings {
            include_namespace: false,
            include_attributes: false,
            hide_namespace: None,
        })
    }
}
impl<A> ToQuoted for A
where
    A: AtomCore,
{
    fn to_quoted(&self) -> String {
        // let mut opts = PrintOptions::file();
        // opts.hide_namespace = Some("gammalooprs");
        self.to_canonically_ordered_string(CanonicalOrderingSettings {
                include_namespace: true,
                include_attributes: false,
                hide_namespace: Some("gammalooprs")
            }).to_string()
        //.printer(opts))
    }
}

pub trait FromStripedStr: Sized {
    fn strip_from(string: &str) -> Result<Self>;
}

pub trait StripParse {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T>;
}

impl StripParse for String {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T> {
        T::strip_from(self.as_str())
    }
}
impl StripParse for &String {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T> {
        T::strip_from(self.as_str())
    }
}
impl StripParse for &str {
    fn strip_parse<T: FromStripedStr>(&self) -> Result<T> {
        T::strip_from(self)
    }
}
impl FromStripedStr for String {
    fn strip_from(string: &str) -> Result<Self> {
        Ok(string
            .strip_prefix('"')
            .unwrap_or(string)
            .strip_suffix('"')
            .unwrap_or(string)
            .into())
    }
}

impl FromStripedStr for Atom {
    fn strip_from(string: &str) -> Result<Self> {
        let a = string
            .strip_prefix('"')
            .unwrap_or(string)
            .strip_suffix('"')
            .unwrap_or(string);
        try_parse!(a).map_err(|e| eyre!("Symbolica parsing error: {e}"))
    }
}

impl FromStripedStr for bool {
    fn strip_from(s: &str) -> Result<Self> {
        Ok(s.strip_prefix('"')
            .unwrap_or(s)
            .strip_suffix('"')
            .unwrap_or(s)
            .parse()?)
    }
}
impl FromStripedStr for i32 {
    fn strip_from(s: &str) -> Result<Self> {
        Ok(s.strip_prefix('"')
            .unwrap_or(s)
            .strip_suffix('"')
            .unwrap_or(s)
            .parse()?)
    }
}
impl FromStripedStr for usize {
    fn strip_from(s: &str) -> Result<Self> {
        Ok(s.strip_prefix('"')
            .unwrap_or(s)
            .strip_suffix('"')
            .unwrap_or(s)
            .parse()?)
    }
}
