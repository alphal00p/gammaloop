use spenso::network::tags::SPENSO_TAG;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, NamespacedSymbol, Symbol, UserData, representation::FunView},
    coefficient::CoefficientView,
    id::MatchSettings,
};

use crate::rep_symbols::RS;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy)]
pub enum CookingError {
    Add,
    Mul,
    Pow,
    RatCoeff,
    FiniteField,
    Float,
    Symbol,
}

impl CookingError {
    pub fn from_view(expr: AtomView) -> Self {
        match expr {
            AtomView::Var(_) => CookingError::Symbol,
            AtomView::Add(_) => CookingError::Add,
            AtomView::Mul(_) => CookingError::Mul,
            AtomView::Pow(_) => CookingError::Pow,
            AtomView::Num(_) => CookingError::RatCoeff,
            _ => CookingError::Symbol,
        }
    }
}

/// Selects how a cooked function payload is represented as a symbol.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CookMode {
    /// Build a readable symbol name from the function name and arguments.
    FlattenedSymbol,
    /// Build a stable symbol carrying the original atom in `UserData::Atom`.
    ReversibleEncoding,
}

/// Settings for cooking functions and representation-slot indices.
///
/// Tags are attached through `Symbol::new(...).with_tags(...)`, so they must be
/// fully namespaced Symbolica tags such as `idenso::cooked`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CookSettings {
    mode: CookMode,
    tags: Vec<String>,
}

impl Default for CookSettings {
    fn default() -> Self {
        Self::flattened()
    }
}

impl CookSettings {
    /// Use readable flattened names such as `f_0` or `f_g_mu`.
    pub fn flattened() -> Self {
        Self {
            mode: CookMode::FlattenedSymbol,
            tags: vec![],
        }
    }

    pub fn indices() -> Self {
        Self {
            mode: CookMode::FlattenedSymbol,
            tags: vec![SPENSO_TAG.index.clone()],
        }
    }

    /// Use a hash-stable symbol that can be restored with matching settings.
    pub fn reversible() -> Self {
        Self {
            mode: CookMode::ReversibleEncoding,
            tags: vec!["idenso::cooked".to_string()],
        }
    }

    /// Switch between flattened and reversible cooking.
    pub fn with_mode(mut self, mode: CookMode) -> Self {
        self.mode = mode;
        self
    }

    /// Replace the tags used for any newly created cooked symbols.
    pub fn with_tags<I, S>(mut self, tags: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.tags = tags
            .into_iter()
            .map(|tag| tag.as_ref().to_string())
            .collect();
        self
    }

    /// Current cooking mode.
    pub fn mode(&self) -> CookMode {
        self.mode
    }

    /// Tags attached to newly created cooked symbols.
    pub fn tags(&self) -> &[String] {
        &self.tags
    }

    /// Cook every function-like subexpression according to these settings.
    pub fn cook(&self, view: AtomView<'_>) -> Atom {
        view.replace_map(|a, _, out| {
            if let AtomView::Fun(f) = a
                && let Ok(cooked) = self.cook_function(f)
            {
                **out = Atom::var(cooked);
            }
        })
    }

    /// Restore reversible cooked symbols matching these settings.
    pub fn uncook(&self, view: AtomView<'_>) -> Atom {
        view.replace_map(|a, _, out| {
            if let AtomView::Var(s) = a {
                let symbol = s.get_symbol();
                if self.should_uncook_symbol(symbol)
                    && let UserData::Atom(a) = symbol.get_data()
                {
                    **out = a.clone();
                }
            }
        })
    }

    /// Cook only the abstract-index payloads of recognized representation slots.
    pub fn cook_indices(&self, view: AtomView<'_>) -> Atom {
        let mut expr = view.to_owned();

        let ipat = SPENSO_TAG.rep_::<0, _>([RS.d_, RS.a_]).to_pattern();
        expr = expr
            .replace(ipat)
            .when(RS.a_.filter_single(|a| a.as_fun_view().is_some()))
            .with_map(move |m| {
                if let Ok(aind) =
                    self.cook_function(m.get(RS.a_).unwrap().to_atom().as_fun_view().unwrap())
                {
                    SPENSO_TAG
                        .rep_::<0, _>([RS.d_, aind])
                        .to_pattern()
                        .replace_wildcards_with_matches(m)
                } else {
                    SPENSO_TAG
                        .rep_::<0, _>([RS.d_, RS.a_])
                        .to_pattern()
                        .replace_wildcards_with_matches(m)
                }
            });

        expr
    }

    /// Cook one atom into a symbol, preserving vars and numbers as-is.
    pub fn cook_function(&self, fun: FunView<'_>) -> Result<Symbol, CookingError> {
        match self.mode {
            CookMode::FlattenedSymbol => {
                let name = CookFunctionName::from_fun(fun)?;
                self.build_symbol(name, None)
            }
            CookMode::ReversibleEncoding => self.build_symbol(
                self.reversible_name(fun.as_view()),
                Some(fun.as_view().to_owned()),
            ),
        }
    }

    fn build_symbol(&self, name: String, source: Option<Atom>) -> Result<Symbol, CookingError> {
        let name = Self::namespaced_name(name);
        let mut builder = Symbol::new(NamespacedSymbol::parse(&name));
        if !self.tags.is_empty() {
            builder = builder.with_tags(&self.tags);
        }
        if let Some(source) = source {
            builder = builder.with_user_data(UserData::Atom(source));
        }
        builder.build().map_err(|_| CookingError::Symbol)
    }

    fn reversible_name(&self, source: AtomView<'_>) -> String {
        let mut hasher = blake3::Hasher::new();
        hasher.update(source.get_data());
        for tag in &self.tags {
            hasher.update(tag.as_bytes());
            hasher.update(&[0]);
        }
        let hash = hasher.finalize();
        let hex = hash.to_hex();
        format!("idenso::cooked_{}", &hex[..24])
    }

    fn namespaced_name(name: String) -> String {
        if name.contains("::") {
            name
        } else {
            format!("idenso::{name}")
        }
    }

    fn should_uncook_symbol(&self, symbol: Symbol) -> bool {
        self.mode == CookMode::ReversibleEncoding && self.tags.iter().all(|tag| symbol.has_tag(tag))
    }
}

/// Convenience methods for cooking Symbolica atoms.
pub trait Cookable {
    /// Reversibly cook function-like subexpressions with the default cooked tag.
    fn cook(&self) -> Atom;
    /// Cook function-like subexpressions with explicit settings.
    fn cook_with_settings(&self, settings: &CookSettings) -> Atom;
    /// Undo default reversible cooking.
    fn uncook(&self) -> Atom;
    /// Undo reversible cooking that used explicit settings.
    fn uncook_with_settings(&self, settings: &CookSettings) -> Atom;
    /// Flatten cook abstract-index payloads in recognized representation slots.
    fn cook_indices(&self) -> Atom;
    /// Cook representation-slot index payloads with explicit settings.
    fn cook_indices_with_settings(&self, settings: &CookSettings) -> Atom;
    /// Flatten cook a single atom into one symbol.
    fn cook_function(&self) -> Result<Atom, CookingError>;
    /// Cook a single atom into one symbol with explicit settings.
    fn cook_function_with_settings(&self, settings: &CookSettings) -> Result<Atom, CookingError>;
}

impl Cookable for Atom {
    fn cook(&self) -> Atom {
        self.as_view().cook()
    }

    fn cook_with_settings(&self, settings: &CookSettings) -> Atom {
        self.as_view().cook_with_settings(settings)
    }

    fn uncook(&self) -> Atom {
        self.as_view().uncook()
    }

    fn uncook_with_settings(&self, settings: &CookSettings) -> Atom {
        self.as_view().uncook_with_settings(settings)
    }

    fn cook_indices(&self) -> Atom {
        self.as_view().cook_indices()
    }

    fn cook_indices_with_settings(&self, settings: &CookSettings) -> Atom {
        self.as_view().cook_indices_with_settings(settings)
    }

    fn cook_function(&self) -> Result<Atom, CookingError> {
        self.as_view().cook_function()
    }

    fn cook_function_with_settings(&self, settings: &CookSettings) -> Result<Atom, CookingError> {
        self.as_view().cook_function_with_settings(settings)
    }
}

impl Cookable for AtomView<'_> {
    fn cook(&self) -> Atom {
        self.cook_with_settings(&CookSettings::reversible())
    }

    fn cook_with_settings(&self, settings: &CookSettings) -> Atom {
        settings.cook(*self)
    }

    fn uncook(&self) -> Atom {
        self.uncook_with_settings(&CookSettings::reversible())
    }

    fn uncook_with_settings(&self, settings: &CookSettings) -> Atom {
        settings.uncook(*self)
    }

    fn cook_indices(&self) -> Atom {
        self.cook_indices_with_settings(&CookSettings::indices())
    }

    fn cook_indices_with_settings(&self, settings: &CookSettings) -> Atom {
        settings.cook_indices(*self)
    }

    fn cook_function(&self) -> Result<Atom, CookingError> {
        self.cook_function_with_settings(&CookSettings::flattened())
    }

    fn cook_function_with_settings(&self, settings: &CookSettings) -> Result<Atom, CookingError> {
        settings
            .cook_function(self.as_fun_view().ok_or(CookingError::from_view(*self))?)
            .map(Atom::var)
    }
}

struct CookFunctionName;

impl CookFunctionName {
    fn from_fun(fun: FunView<'_>) -> Result<String, CookingError> {
        let mut name = fun.get_symbol().get_name().to_string();

        for arg in fun.iter() {
            name.push('_');
            Self::push_arg_name(&mut name, arg)?;
        }

        Ok(name)
    }

    fn push_arg_name(name: &mut String, arg: AtomView<'_>) -> Result<(), CookingError> {
        match arg {
            AtomView::Fun(f) => {
                let arg_name = Self::from_fun(f)?;
                name.push_str(Self::stripped_name(&arg_name));
            }
            AtomView::Num(n) => Self::push_number_name(name, n.get_coeff_view())?,
            AtomView::Var(s) => {
                name.push_str(s.get_symbol().get_stripped_name());
            }
            AtomView::Pow(_) => return Err(CookingError::Pow),
            AtomView::Add(_) => return Err(CookingError::Add),
            AtomView::Mul(_) => return Err(CookingError::Mul),
        }

        Ok(())
    }

    fn push_number_name(name: &mut String, coeff: CoefficientView<'_>) -> Result<(), CookingError> {
        match coeff {
            CoefficientView::Indeterminate => name.push_str("ind"),
            CoefficientView::Infinity(_) => name.push_str("inf"),
            CoefficientView::FiniteField(_, _) => return Err(CookingError::FiniteField),
            CoefficientView::Natural(n, d, imnum, imden) => {
                name.push_str(&n.to_string());
                if d != 1 {
                    name.push(':');
                    name.push_str(&d.to_string());
                }
                if imnum != 0 {
                    name.push('i');
                    name.push_str(&imnum.to_string());
                    if d != 1 {
                        name.push(':');
                        name.push_str(&imden.to_string());
                    }
                }
            }
            CoefficientView::Float(_, _) => return Err(CookingError::Float),
            CoefficientView::Large(r, imr) => {
                let rat = r.to_rat();
                name.push_str(&rat.numerator().to_string());
                if !rat.is_integer() {
                    name.push(':');
                    name.push_str(&rat.denominator().to_string());
                }

                if !imr.is_zero() {
                    let rat = imr.to_rat();
                    name.push('i');
                    name.push_str(&rat.numerator().to_string());
                    if !rat.is_integer() {
                        name.push(':');
                        name.push_str(&rat.denominator().to_string());
                    }
                }
            }
            CoefficientView::RationalPolynomial(_) => return Err(CookingError::RatCoeff),
        }

        Ok(())
    }

    fn stripped_name(name: &str) -> &str {
        name.rsplit_once("::")
            .map_or(name, |(_, stripped)| stripped)
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{
        atom::{Atom, AtomView, Symbol, UserData},
        parse_lit,
    };

    use crate::{Cookable, representations::initialize};

    use super::CookSettings;

    struct CookedIndex;

    impl CookedIndex {
        fn first_index_symbol(expr: &Atom) -> Symbol {
            let AtomView::Fun(tensor) = expr.as_view() else {
                panic!("expected tensor function");
            };
            let Some(slot) = tensor.iter().next() else {
                panic!("expected tensor argument");
            };
            let AtomView::Fun(slot) = slot else {
                panic!("expected representation slot");
            };
            let Some(index) = slot.iter().nth(1) else {
                panic!("expected representation index");
            };
            let AtomView::Var(index) = index else {
                panic!("expected cooked index symbol");
            };
            index.get_symbol()
        }
    }

    #[test]
    fn flattened_index_cooking_uses_custom_tags() {
        initialize();
        let expr = parse_lit!(p(spenso::mink(4, f(0))));
        let settings = CookSettings::flattened().with_tags(["idenso::flat_index"]);

        let cooked = expr.cook_indices_with_settings(&settings);
        let symbol = CookedIndex::first_index_symbol(&cooked);

        assert_eq!(symbol.get_stripped_name(), "f_0");
        assert!(symbol.has_tag("idenso::flat_index"));
    }

    #[test]
    fn reversible_index_cooking_uses_custom_tags_and_uncooks_with_settings() {
        initialize();
        let expr = parse_lit!(p(spenso::mink(4, f(0))));
        let settings = CookSettings::reversible().with_tags(["idenso::reversible_index"]);

        let cooked = expr.cook_indices_with_settings(&settings);
        let symbol = CookedIndex::first_index_symbol(&cooked);

        assert!(symbol.has_tag("idenso::reversible_index"));
        assert_eq!(symbol.get_data(), &UserData::Atom(parse_lit!(f(0))));
        assert_eq!(cooked.uncook_with_settings(&settings), expr);
    }

    #[test]
    fn reversible_expression_cooking_round_trips_with_settings() {
        let expr = parse_lit!(f(g(1)));
        let settings = CookSettings::reversible().with_tags(["idenso::reversible_expr"]);

        let cooked = expr.cook_with_settings(&settings);

        assert_eq!(cooked.uncook_with_settings(&settings), expr);
    }
}
