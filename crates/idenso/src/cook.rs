use crate::tensor::inference::InterfaceInference;
use ahash::HashMap;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex, AbstractIndexError};
use std::sync::{Arc, Mutex};
use symbolica::atom::representation::FunView;
use symbolica::atom::{
    Atom, AtomCore, AtomView, FunctionBuilder, NamespacedSymbol, Symbol, SymbolBuilder, UserData,
};
use symbolica::coefficient::CoefficientView;
use symbolica::parser::{ParseSettings, Token};
use thiserror::Error;

/// Failures when encoding symbolic function, dimension, or index payloads.
///
/// See the [native-index example](CookSettings::indices).
#[derive(Debug, Clone, PartialEq, Eq, Hash, Error)]
pub enum CookingError {
    #[error("cannot flatten a sum into a cooked symbol")]
    Add,
    #[error("cannot flatten a product into a cooked symbol")]
    Mul,
    #[error("cannot flatten a power into a cooked symbol")]
    Pow,
    #[error("cannot flatten this rational coefficient into a cooked symbol")]
    RatCoeff,
    #[error("cannot flatten a finite-field coefficient into a cooked symbol")]
    FiniteField,
    #[error("cannot flatten a floating-point coefficient into a cooked symbol")]
    Float,
    #[error("expected a function to cook")]
    Symbol,
    #[error("cannot create cooked symbol: {0}")]
    Symbolica(String),
    #[error("cooking dimensions requires reversible encoding")]
    DimensionRequiresReversibleEncoding,
    #[error("cannot cook a dimension with a normalization callback")]
    DimensionCallback,
    #[error("unsupported dimension payload for reversible cooking")]
    UnsupportedDimensionPayload,
    #[error("cyclic reversible index payload")]
    CyclicIndexPayload,
    #[error("invalid native index: {0}")]
    InvalidIndex(#[from] AbstractIndexError),
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

#[derive(Debug, Clone)]
struct ArcMutexOption<E> {
    value: Arc<Mutex<Option<E>>>,
}

impl<E> ArcMutexOption<E> {
    fn empty() -> Self {
        Self {
            value: Arc::new(Mutex::new(None)),
        }
    }

    fn set_once(&self, error: E) {
        let mut value = self
            .value
            .lock()
            .expect("replace_map error capture mutex poisoned");
        if value.is_none() {
            *value = Some(error);
        }
    }

    fn into_result<T>(self, value: T) -> Result<T, E> {
        match self
            .value
            .lock()
            .expect("replace_map error capture mutex poisoned")
            .take()
        {
            Some(error) => Err(error),
            None => Ok(value),
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

/// Selects which function occurrences are cooked.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CookSourceFilter {
    /// Cook every function-like subexpression.
    AnyFunction,
    /// Cook only functions whose head symbol passes this tag filter.
    FunctionTags(CookTagFilter),
    /// Cook selected payloads of representation declarations and explicit slots.
    /// The tag filter selects index function heads; dimension selection is explicit.
    /// Scalar-tagged function containers remain opaque metadata.
    RepresentationPayload {
        indices: bool,
        dimensions: bool,
        filter: Option<CookTagFilter>,
    },
}

/// Whether the function itself is an index identity or an ordinary payload.
#[derive(Clone, Copy, Debug)]
enum CookContext {
    Function,
    Index,
}

impl CookSourceFilter {
    fn regular_filter(&self) -> Option<Option<&CookTagFilter>> {
        match self {
            CookSourceFilter::AnyFunction => Some(None),
            CookSourceFilter::FunctionTags(filter) => Some(Some(filter)),
            CookSourceFilter::RepresentationPayload { .. } => None,
        }
    }

    fn index_payload_filter(&self) -> Option<&CookTagFilter> {
        match self {
            CookSourceFilter::AnyFunction => None,
            CookSourceFilter::FunctionTags(filter) => Some(filter),
            CookSourceFilter::RepresentationPayload { filter, .. } => filter.as_ref(),
        }
    }

    fn with_tag_filter(&mut self, filter: CookTagFilter) {
        match self {
            CookSourceFilter::RepresentationPayload {
                filter: index_filter,
                ..
            } => {
                *index_filter = Some(filter);
            }
            _ => {
                *self = CookSourceFilter::FunctionTags(filter);
            }
        }
    }
}

/// Tag predicate used to select input function heads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CookTagFilter {
    /// Match if any listed tag is present.
    Any(Vec<String>),
    /// Match only if all listed tags are present.
    All(Vec<String>),
    /// Match against the settings' explicit output tags.
    MatchedOutputTags,
}

impl CookTagFilter {
    pub fn any<I, S>(tags: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        Self::Any(Self::collect_tags(tags))
    }

    pub fn all<I, S>(tags: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        Self::All(Self::collect_tags(tags))
    }

    fn matches(&self, symbol: Symbol, output_tags: &[String]) -> bool {
        let tags = self.tags(output_tags);
        if tags.is_empty() {
            return false;
        }

        match self {
            CookTagFilter::Any(_) | CookTagFilter::MatchedOutputTags => {
                tags.iter().any(|tag| symbol.has_tag(tag))
            }
            CookTagFilter::All(_) => tags.iter().all(|tag| symbol.has_tag(tag)),
        }
    }

    fn matched_tags(&self, symbol: Symbol, output_tags: &[String]) -> Vec<String> {
        let tags = self.tags(output_tags);
        match self {
            CookTagFilter::All(_) if self.matches(symbol, output_tags) => tags.to_vec(),
            CookTagFilter::All(_) => vec![],
            CookTagFilter::Any(_) | CookTagFilter::MatchedOutputTags => tags
                .iter()
                .filter(|tag| symbol.has_tag(*tag))
                .cloned()
                .collect(),
        }
    }

    fn tags<'a>(&'a self, output_tags: &'a [String]) -> &'a [String] {
        match self {
            CookTagFilter::Any(tags) | CookTagFilter::All(tags) => tags,
            CookTagFilter::MatchedOutputTags => output_tags,
        }
    }

    fn collect_tags<I, S>(tags: I) -> Vec<String>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        tags.into_iter()
            .map(|tag| tag.as_ref().to_string())
            .collect()
    }
}

/// Selects which tags are attached to cooked output symbols.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CookOutputTags {
    /// Attach exactly these tags.
    Explicit(Vec<String>),
    /// Attach the input tags that were used to select the function.
    PreserveTags,
    /// Attach explicit output tags and preserve selected input tags.
    ExplicitAndPreserve(Vec<String>),
}

impl CookOutputTags {
    fn explicit<I, S>(tags: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        Self::Explicit(CookTagFilter::collect_tags(tags))
    }

    fn preserving(self) -> Self {
        match self {
            CookOutputTags::Explicit(tags) if tags.is_empty() => CookOutputTags::PreserveTags,
            CookOutputTags::Explicit(tags) | CookOutputTags::ExplicitAndPreserve(tags) => {
                CookOutputTags::ExplicitAndPreserve(tags)
            }
            CookOutputTags::PreserveTags => CookOutputTags::PreserveTags,
        }
    }

    fn explicit_tags(&self) -> &[String] {
        match self {
            CookOutputTags::Explicit(tags) | CookOutputTags::ExplicitAndPreserve(tags) => tags,
            CookOutputTags::PreserveTags => &[],
        }
    }

    fn tags_for(&self, symbol: Symbol, filter: Option<&CookTagFilter>) -> Vec<String> {
        match self {
            CookOutputTags::Explicit(tags) => tags.clone(),
            CookOutputTags::PreserveTags => {
                Self::preserved_tags(symbol, filter, self.explicit_tags())
            }
            CookOutputTags::ExplicitAndPreserve(tags) => {
                let mut output = tags.clone();
                for tag in Self::preserved_tags(symbol, filter, tags) {
                    Self::push_unique(&mut output, tag);
                }
                output
            }
        }
    }

    fn preserved_tags(
        symbol: Symbol,
        filter: Option<&CookTagFilter>,
        output_tags: &[String],
    ) -> Vec<String> {
        match filter {
            Some(filter) => filter.matched_tags(symbol, output_tags),
            None => symbol.get_tags().to_vec(),
        }
    }

    fn push_unique(tags: &mut Vec<String>, tag: String) {
        if !tags.iter().any(|existing| existing == &tag) {
            tags.push(tag);
        }
    }
}

/// Settings for cooking functions and representation-slot indices.
///
/// Tags are attached through `SymbolBuilder::new(...).with_tags(...)`, so they must be
/// fully namespaced Symbolica tags such as `idenso::cooked`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CookSettings {
    mode: CookMode,
    source: CookSourceFilter,
    output_tags: CookOutputTags,
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
            source: CookSourceFilter::AnyFunction,
            output_tags: CookOutputTags::Explicit(vec![]),
        }
    }

    /// Intern unsupported function labels inside representation slots, preserving
    /// the index forms already understood by `AbstractIndex`.
    ///
    /// # Example
    ///
    /// ```
    /// use idenso::CookSettings;
    /// use idenso::CookingError;
    /// use spenso::structure::abstract_index::AbstractIndex;
    /// use spenso::structure::representation::{Euclidean, RepName};
    /// use spenso::structure::slot::IsAbstractSlot;
    ///
    /// let port = Euclidean {}
    ///     .new_rep(3)
    ///     .to_lib()
    ///     .slot::<AbstractIndex, _>(AbstractIndex::Normal(7))
    ///     .to_atom();
    /// assert_eq!(CookSettings::indices().try_cook_indices(port.as_view())?, port);
    /// # Ok::<(), CookingError>(())
    /// ```
    pub fn indices() -> Self {
        Self {
            mode: CookMode::FlattenedSymbol,
            source: CookSourceFilter::RepresentationPayload {
                indices: true,
                dimensions: false,
                filter: None,
            },
            output_tags: CookOutputTags::PreserveTags,
        }
    }

    /// Use a hash-stable symbol that can be restored with matching settings.
    pub fn reversible() -> Self {
        Self {
            mode: CookMode::ReversibleEncoding,
            source: CookSourceFilter::AnyFunction,
            output_tags: CookOutputTags::Explicit(vec!["idenso::cooked".to_string()]),
        }
    }

    /// Switch between flattened and reversible cooking.
    pub fn with_mode(mut self, mode: CookMode) -> Self {
        self.mode = mode;
        self
    }

    /// Replace the tags used for any newly created cooked symbols.
    pub fn with_output_tags<I, S>(mut self, tags: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.output_tags = CookOutputTags::explicit(tags);
        self
    }

    /// Replace the source filter.
    pub fn with_source_filter(mut self, source: CookSourceFilter) -> Self {
        self.source = source;
        self
    }

    /// Select input function heads that carry any listed tag.
    pub fn with_input_tags<I, S>(mut self, tags: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.source.with_tag_filter(CookTagFilter::any(tags));
        self
    }

    /// Select input function heads that carry all listed tags.
    pub fn with_all_input_tags<I, S>(mut self, tags: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.source.with_tag_filter(CookTagFilter::all(tags));
        self
    }

    /// Select input function heads by matching the explicit output tags.
    pub fn matched_filter(mut self) -> Self {
        self.source
            .with_tag_filter(CookTagFilter::MatchedOutputTags);
        self
    }

    /// Preserve matched input tags on cooked output symbols.
    pub fn preserve_tags(mut self) -> Self {
        self.output_tags = self.output_tags.preserving();
        self
    }

    /// Cook only representation index payloads, with an optional tag filter.
    pub fn with_index_payload_filter(mut self, filter: Option<CookTagFilter>) -> Self {
        self.source = CookSourceFilter::RepresentationPayload {
            indices: true,
            dimensions: false,
            filter,
        };
        self
    }

    /// Select representation payloads without changing the cooking mode or index tag filter.
    ///
    /// Dimension cooking leaves atomic dimensions unchanged and requires reversible mode
    /// for compound dimensions. It supports exact rational arithmetic in ordinary symbols;
    /// functions, callbacks, and symbol metadata are rejected rather than hidden from inference.
    pub fn with_representation_payloads(mut self, indices: bool, dimensions: bool) -> Self {
        let filter = self.source.index_payload_filter().cloned();
        self.source = CookSourceFilter::RepresentationPayload {
            indices,
            dimensions,
            filter,
        };
        self
    }

    /// Current cooking mode.
    pub fn mode(&self) -> CookMode {
        self.mode
    }

    /// Which function occurrences are cooked.
    pub fn source_filter(&self) -> &CookSourceFilter {
        &self.source
    }

    /// Tags attached to newly created cooked symbols.
    pub fn output_tags(&self) -> &[String] {
        self.output_tags.explicit_tags()
    }

    /// Whether selected input tags are preserved on cooked output symbols.
    pub fn preserves_tags(&self) -> bool {
        matches!(
            self.output_tags,
            CookOutputTags::PreserveTags | CookOutputTags::ExplicitAndPreserve(_)
        )
    }

    /// Cook function-like subexpressions according to these settings.
    pub fn cook(&self, view: AtomView<'_>) -> Atom {
        self.try_cook(view)
            .unwrap_or_else(|error| panic!("failed to cook expression: {error:?}"))
    }

    /// Try to cook function-like subexpressions according to these settings.
    pub fn try_cook(&self, view: AtomView<'_>) -> Result<Atom, CookingError> {
        let Some(filter) = self.source.regular_filter() else {
            let CookSourceFilter::RepresentationPayload {
                indices,
                dimensions,
                filter,
            } = &self.source
            else {
                unreachable!();
            };
            return self.cook_representation_payloads(view, *indices, *dimensions, filter.as_ref());
        };

        let error = ArcMutexOption::empty();
        let cooked = view.replace_map(|a, _, out| {
            let AtomView::Fun(f) = a else {
                return;
            };

            match self.cook_function_atom(f, filter, CookContext::Function) {
                Ok(Some(cooked)) => **out = cooked,
                Ok(None) => {}
                Err(e) => error.set_once(e),
            }
        });

        error.into_result(cooked)
    }

    /// Restore reversible cooked symbols matching these settings.
    pub fn uncook(&self, view: AtomView<'_>) -> Atom {
        view.replace_map(|a, _, out| {
            if matches!(&self.source, CookSourceFilter::RepresentationPayload { .. })
                && matches!(a, AtomView::Fun(function) if function.get_symbol().is_scalar())
            {
                **out = a.to_owned();
                return;
            }
            if let AtomView::Var(s) = a {
                let symbol = s.get_symbol();
                if self.mode == CookMode::ReversibleEncoding
                    && let Some(payload) = Self::dimension_payload(symbol)
                {
                    **out = payload;
                } else if self.should_uncook_symbol(symbol)
                    && let UserData::Atom(a) = symbol.get_data()
                {
                    **out = a.clone();
                }
            }
        })
    }

    /// Inspect one public index through reversible payloads without rebuilding it.
    /// Ordinary function-valued labels remain opaque explicit identities. Only
    /// native Open/scope heads and this owner's symbol payloads are followed;
    /// in particular, ordinary labels never construct a failed index parse.
    pub(crate) fn index_has_open_port(
        mut index: spenso::structure::abstract_index::AbstractIndex,
    ) -> Result<bool, CookingError> {
        use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};

        while let AbstractIndex::Scoped(scope) = index {
            index = scope.index();
        }
        let AbstractIndex::Symbol(symbol) = index else {
            return Ok(matches!(index, AbstractIndex::Open { .. }));
        };
        let mut symbol: Symbol = symbol.into();
        let settings = Self::indices().with_mode(CookMode::ReversibleEncoding);
        // Ordinary cooked labels need no allocation. A visited list is needed
        // only when a stored payload refers to another encoded symbol.
        let mut visited = Vec::new();
        loop {
            if !settings.should_uncook_symbol(symbol) {
                return Ok(false);
            }
            let UserData::Atom(payload) = symbol.get_data() else {
                return Ok(false);
            };
            let mut view = payload.as_view();
            loop {
                match view {
                    AtomView::Fun(function) if function.get_symbol() == AIND_SYMBOLS.openind => {
                        return Ok(true);
                    }
                    AtomView::Fun(function)
                        if function.get_symbol() == AbstractIndex::scope_symbol()
                            && function.get_nargs() == 2 =>
                    {
                        let mut arguments = function.iter();
                        if !matches!(arguments.next(), Some(AtomView::Var(_))) {
                            return Ok(false);
                        }
                        view = arguments.next().unwrap();
                    }
                    AtomView::Var(variable) => {
                        let next = variable.get_symbol();
                        if next == symbol || visited.contains(&next) {
                            return Err(CookingError::CyclicIndexPayload);
                        }
                        visited.push(symbol);
                        symbol = next;
                        break;
                    }
                    _ => return Ok(false),
                }
            }
        }
    }

    /// Cook only the abstract-index payloads of recognized representation slots.
    pub fn cook_indices(&self, view: AtomView<'_>) -> Atom {
        self.try_cook_indices(view)
            .unwrap_or_else(|error| panic!("failed to cook representation indices: {error:?}"))
    }

    /// Try to cook only the abstract-index payloads of recognized representation slots.
    pub fn try_cook_indices(&self, view: AtomView<'_>) -> Result<Atom, CookingError> {
        self.cook_representation_payloads(view, true, false, self.source.index_payload_filter())
    }

    /// Cook one function call into an atom, or return it unchanged if skipped.
    /// Index-payload settings preserve validated native index forms.
    /// See [`Self::indices`] for native-index preservation.
    ///
    /// # Errors
    ///
    /// Rejects malformed native index functions, unsupported payload shapes, or
    /// symbol creation failures.
    pub fn cook_function(&self, fun: FunView<'_>) -> Result<Atom, CookingError> {
        let context = match self.source {
            CookSourceFilter::RepresentationPayload { indices: true, .. } => CookContext::Index,
            _ => CookContext::Function,
        };
        self.cook_function_atom(fun, self.source.index_payload_filter(), context)
            .map(|atom| atom.unwrap_or_else(|| fun.as_view().to_owned()))
    }

    fn cook_representation_payloads(
        &self,
        view: AtomView<'_>,
        indices: bool,
        dimensions: bool,
        filter: Option<&CookTagFilter>,
    ) -> Result<Atom, CookingError> {
        let error = ArcMutexOption::empty();
        // A contracted index appears on multiple tensor ports. Reuse its symbol
        // within this expression instead of repeatedly encoding and registering it.
        let mut cooked_indices = HashMap::default();
        let cooked = view.replace_map(|a, _, out| {
            let AtomView::Fun(rep) = a else {
                return;
            };
            if rep.get_symbol().is_scalar() {
                **out = a.to_owned();
                return;
            }
            if !rep.get_symbol().has_tag(&SPENSO_TAG.representation)
                || !(1..=2).contains(&rep.get_nargs())
            {
                return;
            }
            let mut args = rep.iter();
            let dim = args.next().unwrap();
            let index = args.next();
            let cooked_dim = if dimensions {
                match self.cook_dimension(dim) {
                    Ok(value) => value,
                    Err(e) => {
                        error.set_once(e);
                        return;
                    }
                }
            } else {
                None
            };
            let cooked_index = if indices && let Some(AtomView::Fun(index)) = index {
                match cooked_indices
                    .entry(index.as_view().to_owned())
                    .or_insert_with(|| self.cook_function_atom(index, filter, CookContext::Index))
                    .clone()
                {
                    Ok(value) => value,
                    Err(e) => {
                        error.set_once(e);
                        return;
                    }
                }
            } else {
                None
            };
            if cooked_dim.is_some() || cooked_index.is_some() {
                let mut builder = FunctionBuilder::new(rep.get_symbol())
                    .add_arg(cooked_dim.as_ref().map_or(dim, Atom::as_view));
                if let Some(index) = index {
                    builder = builder.add_arg(cooked_index.as_ref().map_or(index, Atom::as_view));
                }
                **out = builder.finish();
            }
        });
        error.into_result(cooked)
    }

    fn cook_dimension(&self, dimension: AtomView<'_>) -> Result<Option<Atom>, CookingError> {
        if matches!(dimension, AtomView::Var(_) | AtomView::Num(_)) {
            return Ok(None);
        }
        if self.mode != CookMode::ReversibleEncoding {
            return Err(CookingError::DimensionRequiresReversibleEncoding);
        }
        if !InterfaceInference::normalization_is_intrinsic(dimension) {
            return Err(CookingError::DimensionCallback);
        }
        if !Self::portable_dimension(dimension) {
            return Err(CookingError::UnsupportedDimensionPayload);
        }
        // Unlike the existing hash/UserData function encoding, this identity survives
        // plain-text export and a fresh process with no previously registered payload.
        let text = dimension.to_canonical_string();
        let mut name = String::from("idenso::cooked_dimension_v1_");
        for byte in text.bytes() {
            use std::fmt::Write;
            write!(&mut name, "{byte:02x}").unwrap();
        }
        self.build_symbol(name, None, self.output_tags.explicit_tags())
            .map(Atom::var)
            .map(Some)
    }

    fn portable_dimension(value: AtomView<'_>) -> bool {
        let mut supported = true;
        value.visitor(&mut |node| {
            supported &= match node {
                AtomView::Var(var) => {
                    let symbol = var.get_symbol();
                    matches!(symbol.get_data(), UserData::None)
                        && symbol.get_tags().is_empty()
                        // Nc is the existing symbolic color parameter. Its registered
                        // numerical value is not evaluated by this algebraic boundary.
                        && (symbol.get_evaluation_info().is_none() || symbol == crate::color::CS.nc)
                }
                AtomView::Num(number) => match number.get_coeff_view() {
                    CoefficientView::Natural(_, _, imaginary, _) => imaginary == 0,
                    CoefficientView::Large(_, imaginary) => imaginary.is_zero(),
                    _ => false,
                },
                AtomView::Fun(_) => false,
                AtomView::Pow(power) => matches!(power.get_exp(), AtomView::Num(exponent)
                    if exponent.get_coeff_view().is_integer()),
                _ => true,
            };
            supported
        });
        supported
    }

    /// Read a portable dimension encoded by this owner, without a process-local registry.
    /// Returns `None` for ordinary symbols or invalid encodings. Consumers can
    /// compare this scalar dimension while retaining native tensor slot syntax.
    pub fn dimension_payload(symbol: Symbol) -> Option<Atom> {
        let hex = symbol
            .get_name()
            .strip_prefix("idenso::cooked_dimension_v1_")?;
        if !hex.len().is_multiple_of(2) {
            return None;
        }
        let bytes = hex
            .as_bytes()
            .as_chunks::<2>()
            .0
            .iter()
            .map(|pair| {
                let text = std::str::from_utf8(pair).ok()?;
                u8::from_str_radix(text, 16).ok()
            })
            .collect::<Option<Vec<_>>>()?;
        let text = String::from_utf8(bytes).ok()?;
        // Tokenize without constructing products. Reject function calls before any
        // normalization can invoke a callback on a forged or unsupported payload.
        let token =
            Token::parse(&text, ParseSettings::default().convert_mul_to_atom(false)).ok()?;
        let mut pending = vec![&token];
        while let Some(token) = pending.pop() {
            match token {
                Token::Number(_, false) | Token::ID(_) => {}
                Token::Op(_, _, _, children) => pending.extend(children),
                _ => return None,
            }
        }
        let atom = Atom::parse(&text, "idenso", ParseSettings::default()).ok()?;
        Self::portable_dimension(atom.as_view()).then_some(atom)
    }

    fn cook_function_atom(
        &self,
        fun: FunView<'_>,
        filter: Option<&CookTagFilter>,
        context: CookContext,
    ) -> Result<Option<Atom>, CookingError> {
        let explicit_output_tags = self.output_tags.explicit_tags();
        if filter.is_some_and(|filter| !filter.matches(fun.get_symbol(), explicit_output_tags)) {
            return Ok(None);
        }

        let head = fun.get_symbol();
        if matches!(context, CookContext::Index) {
            if head == AIND_SYMBOLS.openind || head == AbstractIndex::scope_symbol() {
                return Ok(Some(Atom::from(AbstractIndex::try_from(fun.as_view())?)));
            }
            if head.has_tag(&SPENSO_TAG.index)
                && let Ok(index) = AbstractIndex::try_from(fun.as_view())
            {
                return Ok(Some(Atom::from(index)));
            }
        }

        let output_tags = self.output_tags.tags_for(fun.get_symbol(), filter);
        match self.mode {
            CookMode::FlattenedSymbol => {
                let name = CookFunctionName::from_fun(fun)?;
                self.build_symbol(name, None, &output_tags)
                    .map(|symbol| Some(Atom::var(symbol)))
            }
            CookMode::ReversibleEncoding => self
                .build_symbol(
                    self.reversible_name(fun.as_view(), &output_tags),
                    Some(fun.as_view().to_owned()),
                    &output_tags,
                )
                .map(|symbol| Some(Atom::var(symbol))),
        }
    }

    fn build_symbol(
        &self,
        name: String,
        source: Option<Atom>,
        tags: &[String],
    ) -> Result<Symbol, CookingError> {
        let name = Self::namespaced_name(name);
        let mut builder = SymbolBuilder::new(NamespacedSymbol::parse(&name));
        if !tags.is_empty() {
            builder = builder.with_tags(tags);
        }
        if let Some(source) = source {
            builder = builder.with_user_data(UserData::Atom(source));
        }
        builder
            .build()
            .map_err(|error| CookingError::Symbolica(format!("{error:?}")))
    }

    fn reversible_name(&self, source: AtomView<'_>, tags: &[String]) -> String {
        let mut hasher = blake3::Hasher::new();
        hasher.update(source.get_data());
        for tag in tags {
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
        if self.mode != CookMode::ReversibleEncoding {
            return false;
        }

        let explicit_tags = self.output_tags.explicit_tags();
        if explicit_tags.is_empty() {
            symbol.get_stripped_name().starts_with("cooked_")
        } else {
            explicit_tags.iter().all(|tag| symbol.has_tag(tag))
        }
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
        settings.cook_function(self.as_fun_view().ok_or(CookingError::from_view(*self))?)
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
    use super::{CookMode, CookSettings, CookingError};
    use crate::tensor::SymbolicTensor;
    use crate::test_support::test_initialize;
    use crate::{Cookable, coad, test_support};
    use spenso::network::tags::SPENSO_TAG;
    use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
    use spenso::structure::partial::PartialStructure;
    use spenso::structure::representation::{Minkowski, RepName};
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::structure::slot::ParseableAind;
    use spenso::tensor_symbol;
    use spenso::{index_symbol, mink};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::{
        Atom, AtomCore, AtomView, FunctionBuilder, NamespacedSymbol, Symbol, SymbolBuilder,
        UserData,
    };
    use symbolica::{function, parse_lit, symbol};

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

    struct CookedFunction;

    impl CookedFunction {
        fn tagged_symbol(name: &str, tags: &[&str]) -> Symbol {
            SymbolBuilder::new(NamespacedSymbol::parse(name))
                .with_tags(tags)
                .build()
                .unwrap()
        }

        fn call(symbol: Symbol, arg: i64) -> Atom {
            FunctionBuilder::new(symbol)
                .add_arg(Atom::num(arg))
                .finish()
        }

        fn assert_var_with_tag(expr: &Atom, name: &str, tag: &str) {
            let AtomView::Var(var) = expr.as_view() else {
                panic!("expected cooked variable");
            };

            let symbol = var.get_symbol();
            assert_eq!(symbol.get_stripped_name(), name);
            assert!(symbol.has_tag(tag));
        }
    }

    #[test]
    fn flattened_index_cooking_uses_custom_tags() {
        test_initialize();
        let expr = parse_lit!(p(spenso::mink(4, f(0))));
        let settings = CookSettings::flattened().with_output_tags(["idenso::flat_index"]);

        let cooked = expr.cook_indices_with_settings(&settings);
        let symbol = CookedIndex::first_index_symbol(&cooked);

        assert_eq!(symbol.get_stripped_name(), "f_0");
        assert!(symbol.has_tag("idenso::flat_index"));
    }

    #[test]
    fn reversible_index_cooking_uses_custom_tags_and_uncooks_with_settings() {
        test_initialize();
        let expr = parse_lit!(p(spenso::mink(4, f(0))));
        let settings = CookSettings::reversible().with_output_tags(["idenso::reversible_index"]);

        let cooked = expr.cook_indices_with_settings(&settings);
        let symbol = CookedIndex::first_index_symbol(&cooked);

        assert!(symbol.has_tag("idenso::reversible_index"));
        assert_eq!(symbol.get_data(), &UserData::Atom(parse_lit!(f(0))));
        assert_eq!(cooked.uncook_with_settings(&settings), expr);
    }

    #[test]
    fn native_index_forms_survive_both_interning_modes() {
        test_support::test_initialize();
        let head = index_symbol!("idenso::native_cooking_index");
        let named = AbstractIndex::Named(head.into(), 17, 2);
        let indices = [
            Atom::num(17),
            Atom::num(-17),
            parse_lit!(3 / 5),
            Atom::var(symbol!("native_cooking_mu")),
            named.to_atom(),
            AbstractIndex::Open { owner: 17, axis: 2 }.to_atom(),
            named.scoped(symbol!("native_cooking_scope")).to_atom(),
        ];
        for mode in [CookMode::FlattenedSymbol, CookMode::ReversibleEncoding] {
            let settings = CookSettings::indices().with_mode(mode);
            for index in &indices {
                let slot = mink!(4, index.clone());
                assert_eq!(settings.try_cook(slot.as_view()).unwrap(), slot);
                assert_eq!(settings.try_cook_indices(slot.as_view()).unwrap(), slot);
                if let Some(function) = index.as_fun_view() {
                    assert_eq!(settings.cook_function(function).unwrap(), *index);
                }
            }
        }
    }

    #[test]
    fn interning_mixed_indices_preserves_native_labels_and_shared_unknown_identity() {
        test_support::test_initialize();
        let native = function!(index_symbol!("idenso::mixed_native_index"), 17, 2);
        let opaque = parse_lit!(mixed_opaque_index(17, 2));
        let tensor = symbol!("mixed_index_cooking_tensor");
        let source = function!(
            tensor,
            mink!(4, native.clone()),
            mink!(4, opaque.clone()),
            coad!(8, opaque.clone())
        );
        for mode in [CookMode::FlattenedSymbol, CookMode::ReversibleEncoding] {
            let settings = CookSettings::indices().with_mode(mode);
            let interned = settings
                .cook_function(opaque.as_fun_view().unwrap())
                .unwrap();
            assert!(matches!(interned.as_view(), AtomView::Var(_)));
            let expected = function!(
                tensor,
                mink!(4, native.clone()),
                mink!(4, interned.clone()),
                coad!(8, interned.clone())
            );
            let actual = settings.try_cook_indices(source.as_view()).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(settings.try_cook_indices(actual.as_view()).unwrap(), actual);
            if mode == CookMode::ReversibleEncoding {
                assert_eq!(settings.uncook(actual.as_view()), source);
            }

            // The separate all-function cooking operation retains its contract.
            let generic = CookSettings::flattened().with_mode(mode);
            assert!(matches!(
                generic
                    .cook_function(native.as_fun_view().unwrap())
                    .unwrap()
                    .as_view(),
                AtomView::Var(_)
            ));
            assert!(matches!(
                generic.try_cook(native.as_view()).unwrap().as_view(),
                AtomView::Var(_)
            ));
        }
    }

    #[test]
    fn interning_accepts_opaque_tagged_heads_and_preserves_native_calls() {
        test_support::test_initialize();
        let head = index_symbol!("idenso::malformed_native_index");
        // Register GammaLoop-style heads before cooking, as in a shared process.
        let dummy = index_symbol!("idenso::cooking_gammaloop::dummy");
        let hedge = index_symbol!("idenso::cooking_gammaloop::hedge");
        let opaque = [function!(head, 17, -1), function!(head, 17, 2, 3)];
        let native = [
            function!(head, 17),
            function!(head, 17, 2),
            function!(dummy, 1),
            function!(hedge, 3),
        ];
        for mode in [CookMode::FlattenedSymbol, CookMode::ReversibleEncoding] {
            let settings = CookSettings::indices().with_mode(mode);
            for index in &opaque {
                if mode == CookMode::FlattenedSymbol && *index == function!(head, 17, -1) {
                    // Base flattened naming contains a literal minus sign, which
                    // Symbolica rejects. Preserve that symbol-creation error.
                    assert!(matches!(
                        settings.cook_function(index.as_fun_view().unwrap()),
                        Err(CookingError::Symbolica(_))
                    ));
                    assert!(matches!(
                        settings.try_cook_indices(mink!(4, index.clone()).as_view()),
                        Err(CookingError::Symbolica(_))
                    ));
                    continue;
                }
                let cooked = settings
                    .cook_function(index.as_fun_view().unwrap())
                    .unwrap();
                assert!(matches!(cooked.as_view(), AtomView::Var(_)));
                assert_eq!(
                    settings
                        .try_cook_indices(mink!(4, index.clone()).as_view())
                        .unwrap(),
                    mink!(4, cooked.clone())
                );
                if mode == CookMode::ReversibleEncoding {
                    assert_eq!(settings.uncook(cooked.as_view()), *index);
                }
            }
            for native in &native {
                assert_eq!(
                    settings
                        .cook_function(native.as_fun_view().unwrap())
                        .unwrap(),
                    *native
                );
                let slot = mink!(4, native.clone());
                assert_eq!(settings.try_cook_indices(slot.as_view()).unwrap(), slot);
            }
        }
    }

    #[test]
    fn index_cooking_canonicalizes_local_zero_aliases() {
        test_support::test_initialize();
        let head = index_symbol!("idenso::cooking_local_zero");
        let canonical = function!(head, 17);
        let alias = function!(head, 17, 0);
        for mode in [CookMode::FlattenedSymbol, CookMode::ReversibleEncoding] {
            let settings = CookSettings::indices().with_mode(mode);
            assert_eq!(
                settings
                    .cook_function(alias.as_fun_view().unwrap())
                    .unwrap(),
                canonical
            );
            let mixed = mink!(4, canonical.clone()) * mink!(4, alias.clone());
            let expected = mink!(4, canonical.clone()).pow(2);
            assert_eq!(
                settings.try_cook_indices(mixed.as_view()).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn interning_rejects_malformed_spenso_owned_index_calls() {
        test_support::test_initialize();
        let invalid = [
            function!(AbstractIndex::scope_symbol(), 17, 2),
            function!(AIND_SYMBOLS.openind, 17, -1),
        ];
        for mode in [CookMode::FlattenedSymbol, CookMode::ReversibleEncoding] {
            let settings = CookSettings::indices().with_mode(mode);
            for index in &invalid {
                assert!(matches!(
                    settings.cook_function(index.as_fun_view().unwrap()),
                    Err(CookingError::InvalidIndex(_))
                ));
                assert!(matches!(
                    settings.try_cook_indices(mink!(4, index.clone()).as_view()),
                    Err(CookingError::InvalidIndex(_))
                ));
            }
        }
    }

    #[test]
    fn reversible_expression_cooking_round_trips_with_settings() {
        let expr = parse_lit!(f(g(1)));
        let settings = CookSettings::reversible().with_output_tags(["idenso::reversible_expr"]);

        let cooked = expr.cook_with_settings(&settings);

        assert_eq!(cooked.uncook_with_settings(&settings), expr);
    }

    #[test]
    fn reversible_dimension_cooking_restores_declarations_and_slots_only() {
        test_initialize();
        let settings = CookSettings::reversible().with_representation_payloads(false, true);
        let expression = parse_lit!(trace(
            spenso::coad(Nc ^ 2 - 1),
            spectator(spenso::coad(Nc ^ 2 - 1, a))
        )) * parse_lit!(f(g(1))).cook_with_settings(&settings);
        let cooked = settings
            .try_cook(expression.as_view())
            .expect("formal dimensions should cook");

        assert_ne!(cooked, expression);
        assert_eq!(settings.uncook(cooked.as_view()), expression);
    }

    #[test]
    fn function_cooking_filters_by_input_tags() {
        let cook_tag = "idenso::cook_me";
        let tagged = CookedFunction::tagged_symbol("idenso::tagged_filter_function", &[cook_tag]);
        let plain = symbol!("plain_function");
        let expr = CookedFunction::call(tagged, 0) + CookedFunction::call(plain, 1);
        let settings = CookSettings::flattened().with_input_tags([cook_tag]);

        let cooked = expr.cook_with_settings(&settings);
        let AtomView::Add(sum) = cooked.as_view() else {
            panic!("expected sum");
        };

        let mut cooked_tagged = false;
        let mut kept_plain = false;
        for term in sum.iter() {
            match term {
                AtomView::Var(var)
                    if var.get_symbol().get_stripped_name() == "tagged_filter_function_0" =>
                {
                    cooked_tagged = true;
                }
                AtomView::Fun(fun) if fun.get_symbol() == plain => {
                    kept_plain = true;
                }
                _ => {}
            }
        }

        assert!(cooked_tagged);
        assert!(kept_plain);
    }

    #[test]
    fn matched_filter_uses_output_tags_as_input_tags() {
        let shared_tag = "idenso::shared";
        let tagged =
            CookedFunction::tagged_symbol("idenso::tagged_matched_function", &[shared_tag]);
        let expr = CookedFunction::call(tagged, 0);
        let settings = CookSettings::flattened()
            .with_output_tags([shared_tag])
            .matched_filter();

        let cooked = expr.cook_with_settings(&settings);

        CookedFunction::assert_var_with_tag(&cooked, "tagged_matched_function_0", shared_tag);
    }

    #[test]
    fn preserve_tags_copies_matched_input_tags() {
        let input_tag = "idenso::preserve_me";
        let tagged =
            CookedFunction::tagged_symbol("idenso::tagged_preserved_function", &[input_tag]);
        let expr = CookedFunction::call(tagged, 0);
        let settings = CookSettings::flattened()
            .with_input_tags([input_tag])
            .preserve_tags();

        let cooked = expr.cook_with_settings(&settings);

        CookedFunction::assert_var_with_tag(&cooked, "tagged_preserved_function_0", input_tag);
    }

    #[test]
    fn index_payload_filter_cooks_only_matching_payloads() {
        test_initialize();
        let input_tag = "idenso::index_payload";
        let tagged = CookedFunction::tagged_symbol("idenso::tagged_index", &[input_tag]);
        let payload = CookedFunction::call(tagged, 0);
        let expr = FunctionBuilder::new(symbol!("p"))
            .add_arg(SPENSO_TAG.rep_::<0, _>([Atom::num(4), payload]))
            .finish();
        let settings = CookSettings::indices().with_input_tags([input_tag]);

        let cooked = expr.cook_with_settings(&settings);
        let symbol = CookedIndex::first_index_symbol(&cooked);

        assert_eq!(symbol.get_stripped_name(), "tagged_index_0");
        assert!(symbol.has_tag(input_tag));
    }

    #[test]
    fn index_payload_filter_leaves_unmatched_payloads() {
        test_initialize();
        let expr = parse_lit!(p(spenso::mink(4, f(0))));
        let settings = CookSettings::indices().with_input_tags(["idenso::index_payload"]);

        assert_eq!(expr.cook_with_settings(&settings), expr);
    }

    #[test]
    fn representation_payload_cooking_preserves_shared_index_identity() {
        use crate::CookMode;
        use symbolica::function;

        test_initialize();
        let index = parse_lit!(shared_payload(7));
        let color = crate::coad!(8, index.clone());
        let lorentz = spenso::mink!(4, index.clone());
        let source = function!(
            symbol!("shared_payload_tensor"),
            color.clone(),
            lorentz.clone()
        ) * function!(
            symbol!("shared_payload_spectator"),
            color.clone(),
            index.clone()
        );
        for mode in [CookMode::FlattenedSymbol, CookMode::ReversibleEncoding] {
            let settings = CookSettings::indices().with_mode(mode);
            let cooked_color = settings.cook(color.as_view());
            let cooked_lorentz = settings.cook(lorentz.as_view());
            let expected = function!(
                symbol!("shared_payload_tensor"),
                cooked_color.clone(),
                cooked_lorentz
            ) * function!(
                symbol!("shared_payload_spectator"),
                cooked_color,
                index.clone()
            );
            let cooked = settings.cook(source.as_view());
            assert_eq!(cooked, expected);
            assert_eq!(settings.cook(cooked.as_view()), cooked);
            if mode == CookMode::ReversibleEncoding {
                assert_eq!(settings.uncook(cooked.as_view()), source);
            }
        }
    }

    #[test]
    fn representation_dimensions_are_explicit_reversible_and_atomic_identity() {
        use crate::{CookMode, CookingError, tensor::SymbolicTensor};
        use symbolica::parse;
        test_initialize();
        let tensor = spenso::tensor_symbol!("dimension_cooking_tensor");
        let _ = symbol!("idenso::metadata"; Scalar);
        let source = FunctionBuilder::new(tensor)
            .add_arg(parse!("idenso::metadata(spenso::Nc^2-1)"))
            .add_arg(parse!("spenso::coad(spenso::Nc^2-1,idenso::leg(0))"))
            .finish();
        let indices = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let dimensions = indices.clone().with_representation_payloads(true, true);
        let cooked = dimensions.try_cook(source.as_view()).unwrap();
        assert!(SymbolicTensor::infer(indices.cook(source.as_view())).is_err());
        SymbolicTensor::infer(cooked.clone())
            .expect("cooked dimensions with declared scalar metadata");
        assert_eq!(dimensions.uncook(cooked.as_view()), source);
        assert_eq!(dimensions.cook(cooked.as_view()), cooked);
        assert!(cooked.to_string().contains("metadata"));
        let atomic = parse!("idenso::dimension_cooking_tensor(spenso::coad(8,a))");
        assert_eq!(dimensions.cook(atomic.as_view()), atomic);
        assert_eq!(CookSettings::default().mode(), CookMode::FlattenedSymbol);
        assert_eq!(CookSettings::indices().mode(), CookMode::FlattenedSymbol);
        assert_eq!(
            CookSettings::indices()
                .with_representation_payloads(true, true)
                .try_cook(source.as_view()),
            Err(CookingError::DimensionRequiresReversibleEncoding)
        );
    }

    #[test]
    fn dimension_cooking_refuses_callbacks_and_metadata_without_invoking_them() {
        use crate::{CookMode, CookingError};
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        use symbolica::atom::AtomCore;
        test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = calls.clone();
        let hook = symbol!(
            "dimension_cooking_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let settings = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(false, true);
        let callback = FunctionBuilder::new(hook).add_arg(1).finish();
        let source = FunctionBuilder::new(crate::color::CS.adjoint_rep)
            .add_arg(&callback)
            .finish();
        calls.store(0, Ordering::Relaxed);
        assert_eq!(
            settings.try_cook(source.as_view()),
            Err(CookingError::DimensionCallback)
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        let metadata = symbol!("dimension_metadata", data = UserData::Integer(17));
        let source = FunctionBuilder::new(crate::color::CS.adjoint_rep)
            .add_arg(Atom::var(metadata) + Atom::num(1))
            .finish();
        assert_eq!(
            settings.try_cook(source.as_view()),
            Err(CookingError::UnsupportedDimensionPayload)
        );
        let evaluated = symbol!(
            "dimension_evaluated_parameter",
            eval = symbolica::atom::EvaluationInfo::constant(|_tags, prec| {
                Ok(symbolica::domains::rational::Rational::new(7, 1)
                    .to_multi_prec_float(prec)
                    .into())
            })
        );
        let source = FunctionBuilder::new(crate::color::CS.adjoint_rep)
            .add_arg(Atom::var(evaluated) + Atom::num(1))
            .finish();
        assert_eq!(
            settings.try_cook(source.as_view()),
            Err(CookingError::UnsupportedDimensionPayload)
        );
        for exponent in [
            Atom::num(1) / Atom::num(2),
            Atom::var(symbol!("dimension_exponent")),
        ] {
            let dimension = Atom::var(symbol!("dimension_base")).pow(exponent);
            let source = FunctionBuilder::new(crate::color::CS.adjoint_rep)
                .add_arg(dimension)
                .finish();
            assert_eq!(
                settings.try_cook(source.as_view()),
                Err(CookingError::UnsupportedDimensionPayload)
            );
        }
        let function = parse_lit!(opaque(n));
        let source = FunctionBuilder::new(crate::color::CS.adjoint_rep)
            .add_arg(function)
            .finish();
        assert_eq!(
            settings.try_cook(source.as_view()),
            Err(CookingError::UnsupportedDimensionPayload)
        );
        // A forged payload must be rejected before a function can be normalized.
        let text = callback.to_canonical_string();
        let name = format!(
            "idenso::cooked_dimension_v1_{}",
            text.bytes().map(|x| format!("{x:02x}")).collect::<String>()
        );
        let forged = Atom::var(
            SymbolBuilder::new(NamespacedSymbol::parse(&name))
                .build()
                .unwrap(),
        );
        calls.store(0, Ordering::Relaxed);
        assert_eq!(settings.uncook(forged.as_view()), forged);
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn portable_dimension_plain_text_round_trips_in_a_fresh_process() {
        use crate::CookMode;
        use symbolica::parse;
        const KEY: &str = "IDENSO_DIMENSION_PLAIN_TEST";
        test_initialize();
        let settings = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(false, true);
        let expected = parse!("spenso::coad(spenso::Nc^2-1)");
        if let Ok(text) = std::env::var(KEY) {
            let parsed = Atom::parse(text, "idenso", Default::default()).unwrap();
            let AtomView::Fun(rep) = parsed.as_view() else {
                panic!("representation");
            };
            let AtomView::Var(dim) = rep.iter().next().unwrap() else {
                panic!("atomic dimension");
            };
            assert_eq!(dim.get_symbol().get_data(), &UserData::None);
            assert_eq!(settings.uncook(parsed.as_view()), expected);
            return;
        }
        let cooked = settings.cook(expected.as_view());
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "cook::tests::portable_dimension_plain_text_round_trips_in_a_fresh_process",
                "--nocapture",
            ])
            .env(KEY, cooked.as_view().to_plain_string())
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(String::from_utf8_lossy(&output.stdout).contains("1 passed"));
    }

    #[test]
    fn representation_payload_cooking_keeps_scalar_metadata_and_callback_input() {
        test_initialize();
        let scope = symbol!("cook_scalar_metadata_scope");
        let inner = Minkowski {}
            .new_rep(4)
            .to_symbolic([function!(symbol!("cook_scalar_private_label"), 1701)]);
        let outer = Minkowski {}
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(1702).scoped(scope))
            .to_atom();
        let metadata = function!(
            symbol!("cook_scalar_private_metadata"; Scalar),
            function!(tensor_symbol!("cook_scalar_private_a"), &inner)
                * function!(tensor_symbol!("cook_scalar_private_b"), &inner)
        );
        let bad_calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&bad_calls);
        let expected = metadata.clone();
        let tensor = tensor_symbol!(
            "cook_scalar_metadata_tensor",
            norm = move |value, output| {
                if !matches!(value, AtomView::Fun(function)
                    if function.iter().next() == Some(expected.as_view()))
                {
                    observed.fetch_add(1, Ordering::Relaxed);
                    output.to_num(0);
                }
            }
        );
        let component =
            SymbolicTensor::<PartialStructure>::infer(function!(tensor, &metadata, outer)).unwrap();
        let source = component
            .pow(&SymbolicTensor::infer(Atom::num(2)).unwrap())
            .unwrap();
        assert!(!source.expression().is_zero());
        assert_eq!(bad_calls.load(Ordering::Relaxed), 0);
        for mode in [CookMode::FlattenedSymbol, CookMode::ReversibleEncoding] {
            let settings = CookSettings::indices().with_mode(mode);
            assert_ne!(settings.try_cook(inner.as_view()).unwrap(), inner);
            let cooked = settings.try_cook(source.expression().as_view()).unwrap();
            assert!(!cooked.is_zero());
            assert_eq!(&cooked, source.expression());
            assert_eq!(bad_calls.load(Ordering::Relaxed), 0);
            assert_eq!(settings.cook(cooked.as_view()), cooked);
            if mode == CookMode::ReversibleEncoding {
                assert_eq!(settings.uncook(cooked.as_view()), *source.expression());
                assert_eq!(bad_calls.load(Ordering::Relaxed), 0);
            }
        }
    }

    #[test]
    fn representation_payload_uncooking_preserves_literal_encoded_scalar_metadata() {
        test_initialize();
        let label = function!(symbol!("cook_scalar_explicit_label"; Scalar), 5);
        let original = function!(
            tensor_symbol!("cook_scalar_literal_tensor"),
            Minkowski {}.new_rep(4).to_symbolic([label])
        );
        let settings = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let cooked = settings.cook(original.as_view());
        assert_ne!(cooked, original);
        assert_eq!(settings.uncook(cooked.as_view()), original);
        let index = Atom::var(CookedIndex::first_index_symbol(&cooked));
        let dimensions = settings.with_representation_payloads(true, true);
        let dimension_source = parse_lit!(spenso::coad(cook_scalar_rank ^ 2 - 1));
        let encoded_dimension = dimensions.cook(dimension_source.as_view());
        let AtomView::Fun(rep) = encoded_dimension.as_view() else {
            panic!("cooked representation");
        };
        let dimension = rep.iter().next().unwrap().to_owned();
        let metadata = function!(
            symbol!("cook_scalar_literal_metadata"; Scalar),
            index,
            dimension
        );
        let source = function!(
            tensor_symbol!("cook_scalar_literal_spectator"),
            &metadata,
            &dimension_source
        );
        let expected_cooked = function!(
            tensor_symbol!("cook_scalar_literal_spectator"),
            metadata,
            encoded_dimension
        );
        assert_eq!(dimensions.cook(source.as_view()), expected_cooked);
        assert_eq!(dimensions.uncook(expected_cooked.as_view()), source);
    }

    #[test]
    fn any_function_cooking_and_uncooking_still_enters_scalar_functions() {
        let head = symbol!("cook_scalar_any_function"; Scalar);
        let payload = parse_lit!(cook_scalar_any_payload(7));
        let settings = CookSettings::reversible();
        let original = function!(head, &payload);
        let cooked = settings.cook(original.as_view());
        assert!(matches!(cooked.as_view(), AtomView::Var(_)));
        assert_eq!(settings.uncook(cooked.as_view()), original);
        let encoded_payload = settings.cook(payload.as_view());
        let literal = function!(head, encoded_payload);
        assert_eq!(settings.uncook(literal.as_view()), original);
    }
}
