//! Shared process-card grammar. Runtime/build options remain with their phase owners.
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::LazyLock,
};

use feynkit_model::Model;
use regex::Regex;
use thiserror::Error;

use crate::{
    CutAmplitude, Cuts, GenerationType, Loops, Orders, ParticleSelector, Process, ProcessError,
    Range, Request, Restrict,
};

#[derive(Debug, Error, PartialEq, Eq)]
pub enum RequestParseError {
    #[error("expected 'to' or '>' between initial and final states")]
    MissingArrow,
    #[error("unbalanced brackets/braces")]
    Unbalanced,
    #[error("invalid token '{0}'")]
    InvalidToken(String),
    #[error("unknown particle '{name}'. Valid choices: {choices}")]
    UnknownParticle { name: String, choices: String },
    #[error("unknown coupling '{name}'. Valid choices: {choices}")]
    UnknownCoupling { name: String, choices: String },
    #[error("every generation request requires an explicit loop count")]
    MissingLoops,
    #[error("amplitudes require total loops; cross sections require forward-graph loops")]
    WrongLoopDomain,
    #[error("the process text and explicit loop argument specify different loop counts")]
    ConflictingLoops,
    #[error(transparent)]
    Process(#[from] ProcessError),
}

impl Request {
    /// Parse the same process grammar used by the CLI, including squared coupling orders.
    /// `loops` is an explicit override from a separate argument; it must agree with the text.
    /// An amplitude uses `[{n}]`; a cross section uses `[{{n}}]` for forward loops and
    /// optionally `[{m}]` for the sum of the left and right cut-amplitude loop counts.
    pub fn parse(
        model: &Model,
        kind: GenerationType,
        text: &str,
        loops: Option<Loops>,
    ) -> Result<Self, RequestParseError> {
        let parsed = Syntax::parse(text)?;
        let incoming = parsed
            .incoming
            .iter()
            .map(|name| Syntax::particle(model, name))
            .collect::<Result<Vec<_>, _>>()?;
        let finals = parsed
            .finals
            .iter()
            .map(|states| {
                states
                    .iter()
                    .map(|name| Syntax::particle(model, name))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        if kind == GenerationType::Amplitude && finals.len() != 1 {
            return Err(ProcessError::MultipleAmplitudeFinalStates.into());
        }
        let process = Process::new(incoming, finals.first().cloned().unwrap_or_default())
            .with_final_state_alternatives(finals)?;
        let text_loops = match kind {
            GenerationType::Amplitude => {
                if parsed.forward_loops.is_some() {
                    return Err(RequestParseError::WrongLoopDomain);
                }
                parsed.amplitude_loops.map(|n| Loops::total(n..=n))
            }
            GenerationType::CrossSection => parsed.forward_loops.map(|n| Loops::forward(n..=n)),
        };
        if loops.as_ref().is_some_and(|bound| {
            !matches!(
                (kind, bound),
                (GenerationType::Amplitude, Loops::Total(_))
                    | (GenerationType::CrossSection, Loops::Forward(_))
            )
        }) {
            return Err(RequestParseError::WrongLoopDomain);
        }
        if let (Some(text), Some(argument)) = (&text_loops, &loops) {
            if text != argument {
                return Err(RequestParseError::ConflictingLoops);
            }
        }
        let loops = loops
            .or(text_loops)
            .ok_or(RequestParseError::MissingLoops)?;
        if loops.range().start() > loops.range().end() {
            return Err(ProcessError::InvalidLoopRange {
                minimum: *loops.range().start(),
                maximum: *loops.range().end(),
            }
            .into());
        }
        let mut restrict = Restrict::new().exclude_particles(
            parsed
                .veto
                .iter()
                .map(|name| Syntax::particle(model, name))
                .collect::<Result<Vec<_>, _>>()?,
        );
        for only in parsed.only {
            restrict = restrict.only_particles(
                only.iter()
                    .map(|name| Syntax::particle(model, name))
                    .collect::<Result<Vec<_>, _>>()?,
            );
        }
        let mut orders = Orders::new();
        let mut cut_orders = Orders::new();
        for (name, squared, op, power) in parsed.orders {
            let name = Syntax::coupling(model, &name)?;
            let target = if kind == GenerationType::CrossSection && !squared {
                &mut cut_orders
            } else {
                &mut orders
            };
            *target = match (squared, op.as_str()) {
                (false, "==") => std::mem::take(target).exact(name, power),
                (false, ">=") => std::mem::take(target).min(name, power),
                (false, "<=") => std::mem::take(target).max(name, power),
                (true, "==") => std::mem::take(target).squared_exact(name, power),
                (true, ">=") => std::mem::take(target).squared_min(name, power),
                (true, "<=") => std::mem::take(target).squared_max(name, power),
                _ => unreachable!("the tokenizer checks order operators"),
            };
        }
        for (name, power) in parsed.perturbative {
            orders = orders.perturbative(Syntax::coupling(model, &name)?, power);
        }
        let request = match kind {
            GenerationType::Amplitude => Request::amplitude(process, loops),
            GenerationType::CrossSection => Request::cross_section(process, loops).with_cuts(
                Cuts::default().with_cut_amplitude(CutAmplitude {
                    orders: cut_orders,
                    loops: parsed.amplitude_loops.map_or(Range::Unrestricted, |n| {
                        Range::Inclusive {
                            min: n,
                            max: Some(n),
                        }
                    }),
                }),
            ),
        };
        Ok(request.with_restrict(restrict).with_orders(orders))
    }
}

struct Syntax {
    incoming: Vec<String>,
    finals: Vec<Vec<String>>,
    veto: BTreeSet<String>,
    only: Vec<BTreeSet<String>>,
    orders: Vec<(String, bool, String, usize)>,
    amplitude_loops: Option<usize>,
    forward_loops: Option<usize>,
    perturbative: BTreeMap<String, usize>,
}

impl Syntax {
    fn parse(text: &str) -> Result<Self, RequestParseError> {
        let mut braces = 0_i32;
        let mut brackets = 0_i32;
        for c in text.chars() {
            match c {
                '{' => braces += 1,
                '}' => braces -= 1,
                '[' => brackets += 1,
                ']' => brackets -= 1,
                _ => {}
            }
            if braces < 0 || brackets < 0 {
                return Err(RequestParseError::Unbalanced);
            }
        }
        if braces != 0 || brackets != 0 {
            return Err(RequestParseError::Unbalanced);
        }
        let (lhs, rhs) = Self::arrow(text).ok_or(RequestParseError::MissingArrow)?;
        let mut result = Self {
            incoming: if lhs.trim() == "{}" || lhs.trim().eq_ignore_ascii_case("empty-initial") {
                Vec::new()
            } else {
                split_ws(lhs)
            },
            finals: Vec::new(),
            veto: BTreeSet::new(),
            only: Vec::new(),
            orders: Vec::new(),
            amplitude_loops: None,
            forward_loops: None,
            perturbative: BTreeMap::new(),
        };
        let mut rhs = rhs.to_owned();
        if let Some(start) = rhs.find('[') {
            let end = rhs[start..]
                .find(']')
                .map(|n| n + start)
                .ok_or(RequestParseError::Unbalanced)?;
            for token in rhs[start + 1..end].split_whitespace() {
                let (slot, digits) = if let Some(n) =
                    token.strip_prefix("{{").and_then(|s| s.strip_suffix("}}"))
                {
                    (&mut result.forward_loops, Some(n))
                } else if let Some(n) = token.strip_prefix('{').and_then(|s| s.strip_suffix('}')) {
                    (&mut result.amplitude_loops, Some(n))
                } else {
                    (&mut result.amplitude_loops, None)
                };
                if let Some(digits) = digits {
                    let n = digits
                        .parse::<usize>()
                        .map_err(|_| RequestParseError::InvalidToken(token.into()))?;
                    if slot.is_some_and(|old| old != n) {
                        return Err(RequestParseError::ConflictingLoops);
                    }
                    *slot = Some(n);
                } else {
                    let (name, n) = token.split_once('=').map_or((token, "1"), |p| p);
                    if !NAME.is_match(name) {
                        return Err(RequestParseError::InvalidToken(token.into()));
                    }
                    let n = n
                        .parse::<usize>()
                        .map_err(|_| RequestParseError::InvalidToken(token.into()))?;
                    if result
                        .perturbative
                        .insert(name.into(), n)
                        .is_some_and(|old| old != n)
                    {
                        return Err(RequestParseError::InvalidToken(token.into()));
                    }
                }
            }
            rhs.replace_range(start..=end, "");
            if rhs.contains(['[', ']']) {
                return Err(RequestParseError::Unbalanced);
            }
        }
        let (finals, _, options) = parse_final_states(tokenize(rhs)?)?;
        result.finals = if finals.sets.is_empty() {
            vec![finals.primary]
        } else {
            finals.sets
        };
        let mut rest = options.iter().peekable();
        while let Some(token) = rest.next() {
            match token.as_str() {
                "/" | "veto" | "|" | "only" => {
                    let set = if matches!(token.as_str(), "/" | "veto") {
                        &mut result.veto
                    } else {
                        {
                            result.only.push(BTreeSet::new());
                            result.only.last_mut().unwrap()
                        }
                    };
                    while rest.peek().is_some_and(|s| !is_option_token(s)) {
                        set.insert(rest.next().unwrap().clone());
                    }
                }
                _ => {
                    let c = ORDER
                        .captures(token)
                        .ok_or_else(|| RequestParseError::InvalidToken(token.clone()))?;
                    let squared = match c.get(2).map(|m| m.as_str()) {
                        None | Some("1") => false,
                        Some("2") => true,
                        _ => return Err(RequestParseError::InvalidToken(token.clone())),
                    };
                    let power = c[4]
                        .parse::<usize>()
                        .map_err(|_| RequestParseError::InvalidToken(token.clone()))?;
                    result
                        .orders
                        .push((c[1].into(), squared, c[3].into(), power));
                }
            }
        }
        Ok(result)
    }

    fn arrow(text: &str) -> Option<(&str, &str)> {
        let mut depth = 0;
        for (i, c) in text.char_indices() {
            match c {
                '[' | '{' => depth += 1,
                ']' | '}' => depth -= 1,
                '>' if depth == 0 => return Some((&text[..i], &text[i + 1..])),
                't' | 'T'
                    if depth == 0
                        && text[i..]
                            .get(..2)
                            .is_some_and(|s| s.eq_ignore_ascii_case("to")) =>
                {
                    let before = text[..i].chars().next_back();
                    let after = text[i + 2..].chars().next();
                    if before.is_none_or(char::is_whitespace)
                        && after.is_none_or(char::is_whitespace)
                    {
                        return Some((&text[..i], &text[i + 2..]));
                    }
                }
                _ => {}
            }
        }
        None
    }

    fn particle(model: &Model, name: &str) -> Result<ParticleSelector, RequestParseError> {
        if let Ok(pdg) = name.parse::<i64>() {
            if model.particle_id_by_pdg(pdg).is_ok() {
                return Ok(ParticleSelector::Pdg(pdg));
            }
        } else if let Some(particle) =
            model
                .particles()
                .iter()
                .find(|p| p.name == name)
                .or_else(|| {
                    model
                        .particles()
                        .iter()
                        .find(|p| p.name.eq_ignore_ascii_case(name))
                })
        {
            return Ok(ParticleSelector::Name(particle.name.to_string()));
        }
        Err(RequestParseError::UnknownParticle {
            name: name.into(),
            choices: model
                .particles()
                .iter()
                .map(|p| p.name.to_string())
                .collect::<Vec<_>>()
                .join(", "),
        })
    }

    fn coupling(model: &Model, name: &str) -> Result<String, RequestParseError> {
        if model.orders().is_empty() {
            return Ok(name.into());
        }
        if let Some(order) = model
            .orders()
            .iter()
            .find(|o| o.name.eq_ignore_ascii_case(name))
        {
            return Ok(order.name.to_string());
        }
        Err(RequestParseError::UnknownCoupling {
            name: name.into(),
            choices: model
                .orders()
                .iter()
                .map(|o| o.name.to_string())
                .collect::<Vec<_>>()
                .join(", "),
        })
    }
}

static NAME: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^[A-Za-z_][A-Za-z0-9_]*$").unwrap());
static ORDER: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^([A-Za-z_][A-Za-z0-9_]*)(?:\^(\d+))?(==|>=|<=)(\d+)$").unwrap());

fn is_option_token(token: &str) -> bool {
    matches!(
        token,
        "/" | "|" | "veto" | "only" | "amp" | "xs" | "pert" | "sets" | "empty-final"
    ) || token.starts_with('[')
        || ORDER.is_match(token)
}

fn split_ws(text: &str) -> Vec<String> {
    text.split_whitespace().map(str::to_owned).collect()
}

fn tokenize(s: String) -> Result<Vec<String>, RequestParseError> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut depth = 0i32;
    for ch in s.chars() {
        match ch {
            '{' => {
                depth += 1;
                cur.push(ch);
            }
            '}' => {
                depth -= 1;
                cur.push(ch);
                if depth < 0 {
                    return Err(RequestParseError::Unbalanced);
                }
            }
            c if c.is_whitespace() && depth == 0 => {
                if !cur.trim().is_empty() {
                    out.push(cur.trim().to_string());
                }
                cur.clear();
            }
            _ => cur.push(ch),
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_string());
    }
    if depth != 0 {
        return Err(RequestParseError::Unbalanced);
    }
    Ok(out)
}

fn parse_final_states(
    tokens: Vec<String>,
) -> Result<(FinalSpec, bool, Vec<String>), RequestParseError> {
    if tokens.is_empty() {
        return Ok((FinalSpec::default(), false, vec![]));
    }
    let mut t = tokens;
    let mut finals = Vec::<String>::new();
    let mut empty_final = false;

    if t[0].eq_ignore_ascii_case("empty-final") || t[0] == "{}" {
        empty_final = true;
        t.remove(0);
        return Ok((FinalSpec::default(), empty_final, t));
    }

    if t[0].starts_with('{') {
        let raw = t.remove(0);
        let sets = parse_alt_sets(&raw)?;
        return Ok((
            FinalSpec {
                primary: vec![],
                sets,
            },
            empty_final,
            t,
        ));
    }

    while !t.is_empty() {
        if is_option_token(&t[0]) {
            break;
        }
        finals.push(t.remove(0));
    }
    Ok((
        FinalSpec {
            primary: finals,
            sets: vec![],
        },
        empty_final,
        t,
    ))
}

#[derive(Default)]
struct FinalSpec {
    primary: Vec<String>,
    sets: Vec<Vec<String>>,
}

fn parse_alt_sets(s: &str) -> Result<Vec<Vec<String>>, RequestParseError> {
    let trimmed = s.trim();
    if !(trimmed.starts_with('{') && trimmed.ends_with('}')) {
        return Err(RequestParseError::InvalidToken(s.to_string()));
    }
    let inner = &trimmed[1..trimmed.len() - 1];
    let mut out: Vec<Vec<String>> = Vec::new();
    let mut cur = String::new();
    let mut depth = 0i32;
    for ch in inner.chars() {
        match ch {
            '{' => {
                depth += 1;
                cur.push(ch);
            }
            '}' => {
                depth -= 1;
                if depth < 0 {
                    return Err(RequestParseError::Unbalanced);
                }
                cur.push(ch);
            }
            ',' if depth == 0 => {
                let items = split_ws(cur.trim());
                if !items.is_empty() {
                    out.push(items);
                } else {
                    out.push(vec![]);
                }
                cur.clear();
            }
            _ => cur.push(ch),
        }
    }
    if depth != 0 {
        return Err(RequestParseError::Unbalanced);
    }
    if !cur.trim().is_empty() {
        out.push(split_ws(cur.trim()));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loops_are_explicit_for_every_request() {
        let model = Model::qcd();
        for kind in [GenerationType::Amplitude, GenerationType::CrossSection] {
            for text in ["21 21 > 21 21", "{} > {}"] {
                assert_eq!(
                    Request::parse(&model, kind, text, None).unwrap_err(),
                    RequestParseError::MissingLoops
                );
            }
        }
        let amplitude = Request::parse(
            &model,
            GenerationType::Amplitude,
            "{} > {}",
            Some(Loops::total(2..=3)),
        )
        .unwrap();
        assert_eq!(amplitude.loops(), &Loops::total(2..=3));
        assert_eq!(
            Request::parse(&model, GenerationType::Amplitude, "{} > {} [{{2}}]", None).unwrap_err(),
            RequestParseError::WrongLoopDomain
        );
        assert_eq!(
            Request::parse(
                &model,
                GenerationType::CrossSection,
                "21 21 > 21 21 [{0}]",
                None
            )
            .unwrap_err(),
            RequestParseError::MissingLoops
        );
        assert_eq!(
            Request::parse(
                &model,
                GenerationType::CrossSection,
                "21 21 > 21 21",
                Some(Loops::total(1..=1))
            )
            .unwrap_err(),
            RequestParseError::WrongLoopDomain
        );
    }

    #[test]
    fn explicit_loop_argument_must_agree_with_the_text() {
        let model = Model::qcd();
        let text = "21 21 > 21 21 [{1}]";
        Request::parse(
            &model,
            GenerationType::Amplitude,
            text,
            Some(Loops::total(1..=1)),
        )
        .unwrap();
        assert_eq!(
            Request::parse(
                &model,
                GenerationType::Amplitude,
                text,
                Some(Loops::total(0..=0))
            )
            .unwrap_err(),
            RequestParseError::ConflictingLoops
        );
        assert_eq!(
            Request::parse(
                &model,
                GenerationType::Amplitude,
                "21 21 > 21 21 [{0} {1}]",
                None
            )
            .unwrap_err(),
            RequestParseError::ConflictingLoops
        );
    }

    #[test]
    fn squared_orders_keep_their_power_and_intersect_with_plain_orders() {
        let model = Model::qcd();
        let request = Request::parse(
            &model,
            GenerationType::Amplitude,
            "21 21 > 21 21 [{1}] QCD>=1 QCD<=3 QCD^2==4",
            None,
        )
        .unwrap();
        assert_eq!(
            request.orders().plain["QCD"],
            Range::Inclusive {
                min: 1,
                max: Some(3)
            }
        );
        assert_eq!(request.orders().squared["QCD"], Range::exact(4));
        let empty = Request::parse(
            &model,
            GenerationType::Amplitude,
            "21 21 > 21 21 [{1}] QCD==2 QCD>=3",
            None,
        )
        .unwrap();
        assert!(empty.orders().is_empty());
        assert!(matches!(
            Request::parse(
                &model,
                GenerationType::Amplitude,
                "21 21 > 21 21 [{1}] QCD^3==4",
                None
            ),
            Err(RequestParseError::InvalidToken(_))
        ));
    }

    #[test]
    fn cross_section_scopes_are_preserved() {
        let model = Model::qcd();
        let request = Request::parse(
            &model,
            GenerationType::CrossSection,
            "21 21 > {21 21, 21 21 21} [{0} {{2}} QCD=1] QCD==2 QCD^2<=6",
            None,
        )
        .unwrap();
        assert_eq!(request.loops(), &Loops::forward(2..=2));
        assert_eq!(request.process().outgoing_alternatives().len(), 2);
        assert_eq!(request.cuts().amplitude.loops, Range::exact(0));
        assert_eq!(
            request.cuts().amplitude.orders.plain["QCD"],
            Range::exact(2)
        );
        assert!(request.orders().plain.is_empty());
        assert_eq!(request.orders().squared["QCD"], Range::at_most(6));
        assert_eq!(request.orders().perturbative["QCD"], 1);
    }

    #[test]
    fn amplitude_alternatives_are_an_error() {
        assert_eq!(
            Request::parse(
                &Model::qcd(),
                GenerationType::Amplitude,
                "21 21 > {21 21, 21 21 21} [{0}]",
                None
            )
            .unwrap_err(),
            RequestParseError::Process(ProcessError::MultipleAmplitudeFinalStates)
        );
    }

    #[test]
    fn empty_and_repeated_allow_clauses_keep_their_meaning() {
        let model = Model::qcd();
        let all = Request::parse(
            &model,
            GenerationType::Amplitude,
            "21 21 > 21 21 [{0}]",
            None,
        )
        .unwrap();
        assert!(all.restrictions().particle_allow.is_empty());
        let none = Request::parse(
            &model,
            GenerationType::Amplitude,
            "21 21 > 21 21 [{0}] only",
            None,
        )
        .unwrap();
        assert_eq!(none.restrictions().particle_allow, vec![BTreeSet::new()]);
        let intersect = Request::parse(
            &model,
            GenerationType::Amplitude,
            "21 21 > 21 21 [{0}] only 21 only",
            None,
        )
        .unwrap();
        assert_eq!(intersect.restrictions().particle_allow.len(), 2);
        assert!(
            intersect
                .restrictions()
                .resolve(&model)
                .unwrap()
                .particles
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn invalid_names_delimiters_and_large_numbers_return_errors() {
        let model = Model::qcd();
        for text in [
            "21 > 21 [{0}",
            "21 > 21 [{0}] QCD^2==184467440737095516160",
            "21 > 21 [{184467440737095516160}]",
            "21 > 21 [{{0}}]",
            "21 > 21 [{0}] [QCD]",
        ] {
            assert!(
                Request::parse(&model, GenerationType::Amplitude, text, None).is_err(),
                "{text}"
            );
        }
        assert!(matches!(
            Request::parse(
                &model,
                GenerationType::Amplitude,
                "unknown > 21 [{0}]",
                None
            ),
            Err(RequestParseError::UnknownParticle { .. })
        ));
        assert!(matches!(
            Request::parse(
                &model,
                GenerationType::Amplitude,
                "99999999 > 21 [{0}]",
                None
            ),
            Err(RequestParseError::UnknownParticle { .. })
        ));
        assert!(Syntax::parse("λ > φ [{0}]").is_ok());
    }
}
