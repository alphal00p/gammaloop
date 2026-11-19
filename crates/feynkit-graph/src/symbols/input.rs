use crate::symbols;
use feynkit_model::Model;
use idenso::color::CS;
use idenso::dirac::AGS;
use idenso::epsilon::EPSILON_SYMBOL;
use idenso::representations::{Bispinor, ColorAdjoint, ColorFundamental, SpinFundamental};
use spenso::network::library::symbolic::ETS;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AIND_SYMBOLS;
use spenso::structure::representation::{Euclidean, LibraryRep, Lorentz, Minkowski};
use std::collections::HashMap;
use std::sync::LazyLock;
use symbolica::atom::{Atom, Symbol};
use symbolica::parser::{ParseSettings, Token};
use symbolica::state::{State, Workspace};
use symbolica::with_default_namespace;

static HEADS: LazyLock<HashMap<String, Symbol>> = LazyLock::new(|| {
    [
        symbols::momentum(),
        symbols::external_momentum(),
        symbols::loop_momentum(),
        symbols::denominator(),
        symbols::edge_index(),
        symbols::vertex_index(),
        symbols::hedge_index(),
        symbols::u(),
        symbols::ubar(),
        symbols::v(),
        symbols::vbar(),
        symbols::epsilon(),
        symbols::epsilonbar(),
        ETS.metric,
        SPENSO_TAG.dot,
        SPENSO_TAG.chain,
        SPENSO_TAG.trace,
        ETS.flat,
        ETS.delta,
        AIND_SYMBOLS.aind,
        AIND_SYMBOLS.uind,
        AIND_SYMBOLS.dind,
        AIND_SYMBOLS.selfdualind,
        AIND_SYMBOLS.openind,
        AIND_SYMBOLS.scope,
        AIND_SYMBOLS.cind,
        AIND_SYMBOLS.find,
        LibraryRep::from(Minkowski {}).symbol(),
        LibraryRep::from(Lorentz {}).symbol(),
        LibraryRep::from(Euclidean {}).symbol(),
        LibraryRep::from(Bispinor {}).symbol(),
        LibraryRep::from(SpinFundamental {}).symbol(),
        LibraryRep::from(ColorFundamental {}).symbol(),
        LibraryRep::from(ColorAdjoint {}).symbol(),
        AGS.gamma,
        AGS.gammaconj,
        AGS.gammaadj,
        AGS.gamma0,
        AGS.gamma5,
        AGS.projp,
        AGS.projm,
        AGS.charge_conjugation,
        AGS.sigma,
        *EPSILON_SYMBOL,
        CS.t,
        CS.f,
        CS.d,
        CS.gram,
        CS.cas,
        CS.idx,
    ]
    .into_iter()
    .map(|head| (head.get_stripped_name().to_owned(), head))
    .collect()
});

/// Parse graph input using canonical graph and tensor heads in function position.
/// Bare scalars can name parameters and couplings of the supplied model.
/// Qualified names, Symbolica built-ins and the caller's namespace for unknown
/// symbols retain their ordinary Symbolica meanings. GammaLoop DOT retains its
/// native spatial `dot` head in the `gammalooprs` namespace.
///
/// ```
/// use feynkit_graph::symbols;
/// use feynkit_model::Model;
/// use symbolica::parse;
/// let weight = symbols::parse_expression("lam*Q(0,cind(0))", "example", Some(&Model::phi4()))?;
/// assert_eq!(weight, parse!("UFO::lam*gammalooprs::Q(0,spenso::cind(0))"));
/// # Ok::<(), String>(())
/// ```
///
/// # Errors
/// Returns the Symbolica parse error for invalid input.
pub fn parse_expression(
    expression: &str,
    default_namespace: &str,
    model: Option<&Model>,
) -> Result<Atom, String> {
    // Resolve tokens before normalization. The same spelling may denote a
    // scalar and a tensor head in one expression, for example g*g(i,j).
    let mut token = Token::parse(expression, ParseSettings::default())?;
    let mut pending = vec![&mut token];
    while let Some(node) = pending.pop() {
        match node {
            Token::Fn(_, _, call) => {
                if let Some((Token::ID(name), args)) = call.split_first_mut() {
                    if let Some(&qualified) = HEADS.get(name.as_str())
                        && !State::is_builtin(name.as_str())
                        && !(qualified == SPENSO_TAG.dot && default_namespace == "gammalooprs")
                    {
                        *name = qualified.get_name().into();
                    }
                    pending.extend(args);
                }
            }
            Token::Op(_, _, _, operands) => pending.extend(operands),
            Token::ID(name) if !name.contains("::") && !State::is_builtin(name.as_str()) => {
                if let Some(parameter) = model.and_then(|model| model.scalar_symbol(name)) {
                    *name = parameter.get_name().into();
                } else if name.as_str() == symbols::dimension().get_stripped_name() {
                    *name = symbols::dimension().get_name().into();
                } else if name.as_str() == CS.nc.get_stripped_name() {
                    *name = CS.nc.get_name().into();
                }
            }
            _ => {}
        }
    }
    let namespace = with_default_namespace!(expression, default_namespace.to_owned());
    Workspace::get_local().with(|ws| token.to_atom(&namespace, &mut Default::default(), ws))
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::atom::AtomCore;

    #[test]
    fn resolves_names_by_role_and_preserves_qualified_symbols() {
        let source = "g*g(mink(dim,1),mink(dim,2))+custom::g(dim)+cos(g)+dot(K(0),P(0))";
        let reference = "g*spenso::g(spenso::mink(gammalooprs::dim,1),spenso::mink(gammalooprs::dim,2))+custom::g(gammalooprs::dim)+cos(g)+spenso::dot(gammalooprs::K(0),gammalooprs::P(0))";
        let expected = Atom::parse(reference, "local", ParseSettings::default()).unwrap();
        let parsed = parse_expression(source, "local", None).unwrap();
        assert_eq!(parsed, expected);
        assert_eq!(
            parse_expression("dot(Q(0),Q(1))", "gammalooprs", None).unwrap(),
            Atom::parse(
                "gammalooprs::dot(Q(0),Q(1))",
                "gammalooprs",
                ParseSettings::default()
            )
            .unwrap()
        );
        assert_eq!(
            parse_expression(&parsed.to_canonical_string(), "other", None).unwrap(),
            expected
        );
    }
}
