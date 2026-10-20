//! Explicit, portable adoption of independent index spaces by the native color owner.

use std::collections::HashMap;

use spenso::{
    network::tags::SPENSO_TAG,
    structure::{
        partial::PartialStructure,
        representation::{LibraryRep, RepName, Representation},
    },
};
use symbolica::atom::{
    Atom, AtomCore, AtomView, NamespacedSymbol, Symbol, SymbolBuilder, UserData,
    representation::FunView,
};

use super::{CS, ColorSimplifySettings};
use crate::tensor::{AlgebraSettings, SymbolicTensor};

const TAG: &str = "idenso::color-family:v1";

#[cfg(test)]
mod tests;

/// The native Hermitian-generator convention `[T^a,T^b] = i f^{abc} T^c`.
/// `Tr(t^a t^b) = I2(F) delta^{ab}` remains symbolic by default. The native
/// explicit SU(N) normalization sets `I2(F) = 1/2`. Dual generators are
/// `-(t^a)^T`, and adjoint matrix elements are `-i f^{abc}`.
///
/// An SU(2) fundamental keeps its distinct dual space. Identifying them requires
/// an explicitly supplied epsilon intertwiner; dimension two never selects it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ColorFamilyConvention {
    SymbolicSuN,
}

/// Complete immutable association between scoped spaces and a native color family.
/// Dimensions do not select a family: callers explicitly supply both representation
/// identities and both tensor-head names. Component libraries are a separate capability.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ColorFamilyDeclaration {
    fundamental: LibraryRep,
    adjoint: LibraryRep,
    heads: [String; 2],
    nc: Symbol,
    convention: ColorFamilyConvention,
}

impl ColorFamilyDeclaration {
    /// Declare the native symbolic SU(N) convention on independent scoped spaces.
    /// An independent dimension symbol is mandatory; native `Nc` is never shared
    /// implicitly between two families. Concrete slot dimensions remain explicit.
    pub fn su_n(
        fundamental: LibraryRep,
        adjoint: LibraryRep,
        heads: [String; 2],
        nc: Symbol,
    ) -> Result<Self, String> {
        if fundamental.is_self_dual()
            || !fundamental.is_base()
            || !matches!(adjoint, LibraryRep::SelfDual(_))
        {
            return Err(
                "scoped color requires a base dualizable fundamental and a self-dual adjoint"
                    .into(),
            );
        }
        if fundamental.symbol() == CS.fundamental_rep || adjoint.symbol() == CS.adjoint_rep {
            return Err("scoped color requires independent representation identities; use the built-in owner directly for built-in spaces".into());
        }
        let mut names = std::collections::HashSet::new();
        for name in &heads {
            if !name.contains("::")
                || !names.insert(name)
                || [CS.t, CS.f].iter().any(|s| s.get_name() == name)
            {
                return Err(
                    "scoped color requires two distinct namespaced tensor-head names".into(),
                );
            }
        }
        if !nc.is_real()
            || nc.has_tag(&SPENSO_TAG.tensor)
            || nc.has_tag(&SPENSO_TAG.representation)
            || [CS.nc, CS.t, CS.f, fundamental.symbol(), adjoint.symbol()].contains(&nc)
            || heads.iter().any(|name| name == nc.get_name())
        {
            return Err(
                "scoped color dimension symbol must have an independent real scalar identity"
                    .into(),
            );
        }
        Ok(Self {
            fundamental,
            adjoint,
            heads,
            nc,
            convention: ColorFamilyConvention::SymbolicSuN,
        })
    }

    pub fn fundamental(&self) -> LibraryRep {
        self.fundamental
    }
    pub fn adjoint(&self) -> LibraryRep {
        self.adjoint
    }
    pub fn number_of_colors(&self) -> Symbol {
        self.nc
    }
    pub fn convention(&self) -> ColorFamilyConvention {
        self.convention
    }
    pub fn head_names(&self) -> &[String; 2] {
        &self.heads
    }

    /// Serialize native symbol references as atoms, so Symbolica archives retain
    /// their representation metadata and remap references during state import.
    pub fn to_user_data(&self) -> UserData {
        UserData::List(vec![
            UserData::String(TAG.into()),
            UserData::String("symbolic-su-n".into()),
            UserData::Atom(Atom::var(self.fundamental.symbol())),
            UserData::Atom(Atom::var(self.adjoint.symbol())),
            UserData::Atom(Atom::var(self.nc)),
            UserData::List(self.heads.iter().cloned().map(UserData::String).collect()),
        ])
    }

    pub fn from_user_data(data: &UserData) -> Result<Self, String> {
        let UserData::List(data) = data else {
            return Err("invalid scoped color declaration".into());
        };
        let [
            UserData::String(version),
            UserData::String(convention),
            fundamental,
            adjoint,
            nc,
            UserData::List(heads),
        ] = data.as_slice()
        else {
            return Err("invalid scoped color declaration fields".into());
        };
        if version != TAG || convention != "symbolic-su-n" {
            return Err("unsupported scoped color declaration version or convention".into());
        }
        fn symbol(data: &UserData) -> Result<Symbol, String> {
            if let UserData::Atom(atom) = data
                && let AtomView::Var(v) = atom.as_view()
            {
                return Ok(v.get_symbol());
            }
            Err("scoped color symbol reference must be a native variable atom".into())
        }
        let fundamental =
            LibraryRep::try_from_symbol_coerced(symbol(fundamental)?).map_err(|e| e.to_string())?;
        let adjoint =
            LibraryRep::try_from_symbol_coerced(symbol(adjoint)?).map_err(|e| e.to_string())?;
        let nc = symbol(nc)?;
        let heads = heads
            .iter()
            .map(|h| {
                if let UserData::String(h) = h {
                    Ok(h.clone())
                } else {
                    Err("scoped color head name must be a string".to_string())
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::su_n(
            fundamental,
            adjoint,
            heads
                .try_into()
                .map_err(|_| "scoped color requires exactly two tensor heads")?,
            nc,
        )
    }

    /// Register or restore all family heads from the complete persisted declaration.
    pub fn register(&self) -> Result<ScopedColorFamily, String> {
        let native = [CS.t, CS.f];
        let mut heads = Vec::new();
        for (role, (name, source)) in self.heads.iter().zip(native).enumerate() {
            let mut tags = source.get_tags().to_vec();
            tags.push(TAG.into());
            let data = UserData::List(vec![UserData::Integer(role as i64), self.to_user_data()]);
            let head = SymbolBuilder::new(NamespacedSymbol::parse(name))
                .with_attributes(source.get_attributes().to_vec())
                .with_tags(tags)
                .with_user_data(data)
                .build()
                .map_err(|e| format!("cannot register scoped color head {name}: {e}"))?;
            heads.push(head);
        }
        Ok(ScopedColorFamily {
            declaration: self.clone(),
            heads: heads.try_into().expect("two heads"),
        })
    }

    /// Remap the owned declaration, including its explicit head-name payload.
    /// Generic Atom renaming alone cannot change opaque strings in symbol data.
    pub fn remapped(
        &self,
        symbols: &HashMap<Symbol, Symbol>,
        heads: [String; 2],
    ) -> Result<Self, String> {
        let mapped = |s: Symbol| symbols.get(&s).copied().unwrap_or(s);
        Self::su_n(
            LibraryRep::try_from_symbol_coerced(mapped(self.fundamental.symbol()))
                .map_err(|e| e.to_string())?,
            LibraryRep::try_from_symbol_coerced(mapped(self.adjoint.symbol()))
                .map_err(|e| e.to_string())?,
            heads,
            mapped(self.nc),
        )
    }
}

/// Registered native color action on independent scoped representation identities.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScopedColorFamily {
    declaration: ColorFamilyDeclaration,
    heads: [Symbol; 2],
}

impl ScopedColorFamily {
    pub fn declaration(&self) -> &ColorFamilyDeclaration {
        &self.declaration
    }
    pub fn generator_symbol(&self) -> Symbol {
        self.heads[0]
    }
    pub fn structure_symbol(&self) -> Symbol {
        self.heads[1]
    }

    /// Recover a family's complete declaration from any of its tensor heads.
    /// Missing or inconsistent metadata is diagnosed; names never imply an action.
    pub fn from_symbol(symbol: Symbol) -> Result<Option<Self>, String> {
        if symbol
            .get_tags()
            .iter()
            .any(|tag| tag.starts_with("idenso::color-family:") && tag != TAG)
        {
            return Err("unsupported scoped color metadata version".into());
        }
        if !symbol.has_tag(TAG) {
            return Ok(None);
        }
        let UserData::List(fields) = symbol.get_data() else {
            return Err("scoped color head lacks its declaration".into());
        };
        let [UserData::Integer(role), data] = fields.as_slice() else {
            return Err("invalid scoped color head metadata".into());
        };
        let role = usize::try_from(*role)
            .ok()
            .filter(|r| *r < 2)
            .ok_or("invalid scoped color head role")?;
        let declaration = ColorFamilyDeclaration::from_user_data(data)?;
        if declaration.heads[role] != symbol.get_name() {
            return Err("scoped color head was renamed without remapping its declaration payload; restore a remapped native declaration first".into());
        }
        declaration.register().map(Some)
    }

    /// Recover actions from generator and structure heads in a native expression.
    /// Reduced invariants use shared native heads with explicit representation
    /// arguments; persist their declaration separately when further reduction is needed.
    pub fn in_expression(expression: AtomView<'_>) -> Result<Vec<Self>, String> {
        let mut families = Vec::new();
        let mut error = None;
        expression.visitor(&mut |node| {
            if error.is_some() {
                return false;
            }
            if let AtomView::Fun(function) = node {
                if function.get_symbol().is_scalar() {
                    return false;
                }
                match Self::from_symbol(function.get_symbol()) {
                    Ok(Some(family)) if !families.contains(&family) => families.push(family),
                    Err(e) => {
                        error = Some(e);
                        return false;
                    }
                    _ => {}
                }
            }
            true
        });
        error.map_or(Ok(families), Err)
    }

    /// Validate a declared tensor's native family roles before typed admission.
    pub fn validate_function(function: FunView<'_>) -> Result<(), String> {
        Self::validate_arguments(function.get_symbol(), &function.iter().collect::<Vec<_>>())
    }

    pub(crate) fn validate_arguments(
        symbol: Symbol,
        arguments: &[AtomView<'_>],
    ) -> Result<(), String> {
        let Some(family) = Self::from_symbol(symbol)? else {
            return Ok(());
        };
        // Reuse the native signature owner after one complete symbol adoption.
        // The preliminary role check prevents unrelated built-in spaces from
        // acquiring this family's semantics merely by occupying a valid port.
        let generator = symbol == family.heads[0];
        if arguments.len() != 3 {
            return Err(
                "scoped color generators and structure constants require three slots".into(),
            );
        }
        for (position, argument) in arguments.iter().enumerate() {
            if matches!(argument, AtomView::Var(v) if [SPENSO_TAG.chain_in, SPENSO_TAG.chain_out].contains(&v.get_symbol()))
                && (!generator || position > 0)
            {
                continue;
            }
            let representation = Representation::<LibraryRep>::try_from(*argument).map_err(
                |_| "scoped color tensor arguments must be native slots or matrix placeholders",
            )?;
            let expected = if !generator || position == 0 {
                family.declaration.adjoint
            } else if position == 1 {
                family.declaration.fundamental
            } else {
                family.declaration.fundamental.dual()
            };
            if representation.rep != expected {
                return Err(
                    "scoped color tensor slot does not belong to its declared representation role"
                        .into(),
                );
            }
        }
        Ok(())
    }

    pub(crate) fn native_to_scoped(&self) -> HashMap<Symbol, Symbol> {
        let mut map = HashMap::from([
            (CS.fundamental_rep, self.declaration.fundamental.symbol()),
            (CS.adjoint_rep, self.declaration.adjoint.symbol()),
            (CS.t, self.heads[0]),
            (CS.f, self.heads[1]),
        ]);
        map.insert(CS.nc, self.declaration.nc);
        map
    }

    pub(crate) fn scoped_to_native(&self) -> HashMap<Symbol, Symbol> {
        self.native_to_scoped()
            .into_iter()
            .map(|(a, b)| (b, a))
            .collect()
    }

    /// Produce the complete map needed to move an expression with its declaration.
    pub fn symbol_map_to(&self, target: &Self) -> Result<HashMap<Symbol, Symbol>, String> {
        let target = target.native_to_scoped();
        Ok(self
            .native_to_scoped()
            .into_iter()
            .map(|(native, scoped)| (scoped, target[&native]))
            .collect())
    }

    /// Check that an explicitly requested representation belongs to this action.
    /// Arbitrary representations need a declared native generator template/backend.
    pub fn validate_representation(&self, representation: LibraryRep) -> Result<(), String> {
        if representation.base() == self.declaration.fundamental
            || representation == self.declaration.adjoint
        {
            Ok(())
        } else {
            Err("representation has no action in this scoped color family; supply an explicit native generator template/backend".into())
        }
    }

    /// Construct the fundamental Hermitian generator using the existing native grammar.
    pub fn fundamental_generator(
        &self,
        adjoint: Atom,
        left: Atom,
        right: Atom,
    ) -> Result<SymbolicTensor<PartialStructure>, String> {
        let atom = CS
            .explicit_t(adjoint, left, right)
            .rename_symbols(&HashMap::from([(CS.t, self.heads[0])]));
        SymbolicTensor::infer(atom).map_err(|e| e.to_string())
    }

    /// Construct a generator in the declared fundamental, dual, or adjoint action.
    /// Matrix slots are `(left, right)`; the dual action is `-t^T` and the
    /// adjoint action is `-i f`. Their index identities remain native slots.
    pub fn generator(
        &self,
        representation: LibraryRep,
        adjoint: Atom,
        left: Atom,
        right: Atom,
    ) -> Result<SymbolicTensor<PartialStructure>, String> {
        self.validate_representation(representation)?;
        if representation == self.declaration.fundamental {
            self.fundamental_generator(adjoint, left, right)
        } else if representation == self.declaration.fundamental.dual() {
            Ok(-self.fundamental_generator(adjoint, right, left)?)
        } else {
            let structure = self.structure_constant(adjoint, left, right)?;
            structure
                .with_checked_expression(-Atom::i() * structure.expression())
                .map_err(|e| e.to_string())
        }
    }

    /// Construct the native real antisymmetric structure constant.
    pub fn structure_constant(
        &self,
        a: Atom,
        b: Atom,
        c: Atom,
    ) -> Result<SymbolicTensor<PartialStructure>, String> {
        let atom = CS
            .structure_f(a, b, c)
            .rename_symbols(&HashMap::from([(CS.f, self.heads[1])]));
        SymbolicTensor::infer(atom).map_err(|e| e.to_string())
    }

    /// Reduce one selected family and immediately restore its scoped identities.
    /// Native renaming owns injectivity/collision checks and native algebra owns all identities.
    pub fn reduce(
        &self,
        value: &SymbolicTensor<PartialStructure>,
        settings: ColorSimplifySettings,
    ) -> Result<SymbolicTensor<PartialStructure>, String> {
        value
            .rename_symbols(&self.scoped_to_native())
            .map_err(|e| e.to_string())?
            .simplify_algebra(&AlgebraSettings {
                color: Some(settings),
                ..Default::default()
            })
            .map_err(|e| e.to_string())?
            .rename_symbols(&self.native_to_scoped())
            .map_err(|e| e.to_string())
    }
}
