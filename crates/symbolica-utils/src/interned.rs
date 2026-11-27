//! Process-wide interning of Symbolica atoms behind a small copyable handle.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fmt::Result as FmtResult;
use std::fmt::{Debug, Display, Formatter};
use std::hash::{BuildHasher, Hash, Hasher, RandomState};
use std::sync::{LazyLock, RwLock, RwLockReadGuard, RwLockWriteGuard};

use bincode::de::{BorrowDecoder, Decoder};
use bincode::enc::Encoder;
use bincode::error::{DecodeError, EncodeError};
use bincode::{BorrowDecode, Decode, Encode};
use serde::de::Error as DeError;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::state::HasStateMap;
use symbolica::try_parse;

/// The handles of the interned atoms that share one hash.
enum Bucket {
    One(u32),
    Many(Vec<u32>),
}

#[derive(Default)]
struct AtomInterner {
    /// Interned atoms by handle number. They are never released, so their
    /// references outlive every lock guard.
    atoms: Vec<&'static Atom>,
    buckets: HashMap<u64, Bucket>,
    hasher: RandomState,
}

impl AtomInterner {
    fn read() -> RwLockReadGuard<'static, Self> {
        INTERNER
            .read()
            .expect("the atom interner lock must not be poisoned")
    }

    fn write() -> RwLockWriteGuard<'static, Self> {
        INTERNER
            .write()
            .expect("the atom interner lock must not be poisoned")
    }

    /// The handle number of `view`, whose hash is `hash`, if it is interned.
    fn find(&self, hash: u64, view: AtomView<'_>) -> Option<u32> {
        let interns = |id: &u32| self.atoms[*id as usize].as_view() == view;
        match self.buckets.get(&hash)? {
            Bucket::One(id) => interns(id).then_some(*id),
            Bucket::Many(ids) => ids.iter().copied().find(interns),
        }
    }

    /// Interns `atom`, which is not interned yet and hashes to `hash`.
    ///
    /// # Panics
    ///
    /// Panics if more than `u32::MAX` distinct atoms are interned.
    fn insert(&mut self, hash: u64, atom: Atom) -> u32 {
        let id = u32::try_from(self.atoms.len())
            .expect("the atom interner holds fewer than u32::MAX atoms");
        self.atoms.push(Box::leak(Box::new(atom)));
        match self.buckets.entry(hash) {
            Entry::Vacant(entry) => {
                entry.insert(Bucket::One(id));
            }
            Entry::Occupied(mut entry) => {
                let bucket = entry.get_mut();
                match bucket {
                    Bucket::One(first) => *bucket = Bucket::Many(vec![*first, id]),
                    Bucket::Many(ids) => ids.push(id),
                }
            }
        }
        id
    }
}

static INTERNER: LazyLock<RwLock<AtomInterner>> = LazyLock::new(Default::default);

/// A copyable handle to a Symbolica atom interned for the lifetime of the process.
///
/// Interning is bijective: two handles are equal exactly when their atoms are
/// equal, so equality and hashing use the handle alone. The handle number
/// depends on the order in which threads first intern atoms. `Ord` therefore
/// compares the underlying [`AtomView`]s instead: the order is independent of
/// interning order, consistent with equality, and a deterministic total order
/// within a Symbolica state, since [`AtomView`] ordering follows the symbol
/// order of the active state. Serialization stores the atom, never the number.
///
/// ```
/// use symbolica::parse;
/// use symbolica_utils::InternedAtomId;
///
/// let dimension = parse!("interned_doc::nc^2 - 1");
/// let id = InternedAtomId::intern(dimension.as_view());
///
/// assert_eq!(id, InternedAtomId::intern(dimension.as_view()));
/// assert_eq!(id.to_atom(), dimension);
/// ```
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct InternedAtomId(u32);

impl InternedAtomId {
    /// Intern `view`, returning the handle shared by every equal atom.
    ///
    /// An atom that is already interned is found without copying `view`.
    ///
    /// # Panics
    ///
    /// Panics if the interner lock was poisoned or more than `u32::MAX`
    /// distinct atoms are interned.
    #[must_use]
    pub fn intern(view: AtomView<'_>) -> Self {
        let hash = {
            let interner = AtomInterner::read();
            let hash = interner.hasher.hash_one(view);
            if let Some(id) = interner.find(hash, view) {
                return Self(id);
            }
            hash
        };
        let mut interner = AtomInterner::write();
        // Another thread may have interned the atom between the two locks.
        if let Some(id) = interner.find(hash, view) {
            return Self(id);
        }
        Self(interner.insert(hash, view.to_owned()))
    }

    /// The interned atom. The lock is held only to look up the reference, so
    /// callers never run while holding it.
    fn atom(self) -> &'static Atom {
        AtomInterner::read().atoms[self.0 as usize]
    }

    /// The handle number, unique per interned atom within this process.
    ///
    /// The number depends on interning order; never persist or order by it.
    #[must_use]
    pub fn get_id(self) -> u32 {
        self.0
    }

    /// A copy of the interned atom.
    ///
    /// # Panics
    ///
    /// Panics if the interner lock was poisoned.
    pub fn to_atom(self) -> Atom {
        self.atom().clone()
    }
}

impl Hash for InternedAtomId {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Equal handles name equal atoms, so the number is a consistent key.
        self.0.hash(state);
    }
}

impl Ord for InternedAtomId {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.0 == other.0 {
            return Ordering::Equal;
        }
        self.atom().as_view().cmp(&other.atom().as_view())
    }
}

impl PartialOrd for InternedAtomId {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Debug for InternedAtomId {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.debug_tuple("InternedAtomId").field(self.atom()).finish()
    }
}

impl Display for InternedAtomId {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        Display::fmt(self.atom(), f)
    }
}

impl From<AtomView<'_>> for InternedAtomId {
    fn from(view: AtomView<'_>) -> Self {
        Self::intern(view)
    }
}

impl From<&Atom> for InternedAtomId {
    fn from(atom: &Atom) -> Self {
        Self::intern(atom.as_view())
    }
}

impl From<InternedAtomId> for Atom {
    fn from(id: InternedAtomId) -> Self {
        id.to_atom()
    }
}

impl Serialize for InternedAtomId {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.atom().to_canonical_string().serialize(serializer)
    }
}

impl<'d> Deserialize<'d> for InternedAtomId {
    fn deserialize<D: Deserializer<'d>>(deserializer: D) -> Result<Self, D::Error> {
        let value = String::deserialize(deserializer)?;
        let atom = try_parse!(&value).map_err(DeError::custom)?;
        Ok(Self::intern(atom.as_view()))
    }
}

impl Encode for InternedAtomId {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        self.atom().encode(encoder)
    }
}

impl<C: HasStateMap> Decode<C> for InternedAtomId {
    fn decode<D: Decoder<Context = C>>(decoder: &mut D) -> Result<Self, DecodeError> {
        // The atom decoder translates symbols into the current state.
        let atom = Atom::decode(decoder)?;
        Ok(Self::intern(atom.as_view()))
    }
}

impl<'de, C: HasStateMap> BorrowDecode<'de, C> for InternedAtomId {
    fn borrow_decode<D: BorrowDecoder<'de, Context = C>>(
        decoder: &mut D,
    ) -> Result<Self, DecodeError> {
        Self::decode(decoder)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::env;
    use std::fs;
    use std::io::Cursor;
    use std::path::PathBuf;
    use std::process;
    use std::process::Command;

    use bincode::config;
    use bincode::serde as bincode_serde;
    use symbolica::state::State;
    use symbolica::{parse, symbol};

    use super::*;

    #[test]
    fn interning_is_bijective() {
        let first = parse!("interned_test::bijective_a^2 - 1");
        let second = parse!("interned_test::bijective_b + 4");
        let ids = [&first, &second, &first].map(InternedAtomId::from);

        assert_eq!(ids[0], ids[2]);
        assert_ne!(ids[0], ids[1]);
        assert_eq!(HashSet::from(ids).len(), 2);
        assert_eq!(ids[0].to_atom(), first);
        assert_eq!(Atom::from(ids[1]), second);
    }

    #[test]
    fn ordering_follows_atoms_not_insertion_order() {
        let atoms = ["ordering_a", "ordering_b", "ordering_c", "ordering_d"]
            .map(|name| parse!(&format!("interned_test::{name}^2 - 1")));
        // Intern half in reverse order so handle numbers disagree with atoms.
        let late = [&atoms[3], &atoms[2]].map(InternedAtomId::from);
        let early = [&atoms[0], &atoms[1]].map(InternedAtomId::from);
        let ids = [early[0], early[1], late[1], late[0]];

        for (left, left_atom) in ids.iter().zip(&atoms) {
            for (right, right_atom) in ids.iter().zip(&atoms) {
                assert_eq!(left.cmp(right), left_atom.cmp(right_atom));
            }
        }
    }

    #[test]
    fn serde_round_trips_through_the_atom() {
        let atom = parse!("interned_test::serde_f(interned_test::serde_x, 3)");
        let id = InternedAtomId::intern(atom.as_view());
        let bytes = bincode_serde::encode_to_vec(id, config::standard()).unwrap();
        let (decoded, _): (InternedAtomId, usize) =
            bincode_serde::decode_from_slice(&bytes, config::standard()).unwrap();

        assert_eq!(decoded, id);
    }

    #[test]
    fn bincode_round_trips_with_a_state_map() {
        let atom = parse!("interned_test::state_n^2 - 1");
        let id = InternedAtomId::intern(atom.as_view());
        let bytes = bincode::encode_to_vec(id, config::standard()).unwrap();
        let mut symbols = Vec::new();
        State::export(&mut symbols).unwrap();
        let state_map = State::import(&mut Cursor::new(symbols), None).unwrap();
        let (decoded, consumed): (InternedAtomId, usize) =
            bincode::decode_from_slice_with_context(&bytes, config::standard(), state_map).unwrap();

        assert_eq!(consumed, bytes.len());
        assert_eq!(decoded, id);
    }

    /// Names the directory holding the exported state and encoded handle that
    /// [`bincode_decodes_after_the_symbol_ids_shift`] hands to its child.
    const SHIFTED_STATE_DIR: &str = "SYMBOLICA_UTILS_INTERNED_SHIFTED_STATE_DIR";
    const SHIFTED_SYMBOL: &str = "interned_test::shifted_nc";
    const SHIFTED_DIMENSION: &str = "interned_test::shifted_nc^2 - 1";

    #[test]
    fn bincode_decodes_after_the_symbol_ids_shift() {
        let id = InternedAtomId::intern(parse!(SHIFTED_DIMENSION).as_view());
        let dir = env::temp_dir().join(format!("symbolica-utils-interned-{}", process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut state = Vec::new();
        State::export(&mut state).unwrap();
        fs::write(dir.join("state"), state).unwrap();
        let handle = bincode::encode_to_vec(id, config::standard()).unwrap();
        fs::write(dir.join("handle"), handle).unwrap();
        let symbol_id = symbol!(SHIFTED_SYMBOL).get_id();
        fs::write(dir.join("symbol_id"), symbol_id.to_string()).unwrap();

        let output = Command::new(env::current_exe().unwrap())
            .args([
                "--exact",
                "interned::tests::decodes_in_a_shifted_state",
                "--nocapture",
            ])
            .env(SHIFTED_STATE_DIR, &dir)
            .output()
            .unwrap();
        fs::remove_dir_all(&dir).unwrap();

        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{stdout}\n{stderr}");
        assert!(
            stdout.contains("1 passed"),
            "the child ran no test:\n{stdout}"
        );
    }

    /// The child of [`bincode_decodes_after_the_symbol_ids_shift`], and a no-op
    /// without its directory. It registers fresh symbols before importing the
    /// parent's state, so the imported symbols take other identifiers.
    #[test]
    fn decodes_in_a_shifted_state() {
        let Some(dir) = env::var_os(SHIFTED_STATE_DIR).map(PathBuf::from) else {
            return;
        };
        for index in 0..64 {
            symbol!(&format!("interned_test::shift_{index}"));
        }
        let state = fs::read(dir.join("state")).unwrap();
        let state_map = State::import(&mut Cursor::new(state), None).unwrap();
        let bytes = fs::read(dir.join("handle")).unwrap();
        let (decoded, consumed): (InternedAtomId, usize) =
            bincode::decode_from_slice_with_context(&bytes, config::standard(), state_map).unwrap();
        let parent_symbol_id: u32 = fs::read_to_string(dir.join("symbol_id"))
            .unwrap()
            .parse()
            .unwrap();

        assert_ne!(symbol!(SHIFTED_SYMBOL).get_id(), parent_symbol_id);
        assert_eq!(consumed, bytes.len());
        assert_eq!(decoded.to_atom(), parse!(SHIFTED_DIMENSION));
    }

    #[test]
    fn colliding_hashes_keep_atoms_apart() {
        let atoms = ["collision_a", "collision_b"]
            .map(|name| parse!(&format!("interned_test::{name} + 1")));
        // Force both atoms into one bucket, as a hash collision would.
        let hash = 0;
        let mut interner = AtomInterner::default();
        for atom in &atoms {
            interner.insert(hash, atom.clone());
        }
        let found = atoms
            .each_ref()
            .map(|atom| interner.find(hash, atom.as_view()));

        assert!(found[0].is_some() && found[1].is_some());
        assert_ne!(found[0], found[1]);
    }
}
