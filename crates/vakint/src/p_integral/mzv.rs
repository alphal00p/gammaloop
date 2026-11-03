//! Arbitrary-precision multiple zeta values for `vakint::MZV(n_1, ..., n_k)`.
//!
//! `MZV(n_1,...,n_k) = Σ_(m_1 > ... > m_k >= 1) Π m_i^(-n_i)` (the datamine
//! convention of the Forcer master table) is the iterated integral
//! `I(0; 0^(n_1-1) 1 ... 0^(n_k-1) 1; 1)` with `ω_0 = dt/t`, `ω_1 = dt/(1-t)`.
//! Splitting the path at `1/2` (Hölder convolution, Borwein, Bradley,
//! Broadhurst and Lisoněk, math/9910045) gives
//!
//! ```text
//! I(0; a_1...a_n; 1) = Σ_j I(0; ā_j...ā_1; 1/2) · I(0; a_(j+1)...a_n; 1/2),
//! ```
//!
//! with `ā` exchanging the letters, and every factor is a multiple
//! polylogarithm `Li_(s_1..s_k)(1/2)` whose nested sum converges like `2^-m`.
//! Symbolica's GPL evaluator stops at weight six; the Forcer masters need
//! weight twelve.

use std::collections::BTreeMap;
use std::sync::{LazyLock, Mutex};

use symbolica::atom::{Atom, AtomView, EvaluationInfo, Symbol};
use symbolica::domains::float::{
    Complex, ComplexBall, ErrorQuality, Float, RealBall, TrackedComplex,
};
use symbolica::symbol;

/// The inert tag `vakint::mzv_indices(n_1, ..., n_k)`.
fn indices_symbol() -> Symbol {
    static SYMBOL: LazyLock<Symbol> =
        LazyLock::new(|| symbol!(format!("{}::mzv_indices", crate::NAMESPACE);));
    *SYMBOL
}

/// The numerically evaluable `vakint::MZVvalue(mzv_indices(n_1, ..., n_k))`,
/// the definition of `MZV(n_1, ..., n_k)` for the numerical evaluator.
pub(crate) fn mzv_value_symbol() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        let _ = indices_symbol();
        symbol!(
            format!("{}::MZVvalue", crate::NAMESPACE),
            eval = EvaluationInfo::constant(|tags, bits| {
                mzv_float(&parse_tags(tags)?, bits)
                    .map(|value| Complex::new(value, Float::new(bits)))
            })
            .with_tags(1)
            .register_tracked(|tags, arguments, bits| {
                if !arguments.is_empty() {
                    return Err("MZVvalue expects exactly one symbolic tag".into());
                }
                tracked_mzv(&parse_tags(tags)?, bits)
            })
        )
    });
    *SYMBOL
}

/// `MZVvalue(mzv_indices(indices...))`.
pub(crate) fn mzv_value(indices: &[i64]) -> Atom {
    let tag = indices_symbol().call(&indices.iter().map(|&n| Atom::num(n)).collect::<Vec<_>>());
    symbolica::function!(mzv_value_symbol(), tag)
}

fn parse_tags(tags: &[AtomView]) -> Result<Vec<i64>, String> {
    let [AtomView::Fun(indices)] = tags else {
        return Err("MZVvalue expects mzv_indices(...) as its only argument".into());
    };
    if indices.get_symbol() != indices_symbol() {
        return Err("MZVvalue expects mzv_indices(...) as its only argument".into());
    }
    let indices = indices
        .iter()
        .map(|n| i64::try_from(n).ok())
        .collect::<Option<Vec<_>>>()
        .ok_or("MZV indices must be exact integers")?;
    if indices.is_empty() || indices[0] < 2 || indices.iter().any(|&n| n < 1) {
        return Err(format!("divergent or malformed MZV{indices:?}"));
    }
    Ok(indices)
}

/// Guard bits of the internal summation.
const GUARD_BITS: u32 = 64;

/// `Li_(s_1..s_k)(1/2) = Σ_(m_1 > ... > m_k >= 1) 2^(-m_1) Π m_i^(-s_i)`,
/// truncated at `m_1 <= terms`.
fn polylog_at_half(indices: &[usize], terms: usize, bits: u32) -> Float {
    let one = Float::with_val(bits, 1);
    if indices.is_empty() {
        return one;
    }
    let power = |m: usize, s: usize| {
        (0..s).fold(one.clone(), |value, _| {
            value / Float::with_val(bits, m as u64)
        })
    };
    // inner[m] = Σ over the deeper indices with the deepest sum below m+1.
    let mut inner = vec![one.clone(); terms + 1];
    for &s in indices[1..].iter().rev() {
        let mut cumulative = Float::new(bits);
        let mut next = vec![Float::new(bits); terms + 1];
        for m in 1..=terms {
            cumulative += power(m, s) * &inner[m - 1];
            next[m] = cumulative.clone();
        }
        inner = next;
    }
    let mut total = Float::new(bits);
    let mut half = one.clone();
    for m in 1..=terms {
        half = half / 2;
        total += half.clone() * &power(m, indices[0]) * &inner[m - 1];
    }
    total
}

/// The `{0,1}` letters of a word to its multiple polylogarithm indices
/// (`0^(s-1) 1` is the index `s`). The word must end in `1`.
fn word_indices(letters: &[u8]) -> Vec<usize> {
    let mut indices = Vec::new();
    let mut zeros = 0;
    for &letter in letters {
        if letter == 0 {
            zeros += 1;
        } else {
            indices.push(zeros + 1);
            zeros = 0;
        }
    }
    debug_assert_eq!(zeros, 0, "an MZV word ends in one");
    indices
}

/// The MZV at `bits` of precision (summed with guard bits, the series tail
/// below `2^-(bits + 32)`).
fn mzv_float(indices: &[i64], bits: u32) -> Result<Float, String> {
    type Cache = BTreeMap<(Vec<i64>, u32), Float>;
    static CACHE: LazyLock<Mutex<Cache>> = LazyLock::new(Default::default);
    let key = (indices.to_vec(), bits);
    if let Some(value) = CACHE.lock().map_err(|_| "MZV cache is poisoned")?.get(&key) {
        return Ok(value.clone());
    }
    let working = bits + GUARD_BITS;
    let letters = indices
        .iter()
        .flat_map(|&n| std::iter::repeat_n(0_u8, n as usize - 1).chain([1]))
        .collect::<Vec<_>>();
    // The nested sums are bounded by (1 + ln m)^(k-1) <= 2^24 for the depths
    // and term counts used here, so 2^-terms bounds the tail.
    let terms = (bits + GUARD_BITS / 2 + 32) as usize;
    let mut total = Float::new(working);
    for split in 0..=letters.len() {
        let first = letters[..split]
            .iter()
            .rev()
            .map(|letter| 1 - letter)
            .collect::<Vec<_>>();
        total += polylog_at_half(&word_indices(&first), terms, working)
            * &polylog_at_half(&word_indices(&letters[split..]), terms, working);
    }
    let mut value = total;
    value.set_prec(bits);
    let mut cache = CACHE.lock().map_err(|_| "MZV cache is poisoned")?;
    if cache.len() > 4096 {
        cache.clear();
    }
    cache.insert(key, value.clone());
    Ok(value)
}

fn tracked_mzv(indices: &[i64], bits: u32) -> Result<TrackedComplex, String> {
    let value = mzv_float(indices, bits)?;
    // Rounding of the guarded summation plus the series tail, relative to a
    // value below two: well inside one unit of the requested precision.
    let radius = (0..bits + 4).fold(Float::with_val(bits, 1), |radius, _| radius / 2);
    let ball: ComplexBall = Complex::new(
        RealBall::new(value, radius),
        RealBall::exact(Float::new(bits)),
    );
    Ok(TrackedComplex::from_ball(ball, ErrorQuality::Bounded))
}

#[cfg(test)]
mod tests {
    use symbolica::domains::float::RealLike;

    use super::*;

    fn relative_difference(a: Float, b: &Float) -> f64 {
        ((a - b) / b.clone()).to_f64().abs()
    }

    #[test]
    fn multiple_zeta_values_match_known_digits() {
        let value = |indices: &[i64]| mzv_float(indices, 200).unwrap();
        let known = |digits: &str| Float::parse(digits, Some(200)).unwrap();
        // ζ(2) = π²/6.
        assert!(
            relative_difference(
                value(&[2]),
                &known("1.64493406684822643647241516664602518921894990120679843773555823")
            ) < 1e-55
        );
        // ζ(5,3) (datamine convention, m_1 > m_2).
        assert!(
            relative_difference(
                value(&[5, 3]),
                &known("0.037707672984847544011304782293659")
            ) < 1e-30
        );
        // Stuffle: ζ(3)ζ(5) = ζ(5,3) + ζ(3,5) + ζ(8).
        let stuffle = value(&[5, 3]) + &value(&[3, 5]) + &value(&[8]);
        assert!(relative_difference(stuffle, &(value(&[3]) * &value(&[5]))) < 1e-55);
        // Duality: ζ(2,1) = ζ(3).
        assert!(relative_difference(value(&[2, 1]), &value(&[3])) < 1e-55);
    }
}
