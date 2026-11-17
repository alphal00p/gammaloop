//! Tests of the five-loop p master basis change (lane vk5l-mastermap).
//!
//! Besides the data checks, `inverse_at` is the inverse direction as a test
//! utility (pole-bounded -> original): exact at a rational point `ε0`, where
//! `B(ε0)` is inverted exactly over Q; a round trip original -> pole-bounded
//! -> original is asserted at several points on real expressions.

use std::collections::{BTreeMap, BTreeSet};

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::id::Replacement;

use super::super::master_catalogue::{PMasterCatalogue, PMasterSource};
use super::super::p_master_symbol;
use super::*;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};

fn epsilon() -> Symbol {
    vk_symbol!("ep")
}

fn catalogue() -> PMasterCatalogue {
    crate::Vakint::initialize_vakint_symbols();
    PMasterCatalogue::from_json(
        "built-in p-master-catalogue-5l.json",
        include_str!("../../data/p_master/p-master-catalogue-5l.json"),
    )
    .unwrap()
}

fn master(catalogue: &PMasterCatalogue, id: &str) -> Atom {
    p_master_symbol().call(&[Atom::var(vk_symbol!(catalogue.symbol_of_id(id).unwrap()))][..])
}

/// `atom` at `ε = point` (a rational number), expanded.
fn at(atom: &Atom, point: &Atom) -> Atom {
    atom.replace(Atom::var(epsilon()).to_pattern())
        .with(point.to_pattern())
        .expand()
}

/// `(monomial in the masters, coefficient)` of `expression`, expanded.
fn masters_of(expression: &Atom) -> Vec<(Atom, Atom)> {
    let expanded = expression.expand();
    let mut masters = BTreeSet::new();
    expanded.visitor(&mut |part: AtomView| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == p_master_symbol()
        {
            masters.insert(part.to_owned());
            return false;
        }
        true
    });
    let masters: Vec<Atom> = masters.into_iter().collect();
    expanded.coefficient_list::<i8>(&masters)
}

/// The coefficients of `expression` per master `PMaster(..)` (exact atoms).
fn by_master(expression: &Atom) -> BTreeMap<String, Atom> {
    let expanded = expression.expand();
    let mut masters = BTreeSet::new();
    expanded.visitor(&mut |part: AtomView| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == p_master_symbol()
        {
            masters.insert(part.to_owned());
            return false;
        }
        true
    });
    let masters: Vec<Atom> = masters.into_iter().collect();
    let mut out = BTreeMap::new();
    for (monomial, coefficient) in expanded.coefficient_list::<i8>(&masters) {
        *out.entry(monomial.to_canonical_string())
            .or_insert(Atom::Zero) += coefficient;
    }
    out
}

/// The inverse direction at `ε = point`: every new master of the first layer
/// in old masters, `new = X(old - A rest)`, `X = B(point)^-1` exact over Q.
/// Replaces the new masters of `rotated` (already evaluated at `point`).
fn inverse_at(
    basis: &PMasterBasisData,
    catalogue: &PMasterCatalogue,
    rotated_at_point: &Atom,
    point: &Atom,
) -> Atom {
    let rules = basis.rules(catalogue, epsilon()).unwrap();
    // the layers undone last to first, each through the inverse of its block
    let mut ends = Vec::new();
    let mut end = 0;
    for layer in &basis.layers {
        ends.push((end, end + layer.exchanges.len()));
        end += layer.exchanges.len();
    }
    let mut back = rotated_at_point.clone();
    for (layer, &(start, end)) in basis.layers.iter().zip(&ends).rev() {
        back = inverse_layer_at(
            catalogue,
            &rules[start..end],
            &layer.new_masters,
            &back,
            point,
        );
    }
    back
}

/// One layer undone at `point`: `old = B new + rest` (its rules) solved for
/// its new masters `news` (B square) and substituted into `expression`.
fn inverse_layer_at(
    catalogue: &PMasterCatalogue,
    rules: &[(Atom, Atom)],
    news: &[String],
    expression: &Atom,
    point: &Atom,
) -> Atom {
    let n = rules.len();
    assert_eq!(news.len(), n, "B must be square");
    let new_atoms: Vec<Atom> = news.iter().map(|id| master(catalogue, id)).collect();
    // old_m = sum_n B_mn new_n + rest_m  (all at the point)
    let mut matrix = vec![vec![Atom::Zero; n]; n];
    let mut rest = vec![Atom::Zero; n];
    for (row, (_, value)) in rules.iter().enumerate() {
        for (monomial, coefficient) in masters_of(&at(value, point)) {
            match new_atoms.iter().position(|new| *new == monomial) {
                Some(col) => matrix[row][col] += coefficient,
                None => rest[row] += coefficient * monomial,
            }
        }
    }
    // Gauss-Jordan on [B | I] over Q.
    let mut inverse: Vec<Vec<Atom>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| if i == j { Atom::num(1) } else { Atom::Zero })
                .collect()
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n)
            .find(|&r| !matrix[r][col].is_zero())
            .expect("B(point) singular");
        matrix.swap(col, pivot);
        inverse.swap(col, pivot);
        let p = matrix[col][col].clone();
        for j in 0..n {
            matrix[col][j] = (&matrix[col][j] / &p).expand();
            inverse[col][j] = (&inverse[col][j] / &p).expand();
        }
        for r in 0..n {
            if r != col && !matrix[r][col].is_zero() {
                let f = matrix[r][col].clone();
                for j in 0..n {
                    matrix[r][j] = (&matrix[r][j] - &f * &matrix[col][j]).expand();
                    inverse[r][j] = (&inverse[r][j] - &f * &inverse[col][j]).expand();
                }
            }
        }
    }
    // new_n = sum_m X_nm (old_m - rest_m)
    let olds: Vec<Atom> = rules.iter().map(|(old, _)| old.clone()).collect();
    let replacements: Vec<Replacement> = (0..n)
        .map(|i| {
            let mut value = Atom::Zero;
            for m in 0..n {
                value += &inverse[i][m] * (&olds[m] - &rest[m]);
            }
            Replacement::new(new_atoms[i].to_pattern(), value.expand())
        })
        .collect();
    expression.replace_multiple(&replacements).expand()
}

/// The in-repo fixture: the masters directory restricted to the one
/// exchange T201 meets (int137) and its 8 new masters.
fn fixture_directory() -> std::path::PathBuf {
    std::path::PathBuf::from(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/src/p_integral/master_basis_fixture"
    ))
}

/// The complete masters directory (VAKINT_TEST_P_MASTERS_DIRECTORY, e.g. a
/// five-loop package's `masters/`), if given.
fn full_directory() -> Option<std::path::PathBuf> {
    std::env::var_os("VAKINT_TEST_P_MASTERS_DIRECTORY").map(std::path::PathBuf::from)
}

/// The complete data if given, else the fixture.
fn basis() -> PMasterBasisData {
    PMasterBasisData::from_directory(&full_directory().unwrap_or_else(fixture_directory)).unwrap()
}

/// The fixture plus a second layer of ours (as `our_exchanges` of helper VH 2
/// will add): `ggppss_new(1) = (2 + eps) our(1) + 1/(1 - eps) ggppss_int(1)`,
/// with `our(1)` a new master whose orders are all unknown. The first layer
/// is kept byte-identical (its data version in `layer_versions`); the masters
/// file and the basis version are new. Written under the scratch directory.
fn two_layer_directory() -> std::path::PathBuf {
    use std::io::{Read, Write};
    let gz = |value: &serde_json::Value| -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(value.to_string().as_bytes()).unwrap();
        encoder.finish().unwrap()
    };
    let ungz = |bytes: &[u8]| -> serde_json::Value {
        let mut text = String::new();
        flate2::read::GzDecoder::new(bytes)
            .read_to_string(&mut text)
            .unwrap();
        serde_json::from_str(&text).unwrap()
    };
    let source = fixture_directory();
    let directory = std::env::temp_dir().join(format!(
        "vakint-master-basis-two-layers-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(directory.join("basis-change")).unwrap();
    std::fs::create_dir_all(directory.join("pole-bounded")).unwrap();
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(source.join("MANIFEST.json")).unwrap()).unwrap();
    let first_version = manifest["basis_version"].as_str().unwrap().to_owned();
    let version = format!("{first_version}+ours-test-1");
    let first_layer = "basis-change/00-ggppss-2104.08272-new-basis.json.gz";
    std::fs::copy(source.join(first_layer), directory.join(first_layer)).unwrap();
    let masters_path = manifest["files"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["role"] == "pole-bounded-masters")
        .unwrap()["path"]
        .as_str()
        .unwrap()
        .to_owned();
    let mut masters = ungz(&std::fs::read(source.join(&masters_path)).unwrap());
    masters["data_version"] = serde_json::json!(version);
    masters["masters"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({
            "id": "our(1)", "symbol": "p_our_1", "loops": 5,
            "graph": "test master of a second layer", "leading_order": -5, "coefficients": {}
        }));
    let second = serde_json::json!({
        "schema": P_MASTER_BASIS_SCHEMA, "version": P_MASTER_BASIS_VERSION,
        "data_version": version, "variable": "eps", "dimension": "d = 4 - 2*eps",
        "normalisation": "g-scheme (as the p master catalogue)",
        "layer": {
            "name": "01-ours-test", "source": "test", "direction": "old = sum of terms",
            "new_masters": ["our(1)"],
            "exchanges": [{"old": "ggppss_new(1)", "terms": [
                {"master": "our(1)", "num": ["2", "1"], "den": ["1"]},
                {"master": "ggppss_int(1)", "num": ["1"], "den": ["1", "-1"]}
            ]}]
        }
    });
    let mut files = Vec::new();
    let mut put = |path: &str, bytes: Vec<u8>, role: &str| {
        std::fs::write(directory.join(path), &bytes).unwrap();
        files.push(serde_json::json!({
            "path": path, "role": role,
            "sha256": format!("{:x}", Sha256::digest(&bytes)), "bytes": bytes.len()
        }));
    };
    put(
        first_layer,
        std::fs::read(source.join(first_layer)).unwrap(),
        "basis-change",
    );
    put(
        "basis-change/01-ours-test.json.gz",
        gz(&second),
        "basis-change",
    );
    let new_masters_path = format!("pole-bounded/masters.{version}.json.gz");
    put(&new_masters_path, gz(&masters), "pole-bounded-masters");
    manifest["basis_version"] = serde_json::json!(version);
    manifest["layers"] = serde_json::json!(["00-ggppss-2104.08272-new-basis", "01-ours-test"]);
    manifest["layer_versions"] =
        serde_json::json!({ "00-ggppss-2104.08272-new-basis": first_version });
    manifest["files"] = serde_json::Value::Array(files);
    std::fs::write(directory.join("MANIFEST.json"), manifest.to_string()).unwrap();
    directory
}

#[test]
fn a_second_layer_of_ours_applies_after_the_first_and_round_trips() {
    let directory = two_layer_directory();
    let basis = PMasterBasisData::from_directory(&directory).unwrap();
    assert_eq!(
        basis.record().layers,
        vec!["00-ggppss-2104.08272-new-basis", "01-ours-test"]
    );
    assert!(basis.record().data_version.ends_with("+ours-test-1"));
    let catalogue = basis.extended_catalogue(&catalogue()).unwrap();
    let ep = Atom::var(epsilon());
    let combination = master(&catalogue, "ggppss_int(137)") * (&ep + Atom::num(3))
        + master(&catalogue, "ggppss_int(2)") / (Atom::num(1) + &ep);
    let rotated = basis.rotate(&combination, &catalogue, epsilon()).unwrap();
    let left = by_master(&rotated);
    for gone in ["ggppss_int(137)", "ggppss_new(1)"] {
        let symbol = catalogue.symbol_of_id(gone).unwrap();
        assert!(
            !left.keys().any(|m| m.contains(&format!("{symbol})"))),
            "{gone} survived the rotation"
        );
    }
    assert!(
        left.keys().any(|m| m.contains("p_our_1)")),
        "our(1) missing"
    );
    // undoing the second layer gives the first layer's rotation (the
    // fixture's first layer is a subset, not invertible on its own)
    let first = PMasterBasisData::from_directory(&fixture_directory()).unwrap();
    let first_rotated = first.rotate(&combination, &catalogue, epsilon()).unwrap();
    let rules = basis.rules(&catalogue, epsilon()).unwrap();
    let start = basis.layers[0].exchanges.len();
    for point in ["3/11", "-5/7"] {
        let point = vk_parse!(point).unwrap();
        let back = inverse_layer_at(
            &catalogue,
            &rules[start..],
            &basis.layers[1].new_masters,
            &at(&rotated, &point),
            &point,
        );
        assert_eq!(
            by_master(&back),
            by_master(&at(&first_rotated, &point)),
            "second layer round trip at eps = {point}"
        );
    }
    // the first layer alone keeps its own data version: a wrong one is refused
    let manifest_path = directory.join("MANIFEST.json");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&manifest_path).unwrap()).unwrap();
    manifest["layer_versions"] = serde_json::json!({});
    std::fs::write(&manifest_path, manifest.to_string()).unwrap();
    let refused = PMasterBasisData::from_directory(&directory)
        .unwrap_err()
        .to_string();
    assert!(refused.contains("do not match the manifest"), "{refused}");
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn the_masters_directory_is_bound_by_its_manifest() {
    let basis = PMasterBasisData::from_directory(&fixture_directory()).unwrap();
    let record = basis.record();
    assert_eq!(record.basis, "pole-bounded");
    assert_eq!(
        record.layers,
        vec!["00-ggppss-2104.08272-new-basis".to_owned()]
    );
    assert_eq!(record.sha256.len(), 64);
    // a changed byte in a layer file, or a missing manifest, is refused
    let copy = std::env::temp_dir().join(format!("vakint-master-basis-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&copy);
    for sub in ["basis-change", "pole-bounded"] {
        std::fs::create_dir_all(copy.join(sub)).unwrap();
        for entry in std::fs::read_dir(fixture_directory().join(sub)).unwrap() {
            let entry = entry.unwrap();
            std::fs::copy(entry.path(), copy.join(sub).join(entry.file_name())).unwrap();
        }
    }
    let missing = PMasterBasisData::from_directory(&copy)
        .unwrap_err()
        .to_string();
    assert!(missing.contains("MANIFEST.json"), "{missing}");
    std::fs::copy(
        fixture_directory().join("MANIFEST.json"),
        copy.join("MANIFEST.json"),
    )
    .unwrap();
    let layer = copy.join("basis-change/00-ggppss-2104.08272-new-basis.json.gz");
    let mut bytes = std::fs::read(&layer).unwrap();
    let last = bytes.len() - 1;
    bytes[last] ^= 1;
    std::fs::write(&layer, bytes).unwrap();
    let error = PMasterBasisData::from_directory(&copy)
        .unwrap_err()
        .to_string();
    assert!(error.contains("sha256"), "{error}");
    let _ = std::fs::remove_dir_all(&copy);
    // no package and no variable: an error naming the expected path
    if std::env::var_os(P_MASTERS_DIRECTORY_VARIABLE).is_none() {
        let none = load(None).unwrap_err().to_string();
        assert!(none.contains("masters/MANIFEST.json"), "{none}");
    }
}

#[test]
fn a_layer_never_feeds_its_own_exchanges_and_names_known_masters() {
    let basis = basis();
    let basis = &basis;
    let catalogue = basis.extended_catalogue(&catalogue()).unwrap();
    // masters exchanged by an earlier layer never come back in a later one
    let mut gone: BTreeSet<&str> = BTreeSet::new();
    for layer in &basis.layers {
        let (name, exchanges) = (&layer.name, &layer.exchanges);
        let old: BTreeSet<&str> = exchanges.iter().map(|e| e.old.as_str()).collect();
        assert_eq!(
            old.len(),
            exchanges.len(),
            "{name}: an old master exchanged twice"
        );
        for exchange in exchanges {
            assert!(
                catalogue.symbol_of_id(&exchange.old).is_some(),
                "{}",
                exchange.old
            );
            for term in &exchange.terms {
                let (num, den) = (&term.num, &term.den);
                let Some(master) = &term.master else {
                    assert!(!num.is_empty() && !den.is_empty());
                    continue;
                };
                assert!(
                    !old.contains(master.as_str()),
                    "{name}: {} feeds {master}",
                    exchange.old
                );
                assert!(catalogue.symbol_of_id(master).is_some(), "{master} unknown");
                assert!(
                    !gone.contains(master.as_str()),
                    "{name}: {master} was exchanged by an earlier layer"
                );
                assert!(!num.is_empty() && !den.is_empty());
            }
        }
        // a layer's new masters are its square block's columns
        let used: BTreeSet<&str> = exchanges
            .iter()
            .flat_map(|e| e.terms.iter().filter_map(|t| t.master.as_deref()))
            .filter(|m| layer.new_masters.iter().any(|n| n == m))
            .collect();
        assert_eq!(used.len(), layer.new_masters.len(), "{name}: new masters");
        if full_directory().is_some() || layer.name.starts_with("01-") {
            // the in-repo fixture's first layer is a subset (one exchange)
            assert_eq!(used.len(), exchanges.len(), "{name}: B is not square");
        }
        gone.extend(old);
    }
    if full_directory().is_some() {
        assert_eq!(basis.layers[0].exchanges.len(), 31);
    }
}

/// The new masters' values (exact through their known order) reproduce the
/// exchanges `old = B new + A rest` order by order through eps^0 (every
/// order below the first unknown one vanishes exactly, and no unknown order
/// enters before eps^1): the values after the eps expansion are the same in
/// both bases. `only`: the old masters checked (all when empty).
fn check_new_master_values(only: &[&str]) {
    let basis = basis();
    let basis = &basis;
    let catalogue = basis.extended_catalogue(&catalogue()).unwrap();
    let rules = basis.rules(&catalogue, epsilon()).unwrap();
    let source = PMasterSource::new(Some(&catalogue), "master basis test");
    let olds: Vec<&str> = basis.layers[0]
        .exchanges
        .iter()
        .map(|e| e.old.as_str())
        .collect();
    let mut checked = 0usize;
    for ((old, value), id) in rules.iter().zip(&olds) {
        if !only.is_empty() && !only.contains(id) {
            continue;
        }
        let started = std::time::Instant::now();
        let expanded = super::super::laurent::expand_with(
            &Atom::num(1),
            &(old - value),
            epsilon(),
            0,
            &|name, e| source.expansion(name, e),
        )
        .unwrap();
        let mut first_unknown = i64::MAX;
        let mut by_order = BTreeMap::new();
        for (power, coefficient) in expanded.coefficient_list::<i8>(&[Atom::var(epsilon())]) {
            let order = if power == Atom::num(1) {
                0
            } else if power == Atom::var(epsilon()) {
                1
            } else {
                let AtomView::Pow(p) = power.as_view() else {
                    panic!("{power}")
                };
                crate::get_integer_from_atom(p.get_base_exp().1).unwrap()
            };
            if crate::truncation_guard::has_unknown_tail(coefficient.as_view()) {
                first_unknown = first_unknown.min(order);
            }
            *by_order.entry(order).or_insert(Atom::Zero) += coefficient;
        }
        for (order, coefficient) in by_order {
            if order < first_unknown {
                let value = coefficient.expand();
                assert!(
                    value.is_zero(),
                    "{id}: eps^{order} of old - (B new + A rest) = {value}"
                );
            }
        }
        assert!(
            first_unknown >= 1,
            "{id}: an unknown order enters at eps^{first_unknown} (<= 0)"
        );
        eprintln!(
            "new_master_values: {id}: old - (B new + A rest) = 0 through eps^0 ({:.1} s)",
            started.elapsed().as_secs_f64()
        );
        checked += 1;
    }
    assert_eq!(
        checked,
        if only.is_empty() {
            olds.len()
        } else {
            only.len()
        }
    );
}

/// [`check_new_master_values`] on three exchanges (int137, the one T201
/// meets; int116, the first; int189, T178's corner): under a minute.
#[test]
fn new_master_values_reproduce_exchanges_through_eps0() {
    // the fixture holds int137 only; the complete data three exchanges
    if full_directory().is_some() {
        check_new_master_values(&["ggppss_int(137)", "ggppss_int(116)", "ggppss_int(189)"]);
    } else {
        check_new_master_values(&["ggppss_int(137)"]);
    }
}

/// [`check_new_master_values`] on all 31 exchanges (about 20 minutes on one
/// core: run explicitly, `--run-ignored only`).
#[test]
#[ignore]
fn new_master_values_reproduce_all_31_exchanges_through_eps0() {
    assert!(
        full_directory().is_some(),
        "set VAKINT_TEST_P_MASTERS_DIRECTORY"
    );
    check_new_master_values(&[]);
}

/// Original -> pole-bounded -> original (the inverse at a rational point,
/// exact over Q) on every exchanged master alone and on a combination of all
/// 31 with rational-function coefficients and an unchanged master.
#[test]
fn rotation_round_trips_through_the_inverse_at_points() {
    if full_directory().is_none() {
        eprintln!(
            "VAKINT_TEST_P_MASTERS_DIRECTORY unset: the inverse needs all 31 exchanges; skipped"
        );
        return;
    }
    let basis = basis();
    let basis = &basis;
    let catalogue = basis.extended_catalogue(&catalogue()).unwrap();
    let ep = Atom::var(epsilon());
    let mut combination = master(&catalogue, "ggppss_int(1)") * (&ep + Atom::num(3));
    for (i, layer) in basis.layers.iter().enumerate() {
        for (j, exchange) in layer.exchanges.iter().enumerate() {
            combination += master(&catalogue, &exchange.old)
                * (Atom::num((i + j + 1) as i64) + &ep * &ep / Atom::num((j + 2) as i64));
        }
    }
    let rotated = basis.rotate(&combination, &catalogue, epsilon()).unwrap();
    assert!(
        !by_master(&rotated).keys().any(|m| basis
            .layers
            .iter()
            .flat_map(|layer| &layer.exchanges)
            .any(|e| m.contains(&format!("{})", catalogue.symbol_of_id(&e.old).unwrap())))),
        "an exchanged master survived the rotation"
    );
    for point in ["3/11", "-5/7", "2/13"] {
        let point = vk_parse!(point).unwrap();
        let back = inverse_at(basis, &catalogue, &at(&rotated, &point), &point);
        let original = at(&combination, &point);
        assert_eq!(
            by_master(&back),
            by_master(&original),
            "round trip at eps = {point}"
        );
    }
}

/// The rotation of helper VH 2's composed T201 keys (its exact Python
/// reference, `helper-share/basis/testvectors`, schema
/// helper-vh2.basis-testvectors.v1) equals the expected output at rational
/// points, and in the rotated result no master's coefficient has a pole
/// beyond the master's known depth (need = -valuation at d = 4 <= q). Files:
/// VAKINT_TEST_P_MASTER_BASIS_VECTORS (':'-separated); skipped when unset.
/// With a layer of ours loaded the reference is `expected_output_ours` (T211
/// (0,1): the 15 masters beyond their depth in the GGPPSS basis are all
/// bounded once int72 is exchanged for our(1)).
#[test]
fn rotation_equals_the_reference_on_its_test_vectors() {
    let Some(paths) = std::env::var("VAKINT_TEST_P_MASTER_BASIS_VECTORS").ok() else {
        eprintln!("VAKINT_TEST_P_MASTER_BASIS_VECTORS unset: skipped");
        return;
    };
    let basis = basis();
    let basis = &basis;
    let catalogue = basis.extended_catalogue(&catalogue()).unwrap();
    let ep = Atom::var(epsilon());
    let d = Atom::num(4) - Atom::num(2) * &ep;
    let id = |name: &str| -> Option<String> {
        if name == "CLOSED" {
            None
        } else if let Some(n) = name.strip_prefix("new[") {
            Some(format!("ggppss_new({})", n.trim_end_matches(']')))
        } else if let Some(n) = name.strip_prefix("our[") {
            Some(format!("our({})", n.trim_end_matches(']')))
        } else {
            Some(name.to_owned())
        }
    };
    // master id (None = closed part) -> exact coefficient in ep, tags kept
    let parse = |table: &serde_json::Value| -> Vec<(Option<String>, Atom)> {
        let mut out = Vec::new();
        for (name, entries) in table.as_object().unwrap() {
            let mut coefficient = Atom::Zero;
            for entry in entries.as_array().unwrap() {
                let rational = |key: &str| {
                    vk_parse!(entry[key].as_str().unwrap())
                        .unwrap()
                        .replace(Atom::var(vk_symbol!("d")).to_pattern())
                        .with(d.to_pattern())
                };
                let mut tag = Atom::num(1);
                for pair in entry["tag"].as_array().unwrap() {
                    let m = pair[0].as_i64().unwrap();
                    let e = pair[1].as_i64().unwrap();
                    tag *= super::super::p_gamma_symbol()
                        .call(&[Atom::num(m)][..])
                        .pow(Atom::num(e));
                }
                coefficient += rational("num") / rational("den") * tag;
            }
            out.push((id(name), coefficient));
        }
        out
    };
    let total = |terms: &[(Option<String>, Atom)]| -> Atom {
        let mut sum = Atom::Zero;
        for (id, coefficient) in terms {
            sum += coefficient
                * match id {
                    None => Atom::num(1),
                    Some(id) => master(&catalogue, id),
                };
        }
        sum
    };
    let mut keys_checked = 0;
    for path in paths.split(':').filter(|p| !p.is_empty()) {
        let vectors: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        for key in vectors["keys"].as_array().unwrap() {
            let original = total(&parse(&key["input"]));
            // with a layer of ours: the reference with our exchanges applied
            let field = if basis.layers.len() > 1 {
                "expected_output_ours"
            } else {
                "expected_output"
            };
            assert!(key.get(field).is_some(), "{path}: no {field}");
            let expected_terms = parse(&key[field]);
            let expected = total(&expected_terms);
            let rotated = basis.rotate(&original, &catalogue, epsilon()).unwrap();
            for point in ["3/11", "-5/7"] {
                let point = vk_parse!(point).unwrap();
                let rotated_at = at(&rotated, &point);
                assert_eq!(
                    by_master(&rotated_at),
                    by_master(&at(&expected, &point)),
                    "{path} key {}",
                    key["idx"]
                );
                // round trip on the real result: back to the original basis
                // (the inverse needs all 31 exchanges)
                if full_directory().is_some() {
                    assert_eq!(
                        by_master(&inverse_at(basis, &catalogue, &rotated_at, &point)),
                        by_master(&at(&original, &point)),
                        "{path} key {}: round trip",
                        key["idx"]
                    );
                }
            }
            // pole bound on the rotated result, per master
            for (id, coefficient) in by_master_ids(&rotated, &catalogue) {
                let Some(id) = id else { continue };
                let rational = coefficient
                    .replace(
                        super::super::p_gamma_symbol()
                            .call(&[Atom::var(crate::symbols::S.x_)][..])
                            .to_pattern(),
                    )
                    .with(Atom::num(1).to_pattern());
                // an exactly cancelled coefficient needs no order (key K* of
                // our(1) itself: every other coefficient cancels)
                let Some(valuation) = valuation(&rational) else {
                    continue;
                };
                let symbol = catalogue.symbol_of_id(&id).unwrap();
                let q = catalogue.master(symbol).unwrap().known_through;
                assert!(
                    -valuation <= q,
                    "{path} key {}: {id} needs eps^{} beyond its known depth eps^{q}",
                    key["idx"],
                    -valuation
                );
            }
            keys_checked += 1;
        }
    }
    assert!(keys_checked > 0);
    eprintln!("rotation_equals_the_reference_on_its_test_vectors: {keys_checked} keys");
}

/// The exact valuation in ε of a rational function of ε (its numerator's and
/// denominator's lowest powers), `None` for 0.
fn valuation(rational: &Atom) -> Option<i64> {
    use symbolica::domains::integer::Z;
    use symbolica::domains::rational::Q;
    let r: symbolica::domains::rational_polynomial::RationalPolynomial<_, u16> =
        rational.to_rational_polynomial(&Q, &Z, None);
    if r.numerator.is_zero() {
        return None;
    }
    let lowest = |p: &symbolica::poly::polynomial::MultivariatePolynomial<_, u16>| -> i64 {
        let n = p.variables().len();
        assert!(n <= 1, "valuation: more than one variable");
        if n == 0 {
            0
        } else {
            *p.exponents.iter().min().unwrap() as i64
        }
    };
    Some(lowest(&r.numerator) - lowest(&r.denominator))
}

/// `(master id or None for the closed part, coefficient)` of `expression`.
fn by_master_ids(expression: &Atom, catalogue: &PMasterCatalogue) -> Vec<(Option<String>, Atom)> {
    let mut out: BTreeMap<Option<String>, Atom> = BTreeMap::new();
    for (monomial, coefficient) in masters_of(expression) {
        let id = match monomial.as_view() {
            AtomView::Fun(function) if function.get_symbol() == p_master_symbol() => {
                let name = function.iter().next().unwrap().to_owned();
                let AtomView::Var(var) = name.as_view() else {
                    unreachable!()
                };
                let symbol = var.get_symbol().get_stripped_name().to_owned();
                Some(
                    catalogue
                        .masters_ids()
                        .find(|(s, _)| *s == symbol)
                        .map(|(_, id)| id.to_owned())
                        .unwrap_or(symbol),
                )
            }
            _ => None,
        };
        *out.entry(id).or_insert(Atom::Zero) += coefficient;
    }
    out.into_iter().collect()
}
