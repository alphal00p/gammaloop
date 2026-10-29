from pathlib import Path
import shutil
src=Path('/home/lcnbr/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/kurbo-0.13.1')
shutil.copytree(src,'vendor/kurbo',ignore=shutil.ignore_patterns('.cargo-ok','.cargo_vcs_info.json'))
p=Path('Cargo.toml');s=p.read_text().replace('exclude = ["vendor/clarabel"]','exclude = ["vendor/clarabel", "vendor/kurbo"]');s=s.replace('[patch.crates-io]','[patch.crates-io]\n# Stable finite miter evaluation for native and Wasm drawing coverage.\nkurbo = { path = "vendor/kurbo" }',1);p.write_text(s)
p=Path('examples/notebooks/symbolica-host/Cargo.toml');s=p.read_text().replace('[patch.crates-io]','[patch.crates-io]\n# Use the same authoritative stroke owner as the drawing workspace.\nkurbo = { path = "../../../vendor/kurbo" }',1);p.write_text(s)
for filename in ['Cargo.lock','examples/notebooks/symbolica-host/Cargo.lock']:
 p=Path(filename);s=p.read_text();s=s.replace('name = "kurbo"\nversion = "0.13.1"\nsource = "registry+https://github.com/rust-lang/crates.io-index"\nchecksum = "4b60dfc32f652b926df6192e55525b16d186c69d47876c3ead4da5cc9f8450e2"','name = "kurbo"\nversion = "0.13.1"');p.write_text(s)
p=Path('vendor/kurbo/src/stroke.rs');s=p.read_text();start=s.index('    fn do_join(')
s=s[:start]+'''    // A power-of-two rescaling retains direction bits without underflowing
    // cross products or overflowing dot products of finite tangents.
    fn scaled_tangent(tangent: Vec2) -> Vec2 {
        let largest = tangent.x.abs().max(tangent.y.abs());
        if largest == 0.0 {
            return tangent;
        }
        let bits = largest.to_bits();
        let exponent = (bits >> 52) & 0x7ff;
        let scale = if exponent == 0 {
            f64::from_bits(1_u64 << (63 - bits.leading_zeros()))
        } else {
            f64::from_bits(exponent << 52)
        };
        Vec2::new(tangent.x / scale, tangent.y / scale)
    }

    // Preserve the small turning angle when the two products nearly cancel.
    fn tangent_cross(a: Vec2, b: Vec2) -> f64 {
        let product = a.y * b.x;
        a.x.mul_add(b.y, -product) - a.y.mul_add(b.x, -product)
    }

'''+s[start:]
s=s.replace('let scale = 0.5 * style.width / tan0.hypot();\n        let norm = scale * Vec2::new(-tan0.y, tan0.x);','let cd = Self::scaled_tangent(tan0);\n        let scale = 0.5 * style.width / cd.hypot();\n        let norm = scale * Vec2::new(-cd.y, cd.x);',1)
s=s.replace('let ab = self.last_tan;\n            let cd = tan0;\n            let cross = ab.cross(cd);','let ab = Self::scaled_tangent(self.last_tan);\n            let cross = Self::tangent_cross(ab, cd);',1)
s=s.replace('''                        if 2.0 * hypot < (hypot + dot) * style.miter_limit.powi(2) {''','''                        // Rationalize 1 + cos(theta) near a half-turn. The
                        // direct sum hypot + dot can round to zero even when
                        // the finite miter is within the requested limit.
                        let miter_ratio = if dot < 0.0 {
                            (2.0 * hypot).sqrt() * (hypot - dot).sqrt() / cross.abs()
                        } else {
                            (2.0 * hypot / (hypot + dot)).sqrt()
                        };
                        if miter_ratio < style.miter_limit {''',1)
s=s.replace('''                                let fp_last = p0 - last_norm;
                                let fp_this = p0 - norm;
                                let h = ab.cross(fp_this - fp_last) / cross;
                                let miter_pt = fp_this - cd * h;''','''                                // Intersect offsets in the vertex's local
                                // frame before translating: subtracting two
                                // already translated points loses the offset
                                // at large coordinates and deep notches.
                                let h = Self::tangent_cross(ab, last_norm - norm) / cross;
                                let miter_pt = p0 + (-norm - cd * h);''',1)
s=s.replace('''                                let fp_last = p0 + last_norm;
                                let fp_this = p0 + norm;
                                let h = ab.cross(fp_this - fp_last) / cross;
                                let miter_pt = fp_this - cd * h;''','''                                let h = Self::tangent_cross(ab, norm - last_norm) / cross;
                                let miter_pt = p0 + (norm - cd * h);''',1)
p.write_text(s)
