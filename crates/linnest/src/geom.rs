use cgmath::{Angle, InnerSpace, Point2, Rad, Vector2};
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum GeomError {
    #[error("degenerate triangle: two points coincide")]
    CoincidentPoints,
    #[error("degenerate triangle: points are collinear")]
    Collinear,
    #[error("Not computed")]
    NotComputed,
}

#[inline]
fn cross_z(u: Vector2<f64>, v: Vector2<f64>) -> f64 {
    u.x * v.y - u.y * v.x
}

/// Oriented angle from AB to the tangent ray at A that lies on C's side of AB.
/// Returns `Rad<f64>` in (-π, π]. Positive = CCW from AB.
/// This equals `sign(orient(A,B,C)) * (π - ∠ACB)`.
pub fn tangent_angle_toward_c_side(
    a: Point2<f64>,
    b: Point2<f64>,
    c: Point2<f64>,
) -> Result<Rad<f64>, GeomError> {
    // Exact coincidence check
    if a == b || a == c || b == c {
        return Err(GeomError::CoincidentPoints);
    }

    let ab = b - a;
    let ac = c - a;

    // Collinearity check with relative tolerance
    const COL_EPS: f64 = 1e-14;
    let area2 = cross_z(ab, ac);
    let ab_len = ab.magnitude();
    let ac_len = ac.magnitude();
    if area2.abs() <= COL_EPS * (ab_len * ac_len) {
        return Err(GeomError::Collinear);
    }

    // cgmath's angle() returns a signed angle that already encodes the orientation
    // of the triangle, so we don't need to compute orientation separately
    let ca = a - c;
    let cb = b - c;
    let angle_c = ca.angle(cb);

    // The tangent angle is the supplement of ∠ACB
    Ok((Rad::<f64>::turn_div_2() - angle_c).normalize_signed())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[inline]
    fn turns(t: f64) -> Rad<f64> {
        Rad::<f64>::full_turn() * t
    }

    #[test]
    fn basic_right_triangle_toward_c_side() {
        // A=(0,0), B=(1,0), C=(0,1): ∠ACB=π/4 ⇒ θ=+3π/4 (0.375 turns)
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.0, 1.0);
        let ang = tangent_angle_toward_c_side(a, b, c).unwrap();

        assert!((ang - turns(3.0 / 8.0)).0.abs() < 1e-12);
    }

    #[test]
    fn mirrored_right_triangle_toward_c_side() {
        // Reflect C below AB ⇒ θ=-3π/4
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.0, -1.0);
        let ang = tangent_angle_toward_c_side(a, b, c).unwrap();

        assert!((ang + turns(3.0 / 8.0)).0.abs() < 1e-12);
    }

    #[test]
    fn antisymmetry_when_swapping_ab() {
        // Some non-degenerate, non-right triangle
        let a = Point2::new(0.2, -0.7);
        let b = Point2::new(2.5, 0.1);
        let c = Point2::new(-0.4, 1.3);

        let theta_ab = tangent_angle_toward_c_side(a, b, c).unwrap();
        let theta_ba = tangent_angle_toward_c_side(b, a, c).unwrap();

        // Oriented angle must flip when the base ray is reversed
        assert!(
            (theta_ab.0 + theta_ba.0).abs() < 1e-12,
            "Expected θ(A→B,C) = -θ(B→A,C), got {} vs {}",
            theta_ab.0,
            theta_ba.0
        );
    }

    #[test]
    fn tangent_really_points_to_cs_side() {
        // A few diverse triangles (randomize more if you like)
        let cases = [
            (
                Point2::new(0.2, -0.7),
                Point2::new(2.5, 0.1),
                Point2::new(-0.4, 1.3),
            ),
            (
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(0.3, 0.9),
            ),
            (
                Point2::new(-1.2, 2.0),
                Point2::new(0.4, 3.1),
                Point2::new(-2.0, -0.8),
            ),
            (
                Point2::new(5.0, -3.0),
                Point2::new(8.0, 1.0),
                Point2::new(6.0, -6.0),
            ),
        ];

        for (a, b, c) in cases {
            let theta = tangent_angle_toward_c_side(a, b, c).unwrap();

            // Rotate AB by the reported angle and check the side:
            let ab = b - a;
            let (s, c0) = theta.0.sin_cos();
            let t = Vector2::new(c0 * ab.x - s * ab.y, s * ab.x + c0 * ab.y);

            let orient = cross_z(b - a, c - a);
            let side_rotated = cross_z(ab, t);
            assert!(
                orient * side_rotated > 0.0,
                "tangent not on C's side: orient={}, side_rotated={}, θ={}",
                orient,
                side_rotated,
                theta.0
            );
        }
    }

    #[test]
    fn coincident_points_error() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(0.0, 0.0);
        let c = Point2::new(1.0, 0.0);
        assert_eq!(
            tangent_angle_toward_c_side(a, b, c).unwrap_err(),
            GeomError::CoincidentPoints
        );
    }

    #[test]
    fn collinear_error() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(2.0, 0.0);
        let c = Point2::new(5e-16, 0.0);
        assert_eq!(
            tangent_angle_toward_c_side(a, b, c).unwrap_err(),
            GeomError::Collinear
        );
    }

    #[test]
    fn invariance_under_rigid_transformations() {
        // Original triangle
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.0, 1.0);
        let original_angle = tangent_angle_toward_c_side(a, b, c).unwrap();

        // Apply translation
        let translation = Vector2::new(5.0, -3.0);
        let a_translated = a + translation;
        let b_translated = b + translation;
        let c_translated = c + translation;
        let translated_angle =
            tangent_angle_toward_c_side(a_translated, b_translated, c_translated).unwrap();

        // Apply rotation (45 degrees)
        let cos_theta = (std::f64::consts::PI / 4.0).cos();
        let sin_theta = (std::f64::consts::PI / 4.0).sin();

        let rotate = |p: Point2<f64>| -> Point2<f64> {
            Point2::new(
                p.x * cos_theta - p.y * sin_theta,
                p.x * sin_theta + p.y * cos_theta,
            )
        };

        let a_rotated = rotate(a);
        let b_rotated = rotate(b);
        let c_rotated = rotate(c);
        let rotated_angle = tangent_angle_toward_c_side(a_rotated, b_rotated, c_rotated).unwrap();

        // Apply both translation and rotation
        let a_transformed = rotate(a) + translation;
        let b_transformed = rotate(b) + translation;
        let c_transformed = rotate(c) + translation;
        let transformed_angle =
            tangent_angle_toward_c_side(a_transformed, b_transformed, c_transformed).unwrap();

        // All angles should be the same (within numerical tolerance)
        assert!(
            (original_angle.0 - translated_angle.0).abs() < 1e-12,
            "Translation changed angle: {} vs {}",
            original_angle.0,
            translated_angle.0
        );
        assert!(
            (original_angle.0 - rotated_angle.0).abs() < 1e-12,
            "Rotation changed angle: {} vs {}",
            original_angle.0,
            rotated_angle.0
        );
        assert!(
            (original_angle.0 - transformed_angle.0).abs() < 1e-12,
            "Combined transformation changed angle: {} vs {}",
            original_angle.0,
            transformed_angle.0
        );
    }

    #[test]
    fn invariance_under_rigid_transformations_mirrored() {
        // Original mirrored triangle (negative angle case)
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.0, -1.0);
        let original_angle = tangent_angle_toward_c_side(a, b, c).unwrap();

        // Apply translation
        let translation = Vector2::new(-2.5, 7.3);
        let a_translated = a + translation;
        let b_translated = b + translation;
        let c_translated = c + translation;
        let translated_angle =
            tangent_angle_toward_c_side(a_translated, b_translated, c_translated).unwrap();

        // Apply rotation (120 degrees)
        let cos_theta = (2.0 * std::f64::consts::PI / 3.0).cos();
        let sin_theta = (2.0 * std::f64::consts::PI / 3.0).sin();

        let rotate = |p: Point2<f64>| -> Point2<f64> {
            Point2::new(
                p.x * cos_theta - p.y * sin_theta,
                p.x * sin_theta + p.y * cos_theta,
            )
        };

        let a_rotated = rotate(a);
        let b_rotated = rotate(b);
        let c_rotated = rotate(c);
        let rotated_angle = tangent_angle_toward_c_side(a_rotated, b_rotated, c_rotated).unwrap();

        // Apply both translation and rotation
        let a_transformed = rotate(a) + translation;
        let b_transformed = rotate(b) + translation;
        let c_transformed = rotate(c) + translation;
        let transformed_angle =
            tangent_angle_toward_c_side(a_transformed, b_transformed, c_transformed).unwrap();

        // All angles should be the same (within numerical tolerance)
        assert!(
            (original_angle.0 - translated_angle.0).abs() < 1e-12,
            "Translation changed mirrored angle: {} vs {}",
            original_angle.0,
            translated_angle.0
        );
        assert!(
            (original_angle.0 - rotated_angle.0).abs() < 1e-12,
            "Rotation changed mirrored angle: {} vs {}",
            original_angle.0,
            rotated_angle.0
        );
        assert!(
            (original_angle.0 - transformed_angle.0).abs() < 1e-12,
            "Combined transformation changed mirrored angle: {} vs {}",
            original_angle.0,
            transformed_angle.0
        );
    }
}
