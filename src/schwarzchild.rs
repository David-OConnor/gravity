//! Code related to the Scharzchild geometry. Note that some functionality is in methods of
//! `Worldline` and `MetricTensor`.

use crate::tensors::Vec4Minkowski;

use lin_alg2::f64::Vec3;

/// Helper fn for generating neighboring points for the Schwarzchild metric
/// This is a conversion from cartesian to spherical coordinates.
pub fn find_params(posit_sample: Vec4Minkowski, posit_mass: Vec3) -> (f64, f64) {
    let diff = Vec3::new(
        posit_sample.x() - posit_mass.x,
        posit_sample.y() - posit_mass.y,
        posit_sample.z() - posit_mass.z,
    );

    let r = diff.magnitude();

    // todo: phi or theta??
    let θ = (diff.x.powi(2) + diff.y.powi(2)).sqrt().atan2(diff.z);

    (r, θ)
}
