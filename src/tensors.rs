//! This file contains definitions for the most important tensors used
//! in GR.

use std::f64::consts::TAU;

use lin_alg2::f64::Mat4;

use crate::christoffel::Christoffel;
use crate::{
    metric::MetricTensor, Arr4dChristoffel, Arr4dMetric, Arr4dReal, Worldline, C, C_SQ, G,
};

pub const COMPS: [V4Component; 4] = [
    V4Component::T,
    V4Component::X,
    V4Component::y,
    V4Component::Z,
];

#[rustfmt::skip]
const ETA_MINKOWSKI: Mat4 = Mat4 { data: [
    -1., 0., 0., 0.,
    0., 1., 0., 0.,
    0., 0., 1., 0.,
    0., 0., 0., 1.
]};

#[derive(Clone, Copy)]
pub enum Tensor1Config {
    U,
    L,
}

#[derive(Clone, Copy)]
pub enum Tensor2Config {
    Uu,
    Ul,
    Lu,
    Ll,
}

// todo: deprecate?
/// Numerical components of a vector that may be contravariant, or covaraint.
#[derive(Clone, Copy, Default)]
pub struct Vec4 {
    pub t: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec4 {
    pub fn dot(self, other: Self) -> f64 {
        -(self.t * other.t) + self.x * other.x + self.y * other.y + self.z * other.z
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
/// Represents a component of a 4-vector in Minkowski space. Compact syntax sacrifices explicitness
/// for code readability since it appears in groups.
pub enum V4Component {
    T,
    X,
    Y,
    Z,
}

// Our equations make heavy use of this; keep things terse.
pub type C = V4Component;

// /// Group Vec4s, for calculating partial derivatives
// struct Vec4Set {
//     on_pt: Vec4Minkowski
// }

/// 4-vector in Minkowski coordinates. Capable of representations in covariant and
/// contravariant forms. Invariant.
#[derive(Clone, Copy, Default)]
pub struct Vec4Minkowski {
    /// We store the numerical representation internally in its contravariant form (Upper index).
    /// We use the `as_lower` method to get its value in covariant form.
    pub components_u: [f64; 4],
}

impl Vec4Minkowski {
    pub fn new(t: f64, x: f64, y: f64, z: f64) -> Self {
        // Self { value_upper: Vec4 { t, x, y, z } }
        Self {
            components_u: [t, x, y, z],
        }
    }

    // /// Get the vector's contravariant (upper-index) numerical form.
    // pub fn as_upper(&self) -> Vec4 {
    //     self.value_upper
    // }
    //
    // /// Get the vector's covariant (lower-index) numerical form.
    // /// V_μ = g_{μ ν} V^ν
    // pub fn as_lower(&self, g: &MetricTensor) -> Vec4 {
    //     let mut result = Self::default();
    //
    //     for μ in &COMPS {
    //         for ν in &COMPS {
    //             *result.val_mut(*μ) = g.val(*μ, *ν, Tensor2Config::Ll) * self.val(*ν, Tensor1Config::U);
    //         }
    //     }
    //
    //     result
    // }

    // pub fn dot(&self, other: &Self) -> f64 {
    //     let this = self.as_upper();
    //     let other = other.as_upper();
    //
    //     -(this.t * other.t) + this.x * other.x + this.y * other.y + this.z * other.z
    // }

    // // todo
    //
    // /// Perform a coordinate-system transform of the upper-indexed (contravariant) form.
    // /// todo: Should these var names be `ds`, or `dx`? A convention question.
    // pub fn transform_upper(&self, dx: Vec4, dx_p: Vec4) -> Vec4 {
    //     let A = self.as_upper();
    //     let t =
    //         dx_p.t / dx.t * A.t + dx_p.t / dx.x * A.x + dx_p.t / dx.y * A.y + dx_p.t / dx.z * A.z;
    //     let x =
    //         dx_p.x / dx.t * A.t + dx_p.x / dx.x * A.x + dx_p.x / dx.y * A.y + dx_p.x / dx.z * A.z;
    //     let y =
    //         dx_p.y / dx.t * A.t + dx_p.y / dx.x * A.x + dx_p.y / dx.y * A.y + dx_p.y / dx.z * A.z;
    //     let z =
    //         dx_p.z / dx.t * A.t + dx_p.z / dx.x * A.x + dx_p.z / dx.y * A.y + dx_p.z / dx.z * A.z;
    //
    //     Vec4 { t, x, y, z }
    // }
    //
    // // todo
    //
    // pub fn transform_lower(&self, metric: &MetricTensor, dx: Vec4, dx_p: Vec4) -> Vec4 {
    //     let A = self.as_lower(metric);
    //
    //     let t =
    //         dx.t / dx_p.t * A.t + dx.x / dx_p.t * A.x + dx.y / dx_p.t * A.y + dx.z / dx_p.t * A.z;
    //     let x =
    //         dx.t / dx_p.x * A.t + dx.x / dx_p.x * A.x + dx.y / dx_p.x * A.y + dx.z / dx_p.x * A.z;
    //     let y =
    //         dx.t / dx_p.y * A.t + dx.x / dx_p.y * A.x + dx.y / dx_p.y * A.y + dx.z / dx_p.y * A.z;
    //     let z =
    //         dx.t / dx_p.z * A.t + dx.x / dx_p.z * A.x + dx.y / dx_p.z * A.y + dx.z / dx_p.z * A.z;
    //
    //     Vec4 { t, x, y, z }
    // }

    pub fn mag_sq(&self) -> f64 {
        -(C_SQ * self.t.powi(2)) + self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    // todo: Sort out how you get covariant config.

    // pub fn val(&self, comp: V4Component, config: Tensor1Config) -> f64 {
    pub fn val(&self, comp: V4Component) -> f64 {
        match comp {
            C::T => self.components_u[0],
            C::X => self.components_u[1],
            C::Y => self.components_u[2],
            C::Z => self.components_u[3],
        }
    }

    // pub fn val_mut(&mut self, comp: V4Component, config: Tensor1Config) -> &mut f64 {
    pub fn val_mut(&mut self, comp: V4Component) -> &mut f64 {
        &mut match comp {
            C::T => self.components_u[0],
            C::X => self.components_u[1],
            C::Y => self.components_u[2],
            C::Z => self.components_u[3],
        }
    }

    // todo
    // pub fn lortenz_transform(&self, v: f64) -> Self {
    //     let β = v.powi(2) / C_SQ;
    //     let γ = 1. / (1. - β).sqrt();
    //
    //     #[rustfmt::skip]
    //     let mat = Mat4::new([
    //        γ, -γ * β, 0., 0.,
    //        -γ * β, γ,  0., 0.,
    //        0., 0., 1., 0.,
    //        0., 0., 0., 1.,
    //     ]);
    //
    //     // mat.dot(&self.as_upper())
    // }
}

/// C+P from Metric Tensor, with changes A/R. Combine with other tensors like Metric A/R.
/// Can be used to generate the Einstein tensor.
pub struct StressEnergyTensor {
    // /// Configuration of indices
    // config: Tensor2Config,
    /// Matrix representation: A 4x4 matrix of minkowski space. This is in UU config. ie g^{m n}.
    /// We get other configs using our `as_config` method.
    matrix_uu: Mat4,
    // Dependent on uu.
    matrix_ll: Mat4,
}

impl StressEnergyTensor {
    // pub fn as_config(&self, config: Tensor2Config) -> &Mat4 {
    //     // todo
    //     &match config {
    //         Tensor2Config::Uu => self.matrix_uu,
    //         Tensor2Config::Ul => {
    //             // let a = ETA_MINKOWSKI * self.matrix_uu;
    //             unimplemented!()
    //         }
    //         Tensor2Config::Lu => {
    //             unimplemented!()
    //         }
    //         Tensor2Config::Ll => self.matrix_ll,
    //     }
    // }

    /// Return the matrix rep of the Einstein tensor, as a UU config.
    pub fn get_einstein_tensor(&self) -> Mat4 {
        self.matrix_uu.clone() * 4. * TAU * G
    }
}

pub struct EinsteinTensor {
    pub coords: [f64; 1],
    pub con: i8,
    pub components: [f64; 1],
    pub index_values: [i8; 1],
}

impl EinsteinTensor {
    pub fn new(ricci: &RicciTensor) {}
}

/// The path a particle takes through spacetime. Can be a geodesic.
pub struct Worldline {
    /// A list of events, with constant proper-time spacing.
    pub events: Vec<Vec4Minkowski>,
    // todo: v and accel too?
    /// Defines the spacing between events. A smaller spacing is more precise.
    /// todo: Maybe this isn't general, ie doesn't apply if mass is 0, like for photons?
    /// todo In that case, you need something else. (Proper time is 0 along a photon's worldline)
    /// todo: Is this unit meters, seconds, or something else?
    pub dτ: f64,
}

impl Default for Worldline {
    fn default() -> Self {
        Self {
            events: Vec::new(),
            dτ: 1.,
        }
    }
}

impl Worldline {
    /// Create a geodesic, given a Christoffel grid, an initial position, and initial velocity.
    /// todo: Currently uses Euler integration; improve this.
    pub fn new_geodesic(
        Γs: Arr4dChristoffel,
        grid_posits: &Arr4dReal,
        posit_init_i: &[usize; 4],
        v_init: Vec4Minkowski,
        num_events: usize,
        dτ: f64,
    ) -> Self {
        let mut result = Self {
            events: Vec::new(),
            dτ,
        };

        let mut a = Vec4Minkowski::default();
        let mut v = v_init;
        let mut s = grid_posits[posit_init_i][0][posit_init_i][1][posit_init_i][2][posit_init_i][3];

        for _ in num_events {
            result.events.push(s);

            // todo: Big problem: Your Chrisftoffels are calcualted at discrete points, so how do
            // todo you get them on the arbitrary geodesic points? Sure, you can do this if you
            // todo have an analytic metric like Schw.
            a = Self::calc_geodesic_accel(Γs[])
        }

        result
    }

    /// τ_AB = \int (0, 1) sqrt(-g_μv(x^alpha(\sigma)) dx^u/dsigma dx^v / dsigma) dsigma
    pub fn calc_proper_time(&self) -> f64 {
        0.
    }

    /// Calculate the acceleration (relative to τ) of a geodesic, at a given spacetime point (encompassed by
    /// the Christoffel symbol), and velocity (relative to τ) at that point.
    /// d^2 x^μ / dτ^2 + Γ^μ_ρσ dx^ρ /dτ dx^σ /dτ = 0
    pub fn calc_geodesic_accel(Γ: &Christoffel, v: Vec4Minkowski) -> Vec4Minkowski {
        let mut result = Vec4Minkowski::default();

        for μ in &COMPS {
            for ρ in &COMPS {
                for σ in &COMPS {
                    *result.val_mut(*μ) = -Γ.val(*μ, *ρ, *σ) * v.val(*ρ) * v.val(*σ);
                }
            }
        }

        result
    }
}
