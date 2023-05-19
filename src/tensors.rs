//! This file contains definitions for the most important tensors used
//! in GR.

use std::f64::consts::TAU;

use lin_alg2::f64::Mat4;

use crate::{C, G, Arr4dMetric, C_SQ};

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

/// Numerical components of a vector that may be contravariant, or covaraint.
#[derive(Clone, Copy)]
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

/// 4-vector in Minkowski coordinates. Capable of representations in covariant and
/// contravariant forms. Invariant.
#[derive(Clone, Copy)]
pub struct Vec4Minkowski {
    /// We store the numerical representation internally in its contravariant form (Upper index).
    /// We use the `as_lower` method to get its value in covariant form.
    pub value_upper: Vec4,
}

impl Vec4Minkowski {
    pub fn new(t: f64, x: f64, y: f64, z: f64) -> Self {
        Self { value_upper: Vec4 { t, x, y, z } }
    }

    /// Get the vector's contravariant (upper-index) numerical form.
    pub fn as_upper(&self) -> Vec4 {
        self.value_upper
    }

    /// Get the vector's covariant (lower-index) numerical form.
    /// V_m = g_{m n} V^n
    pub fn as_lower(&self, g: &MetricTensor) -> Vec4 {
        let V = &self.value_upper;

        let l = Tensor2Config::Ll;

        // m = t
        let t = g.val(C::T, C::T, l) * V.t
            + g.val(C::T, C::X, l) * V.x
            + g.val(C::T, C::Y, l) * V.y
            + g.val(C::T, C::Z, l) * V.z;

        // m = x
        let x = g.val(C::X, C::T, l) * V.t
            + g.val(C::X, C::X, l) * V.x
            + g.val(C::X, C::Y, l) * V.y
            + g.val(C::X, C::Z, l) * V.z;

        // m = y
        let y = g.val(C::Y, C::T, l) * V.t
            + g.val(C::Y, C::X, l) * V.x
            + g.val(C::Y, C::Y, l) * V.y
            + g.val(C::Y, C::Z, l) * V.z;

        // m = z
        let z = g.val(C::Z, C::T, l) * V.t
            + g.val(C::Z, C::X, l) * V.x
            + g.val(C::Z, C::Y, l) * V.y
            + g.val(C::Z, C::Z, l) * V.z;

        Vec4 { t, x, y, z }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        let this = self.as_upper();
        let other = other.as_upper();

        -(this.t * other.t) + this.x * other.x + this.y * other.y + this.z * other.z
    }

    /// Perform a coordinate-system transform of the upper-indexed (contravariant) form.
    /// todo: Should these var names be `ds`, or `dx`? A convention question.
    pub fn transform_upper(&self, dx: Vec4, dx_p: Vec4) -> Vec4 {
        let A = self.as_upper();

        let t =
            dx_p.t / dx.t * A.t + dx_p.t / dx.x * A.x + dx_p.t / dx.y * A.y + dx_p.t / dx.z * A.z;
        let x =
            dx_p.x / dx.t * A.t + dx_p.x / dx.x * A.x + dx_p.x / dx.y * A.y + dx_p.x / dx.z * A.z;
        let y =
            dx_p.y / dx.t * A.t + dx_p.y / dx.x * A.x + dx_p.y / dx.y * A.y + dx_p.y / dx.z * A.z;
        let z =
            dx_p.z / dx.t * A.t + dx_p.z / dx.x * A.x + dx_p.z / dx.y * A.y + dx_p.z / dx.z * A.z;

        Vec4 { t, x, y, z }
    }

    pub fn transform_lower(&self, metric: &MetricTensor, dx: Vec4, dx_p: Vec4) -> Vec4 {
        let A = self.as_lower(metric);

        let t =
            dx.t / dx_p.t * A.t + dx.x / dx_p.t * A.x + dx.y / dx_p.t * A.y + dx.z / dx_p.t * A.z;
        let x =
            dx.t / dx_p.x * A.t + dx.x / dx_p.x * A.x + dx.y / dx_p.x * A.y + dx.z / dx_p.x * A.z;
        let y =
            dx.t / dx_p.y * A.t + dx.x / dx_p.y * A.x + dx.y / dx_p.y * A.y + dx.z / dx_p.y * A.z;
        let z =
            dx.t / dx_p.z * A.t + dx.x / dx_p.z * A.x + dx.y / dx_p.z * A.y + dx.z / dx_p.z * A.z;

        Vec4 { t, x, y, z }
    }

    pub fn mag_sq(&self) -> f64 {
        -(crate::C_SQ * self.t.powi(2)) + self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
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

/// https://github.com/einsteinpy/einsteinpy/blob/main/src/einsteinpy/metric/base_metric.py
/// A metric tensor of order (2? 3?)
/// An example: https://docs.einsteinpy.org/en/stable/_modules/einsteinpy/symbolic/metric.html#MetricTensor
///
/// See the accepted answer here for a nice explanation: https://math.stackexchange.com/questions/1410793/matrix-representations-of-tensors
#[derive(Clone)]
pub struct MetricTensor {
    // /// Configuration of indices
    // config: Tensor2Config,
    /// Matrix representation: A 4x4 matrix of minkowski space. This is in LL config. ie g_{m n}.
    /// We get other configs using our `as_config` method.
    // We use a Mat4 here to take advantage of inverse, and multiplying by a constant.
    matrix_ll: Mat4,
    /// Ie the inverse metric; it's upper config.
    matrix_uu: Mat4,
}

/// Eg initialization of a metric field
impl Default for MetricTensor {
    fn default() -> Self {
        Self {
            matrix_ll: Mat4::new_identity(),
            matrix_uu: Mat4::new_identity(),
        }
    }
}

impl MetricTensor {
    /// Compute the metric tensor from coordinate bases. Uses the defintion `g_{m n} = e_m · e_n
    /// These bases are all lower-index.
    pub fn from_coordinate_bases(e_t: Vec4, e_x: Vec4, e_y: Vec4, e_z: Vec4) -> Self {
        // This matrix array layout is column-major, and lower-indexed.
        let mut g = Self {
            matrix_ll: Mat4 { data: [0.; 16] },
            matrix_uu: Mat4 { data: [0.; 16] },
        };

        let l = Tensor2Config::Ll;

        // For g_{m n}, m = t
        *g.val_mut(C::T, C::T, l) = e_t.dot(e_t); // n = t
        *g.val_mut(C::T, C::X, l) = e_t.dot(e_x); // n = x
        *g.val_mut(C::T, C::Y, l) = e_t.dot(e_y); // n = y
        *g.val_mut(C::T, C::Z, l) = e_t.dot(e_z); // n = z

        // For g_{m n}, m = x
        *g.val_mut(C::X, C::T, l) = e_x.dot(e_t); // n = x ( etc for the rest; follow the pattern)
        *g.val_mut(C::X, C::X, l) = e_x.dot(e_x);
        *g.val_mut(C::X, C::Y, l) = e_x.dot(e_y);
        *g.val_mut(C::X, C::Z, l) = e_x.dot(e_z);

        // For g_{m n}, m = y
        *g.val_mut(C::Y, C::T, l) = e_y.dot(e_t);
        *g.val_mut(C::Y, C::X, l) = e_y.dot(e_x);
        *g.val_mut(C::Y, C::Y, l) = e_y.dot(e_y);
        *g.val_mut(C::Y, C::Z, l) = e_y.dot(e_z);

        // For g_{m n}, m = z
        *g.val_mut(C::Z, C::T, l) = e_z.dot(e_t);
        *g.val_mut(C::Z, C::X, l) = e_z.dot(e_x);
        *g.val_mut(C::Z, C::Y, l) = e_z.dot(e_y);
        *g.val_mut(C::Z, C::Z, l) = e_z.dot(e_z);

        g.generate_inverse();

        g
    }

    /// Create a new metric tensor given Scharzchild gemoetry, ie a given distance from a
    /// given black hole (or spherical object in general?) with Scharzchild radius r_s.
    pub fn new_schwarzchild(M: f64, r: f64, θ: f64) -> Self {
        let r_s = 2. * M * G / C_SQ;

        let mut result = Self {
            matrix_ll: Mat4 { data: [0.; 16] },
            matrix_uu: Mat4 { data: [0.; 16] },
        };

        let a = 1. - r_s / r;
        let r_sq = r.powi(2);

        result.matrix_ll[0] = -a;
        result.matrix_ll[5] = 1. / a;
        result.matrix_ll[10] = r_sq;
        result.matrix_ll[15] = r_sq * (θ.sin()).powi(2);

        result.generate_inverse();

        result
    }

    /// Generate the inverse metric (upper-upper config).
    fn generate_inverse(&mut self) {
        self.matrix_uu = self.matrix_ll.inverse().unwrap();
    }

    // // todo :You need some method of getting by named index *while in a different config*; ie you have an immediate
    // // todo need to do that with a UU config now for finding christoffel symbols.
    // pub fn as_config(&self, config: Tensor2Config) -> &Mat4 {
    //     &match config {
    //         Tensor2Config::Uu => self.matrix_uu,
    //         Tensor2Config::Ul => {
    //             // let a = ETA_MINKOWSKI * self.matrix_ll;
    //             unimplemented!()
    //         }
    //         Tensor2Config::Lu => {
    //             unimplemented!()
    //         }
    //         Tensor2Config::Ll => self.matrix_ll,
    //     }
    // }

    /// Take the inner product of 2 vectors, using this metric; V·V = g_{mn} V^m V^n
    pub fn inner_product(&self, vec_a: Vec4Minkowski, vec_b: Vec4Minkowski) -> f64 {
        let a = vec_a.as_upper();
        let b = vec_b.as_upper();

        let l = Tensor2Config::Ll;

        self.val(C::T, C::T, l) * a.t * b.t
            + self.val(C::T, C::X, l) * a.t * b.x
            + self.val(C::T, C::Y, l) * a.t * b.y
            + self.val(C::T, C::Z, l) * a.t * b.z
            + self.val(C::X, C::T, l) * a.x * b.t
            + self.val(C::X, C::X, l) * a.x * b.x
            + self.val(C::X, C::Y, l) * a.x * b.y
            + self.val(C::X, C::Z, l) * a.x * b.z
            + self.val(C::Y, C::T, l) * a.y * b.t
            + self.val(C::Y, C::X, l) * a.y * b.x
            + self.val(C::Y, C::Y, l) * a.y * b.y
            + self.val(C::Y, C::Z, l) * a.y * b.z
            + self.val(C::Z, C::T, l) * a.z * b.t
            + self.val(C::Z, C::X, l) * a.z * b.x
            + self.val(C::Z, C::Y, l) * a.z * b.y
            + self.val(C::Z, C::Z, l) * a.z * b.z
    }

    /// Get a specific matrix component. This keeps code that uses this matrix consistent
    /// with conventions. Assumptions: Column-major array representation of the matrix, and
    /// the internal matrix is in lower-lower representation.
    pub fn val(&self, μ: V4Component, ν: V4Component, config: Tensor2Config) -> f64 {
        let g = &match config {
            Tensor2Config::Uu => &self.matrix_uu,
            Tensor2Config::Ul => {
                // let a = ETA_MINKOWSKI * self.matrix_ll;
                unimplemented!()
            }
            Tensor2Config::Lu => {
                unimplemented!()
            }
            Tensor2Config::Ll => &self.matrix_ll,
        }
        .data;

        match μ {
            C::T => match ν {
                C::T => g[0],
                C::X => g[1],
                C::Y => g[2],
                C::Z => g[3],
            },

            C::X => match ν {
                C::T => g[4],
                C::X => g[5],
                C::Y => g[6],
                C::Z => g[7],
            },

            C::Y => match ν {
                C::T => g[8],
                C::X => g[9],
                C::Y => g[10],
                C::Z => g[11],
            },

            C::Z => match ν {
                C::T => g[12],
                C::X => g[13],
                C::Y => g[14],
                C::Z => g[15],
            },
        }
    }

    /// Similar to `val`, but mutable.
    /// todo: DRY with `.val`
    pub fn val_mut(&mut self, μ: V4Component, ν: V4Component, config: Tensor2Config) -> &mut f64 {
        let g = &mut match config {
            Tensor2Config::Uu => &mut self.matrix_uu.data,
            Tensor2Config::Ul => {
                // let a = ETA_MINKOWSKI * self.matrix_ll;
                unimplemented!()
            }
            Tensor2Config::Lu => {
                unimplemented!()
            }
            Tensor2Config::Ll => &mut self.matrix_ll.data,
        };

        &mut match μ {
            C::T => match ν {
                C::T => g[0],
                C::X => g[1],
                C::Y => g[2],
                C::Z => g[3],
            },

            C::X => match ν {
                C::T => g[4],
                C::X => g[5],
                C::Y => g[6],
                C::Z => g[7],
            },

            C::Y => match ν {
                C::T => g[8],
                C::X => g[9],
                C::Y => g[10],
                C::Z => g[11],
            },

            C::Z => match ν {
                C::T => g[12],
                C::X => g[13],
                C::Y => g[14],
                C::Z => g[15],
            },
        }
    }
}

/// Group that includes the metric at a point, and at points surrounding it, an infinetesimal difference
/// in both directions along each spacial axis.
/// This is useful for when we have a Metric tensor we can find analytically; allows for accurate
/// derivatives.
#[derive(Clone)]
pub struct MetricWDiffs {
    pub on_pt: Arr4dMetric,
    pub t_prev: Arr4dMetric,
    pub t_next: Arr4dMetric,
    pub x_prev: Arr4dMetric,
    pub x_next: Arr4dMetric,
    pub y_prev: Arr4dMetric,
    pub y_next: Arr4dMetric,
    pub z_prev: Arr4dMetric,
    pub z_next: Arr4dMetric,
}

impl MetricWDiffs {
    pub fn new(grid_n: usize) -> Self {
        let data = crate::new_data_metric(grid_n);

        Self{
            on_pt: data.clone(),
            t_prev: data.clone(),
            t_next: data.clone(),
            x_prev: data.clone(),
            x_next: data.clone(),
            y_prev: data.clone(),
            y_next: data.clone(),
            z_prev: data.clone(),
            z_next: data.clone(),
        }
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
