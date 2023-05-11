//! This file contains definitions for the most important tensors used
//! in GR.

use std::f64::consts::TAU;

use lin_alg2::f64::Mat4;

use crate::G;

#[rustfmt::skip]
const ETA_MINKOWSKI: Mat4 = Mat4 { data: [
    -1., 0., 0., 0.,
    0., 1., 0., 0.,
    0., 0., 1., 0.,
    0., 0., 0., 1.
]};

#[derive(Clone, Copy)]
enum Tensor1Config {
    U,
    L,
}

#[derive(Clone, Copy)]
enum Tensor2Config {
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

#[derie(Clone, Copy)]
enum V4Comp {
    T,
    X,
    Y,
    Z,
}

/// Christoffel symbol.
pub struct Christoffel {}

/// 4-vector in Minkowski coordinates. Capable of representations in covariant and
/// contravariant forms. Invariant.
#[derive(Clone, Copy)]
pub struct Vec4Minkowski {
    /// We store the numerical representation internally in its contravariant form (Upper index).
    /// We use the `as_lower` method to get its value in covariant form.
    pub value_upper: Vec4,
}

impl Vec4Minkowski {
    /// Get the vector's contravariant (upper-index) numerical form.
    pub fn as_upper(&self) -> Vec4 {
        self.value_upper
    }

    /// Get the vector's covariant (lower-index) numerical form.
    /// V_m = g_{m n} V^n
    pub fn as_lower(&self, g: &MetricTensor) -> Vec4 {
        let V = &self.value_upper;

        // m = t
        let t = g.val(V4Comp::T, V4Comp::T) * V.t
            + g.val(V4Comp::T, V4Comp::X) * V.x
            + g.val(V4Comp::T, V4Comp::Y) * V.y
            + g.val(V4Comp::T, V4Comp::Z) * V.z;

        let x = g.val(V4Comp::X, V4Comp::T) * V.t
            + g.val(V4Comp::T, V4Comp::X) * V.x
            + g.val(V4Comp::T, V4Comp::Y) * V.y
            + g.val(V4Comp::T, V4Comp::Z) * V.z;

        let y = g.val(V4Comp::Y, V4Comp::T) * V.t
            + g.val(V4Comp::T, V4Comp::X) * V.x
            + g.val(V4Comp::T, V4Comp::Y) * V.y
            + g.val(V4Comp::T, V4Comp::Z) * V.z;

        let z = g.val(V4Comp::Z, V4Comp::T) * V.t
            + g.val(V4Comp::T, V4Comp::X) * V.x
            + g.val(V4Comp::T, V4Comp::Y) * V.y
            + g.val(V4Comp::T, V4Comp::Z) * V.z;

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
}

impl StressEnergyTensor {
    pub fn as_config(&self, config: Tensor2Config) -> Mat4 {
        // todo
        match config {
            // todo: Unnecessary clone; find another way.
            Tensor2Config::Uu => self.matrix_uu.clone(),
            Tensor2Config::Ul => {
                let a = ETA_MINKOWSKI * self.matrix_uu;
                unimplemented!()
            }
            Tensor2Config::Lu => {
                unimplemented!()
            }
            Tensor2Config::Ll => self.matrix_uu.inverse().unwrap(),
        }
    }

    /// Return the matrix rep of the Einstein tensor, as a UU config.
    pub fn get_einstein_tensor(&self) -> Mat4 {
        self.matrix_uu * 4. * TAU * G
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
    pub matrix_ll: Mat4,
}

impl MetricTensor {
    /// Compute the metric tensor from coordinate bases. Uses the defintion `g_{m n} = e_m · e_n
    /// These bases are all lower-index.
    pub fn from_coordinate_bases(e_t: Vec4, e_x: Vec4, e_y: Vec4, e_z: Vec4) -> Self {
        // This matrix array layout is column-major, and lower-indexed.
        let mut g = Self {
            matrix_ll: Mat4 { data: [0.; 16] },
        };

        // For g_{m n}, m = t
        *g.val_mut(V4Comp::T, V4Comp::T) = e_t.dot(e_t); // n = t
        *g.val_mut(V4Comp::T, V4Comp::X) = e_t.dot(e_x); // n = x
        *g.val_mut(V4Comp::T, V4Comp::Y) = e_t.dot(e_y); // n = y
        *g.val_mut(V4Comp::T, V4Comp::Z) = e_t.dot(e_z); // n = z

        // For g_{m n}, m = x
        *g.val_mut(V4Comp::X, V4Comp::T) = e_x.dot(e_t); // n = x ( etc for the rest; follow the pattern)
        *g.val_mut(V4Comp::X, V4Comp::X) = e_x.dot(e_x);
        *g.val_mut(V4Comp::X, V4Comp::Y) = e_x.dot(e_y);
        *g.val_mut(V4Comp::X, V4Comp::Z) = e_x.dot(e_z);

        // For g_{m n}, m = y
        *g.val_mut(V4Comp::Y, V4Comp::T) = e_y.dot(e_t);
        *g.val_mut(V4Comp::Y, V4Comp::X) = e_y.dot(e_x);
        *g.val_mut(V4Comp::Y, V4Comp::Y) = e_y.dot(e_y);
        *g.val_mut(V4Comp::Y, V4Comp::Z) = e_y.dot(e_z);

        // For g_{m n}, m = z
        *g.val_mut(V4Comp::Z, V4Comp::T) = e_z.dot(e_t);
        *g.val_mut(V4Comp::Z, V4Comp::X) = e_z.dot(e_x);
        *g.val_mut(V4Comp::Z, V4Comp::Y) = e_z.dot(e_y);
        *g.val_mut(V4Comp::Z, V4Comp::Z) = e_z.dot(e_z);

        g
    }

    pub fn as_config(&self, config: Tensor2Config) -> Mat4 {
        match config {
            Tensor2Config::Uu => self.matrix_ll.inverse().unwrap(),
            Tensor2Config::Ul => {
                // let a = ETA_MINKOWSKI * self.matrix_ll;
                unimplemented!()
            }
            Tensor2Config::Lu => {
                unimplemented!()
            }
            // todo: Unnecessary clone; find another way.
            Tensor2Config::Ll => self.matrix_ll.clone(),
        }
    }

    /// Take the inner product of 2 vectors, using this metric; V·V = g_{mn} V^m V^n
    pub fn inner_product(&self, vec_a: Vec4Minkowski, vec_b: Vec4Minkowski) -> f64 {
        // todo: Maybe methods on *this* struct?
        let a = vec_a.as_upper();
        let b = vec_b.as_upper();

        self.val(V4Comp::T, V4Comp::T) * a.t * b.t
            + self.val(V4Comp::T, V4Comp::X) * a.t * b.x
            + self.val(V4Comp::T, V4Comp::Y) * a.t * b.y
            + self.val(V4Comp::T, V4Comp::Z) * a.t * b.z
            + self.val(V4Comp::X, V4Comp::T) * a.x * b.t
            + self.val(V4Comp::X, V4Comp::X) * a.x * b.x
            + self.val(V4Comp::X, V4Comp::Y) * a.x * b.y
            + self.val(V4Comp::X, V4Comp::Z) * a.x * b.z
            + self.val(V4Comp::Y, V4Comp::T) * a.y * b.t
            + self.val(V4Comp::Y, V4Comp::X) * a.y * b.x
            + self.val(V4Comp::Y, V4Comp::Y) * a.y * b.y
            + self.val(V4Comp::Y, V4Comp::Z) * a.y * b.z
            + self.val(V4Comp::Z, V4Comp::T) * a.z * b.t
            + self.val(V4Comp::Z, V4Comp::X) * a.z * b.x
            + self.val(V4Comp::Z, V4Comp::Y) * a.z * b.y
            + self.val(V4Comp::Z, V4Comp::Z) * a.z * b.z
    }

    /// Get a specific matrix component. This keeps code that uses this matrix consistent
    /// with conventions. Assumptions: Column-major array representation of the matrix, and
    /// the internal matrix is in lower-lower representation.
    pub fn val(&self, m: V4Comp, n: V4Comp) -> f64 {
        // todo: Order? Is this reversed?
        match m {
            V4Comp::T => match n {
                V4Comp::T => self.matrix_ll.data[0],
                V4Comp::X => self.matrix_ll.data[1],
                V4Comp::Y => self.matrix_ll.data[2],
                V4Comp::Z => self.matrix_ll.data[3],
            },

            V4Comp::X => match n {
                V4Comp::T => self.matrix_ll.data[4],
                V4Comp::X => self.matrix_ll.data[5],
                V4Comp::Y => self.matrix_ll.data[6],
                V4Comp::Z => self.matrix_ll.data[7],
            },

            V4Comp::Y => match n {
                V4Comp::T => self.matrix_ll.data[8],
                V4Comp::X => self.matrix_ll.data[9],
                V4Comp::Y => self.matrix_ll.data[10],
                V4Comp::Z => self.matrix_ll.data[11],
            },

            V4Comp::Z => match n {
                V4Comp::T => self.matrix_ll.data[12],
                V4Comp::X => self.matrix_ll.data[13],
                V4Comp::Y => self.matrix_ll.data[14],
                V4Comp::Z => self.matrix_ll.data[15],
            },
        }
    }

    /// Similar to `val`, but mutable.
    /// todo: DRY with `.val`
    pub fn val_mut(&mut self, m: V4Comp, n: V4Comp) -> &mut f64 {
        // todo: Same order caveat as with Val.
        &mut match m {
            V4Comp::T => match n {
                V4Comp::T => self.matrix_ll.data[0],
                V4Comp::X => self.matrix_ll.data[1],
                V4Comp::Y => self.matrix_ll.data[2],
                V4Comp::Z => self.matrix_ll.data[3],
            },

            V4Comp::X => match n {
                V4Comp::T => self.matrix_ll.data[4],
                V4Comp::X => self.matrix_ll.data[5],
                V4Comp::Y => self.matrix_ll.data[6],
                V4Comp::Z => self.matrix_ll.data[7],
            },

            V4Comp::Y => match n {
                V4Comp::T => self.matrix_ll.data[8],
                V4Comp::X => self.matrix_ll.data[9],
                V4Comp::Y => self.matrix_ll.data[10],
                V4Comp::Z => self.matrix_ll.data[11],
            },

            V4Comp::Z => match n {
                V4Comp::T => self.matrix_ll.data[12],
                V4Comp::X => self.matrix_ll.data[13],
                V4Comp::Y => self.matrix_ll.data[14],
                V4Comp::Z => self.matrix_ll.data[15],
            },
        }
    }
}

// /// todo: As method of Christoffel instead?
// pub fn to_christoffel(&self) -> Christoffel {
//     let term1 = self(a).diff(b, c);
//     let term2 = self(b).diff(a, c);
//     let term3 = self(c).diff(a, b);
//
//     0.5 * self.matrix * (term1 + term2 - term3)
// }

/// https://github.com/wojciechczaja/GraviPy/blob/master/gravipy/tensorial.py
pub struct ReimannTensor {
    pub coords: [f64; 1],
    pub con: i8,
    pub components: [f64; 1],
    pub index_values: [i8; 1],
}

impl ReimannTensor {}

/// https://github.com/wojciechczaja/GraviPy/blob/master/gravipy/tensorial.py
pub struct RicciTensor {
    pub coords: [f64; 1],
    pub con: i8,
    pub components: [f64; 20],
    pub index_values: [i8; 1],
}

impl RicciTensor {
    // pub fn new() -> Self {}
    //
    // pub fn compute_covaraiant_component(&self) -> f64 {}
    //
    // pub fn scaler(&self) -> f64 {}
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
