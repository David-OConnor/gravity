//! Modeling a galaxy, with the hope of learning something about Dark matter, or perhaps
//! physics/chem modelling and rendering skills in general.
//!
//! https://docs.einsteinpy.org/en/stable/getting_started.html
//!
//!
//!
//!
//! Thought: Can you model 3D space and/or 4D spacetime as a being a manifold in a higher-dimension
//! space? perhaps, just one higher dimension? Perhaps each point in space has a mass associated
//! with it. Figure out if you can then, using the curvature (eg second derivatives) of the space
//! (in 4d and/or 5d?) model gravity.

use lin_alg2::f64::{Mat4, Vec3, Vec4};

mod render;

pub type Arr3dReal = Vec<Vec<Vec<f64>>>;

pub type Arr3d = Vec<Vec<Vec<Cplx>>>;
pub type Arr3dMetric = Vec<Vec<Vec<MetricTensor>>>;

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

/// Christoffel symbol.
pub struct Christoffel {

}

/// 4-vector in Minkowski coordinates. Capable of representations in covariant and
/// contravariant forms.
struct Vec4Minkowski {
    /// We store the numerical representation internally in its contravariant form (Upper index).
    /// We use the `as_lower` method to get its value in covariant form.
    xyz_u: Vec3,
    t_u: f64,
}

impl Vec4Minkowski {
    /// Get the vector's contravariant (upper-index) numerical form.
    pub fn as_upper(&self) -> Vec4 {
        // todo: Should t be the 0th index?
        Vec4::new(self.xyz_u.x, self.xyz_u.y, self.xyz_u.z, self.t_u)
    }

    /// Get the vector's covariant (lower-index) numerical form.
    pub fn as_lower(&self, metric: MetricTensor) -> Vec4 {
        metric.matrix_uu.as_config(Tensor2Config::Ll) * self.as_upper()
    }
}

/// https://github.com/einsteinpy/einsteinpy/blob/main/src/einsteinpy/metric/base_metric.py
/// A metric tensor of order (2? 3?)
/// An example: https://docs.einsteinpy.org/en/stable/_modules/einsteinpy/symbolic/metric.html#MetricTensor
///
/// See the accepted answer here for a nice explanation: https://math.stackexchange.com/questions/1410793/matrix-representations-of-tensors
pub struct MetricTensor {
    // /// Configuration of indices
    // config: Tensor2Config,
    /// Matrix representation: A 4x4 matrix of minkowski space. This is in UU config. ie g^{m n}.
    /// We get other configs using our `as_config` method.
    matrix_uu: Mat4,
}

impl MetricTensor {
    pub fn as_config(&self, config: Tensor2Config) -> Mat4 {
        // todo
        match config {
            Tensor2Config::Uu => {
                self.matrix_uu
            }
            Tensor2Config::Ul => {
                let a = ETA_MINKOWSKI * self.matrix;

            }
            Tensor2Config::Lu => {

            }
            Tensor2Config::Ll => {
                self.matrix_uu.inverse().unwrap()
            }
        }
    }

    /// Take the inner product of 2 vectors, using this metric; V·V = g_{mn} V^m V^n
    pub fn inner_product(&self, vec_a: Vec4Minkowski, vec_b: Vec4Minkowski) -> Vec4Minkowski {
        self.matrix_uu *
            self.change_vec_config(vec_a, Tensor1Config::L).to_vec4() *
            self.change_vec_config(vec_b, Tensor1Config::L).to_vec4()
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

// /// Returns the inverse of the Metric.
// /// Returns contravariant Metric if it is originally covariant or vice-versa.
// pub fn inverse(&self) -> Option<Mat4> {
//     self.matrix.inverse()
// }
//
// /// Changes the index configuration(contravariant/covariant)
// /// todo: N/A for metric tensor.
// pub fn _change_config(&mut self, config: TensorConfig) {
//
// }
//
// /// Returns a covariant instance of the given metric tensor.
// ///same instance if the configuration is already lower or
// /// inverse of given metric if configuration is upper
// pub fn lower_config(&self) -> Self {
//
// }
//
// /// Performs a Lorentz transform on the tensor.
// pub fn lorentz_transform(&self, transformation_matrix: i8) -> Self {
//
// }

/// https://github.com/wojciechczaja/GraviPy/blob/master/gravipy/tensorial.py
pub struct ReimannTensor {
    pub coords: [f64; 1],
    pub con: i8,
    pub components: [f64; 1],
    pub index_values: [i8; 1],
}

impl ReimannTensor {

}

/// https://github.com/wojciechczaja/GraviPy/blob/master/gravipy/tensorial.py
pub struct RicciTensor {
    pub coords: [f64; 1],
    pub con: i8,
    pub components: [f64; 20],
    pub index_values: [i8; 1],
}

impl RicciTensor {
    pub fn new() -> Self {

    }

    pub fn compute_covaraiant_component(&self) -> f64 {

    }

    pub fn scaler(&self) -> f64 {
    }
}

pub struct EinsteinTensor {
    pub coords: [f64; 1],
    pub con: i8,
    pub components: [f64; 1],
    pub index_values: [i8; 1],
}

impl EinsteinTensor {
    pub fn new(ricci: &RicciTensor) {

    }
}

pub struct GeodesicTensor {

}

struct Particle {
    pub posit: Vec3,
    pub mass: f64,
}

pub struct State {
    particles: Vec<Particle>,
    field_grav: Arr3dReal,
    field_metric: Arr3dMetric,
}


/// Make a new 3D grid of position, time 4-vecs in Minknowski space
pub fn new_data_vec(n: usize) -> Arr3d4Vec {
    let mut z = Vec::new();
    z.resize(n, Vec4Minkowski { t: 0., x_y_z: Vec3::new_zero()});

    let mut y = Vec::new();
    y.resize(n, z);

    let mut x = Vec::new();
    x.resize(n, y);

    x
}

fn main() {
    let state = State {
        particles: Vec::new(),
    };

    render::render(state);
}
