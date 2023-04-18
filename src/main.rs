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

use lin_alg2::f64::Vec3;

mod render;

pub type Arr3dReal = Vec<Vec<Vec<f64>>>;

pub type Arr3d = Vec<Vec<Vec<Cplx>>>;
pub type Arr3dMetric = Vec<Vec<Vec<MetricTensor>>>;

#[derive(Clone, Copy)]
enum TensorConfig {
    Uu,
    Ul,
    Lu,
    Ll,
}

/// https://github.com/einsteinpy/einsteinpy/blob/main/src/einsteinpy/metric/base_metric.py
/// A metric tensor of order (2? 3?)
/// An example: https://docs.einsteinpy.org/en/stable/_modules/einsteinpy/symbolic/metric.html#MetricTensor
///
/// See the accepted answer here for a nice explanation: https://math.stackexchange.com/questions/1410793/matrix-representations-of-tensors
pub struct MetricTensor {
    // R1: L to R. R2: L to R. (etc)
    // Noting that this is a symmetric matrix.
    matrix_vals: [f64; 10],
    pub coords: [f64; 1],
    pub con: i8,
    pub components: [f64; 1],
    pub index_values: [i8; 1],
}

impl MetricTensor {
    /// Returns the inverse of the Metric.
    /// Returns contravariant Metric if it is originally covariant or vice-versa.
    pub fn inverse(&self) -> Self {

    }

    /// Changes the index configuration(contravariant/covariant)
    /// todo: N/A for metric tensor.
    pub fn _change_config(&mut self, config: TensorConfig) {

    }

    /// Returns a covariant instance of the given metric tensor.
    ///same instance if the configuration is already lower or
    /// inverse of given metric if configuration is upper
    pub fn lower_config(&self) -> Self {

    }

    /// Performs a Lorentz transform on the tensor.
    pub fn lorentz_transform(&self, transformation_matrix: i8) -> Self {

    }
}

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

struct Vec4 {
    pub x_y_z: Vec3,
    pub t: f64,
}



/// Make a new 3D grid of position, time 4-vecs in Minknowski space
pub fn new_data_vec(n: usize) -> Arr3d4Vec {
    let mut z = Vec::new();
    z.resize(n, Vec4 { t: 0., x_y_z: Vec3::new_zero()});

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
