#![allow(non_ascii_idents)]
#![allow(non_snake_case)]

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

use lin_alg2::f64::{Mat4, Vec3, Vec4 as Vec4Linalg};
use std::f64::consts::TAU;

mod Christoffel;
mod render;
mod tensors;
mod ui;

use crate::tensors::{Christoffel, MetricTensor, Vec4, Vec4Minkowski};

// Gravitational constant.
const C: f64 = 1.;
const C_SQ: f64 = 1.;
const G: f64 = 1.; // todo: Change A/R.

// pub type Arr3dReal = Vec<Vec<Vec<f64>>>;
// pub type Arr3d = Vec<Vec<Vec<Cplx>>>;
// pub type Arr3dVec = Vec<Vec<Vec<Vec3>>>;
pub type Arr3dVec4 = Vec<Vec<Vec<Vec4Minkowski>>>;
pub type Arr3dMetric = Vec<Vec<Vec<MetricTensor>>>;

// todo: Consider the form of this.
fn gamma(v: f64) -> f64 {
    1. / (1. - v.powi(2) / C_SQ).sqrt()
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

// todo: Geodesic struct?

/// Find a geodesic
pub fn compute_geodesic(
    grid_posits: &Arr3dVec4Minkowski,
    metrics: &Arr3dMetric,
    posit: &Vec4Minkowski,
) -> Vec<Vec4Minkowski> {
    // todo:
    let Γ = Christoffel::from_metric(metrics, posit);
}

/// Make a new 3D grid of position, time 4-vecs in Minknowski space
pub fn new_data_vec(n: usize) -> Arr3dVec4 {
    let mut z = Vec::new();
    z.resize(
        n,
        Vec4Minkowski {
            value_upper: Vec4 {
                t: 1.,
                x: 0.,
                y: 0.,
                z: 0.,
            },
        },
    );

    let mut y = Vec::new();
    y.resize(n, z);

    let mut x = Vec::new();
    x.resize(n, y);

    x
}

/// Make a new 3D grid of Metric tensors. Used to model the curvature of space.
/// Defaults to the identity transform, ie flat space.
pub fn new_data_metric(n: usize) -> Arr3dMetric {
    let mut z = Vec::new();

    z.resize(
        n,
        MetricTensor {
            matrix_ll: Mat4::new_identity(),
        },
    );

    let mut y = Vec::new();
    y.resize(n, z);

    let mut x = Vec::new();
    x.resize(n, y);

    x
}

fn main() {
    // let state = State {
    //     particles: Vec::new(),
    // };

    // render::render(state);
}
