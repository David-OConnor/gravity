#![allow(non_ascii_idents)]
#![allow(non_snake_case)]
#![allow(confusable_idents)] // nu and v

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

mod christoffel;
mod render;
mod tensors;
mod ui;
mod reimann;
mod util;

use crate::{
    tensors::{MetricTensor, Vec4, Vec4Minkowski},
    christoffel::Christoffel,
};

// Gravitational constant.
const C: f64 = 1.;
const C_SQ: f64 = 1.;
const G: f64 = 1.;

// Used for numerical derivatives
pub const H: f64 = 0.0001;

pub type Arr4dReal = Vec<Vec<Vec<Vec<f64>>>>;
// pub type Arr3d = Vec<Vec<Vec<Cplx>>>;
// pub type Arr3dVec = Vec<Vec<Vec<Vec3>>>;
// pub type Arr3dVec4 = Vec<Vec<Vec<Vec4Minkowski>>>;
// pub type Arr3dMetric = Vec<Vec<Vec<MetricTensor>>>;
pub type Arr4dVec4 = Vec<Vec<Vec<Vec<Vec4Minkowski>>>>;
pub type Arr4dMetric = Vec<Vec<Vec<Vec<MetricTensor>>>>;
pub type Arr4dChristoffel = Vec<Vec<Vec<Vec<Christoffel>>>>;

// todo: Consider the form of this.
fn gamma(v: f64) -> f64 {
    1. / (1. - v.powi(2) / C_SQ).sqrt()
}

// Our equations make heavy use of this; keep things terse.
pub type C = tensors::V4Component;

/// Used for indexing into spacetime grids, eg for metric tensors.
#[derive(Clone, Copy)]
pub struct PositIndex {
    pub t: usize,
    pub x: usize,
    pub y: usize,
    pub z: usize,
}

struct Particle {
    pub posit: Vec3,
    pub mass: f64,
}

pub struct State {
    // A single swarzchild point
    pub schwarzchild_coords: Vec3, // interesting how this is a 3-vec
    pub schwarzchild_rad: f64,
    pub particles: Vec<Particle>,
    // pub field_metric: Arr4dMetric,
    pub field_metric: tensors::MetricWDiffs,
    pub field_christoffel: Arr4dChristoffel,
    pub grid_posits: Arr4dVec,
    pub grid_n: usize,

}

/// The path a particle takes through spacetime. Can be a geodesic.
pub struct Worldline {
    /// A list of events, with constant proper-time spacing.
    pub events: Vec<Vec4Minkowski>,
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
    /// τ_AB = \int (0, 1) sqrt(-g_μv(x^alpha(\sigma)) dx^u/dsigma dx^v / dsigma) dsigma
    pub fn calc_proper_time(&self) -> f64 {
        0.
    }

    // /// Find a geodesic:
    // /// d^2 x^μ / dλ^2 + Γ^μ_ρσ dx^ρ /dλ dx^σ /dλ = 0
    // pub fn find_geodesic(
    //     grid_posits: &Arr3dVec4Minkowski,
    //     metrics: &Arr3dMetric,
    //     posit: &Vec4Minkowski,
    // ) -> Self {
    //     // todo:
    //     let Γ = Christoffel::from_metric(metrics, posit);
    //     // This is a differential equation; solving it may not be a simple fn...
    //
    // }
}

/// Make a new 3D grid of position, time 4-vecs in Minknowski space
pub fn new_data_vec(n: usize) -> Arr4dVec4 {
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

    let mut t = Vec::new();
    t.resize(n, x);

    t
}

/// Make a new 3D grid of Metric tensors. Used to model the curvature of space.
/// Defaults to the identity transform, ie flat space.
pub fn new_data_metric(n: usize) -> Arr4dMetric {
    let mut z = Vec::new();

    z.resize(n, MetricTensor::default());

    let mut y = Vec::new();
    y.resize(n, z);

    let mut x = Vec::new();
    x.resize(n, y);

    let mut t = Vec::new();
    t.resize(n, x);

    t
}

pub fn new_data_christoffel(n: usize) -> Arr4dChristoffel {
    let mut z = Vec::new();

    z.resize(n, christoffel::default());

    let mut y = Vec::new();
    y.resize(n, z);

    let mut x = Vec::new();
    x.resize(n, y);

    let mut t = Vec::new();
    t.resize(n, x);

    t
}

/// Calcaulte teh scharzchild spacetime interal.
pub fn swarzchild_interval(r_s: f64, r: f64, dr: f64, dθ: f64, dφ: f64) -> f64 {
    let r_const = 1. - r_s / r;
    let g_omega = dθ.powi(2) + ((dφ.powi(2)).sin()).powi(2);

    -r_const * C_SQ * dt.powi(2) + 1. / r_const * dr.powi(2) + r.powi(2) * g_omega
}

/// Update our grid positions. Run this when we change grid bounds or spacing.
pub fn update_grid_posits(
    grid_posits: &mut Arr3dVec,
    grid_min: f64,
    grid_max: f64,
    n: usize,
) {
    let grid_lin = util::linspace((grid_min, grid_max), n);

    // Set up a grid with values that increase in distance the farther we are from the center.
    let mut grid_1d = vec![0.; n];


    for (i, x) in grid_lin.iter().enumerate() {
        for (j, y) in grid_lin.iter().enumerate() {
            for (k, z) in grid_lin.iter().enumerate() {
                grid_posits[i][j][k] = Vec3::new(*x, *y, *z);
            }
        }
    }
}

/// Helper fn for generating neighboring points for the Schwarzchild metric
fn find_r_θ(posit_sample: Vec4Minkowski, posit_mass: Vec3) -> (f64, f64) {
    let diff = posit_sample - posit_mass; // todo: Invalid op...
    // todo: Mass is at a Vec3, sample is at a vec4...
    let r = posit.mag_sq().sqrt();
// todo: QC this
    let θ = (-posit.t.powi(2) + posit.x.powi(2) + posit.y.powi(2)).sqrt().atan2(posit.z);

    (r, θ)
}

fn main() {
    let schwarzchild_coords = Vec3::new_zero();
    let schwarzchild_rad = 2.;

    let grid_n = 30;

    let mut field_metric = tensors::MetricWDiffs::new(grid_n);
    let mut field_christoffel = new_data_christoffel(grid_n);

    let mut grid_posits = new_data_vec(grid_n);
    let grid_max = 10.;
    let grid_min = -grid_max;
    update_grid_posits(&mut grid_posits, grid_min, grid_max, grid_n);

    // For now, we take advantage of the Analytic Swarzchild metric to take
    // accurate derivatives a small distance away.

    // todo: Try to include a grid-based derivative later for when moving away from Swarchild.

    for i in 0..grid_n {
        for j in 0..grid_n {
            for k in 0..grid_n {
                for l in 0..grid_n {
                    let posit = grid_posits[i][j][k][l];

                    let posit_t_prev = Vec4Minkowski::new(posit.t - H, posit_sample.x, posit_sample.y, posit_sample.z);
                    let posit_t_next = Vec4Minkowski::new(posit.t + H, posit_sample.x, posit_sample.y, posit_sample.z);

                    let posit_x_prev = Vec4Minkowski::new(posit.t, posit_sample.x - H, posit_sample.y, posit_sample.z);
                    let posit_x_next = Vec4Minkowski::new(posit.t, posit_sample.x + H, posit_sample.y, posit_sample.z);

                    let posit_y_prev = Vec4Minkowski::new(posit.t, posit_sample.x, posit_sample.y - H, posit_sample.z);
                    let posit_y_next = Vec4Minkowski::new(posit.t, posit_sample.x, posit_sample.y + H, posit_sample.z);

                    let posit_z_prev = Vec4Minkowski::new(posit.t, posit_sample.x, posit_sample.y, posit_sample.z - H);
                    let posit_z_next = Vec4Minkowski::new(posit.t, posit_sample.x, posit_sample.y, posit_sample.z + H);


                    let (r_on_pt, θ_on_pt ) = find_r_θ(posit, schwarzchild_coords);

                    let (r_t_prev, θ_t_prev ) = find_r_θ(posit_t_prev, schwarzchild_coords);
                    let (r_t_next, θ_t_next ) = find_r_θ(posit_t_next, schwarzchild_coords);

                    let (r_x_prev, θ_x_prev ) = find_r_θ(posit_x_prev, schwarzchild_coords);
                    let (r_x_next, θ_x_next ) = find_r_θ(posit_x_next, schwarzchild_coords);

                    let (r_y_prev, θ_y_prev ) = find_r_θ(posit_y_prev, schwarzchild_coords);
                    let (r_y_next, θ_y_next ) = find_r_θ(posit_y_next, schwarzchild_coords);

                    let (r_z_prev, θ_z_prev ) = find_r_θ(posit_z_prev, schwarzchild_coords);
                    let (r_z_next, θ_z_next ) = find_r_θ(posit_z_next, schwarzchild_coords);


                    field_metric.on_pt[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_on_pt, θ_on_pt
                    );

                    field_metric.t_prev[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_t_prev, θ_t_prev
                    );
                    field_metric.t_next[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_t_next, θ_t_next
                    );
                    field_metric.x_prev[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_x_prev, θ_x_prev
                    );
                    field_metric.x_next[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_x_next, θ_x_next
                    );
                    field_metric.y_prev[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_y_prev, θ_y_prev
                    );
                    field_metric.y_next[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_y_next, θ_y_next
                    );
                    field_metric.z_prev[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_z_prev, θ_z_prev
                    );
                    field_metric.z_next[i][j][k][l] = MetricTensor::new_schwarzchild(
                        schwarzchild_rad, r_z_next, θ_z_next
                    );
                }
            }
        }
    }


    // // Second iteration, since we need the metric's neighbors.
    // for i in 0..grid_n {
    //     for j in 0..grid_n {
    //         for k in 0..grid_n {
    //             for l in 0..grid_n {
    //                 let posit_sample = grid_posits[i][j][k][l];
    //
    //                 // field_christoffel[i][j][k][l] = Christoffel::from_metrics()
    //             }
    //         }
    //     }
    // }


    let state = State {
        schwarzchild_coords,
        schwarzchild_rad,
        particles: Vec::new(),
        field_metric,
        field_christoffel,
        grid_posits,
        grid_n,
    };

    render::render(state);
}
