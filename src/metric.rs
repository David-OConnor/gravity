//! The Metric tensor and related data.

use crate::{
    tensors::{PositIndex, PrevNext, Tensor2Config, V4Component, Vec4, Vec4Minkowski, C, COMPS},
    Arr4dMetric, Arr4dReal, C_SQ, G, H,
};

use lin_alg2::f64::{Mat4, Vec3};

/// https://github.com/einsteinpy/einsteinpy/blob/main/src/einsteinpy/metric/base_metric.py
/// A metric tensor of order (2? 3?)
/// An example: https://docs.einsteinpy.org/en/stable/_modules/einsteinpy/symbolic/metric.html#MetricTensor
///
/// See the accepted answer here for a nice explanation: https://math.stackexchange.com/questions/1410793/matrix-representations-of-tensors
#[derive(Clone)]
pub struct MetricTensor {
    /// These component identities are along diagonals; first the top left-bottom right diagonal.
    /// 0 4 7 9
    /// 4 1 5 8
    /// 7 5 2 6
    /// 9 8 6 3
    components_ll: [f64; 10],
    components_uu: [f64; 10],
    // matrix_ll: Mat4,
    // /// Ie the inverse metric; it's upper config.
    // matrix_uu: Mat4,
}

/// Eg initialization of a metric field
impl Default for MetricTensor {
    fn default() -> Self {
        let mut identity = [0.; 10];
        identity[0..4].copy_from_slice(&[1., 1., 1., 1.]);

        Self {
            components_ll: identity.clone(),
            components_uu: identity,
        }
    }
}

impl MetricTensor {
    /// Compute the metric tensor from coordinate bases. Uses the defintion `g_{m n} = e_m · e_n
    /// These bases are all lower-index.
    pub fn from_coordinate_bases(e_t: Vec4, e_x: Vec4, e_y: Vec4, e_z: Vec4) -> Self {
        // This matrix array layout is column-major, and lower-indexed.
        let mut g = Self::default();

        let l = Tensor2Config::Ll;

        // For g_{m n}, m = t
        *g.val_mut(C::T, C::T, l) = e_t.dot(e_t); // n = t
        *g.val_mut(C::T, C::X, l) = e_t.dot(e_x); // n = x
        *g.val_mut(C::T, C::Y, l) = e_t.dot(e_y); // n = y
        *g.val_mut(C::T, C::Z, l) = e_t.dot(e_z); // n = z

        // For g_{m n}, m = x
        *g.val_mut(C::X, C::X, l) = e_x.dot(e_x);
        *g.val_mut(C::X, C::Y, l) = e_x.dot(e_y);
        *g.val_mut(C::X, C::Z, l) = e_x.dot(e_z);

        // For g_{m n}, m = y
        *g.val_mut(C::Y, C::Y, l) = e_y.dot(e_y);
        *g.val_mut(C::Y, C::Z, l) = e_y.dot(e_z);

        // For g_{m n}, m = z
        *g.val_mut(C::Z, C::Z, l) = e_z.dot(e_z);

        g.generate_inverse();

        g
    }

    /// Create a new metric tensor given Scharzchild gemoetry, ie a given distance from a
    /// given black hole (or spherical object in general?) with Scharzchild radius r_s.
    pub fn new_schwarzchild(M: f64, r: f64, θ: f64) -> Self {
        // todo: phi instead of theta?
        let mut result = Self::default();

        let r_s = 2. * M * G / C_SQ;

        let a = 1. - r_s / r;
        let r_sq = r.powi(2);

        result.components_ll[0] = -a;
        result.components_ll[1] = 1. / a;
        result.components_ll[2] = r_sq;
        result.components_ll[3] = r_sq * (θ.sin()).powi(2);

        result.generate_inverse();

        result
    }

    /// Generate the inverse metric (upper-upper config).
    fn generate_inverse(&mut self) {
        let c = &self.components_ll;

        #[rustfmt::skip]
            let matrix_ll = Mat4::new([
            c[0], c[4], c[7], c[9],
            c[4], c[1], c[2], c[8],
            c[7], c[5], c[2], c[6],
            c[9], c[8], c[6], c[3]
        ]);

        let mu = matrix_ll.inverse().unwrap();
        self.components_uu = [
            mu.data[0],
            mu.data[5],
            mu.data[10],
            mu.data[15],
            mu.data[1],
            mu.data[6],
            mu.data[11],
            mu.data[2],
            mu.data[7],
            mu.data[3],
        ];
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

    /// Take the inner product of 2 vectors, using this metric; V·V = g_{μν} V^μ V^ν
    pub fn inner_product(&self, vec_a: Vec4Minkowski, vec_b: Vec4Minkowski) -> f64 {
        // let a = vec_a.as_upper();
        // let b = vec_b.as_upper();

        let mut result = 0.;

        for μ in &COMPS {
            for ν in &COMPS {
                result += self.val(*μ, *ν, Tensor2Config::Ll) * vec_a.val(*μ) * vec_b.val(*ν);
            }
        }

        result
    }

    /// Get a specific matrix component. This keeps code that uses this matrix consistent
    /// with conventions. Assumptions: Column-major array representation of the matrix, and
    /// the internal matrix is in lower-lower representation.
    pub fn val(&self, μ: V4Component, ν: V4Component, config: Tensor2Config) -> f64 {
        let g = &match config {
            Tensor2Config::Uu => &self.components_uu,
            Tensor2Config::Ul => {
                // let a = ETA_MINKOWSKI * self.matrix_ll;
                unimplemented!()
            }
            Tensor2Config::Lu => {
                unimplemented!()
            }
            Tensor2Config::Ll => &self.components_ll,
        };

        // 0 4 7 9
        // 4 1 5 8
        // 7 5 2 6
        // 9 8 6 3

        match μ {
            C::T => match ν {
                C::T => g[0],
                C::X => g[4],
                C::Y => g[7],
                C::Z => g[9],
            },

            C::X => match ν {
                C::T => g[4],
                C::X => g[1],
                C::Y => g[5],
                C::Z => g[8],
            },

            C::Y => match ν {
                C::T => g[7],
                C::X => g[5],
                C::Y => g[2],
                C::Z => g[6],
            },

            C::Z => match ν {
                C::T => g[9],
                C::X => g[8],
                C::Y => g[6],
                C::Z => g[3],
            },
        }
    }

    /// Similar to `val`, but mutable.
    /// todo: DRY with `.val`
    pub fn val_mut(&mut self, μ: V4Component, ν: V4Component, config: Tensor2Config) -> &mut f64 {
        let g = match config {
            Tensor2Config::Uu => &mut self.components_uu,
            Tensor2Config::Ul => {
                // let a = ETA_MINKOWSKI * self.matrix_ll;
                unimplemented!()
            }
            Tensor2Config::Lu => {
                unimplemented!()
            }
            Tensor2Config::Ll => &mut self.components_ll,
        };

        match μ {
            C::T => match ν {
                C::T => &mut g[0],
                C::X => &mut g[4],
                C::Y => &mut g[7],
                C::Z => &mut g[9],
            },

            C::X => match ν {
                C::T => &mut g[4],
                C::X => &mut g[1],
                C::Y => &mut g[5],
                C::Z => &mut g[8],
            },

            C::Y => match ν {
                C::T => &mut g[7],
                C::X => &mut g[5],
                C::Y => &mut g[2],
                C::Z => &mut g[6],
            },

            C::Z => match ν {
                C::T => &mut g[9],
                C::X => &mut g[8],
                C::Y => &mut g[6],
                C::Z => &mut g[3],
            },
        }
    }
}

/// Group that includes the metric at a point, and at points surrounding it, an infinetesimal difference
/// in both directions along each spacial axis.
/// This is useful for when we have a Metric tensor we can find analytically; allows for accurate
/// derivatives.
#[derive(Clone)]
pub struct MetricGridWDiffs {
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

impl MetricGridWDiffs {
    pub fn new(grid_n: usize) -> Self {
        let data = crate::new_data_metric(grid_n);

        Self {
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

#[derive(Clone, Default)]
pub struct MetricWDiffs {
    pub on_pt: MetricTensor,
    pub t_prev: MetricTensor,
    pub t_next: MetricTensor,
    pub x_prev: MetricTensor,
    pub x_next: MetricTensor,
    pub y_prev: MetricTensor,
    pub y_next: MetricTensor,
    pub z_prev: MetricTensor,
    pub z_next: MetricTensor,
}

impl MetricWDiffs {
    pub fn from_grid(metrics: &Arr4dMetric, p_i: &PositIndex, grid_n: usize) -> Self {
        // todo: Make sure is aren't at the edges. If so, return etc.
        if p_i.t == 0
            || p_i.t == grid_n - 1
            || p_i.x == 0
            || p_i.x == grid_n - 1
            || p_i.y == 0
            || p_i.y == grid_n - 1
            || p_i.z == 0
            || p_i.z == grid_n - 1
        {
            println!("On edge; return an error or some other workaround.")
        }

        // todo: Refs instead of cloning, to make this more performant

        Self {
            on_pt: metrics[p_i.t][p_i.x][p_i.y][p_i.z].clone(),
            t_prev: metrics[p_i.t - 1][p_i.x][p_i.y][p_i.z].clone(),
            t_next: metrics[p_i.t - 1][p_i.x][p_i.y][p_i.z].clone(),
            x_prev: metrics[p_i.t][p_i.x - 1][p_i.y][p_i.z].clone(),
            x_next: metrics[p_i.t][p_i.x + 1][p_i.y][p_i.z].clone(),
            y_prev: metrics[p_i.t][p_i.x][p_i.y - 1][p_i.z].clone(),
            y_next: metrics[p_i.t][p_i.x][p_i.y + 1][p_i.z].clone(),
            z_prev: metrics[p_i.t][p_i.x][p_i.y][p_i.z - 1].clone(),
            z_next: metrics[p_i.t][p_i.x][p_i.y][p_i.z + 1].clone(),
        }
    }

    /// todo: θ or phi?
    pub fn new_schwarz(M: f64, posit_sample: Vec4Minkowski, posit_mass: Vec3) -> Self {
        let posit_t_prev = Vec4Minkowski::new(
            posit_sample.t() - H,
            posit_sample.x(),
            posit_sample.y(),
            posit_sample.z(),
        );
        let posit_t_next = Vec4Minkowski::new(
            posit_sample.t() + H,
            posit_sample.x(),
            posit_sample.y(),
            posit_sample.z(),
        );

        let posit_x_prev = Vec4Minkowski::new(
            posit_sample.t(),
            posit_sample.x() - H,
            posit_sample.y(),
            posit_sample.z(),
        );
        let posit_x_next = Vec4Minkowski::new(
            posit_sample.t(),
            posit_sample.x() + H,
            posit_sample.y(),
            posit_sample.z(),
        );

        let posit_y_prev = Vec4Minkowski::new(
            posit_sample.t(),
            posit_sample.x(),
            posit_sample.y() - H,
            posit_sample.z(),
        );
        let posit_y_next = Vec4Minkowski::new(
            posit_sample.t(),
            posit_sample.x(),
            posit_sample.y() + H,
            posit_sample.z(),
        );

        let posit_z_prev = Vec4Minkowski::new(
            posit_sample.t(),
            posit_sample.x(),
            posit_sample.y(),
            posit_sample.z() - H,
        );
        let posit_z_next = Vec4Minkowski::new(
            posit_sample.t(),
            posit_sample.x(),
            posit_sample.y(),
            posit_sample.z() + H,
        );

        //todo: theta or phi?
        let (r_on_pt, θ_on_pt) = crate::find_schwarz_params(posit_sample, posit_mass);

        let (r_t_prev, θ_t_prev) = crate::find_schwarz_params(posit_t_prev, posit_mass);
        let (r_t_next, θ_t_next) = crate::find_schwarz_params(posit_t_next, posit_mass);
        let (r_x_prev, θ_x_prev) = crate::find_schwarz_params(posit_x_prev, posit_mass);
        let (r_x_next, θ_x_next) = crate::find_schwarz_params(posit_x_next, posit_mass);
        let (r_y_prev, θ_y_prev) = crate::find_schwarz_params(posit_y_prev, posit_mass);
        let (r_y_next, θ_y_next) = crate::find_schwarz_params(posit_y_next, posit_mass);
        let (r_z_prev, θ_z_prev) = crate::find_schwarz_params(posit_z_prev, posit_mass);
        let (r_z_next, θ_z_next) = crate::find_schwarz_params(posit_z_next, posit_mass);

        Self {
            on_pt: MetricTensor::new_schwarzchild(M, r_on_pt, θ_on_pt),
            t_prev: MetricTensor::new_schwarzchild(M, r_t_prev, θ_t_prev),
            t_next: MetricTensor::new_schwarzchild(M, r_t_next, θ_t_next),
            x_prev: MetricTensor::new_schwarzchild(M, r_x_prev, θ_x_prev),
            x_next: MetricTensor::new_schwarzchild(M, r_x_next, θ_x_next),
            y_prev: MetricTensor::new_schwarzchild(M, r_y_prev, θ_y_prev),
            y_next: MetricTensor::new_schwarzchild(M, r_y_next, θ_y_next),
            z_prev: MetricTensor::new_schwarzchild(M, r_z_prev, θ_z_prev),
            z_next: MetricTensor::new_schwarzchild(M, r_z_next, θ_z_next),
        }
    }

    // pub fn val(comp: V4Component, prev_next: PrevNext) -> Arr4dMetric {
    //     match comp {
    //         C::T => {
    //             match prev_next {
    //                 PrevNext::Prev => self,
    //             }
    //         }
    //         C::X => {
    //
    //         }
    //         C::Y => {
    //
    //         }
    //         C::Z => {
    //
    //         }
    //     }
    // }
}
