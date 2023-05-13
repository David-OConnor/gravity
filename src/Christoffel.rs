//! Separated from tensors since it's so verbose. (And not a tensor)

use std::collections::HashMap;

use crate::{
    tensors::{V4Component, Vec4Minkowski},
    Arr3dMetric,
    C,
};
use crate::tensors::MetricTensor;

/// We use this for indexing into metric tensor collections used in derivatives.
#[derive(Clone, Copy)]
enum PrevNext {
    P,
    N,
}

/// Christoffel symbol. (Not a tensor)
/// D = 4 and n = 3, so this has 4^3 = 64 components.
pub struct Christoffel {
    /// This can be thought of a a col-major 4x4x4 matrix.
    components: [f64; 64],
}

impl Christoffel {
    /// Helper function to calculate a single component. Uses the index convention shown in `from_metric`'s
    /// function description.
    fn calc_component(a: C, b: C, c: C, g_this: &MetricTensor, g_diffs: &HashMap<(C, PrevNext), &MetricTensor>) -> f64 {
        let mut result = 0.;

        for d in &[C::T, C::X, C::Y, C::Z] {
            result += g_this.get(a, d) * (
                (g_diffs.get(&(b, PrevNext::N)).get(c, d) - g_diffs.get(&(b, PrevNext::P)).get(c, d)) / 2. +
                    (g_diffs.get(&(c, PrevNext::N)).get(b, d) - g_diffs.get(&(c, PrevNext::P)).get(b, d)) / 2. +
                    (g_diffs.get(&(*d, PrevNext::N)).get(b, c) - g_diffs.get(&(*d, PrevNext::P)).get(b, c)) / 2.
            )
        }

        result
    }

    /// Compute from the metric tensor and its numerical derivatives.
    /// Γ^a _{b c} = g^{a d} ( d_b g_{c d} + d_c g_{b d} - d_d g_{b c})
    pub fn from_metric(metrics: &crate::Arr4dMetric, posit: &Vec4Minkowski, p_i: crate::PositIndex) -> Self {
        let mut result = Self { components: [0.; 64] };
        
        // todo: Make sure is aren't at the edges. If so, return etc.

        // todo: DO you want to do a midpoint, or just one-side? It's a first-deriv, so you
        // todo could pick either.
        let g_this = &metrics[p_i.t][p_i.x][p_i.y][p_i.z];
        let g_t_prev = &metrics[p_i.t - 1][p_i.x][p_i.y][p_i.z];
        let g_t_next = &metrics[p_i.t - 1][p_i.x][p_i.y][p_i.z];
        let g_x_prev = &metrics[p_i.t][p_i.x - 1][p_i.y][p_i.z];
        let g_x_next = &metrics[p_i.t][p_i.x + 1][p_i.y][p_i.z];
        let g_y_prev = &metrics[p_i.t][p_i.x][p_i.y - 1][p_i.z];
        let g_y_next = &metrics[p_i.t][p_i.x][p_i.y + 1][p_i.z];
        let g_z_prev = &metrics[p_i.t][p_i.x][p_i.y][p_i.z - 1];
        let g_z_next = &metrics[p_i.t][p_i.x][p_i.y][p_i.z + 1];

        let mut metrics = HashMap::new();
        metrics.insert((C::T, PrevNext::P), g_t_prev);
        metrics.insert((C::T, PrevNext::N), g_t_next);
        metrics.insert((C::X, PrevNext::P), g_x_prev);
        metrics.insert((C::X, PrevNext::N), g_x_next);
        metrics.insert((C::Y, PrevNext::P), g_y_prev);
        metrics.insert((C::Y, PrevNext::N), g_y_next);
        metrics.insert((C::Z, PrevNext::P), g_z_prev);
        metrics.insert((C::Z, PrevNext::N), g_z_next);

        // todo: You need a way to get metric's upper config by Comp label!! You're currently using
        // todo the lower-rep version!

        let comps = [C::T, C::X, C::Y, C::Z];

        for a in &comps {
            for b in &comps {
                for c in &comps {
                    result.get_mut(a, b, c) += Self::calc_component(*a, *b, *c, &g_this, &metrics);
                }
            }
        }

        result
    }

    /// Same idea as for the metric tensor.
    /// todo boy: This is a bear.
    /// λ is the upper index.
    pub fn val(&self, λ: V4Component, μ: V4Component, ν: V4Component) -> f64 {
        let d = &self.components.data;

        // todo: Order? Is this reversed?
        match λ {
            C::T => match μ {
                C::T => match ν {
                    C::T => d[0],
                    C::X => d[1],
                    C::Y => d[2],
                    C::Z => d[3],
                },

                C::X => match ν {
                    C::T => d[4],
                    C::X => d[5],
                    C::Y => d[6],
                    C::Z => d[7],
                },

                C::Y => match ν {
                    C::T => d[8],
                    C::X => d[9],
                    C::Y => d[10],
                    C::Z => d[11],
                },

                C::Z => match ν {
                    C::T => d[12],
                    C::X => d[13],
                    C::Y => d[14],
                    C::Z => d[15],
                },
            },
            C::X => match μ {
                C::T => match ν {
                    C::T => d[16],
                    C::X => d[17],
                    C::Y => d[18],
                    C::Z => d[19],
                },

                C::X => match ν {
                    C::T => d[20],
                    C::X => d[21],
                    C::Y => d[22],
                    C::Z => d[23],
                },

                C::Y => match ν {
                    C::T => d[24],
                    C::X => d[25],
                    C::Y => d[26],
                    C::Z => d[27],
                },

                C::Z => match ν {
                    C::T => d[28],
                    C::X => d[29],
                    C::Y => d[30],
                    C::Z => d[31],
                },
            },
            C::Y => match μ {
                C::T => match ν {
                    C::T => d[32],
                    C::X => d[33],
                    C::Y => d[34],
                    C::Z => d[35],
                },

                C::X => match ν {
                    C::T => d[36],
                    C::X => d[37],
                    C::Y => d[38],
                    C::Z => d[39],
                },

                C::Y => match ν {
                    C::T => d[40],
                    C::X => d[41],
                    C::Y => d[42],
                    C::Z => d[43],
                },

                C::Z => match ν {
                    C::T => d[44],
                    C::X => d[45],
                    C::Y => d[46],
                    C::Z => d[47],
                },
            },
            C::Z => match μ {
                C::T => match ν {
                    C::T => d[48],
                    C::X => d[49],
                    C::Y => d[50],
                    C::Z => d[51],
                },

                C::X => match ν {
                    C::T => d[52],
                    C::X => d[53],
                    C::Y => d[54],
                    C::Z => d[55],
                },

                C::Y => match ν {
                    C::T => d[56],
                    C::X => d[57],
                    C::Y => d[58],
                    C::Z => d[59],
                },

                C::Z => match ν {
                    C::T => d[60],
                    C::X => d[61],
                    C::Y => d[62],
                    C::Z => d[63],
                },
            },
        }
    }
    // pub fn val_mut(&mut self, m: Comp, n: Comp) -> &mut f64 {} // todo
}
