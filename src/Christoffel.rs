//! Separated from tensors since it's so verbose. (And not a tensor)

use crate::{
    tensors::{Comp, Vec4Minkowski},
    Arr3dMetric,
};

/// Christoffel symbol. (Not a tensor)
/// D = 4 and n = 3, so this has 4^3 = 64 components.
pub struct Christoffel {
    /// This can be thought of a a col-major 4x4x4 matrix.
    components: [f64; 64],
}

impl Christoffel {
    /// Compute from the metric tensor and its numerical derivatives.
    /// Γ^a _{b c} = g^{a d} ( d_b g_{c d} + d_c g_{b d} - d_d g_{b c})
    pub fn from_metric(metrics: &crate::Arr3dMetric, posit: &Vec4Minkowski) {}

    /// Same idea as for the metric tensor.
    /// todo boy: This is a bear.
    /// λ is the upper index.
    pub fn val(&self, λ: Comp, μ: Comp, ν: Comp) -> f64 {
        let d = &self.components.data;

        // todo: Order? Is this reversed?
        match λ {
            Comp::T => match μ {
                Comp::T => match ν {
                    Comp::T => d[0],
                    Comp::X => d[1],
                    Comp::Y => d[2],
                    Comp::Z => d[3],
                },

                Comp::X => match ν {
                    Comp::T => d[4],
                    Comp::X => d[5],
                    Comp::Y => d[6],
                    Comp::Z => d[7],
                },

                Comp::Y => match ν {
                    Comp::T => d[8],
                    Comp::X => d[9],
                    Comp::Y => d[10],
                    Comp::Z => d[11],
                },

                Comp::Z => match ν {
                    Comp::T => d[12],
                    Comp::X => d[13],
                    Comp::Y => d[14],
                    Comp::Z => d[15],
                },
            },
            Comp::X => match μ {
                Comp::T => match ν {
                    Comp::T => d[16],
                    Comp::X => d[17],
                    Comp::Y => d[18],
                    Comp::Z => d[19],
                },

                Comp::X => match ν {
                    Comp::T => d[20],
                    Comp::X => d[21],
                    Comp::Y => d[22],
                    Comp::Z => d[23],
                },

                Comp::Y => match ν {
                    Comp::T => d[24],
                    Comp::X => d[25],
                    Comp::Y => d[26],
                    Comp::Z => d[27],
                },

                Comp::Z => match ν {
                    Comp::T => d[28],
                    Comp::X => d[29],
                    Comp::Y => d[30],
                    Comp::Z => d[31],
                },
            },
            Comp::Y => match μ {
                Comp::T => match ν {
                    Comp::T => d[32],
                    Comp::X => d[33],
                    Comp::Y => d[34],
                    Comp::Z => d[35],
                },

                Comp::X => match ν {
                    Comp::T => d[36],
                    Comp::X => d[37],
                    Comp::Y => d[38],
                    Comp::Z => d[39],
                },

                Comp::Y => match ν {
                    Comp::T => d[40],
                    Comp::X => d[41],
                    Comp::Y => d[42],
                    Comp::Z => d[43],
                },

                Comp::Z => match ν {
                    Comp::T => d[44],
                    Comp::X => d[45],
                    Comp::Y => d[46],
                    Comp::Z => d[47],
                },
            },
            Comp::Z => match μ {
                Comp::T => match ν {
                    Comp::T => d[48],
                    Comp::X => d[49],
                    Comp::Y => d[50],
                    Comp::Z => d[51],
                },

                Comp::X => match ν {
                    Comp::T => d[52],
                    Comp::X => d[53],
                    Comp::Y => d[54],
                    Comp::Z => d[55],
                },

                Comp::Y => match ν {
                    Comp::T => d[56],
                    Comp::X => d[57],
                    Comp::Y => d[58],
                    Comp::Z => d[59],
                },

                Comp::Z => match ν {
                    Comp::T => d[60],
                    Comp::X => d[61],
                    Comp::Y => d[62],
                    Comp::Z => d[63],
                },
            },
        }
    }
    pub fn val_mut(&mut self, m: Comp, n: Comp) -> &mut f64 {}
}
