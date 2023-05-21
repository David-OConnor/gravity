//! Separated from tensors since it's so verbose. (And not a tensor)

use std::collections::HashMap;

use crate::{
    metric::{MetricTensor, MetricWDiffs},
    tensors::{PositIndex, PrevNext, Tensor2Config, V4Component, Vec4, Vec4Minkowski, C, COMPS},
    Arr4dMetric,
};

/// Christoffel symbol. (Not a tensor)
/// D = 4 and n = 3, so this has 4^3 = 64 components.
#[derive(Clone)]
pub struct Christoffel {
    /// This can be thought of a a col-major 4x4x4 matrix.
    /// There are 40 independent components vice 64, since lower indices are swappable.
    /// We are keeping the full representation for now for simplicity and clarity with matrix
    /// index conventions, but our code takes advantage of the degeneracy; we don't use all indices.
    components: [f64; 40],
}

/// Eg initialization of a Christoffel field
impl Default for Christoffel {
    fn default() -> Self {
        Self {
            components: [0.; 40],
        }
    }
}

impl Christoffel {
    /// Helper function to calculate a single component. Uses the index convention shown in `from_metric`'s
    /// function description.
    fn calc_component(
        λ: C,
        μ: C,
        ν: C,
        g_this: &MetricTensor,
        g_diffs: &HashMap<(C, PrevNext), &MetricTensor>,
        ds: f64,
    ) -> f64 {
        let mut result = 0.;

        // todo: YOu need to divide by d_metric... whatever sort of quantity that is.
        // todo: Woudl it be computationally-simpler to use a 2-pt formulta instead of midpoint?

        let factor = 1. / (2. * ds);
        let u = Tensor2Config::Uu;

        for σ in &COMPS {
            result += g_this.val(λ, *σ, u)
                * ((g_diffs.get(&(μ, PrevNext::N)).unwrap().val(ν, *σ, u)
                    - g_diffs.get(&(μ, PrevNext::P)).unwrap().val(ν, *σ, u))
                    * factor
                    + (g_diffs.get(&(ν, PrevNext::N)).unwrap().val(μ, *σ, u)
                        - g_diffs.get(&(ν, PrevNext::P)).unwrap().val(μ, *σ, u))
                        * factor
                    - (g_diffs.get(&(*σ, PrevNext::N)).unwrap().val(μ, ν, u)
                        - g_diffs.get(&(*σ, PrevNext::P)).unwrap().val(μ, ν, u))
                        * factor)
        }

        result
    }

    /// Compute from the metric tensor and its numerical derivatives.
    /// Γ^λ _{μ ν} = g^{λ d} ( d_μ g_{ν d} + d_ν g_{b d} - d_d g_{μ ν})
    /// ds is the diff in (what? Spacetime interval? Proper time?) between points. I suppose
    /// it's just the spacing of the grid.
    // pub fn from_metric(metrics: &crate::Arr4dMetric, posit: &Vec4Minkowski, p_i: crate::PositIndex) -> Self {
    // pub fn from_metric(metrics: &Arr4dMetric, p_i: crate::PositIndex, ds: f64) -> Self {
    pub fn from_metric(metrics: &MetricWDiffs, ds: f64) -> Self {
        let mut result = Self {
            components: [0.; 40],
        };

        // todo: Instead of a Hashmap, you could use Val etc methods on MetricWDiffs, but this is fine too.
        // todo: Actually, you should; this HM is redundant.
        let mut metric_map = HashMap::new();
        metric_map.insert((C::T, PrevNext::P), &metrics.t_prev);
        metric_map.insert((C::T, PrevNext::N), &metrics.t_next);
        metric_map.insert((C::X, PrevNext::P), &metrics.x_prev);
        metric_map.insert((C::X, PrevNext::N), &metrics.x_next);
        metric_map.insert((C::Y, PrevNext::P), &metrics.y_prev);
        metric_map.insert((C::Y, PrevNext::N), &metrics.y_next);
        metric_map.insert((C::Z, PrevNext::P), &metrics.z_prev);
        metric_map.insert((C::Z, PrevNext::N), &metrics.z_next);

        // todo: You need a way to get metric's upper config by Comp label!! You're currently using
        // todo the lower-rep version!

        let comps = [C::T, C::X, C::Y, C::Z];

        for λ in &comps {
            // We use this to skip degenerate values, since the lower indices are interchangeable.
            // Let's say we tried μ = t, ν = x
            // Then when this comes up: μ = x, ν = t, we skip it.
            let mut complete_lower_combos = Vec::new();

            for μ in &comps {
                for ν in &comps {
                    // Todo: Is there a more elegant way to do this?
                    let mut degen = false;
                    for complete in &complete_lower_combos {
                        if (ν, μ) == *complete {
                            // Notice the swapped indices compared.
                            degen = true;
                            break;
                        }
                    }
                    if degen {
                        continue;
                    }

                    *result.val_mut(*λ, *μ, *ν) +=
                        Self::calc_component(*λ, *μ, *ν, &metrics.on_pt, &metric_map, ds);

                    complete_lower_combos.push((μ, ν));
                }
            }
        }

        result
    }

    /// Same idea as for the metric tensor.
    /// todo boy: This is a bear.
    /// λ is the upper index.
    /// Note that these indices are arbitrary; they no longer can be composed into a matrix.
    pub fn val(&self, λ: V4Component, μ: V4Component, ν: V4Component) -> f64 {
        let d = &self.components;

        match λ {
            C::T => match μ {
                C::T => match ν {
                    C::T => d[0],
                    C::X => d[1],
                    C::Y => d[2],
                    C::Z => d[3],
                },

                C::X => match ν {
                    C::T => d[1], //  Degen
                    C::X => d[4],
                    C::Y => d[5],
                    C::Z => d[6],
                },

                C::Y => match ν {
                    C::T => d[2], // Degen
                    C::X => d[6], // Degen
                    C::Y => d[7],
                    C::Z => d[8],
                },

                C::Z => match ν {
                    C::T => d[3], // Degen
                    C::X => d[7], // Degen
                    C::Y => d[8], // Degen
                    C::Z => d[9],
                },
            },
            C::X => match μ {
                C::T => match ν {
                    C::T => d[10],
                    C::X => d[11],
                    C::Y => d[12],
                    C::Z => d[13],
                },

                C::X => match ν {
                    C::T => d[11], // Degen
                    C::X => d[14],
                    C::Y => d[15],
                    C::Z => d[16],
                },

                C::Y => match ν {
                    C::T => d[12], // Degen
                    C::X => d[15], // Degen
                    C::Y => d[17],
                    C::Z => d[18],
                },

                C::Z => match ν {
                    C::T => d[13], // Degen
                    C::X => d[16], // Degen
                    C::Y => d[17], // Degen
                    C::Z => d[19],
                },
            },
            C::Y => match μ {
                C::T => match ν {
                    C::T => d[20],
                    C::X => d[21],
                    C::Y => d[22],
                    C::Z => d[23],
                },

                C::X => match ν {
                    C::T => d[21], // Degen
                    C::X => d[24],
                    C::Y => d[25],
                    C::Z => d[26],
                },

                C::Y => match ν {
                    C::T => d[22], // Degen
                    C::X => d[25], // Degen
                    C::Y => d[27],
                    C::Z => d[28],
                },

                C::Z => match ν {
                    C::T => d[23], // Degen
                    C::X => d[26], // Degen
                    C::Y => d[28], // Degen
                    C::Z => d[29],
                },
            },
            C::Z => match μ {
                C::T => match ν {
                    C::T => d[30],
                    C::X => d[31],
                    C::Y => d[32],
                    C::Z => d[33],
                },

                C::X => match ν {
                    C::T => d[31], // Degen
                    C::X => d[34],
                    C::Y => d[35],
                    C::Z => d[36],
                },

                C::Y => match ν {
                    C::T => d[32], // Degen
                    C::X => d[35], // Degen
                    C::Y => d[37],
                    C::Z => d[38],
                },

                C::Z => match ν {
                    C::T => d[33], // Degen
                    C::X => d[36], // Degen
                    C::Y => d[38], // Degen
                    C::Z => d[39],
                },
            },
        }
    }

    /// DRY!
    pub fn val_mut(&mut self, λ: V4Component, μ: V4Component, ν: V4Component) -> &mut f64 {
        // todo: Seeing if this simplifier works. If it does, apply to metric and other tensors a/r
        // &mut self.val(λ, μ, ν)

        // todo: Massive DRY
        let d = &mut self.components;

        match λ {
            C::T => match μ {
                C::T => match ν {
                    C::T => &mut d[0],
                    C::X => &mut d[1],
                    C::Y => &mut d[2],
                    C::Z => &mut d[3],
                },

                C::X => match ν {
                    C::T => &mut d[1], //  Degen
                    C::X => &mut d[4],
                    C::Y => &mut d[5],
                    C::Z => &mut d[6],
                },

                C::Y => match ν {
                    C::T => &mut d[2], //  Degen wit
                    C::X => &mut d[6], // Degen wit
                    C::Y => &mut d[7],
                    C::Z => &mut d[8],
                },

                C::Z => match ν {
                    C::T => &mut d[3], // Degen
                    C::X => &mut d[7], // Degen
                    C::Y => &mut d[8], // Degen
                    C::Z => &mut d[9],
                },
            },
            C::X => match μ {
                C::T => match ν {
                    C::T => &mut d[10],
                    C::X => &mut d[11],
                    C::Y => &mut d[12],
                    C::Z => &mut d[13],
                },

                C::X => match ν {
                    C::T => &mut d[11], // Degen
                    C::X => &mut d[14],
                    C::Y => &mut d[15],
                    C::Z => &mut d[16],
                },

                C::Y => match ν {
                    C::T => &mut d[12], // Degen
                    C::X => &mut d[15], // Degen
                    C::Y => &mut d[17],
                    C::Z => &mut d[18],
                },

                C::Z => match ν {
                    C::T => &mut d[13], // Degen
                    C::X => &mut d[16], // Degen
                    C::Y => &mut d[17], // Degen
                    C::Z => &mut d[19],
                },
            },
            C::Y => match μ {
                C::T => match ν {
                    C::T => &mut d[20],
                    C::X => &mut d[21],
                    C::Y => &mut d[22],
                    C::Z => &mut d[23],
                },

                C::X => match ν {
                    C::T => &mut d[21], // Degen
                    C::X => &mut d[24],
                    C::Y => &mut d[25],
                    C::Z => &mut d[26],
                },

                C::Y => match ν {
                    C::T => &mut d[22], // Degen
                    C::X => &mut d[25], // Degen
                    C::Y => &mut d[27],
                    C::Z => &mut d[28],
                },

                C::Z => match ν {
                    C::T => &mut d[23], // Degen
                    C::X => &mut d[26], // Degen
                    C::Y => &mut d[28], // Degen
                    C::Z => &mut d[29],
                },
            },
            C::Z => match μ {
                C::T => match ν {
                    C::T => &mut d[30],
                    C::X => &mut d[31],
                    C::Y => &mut d[32],
                    C::Z => &mut d[33],
                },

                C::X => match ν {
                    C::T => &mut d[31], // Degen
                    C::X => &mut d[34],
                    C::Y => &mut d[35],
                    C::Z => &mut d[36],
                },

                C::Y => match ν {
                    C::T => &mut d[32], // Degen
                    C::X => &mut d[35], // Degen
                    C::Y => &mut d[37],
                    C::Z => &mut d[38],
                },

                C::Z => match ν {
                    C::T => &mut d[33], // Degen
                    C::X => &mut d[36], // Degen
                    C::Y => &mut d[38], // Degen
                    C::Z => &mut d[39],
                },
            },
        }
    }
}
