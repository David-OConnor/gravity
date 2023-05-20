//! This module contains code for the Reimann and Ricci tensors.

use std::collections::HashMap;

use crate::tensors::C;
use crate::{
    christoffel::Christoffel,
    tensors::{MetricTensor, Tensor2Config, V4Component, COMPS},
    Arr4dChristoffel, Arr4dMetric, C,
};

/// Reimann tensor. ρ
pub struct Riemann {
    /// Despite this being a (class?)-4 tensor, it's only 20 components long due to most being
    /// degenerate.
    components: [f64; 20],
}

impl Riemann {
    /// Helper function to calculate a single component. Uses the index convention shown in `from_metric`'s
    /// function description.
    fn calc_component(
        ρ: C,
        σ: C,
        μ: C,
        ν: C,
        Γ_this: &Christoffel,
        Γ_diffs: &HashMap<(C, PrevNext), &MetricTensor>,
        ds: f64,
    ) -> f64 {
        let mut result = 0.;

        // todo: Woudl it be computationally-simpler to use a 2-pt formulta instead of midpoint?

        let factor = 1. / (2. * ds);

        for λ in &COMPS {
            let part_1 = (Γ_diffs.get(&(μ, PrevNext::N)).unwrap().val(ρ, *ν, *σ)
                - Γ_diffs.get(&(μ, PrevNext::P)).unwrap().val(ρ, *ν, *σ))
                * factor;

            let part_2 = (Γ_diffs.get(&(ν, PrevNext::N)).unwrap().val(ρ, *μ, *σ)
                - Γ_diffs.get(&(ν, PrevNext::P)).unwrap().val(ρ, *μ, *σ))
                * factor;

            let part_3 = Γ.val(ρ, μ, λ) * Γ.val(λ, ν, σ) * factor;
            let part_4 = Γ.val(ρ, ν, λ) * Γ.val(λ, μ, σ) * factor;

            result += part_1 - part_2 + part_3 - part_4;
        }

        result
    }

    ///R^ρ_σμν = d_μ Γ^ρ_νσ - d_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
    /// This assumes we've pre-calculated Christoffels from the metric.
    pub fn from_christoffel(Γs: &Arr4dChristoffel, p_i: crate::PositIndex, ds: f64) -> Self {
        // todo: Make sure is aren't at the edges. If so, return etc.

        // todo: Figure out which components you can eliminate or combine.

        let Γ_t_prev = &Γs[p_i.t - 1][p_i.x][p_i.y][p_i.z];
        let Γ_t_next = &Γs[p_i.t - 1][p_i.x][p_i.y][p_i.z];
        let Γ_x_prev = &Γs[p_i.t][p_i.x - 1][p_i.y][p_i.z];
        let Γ_x_next = &Γs[p_i.t][p_i.x + 1][p_i.y][p_i.z];
        let Γ_y_prev = &Γs[p_i.t][p_i.x][p_i.y - 1][p_i.z];
        let Γ_y_next = &Γs[p_i.t][p_i.x][p_i.y + 1][p_i.z];
        let Γ_z_prev = &Γs[p_i.t][p_i.x][p_i.y][p_i.z - 1];
        let Γ_z_next = &Γs[p_i.t][p_i.x][p_i.y][p_i.z + 1];

        let mut christoffels = HashMap::new();
        christoffels.insert((C::T, PrevNext::P), Γ_t_prev);
        christoffels.insert((C::T, PrevNext::N), Γ_t_next);
        christoffels.insert((C::X, PrevNext::P), Γ_x_prev);
        christoffels.insert((C::X, PrevNext::N), Γ_x_next);
        christoffels.insert((C::Y, PrevNext::P), Γ_y_prev);
        christoffels.insert((C::Y, PrevNext::N), Γ_y_next);
        christoffels.insert((C::Z, PrevNext::P), Γ_z_prev);
        christoffels.insert((C::Z, PrevNext::N), Γ_z_next);

        for ρ in &comps {
            // We use this to skip degenerate values, since the lower indices are interchangeable.
            // Let's say we tried μ = t, ν = x
            // Then when this comes up: μ = x, ν = t, we skip it.
            let mut complete_lower_combos = Vec::new();

            for σ in &comps {
                for μ in &comps {
                    for ν in &comps {
                        let mut degen = false;
                        // for complete in &complete_lower_combos {
                        //     if (ν, μ) == *complete {
                        //         // Notice the swapped indices compared.
                        //         degen = true;
                        //         break;
                        //     }
                        // }
                        // if degen {
                        //     continue;
                        // }

                        *result.val_mut(*ρ, *σ, *μ, *ν) +=
                            Self::calc_component(*λ, *σ, *μ, *ν, &Γ_this, &christoffels, ds);

                        // complete_lower_combos.push((μ, ν));
                    }
                }
            }
        }
    }

    /// This calculates Christoffels from the metric, then calls our other constructor using them.
    pub fn from_metric(metrics: &Arr4dMetric, p_i: crate::PositIndex, ds: f64) -> Self {
        let Γs = Christoffel::from_metric(metrics, p_i, ds);

        Self::from_christoffel(&Γs, p_i, ds)
    }

    /// ρ is the upper index.
    /// Note that these indices are arbitrary; they no longer can be composed into a matrix.
    pub fn val(&self, ρ: V4Component, σ: V4Component, μ: V4Component, ν: V4Component) -> f64 {
        let d = &self.components;

        match λ {
            C::T => match μ {
                C::T => match ν {
                    C::T => d[0],
                    C::X => d[1],
                    C::Y => d[2],
                    C::Z => d[3],
                },
            },
        }
    }
}

/// Ricci tensor
pub struct Ricci {
    /// 16-components due to order-2 tensor; only 6 are independent.
    components: [f64; 6],
}

impl Ricci {
    /// Note that these indices are arbitrary; they no longer can be composed into a matrix.
    pub fn val(&self, μ: V4Component, ν: V4Component) -> f64 {
        let d = &self.components;

        match λ {
            C::T => match μ {
                C::T => match ν {
                    C::T => d[0],
                    C::X => d[1],
                    C::Y => d[2],
                    C::Z => d[3],
                },
            },
        }
    }
}
