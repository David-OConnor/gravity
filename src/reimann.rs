//! This module contains code for the Reimann and Ricci tensors.

use crate::{
    christoffel::{Christoffel, ChristoffelWDiffs},
    metric::MetricTensor,
    tensors::{MetricTensor, PrevNext, ReimannConfig, Tensor2Config, V4Component, C, COMPS},
    Arr4dChristoffel, Arr4dMetric,
};

/// Reimann tensor. ρ
pub struct Riemann {
    /// Despite this being a (class?)-4 tensor, it's only 20 components long due to most being
    /// degenerate.
    /// Order is, starting top left, main diag, diag shifted down/left one, ripple.
    components_ulll: [f64; 20],
    components_llll: [f64; 20],
}

impl Riemann {
    /// Helper function to calculate a single component. Uses the index convention shown in `from_metric`'s
    /// function description.
    fn calc_component(ρ: C, σ: C, μ: C, ν: C, Γs: &ChristoffelWDiffs, ds: f64) -> f64 {
        let mut result = 0.;

        // todo: Woudl it be computationally-simpler to use a 2-pt formulta instead of midpoint?

        let factor = 1. / (2. * ds);

        for λ in &COMPS {
            let part_1 = (Γs.val(μ, PrevNext::N).val(ρ, ν, σ)
                - Γs.val(μ, PrevNext::P).val(ρ, ν, σ))
                * factor;

            let part_2 = (Γs.val(ν, PrevNext::N).val(ρ, μ, σ)
                - Γs.val(ν, PrevNext::P).val(ρ, μ, σ))
                * factor;

            // C::T is a dummy val here.
            let part_3 = Γs.val(C::T, PrevNext::OnPt).val(ρ, μ, λ) * Γ.val(λ, ν, σ) * factor;
            let part_4 = Γs.val(C::T, PrevNext::OnPt).val(ρ, ν, λ) * Γ.val(λ, μ, σ) * factor;

            result += part_1 - part_2 + part_3 - part_4;
        }

        result
    }

    ///R^ρ_σμν = d_μ Γ^ρ_νσ - d_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
    /// This assumes we've pre-calculated Christoffels from the metric.
    pub fn from_christoffel(Γs: &ChristoffelWDiffs, ds: f64) -> Self {
        // todo: Make sure is aren't at the edges. If so, return etc.
        let mut result = Self {
            components_ulll: [0.; 20],
            components_llll: [0.; 20],
        };

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
                            Self::calc_component(*λ, *σ, *μ, *ν, Γs, ds);

                        // complete_lower_combos.push((μ, ν));
                    }
                }
            }
        }

        result
    }

    /// This calculates Christoffels from the metric, then calls our other constructor using them.
    pub fn from_metric(metrics: &Arr4dMetric, p_i: crate::PositIndex, ds: f64) -> Self {
        let Γs = Christoffel::from_metric(metrics, p_i, ds);

        Self::from_christoffel(&Γs, p_i, ds)
    }

    /// ρ is the upper index.
    /// Note that these indices are arbitrary; they no longer can be composed into a matrix.
    /// See *A General Relativity Workbook*, page 226.
    pub fn val(
        &self,
        ρ: V4Component,
        σ: V4Component,
        μ: V4Component,
        ν: V4Component,
        config: ReimannConfig,
    ) -> f64 {
        let d = &self.components_llll;

        // todo: This is likely valid for llll only.

        // todo: Is this right?
        // if ρ == σ || μ == ν {
        //     return 0.;
        // }

        // Before this, you should probably ahve some manipulation code to coerce from 2 forms to
        // the ones we're comparing to below.

        match (ρ, σ, μ, ν) {
            // Main diagonal
            (C::T, C::X, C::T, C::X) => d[0],
            (C::T, C::Y, C::T, C::Y) => d[1],
            (C::T, C::Z, C::T, C::Z) => d[2],
            (C::X, C::Y, C::X, C::Y) => d[3],
            (C::X, C::Z, C::X, C::Z) => d[4],
            (C::Y, C::Z, C::Y, C::Z) => d[5],

            // Shifted down and left 1
            (C::T, C::Y, C::T, C::X) => d[6],
            (C::T, C::Z, C::T, C::Y) => d[7],
            (C::X, C::Y, C::T, C::Z) => d[8],
            (C::X, C::Z, C::X, C::Y) => d[9],
            (C::Y, C::Z, C::X, C::Z) => d[10],

            (C::T, C::Z, C::T, C::X) => d[11],
            (C::X, C::Y, C::T, C::Y) => d[12],
            (C::X, C::Z, C::T, C::Z) => d[13],
            (C::Y, C::Z, C::X, C::Y) => d[14],

            (C::X, C::Y, C::T, C::X) => d[15],
            (C::X, C::Z, C::T, C::Y) => d[16],
            (C::Y, C::Z, C::T, C::Z) => d[17],

            (C::X, C::Z, C::T, C::X) => d[18],
            (C::Y, C::Z, C::T, C::Y) => d[19],

            // todo: The cyclic relation means we can elimintate the 21st component; fix this;
            // todo there is no index 20.
            // (C::Y, C::Z, C::T, C::X) => d[20],
            (C::Y, C::Z, C::T, C::X) => 69., //??

            // todo: Is there a less-repetative way of showing this symmetry than manually
            // todo assigning the upper indices?

            // Shifted up and right 1.
            (C::T, C::X, C::T, C::Y) => d[6],
            (C::T, C::Y, C::T, C::Z) => d[7],
            (C::T, C::Z, C::X, C::Y) => d[8],
            (C::X, C::Y, C::X, C::Z) => d[9],
            (C::X, C::Z, C::Y, C::Z) => d[10],

            (C::T, C::X, C::T, C::Z) => d[11],
            (C::T, C::Y, C::X, C::Y) => d[12],
            (C::T, C::Z, C::X, C::Z) => d[13],
            (C::X, C::Y, C::Y, C::Z) => d[14],

            (C::T, C::X, C::X, C::Y) => d[15],
            (C::T, C::Y, C::X, C::Z) => d[16],
            (C::T, C::Z, C::Y, C::Z) => d[17],

            (C::T, C::X, C::X, C::Z) => d[18],
            (C::T, C::Y, C::Y, C::Z) => d[19],

            (C::T, C::X, C::Y, C::Z) => 69., //??

            _ => 0.,
        }
    }
}

/// Ricci tensor
#[derive(Default)]
pub struct Ricci {
    /// 16-components due to order-2 tensor; only 6 are independent.
    /// todo: 6, or 10?
    /// todo: We currently use 16 components so we can re-use the Metric tensor's operaations. This should
    /// todo be fine, but long-term, switch to 6 or 10. Either re-implement the get/set, or use a general get/set algo called
    /// todo by both.
    components_ll: [f64; 16],
}

impl Ricci {
    /// Contract over the final covariant index to create the Ricci tensor from the Reimann one.
    /// http://astro.dur.ac.uk/~done/gr/l11.pdf
    /// Note: This contracts over the third index. An alternate convention is contracting over the fourth.
    pub fn from_reimann(reimann: &Reimann) -> Self {
        let mut result = Self::default();

        for ρ in &COMPS {
            for σ in &COMPS {
                for ν in &COMPS {
                    result.val_mut(*σ, *μ) += reimann.val(ρ, ρ, μ, ν);
                }
            }
        }

        result
    }

    pub fn make_scaler(&self, metric: &MetricTensor) -> f64 {
        let mut result = 0.;

        for β in &COMPS {
            for ν in &COMPS {
                result += metric.val(*β, *ν, Tensor2Config::Uu) * self.val(*β, *ν);
            }
        }
        result
    }

    /// We piggyback on `Metric`'s getters and setters for now. Note that Ricci is symmetric,
    /// while Metric isn't always symmetric.
    /// todo: How does it reduce to 6 components? Maybe that's not true.
    pub fn val(&self, μ: V4Component, ν: V4Component) -> f64 {
        let m = MetricTensor {
            components_ll: self.components_ll,
            components_uu: [0.; 16],
        };

        m.val(μ, ν, Tensor2Config::Ll)
    }

    pub fn val_mut(&mut self, μ: V4Component, ν: V4Component) -> &mut f64 {
        let mut m = MetricTensor {
            components_ll: self.components_ll,
            components_uu: [0.; 16],
        };
        m.val_mut(μ, ν, Tensor2Config::Ll)
    }
}
