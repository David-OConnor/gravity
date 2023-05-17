//! This module contains code for the Reimann and Ricci tensors.

use crate::{
    tensors::V4Component,
    C,
};


/// Reimann tensor. ρ
pub struct Riemann {
    /// Despite this being a (class?)-4 tensor, it's only 20 components long due to most being
    /// degenerate.
    components: [f64; 20],
}

impl Riemann {
    /// ρ is the upper index.
    /// Note that these indices are arbitrary; they no longer can be composed into a matrix.
    pub fn val(&self, ρ: V4Component, σ: V4Component, μ: V4Component, ν: V4Component, ) -> f64 {
        let d = &self.components;

        match λ {
            C::T => match μ {
                C::T => match ν {
                    C::T => d[0],
                    C::X => d[1],
                    C::Y => d[2],
                    C::Z => d[3],
                },
            }
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
    pub fn val(&self, μ: V4Component, ν: V4Component, ) -> f64 {
        let d = &self.components;

        match λ {
            C::T => match μ {
                C::T => match ν {
                    C::T => d[0],
                    C::X => d[1],
                    C::Y => d[2],
                    C::Z => d[3],
                },
            }
        }
    }

}
