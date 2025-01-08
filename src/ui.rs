use egui::{self, Color32, RichText};
use graphics::{EngineUpdates, Scene};

use crate::State;

const UI_WIDTH: f32 = 300.;
const SIDE_PANEL_SIZE: f32 = 400.;

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn ui_handler(state: &mut State, cx: &egui::Context, scene: &mut Scene) -> EngineUpdates {
    let mut engine_updates = EngineUpdates::default();

    let panel = egui::SidePanel::left(0) // ID must be unique among panels.
        .default_width(SIDE_PANEL_SIZE);

    panel.show(cx, |ui| {
        engine_updates.ui_size = ui.available_width();
        ui.spacing_mut().item_spacing = egui::vec2(10.0, 12.0);

        ui.set_max_width(UI_WIDTH);
    });

    engine_updates
}
