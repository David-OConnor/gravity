use core::f32::consts::TAU;

use crate::State;

use graphics::{
    self, Camera, ControlScheme, DeviceEvent, EngineUpdates, Entity, InputSettings, LightType,
    Lighting, Mesh, PointLight, Scene, UiLayout, UiSettings,
};

use lin_alg2::{
    f32::{Quaternion, Vec3},
    f64::Vec3 as Vec3F64,
};

type Color = (f32, f32, f32);

const WINDOW_TITLE: &str = "ψ lab";
const WINDOW_SIZE_X: f32 = 1_600.;
const WINDOW_SIZE_Y: f32 = 1_200.;
const RENDER_DIST: f32 = 100.;
const BACKGROUND_COLOR: Color = (0.5, 0.5, 0.5);

/// Entry point to our render and event loop.
pub fn render(state: State) {
    let mut scene = Scene {
        meshes: Vec::new(),   // updated below.
        entities: Vec::new(), // updated below.
        camera: Camera {
            fov_y: TAU / 8.,
            position: Vec3::new(0., 6., -15.),
            far: RENDER_DIST,
            orientation: Quaternion::from_axis_angle(Vec3::new(1., 0., 0.), TAU / 16.),
            ..Default::default()
        },
        lighting: Lighting {
            ambient_color: [-1., 1., 1., 0.5],
            ambient_intensity: 0.03,
            point_lights: vec![
                // Light from above. The sun?
                PointLight {
                    type_: LightType::Omnidirectional,
                    position: Vec3::new(0., 100., 0.),
                    diffuse_color: [0.6, 0.4, 0.3, 1.],
                    specular_color: [0.6, 0.4, 0.3, 1.],
                    diffuse_intensity: 10_000.,
                    specular_intensity: 10_000.,
                },
                PointLight {
                    type_: LightType::Omnidirectional,
                    position: Vec3::new(30., 100., 30.),
                    diffuse_color: [0.3, 0.4, 0.5, 1.],
                    specular_color: [0.3, 0.4, 0.5, 1.],
                    diffuse_intensity: 10_000.,
                    specular_intensity: 10_000.,
                },
            ],
        },
        background_color: BACKGROUND_COLOR,
        window_size: (WINDOW_SIZE_X, WINDOW_SIZE_Y),
        window_title: WINDOW_TITLE.to_owned(),
    };

    // update_meshes(
    //     &state.surfaces,
    //     state.z_displayed,
    //     &mut scene,
    //     // state.grid_min,
    //     // state.grid_max,
    // );
    // update_entities(&state.charges, &state.show_surfaces, &mut scene);

    let input_settings = InputSettings {
        initial_controls: ControlScheme::FreeCamera,
        ..Default::default()
    };
    let ui_settings = UiSettings {
        layout: UiLayout::Left,
        // todo: How to handle this? For blocking keyboard and moues inputs when over the UI.
        // width: gui::UI_WIDTH as f64, // todo: Not working correctly.
        size: 0., // todo: Bad API here.
        icon_path: None,
    };

    graphics::run(
        state,
        scene,
        input_settings,
        ui_settings,
        render_handler,
        event_handler,
        crate::ui::ui_handler,
        include_str!("shader_compute.wgsl").into(),
    );
}
