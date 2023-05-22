use core::f32::consts::TAU;

use crate::{
    tensors::{Event, Vec4Minkowski},
    State,
};

use graphics::{
    self, Camera, ControlScheme, DeviceEvent, EngineUpdates, Entity, InputSettings, LightType,
    Lighting, Mesh, PointLight, Scene, UiLayout, UiSettings,
};

use crate::tensors::Worldline;
use lin_alg2::{
    f32::{Quaternion, Vec3},
    f64::Vec3 as Vec3F64,
};

type Color = (f32, f32, f32);

const WINDOW_TITLE: &str = "Gravity";
const WINDOW_SIZE_X: f32 = 1_600.;
const WINDOW_SIZE_Y: f32 = 1_200.;
const RENDER_DIST: f32 = 100.;
const BACKGROUND_COLOR: Color = (0.5, 0.5, 0.5);

const SCHWARZ_SPHERE_SIZE: f32 = 0.5;
const GEODESIC_SPHERE_SIZE: f32 = 0.1;

const COLOR_MASS: Color = (0., 0., 0.);

const WORLDLINE_COLORS: [Color; 10] = [
    (0., 0., 1.),
    (0., 0.5, 0.2),
    (1., 0., 0.),
    (0., 0.5, 0.5),
    (0., 0.3, 0.8),
    (0.5, 0.5, 0.),
    (0.33, 0.33, 0.333),
    (0.5, 0., 0.5),
    (0.5, 0.4, 0.2),
    (0.5, 0.4, 0.3),
];

const SHINYNESS: f32 = 10.;

fn event_handler(
    _state: &mut State,
    _event: DeviceEvent,
    _scene: &mut Scene,
    _dt: f32,
) -> EngineUpdates {
    // match event {
    //     DeviceEvent::Key(key) => {}
    //     _ => (),
    // }
    EngineUpdates::default()
}

/// This runs each frame. Currently, no updates.
fn render_handler(_state: &mut State, _scene: &mut Scene, _dt: f32) -> EngineUpdates {
    // EngineUpdates::default()

    EngineUpdates {
        // compute: true,
        ..Default::default()
    }
}

pub fn update_meshes(z_displayed: f64, scene: &mut Scene) {
    // Our meshes are defined in terms of a start point,
    // and a step. Adjust the step to center the grid at
    // the renderer's center.
    // const SFC_MESH_START: f32 = -4.;
    // let sfc_mesh_start = grid_min as f32; // todo: Sync graphics and atomic coords?
    // let sfc_mesh_step: f32 = -2. * sfc_mesh_start / N as f32;

    // `z_displayed` is a value float. Convert this to an index. Rounds to the nearest integer.
    // let z_i = map_linear(z_displayed, (grid_min, grid_max), (0., N as f64)) as usize;
    // todo: Using your new system, you can't use a linear map here!

    // let mut z_i = 0;
    // for i in 0..grid_posits.len() {
    //     if grid_posits[0][0][i].z > z_displayed {
    //         z_i = i;
    //         break;
    //     }
    // }
    //
    let mut meshes = Vec::new();
    //
    // meshes.push(Mesh::new_surface(
    //     // todo: Be able to show shared V and per-elec. Currently hard-coded.
    //     &prepare_2d_mesh_real(grid_posits, &surfaces.V, z_i, 1., grid_n),
    //     true,
    // ));
    //
    // meshes.push(Mesh::new_surface(
    //     &prepare_2d_mesh(
    //         grid_posits,
    //         &surfaces.psi.on_pt,
    //         z_i,
    //         PSI_SCALER,
    //         mag_phase,
    //         false,
    //         grid_n,
    //     ),
    //     true,
    // ));
    //
    // meshes.push(Mesh::new_surface(
    //     &prepare_2d_mesh(
    //         grid_posits,
    //         &surfaces.psi.on_pt,
    //         z_i,
    //         PSI_SCALER,
    //         mag_phase,
    //         true,
    //         grid_n,
    //     ),
    //     true,
    // ));
    //
    // let mut psi_sq = new_data_real(grid_n);
    // util::norm_sq(&mut psi_sq, &surfaces.psi.on_pt, grid_n);
    //
    // meshes.push(Mesh::new_surface(
    //     &prepare_2d_mesh_real(grid_posits, &psi_sq, z_i, PSI_SQ_SCALER, grid_n),
    //     true,
    // ));
    //
    // for (scaler, sfc) in [
    //     (PSI_PP_SCALER, &surfaces.psi_pp_calculated),
    //     (PSI_PP_SCALER, &surfaces.psi_pp_measured),
    // ] {
    //     meshes.push(Mesh::new_surface(
    //         &prepare_2d_mesh(grid_posits, sfc, z_i, scaler, mag_phase, false, grid_n),
    //         true,
    //     ));
    //
    //     meshes.push(Mesh::new_surface(
    //         &prepare_2d_mesh(grid_posits, sfc, z_i, scaler, mag_phase, true, grid_n),
    //         true,
    //     ));
    // }

    // meshes.push(Mesh::new_surface(
    //     // &prepare_2d_mesh(&surfaces.grid_posits, &surfaces.aux2, z_i, ELEC_V_SCALER),
    //     &prepare_2d_mesh(&surfaces.grid_posits, &surfaces.aux2, z_i, 10.),
    //     // todo: Center! Maybe offset in entities.
    //     // sfc_mesh_start,
    //     // sfc_mesh_step,
    //     true,
    // ));

    meshes.push(Mesh::new_sphere(SCHWARZ_SPHERE_SIZE, 14, 14));
    meshes.push(Mesh::new_sphere(GEODESIC_SPHERE_SIZE, 8, 8));

    scene.meshes = meshes;
}

pub fn update_entities(schwarz_posit: Vec3F64, worldlines: &[Worldline], scene: &mut Scene) {
    let mut entities = Vec::new();

    entities.push(Entity::new(
        0,
        Vec3::new(
            schwarz_posit.x as f32,
            // We invert Y and Z due to diff coord systems
            // between the meshes and the renderer.
            schwarz_posit.z as f32,
            schwarz_posit.y as f32,
        ),
        Quaternion::new_identity(),
        1.,
        COLOR_MASS,
        SHINYNESS,
    ));

    for (i, worldline) in worldlines.iter().enumerate() {
        for event in &worldline.events {
            entities.push(Entity::new(
                0,
                Vec3::new(
                    event.posit.t() as f32,
                    // We invert Y and Z due to diff coord systems
                    // between the meshes and the renderer.
                    event.posit.x() as f32,
                    event.posit.y() as f32,
                    // z left out for now
                ),
                Quaternion::new_identity(),
                1.,
                WORLDLINE_COLORS[i],
                SHINYNESS,
            ));
        }
    }

    scene.entities = entities;
}

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

    update_meshes(state.z_displayed, &mut scene);
    update_entities(state.schwarz_mass.0, &state.worldlines, &mut scene);

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
