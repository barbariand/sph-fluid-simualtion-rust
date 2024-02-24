use crate::fluid_simulation::particle::Particle;
use graphics::*;
use itertools::Itertools;
use opengl_graphics::GlGraphics;
use piston::RenderArgs;

pub struct RenderManager {
    gl: GlGraphics,
}
const BLACK_COLOR: [f32; 4] = [0.0, 0.0, 0.0, 1.0];
impl RenderManager {
    pub fn new(gl: GlGraphics) -> Self {
        RenderManager { gl }
    }

    pub fn render(&mut self, args: &RenderArgs, particles: Vec<Particle>) {
        self.gl.draw(args.viewport(), |c, gl| {
            clear(BLACK_COLOR, gl);

            let draw_state = Default::default();
            particles
                .into_iter()
                .group_by(|v| speed_to_color_gradient(v.speed() as f32))
                .into_iter()
                .map(|(color, particle_group)| {
                    let len = particle_group.size_hint().1.unwrap_or(100) * 10;
                    (
                        color,
                        particle_group.fold(Vec::with_capacity(len), |mut vec, particle| {
                            triangulation::with_ellipse_tri_list(
                                10,
                                c.transform,
                                [
                                    particle.position.x as f64,
                                    particle.position.y as f64,
                                    5.0,
                                    5.0,
                                ],
                                |vertices: &[[f32; 2]]| vec.extend_from_slice(vertices),
                            );

                            vec
                        }),
                    )
                })
                .for_each(|(color, particle_group)| {
                    gl.tri_list(&draw_state, &color, |f| {
                        let (v, e): (&[[[f32; 2]; BACK_END_MAX_VERTEX_COUNT]], &[[f32; 2]]) =
                            particle_group.as_chunks();
                        f(e);

                        v.into_iter().for_each(|val| f(val));
                    });
                });
        });
    }
}
const INVERSED_MAX_SPEED: f32 = 1.0 / 250.0;
const INVERSE_255: f32 = 1.0 / 255.0;
fn speed_to_color_gradient(speed: f32) -> [f32; 4] {
    let ratio: f32 = speed * INVERSED_MAX_SPEED;
    let region = (ratio * 4.0) as i32;

    let x = (region % 256) as f32;
    return match region {
        3 => [1.0, (255.0 - x) * INVERSE_255, 0.0, 1.0],
        2 => [(x as f32) * INVERSE_255, 1.0, 0.0, 1.0],
        1 => [0.0, 1.0, (255.0 - x) * INVERSE_255, 1.0],
        0 => [0.0, x * INVERSE_255, 1.0, 1.0],
        _ => [1.0, 0.0, 0.0, 1.0],
    };
}
