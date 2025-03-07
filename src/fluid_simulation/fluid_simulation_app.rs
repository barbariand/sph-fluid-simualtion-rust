use crate::fluid_simulation::cell_manager::CellManager;
use crate::fluid_simulation::collision_manager::CollisionManager;
use crate::fluid_simulation::external_attractor::ExternalAttractor;
use crate::fluid_simulation::particle::Particle;
use crate::fluid_simulation::particle_dynamics_manager::ParticleDynamicsManager;
use crate::fluid_simulation::smoothed_interaction::SmoothedInteraction;
use crate::piston::PressEvent;
use piston::Button;
use piston::Event;
use piston::Key;
use piston::MouseButton;
use piston::ReleaseEvent;
use piston::UpdateArgs;
use rand::Rng;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;

use vector2d::Vector2D;

pub struct FluidSimulationApp {
    pub particles: Vec<Particle>,
    dynamics_manager: ParticleDynamicsManager,
    smoothed_interaction: SmoothedInteraction,
    external_attractor: ExternalAttractor,
    collision_manager: CollisionManager,
    cell_manager: CellManager,
}

impl FluidSimulationApp {
    pub fn new(box_dimensions: [usize; 2]) -> Self {
        let mut rng = rand::thread_rng();
        let particle_count = 6000;
        let delta_time = 1.0 / 60.0;
        let pressure_multiplier: f32 = 220000.0;
        let target_density: f32 = 0.00002;
        let smoothing_radius: f32 = 14.0;
        let viscosity: f32 = 0.008;
        let particles = (0..particle_count)
            .map(|index| {
                Particle::new(
                    index,
                    Vector2D::new(
                        rng.gen_range(0.0..300_f32),
                        rng.gen_range(0.0..(box_dimensions[1] as f32)),
                    ),
                )
            })
            .collect();
        FluidSimulationApp {
            particles,
            dynamics_manager: ParticleDynamicsManager::new(true, delta_time),
            smoothed_interaction: SmoothedInteraction::new(
                pressure_multiplier,
                target_density,
                smoothing_radius,
                viscosity,
            ),
            external_attractor: ExternalAttractor::new(),
            collision_manager: CollisionManager::new(box_dimensions),
            cell_manager: CellManager::new(particle_count, box_dimensions, smoothing_radius),
        }
    }

    pub fn update(&mut self, _args: &UpdateArgs) {
        let len = self.particles.len();
        self.particles.par_iter_mut().for_each(|particle| {
            self.dynamics_manager.update_position(particle);
            self.collision_manager.apply_boundary_conditions(particle);
        });
        self.cell_manager.update(&mut self.particles);

        let denacc = (0..len)
            .into_par_iter() //rayon
            .map(|particle_index| {
                let adjacente_particles_indices = self
                    .cell_manager
                    .get_adjacent_particles(self.particles[particle_index].position);

                let mut acceleration = self.smoothed_interaction.calculate_pressure(
                    particle_index,
                    adjacente_particles_indices,
                    &self.particles,
                );
                acceleration += self.smoothed_interaction.calculate_viscosity(
                    particle_index,
                    adjacente_particles_indices,
                    &self.particles,
                );

                (
                    acceleration,
                    self.smoothed_interaction.calculate_density(
                        particle_index,
                        adjacente_particles_indices,
                        &self.particles,
                    ),
                ) // calulate as before but dont apply it
            })
            .collect::<Vec<_>>();

        self.particles
            .par_iter_mut()
            .zip(denacc.into_par_iter())
            .for_each(|(particle, (acceleration, dencity))| {
                particle.acceleration = acceleration
                    + self
                        .external_attractor
                        .get_external_attraction_acceleration(particle);
                particle.local_density = dencity;
                self.dynamics_manager.update_velocity(particle);
            }); // apply it
    }

    pub fn handle_event(&mut self, event: Event) {
        if let Some(Button::Keyboard(Key::G)) = event.press_args() {
            self.dynamics_manager.toggle_gravity();
        }
        if let Some(Button::Keyboard(Key::D)) = event.press_args() {
            self.collision_manager.break_dam();
        }
        if let Some(Button::Mouse(MouseButton::Left)) = event.press_args() {
            self.external_attractor
                .activate(Vector2D::new(400.0_f32, 100.0_f32));
        }
        if let Some(Button::Mouse(MouseButton::Left)) = event.release_args() {
            self.external_attractor.active = false;
        }
    }
}
