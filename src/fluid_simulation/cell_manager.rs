use crate::fluid_simulation::particle::Particle;
use vector2d::Vector2D;

pub struct CellManager {
    particle_count: usize,
    pub spatial_lookup: Vec<(usize, usize)>,
    pub starting_indices: Vec<usize>,
    number_of_columns: usize,
    number_of_rows: usize,
    cell_size: f32,
    number_of_cells: usize,
}

impl CellManager {
    pub fn new(particle_count: usize, box_dimensions: [usize; 2], smoothing_radius: f32) -> Self {
        let cell_size = 2.0 * smoothing_radius;
        let number_of_columns = (box_dimensions[0] as f32 / cell_size).ceil() as usize;
        let number_of_rows = (box_dimensions[1] as f32 / cell_size).ceil() as usize;
        let number_of_cells = (number_of_columns * number_of_rows) as usize;
        CellManager {
            particle_count,
            spatial_lookup: (0..particle_count)
                .map(|_| (number_of_cells as usize, 0))
                .collect(),
            starting_indices: (0..number_of_cells)
                .map(|_| number_of_cells as usize)
                .collect(),
            number_of_columns,
            number_of_rows,
            cell_size,
            number_of_cells,
        }
    }

    pub fn update(&mut self, particles: &mut Vec<Particle>) {
        self.spatial_lookup = (0..self.particle_count)
            .map(|_| (self.number_of_cells as usize, 0))
            .collect();
        for particle in particles {
            self.to_spacial_lookup(particle)
        }
        self.spatial_lookup.sort_by(|s_a, s_b| s_a.0.cmp(&s_b.0));
        self.generate_start_indices();
    }

    pub fn get_adjacent_particles(&self, particle_position: Vector2D<f32>) -> Vec<usize> {
        let mut adjacent_particle_indexes: Vec<usize> = vec![];
        let adjacent_cell_keys = self.get_adjacent_cell_keys_from_position(particle_position);
        for adjacent_cell_key in adjacent_cell_keys {
            let particle_indexes = self.get_particle_indexes_from_cell(adjacent_cell_key);
            adjacent_particle_indexes.extend(particle_indexes)
        }
        adjacent_particle_indexes
    }

    fn to_spacial_lookup(&mut self, particle: &mut Particle) {
        let cell_coord = self.particle_position_to_cell_coord(particle.position);
        let cell_key = self.cell_coord_to_cell_key(cell_coord);
        particle.cell_key = cell_key;
        self.spatial_lookup[particle.id as usize] = (cell_key, particle.id)
    }

    fn generate_start_indices(&mut self) {
        let mut starting_indices: Vec<usize> =
            vec![self.particle_count as usize; self.number_of_cells as usize];
        for (sl_index, &(cell_key, _)) in self.spatial_lookup.iter().enumerate() {
            if starting_indices[cell_key] == self.particle_count as usize {
                starting_indices[cell_key] = sl_index;
            }
        }
        self.starting_indices = starting_indices;
    }

    pub fn get_adjacent_cell_keys_from_position(&self, position: Vector2D<f32>) -> Vec<usize> {
        let current_cell_coord = self.particle_position_to_cell_coord(position);

        let mut adjacent_cell_coords = vec![
            current_cell_coord + Vector2D::new(-1, -1),
            current_cell_coord + Vector2D::new(-1, 0),
            current_cell_coord + Vector2D::new(-1, 1),
            current_cell_coord + Vector2D::new(0, -1),
            current_cell_coord,
            current_cell_coord + Vector2D::new(0, 1),
            current_cell_coord + Vector2D::new(1, -1),
            current_cell_coord + Vector2D::new(1, 0),
            current_cell_coord + Vector2D::new(1, 1),
        ];

        adjacent_cell_coords = adjacent_cell_coords
            .iter()
            .cloned()
            .filter(|coord| {
                coord.x >= 0
                    && coord.x < self.number_of_columns
                    && coord.y >= 0
                    && coord.y < self.number_of_rows
            })
            .collect();
        let adjacent_cell_keys = adjacent_cell_coords
            .iter()
            .map(|coord| self.cell_coord_to_cell_key(*coord))
            .collect();
        adjacent_cell_keys
    }

    pub fn get_particle_indexes_from_cell(&self, cell_key: usize) -> Vec<usize> {
        let mut particle_indexes: Vec<usize> = Vec::new();
        let mut spatial_lookup_cell: usize = cell_key;
        let mut spatial_lookup_index = self.starting_indices[cell_key];
        if spatial_lookup_index >= self.particle_count as usize {
            return particle_indexes;
        }
        while cell_key == spatial_lookup_cell {
            let particle_index = self.spatial_lookup[spatial_lookup_index].1;
            particle_indexes.push(particle_index);
            spatial_lookup_index += 1;
            if spatial_lookup_index >= self.particle_count as usize {
                break;
            }
            spatial_lookup_cell = self.spatial_lookup[spatial_lookup_index].0;
        }
        particle_indexes
    }

    pub fn particle_position_to_cell_coord(&self, position: Vector2D<f32>) -> Vector2D<usize> {
        let x = (position.x / self.cell_size).floor() as usize;
        let y = (position.y / self.cell_size).floor() as usize;
        Vector2D::new(x, y)
    }

    pub fn cell_coord_to_cell_key(&self, coord: Vector2D<usize>) -> usize {
        ((coord.x * self.number_of_rows) + coord.y) as usize
    }
}
