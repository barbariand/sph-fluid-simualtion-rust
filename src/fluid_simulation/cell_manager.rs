use crate::fluid_simulation::particle::Particle;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use vector2d::Vector2D;

pub struct CellManager {
    particle_count: usize,
    pub spatial_lookup: Vec<(usize, usize)>,
    pub starting_indices: Vec<usize>,
    number_of_columns: usize,
    number_of_rows: usize,
    cell_size: f32,
    number_of_cells: usize,
    adjacencies_cache: Vec<Vec<usize>>,
}

impl CellManager {
    pub fn new(particle_count: usize, box_dimensions: [usize; 2], smoothing_radius: f32) -> Self {
        let cell_size = 2.0 * smoothing_radius;
        let number_of_columns = (box_dimensions[0] as f32 / cell_size).ceil() as usize;
        let number_of_rows = (box_dimensions[1] as f32 / cell_size).ceil() as usize;
        let number_of_cells = number_of_columns * number_of_rows;
        CellManager {
            particle_count,
            spatial_lookup: (0..particle_count).map(|_| (number_of_cells, 0)).collect(),
            starting_indices: (0..number_of_cells).collect(),
            number_of_columns,
            number_of_rows,
            cell_size,
            adjacencies_cache: Vec::with_capacity(number_of_cells),
            number_of_cells,
        }
    }

    pub fn update(&mut self, particles: &mut Vec<Particle>) {
        self.spatial_lookup = vec![(self.number_of_cells, 0); self.particle_count];

        for particle in particles.iter_mut() {
            self.to_spacial_lookup(particle)
        }
        self.rebuild_adjacencies_cache(particles);
        self.spatial_lookup.sort_by(|s_a, s_b| s_a.0.cmp(&s_b.0));
        self.generate_start_indices();
    }

    pub fn get_adjacent_particles(&self, particle_position: Vector2D<f32>) -> &Vec<usize> {
        let cell_key =
            self.cell_coord_to_cell_key(self.particle_position_to_cell_coord(particle_position));
        &self.adjacencies_cache[cell_key]
    }

    fn to_spacial_lookup(&mut self, particle: &mut Particle) {
        let cell_coord = self.particle_position_to_cell_coord(particle.position);
        let cell_key = self.cell_coord_to_cell_key(cell_coord);
        particle.cell_key = cell_key;
        self.spatial_lookup[particle.id] = (cell_key, particle.id)
    }

    fn generate_start_indices(&mut self) {
        let mut starting_indices: Vec<usize> = vec![self.particle_count; self.number_of_cells];
        for (sl_index, &(cell_key, _)) in self.spatial_lookup.iter().enumerate() {
            if starting_indices[cell_key] == self.particle_count {
                starting_indices[cell_key] = sl_index;
            }
        }
        self.starting_indices = starting_indices;
    }
    fn rebuild_adjacencies_cache(&mut self, particles: &[Particle]) {
        let mut temp_cache = vec![vec![]; self.number_of_cells];

        // Step 1: Assign particles to their cells
        for particle in particles {
            let cell_key = self
                .cell_coord_to_cell_key(self.particle_position_to_cell_coord(particle.position));
            temp_cache[cell_key].push(particle.id);
        }

        // Step 2: Create a comprehensive adjacencies cache
        self.adjacencies_cache = vec![vec![]; self.number_of_cells];
        for cell_key in 0..self.number_of_cells {
            let adj_keys = self.get_adjacent_cell_keys(cell_key); // Implement this method
            let mut all_adjacent_particles = Vec::new();

            // Include current cell's particles
            all_adjacent_particles.extend(&temp_cache[cell_key]);

            // Include adjacent cells' particles
            for adj_key in adj_keys {
                if let Some(particles) = temp_cache.get(adj_key) {
                    all_adjacent_particles.extend(particles);
                }
            }

            // Deduplicate the particles IDs for each cell
            all_adjacent_particles.sort_unstable();
            all_adjacent_particles.dedup();

            // Update the cache with particles from adjacent cells included
            self.adjacencies_cache[cell_key] = all_adjacent_particles;
        }
    }

    // Helper method to get cell keys of all adjacent cells including itself
    fn get_adjacent_cell_keys(&self, cell_key: usize) -> Vec<usize> {
        let (cx, cy) = self.cell_key_to_coord(cell_key);
        let mut keys = Vec::new();

        for dx in -1..=1 {
            for dy in -1..=1 {
                let nx = cx as isize + dx;
                let ny = cy as isize + dy;
                if nx >= 0
                    && nx < self.number_of_columns as isize
                    && ny >= 0
                    && ny < self.number_of_rows as isize
                {
                    keys.push(self.cell_coord_to_cell_key(Vector2D::new(nx as usize, ny as usize)));
                }
            }
        }

        keys
    }
    pub fn get_adjacent_cell_keys_from_position(&self, position: Vector2D<f32>) -> Vec<usize> {
        let current_cell_coord = self.particle_position_to_cell_coord(position).as_isizes();
        // use iterators beacuse faster to use
        vec![
            current_cell_coord + Vector2D::new(-1, -1),
            current_cell_coord + Vector2D::new(-1, 0),
            current_cell_coord + Vector2D::new(-1, 1),
            current_cell_coord + Vector2D::new(0, -1),
            current_cell_coord,
            current_cell_coord + Vector2D::new(0, 1),
            current_cell_coord + Vector2D::new(1, -1),
            current_cell_coord + Vector2D::new(1, 0),
            current_cell_coord + Vector2D::new(1, 1),
        ] //use into_iter when you wil not be using the values in other places skipss the need for cloning
        .into_par_iter()
        .filter_map(|coord: Vector2D<isize>| {
            //filter and map usng a Option/None where None is filterd out
            (coord.x >= 0
                && coord.x < self.number_of_columns as isize
                && coord.y >= 0
                && coord.y < self.number_of_rows as isize)
                .then_some(self.cell_coord_to_cell_key(coord.as_usizes()))
        })
        .collect()
    }

    pub fn get_particle_indexes_from_cell(&self, cell_key: usize) -> Vec<usize> {
        // needs optimizing here
        let mut particle_indexes: Vec<usize> = Vec::new();
        let mut spatial_lookup_cell: usize = cell_key;
        let mut spatial_lookup_index = self.starting_indices[cell_key];
        if spatial_lookup_index >= self.particle_count {
            return particle_indexes;
        }
        while cell_key == spatial_lookup_cell {
            let particle_index = self.spatial_lookup[spatial_lookup_index].1;
            particle_indexes.push(particle_index);
            spatial_lookup_index += 1;
            if spatial_lookup_index >= self.particle_count {
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
        (coord.x * self.number_of_rows) + coord.y
    }
    pub fn cell_key_to_coord(&self, cell_key: usize) -> (usize, usize) {
        let x = cell_key / self.number_of_rows; // Division gives the x coordinate
        let y = cell_key % self.number_of_rows; // Modulus gives the y coordinate

        (x, y)
    }
}
