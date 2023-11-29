# ----------------------------------------------------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Particle injection from the Xmin and Xmax boundaries
#
# ----------------------------------------------------------------------------------------

import math
import numpy as np

# Mean velocity
mean_velocity = 0.3
# Electron temperature
Te = 1
# Ion temperature
Ti = 1
# Ion charge
Zi = 1
# Density
n0 = 1
# Debye length
Debye_length = 1. / np.sqrt( n0 / Te + Zi * n0 / Ti )
# Cell length
cell_length = [Debye_length*0.5, Debye_length*0.5]
# Number of patches
number_of_patches =[4, 16]
# Cells per patches (patch shape)
cells_per_patch = [32., 8.]
# Grid length
grid_length = [0.,0.]
for i in range(2):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 32
# Position init
position_initialization = 'random'
# Time step
timestep = 0.95/np.sqrt(1./ cell_length[0]**2 + 1./ cell_length[1]**2 )
# Total simulation time
simulation_time = 3*((1.5 - 0.125)*grid_length[0])/mean_velocity          # duration of the simulation
# Period of output for the diags
diag_every = int(simulation_time / timestep)
# Boundary conditions for particles
particle_boundary_conditions = [["remove", "remove"],["remove", "remove"]]
# Boundary conditions for fields
field_boundary_conditions = [['silver-muller'],['silver-muller']]

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2 ,
    cell_length = cell_length,
    grid_length  = grid_length,
    number_of_patches = number_of_patches,
    #cell_sorting = True,
    timestep = timestep,
    simulation_time = simulation_time,
    EM_boundary_conditions = field_boundary_conditions,
    random_seed = smilei_mpi_rank,
)

LoadBalancing(
	every = 100
)

# Initial plasma shape
ratio = 32.
fp = trapezoidal(1., xvacuum=0.                 ,xplateau=grid_length[0]/ratio)

Species(
	name = 'pon1',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1836.0,
	charge = 1.0,
	number_density = fp,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [Ti],
	time_frozen = 0.0,
	boundary_conditions = particle_boundary_conditions,
)

ParticleInjector(
    species = 'pon1',
    box_side = 'xmin',
)

Species(
	name = 'eon1',
	position_initialization = position_initialization,
	momentum_initialization = 'mj',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	c_part_max = 1.0,
	mass = 1.0,
	charge = -1.0,
	number_density = fp,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [Te],
	time_frozen = 0.0,
	boundary_conditions = particle_boundary_conditions,
)
ParticleInjector(
    species = 'eon1',
    box_side = 'xmin',
)

# fd = trapezoidal(1., yvacuum=0.					,yplateau=grid_length[0]/ratio)
# fu = trapezoidal(1., yvacuum=(ratio-1)*grid_length[0]/ratio,yplateau=grid_length[0])

# Species(
# 	name = 'pon2',
# 	position_initialization = position_initialization,
# 	momentum_initialization = 'mj',
# 	ionization_model = 'none',
# 	particles_per_cell = particles_per_cell,
# 	c_part_max = 1.0,
# 	mass = 1836.0,
# 	charge = 1.0,
# 	number_density = fd,
# 	mean_velocity = [0.,mean_velocity,0.],
# 	temperature = [Ti],
# 	time_frozen = 0.0,
# 	boundary_conditions = particle_boundary_conditions,
# )

# ParticleInjector(
#     species = 'pon2',
#     box_side = 'ymin',
# )

# Species(
# 	name = 'eon2',
# 	position_initialization = position_initialization,
# 	momentum_initialization = 'mj',
# 	ionization_model = 'none',
# 	particles_per_cell = particles_per_cell,
# 	c_part_max = 1.0,
# 	mass = 1.0,
# 	charge = -1.0,
# 	number_density = fd,
# 	mean_velocity = [0.,mean_velocity,0.],
# 	temperature = [Te],
# 	time_frozen = 0.0,
# 	boundary_conditions = particle_boundary_conditions,
# )
# ParticleInjector(
#     species = 'eon2',
#     box_side = 'ymin',
# )


DiagScalar(every=1)

DiagFields(
    every = 10, 
    fields = ['Rho_pon1','Rho_eon1']
)
