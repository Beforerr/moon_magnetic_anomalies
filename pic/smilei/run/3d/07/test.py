# PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
# Reference: Deca, J., et.al. (2015). General mechanism and dynamics of the solar wind interaction with lunar magnetic anomalies from 3-D particle-in-cell simulations. Journal of Geophysical Research: Space Physics.
# ----------------------------------------------------------------------------------------
from math import pi
import numpy as np
import unyt
from unyt import m, cm, s, A, c, me, qp, mu_0, eps_0, eV, nT
# ----------------------------------------------------------------------------------------
# A: ampere
# me: mass_electron, electron_mass;
# qp: elementary_charge
# ----------------------------------------------------------------------------------------

# Reference quantities
## Basic reference quantities
V_r = c # reference velocity
M_r = me    # reference mass
Q_r = qp   # reference electric charge
K_r = me * c**2    # reference energy
P_r = me * c       # reference momentum
## Arbitrary reference quantities
L_r = d_i = 1.3e5 * m # reference length: ion inertial length (m)
w_r = c/L_r # reference frequency
T_r = 1/w_r
E_r = me*c*w_r/qp
B_r = me*w_r/qp
N_r = eps_0*me*w_r**2/qp**2
J_r = c*qp*N_r

# %%
# Solar wind
kT_sw = 35 * eV # ion and electron temperature
v_sw = 3.5e5 * m/s # bulk flow velocity (m/s)
n_solar = 3 * cm**-3 # (m^-3)
B_IMF = 3 * nT

# %%
# Simulation setting
L_x = 0.625*d_i; L_y = 1.25*d_i; L_z = 1.25*d_i;
grid_length = [L_x, L_y, L_z]
patch_arrangement = "hilbertian"
number_of_patches = [8, 16, 16] # Number of patches
cells_per_patch = [8, 8, 8]
cell_length = [0.,0.,0.]
number_of_cells = [0,0,0]
for i in range(3):
    cell_length[i] = grid_length[i] / number_of_patches[i] /cells_per_patch[i]
    number_of_cells[i] = number_of_patches[i]*cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 64
# Time step
simulation_time =  (1.0*grid_length[0])/v_sw

# --------------------------------------
# Implementation of the LMA field
# --------------------------------------
M_d = 1.12e12 
m_d = [0, M_d, 0] * A * m**2 # dipolar moment
r_0 = [-0.1*d_i, L_y/2, L_z/2]

# normalized magnetic field
def B(x, y, z):
    x = np.array(x);y = np.array(y);z = np.array(z);
    arr_length = x.size
    r = np.stack((x, y, z), axis=-1) * L_r
    return ((mu_0/(4*pi)*( 3*(r-r_0)*unyt.udot(r-r_0, m_d).reshape(arr_length, 1)/(unyt.unorm(r-r_0, axis=-1)**5).reshape(arr_length, 1)  - m_d / (unyt.unorm(r-r_0, axis=-1)**3).reshape(arr_length, 1))/B_r)).value
    # `in_units('dimensionless'))`
    # `convert_to_units('dimensionless')`

def Bx(x, y, z):
    return B(x, y, z)[::,0]

def By(x, y, z):
    return B(x, y, z)[::,1]

def Bz(x, y, z):
    return B(x, y, z)[::,2]

# --------------------------------------
# Normalization and Stripping units off of data
# --------------------------------------
m_ion   =   float(256*me/M_r) # ion-to-electron mass ratio

## Solar wind
kT_sw   =   float(kT_sw / K_r)
v_sw    =   float(v_sw / V_r)
n_solar =   float(n_solar / N_r)
B_IMF   =   float(B_IMF / B_r) # normalized magnetic field

## Simulation setting
grid_length = [float(i/L_r) for i in grid_length]
cell_length = [float(i/L_r) for i in cell_length]
simulation_time = float(simulation_time/T_r)
timestep = 0.95/np.sqrt(1./ cell_length[0]**2 + 1./ cell_length[1]**2 + 1./ cell_length[2]**2)
number_of_timesteps = int(simulation_time / timestep)
diagEvery  =  int(number_of_timesteps / 20)

# Position init
position_initialization = "random"
momentum_initialization = "maxwell-juettner"

EM_boundary_conditions = [["silver-muller"]]
# boundary_conditions = [["remove"]]
boundary_conditions = [["stop", "remove"], ["remove","remove"], ["remove","remove"]]
boundary_conditions = [["remove", "remove"], ["remove","remove"], ["remove","remove"]]

mean_velocity = [-v_sw, 0., 0.]

# --------------------------------------
# SMILEI's VARIABLES (defined in blocks)
# --------------------------------------
# MAIN BLOCK
Main(
    geometry = "3Dcartesian",
    cell_length = cell_length,
    grid_length = grid_length,
    number_of_patches = number_of_patches,
    patch_arrangement = patch_arrangement,
    timestep = timestep,
    simulation_time = simulation_time,
    # number_of_timesteps = number_of_timesteps,    # `simulation_time` is more fundamental in c++ code implementation 
    EM_boundary_conditions = EM_boundary_conditions,
    random_seed = smilei_mpi_rank,
)

# SPECIES BLOCKS
## solar wind
# f_solar = trapezoidal(n_solar, xvacuum=grid_length[0]-cell_length[0] ,xplateau=grid_length[0])
f_solar = trapezoidal(n_solar, xvacuum=31*grid_length[0]/32,xplateau=grid_length[0])

Species(
    name = 'solar_ion',
    position_initialization = position_initialization,
    momentum_initialization = momentum_initialization,
    ionization_model = 'none',
    particles_per_cell = particles_per_cell,
    mass = m_ion,
    charge = 1.,
    number_density = f_solar,
    mean_velocity = mean_velocity,
    temperature = [kT_sw],
    time_frozen = 0.0,
    boundary_conditions = boundary_conditions, 
)
ParticleInjector(
    species = 'solar_ion',
    box_side = 'xmax',
)

Species(
    name = 'solar_electron',
    position_initialization = position_initialization,
    momentum_initialization = momentum_initialization,
    ionization_model = 'none',
    particles_per_cell = particles_per_cell,
    mass = 1.,
    charge = -1.,
    number_density = f_solar,
    mean_velocity = mean_velocity,
    temperature = [kT_sw],
    time_frozen = 0.0,
    boundary_conditions = boundary_conditions,
)
ParticleInjector(
    species = 'solar_electron',
    box_side = 'xmax',
)

# EXTERNAL FIELDS BLOCKS
ExternalField(
    field = "Bx",
    profile = Bx
)
ExternalField(
    field = "By",
    profile = By
)
ExternalField(
    field = "Bz",
    profile = Bz
)

## Implementation of the interplanetary magnetic field
ExternalField(
    field = "By",
    profile = constant(B_IMF)
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------
DiagScalar(every=1)

DiagFields(
    every = diagEvery
)

Checkpoints(
    dump_step = number_of_timesteps,
    exit_after_dump = True,
)