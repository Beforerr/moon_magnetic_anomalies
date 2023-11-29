# PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
# Reference: Deca, J., et.al. (2015). General mechanism and dynamics of the solar wind interaction with lunar magnetic anomalies from 3-D particle-in-cell simulations. Journal of Geophysical Research: Space Physics.
# ----------------------------------------------------------------------------------------
from math import pi
from scipy.constants import e, m_e, c, epsilon_0, mu_0, eV
import numpy as np

# Reference quantities
## Basic reference quantities
V_r = c             # reference velocity
M_r = m_e           # reference mass
Q_r = e             # reference electric charge
K_r = m_e * c**2    # reference energy
P_r = m_e * c       # reference momentum
## Arbitrary reference quantities
L_r = d_i = 1.3e5 # reference length: ion inertial length (m)
w_r = c/L_r # reference frequency
N_r = epsilon_0*m_e*w_r**2/e**2
B_r = m_e*w_r/e
T_r = 1/w_r

# %%
# Solar wind
kT_sw = 35 * eV # ion and electron temperature
v_sw = 3.5e5 # bulk flow velocity (m/s)
n_solar = 3e6 # (m^-3)
b_IMF = 3e-9
# v_th_e = sqrt(2*kT_sw/m_e)  # 9.3e5 (m/s)
# v_th_i = sqrt(2*kT_sw/m_ion)  # 6.2e4 (m/s) 

# %%
# Simulation setting
ratio = 10 # Decrease the cell precision
grid_length = [0.625*d_i, 1.25*d_i, 1.25*d_i]
patch_arrangement = "hilbertian"
number_of_patches = [8, 16, 16] # Number of patches
cells_per_patch = [8, 8, 8]
cell_length = [0.,0.,0.]
number_of_cells = [0,0,0]
for i in range(3):
    cell_length[i] = grid_length[i] / number_of_patches[i] /cells_per_patch[i]
    number_of_cells[i] = number_of_patches[i]*cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 16 # Number of patches
# Time step
simulation_time =  (1.0*grid_length[0])/v_sw

# %%
# Normalization
## ion-to-electron mass ratio
m_ion = 256*m_e/M_r
## Solar wind
kT_sw = kT_sw / K_r
v_sw = v_sw / V_r
n_solar = n_solar / N_r
## Simulation setting
grid_length = [i/L_r for i in grid_length]
cell_length = [i/L_r for i in cell_length]
simulation_time = simulation_time/T_r
timestep = 0.95/np.sqrt(1./ cell_length[0]**2 + 1./ cell_length[1]**2 + 1./ cell_length[2]**2)
diagEvery  =  int(simulation_time / timestep / 20)

# Position init
position_initialization = "random"
momentum_initialization = "maxwell-juettner"

EM_boundary_conditions = [["silver-muller"]]
boundary_conditions = [["remove"]]

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
    EM_boundary_conditions = EM_boundary_conditions,
    random_seed = smilei_mpi_rank,
)

# SPECIES BLOCKS
## solar wind
# f_solar = trapezoidal(n_solar, xvacuum=grid_length[0]-cell_length[0] ,xplateau=grid_length[0])
f_solar = trapezoidal(n_solar, xvacuum=7*grid_length[0]/8,xplateau=grid_length[0])

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
## Implementation of the LMA field
M_d = 1.12e12
m = np.array([0, M_d, 0])
r_0 = np.array([-0.1*d_i, 0, 0])

def B(x, y, z): # normalized magnetic field
    r =np.array([x, y, z])
    return mu_0/(4*pi)*(3*(np.dot(m, r-r_0))*(r-r_0)/np.linalg.norm(r-r_0)**5-m/np.linalg.norm(r-r_0)**3)/B_r

def Bx(x, y, z):
    return B(x, y, z)[0]

def By(x, y, z):
    return B(x, y, z)[1]

def Bz(x, y, z):
    return B(x, y, z)[2]

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
B_IMF =  3e-9/B_r # normalized magnetic field
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