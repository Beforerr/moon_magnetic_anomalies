import numpy as np
from astropy.units import m, cm, s, A, eV, nT
from astropy.constants import m_e, e, c, eps0, mu0
from plasmapy.particles import Particle
from plasmapy.formulary.parameters import thermal_speed

e = e.si # TypeError: Constant 'e' does not have physically compatible units across all systems of units and cannot be combined with other values without specifying a system (eg. e.emu)

# Reference quantities
# ------------------

# --- basic reference quantities
V_r = c  # reference velocity
M_r = m_e  # reference mass
Q_r = e  # reference electric charge
K_r = m_e * c**2  # reference energy
P_r = m_e * c  # reference momentum

# --- arbitrary reference quantities
L_r = d_i = 1.3e5 * m  # reference length: ion inertial length
w_r = c / L_r  # reference frequency
T_r = 1 / w_r
E_r = m_e * c * w_r / e
B_r = m_e * w_r / e
N_r = eps0 * m_e * w_r**2 / e**2
J_r = c * e * N_r


# Physics parameters
# ------------------

# --- plasma
kT_sw = 35 * eV  # ion and electron temperature
v_sw_value = 3.5e5
v_sw = v_sw_value * m / s
mean_velocity = [-v_sw_value, 0.0, 0.0] * m / s  # bulk flow velocity (m/s)
n_solar = 3 * cm**-3
position_initialization = "random"
momentum_initialization = "maxwell-juettner"

solar_electron =  Particle('e-')
solar_electron.name = "solar_electron"
solar_proton = Particle('p+')
solar_proton.name = "solar_proton"
solar_species = [solar_electron, solar_proton]


# --- LMA field
M_d = 1.12e12
m_d = [0, M_d, 0] * A * m**2  # dipolar moment
s_d = -0.1 * d_i  # source depth
LMA_filename = 'LMA_fields_input.h5'

# --- interplanetary magnetic field
B_IMF = 3 * nT

# --- normalization
kT_sw_norm = kT_sw / K_r
n_solar_norm = n_solar / N_r
B_IMF_norm = B_IMF / B_r
mean_velocity_norm = (mean_velocity / c).value
v_abs = np.linalg.norm(mean_velocity_norm)
for species in solar_species:
    species.mass_norm = species.mass/M_r
    species.charge_norm = species.charge/Q_r


# Numerics parameters (simulation setting)
# ----------------------------------------

# --- geometry
geometry = "3Dcartesian"
dimensions = 3

# --- boundary conditions
EM_boundary_conditions = [["silver-muller"]]
boundary_conditions = [["remove"]]

# --- grid
L_x = 0.625 * d_i
L_y = 1.25 * d_i
L_z = 1.25 * d_i
grid_length = [L_x, L_y, L_z]
number_of_patches = [8, 16, 16]
cells_per_patch = [8] * dimensions
patch_arrangement = "hilbertian"
cell_length = [0.0] * dimensions
number_of_cells = [0.0] * dimensions
for i in range(dimensions):
    cell_length[i] = grid_length[i] / number_of_patches[i] / cells_per_patch[i]
    number_of_cells[i] = number_of_patches[i] * cells_per_patch[i]

particles_per_cell = 64

# --- time
time_ratio = 10
simulation_time = time_ratio * (1.0 * grid_length[0]) / v_sw

# --- normalization
grid_length = [float(i / L_r) for i in grid_length]
cell_length = [float(i / L_r) for i in cell_length]
simulation_time = float(simulation_time / T_r)
timestep = 0.95 / np.sqrt(
    1.0 / cell_length[0] ** 2 + 1.0 / cell_length[1] ** 2 + 1.0 / cell_length[2] ** 2
)

# Diagnostics
# -----------

# --- DiagParticleBinning
for species in solar_species:
    species.v_rms = thermal_speed(T=kT_sw, particle=species, method="rms", ndim=dimensions)
    species.v_rms_norm = species.v_rms/V_r
    species.v_diag = 3 * (species.v_rms_norm+v_abs)


if __name__ == '__main__':
    from icecream import ic
    for species in solar_species:
        ic(species.name, species.mass_norm, species.charge_norm, species.v_diag)