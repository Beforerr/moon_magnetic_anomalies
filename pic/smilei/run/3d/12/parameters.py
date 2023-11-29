import numpy as np
from unyt import m, cm, s, A, c, me, qp, eps_0, mu_0, eV, nT
from unyt import unyt_quantity
from plasmapy.particles import Particle, proton, electron
from plasmapy.formulary.parameters import thermal_speed

# Reference quantities
# ------------------

# --- basic reference quantities
V_r = c  # reference velocity
M_r = me  # reference mass
Q_r = qp  # reference electric charge
K_r = me * c**2  # reference energy
P_r = me * c  # reference momentum

# --- arbitrary reference quantities
L_r = d_i = 1.3e5 * m  # reference length: ion inertial length
w_r = c / L_r  # reference frequency
T_r = 1 / w_r
E_r = me * c * w_r / qp
B_r = me * w_r / qp
N_r = eps_0 * me * w_r**2 / qp**2
J_r = c * qp * N_r


# Physics parameters
# ------------------

# --- plasma
kT_sw = 35 * eV  # ion and electron temperature
v_sw_value = 3.5e5
v_sw = v_sw_value * m / s
mean_velocity = [-v_sw_value, 0.0, 0.0] * m / s  # bulk flow velocity (m/s)
n_solar = 3 * cm**-3
particle_types = ["electron", "ion"]
species_name = ["solar_electron", "solar_ion"]
m_ion = 256 * me
position_initialization = "random"
momentum_initialization = "maxwell-juettner"

solar_electron =  Particle('e-')
solar_proton = Particle('p+')
solar_species = [solar_electron, solar_proton]

# --- LMA field
M_d = 1.12e12
m_d = [0, M_d, 0] * A * m**2  # dipolar moment
s_d = -0.1 * d_i  # source depth

# --- interplanetary magnetic field
B_IMF = 3 * nT

# --- normalization
kT_sw_norm = float(kT_sw / K_r)
n_solar = float(n_solar / N_r)
m_ion = float(m_ion/M_r)
B_IMF = float(B_IMF / B_r)
mean_velocity_norm = (mean_velocity / c).value
v_abs = np.linalg.norm(mean_velocity_norm)


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
v_electron_rms = thermal_speed(T=kT_sw.to_astropy(), particle=electron, method="rms", ndim=dimensions)
v_proton_rms = thermal_speed(T=kT_sw.to_astropy(), particle=proton, method="rms", ndim=dimensions)
v_electron_rms_norm = float(unyt_quantity.from_astropy(v_electron_rms)/V_r)
v_proton_rms_norm = float(unyt_quantity.from_astropy(v_proton_rms)/V_r)
v_electron_diag = 3 * (v_electron_rms_norm+v_abs)
v_proton_diag = 3 * (v_proton_rms_norm+v_abs)
v_diags = [v_electron_diag, v_proton_diag]

if __name__ == '__main__':
    from icecream import ic
    solar_electron =  Particle('e-')
    solar_proton = Particle('p+')
    solar_species = [solar_electron, solar_proton]
    for species in solar_species:
        species.v_rms = thermal_speed(T=kT_sw.to_astropy(), particle=species, method="rms", ndim=dimensions)
        species.v_rms_norm = float(unyt_quantity.from_astropy(species.v_rms)/V_r)
        species.v_diag = 3 * (species.v_rms_norm+v_abs)
        ic(species.v_diag)
