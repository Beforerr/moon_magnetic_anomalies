import numpy as np
from astropy.units import m, cm, s, A, eV, V, nT
from astropy.constants import m_e, e, c, eps0, mu0
from plasmapy.particles import Particle
from plasmapy.formulary import thermal_speed, inertial_length

e = (
    e.si
)  # TypeError: Constant 'e' does not have physically compatible units across all systems of units and cannot be combined with other values without specifying a system (eg. e.emu)

# Reference quantities
# ------------------

# --- basic reference quantities
V_r = c  # reference velocity
M_r = m_e  # reference mass
Q_r = e  # reference electric charge
K_r = m_e * c**2  # reference energy
P_r = m_e * c  # reference momentum

# --- arbitrary reference quantities
N_r = n_solar = 5 * cm**-3 # independent variable
w_r = np.sqrt((e**2*N_r) / (eps0*m_e)) # reference frequency
L_r = c / w_r  # reference length
T_r = 1 / w_r
E_r = m_e * c * w_r / e
B_r = m_e * w_r / e
J_r = c * e * N_r


# Physics parameters
# ------------------

# --- plasma
kT_sw = 10 * eV  # ion and electron temperature
v_sw_value = 4e5
v_sw = v_sw_value * m / s
mean_velocity = [-v_sw_value, 0.0, 0.0] * m / s  # bulk flow velocity (m/s)
position_initialization = "random"
momentum_initialization = "maxwell-juettner"

solar_electron = Particle("e-")
solar_electron.name = "solar_electron"
solar_proton = Particle("p+")
solar_proton.name = "solar_proton"
solar_species = [solar_electron, solar_proton]

d_i = inertial_length(n= n_solar, particle=solar_proton)


# --- LMA field
M_d_value = 1.12e13
M_d = M_d_value * A * m**2
m_d = [0, M_d_value, 0] * A * m**2  # dipolar moment
s_d = -0.01 * d_i  # source depth
LMA_filename = "LMA_fields_input.h5"

# --- interplanetary magnetic field
B_IMF = 3 * nT

# --- normalization
kT_sw_norm = kT_sw / K_r
n_solar_norm = n_solar / N_r
B_IMF_norm = B_IMF / B_r
mean_velocity_norm = (mean_velocity / c).value
v_abs = np.linalg.norm(mean_velocity_norm)
for species in solar_species:
    species.mass_norm = species.mass / M_r
    species.charge_norm = species.charge / Q_r


# Numerics parameters (simulation setting)
# ----------------------------------------

# --- geometry
geometry = "2Dcartesian"
dimensions = 2

# --- boundary conditions
EM_boundary_conditions = [["silver-muller"]]
boundary_conditions = [["remove"]]

# --- grid
L_x = 0.625 * d_i
L_y = 1.25 * d_i
L_z = 1.25 * d_i

grid_length = [L_x, L_y, L_z]
number_of_patches = [16, 16, 16]
cells_per_patch = [4] * dimensions
patch_arrangement = "hilbertian"
cell_length = [0.0] * dimensions
number_of_cells = [0.0] * dimensions
for i in range(dimensions):
    cell_length[i] = grid_length[i] / number_of_patches[i] / cells_per_patch[i]
    number_of_cells[i] = number_of_patches[i] * cells_per_patch[i]

particles_per_cell = 64

# --- time
time_ratio = 2
simulation_time = time_ratio * (1.0 * grid_length[0]) / v_sw

# --- normalization
grid_length_norm = [float(i / L_r) for i in grid_length]
cell_length_norm = [float(i / L_r) for i in cell_length]
simulation_time_norm = float(simulation_time / T_r)
timestep_norm = 0.999 / np.sqrt(sum(1 / np.square(cell_length_norm)))
number_of_timesteps = int(simulation_time_norm / timestep_norm)

timestep = timestep_norm * T_r

# Diagnostics
# -----------

# --- DiagParticleBinning
for species in solar_species:
    species.v_th = thermal_speed(T=kT_sw, particle=species, ndim=dimensions)
    species.v_rms = thermal_speed(
        T=kT_sw, particle=species, method="rms", ndim=dimensions
    )
    species.v_rms_norm = species.v_rms / V_r
    species.v_diag = 5 * (species.v_rms_norm) + v_abs

if __name__ == "__main__":
    from icecream import ic
    from plasmapy.formulary import (
        beta,
        gyrofrequency,
        gyroradius,
        plasma_frequency,
        Debye_length,
        inertial_length,
        Mag_Reynolds,
    )

    # ic(thermal_speed(T=10 * eV, particle="e-"), thermal_speed(T=10 * eV, particle="p+"))
    ic(grid_length, cell_length, grid_length_norm,  cell_length_norm)
    ic(timestep, simulation_time,timestep_norm, simulation_time_norm)
    ic(v_sw, v_sw / V_r)
    ic(B_IMF, B_IMF_norm.si)
    ic(kT_sw,n_solar)
    ic(beta(T= kT_sw, n=n_solar, B=B_IMF),Debye_length(kT_sw, n_solar),)
    for species in solar_species:
        ic(
            species.name,
            species.mass_norm,
            species.charge_norm,
            species.v_th,
            species.v_rms_norm,
            species.v_rms,
            species.v_diag,
            species.v_th / v_sw,
            gyrofrequency(B_IMF, species, to_hz=True),
            gyroradius(B_IMF, species, T=kT_sw),
            gyroradius(B_IMF, species, Vperp=v_sw),
            inertial_length(n = n_solar,particle =species),
            plasma_frequency(n = n_solar, particle = species, to_hz=True),
        )
    from numpy import pi
    r_x = 0; r_y = L_y / 2;r_z = L_z / 2
    r_0_x = s_d; r_0_y = L_y / 2;r_0_z = L_z / 2
    m_d_x, m_d_y, m_d_z = m_d
    r = np.sqrt((r_x-r_0_x)**2+(r_y- r_0_y)**2+(r_z-r_0_z)**2)
    mr_dot = (r_x - r_0_x)*m_d_x + (r_y - r_0_y)*m_d_y +  (r_z - r_0_z)*m_d_z
    Bx = mu0 / (4 * pi)*(3*(r_x-r_0_x)*(mr_dot)/r**5 - m_d_x/r**3)
    By = mu0 / (4 * pi)*(3*(r_y-r_0_y)*(mr_dot)/r**5 - m_d_y/r**3)
    Bz = mu0 / (4 * pi)*(3*(r_z-r_0_z)*(mr_dot)/r**5 - m_d_z/r**3)
    ic(Bx,By,Bz)
    B_max = By.si + B_IMF
    ic(B_max)
    for species in solar_species:
        ic("Gyro-radius at B_max {}: {}".format(B_max,gyroradius(B_max, species, T=kT_sw)))
        L = ((mu0 * M_d**2) / (16 * pi **2 * n_solar * species.mass * v_sw**2))**(1/6)
        L = L.si
        ic(species.name, L, s_d + L)