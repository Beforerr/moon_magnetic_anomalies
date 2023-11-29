import numpy as np
from astropy.units import m, cm, s, A, eV, V, nT
from astropy.constants import m_e, e, c, eps0, mu0
from plasmapy.particles import Particle
from plasmapy.formulary import thermal_speed, inertial_length, MagneticDipole

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
N_r = n_solar = 5 * cm**-3  # independent variable
w_r = np.sqrt((e**2 * N_r) / (eps0 * m_e))  # reference frequency
L_r = c / w_r  # reference length
T_r = 1 / w_r
E_r = m_e * c * w_r / e
B_r = m_e * w_r / e
J_r = c * e * N_r


# Physics parameters
# ------------------

# --- plasma
kT_sw = 10 * eV  # ion and electron temperature
v_sw_x = -4e5
v_sw_y = 0
v_sw_z = 0
directed_velocity = (
    [v_sw_x, v_sw_y, v_sw_z] * m / s
)  # bulk flow velocity (m/s) following `xyz` directions
v_sw = np.linalg.norm(directed_velocity)
flux_density = n_solar * v_sw

solar_electron = Particle("e-")
solar_electron.name = "electron"
solar_proton = Particle("p+")
solar_proton.name = "proton"
solar_species = [solar_electron, solar_proton]


d_i = inertial_length(n=n_solar, particle=solar_proton)


# --- LMA field
m_d_x = 0
m_d_y = 2e13
m_d_z = 0
m_d = [m_d_x, m_d_y, m_d_z] * A * m**2  # dipolar moment
M_d = np.linalg.norm(m_d)
s_d = -0.1 * d_i  # source depth

# --- interplanetary magnetic field
b_IMF_x = 0
b_IMF_y = 3
b_IMF_z = 0
b_IMF = [b_IMF_x, b_IMF_y, b_IMF_z] * nT
B_IMF = np.linalg.norm(b_IMF)


# Numerics parameters (simulation setting)
# ----------------------------------------

# --- geometry
geometry = "2Dcartesian"
dimensions = 2


# --- grid
simulation_length = 0.5 * d_i

xmin = 0 * d_i
xmax = simulation_length * 1
ymin = -0.5 * d_i
ymax = 0.5 * d_i
zmin = -0.5 * d_i
zmax = 0.5 * d_i
L_x = xmax - xmin
L_z = zmax - zmin
L_y = ymax - ymin

space_ratio = 64
if dimensions == 2:
    grid_length = [L_x, L_z]
    number_of_cells = [_ * space_ratio for _ in [1, 2]]
if dimensions == 3:
    grid_length = [L_x, L_y, L_z]
    number_of_cells = [_ * space_ratio for _ in [2, 1, 1]]

cell_length = [0.0] * dimensions
for i in range(dimensions):
    cell_length[i] = grid_length[i] / number_of_cells[i]

# --- time
time_ratio = 1.5
simulation_time = time_ratio * (1.0 * simulation_length) / v_sw
cfl = 0.999  # default 0.999 in `WarpX`

# --- normalization
cell_length_norm = [float(i / L_r) for i in cell_length]
simulation_time_norm = float(simulation_time / T_r)
timestep_norm = cfl / np.sqrt(sum(1 / np.square(cell_length_norm)))
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

# Others
# -----------
L = ((mu0 * M_d**2) / (16 * np.pi**2 * n_solar * species.mass * v_sw**2)) ** (
    1 / 6
)
L = L.si


def magnetic_field(r_x, r_y, r_z):
    from numpy import pi

    r_0_x = s_d.si.value
    r_0_y = 0
    r_0_z = 0
    p0 = [r_0_x, r_0_y, r_0_z] * m
    _md = MagneticDipole(moment=m_d, p0=p0)
    return _md.magnetic_field([r_x, r_y, r_z])


def main():
    from icecream import ic
    from plasmapy.formulary import (
        beta,
        gyrofrequency,
        gyroradius,
        plasma_frequency,
        Debye_length,
        inertial_length,
    )

    ic(grid_length, cell_length)
    ic(timestep.si, simulation_time, timestep_norm, simulation_time_norm)
    ic(v_sw, v_sw / V_r)
    ic(B_IMF)
    ic(kT_sw, n_solar)
    ic(
        beta(T=kT_sw, n=n_solar, B=B_IMF),
        Debye_length(kT_sw, n_solar),
    )
    for species in solar_species:
        ic(
            species.name,
            species.v_th,
            species.v_rms,
            species.v_th / v_sw,
            gyrofrequency(B_IMF, species, to_hz=True),
            gyroradius(B_IMF, species, T=kT_sw),
            gyroradius(B_IMF, species, Vperp=v_sw),
            inertial_length(n=n_solar, particle=species),
            plasma_frequency(n=n_solar, particle=species, to_hz=True),
        )

    r_x = 0
    r_y = 0
    r_z = 0
    B = magnetic_field(r_x, r_y, r_z)
    B_max = B + b_IMF
    ic(B_max)
    B_max = np.linalg.norm(B_max)

    r_x = (s_d + L).si.value
    r_y = 0
    r_z = 0
    B = magnetic_field(r_x, r_y, r_z)
    B_pause = B + b_IMF
    ic(B_pause)
    B_pause = np.linalg.norm(B_pause)

    for species in solar_species:
        ic(
            species.name,
            L,
            s_d + L,
            "Gyro-radius at B_max {}: {} , {}".format(
                B_max,
                gyroradius(B_max, species, T=kT_sw),
                gyroradius(B_max, species, Vperp=v_sw),
            ),
            "Gyro-radius at B_pause {}: {}, {}".format(
                B_pause,
                gyroradius(B_pause, species, T=kT_sw),
                gyroradius(B_pause, species, Vperp=v_sw),
            ),
        )


if __name__ == "__main__":
    main()
