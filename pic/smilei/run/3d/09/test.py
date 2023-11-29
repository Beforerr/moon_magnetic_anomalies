# PIC-CODE SMILEI
from math import pi
import numpy as np
import unyt
from unyt import m, cm, s, A, c, me, qp, mu_0, eps_0, eV, nT


# Reference quantities
# ------------------

# --- basic reference quantities
V_r = c  # reference velocity
M_r = me  # reference mass
Q_r = qp  # reference electric charge
K_r = me * c**2  # reference energy
P_r = me * c  # reference momentum

# --- arbitrary reference quantities
L_r = d_i = 1.3e5 * m  # reference length: ion inertial length (m)
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
mean_velocity = [-v_sw, 0.0, 0.0] * m / s  # bulk flow velocity (m/s)
n_solar = 3 * cm**-3  # (m^-3)
particle_types = ["electron", "ion"]
species_name = ["solar_electron", "solar_ion"]
m_ion = 256 * me
position_initialization = "random"
momentum_initialization = "maxwell-juettner"

# --- LMA field
M_d = 1.12e12
m_d = [0, M_d, 0] * A * m**2  # dipolar moment
s_d = -0.1 * d_i  # source depth

# --- interplanetary magnetic field
B_IMF = 3 * nT

# --- normalization
kT_sw = float(kT_sw / K_r)
n_solar = float(n_solar / N_r)
m_ion = float(m_ion/M_r)
B_IMF = float(B_IMF / B_r)
mean_velocity = (mean_velocity / c).value


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
simulation_time = (1.0 * grid_length[0]) / v_sw

# --- normalization
grid_length = [float(i / L_r) for i in grid_length]
cell_length = [float(i / L_r) for i in cell_length]
simulation_time = float(simulation_time / T_r)
timestep = 0.95 / np.sqrt(
    1.0 / cell_length[0] ** 2 + 1.0 / cell_length[1] ** 2 + 1.0 / cell_length[2] ** 2
)


# Numerics components
# -------------------
Main(
    geometry=geometry,
    cell_length=cell_length,
    grid_length=grid_length,
    number_of_patches=number_of_patches,
    patch_arrangement=patch_arrangement,
    timestep=timestep,
    simulation_time=simulation_time,
    # number_of_timesteps = number_of_timesteps,    # `simulation_time` is more fundamental in c++ code implementation
    EM_boundary_conditions=EM_boundary_conditions,
    random_seed=smilei_mpi_rank,
)

# Physics components
# ------------------

# --- solar wind plasma
f_solar = trapezoidal(
    n_solar, xvacuum=31 * grid_length[0] / 32, xplateau=grid_length[0], yplateau=1.2*grid_length[1], zplateau=1.2*grid_length[2]
)

Species(
    name="solar_electron",
    position_initialization=position_initialization,
    momentum_initialization=momentum_initialization,
    ionization_model="none",
    particles_per_cell=particles_per_cell,
    mass=1.0,
    charge=-1.0,
    number_density=f_solar,
    mean_velocity=mean_velocity,
    temperature=[kT_sw],
    time_frozen=0.0,
    boundary_conditions=boundary_conditions,
)
ParticleInjector(
    species="solar_electron",
    box_side="xmax",
)

Species(
    name="solar_ion",
    position_initialization=position_initialization,
    momentum_initialization=momentum_initialization,
    ionization_model="none",
    particles_per_cell=particles_per_cell,
    mass=m_ion,
    charge=1.0,
    number_density=f_solar,
    mean_velocity=mean_velocity,
    temperature=[kT_sw],
    time_frozen=0.0,
    boundary_conditions=boundary_conditions,
)
ParticleInjector(
    species="solar_ion",
    box_side="xmax",
)

# --- LMA field
def B(x, y, z):
    # normalized magnetic field
    r_0 = [s_d, L_y / 2, L_z / 2]
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    arr_length = x.size
    r = np.stack((x, y, z), axis=-1) * L_r
    return (
        (
            mu_0
            / (4 * pi)
            * (
                3
                * (r - r_0)
                * unyt.udot(r - r_0, m_d).reshape(arr_length, 1)
                / (unyt.unorm(r - r_0, axis=-1) ** 5).reshape(arr_length, 1)
                - m_d / (unyt.unorm(r - r_0, axis=-1) ** 3).reshape(arr_length, 1)
            )
            / B_r
        )
    ).value


def Bx(x, y, z):
    return B(x, y, z)[::, 0]


def By(x, y, z):
    return B(x, y, z)[::, 1]


def Bz(x, y, z):
    return B(x, y, z)[::, 2]


# ExternalField(field="Bx", profile=Bx)
# ExternalField(field="By", profile=By)
# ExternalField(field="Bz", profile=Bz)

ExternalField(field="Bx", profile="Fields_input.h5/Bx")
ExternalField(field="By", profile="Fields_input.h5/By")
ExternalField(field="Bz", profile="Fields_input.h5/Bz")


# --- interplanetary magnetic field
ExternalField(field="By", profile=constant(B_IMF))


# Diagnostics
# -----------
number_of_timesteps = int(Main.simulation_time / Main.timestep)
diagEvery = int(number_of_timesteps / 20)
flush_every = diagEvery*4
DiagScalar(every=10, precision=3)

DiagFields(every=diagEvery)

for species in species_name:
    velocity_list = ["vx", "vy", "vz"]
    for velocity in velocity_list:
        DiagParticleBinning(
            deposited_quantity="weight",
            every=diagEvery,
            flush_every=flush_every,
            time_average=1,
            species=[species],
            axes=[
                ["x", "auto", "auto", 64],
                [
                    "y",
                    Main.grid_length[1] - Main.cell_length[1],
                    Main.grid_length[1] + Main.cell_length[1],
                    1,
                ],
                [
                    "z",
                    Main.grid_length[2] - Main.cell_length[2],
                    Main.grid_length[2] + Main.cell_length[2],
                    1,
                ],
                [velocity, -0.1, 0.1, 20],
            ],
        )
    DiagTrackParticles(
        species=species,
        every=diagEvery,
        flush_every=flush_every,
        attributes=["x", "y", "z", "px", "py", "pz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"],
    )

Checkpoints(
    dump_step=number_of_timesteps,
    exit_after_dump=True, 
)
