# PIC-CODE SMILEI
import sys

sys.path.append(".")
from utilities.parameters import *


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
for species in solar_species:
    Species(
        name=species.name,
        position_initialization=position_initialization,
        momentum_initialization=momentum_initialization,
        particles_per_cell=particles_per_cell,
        mass=species.mass_norm,
        charge=species.charge_norm,
        number_density=n_solar_norm,
        boundary_conditions=boundary_conditions,
        mean_velocity=mean_velocity_norm,
        temperature=[kT_sw_norm],
    )

# inj_ratio = 1 / 8

# f_solar = trapezoidal(
#     n_solar_norm,
#     xvacuum=grid_length[0] - grid_length[0] * inj_ratio,
#     xplateau=grid_length[0],
# )

ParticleInjector(
    name="injector_solar_electron",
    species="solar_electron",
    box_side="xmax",
    number_density=n_solar_norm,
)
ParticleInjector(
    name="injector_solar_proton",
    species="solar_proton",
    box_side="xmax",
    position_initialization="injector_solar_electron",
    number_density=n_solar_norm,
)

# f_solar_y_max = trapezoidal(
#     n_solar_norm,
#     yvacuum=grid_length[1] - grid_length[1] * inj_ratio,
#     yplateau=grid_length[1],
# )
# ParticleInjector(
#     name="injector_solar_electron_y_max",
#     species="solar_electron",
#     box_side="ymax",
#     number_density=n_solar_norm,
# )
# ParticleInjector(
#     name="injector_solar_proton_y_max",
#     species="solar_proton",
#     box_side="ymax",
#     position_initialization="injector_solar_electron_y_max",
#     number_density=n_solar_norm,
# )

# f_solar_y_min = trapezoidal(
#     n_solar_norm,
#     yplateau=grid_length[1] * inj_ratio,
# )
# ParticleInjector(
#     name="injector_solar_electron_y_min",
#     species="solar_electron",
#     box_side="ymin",
#     number_density=n_solar_norm,
# )
# ParticleInjector(
#     name="injector_solar_proton_y_min",
#     species="solar_proton",
#     box_side="ymin",
#     position_initialization="injector_solar_electron_y_min",
#     number_density=n_solar_norm,
# )


# --- interplanetary magnetic field
ExternalField(field="By", profile=constant(B_IMF_norm))


# Diagnostics
# -----------
number_of_timesteps = int(Main.simulation_time / Main.timestep)
diagEvery = max(int(number_of_timesteps / 10), 1)
flush_every = diagEvery * 4
scalarEvery = max(number_of_timesteps / 1000, 1)
fieldEvery = diagEvery
probeEvery = max(number_of_timesteps / 1000, 1)

dump_step = number_of_timesteps


def my_filter(particles):
    _subset_ratio = 1e-9
    random_subset = np.random.rand(len(particles.id)) < _subset_ratio
    if dimensions == 3:
        # threshold = (
        #     (particles.y > grid_length[1] * 1 / 2 - cell_length[1] * 2)
        #     * (particles.y < grid_length[1] * 1 / 2 + cell_length[1] * 2)
        #     * (particles.z > grid_length[2] * 1 / 2 - cell_length[2] * 2)
        #     * (particles.z < grid_length[2] * 1 / 2 + cell_length[2] * 2)
        # )
        threshold = True
    if dimensions == 2:
        # threshold = (particles.y > grid_length[1] * 1 / 2 - cell_length[1] * 2) * (
        #     particles.y < grid_length[1] * 1 / 2 + cell_length[1] * 2
        # )
        threshold = True

    return (threshold) * random_subset + (particles.id > 0)


DiagScalar(every=scalarEvery, precision=3)

DiagFields(every=fieldEvery, flush_every=flush_every)


for species in solar_species:
    x_axes = ["x", 0, Main.grid_length[0], Main.number_of_cells[0]]

    y_axes = [
        "y",
        Main.grid_length[1] / 2 - Main.cell_length[1] * 11 / 2,
        Main.grid_length[1] / 2 + Main.cell_length[1] * 11 / 2,
        11,
    ]

    velocity_list = ["vx", "vy", "vz"]
    for velocity in velocity_list:
        DiagParticleBinning(
            name=species.name + "_" + velocity,
            deposited_quantity="weight",
            every=diagEvery,
            flush_every=flush_every,
            time_average=1,
            species=[species.name],
            axes=[
                x_axes,
                y_axes,
                [velocity, -species.v_diag, species.v_diag, 100],
            ],
        )


for species in solar_species:
    DiagTrackParticles(
        species=species.name,
        every=[0, number_of_timesteps, 100],
        flush_every=flush_every,
        filter=my_filter,
        attributes=[
            "x",
            "y",
            "px",
            "py",
            "pz",
            "Ex",
            "Ey",
            "Ez",
            "Bx",
            "By",
            "Bz",
            "w",
        ],
    )

Checkpoints(
    dump_step=dump_step,
    exit_after_dump=True,
)
