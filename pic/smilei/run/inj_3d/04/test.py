# PIC-CODE SMILEI
import sys

sys.path.append(".")
from parameters import *


# Numerics components
# -------------------
Main(
    geometry=geometry,
    cell_length=cell_length_norm,
    grid_length=grid_length_norm,
    number_of_patches=number_of_patches,
    patch_arrangement=patch_arrangement,
    timestep=timestep_norm,
    simulation_time=simulation_time,
    # number_of_timesteps = number_of_timesteps,    # `simulation_time` is more fundamental in c++ code implementation
    EM_boundary_conditions=EM_boundary_conditions,
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
        number_density=0,
        boundary_conditions=boundary_conditions,
        mean_velocity=mean_velocity_norm,
        temperature=[kT_sw_norm],
    )

    Species(
        name=species.name + "_test",
        position_initialization=position_initialization,
        momentum_initialization=momentum_initialization,
        particles_per_cell=1,
        mass=species.mass_norm,
        charge=species.charge_norm,
        number_density=0,
        boundary_conditions=boundary_conditions,
        mean_velocity=mean_velocity_norm,
        temperature=[kT_sw_norm],
        is_test=True,
    )

for species in Species:
    ParticleInjector(
        name="injector_" + species.name,
        species=species.name,
        box_side="xmax",
        # time_envelope = ttrapezoidal(start=0.0, plateau=Main.timestep),
        number_density=n_solar_norm,
    )

# --- LMA field

from pathlib import Path

if not Path(LMA_filename).is_file():
    import fields

    fields.main()

for field in ["Bx", "By", "Bz"]:
    ExternalField(field=field, profile=LMA_filename + "/" + field)


# --- interplanetary magnetic field
ExternalField(field="By", profile=constant(B_IMF_norm))


# Diagnostics
# -----------
_DiagFields = True
_DiagProbe = True
_DiagParticleBinning = True
_DiagTrackParticles = True
_Checkpoints = True


number_of_timesteps = int(Main.simulation_time / Main.timestep)
diagEvery = max(int(number_of_timesteps / 10), 1)
flush_every = diagEvery * 4
scalarEvery = max(number_of_timesteps / 1000, 1)
fieldEvery = diagEvery
probeEvery = max(number_of_timesteps / 1000, 1)
trackEvery = 10

dump_step = number_of_timesteps


def my_filter(particles):
    from math import prod

    _subset_ratio = 1 / prod(Main.number_of_cells) / particles_per_cell
    random_subset = np.random.random(len(particles.id)) < _subset_ratio
    if dimensions == 3:
        threshold = True
    if dimensions == 2:
        threshold = True

    return (threshold) * random_subset + (particles.id > 0)


DiagScalar(every=scalarEvery, precision=3)

if _DiagFields:
    DiagFields(every=fieldEvery, flush_every=flush_every)

if _DiagProbe:
    number_of_probes = 9
    rng = np.random.default_rng(0)
    origins = rng.random((number_of_probes, dimensions)) * Main.grid_length
    for origin in origins:
        DiagProbe(
            every=probeEvery,
            flush_every=flush_every,
            origin=origin,
        )

if _DiagParticleBinning:
    number_of_probes = 9
    rng = np.random.default_rng(0)
    origins = rng.random((number_of_probes, dimensions)) * Main.grid_length
    for index, origin in enumerate(origins):
        for species in solar_species:
            x_axes = ["x", origin[0], origin[0], 1]
            y_axes = ["y", origin[1], origin[1], 1]
            z_axes = ["z", origin[2], origin[2], 1]
            vx_axes = ["vx", -species.v_diag, species.v_diag, 100]
            vy_axes = ["vy", -species.v_diag, species.v_diag, 100]
            vz_axes = ["vz", -species.v_diag, species.v_diag, 100]
            DiagParticleBinning(
                name=species.name + str(index),
                deposited_quantity="weight",
                every=diagEvery,
                flush_every=flush_every,
                species=[species.name],
                axes=[x_axes, y_axes, z_axes, vx_axes, vy_axes, vz_axes],
            )

if _DiagTrackParticles:
    for species in Species:
        DiagTrackParticles(
            species=species.name,
            every=[0, number_of_timesteps, trackEvery],
            flush_every=100 * trackEvery,
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

if _Checkpoints:
    Checkpoints(
        dump_step=dump_step,
        exit_after_dump=True,
    )
