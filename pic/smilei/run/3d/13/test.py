# PIC-CODE SMILEI
import sys
sys.path.append('.')
from parameters import *


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
    n_solar_norm, xvacuum=31 * grid_length[0] / 32, xplateau=grid_length[0], yplateau=1.2*grid_length[1], zplateau=1.2*grid_length[2]
)

for species in solar_species:
    Species(
        name=species.name,
        position_initialization=position_initialization,
        momentum_initialization=momentum_initialization,
        ionization_model="none",
        particles_per_cell=particles_per_cell,
        mass=species.mass_norm,
        charge=species.charge_norm,
        number_density=f_solar,
        mean_velocity=mean_velocity_norm,
        temperature=[kT_sw_norm],
        time_frozen=0.0,
        boundary_conditions=boundary_conditions,
    )

ParticleInjector(
    name= "injector_solar_electron",
    species="solar_electron",
    box_side="xmax",
)
ParticleInjector(
    name= "injector_solar_proton",
    position_initialization = 'injector_solar_electron',
    species="solar_proton",
    box_side="xmax",
)

# --- LMA field
for field in ["Bx", "By", "Bz"]:
    ExternalField(field=field, profile= LMA_filename +"/" + field)


# --- interplanetary magnetic field
ExternalField(field="By", profile=constant(B_IMF_norm))


# Diagnostics
# -----------
number_of_timesteps = int(Main.simulation_time / Main.timestep)
diagEvery = int(number_of_timesteps / 20)
flush_every = diagEvery*4
DiagScalar(every=10, precision=3)

DiagFields(every=diagEvery, flush_every=flush_every)

# PIC-CODE SMILEI
import sys
sys.path.append('.')
from parameters import *


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
    n_solar_norm, xvacuum=31 * grid_length[0] / 32, xplateau=grid_length[0], yplateau=1.2*grid_length[1], zplateau=1.2*grid_length[2]
)

for species in solar_species:
    Species(
        name=species.name,
        position_initialization=position_initialization,
        momentum_initialization=momentum_initialization,
        ionization_model="none",
        particles_per_cell=particles_per_cell,
        mass=species.mass_norm,
        charge=species.charge_norm,
        number_density=f_solar,
        mean_velocity=mean_velocity_norm,
        temperature=[kT_sw_norm],
        time_frozen=0.0,
        boundary_conditions=boundary_conditions,
    )

ParticleInjector(
    name= "injector_solar_electron",
    species="solar_electron",
    box_side="xmax",
)
ParticleInjector(
    name= "injector_solar_proton",
    position_initialization = 'injector_solar_electron',
    species="solar_proton",
    box_side="xmax",
)

# --- LMA field
for field in ["Bx", "By", "Bz"]:
    ExternalField(field=field, profile= LMA_filename +"/" + field)


# --- interplanetary magnetic field
ExternalField(field="By", profile=constant(B_IMF_norm))


# Diagnostics
# -----------
number_of_timesteps = int(Main.simulation_time / Main.timestep)
diagEvery = int(number_of_timesteps / time_ratio/10)
flush_every = diagEvery*4
dump_step=int(number_of_timesteps/time_ratio/10)
DiagScalar(every=10, precision=3)

DiagFields(every=diagEvery, flush_every=flush_every)

for species in solar_species:
    velocity_list = ["vx", "vy", "vz"]
    for velocity in velocity_list:
        DiagParticleBinning(
            deposited_quantity="weight",
            every=diagEvery,
            flush_every=flush_every,
            time_average=1,
            species=[species.name],
            axes=[
                ["x", "auto", "auto", Main.number_of_cells[0]*10],
                [
                    "y",
                    Main.grid_length[1] - Main.cell_length[1]*1/2-Main.cell_length[1]*5,
                    Main.grid_length[1] + Main.cell_length[1]*1/2+Main.cell_length[1]*5,
                    11,
                ],
                [
                    "z",
                    Main.grid_length[2] - Main.cell_length[2]*1/2-Main.cell_length[2]*5,
                    Main.grid_length[2] + Main.cell_length[2]*1/2+Main.cell_length[2]*5,
                    11,
                ],
                [velocity, -species.v_diag, species.v_diag, 100],
            ],
        )

Checkpoints(
    dump_step=dump_step,
    exit_after_dump=True, 
)


Checkpoints(
    dump_step=int(number_of_timesteps/time_ratio/10),
    exit_after_dump=True, 
)
