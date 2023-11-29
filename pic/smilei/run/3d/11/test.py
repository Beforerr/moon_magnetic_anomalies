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
    name= "injector_solar_electron",
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
    name= "injector_solar_ion",
    position_initialization = 'injector_solar_electron',
    species="solar_ion",
    box_side="xmax",
)

# --- LMA field
filename = 'LMA_fields_input.h5'

for field in ["Bx", "By", "Bz"]:
    ExternalField(field=field, profile= filename +"/" + field)


# --- interplanetary magnetic field
ExternalField(field="By", profile=constant(B_IMF))


# Diagnostics
# -----------
number_of_timesteps = int(Main.simulation_time / Main.timestep)
diagEvery = int(number_of_timesteps / 20)
flush_every = diagEvery*4
DiagScalar(every=10, precision=3)

DiagFields(every=diagEvery, flush_every=flush_every)

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
                [velocity, -2*v_abs, 2*v_abs, 200],
            ],
        )
    DiagTrackParticles(
        species=species,
        every=diagEvery,
        flush_every=flush_every,
        attributes=["x", "y", "z", "px", "py", "pz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"],
    )

Checkpoints(
    dump_step=number_of_timesteps/time_ratio,
    exit_after_dump=True, 
)
