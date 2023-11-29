# PIC-CODE SMILEI
import sys
import unyt
from unyt import unyt_quantity

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
f_solar = trapezoidal(
    n_solar_norm,
    xvacuum=31 * grid_length[0] / 32,
    xplateau=grid_length[0],
    yplateau=1.2 * grid_length[1],
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
    name="injector_solar_electron",
    species="solar_electron",
    box_side="xmax",
)
ParticleInjector(
    name="injector_solar_proton",
    position_initialization="injector_solar_electron",
    species="solar_proton",
    box_side="xmax",
)

# --- LMA field


def B(x, y, z=L_z / 2 / L_r):
    # normalized magnetic field
    r_0 = [s_d, L_y / 2, L_z / 2]
    r_0 = unyt.unyt_array([unyt_quantity.from_astropy(_) for _ in r_0])
    # r = [x * L_r, y* L_r, z* L_r]
    r = unyt.unyt_array([unyt_quantity.from_astropy(_*L_r) for _ in [x, y, z]])
    m_d_yt = unyt.unyt_array(m_d,unyt_quantity.from_astropy(m_d.unit))
    return (
        (
            unyt_quantity.from_astropy(mu0)
            / (4 * np.pi)
            * (
                3 * (r - r_0) * unyt.udot(r - r_0, m_d_yt) / (unyt.unorm(r - r_0) ** 5)
                - m_d_yt / (unyt.unorm(r - r_0) ** 3)
            )
            / unyt_quantity.from_astropy(B_r)
        )
    ).value


def Bx(x, y):
    return B(x, y)[0]


def By(x, y):
    return B(x, y)[1]


ExternalField(field="Bx", profile=Bx)
ExternalField(field="By", profile=By)


# --- interplanetary magnetic field
ExternalField(field="By", profile=constant(B_IMF_norm))


# Diagnostics
# -----------
number_of_timesteps = int(Main.simulation_time / Main.timestep)
diagEvery = int(number_of_timesteps / 10)
probeEvery = 10
flush_every = diagEvery * 4
dump_step = int(number_of_timesteps)


def my_filter(particles):
    if dimensions == 3:
        return (
            (particles.y > grid_length[1] * 1 / 2 - cell_length[1] * 2)
            * (particles.y < grid_length[1] * 1 / 2 + cell_length[1] * 2)
            * (particles.z > grid_length[2] * 1 / 2 - cell_length[2] * 2)
            * (particles.z < grid_length[2] * 1 / 2 + cell_length[2] * 2)
        )
    if dimensions == 2:
        return (particles.y > grid_length[1] * 1 / 2 - cell_length[1] * 2) * (
            particles.y < grid_length[1] * 1 / 2 + cell_length[1] * 2
        )


DiagScalar(every=10, precision=3)

DiagFields(every=diagEvery, flush_every=flush_every)

DiagProbe(
    every=probeEvery,
    flush_every=flush_every,
    origin=[0, Main.grid_length[1] * 1 / 2],
    corners=[
        [Main.grid_length[0], Main.grid_length[1] * 1 / 2]
    ],
    number=[Main.number_of_cells[0]],
)


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
    x_axes = ["x", 0, Main.grid_length[0], Main.number_of_cells[0]]
    y_axes = [
        "y",
        Main.grid_length[1] / 2 - Main.cell_length[1] / 2,
        Main.grid_length[1] / 2 + Main.cell_length[1] / 2,
        1,
    ]
    every = probeEvery
    DiagParticleBinning(
        name=species.name + "_weight",
        deposited_quantity="weight",
        every=every,
        flush_every=flush_every,
        species=[species.name],
        axes=[
            x_axes,
            y_axes,
        ],
    )
    for weight in ["weight_vx_px", "weight_vy_py", "weight_vz_pz"]:
        DiagParticleBinning(
            name=species.name + "_" + weight,
            deposited_quantity="weight_vx_px",
            every=every,
            flush_every=flush_every,
            species=[species.name],
            axes=[
                x_axes,
                y_axes,
            ],
        )

for species in solar_species:
    DiagTrackParticles(
        species=species.name,
        every=diagEvery,
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
