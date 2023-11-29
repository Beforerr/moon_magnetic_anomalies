from pywarpx import picmi
from parameters import *

constants = picmi.constants

_ParticleDiagnostic = False
##########################
# physics parameters
##########################

plasma_density = 3e-6
plasma_x_velocity = -3.5e5

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = number_of_timesteps
diagnostic_intervals = "::{}".format(int(number_of_timesteps / 10))

# --- Grid
d_i = 1.3e5
xmin = 0
ymin = -0.625 * d_i
zmin = -0.625 * d_i
xmax = +0.625 * d_i
ymax = +0.625 * d_i
zmax = +0.625 * d_i

n_macroparticle_per_cell = [2, 2, 2]

# --- geometry and solver
# em_solver_method = 'CKC'  # Cole-Karkkainen-Cowan stencil
# geometry = '3D'


##########################
# physics components
##########################

ratio = 7 / 8

electron_dist = picmi.UniformDistribution(
    density=plasma_density,
    lower_bound=[ratio * xmax, ymin, zmin],
    upper_bound=[xmax, ymax, zmax],
    directed_velocity=[plasma_x_velocity, 0.0, 0.0],
    rms_velocity=[solar_electron.v_rms.value, solar_electron.v_rms.value, solar_electron.v_rms.value],
    fill_in=True,
)

proton_dist = picmi.UniformDistribution(
    density=plasma_density,
    lower_bound=[ratio * xmax, ymin, zmin],
    upper_bound=[xmax, ymax, zmax],
    directed_velocity=[plasma_x_velocity, 0.0, 0.0],
    rms_velocity=[solar_proton.v_rms.value, solar_proton.v_rms.value, solar_proton.v_rms.value],
    fill_in=True,
)


electrons = picmi.Species(
    particle_type="electron", name="e-", initial_distribution=electron_dist
)
protons = picmi.Species(
    particle_type="proton", name="p+", initial_distribution=proton_dist
)

# plasma = picmi.MultiSpecies(
#                 particle_types = ['electron', 'proton'],
#                 names          = ['e-','p+'],
#                 initial_distribution=[electron_dist,proton_dist])

B_IMF = picmi.ConstantAppliedField(By=3e-9)

##########################
# numerics components
##########################
if geometry == "2Dcartesian":
    grid = picmi.Cartesian2DGrid(
        number_of_cells=number_of_cells,
        lower_bound=[xmin, ymin],
        upper_bound=[xmax, ymax],
        lower_boundary_conditions=["open"] * dimensions,
        upper_boundary_conditions=["open"] * dimensions,
        lower_boundary_conditions_particles=["absorbing"] * dimensions,
        upper_boundary_conditions_particles=["absorbing"] * dimensions,
    )
elif geometry == "3Dcartesian":
    grid = picmi.Cartesian3DGrid(
        number_of_cells=number_of_cells,
        lower_bound=[xmin, ymin, zmin],
        upper_bound=[xmax, ymax, zmax],
        lower_boundary_conditions=["open"] * dimensions,
        upper_boundary_conditions=["open"] * dimensions,
        lower_boundary_conditions_particles=["absorbing"] * dimensions,
        upper_boundary_conditions_particles=["absorbing"] * dimensions,
    )

solver = picmi.ElectromagneticSolver(grid=grid)

##########################
# diagnostics
##########################

field_diag1 = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=diagnostic_intervals,
    data_list=["rho", "rho_e-", "rho_p+", "E", "B", "J"],
    write_dir="./diags/plotfiles/",
    warpx_file_prefix="plt",
)
if _ParticleDiagnostic:
    part_diag1 = picmi.ParticleDiagnostic(
        name="diag1", period=diagnostic_intervals, species=[electrons, protons]
    )

##########################
# simulation setup
##########################

sim = picmi.Simulation(solver=solver, max_steps=max_steps)

sim.add_species(
    protons,
    layout=picmi.PseudoRandomLayout(
        grid=grid, n_macroparticles_per_cell=n_macroparticle_per_cell
    ),
)

sim.add_species(
    electrons,
    layout=picmi.PseudoRandomLayout(
        grid=grid, n_macroparticles_per_cell=n_macroparticle_per_cell
    ),
)

sim.add_applied_field(B_IMF)

sim.add_diagnostic(field_diag1)

if _ParticleDiagnostic:
    sim.add_diagnostic(part_diag1)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name="inputs{}d_from_PICMI".format(dimensions))

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()
