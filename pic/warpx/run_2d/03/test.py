from pywarpx import picmi
from parameters import *


_FieldDiagnostic = True
_ParticleDiagnostic = True

##########################
# physics parameters
##########################

density = n_solar.si.value
directed_velocity = mean_velocity.si.value

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = number_of_timesteps
diagnostic_intervals = "::{}".format(int(number_of_timesteps / 20))

# --- Grid
# d_i is the ion inertial length
xmin = 0
ymin = -0.625 * d_i.si.value
zmin = -0.625 * d_i.si.value
xmax = +0.625 * d_i.si.value
ymax = +0.625 * d_i.si.value
zmax = +0.625 * d_i.si.value

n_macroparticle_per_cell = [8, 8, 8]

##########################
# physics components
##########################

ratio = 15 / 16
lower_bound=[ratio * xmax, None, None]
upper_bound=[xmax, None, None]

electron_dist = picmi.UniformDistribution(
    density=density,
    lower_bound=lower_bound,
    upper_bound=upper_bound,
    directed_velocity=directed_velocity,
    rms_velocity=[
        solar_electron.v_rms.value,
        solar_electron.v_rms.value,
        solar_electron.v_rms.value,
    ],
    fill_in=True,
)

proton_dist = picmi.UniformDistribution(
    density=density,
    lower_bound=lower_bound,
    upper_bound=upper_bound,
    directed_velocity=directed_velocity,
    rms_velocity=[
        solar_proton.v_rms.value,
        solar_proton.v_rms.value,
        solar_proton.v_rms.value,
    ],
    fill_in=True,
)


electrons = picmi.Species(
    particle_type="electron", name="electron", initial_distribution=electron_dist
)
protons = picmi.Species(
    particle_type="proton", name="proton", initial_distribution=proton_dist
)

Bx_expression = "mu0 / 4 / pi * ( 3 * (x-p0_x) * ((x - p0_x)*m_d_x + (y - p0_y)*m_d_y +  (z - p0_z)*m_d_z )/ ((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(5/2) - m_d_x/((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(3/2))"
By_expression = "mu0 / 4 / pi * ( 3 * (y-p0_y) * ((x - p0_x)*m_d_x + (y - p0_y)*m_d_y +  (z - p0_z)*m_d_z )/ ((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(5/2) - m_d_y/((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(3/2))"
Bz_expression = "B_IMF +  mu0 / 4 / pi * ( 3 * (z-p0_z) * ((x - p0_x)*m_d_x + (y - p0_y)*m_d_y +  (z - p0_z)*m_d_z )/ ((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(5/2) - m_d_z/((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(3/2))"

B_dipolar = picmi.AnalyticAppliedField(
    Bx_expression=Bx_expression,
    By_expression=By_expression,
    Bz_expression=Bz_expression,
    p0_x=s_d.si.value,
    p0_y=0.0,
    p0_z=0.0,
    m_d_x=0.0,
    m_d_y=0,
    m_d_z=M_d_value,
    B_IMF=B_IMF.si.value,
)


##########################
# numerics components
##########################
if geometry == "2Dcartesian":
    grid = picmi.Cartesian2DGrid(
        number_of_cells=number_of_cells,
        lower_bound=[xmin, zmin],
        upper_bound=[xmax, zmax],
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

solver = picmi.ElectromagneticSolver(grid=grid,cfl=cfl)

##########################
# diagnostics
##########################

field_diag1 = picmi.FieldDiagnostic(
    name='diag1',
    grid=grid,
    period=diagnostic_intervals,
    data_list=["rho", "rho_electron", "rho_proton", "E", "B", "J"],
    write_dir="./diags/plotfiles/",
    warpx_file_prefix="plt",
)
part_diag1 = picmi.ParticleDiagnostic(
    name='diag1',
    period=diagnostic_intervals,
    species=[electrons, protons],
    data_list=["position", "momentum", "weighting"],
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(solver=solver, max_steps=max_steps,particle_shape='cubic')

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

sim.add_applied_field(B_dipolar)

if _FieldDiagnostic:
    sim.add_diagnostic(field_diag1)

if _ParticleDiagnostic:
    sim.add_diagnostic(part_diag1)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file()


def main():
    sim.step()
    pass


# Alternatively, sim.step will run WarpX, controlling it from Python
if __name__ == "__main__":
    main()
