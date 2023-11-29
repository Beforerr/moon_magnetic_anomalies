from pywarpx import picmi
from parameters import *
import numpy as np

_FieldDiagnostic = True
_ParticleDiagnostic = True

##########################
# physics parameters
##########################

density = n_solar.si.value
directed_velocity = directed_velocity.si.value

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = number_of_timesteps

# --- Grid
# d_i is the ion inertial length
xmin = xmin.si.value
ymin = ymin.si.value
zmin = zmin.si.value
xmax = xmax.si.value
ymax = ymax.si.value
zmax = zmax.si.value

n_macroparticle_per_cell = [1, 1, 1]

##########################
# physics components
##########################

ratio = 15/16
lower_bound = [ratio * simulation_length.si.value, None, None]
upper_bound = [None, None, None]

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

Bx_expression = "b_IMF_x + mu0 / 4 / pi * ( 3 * (x-p0_x) * ((x - p0_x)*m_d_x + (y - p0_y)*m_d_y +  (z - p0_z)*m_d_z )/ ((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(5/2) - m_d_x/((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(3/2))"
By_expression = "b_IMF_y + mu0 / 4 / pi * ( 3 * (y-p0_y) * ((x - p0_x)*m_d_x + (y - p0_y)*m_d_y +  (z - p0_z)*m_d_z )/ ((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(5/2) - m_d_y/((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(3/2))"
Bz_expression = "b_IMF_z + mu0 / 4 / pi * ( 3 * (z-p0_z) * ((x - p0_x)*m_d_x + (y - p0_y)*m_d_y +  (z - p0_z)*m_d_z )/ ((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(5/2) - m_d_z/((x-p0_x)**2+(y- p0_y)**2+(z-p0_z)**2)**(3/2))"

B_dipolar = picmi.AnalyticAppliedField(
    Bx_expression=Bx_expression,
    By_expression=By_expression,
    Bz_expression=Bz_expression,
    p0_x=s_d.si.value,
    p0_y=0,
    p0_z=0,
    m_d_x=m_d[0].si.value,
    m_d_y=m_d[1].si.value,
    m_d_z=m_d[2].si.value,
    b_IMF_x=b_IMF[0].si.value,
    b_IMF_y=b_IMF[1].si.value,
    b_IMF_z=b_IMF[2].si.value,
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
        lower_bound=[xmin, zmin, ymin],
        upper_bound=[xmax, zmax, ymax],
        lower_boundary_conditions=["open"] * dimensions,
        upper_boundary_conditions=["open"] * dimensions,
        lower_boundary_conditions_particles=["absorbing"] * dimensions,
        upper_boundary_conditions_particles=["absorbing"] * dimensions,
    )

solver = picmi.ElectromagneticSolver(grid=grid, cfl=cfl)

##########################
# diagnostics
##########################

diagnostic_intervals = "::{}".format(int(max_steps / 10))
field_diag1 = picmi.FieldDiagnostic(
    grid=grid,
    period=diagnostic_intervals,
    data_list=["rho", "rho_electron", "rho_proton", "E", "B", "J"],
    write_dir="./diags/diag1/",
    upper_bound=[simulation_length.si.value,zmax],
    warpx_file_prefix="plt",
)
part_diag1 = picmi.ParticleDiagnostic(
    period=diagnostic_intervals,
    species=[electrons, protons],
    data_list=["position", "momentum", "weighting"],
    write_dir="./diags/diag2/",
    warpx_plot_filter_function =  "(x<warpx_simulation_length)", 
    warpx_simulation_length = simulation_length.si.value,
    warpx_file_prefix="plt",
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    particle_shape="cubic",
    warpx_load_balance_intervals=128,
)

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


def append_flux():
    with open("inputs", mode="a", encoding="utf-8") as f:
        for species in solar_species:
            f.write("{}.injection_style = nfluxpercell\n".format(species.name))
            f.write("{}.surface_flux_pos = {}\n".format(species.name, xmax))
            f.write("{}.flux_normal_axis = x\n".format(species.name))
            f.write("{}.flux_direction = -1\n".format(species.name))
            f.write("{}.density = {}\n".format(species.name, flux_density.si.value))


reduced_diags = """
warpx.reduced_diags_names = probe_fp
probe_fp.type = FieldProbe
probe_fp.probe_geometry = Point
probe_fp.x_probe = {}
probe_fp.y_probe = 0
probe_fp.z_probe = 0
probe_fp.intervals = ::10
probe_fp.raw_fields = true""".format(
    (s_d + L).si.value
)


def append_reduced_diags():
    with open("inputs", mode="a", encoding="utf-8") as f:
        f.write(reduced_diags)


def main():
    append_reduced_diags()
    # sim.step()
    pass


# Alternatively, sim.step will run WarpX, controlling it from Python
if __name__ == "__main__":
    main()
