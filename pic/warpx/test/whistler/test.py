from pywarpx import picmi
import numpy as np

from astropy.units import m, cm, s, A, eV, V, nT, Hz
from astropy.constants import m_e, e, c, eps0, mu0

from plasmapy.formulary import thermal_speed

e = e.si


geometry = "2Dcartesian"
dimensions = 2

_FieldDiagnostic = True
_ParticleDiagnostic = True

##########################
# physical parameters
##########################

omega_r = omega_pe = 1 * Hz # deduced from B_r
density = n_e = n_p = eps0 * m_e * omega_pe**2 / e**2

Omega_e = omega_pe / 8.0  # cyclotron frequency
beta_par_e = 1.0
B_0 = m_e*Omega_e/e

T_par_e = beta_par_e * B_0**2 / (2*mu0*n_e)
T_ver_e = 3.0 * T_par_e

v_par_e = thermal_speed(T=T_par_e, particle='electron', ndim=dimensions)
v_ver_e = thermal_speed(T=T_ver_e, particle='electron', ndim=dimensions)

##########################
# numerics parameters
##########################

# --- Number of time steps
dt = 0.125/omega_pe
max_steps = int(400/omega_pe/dt)

# --- Grid
# d_i is the ion inertial length

dx = 0.25*c/omega_pe
dy = 0.25*c/omega_pe
dz = 0.25*c/omega_pe
nx = 96
ny = 8
nz = 1

xmin = 0 * m
ymin = 0 * m
zmin = 0 * m
xmax = nx * dx
ymax = ny * dy
zmax = nz * dz

n_macroparticle_per_cell = [1, 1, 1]

##########################
# physics components
##########################
electron_dist = picmi.UniformDistribution(
    density=density.si.value,
    rms_velocity=[
        v_par_e.si.value,
        v_ver_e.si.value,
        v_ver_e.si.value,
    ],
)

proton_dist = picmi.UniformDistribution(
    density=density.si.value,
)


electrons = picmi.Species(
    particle_type="electron", name="electron", initial_distribution=electron_dist
)
protons = picmi.Species(
    particle_type="proton", name="proton", initial_distribution=proton_dist
)

field = picmi.ConstantAppliedField(Bx=B_0.si.value)

##########################
# numerics components
##########################

if geometry == "2Dcartesian":
    number_of_cells = [nx, ny]
    grid = picmi.Cartesian2DGrid(
        number_of_cells = number_of_cells,
        lower_bound=[xmin.si.value, ymin.si.value],
        upper_bound=[xmax.si.value, ymax.si.value],
        lower_boundary_conditions=['periodic'] * dimensions,
        upper_boundary_conditions=['periodic'] * dimensions,
        lower_boundary_conditions_particles=['periodic'] * dimensions,
        upper_boundary_conditions_particles=['periodic'] * dimensions,
    )


elif geometry == "3Dcartesian":
    number_of_cells = [nx, ny, nz]
    grid = picmi.Cartesian3DGrid(
        number_of_cells=number_of_cells,
        lower_bound=[xmin.si.value, ymin.si.value, zmin.si.value],
        upper_bound=[xmax.si.value, ymax.si.value, zmax.si.value],
        lower_boundary_conditions=['periodic'] * dimensions,
        upper_boundary_conditions=['periodic'] * dimensions,
        lower_boundary_conditions_particles=['periodic'] * dimensions,
        upper_boundary_conditions_particles=['periodic'] * dimensions,
    )

solver = picmi.ElectromagneticSolver(grid=grid)

##########################
# diagnostics
##########################
diagnostic_intervals = "::{}".format(int(max_steps / 10))

field_diag = picmi.FieldDiagnostic(
    grid=grid,
    period="::{}".format(int(max_steps / 50)),
    data_list=["rho", "rho_electron", "rho_proton", "E", "B", "J"],
    warpx_format="openpmd",
    warpx_openpmd_backend = "h5",
)

part_diag = picmi.ParticleDiagnostic(
    period="::{}".format(int(max_steps / 10)),
    species=[electrons, protons],
    data_list=["position", "momentum", "weighting"],
    warpx_format="openpmd",
    warpx_openpmd_backend = "h5",
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    # particle_shape="cubic",
    warpx_load_balance_intervals="::32",
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

sim.add_applied_field(field)

if _FieldDiagnostic:
    sim.add_diagnostic(field_diag)

if _ParticleDiagnostic:
    sim.add_diagnostic(part_diag)

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
