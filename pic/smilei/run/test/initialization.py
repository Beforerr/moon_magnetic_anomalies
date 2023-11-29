# PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
# Reference: Deca, J., et.al. (2015). General mechanism and dynamics of the solar wind interaction with lunar magnetic anomalies from 3-D particle-in-cell simulations. Journal of Geophysical Research: Space Physics.
# ----------------------------------------------------------------------------------------
from math import pi
import numpy as np
import unyt
from unyt import m, cm, s, A, c, me, qp, mu_0, eps_0, eV, nT
# ----------------------------------------------------------------------------------------
# A: ampere
# me: mass_electron, electron_mass;
# qp: elementary_charge
# ----------------------------------------------------------------------------------------

# Reference quantities
## Basic reference quantities
V_r = c # reference velocity
M_r = me    # reference mass
Q_r = qp   # reference electric charge
## Arbitrary reference quantities
L_r = d_i = 1.3e5 * m # reference length: ion inertial length (m)
w_r = c/L_r # reference frequency
B_r = me*w_r/qp

# %%
# Simulation setting
L_x = 0.625*d_i; L_y = 1.25*d_i; L_z = 1.25*d_i;
grid_length = [L_x, L_y, L_z]
number_of_patches = [8, 16, 16] # Number of patches
cells_per_patch = [8, 8, 8]
cell_length = [0.,0.,0.]
number_of_cells = [0,0,0]
for i in range(3):
    cell_length[i] = grid_length[i] / number_of_patches[i] /cells_per_patch[i]
    number_of_cells[i] = number_of_patches[i]*cells_per_patch[i]

# --------------------------------------
# Implementation of the LMA field
# --------------------------------------
M_d = 1.12e12 
m_d = [0, M_d, 0] * A * m**2 # dipolar moment
r_0 = [-0.1*d_i, L_y/2, L_z/2]

# normalized magnetic field
def B(x, y, z):
    x = np.array(x);y = np.array(y);z = np.array(z);
    arr_length = x.size
    r = np.stack((x, y, z), axis=-1) * L_r
    return (mu_0/(4*pi)*( 3*(r-r_0)*unyt.udot(r-r_0, m_d).reshape(arr_length, 1)/(unyt.unorm(r-r_0, axis=-1)**5).reshape(arr_length, 1)  - m_d / (unyt.unorm(r-r_0, axis=-1)**3).reshape(arr_length, 1))/B_r).in_units('dimensionless').value

def Bx(x, y, z):
    return B(x, y, z)[::,0]

def By(x, y, z):
    return B(x, y, z)[::,1]

def Bz(x, y, z):
    return B(x, y, z)[::,2]

