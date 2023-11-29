# -------------------
# MY PYTHON VARIABLES
# ------------------
import numpy as np
import math
import matplotlib.pyplot as plt
pi = math.pi
# --------------------------------------
# Define the function
# --------------------------------------
def deltaRing(N,v):
    theta =np.linspace(0, 2*pi, N)
    vx = v*np.cos(theta)
    vy = v*np.sin(theta)
    vz = v*np.zeros(len(theta))
    p_init = np.array([vx, vy, vz])
    return p_init

def position_init(N, nppc, dx):
    s = 0
    x = []
    ddx = dx/nppc
    for i in range(int(N)):
        x_dx = np.linspace(s+dx*i+ddx, s+dx*(i+1)-ddx, nppc)
        x = np.append(x, x_dx)
    return x

#basic parameters_
n = 1.          #plasma number density
#wr = wp,
wp = 1.       #proton plasma frequency
c = 1.
me = 1.
e = 1.
mp = 1600.*me

vA = c/20
wcp = vA*wp/c    #proton gyrofrequency

#Spatial resolution
lamda_i = c/wp
dx = 7.5*10e-4*lamda_i
#Time resolution
dt = 2.5e-5/wcp
tsim = dt*4e4
B0 = wcp*mp/e

#grid setting
Nx = 8192
Lsim = Nx*dx
nppc = 10
number_cell = nppc*Nx

#ring proton parameter
vr = 1.0*vA

#set temperture
vcth = 0.00635*vA             #thermal velocity of cold protons
Tp = 0.5*(vcth**2)*mp         #cold proton temperture
Te = Tp                       #cold electron temperture

# --------------------------------------
# Initialize the momentum and position of the ion
# --------------------------------------
mo_init = deltaRing(number_cell, vA)
pos_init = position_init(Nx, nppc, dx)
w_ion = np.full(int(number_cell), n/nppc * dx)  #w_ion = n/(nppc) * dx
po_init_w = np.vstack((pos_init, w_ion))

diagEvery  = 200  # frequency of outputs for DiagField

Main(
    geometry="1Dcartesian",

    interpolation_order=2,

    timestep=dt,
    simulation_time=tsim,

    cell_length=[dx],
    grid_length=[Lsim],

    number_of_patches=[8],

    EM_boundary_conditions=[
                               ['periodic'],
                               ],

    random_seed=smilei_mpi_rank
)

# APPLY EXTERNAL FIELD

ExternalField(
	field='Bz',
    profile=constant(B0)
)

Species(
    name='electron',
    position_initialization='regular',
    momentum_initialization='maxwell-juettner',
    particles_per_cell=nppc,
    mass=me,
    charge=-e,
    number_density=n,
    temperature=[Te],
    boundary_conditions=[
        ['periodic'],
    ],
)

Species(
    name='cold_ion',
    position_initialization='regular',
    momentum_initialization='maxwell-juettner',
    particles_per_cell = nppc,
    number_density=0.98*n,
    temperature=[Tp],
    mass=mp,
    charge=e,
    boundary_conditions=[
        ['periodic'],
    ],
)

Species(
    name='ring_ion',
    position_initialization=po_init_w,
    momentum_initialization=mo_init,
   #particles_per_cell = nppc,
   #temperature
    number_density=0.02*n,
    mass=mp,
    charge=e,
    boundary_conditions=[
        ['periodic'],
    ],
)

DiagScalar(
    every=100
)


DiagFields(
    every=diagEvery,
    fields=['Ex','Ey','Ez']
)

# DiagTrackParticles(
#     species = "ion",
#     every = diagEvery
# )

