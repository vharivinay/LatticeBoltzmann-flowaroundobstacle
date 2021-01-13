'''
Lattice Boltzmann methods (LBM), originated from the lattice gas 
automata (LGA) method (Hardy-Pomeau-Pazzis and Frisch-Hasslacher-Pomeau 
models), is a class of computational fluid dynamics (CFD) methods for
fluid simulation. Instead of solving the Navier–Stokes equations directly,
a fluid density on a lattice is simulated with streaming and collision 
(relaxation) processes. The method is versatile as the model fluid 
can straightforwardly be made to mimic common fluid behaviour like 
vapour/liquid coexistence, and so fluid systems such as liquid droplets 
can be simulated. 

More Information: https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods

Based on lectures: 
Simulation and modeling of natural processes
Week 5: Lattice Boltzmann modeling of fluid flow
Authoor: Jonas Latt
Universite de Geneve, Switzerland

Course link: 
https://www.coursera.org/learn/modeling-simulation-natural-processes

Modified by: Harivinay V
Github: https://www.github.com/M87K452b

--- Code restructured only a bit. Contains additional 
--- Added functionality to export a video directly from code

Dependencies: 
    ffmpeg - install using pip install ffmpeg-python
    graphviz - install using pip install graphviz

'''

from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
from numba import jit

## ADD PATH TO FOLDER FOR SAVING IMAGES 
import os
my_path = os.path.dirname(__file__)

## Flow definition 
maxIter = 30000  # Total number of time iterations.
Re = 220.0         # Reynolds number.
nx, ny = 420, 180 # Numer of lattice nodes.
ly = ny-1         # Height of the domain in lattice units.
cx, cy, r = nx//4, ny//2, ny//9 # Coordinates of the cylinder.
uLB     = 0.04                  # Velocity in lattice units.
nulb    = uLB*r/Re;             # Viscoscity in lattice units.
omega = 1 / (3*nulb+0.5);    # Relaxation parameter.

## Lattice Constants 
v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])
t = array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

col1 = array([0, 1, 2])
col2 = array([3, 4, 5])
col3 = array([6, 7, 8])

## Function Definitions 
@jit(nopython=True)
def macroscopic(fin):
    '''
    Computes density and velocity distribution
    '''
    rho = sum(fin, axis=0)
    u = zeros((2, nx, ny))
    for i in range(9):
        u[0,:,:] += v[i,0] * fin[i,:,:]
        u[1,:,:] += v[i,1] * fin[i,:,:]
    u /= rho
    return rho, u

def equilibrium(rho, u):
    '''
    Computes Equilibrium distribution.
    '''            
    usqr = 3/2 * (u[0]**2 + u[1]**2)
    feq = zeros((9,nx,ny))
    for i in range(9):
        cu = 3 * (v[i,0]*u[0,:,:] + v[i,1]*u[1,:,:])
        feq[i,:,:] = rho*t[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq

## Setup: cylindrical obstacle and velocity inlet with perturbation ########
# Creation of a mask with 1/0 values, defining the shape of the obstacle.

def obstacle_fun(x, y):
    '''
    This function creates the obstacle.
    This can be modified to get any shape desired
    
    Default: 
        -cylinder 
        (x-xc)**2  + (yc-y)**2 < r*r
        -ellipse
        (x-cx)**2/8 + (cy-y)**2/3 < r*(11/2)
        
    '''
    return (x-cx)**2/8 + (cy-y)**2/3 < r*(11/2) 

obstacle = fromfunction(obstacle_fun, (nx,ny))

# Initial velocity profile: almost zero, with a slight perturbation to trigger
# the instability.
def inivel(d, x, y):
    return (1-d) * uLB * (1 + 1e-4*sin(y/ly*2*pi))

vel = fromfunction(inivel, (2,nx,ny))

# Initialization of the populations at equilibrium with the given velocity.
fin = equilibrium(1, vel)

## Video Visuals - Saving image at desired timesteps
fig = plt.figure()
#list_im = [] # Empty list to store images

## Main time loop
for time in range(maxIter+1):
    # Right wall: outflow condition.
    fin[col3,-1,:] = fin[col3,-2,:] 

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Left wall: inflow condition.
    u[:,0,:] = vel[:,0,:]
    rho[0,:] = 1/(1-u[0,0,:]) * ( sum(fin[col2,0,:], axis=0) +
                                  2*sum(fin[col3,0,:], axis=0) )
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # Collision step.
    fout = fin - omega * (fin - feq)

    # Bounce-back condition for obstacle.
    for i in range(9):
        fout[i, obstacle] = fin[8-i, obstacle]

    # Streaming step.
    for i in range(9):
        fin[i,:,:] = roll(
                            roll(fout[i,:,:], v[i,0], axis=0),
                            v[i,1], axis=1 )
 
    # Visualization of the velocity.
    if (time%100==0):
        print('time = '+ str(time))
        plt.clf()
        plt.tight_layout()
        plt.title('time = {}s'.format(time))
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.jet, animated=True)
        plt.savefig(my_path + "/Re_test/vel.{0:04d}.png".format(time//100))
        #list_im.append([image_frame])

import ffmpeg
(
    ffmpeg
    .input('Re_test/*.png', pattern_type='glob', framerate=24)
    .filter('deflicker', mode='pm', size=10)
    .filter('scale', size='hd1080', force_original_aspect_ratio='increase')
    .output('ellipse_1080.mp4', crf=20, preset='slower', movflags='faststart', pix_fmt='yuv420p')
    .run()
)

