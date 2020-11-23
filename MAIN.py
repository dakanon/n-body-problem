import numpy as np
from Nbody_Problem import nBody
from Display_Results import Display

# Set-up the model
# Units : N, m, kg
dim = 3 
dof = 3
n_body = 3 
time = 0.02
dt = 0.00001

node_crd = np.array([[-0.0555, 0.0, 0.0],            # body 1 - [x y z]
                     [0.0, 0.0, 0.0],            # body 2 - [x y z]
                     [10.0, -1, 10.0]])           # body 3 - [x y z]

node_mass = np.array([3.4e2, 3.4e10, 6.8e3])       # body mass ([1]] [2]  [3])
 
node_bc = np.array([[0, 0.0, 0, 0, 0, 0, 0, 0, 0],     # body 1 - [dof 1 - acc, vel, disp dof 2 ...]
                    [0.0, 0, 0, 0, 0, 0, 0, 0, 0],     # body 2 - [dof 1 - acc, vel, disp dof 2 ...] 
                    [0, 0.0, 0, 0, 0, 0, 0, 0, 0]])    # body 3 - [dof 1 - acc, vel, disp dof 2 ...]


#node_propulsion = np.loadtxt('prop_profile.txt')
node_propulsion = np.zeros(((int(time/dt)),(int(n_body*dof))))
# Solve the N-body Problem
#step 1: construct the model inside the domain
FEM = nBody(dim, dof, n_body, time, dt, node_crd, node_mass, node_bc, node_propulsion)
DISP = Display(dim, dof, n_body, time, dt, node_crd.copy(), node_mass.copy())
del dim, dof, n_body, dt, node_crd, node_mass, node_bc, node_propulsion
#step 2: assemble initial matrix-vector equation
FEM.assemble()
#analysis 
while (FEM.t < time): #while current time is less than analysis time
    #step 3: solve the equation of motion
    FEM.solve(True) # monitor progress
    #step 4: update internal parameters
    FEM.update()
       
# Return Solution (over time)
#acceleration = FEM.node.arec
#velocity = FEM.node.vrec
#displacement = FEM.node.drec
#step 5: plot the current domain
DISP.play_domain(FEM.node.drec, FEM.node.dmstore, FEM.nstep, FEM.dt)
