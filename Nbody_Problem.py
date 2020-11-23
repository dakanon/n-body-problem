import math
import numpy as np

class nBody:
    def __init__(self, dim, dof, n_body, time, dt, node_crd, node_mass, node_bc, node_propulsion):
        ## User input array defintions
        ## node_crd = np.zeros((N_body,DOF))        # x, y, z global
        ## node_bc = np.zeros((3*DOf,N_body))       # acc, vel, disp
        ## node_mass = np.zeros((N_body))           # m1, m2, m3, ...
        ## node_propulsion = np.zeros(((int(time/dt)),(int(n_body*dof)))) or .txt input
        ## node mass can also be defined dynamically by using the variable node.dmstore
        ## change in nodal mass at each step can be passed to dmstore (leaving the first row as 0)
        print("Model building started...")
        # nBody Class variables
        self.DIM = dim                          # 3D space
        self.DOF = dof                          # only translation
        self.N_body = n_body                    # number of bodies
        self.N_dof = dof*n_body                 # number of d.o.f.s
        self.time = time                        # analysis time
        self.dt = dt                            # time step
        self.t = 0                              # current time
        self.step = int(0)                      # current step
        self.nstep = int(time/dt)               # number of steps
        self.node = self.Node(n_body, dof, node_crd, node_mass, node_bc, node_propulsion, self.nstep)
        self.tensor = self.Tensor(self.DOF, self.N_body, self.N_dof, self.node.bc)
        self.determine = self.State_Determination()
        self.method = self.Numerical_Methods()
        print("Model built...")
        
                
    # nBody Class sub-classes
    class Node:
        def __init__(self, N_body, DOF, node_crd, node_mass, node_bc, node_propulsion, nstep):
            #nodal info
            self.crd = node_crd.copy()                                 # node coordinates x, y, z
            self.bc = node_bc.copy()                                   # node boundary conditions
            self.m = node_mass.copy()                                  # body mass vector
            #recorders & storage
            self.dmstore = np.zeros((N_body, nstep + 2))        # change in mass dm1, dm2, dm3, ...
            self.pstore = np.zeros((N_body, DOF, nstep + 2))    # propulsion px, py, pz
            self.grec = np.zeros((N_body, DOF, nstep+1))          # pull gx, gy, gz
            self.arec = np.zeros((N_body, DOF, nstep+1))          # acceleration ax, ay, az
            self.vrec = np.zeros((N_body, DOF, nstep+1))          # velocity ax, ay, az
            self.drec = np.zeros((N_body, DOF, nstep+1))          # displacement ax, ay, az
            #fill propulsion recorder
            node_propulsion = node_propulsion.reshape(nstep, N_body, DOF)
            for i in range(nstep):
                self.pstore[:,:,i] = node_propulsion[i,:,:]
                   
                
    class Tensor:
        def __init__(self, DOF, N_body, N_dof, BC):
            #matrix-vector equation
            self.m = np.zeros((N_dof, N_dof))       # mass matrix mx, my, mz
            self.a = np.zeros((N_dof))              # acceleration vector
            self.a_commit = np.zeros((N_dof))          # converged acceleration vector
            self.v = np.zeros((N_dof))              # velocity vector
            self.v_commit = np.zeros((N_dof))          # converged velocity vector
            self.d = np.zeros((N_dof))              # displacement vector
            self.d_commit = np.zeros((N_dof))          # converged displacement vector
            self.g = np.zeros((N_dof))              # gravitational acceleration vector
            self.p = np.zeros((N_dof))              # propulsion acceleration vector
            #fill-in boundary conditions
            for i in range(N_body): #for each body
                for j in range(DOF): #for each d.o.f.
                    self.a_commit[j+(i*DOF)] = BC[i, (j*DOF)]
                    self.v_commit[j+(i*DOF)] = BC[i, (j*DOF+1)]
                    self.d_commit[j+(i*DOF)] = BC[i, (j*DOF+2)]
    
                    
    class State_Determination:
        def __init__(self):
            self.gravity_field = self.Gravity_Field()
            
        # State_Determination Class sub-class
        class Gravity_Field:
            def radius(self, DIM, nd1, nd2):
                dist = 0
                dir_vect = np.zeros(DIM)
                for i in range(DIM):
                    dir_vect[i] = nd2[i] - nd1[i]
                    dist = dist + ((dir_vect[i])**2)
                
                dist = math.sqrt(dist)
                if dist < 1e-8:
                    dist = 0
                    dir_vect = np.zeros(DIM)
                else:
                    dir_vect = dir_vect/dist
                
                return dist, dir_vect
        
        
            def pull(self, mass, radius):
                # Constants
                G = 6.67*10**(-11)                 # Gravitational constant (Nm2/kg2)
                #initialize
                gacc = 0
                if radius < 1e-8:
                    pass
                else:
                    gacc = (-1*G)*(mass/(radius**2))
            
                return gacc
        

        # State_Determination Class functions
        def pull_state(self, N_body, DOF, DIM, nd_crd, nd_m):
            #initialize
            ps = np.zeros((N_body*DOF))
            #determine
            for i in range(N_body): #for each body
                PULL = np.zeros(DOF)
                #integrate over the gravity field
                for j in range(N_body): #due to each body
                    r_mag, r_hat = self.gravity_field.radius(DIM, (nd_crd[i]), (nd_crd[j]))
                    PULL = PULL + self.gravity_field.pull((nd_m[j]), r_mag) * r_hat
                
                #put it in a vector
                for j in range(DOF):
                    ps[j+(i*DOF)] = PULL[j]
                    
            #return
            return ps.copy()
            
            
        def mass_state(self, N_body, DOF, nd_m):
            #initialize
            ms = np.zeros(((N_body*DOF),(N_body*DOF)))
            #determine
            for i in range(N_body): #for each body
                for j in range(DOF): #for each d.o.f.
                    ms[(j+(i*DOF)), (j+(i*DOF))] = nd_m[i] #diagonal
                            
            #return
            return ms.copy()
            
            
        def prop_state(self, N_body, DOF, nd_pstore, step):
            #initialize
            prs = np.zeros((N_body*DOF))
            #determine
            for i in range(N_body): #for each body
                for j in range(DOF): #for each d.o.f.
                    prs[j+(i*DOF)] = nd_pstore[i, j, step]
                            
            #return
            return prs.copy()
    

    class Numerical_Methods:
        def Newmark(self, algorithm):
            #pick integration algorithm
            if ((algorithm == 'Central Difference') or (algorithm == 'CD')):
                beta = 0.0
                gamma = 0.5
            elif ((algorithm == 'Linear Acceleration') or (algorithm == 'LA')):
                beta = 0.1667
                gamma = 0.5 
            elif ((algorithm == 'Constant Acceleration') or (algorithm == 'CA')):
                beta = 0.25
                gamma = 0.5  
            else: #use constant acceleration -> uncoditionally stable
                beta = 0.25
                gamma = 0.5
                   
            return beta, gamma
            
        
        def L2Norm(self, new, old):
            #initialize
            norm = 0.0
            dof = np.size(new)
            error = np.zeros(dof)
            #compute error
            for i in range(dof):
                error[i] = abs((new[i]-old[i])/(new[i]))
                norm += (error[i]**2)
            
            #compute L2-norm
            norm = math.sqrt(norm)
            #retun norm
            return norm
    
    
    # nBody Class methods      
    def solve(self, progress):
        #numerical solution parameters
        #max_it = 10
        #convergence_criteria = 1e-3
        norm = 0
        n_it = 1
        #no_convergence = True
        #initiate iterative parameters
        acc_it = np.zeros_like(self.tensor.a)
        vel_it = np.zeros_like(self.tensor.v)
        disp_it = np.zeros_like(self.tensor.d)
        vel_tilda = np.zeros_like(self.tensor.v)
        disp_tilda = np.zeros_like(self.tensor.d)
        
        # ---Start Solution---
        #for acc_it computation later
        #advance virtual time to new point
        virtual_step = self.step + 1    #step n+1 (for backward (implicit) mapping)
        #keep track of change in mass with time
        virtual_mass = self.node.m.copy()
        virtual_mass += self.node.dmstore[:, virtual_step].copy()
        #Now,
        #solve Mü(t) = Mg(d(t)) + Ma(t) with Newmark
        algorithm = 'Constant Acceleration'
        beta, gamma = self.method.Newmark(algorithm)
        #compute predictors
        vel_tilda = self.tensor.v_commit.copy() + ((1-gamma)*self.dt*self.tensor.a_commit.copy())
        disp_tilda = self.tensor.d_commit.copy() + (self.dt*self.tensor.v_commit.copy()) + \
            ((0.5*self.dt*self.dt)*(1-(2*beta))*self.tensor.a_commit.copy())
        #initialize virtual position
        virtual_crd = self.node.crd.copy()
        #advance virtual position
        for i in range(self.N_body):
            for j in range(self.DOF):
                virtual_crd[i, j] += disp_tilda[j+(i*self.DOF)].copy()
                
        #compute acceleration at new point
        self.tensor.g = self.determine.pull_state(self.N_body, self.DOF, self.DIM, virtual_crd, virtual_mass)
        self.tensor.p = self.determine.prop_state(self.N_body, self.DOF, \
                                                                   self.node.pstore, virtual_step)
        acc_tilda = -1*(self.tensor.g + self.tensor.p)
        #compute predictor + corrector
        vel_tilda = vel_tilda + (gamma*self.dt*acc_tilda)
        disp_tilda = disp_tilda + (beta*self.dt*self.dt*acc_tilda)    
        #iterative solution
        # while no_convergence:
        #     #reset virtual position
        #     virtual_crd = self.node.crd.copy()
        #     #advance virtual position
        #     for i in range(self.N_body):
        #         for j in range(self.DOF):
        #             virtual_crd[i, j] += disp_tilda[j+(i*self.DOF)].copy()
                    
        #     #compute acceleration at new point
        #     self.tensor.g = self.determine.pull_state(self.N_body, self.DOF, self.DIM, virtual_crd, virtual_mass)
        #     self.tensor.p = self.determine.prop_state(self.N_body, self.DOF, \
        #                                                                        self.node.pstore, virtual_step)
        #     acc_it = -1*(self.tensor.g + self.tensor.p)
        #     #solve again with trapezoidal
        #     vel_it = self.tensor.v_commit.copy() + 0.5*(self.tensor.a_commit.copy() + acc_it)*self.dt
        #     disp_it = self.tensor.d_commit.copy() + 0.5*(self.tensor.v_commit.copy() + vel_it)*self.dt
        #     #compute L2-norm of displacement
        #     norm = self.method.L2Norm(disp_it.copy(), disp_tilda.copy())
        #     #check convergence & loop control
        #     if (norm < convergence_criteria):
        #         no_convergence = False
        #     elif (n_it >= max_it):
        #         raise Exception('''solve() did not converge with norm: {:.2e} \
        #                         at iteration: {:d}. Current time = {:.4f}.'''. format(norm, n_it, self.t))
        #     else:
        #         n_it += 1
        #         disp_tilda = disp_it
        
        acc_it = acc_tilda
        vel_it = vel_tilda
        disp_it = disp_tilda
        
        self.tensor.a = gamma*(self.tensor.a_commit + acc_it)
        self.tensor.v = vel_it
        self.tensor.d = disp_it
        
        if progress:
            self.monitor(virtual_step, n_it, norm, (virtual_step*self.dt))
        else:
            pass
                 
            
    def update(self):
        #operate in the current step
        #for each body
        for i in range(self.N_body):
            for j in range(self.DOF):
                #update coordinates
                self.node.crd[i, j] += self.tensor.d[j+(i*self.DOF)]    
                #update recorders
                self.node.arec[i, j, self.step + 1] = self.tensor.a[j+(i*self.DOF)]
                self.node.vrec[i, j, self.step + 1] = self.tensor.v[j+(i*self.DOF)]
                self.node.drec[i, j, self.step + 1] = self.tensor.d[j+(i*self.DOF)]
                self.node.grec[i, j, self.step + 1] = self.tensor.g[j+(i*self.DOF)]
        
        #operate in the following time step        
        #update time
        self.step += 1 
        self.t = self.dt * self.step
        #update acc, vel, disp tensors
        self.tensor.a_commit = self.tensor.a
        self.tensor.v_commit = self.tensor.v
        self.tensor.d_commit = self.tensor.d
        #update mass
        self.node.m += self.node.dmstore[:, self.step]                
            
    
    def assemble(self):            
        #compute acceleration bc
        self.tensor.g_commit = self.determine.pull_state(self.N_body, self.DOF, self.DIM, self.node.crd, self.node.m)
        self.tensor.p_commit = self.determine.prop_state(self.N_body, self.DOF, \
                                                                               self.node.pstore, self.step)
        self.tensor.a_commit = -1*(self.tensor.g_commit + self.tensor.p_commit)
        #fill recorders with the b.c.
        for i in range(self.N_body):
            for j in range(self.DOF):    
                self.node.arec[i, j, 0] = self.tensor.a_commit[j+(i*self.DOF)]
                self.node.vrec[i, j, 0] = self.tensor.v_commit[j+(i*self.DOF)]
                self.node.drec[i, j, 0] = self.tensor.d_commit[j+(i*self.DOF)]
                self.node.grec[i, j, 0] = self.tensor.g_commit[j+(i*self.DOF)]
           
    
    def monitor(self, step, n_it, norm, t):     
        print("Step No. {:d} converged with norm: {:.2e} at it.: {:d}. Time: {:.4f}."\
              .format(step, norm, n_it, t))     
        