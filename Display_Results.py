import numpy as np
import matplotlib.pyplot as plt

class Display:
    def __init__(self, dim, dof, n_body, time, dt, node_crd, node_mass):
        self.DIM = dim                          # 3D space
        self.DOF = dof                          # only translation
        self.N_body = n_body                    # number of bodies
        self.N_dof = dof*n_body                 # number of d.o.f.s
        self.time = time                        # analysis time
        self.dt = dt                            # time step
        self.t = int(0)                         # current step
        self.nstep = int(time/dt)               # number of steps
        self.crd = node_crd                     # node coordinates x, y, z
        self.m = node_mass                      # body mass vector
        self.locus = np.zeros((n_body, dof, self.nstep + 1))
        self.mass = np.zeros((n_body, self.nstep + 2))
        	                 
    
    def play_domain(self, drec, dmass, nstep, dtime):
        
        def update(i):
            xdata = self.locus[:, 0, i]
            ydata = self.locus[:, 1, i]
            zdata = self.locus[:, 2, i]
            
            domain.set_xdata(xdata)
            domain.set_ydata(ydata)
            domain.set_3d_properties(zdata)
            return domain,
        
        
        #get size of recorder
        size = np.size(drec, axis = 2)
        #update mass and node coordinates
        for i in range(self.nstep + 1):
            for j in range(self.N_body):
                self.mass[j, i] = self.m[j] + dmass[j, i]
                for k in range(self.DIM):
                    self.locus[j, k, i] = self.crd[j, k] + drec[j, k, i]
                    
        #add figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        ax.set_xlim3d(-.5, .5)
        ax.set_ylim3d(-.5, .5)
        ax.set_zlim3d(-.5, .5)
        
        domain, = ax.plot(self.locus[:, 0, 0], self.locus[:, 1, 0], self.locus[:, 2, 0], color = 'red', marker = 'o', linestyle = 'None')
                
        for i in range(size):
            update(i)
            plt.draw()
            plt.pause(0.01)
        
        plt.show()     
            
            
            
            
            
            
            







            

        
        
        
        
        
        
        
        
        
        
        
        
        