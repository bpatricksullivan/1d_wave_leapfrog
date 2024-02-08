import numpy as np
import matplotlib.pyplot as plt

#Initialize physical and numerical constants

a = 1.0  #set the length of the system
dx = 0.1 #set the discrete spatial stepsize
c = 1.0  #define the wave speed

dt = dx/c   #choose a time step to satisfy the CFL condition. Information can't
            #travel further than dx during a time dt or the system will be
            #numerically unstable. 

x = np.arange(0,a*(1+dx),dx) #define an array to store x position data
npts = len(x) # this is the number of spatial points along x
nsteps = 100  #set the number of time steps

f = np.zeros((npts,3))

xc = 0.5 #define the center of the system to locate a Gaussian pulse (see below)
w = 0.05 #define the width of the Gaussian wave pulse

#f[:,0] = np.exp(-(x-xc)**2/(w**2)) #use this initial condition for a Gaussian
f[:,0] = np.sin(2*np.pi*x/a) #use this initial condition for a standing wave

#first time step in the leap frog algorithm
f[1:-1,1] = f[1:-1,0] + .5*c**2*(f[:-2,0]+f[2:,0]-2.*f[1:-1,0])*(dt**2/dx**2) 

#for all additional time steps
for k in range(0,nsteps):
    #for i in range(1,npts-1):
    f[1:-1,2] = -f[1:-1,0] + 2*f[1:-1,1] \
        + c**2*(f[:-2,1]+f[2:,1]-2.*f[1:-1,1])*(dt**2/dx**2) 
    
    #push the data back for the leapfrogging
    f[:,0] = f[:,1]
    f[:,1] = f[:,2]

    f[0,:]=0   #enforce the boundary conditions (if necessary)
    f[-1,:]=0

    if k % 1== 0:#use this line to change how often frames are plotted.
                 #e.g. if k%10 ==0 will plot only every 10th timestep.
        
        plt.figure(1)
        plt.clf()
        plt.plot(x,f[:,2], 'b')

        plt.title('t= %0.2f ' % (k*dt) )
        plt.xlim(0,1)
        plt.ylim(-1.5,1.5)        
        #plt.savefig('frame_%03d.png' % k) #use this line to save frames for animation
        plt.show()
        plt.pause(.01)