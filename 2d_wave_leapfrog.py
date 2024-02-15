import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib import cm

#Initialize physical and numerical constants
L = 1.0  #set the length of the system
dx = 0.01 #set the discrete spatial stepsize
c = 1.0  #define the wave speed

dt = .707*dx/c   #choose a time step to satisfy the CFL condition. Information can't
            #travel further than dx during a time dt or the system will be
            #numerically unstable. 

x = np.arange(0,L*(1+dx),dx) #define an array to store x position data
y = np.arange(0,L*(1+dx),dx) #define an array to store x position data

xx,yy = np.meshgrid(x,y) 

npts = len(x) # this is the number of spatial points along x
nsteps = 199  #set the number of time steps

f = np.zeros((npts,npts, 3))

xc = 0.5 #define the center of the system to locate a Gaussian pulse (see below)
w = 0.05 #define the width of the Gaussian wave pulse

#f[:,0] = np.sin(2*np.pi*x/L) #use this initial condition for a standing wave
#f[:,0] = 0.5*np.sin(2*np.pi*x/L) + 0.5*np.sin(3*np.pi*x/L)
f[:,:,0] = np.exp(-(xx-xc)**2/(w**2))*np.exp(-(yy-xc)**2/(w**2)) #use this initial condition for a Gaussian
#f[:int(npts/4),0] = 4*x[:int(npts/4)]
#f[int(npts/4):,0] = -4/3*x[int(npts/4):]+4/3
orig = copy.deepcopy(f[:,0])

#first time step in the leap frog algorithm
f[1:-1,1:-1, 1] = f[1:-1,1:-1, 0] + .5*c**2*(f[:-2,1:-1, 0]+f[2:,1:-1,0]-2.*f[1:-1,1:-1,0])*(dt**2/dx**2)\
                                  + .5*c**2*(f[1:-1,:-2, 0]+f[1:-1,2:,0]-2.*f[1:-1,1:-1,0])*(dt**2/dx**2)

#for all additional time steps
for k in range(0,nsteps):
    #for i in range(1,npts-1):
    f[1:-1,1:-1,2] = -f[1:-1,1:-1,0] + 2*f[1:-1,1:-1,1] \
        + c**2*(f[:-2,1:-1, 1]+f[2:,1:-1,1]-2.*f[1:-1,1:-1,1])*(dt**2/dx**2)\
        + c**2*(f[1:-1,:-2, 1]+f[1:-1,2:,1]-2.*f[1:-1,1:-1,1])*(dt**2/dx**2)
    
    #push the data back for the leapfrogging
    f[:,:,0] = f[:,:,1]
    f[:,:,1] = f[:,:,2]


    if k % 1 == 0:#use this line to change how often frames are plotted.
                 #e.g. if k%10 ==0 will plot only every 10th timestep.
        
        fig = plt.figure(1)
        fig.clf()
        ax = fig.add_subplot(projection = '3d')
        ax.plot_surface(xx,yy, f[:,:,2], rstride=1, cstride=1, cmap = cm.coolwarm)
        ax.plot_wireframe(xx, yy, f[:,:,2], rstride=10, cstride=10, color='green')
        plt.title('t= %0.2f ' % (k*dt) )
        ax.set_zlim3d(-.25,1)       
        #plt.savefig('frame_%03d.png' % k) #use this line to save frames for animation
        plt.show()
        plt.pause(.001)