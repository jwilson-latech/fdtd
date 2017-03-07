import fdtdcl as fd
import numpy as np
import os
import pyopencl as cl
from numpy import pi, sqrt, cosh, exp, cos, sin
import matplotlib.pyplot as plt
from fullviridis import full_viridis
import matplotlib.colors as colors



def main(width,steps):
    WIDTH = width

    ''' 
    Create space.
    '''
    a = 20
    xx = np.linspace(-a,a,WIDTH) 
    dx = xx[1]-xx[0]
    sigma = 1/10.0
    dt=sigma*dx**2
    X, Y = np.meshgrid(xx, xx, sparse=True)

    ''' 
    Set intial conditions.
    '''
    d1=45 # direction normal to curve in deg
    d2=0 # direction of propagation w.r.t curve in deg
    d2=d2+d1 # direction of propagation in deg

    t1=d1/180.0*pi # direction normal to curve in rad
    t2=d2/180.0*pi # direction normal to curve in rad
    s = 2.0 # parameter. Simulation unstable for p<2
    lamb = -2 # nonlinear coupling constant
    k1 = np.sqrt(-lamb/2.0) # pulse width
    k2 = s*k1 # wave number
    w1 = s*k1**2*abs(np.cos(t1-t2)) # pulse 'frequency'
    w2 = k1**2*(s**2-1) #frequency of phase
    A0 = sqrt(2/(-lamb))*k1 #amplitude
    z0=10 

    def exact(x,y,t):
        ''' This is the function for the initial conditions. 
        '''
        return A0*exp(-1j*(k2*(x*cos(t2)+y*sin(t2)-z0)-w2*t*dt/2.0))/cosh(k1*(x*cos(t1)+y*sin(t1)-z0)-w1*t*dt/2.0)
    
    def exactr(x,y,t):
        ''' Take the real part. 
        '''
        return exact(x,y,t).real
    def exacti(x,y,t):
        ''' Take the imag part. 
        '''
        return exact(x,y,t).imag
    
    u0 = np.float64(exactr(X,Y,0))
    v0 = np.float64(exacti(X,Y,0))
    u1 = np.float64(exactr(X,Y,1))
    v1 = np.float64(exacti(X,Y,1))

    '''
    Setting up the environment. 
    '''
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1' #supress outputs
    env = fd.ComputeNode(0,2) # specifies which GPU to use
    ctx = env.context #shorten this name a bit

    '''
    ABC parameters.
    '''
    cx = (2.0*k2*cos(t2)*dt/dx)
    cy = (2.0*k2*sin(t2)*dt/dx)
    v2x = (k2*cos(t2))**2
    v2y = (k2*sin(t2))**2

    ''' 
    Need to include some global variables in the kernal program.
    '''
    global_vars = ""
    global_vars = global_vars + "\n __constant int WIDTH = %s;" %(WIDTH)
    global_vars = global_vars + "\n __constant int bnd1 = %s;" %(2)
    global_vars = global_vars + "\n __constant int bnd2 = %s;" %(4)
    global_vars = global_vars + "\n __constant double dt = %s;" %(dt)
    global_vars = global_vars + "\n __constant double dx = %s;" %(dx)
    global_vars = global_vars + "\n __constant double dy = %s;" %(dx)
    global_vars = global_vars + "\n __constant double sigma = %s;" %(sigma)
    global_vars = global_vars + "\n __constant double lamb = %s;" %(lamb)
    global_vars = global_vars + "\n __constant double cx = %s;" %(cx)
    global_vars = global_vars + "\n __constant double cy = %s;" %(cy)
    global_vars = global_vars + "\n __constant double v2x = %s;" %(v2x)
    global_vars = global_vars + "\n __constant double v2y = %s;" %(v2y)


    filenames = ['fdm.h', 'regions.h', 'operators.h', 'functions.c']
    # Catenate all the files into a single peice of code.
    with open('kernal_file.c', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
                
    program = fd.CLProgram(ctx, 'kernal_file.c', global_vars=global_vars).program
    
    '''
    Creating simulation.
    ''' 
    wave = fd.Field(env,[u0,u1],[v0,v1],program) #create the wave object

    wave.field_from_device() #take the field from the device
    field_list = [wave] #Simulation needs list of fields
    sim = fd.Simulation(field_list,dt=0.1,dx=1) #make the simluation


    #prg = program.build()
    prg = program.build()

    queue = env.queue

    for i in range(steps):
        sim.update()
    #prg.test(queue, wave.shape, None, wave.values_real_dev, wave.values_imag_dev)
    wave.field_from_device()
    thewave = wave.values_real[1]**2+wave.values_imag[1]**2
    plt.imshow(thewave,interpolation='nearest',cmap=full_viridis,norm=colors.LogNorm(vmin=10**(-5),vmax=1))
    plt.title('width = %s, steps = %s'%(width,steps))
    plt.show()
    

#main(200,4000)




    
    
    
    
    
    
    
    
    
    
    