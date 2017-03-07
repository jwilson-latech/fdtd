import fdtdcl as fd
import numpy as np
import os
import pyopencl as cl
from numpy import pi, sqrt, cosh, exp, cos, sin
import matplotlib.pyplot as plt


def run_soliton(width,steps,snapshot,name):
    WIDTH = width
    ''' 
    Create space.
    '''
    a = 20
    xx = np.linspace(-a,a,WIDTH+1)
    dx = xx[1]-xx[0]
    dy = xx[1]-xx[0]
    xx = xx[:-1]
    sigma = 1/10.0
    dt=sigma*dx**2
    X, Y = np.meshgrid(xx, xx, sparse=True)
    ''' 
  49  Set intial conditions.
    '''
    d1=45 # direction normal to curve in deg
    d2=0 # direction of propagation w.r.t curve in deg
    d2=d2+d1 # direction of propagation in deg

    t1=d1/180.0*pi # direction normal to curve in rad
    t2=d2/180.0*pi # direction normal to curve in rad
    s = 2.0 # parameter. Simulation unstable for p<1
    lamb = -4.0 # nonlinear coupling constant
    k1 = np.sqrt(-lamb/2.0) # pulse width
    k2 = s*k1 # wave number
    w1 = s*k1**2*abs(np.cos(t1-t2)) # pulse 'frequency'
    w2 = k1**2*(s**2-1) #frequency of phase
    A0 = sqrt(2/(-lamb))*k1 #amplitude
    z0=50/sqrt(2)
    k1x = k1*cos(t1)
    k1y = k1*sin(t1)
    k2x = k2*cos(t2)
    k2y = k2*sin(t2)
    k1z0 = k1*z0
    k2z0 = k2*z0

    def exact(x,y,t):
        ''' 
        This is the function for the initial conditions. 
        '''
        z0=10
        return exp(1j*(2*(x-z0))-w1*t)*exp(-((x-z0)*(x-z0)+y*y)/9)+exp(-1j*(2*(x+z0))-w1*t)*exp(-((x+z0)*(x+z0)+y*y)/9)
        
    def mox1(x,y):
        ''' 
        This is the function for the initial conditions. 
        '''
        return exp(-1j*(2*x+2*y-2*z0))/cosh((x+y-z0))
        
    def mox2(x,y):
        ''' 
        This is the function for the initial conditions. 
        '''
        return exp(-1j*(2*x+2*y-2*z0-3*dt))/cosh((x+y-z0)-4*dt)
        
    
    sScale = 1
    rScale = 0
    u0 = sScale*np.float64(exact(X,Y,0).real      + rScale*(np.random.rand(WIDTH,WIDTH)*exp(1j*2*pi*np.random.rand(WIDTH,WIDTH)-0.5)))
    v0 = sScale*np.float64(exact(X,Y,0).imag      + rScale*(np.random.rand(WIDTH,WIDTH)*exp(1j*2*pi*np.random.rand(WIDTH,WIDTH)-0.5)))
    u1 = sScale*np.float64(exact(X,Y,0.5*dt).real + rScale*(np.random.rand(WIDTH,WIDTH)*exp(1j*2*pi*np.random.rand(WIDTH,WIDTH)-0.5)))
    v1 = sScale*np.float64(exact(X,Y,0.5*dt).imag + rScale*(np.random.rand(WIDTH,WIDTH)*exp(1j*2*pi*np.random.rand(WIDTH,WIDTH)-0.5)))
    
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
    global_vars = global_vars + "\n __constant double A0 = %s;" %(A0)
    global_vars = global_vars + "\n __constant double z0 = %s;" %(z0)
    global_vars = global_vars + "\n __constant double k1x = %s;" %(k1x)
    global_vars = global_vars + "\n __constant double k1y = %s;" %(k1y)
    global_vars = global_vars + "\n __constant double k2x = %s;" %(k2x)
    global_vars = global_vars + "\n __constant double k2y = %s;" %(k2y)
    global_vars = global_vars + "\n __constant double k1 = %s;" %(k1)
    global_vars = global_vars + "\n __constant double k2 = %s;" %(k2)
    global_vars = global_vars + "\n __constant double w1 = %s;" %(w1)
    global_vars = global_vars + "\n __constant double w2 = %s;" %(w2)
    global_vars = global_vars + "\n __constant double k1z0 = %s;" %(k1z0)
    global_vars = global_vars + "\n __constant double k2z0 = %s;" %(k2z0)
    
    print global_vars

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
    #create the wave objects
    wave_abc = fd.Field(env,[u0,u1],[v0,v1],program, update_function = gfdtd_abc,multistep_modifier=4)  
    
    field_list = [wave_abc] #Simulation needs list of fields
    for wave in field_list:
        wave.field_from_device()#take the field from the device 
        
    sim = fd.Simulation(field_list,dt=0.1,dx=1) #make the simluation
    
    prg = program.build()

    queue = env.queue
    
    data_abc  = np.array( [ wave_abc.values_real[1] +  1j * wave_abc.values_imag[1] ] )

    for iter in range(steps):
        if (iter % snapshot == 0):
            for wave in field_list:
                wave.field_from_device() #take the field from the device
            thewave_abc = np.array( [ wave_abc.values_real[1] + 1j*wave_abc.values_imag[1] ] )
            data_abc = np.append(data_abc,thewave_abc, axis=0)
            print iter*dt
            None
        sim.update()
    data_abc  = np.delete(data_abc,0,0)
    file_abc  = name + "abc"+ ".npy"
    np.save(file_abc,data_abc)
    return data_abc

def run_gaussian(width,steps,snapshot,name):
    WIDTH = width
    ''' 
    Create space.
    '''
    a = 40
    xx = np.linspace(0,a,WIDTH+1)
    dx = xx[1]-xx[0]
    xx = xx[:-1]
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
    s = 2.0 # parameter. Simulation unstable for p<1
    lamb = -4.0 # nonlinear coupling constant
    k1 = np.sqrt(-lamb/2.0) # pulse width
    k2 = s*k1 # wave number
    w1 = s*k1**2*abs(np.cos(t1-t2)) # pulse 'frequency'
    w2 = k1**2*(s**2-1) #frequency of phase
    A0 = sqrt(2/(-lamb))*k1 #amplitude
    z0=50
    print w1,w2

    def exact(x,y,t):
        ''' 
        This is the function for the initial conditions. 
        '''
        return A0*exp(-1j*(k2*(x*cos(t2)+y*sin(t2)-z0)-w2*t))/cosh(k1*(x*cos(t1)+y*sin(t1)-z0)-w1*t)
        
    def mox(x,y,t):
        ''' 
        This is the function for the initial conditions. 
        '''
        return exp(-1j*(2*(x+y-z0)-w2*t))/cosh((x+y-z0)-w1*t)
        
    
    u0 = np.float64(mox(X,Y,0).real)
    v0 = np.float64(mox(X,Y,0).imag)
    u1 = np.float64(mox(X,Y,0.5*dt).real)
    v1 = np.float64(mox(X,Y,0.5*dt).imag)
    
    '''
    Setting up the environment. 
    '''
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1' #supress outputs
    env = fd.ComputeNode(0,1) # specifies which GPU to use
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
    
    print global_vars

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
    #create the wave objects
    wave  = fd.Field(env,[u0,u1],[v0,v1],program, update_function = gfdtd_abc,multistep_modifier=3) 
    wave2 = fd.Field(env,[u0,u1],[v0,v1],program, update_function = gfdtd_ebc, multistep_modifier=3) 
    wave3 = fd.Field(env,[u0,u1],[v0,v1],program, update_function = gfdtd_nbc, multistep_modifier=3) 

    wave.field_from_device() #take the field from the device
    field_list = [wave,wave2,wave3] #Simulation needs list of fields
    sim = fd.Simulation(field_list,dt=0.1,dx=1) #make the simluation
    
    prg = program.build()

    queue = env.queue

    data  = np.array( [ wave.values_real[1] +  1j * wave.values_imag[1] ] )
    data2 = np.array( [ wave2.values_real[1] + 1j *wave2.values_imag[1] ] )
    data3 = np.array( [ wave3.values_real[1] + 1j *wave3.values_imag[1] ] )

    for iter in range(steps):
        if (iter % snapshot == 0):
            #data = np.stack(data)
            #data2 = np.stack(data2)
            wave.field_from_device()
            wave2.field_from_device()
            wave3.field_from_device()
            thewave = np.array( [ wave.values_real[1] + 1j*wave.values_imag[1] ] )
            thewave2 = np.array( [ wave2.values_real[1] + 1j*wave2.values_imag[1] ] )
            thewave3 = np.array( [ wave3.values_real[1] + 1j*wave3.values_imag[1] ] )
            data = np.append(data,thewave, axis=0)
            data2 = np.append(data2,thewave2, axis=0)
            data3 = np.append(data3,thewave3, axis=0)
            print iter*dt
            None
        sim.update()
    data  = np.delete(data,0,0)
    data2 = np.delete(data2,0,0)
    data3 = np.delete(data3,0,0)
    file  = name + "abc"+ ".npy"
    file2 = name + "ebc"+ ".npy"
    file3 = name + "ebc"+ ".npy"
    np.save(file,data)
    np.save(file2,data2)
    np.save(file3,data3)
    
    return data, data2, data3
    
def gfdtd_abc(object):
    """
    This is the update function for the G-FDTD scheme with absorbing boundary conditions. 
    """
    op_A = object.prg.operator_A
    mover = object.prg.operator_mover
    gfdtd = object.prg.gfdtd
    bc = object.prg.abc
    restore = object.prg.restore
    
    U = object.values_real_dev
    V = object.values_imag_dev
    ctx = object.compute_node.context
    queue = object.compute_node.queue
    # need to calculate A1_ and A3_. A_2 is a dummy
    for substep in [2,3]:
        substep = np.int64(substep)
        mover( queue, object.shape, None, U, V, substep) 
        op_A( queue, object.shape, None, U, V, substep)
        # perform GFDTD calculations
        gfdtd( queue, object.shape, None, U, V,substep)
        # reconstruct boundaries
        bc(queue, object.shape,None,U,V,substep)
        
    # set up for next iteration
    restore(queue, object.shape, None,U,V)
    return None
    
    
def gfdtd_ebc(object):
    """
    This is the update function for the G-FDTD scheme with absorbing boundary conditions. 
    """
    op_A = object.prg.operator_A
    mover = object.prg.operator_mover
    gfdtd = object.prg.gfdtd
    bc = object.prg.ebc
    restore = object.prg.restore
    
    U = object.values_real_dev
    V = object.values_imag_dev
    ctx = object.compute_node.context
    queue = object.compute_node.queue
    # need to calculate A1_ and A3_. A_2 is a dummy
    for substep in [2,3]:
        substep = np.int64(substep)
        thestep = np.int64(object.thestep)
        mover( queue, object.shape, None, U, V, substep) 
        op_A( queue, object.shape, None, U, V, substep)
        # perform GFDTD calculations
        gfdtd( queue, object.shape, None, U, V,substep)
        # reconstruct boundaries
        bc(queue, object.shape,None,U,V,substep,thestep)
        
    # set up for next iteration
    restore(queue, object.shape, None,U,V)
    return None
    
def gfdtd_nbc(object):
    """
    This is the update function for the G-FDTD scheme with absorbing boundary conditions. 
    """
    op_A = object.prg.operator_A
    mover = object.prg.operator_mover
    gfdtd = object.prg.gfdtd
    bc = object.prg.nbc
    restore = object.prg.restore
    
    U = object.values_real_dev
    V = object.values_imag_dev
    ctx = object.compute_node.context
    queue = object.compute_node.queue
    # need to calculate A1_ and A3_. A_2 is a dummy
    for substep in [2,3]:
        substep = np.int64(substep)
        thestep = np.int64(object.thestep)
        mover( queue, object.shape, None, U, V, substep) 
        op_A( queue, object.shape, None, U, V, substep)
        # perform GFDTD calculations
        gfdtd( queue, object.shape, None, U, V,substep)
        # reconstruct boundaries
        bc(queue, object.shape,None,U,V,substep,thestep)
        
    # set up for next iteration
    restore(queue, object.shape, None,U,V)
    return None

    
    




    
    
    
    
    
    
    
    
    
    
    
    
    