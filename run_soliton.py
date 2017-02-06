import fdtdcl as fd
import numpy as np
import os
import pyopencl as cl
from numpy import pi, sqrt, cosh, exp, cos, sin
import matplotlib.pyplot as plt


def runSim(width,steps,snapshot,name):
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
    wave  = fd.Field(env,[u0,u1],[v0,v1],program, update_function = gfdtd_abc,multistep_modifier=3) #create the wave object
    wave2 = fd.Field(env,[u0,u1],[v0,v1],program, update_function = gfdtd_ebc, multistep_modifier=3) 
    wave3 = fd.Field(env,[u0,u1],[v0,v1],program, update_function = gfdtd_nbc, multistep_modifier=3) 
    #print wave.update_function

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
        op_A( queue, object.shape, None, U, V, substep)
        #op_A(queue, object.shape, None, A1_U, A1_V, A2_U, A2_V,substep)
        #op_A(queue, object.shape, None, A2_U, A2_V, A3_U, A3_V,substep)
        # perform GFDTD calculations
        gfdtd( queue, object.shape, None, U, V,substep)
        # reconstruct boundaries
        bc(queue, object.shape,None,U,V,substep)
        
    # set up for next iteration
    restore(queue, object.shape, None,U,V)
    return None
    
def gfdtd_ebc(object):
    """
    This is the update function for the G-FDTD scheme with exact boundary conditions. 
    """
    op_A = object.prg.operator_A
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
        op_A( queue, object.shape, None, U, V, substep)
        #op_A(queue, object.shape, None, A1_U, A1_V, A2_U, A2_V,substep)
        #op_A(queue, object.shape, None, A2_U, A2_V, A3_U, A3_V,substep)
        # perform GFDTD calculations
        gfdtd( queue, object.shape, None, U, V,substep)
        # reconstruct boundaries
        bc(queue, object.shape,None,U,V,substep,np.int64(object.thestep))
        
    # set up for next iteration
    restore(queue, object.shape, None,U,V)
    return None

def gfdtd_nbc(object):
    """
    This is the update function for the G-FDTD scheme with exact boundary conditions. 
    """
    op_A = object.prg.operator_A
    gfdtd = object.prg.gfdtd
    #bc = object.prg.ebc
    restore = object.prg.restore
    
    U = object.values_real_dev
    V = object.values_imag_dev
    ctx = object.compute_node.context
    queue = object.compute_node.queue
    # need to calculate A1_ and A3_. A_2 is a dummy
    for substep in [2,3]:
        substep = np.int64(substep)
        op_A( queue, object.shape, None, U, V, substep)
        #op_A(queue, object.shape, None, A1_U, A1_V, A2_U, A2_V,substep)
        #op_A(queue, object.shape, None, A2_U, A2_V, A3_U, A3_V,substep)
        # perform GFDTD calculations
        gfdtd( queue, object.shape, None, U, V,substep)
        # reconstruct boundaries
        #bc(queue, object.shape,None,U,V,substep,np.int64(object.thestep))
        
    # set up for next iteration
    restore(queue, object.shape, None,U,V)
    return None




    
    
    
    
    
    
    
    
    
    
    