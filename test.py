import operators as op
import numpy as np
import os
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
env = op.ComputeNode(0,2)

x1=np.zeros([3,3])+1 # 3x3 array of complex numbers. Complex numbers are stored as 2-tuple. 
x2=x1+1
wave = op.Field(env,[x1,x2],[x2,x1])
wave.field_from_device()
field_list = [wave]
step_size_list = [0.1,0.2]
sim = op.Simulation(field_list,dt=0.1,dx=1)

ctx = env.context
global_vars = "\n __constant int WIDTH = 10;"
program = op.CL_Program(ctx, 'functions.cpp', global_vars=global_vars).program
prg = program.build()

prg.solve1(x1,x2)

#can access programs now. Need to work out complex variable type. 

 


