from run_soliton import runSim
from plotter import plotter
import os
import numpy as np

mult = 2
NAME = 'soln2d'



############################
# NOTHING TO DO BELOW HERE #
############################

print 'Initializing...'
SIZE = mult*100
MAXITER = mult**2*(100)**2/32
print MAXITER
PICS = 5
IMAGENO = MAXITER/PICS

print 'Running...'
data ,data2, data3 = runSim(SIZE,MAXITER+1,IMAGENO,NAME)

#print 'Saving...'
#file  = NAME + "abc"+ ".npy"
#file2 = NAME + "ebc"+ ".npy"
#data = np.load(file)
#data2 = np.load(file2)

print 'Plotting...'
plotter(data ,NAME+ "abc")
plotter(data2,NAME+ "ebc")
plotter(data3,NAME+ "nbc")

print 'Plots Open'
os.system('open *.eps')

print 'Plots Closed'
