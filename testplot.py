from test import runSim
from plotter import plotter
import os
import numpy as np


NAME = 'soln2dabc'
SIZE = 100
MAXITER = 100
PICS = 5
IMAGENO = MAXITER/PICS

runSim(SIZE,MAXITER,IMAGENO,NAME)
data = np.load(NAME+".npy")
plotter(data,NAME)
os.system('open *.eps')