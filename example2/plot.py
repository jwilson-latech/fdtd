import numpy as np
import os
from plotter import plotter
NAME = 'soln2d' # what to name the files
#print 'Saving...'
file_abc = NAME + "abc"+ ".npy"
file_ebc = NAME + "ebc"+ ".npy"
file_nbc = NAME + "nbc"+ ".npy"

data_abc = np.load(file_abc)
data_ebc = np.load(file_ebc)
data_nbc = np.load(file_nbc)

print 'Plotting...'
plotter(data_abc, NAME+ "_abc")
plotter(data_ebc, NAME+ "_ebc")
plotter(data_nbc, NAME+ "_nbc")

print 'Plots Open'
os.system('open *.eps')