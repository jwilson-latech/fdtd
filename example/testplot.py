from test import run_soliton
from plotter import plotter
import os
import numpy as np
import sys
import time

mult = 2 # determines how many hundreds points to use in the simulation 
iter_mult = 1
NAME = 'soln2d' # what to name the files



################################
### NOTHING TO DO BELOW HERE ###
################################

print 'Initializing...'
SIZE = np.int64(mult*100)
MAXITER = np.int64(iter_mult*mult**2*(100)**2/32/5*3)
print MAXITER
PICS = 5
IMAGENO = MAXITER/PICS

maxTime = (40.0/SIZE)**2/10*MAXITER
# raw_input returns the empty string for "enter"
yes = set(['yes','y', 'ye', ''])
no = set(['no','n'])


timescalar=0.00000001
approxtime = timescalar*MAXITER*SIZE*SIZE
choice = None
while approxtime>1 and choice not in yes:
    outstr = "System will compute %s timesteps (up to t=%s) on a %s x %s grid.\n This will take approximately %s minutes. \n Continue (yes/no)? " %(MAXITER, maxTime,SIZE,SIZE,approxtime)
    sys.stdout.write(outstr) 
    choice = raw_input().lower()
    if choice in yes:
        sys.stdout.write("Continuing...")
    elif choice in no:
        sys.exit()
    else:
        sys.stdout.write("Please answer yes or no.")
        
   
sys.stdout.write("Running") 
start = time.time()
run_soliton(SIZE,MAXITER+1,IMAGENO,NAME)
end = time.time()

totaltime = (end-start)/60
outstr = "Simulation finished! Took %s minutes.\n" %(totaltime)
sys.stdout.write(outstr) 


#print 'Saving...'
#file_abc = NAME + "_abc"+ ".npy"

#print 'Plotting...'
#plotter(data_abc, NAME+ "_abc")

#print 'Plots Open'
#os.system('open *.eps')

#print 'Plots Closed'
