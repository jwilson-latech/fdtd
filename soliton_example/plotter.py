import numpy as np
from numpy import pi, sqrt, cosh, exp, cos, sin, absolute
import matplotlib.pyplot as plt
from fullviridis import full_viridis

import string
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import cmath
from matplotlib import mlab, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import AxesGrid
import os
from matplotlib.ticker import LogFormatter 


def Timestep(i):
    return i*0.6

cmap = full_viridis
golden_ratio = 1.61818
ratio = 1
def myplots(fig,data):
    """
    A grid of 2x2 images with a single colorbar
    """
    grid = AxesGrid(fig, 111,  # similar to subplot(142)
                    nrows_ncols=(1, 6),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="none",
                    )

     
    for i in range(6):
        Z = absolute(data[i])**2
        grid[i].set_title(r"$t={:03.1f}$".format(Timestep(i)),color='black',horizontalalignment='center',verticalalignment='bottom')
        im = grid[i].imshow(Z, extent=(-2, 2, -2, 2), interpolation="gaussian",origin="lower",cmap=cmap,norm=LogNorm(vmin=1e-5,vmax=1))
        grid[i].set_aspect(ratio)
        grid[i].set_xlabel("$x/10$",size=16)
    #plt.colorbar(im, cax = grid.cbar_axes[0])
    #ticks = np.logspace(1e-6,1,7)
    #lf = LogFormatter(10, labelOnlyBase=False)
    grid[0].set_ylabel("$y/10$",size=16)
    pos2 = [0.905,0.25,0.01,0.5]
    position = fig.add_axes(pos2)
    ticks=np.logspace(1e-6,1e-1,6)
    fig.colorbar(im, ax=grid[5],cax=position,extend="both")
    
    for cax in grid.cbar_axes:
        cax.toggle_label(True)
 
    # This affects all axes as share_all = True.
    grid.axes_llc.set_xticks([-2,-1, 0,1])
    #grid[0].set_xticks([-20,-10, 0,10, 20])
    grid.axes_llc.set_yticks([-2, -1, 0, 1,2])
    
def debugmyplots(fig,data):
    """
    A grid of 2x2 images with a single colorbar
    """
    grid = AxesGrid(fig, 111,  # similar to subplot(142)
                    nrows_ncols=(1, 6),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="none",
                    )

     
    for i in range(6):
        Z = np.real(data[i])
        grid[i].set_title(r"$t={:03.1f}$".format(Timestep(i)),color='black',horizontalalignment='center',verticalalignment='bottom')
        im = grid[i].imshow(Z, extent=(-2, 2, -2, 2), interpolation="Gaussian",origin="lower",cmap='seismic',vmin=-1,vmax=1)
        grid[i].set_aspect(ratio)
        grid[i].set_xlabel("$x/10$",size=16)
    #plt.colorbar(im, cax = grid.cbar_axes[0])
    #ticks = np.logspace(1e-6,1,7)
    #lf = LogFormatter(10, labelOnlyBase=False)
    grid[0].set_ylabel("$y/10$",size=16)
    pos2 = [0.905,0.25,0.01,0.5]
    position = fig.add_axes(pos2)
    ticks=np.logspace(1e-6,1e-1,6)
    fig.colorbar(im, ax=grid[5],cax=position,extend="both")
    
    for cax in grid.cbar_axes:
        cax.toggle_label(True)
 
    # This affects all axes as share_all = True.
    grid.axes_llc.set_xticks([-2,-1, 0,1])
    #grid[0].set_xticks([-20,-10, 0,10, 20])
    grid.axes_llc.set_yticks([-2, -1, 0, 1,2])
    
def debugplots(fig,data):
    """
    A grid of 2x2 images with a single colorbar
    """
    grid = AxesGrid(fig, 111,  # similar to subplot(142)
                    nrows_ncols=(1, 6),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="none",
                    )

    Z0=data[0].real
    Z1=data[1].real
    Z2=data[2].real
    Z3=data[3].real
    Z4=data[4].real
    Z5=data[5].real
    
    Z=[Z0,Z1,Z2,Z3,Z4,Z5]
    
    for i in range(6):
        grid[i].set_title(r"$t=%u\Delta t$"%(Timestep(i)),color='black',horizontalalignment='center',verticalalignment='bottom')
        im = grid[i].imshow(Z[i], extent=(-2, 2, -2, 2), interpolation="Nearest",origin="lower",cmap='seismic',vmin=-1,vmax=1)
        grid[i].set_aspect(ratio)
        grid[i].set_xlabel("$x/10$",size=16)
    grid[0].set_ylabel("$y/10$",size=16)
    pos2 = [0.905,0.25,0.01,0.5]
    position = fig.add_axes(pos2)
    fig.colorbar(im, ax=grid[2],cax=position,extend="both")
    
    for cax in grid.cbar_axes:
        cax.toggle_label(True)
 
    # This affects all axes as share_all = True.
    grid.axes_llc.set_xticks([-2,-1, 0,1])
    #grid[0].set_xticks([-20,-10, 0,10, 20])
    grid.axes_llc.set_yticks([-2, -1, 0, 1,2])
    
def plotter(data, fname):
    F = plt.figure(1, (10, 4))
    #F.subplots_adjust(left=0.05, right=0.95)
    debugplots(F,data)
    #myplots(F,data)
    plt.savefig(fname+'.eps')
    
#data = np.load('soln2dlbc.npy')  
#plotter(data,'soln2dlbc')


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    