import numpy as np
from numpy import pi, cosh, sin, cos, exp, sqrt

def gaussian(x,y,stdx=1,stdy=1):
    return exp(-(x**2/stdx**2+y**2/stdy**2)/2.0)
    
def planewave(arg):
    return exp(-1j*arg)

class BrightSoliton: 
    
    def __init__(self, nlcoup, gvel, rot, displacement = 0, amp=1, pdir=0:
        self.nlcoup = nlcoup
        self.amp = amp
        self.gvel = gvel
        self.rot = rot
        self.pdir = rot+pdir
        self.spread = amp*sqrt(nlcoup/2)
        self.mass = 2.0*amp*sqrt(2/nlcoup)
        self.momk = gvel/2.0
        self.pvel = gvel/2.0 - amp**2*nlcoup/gvel
        self.displacement = displacement
        
    def __str__(self):
        ''' 
        This makes printing the soliton nice. 
        
        '''
        A0 = self.amp
        k1x = self.spread*cos(self.rot*pi/180)
        k1y = self.spread*sin(self.rot*pi/180)
        k2x = self.momk*cos(self.pdir*pi/180)
        k2y = self.momk*sin(self.pdir*pi/180)
        w1 = self.spread*self.gvel
        w2 = self.momk*self.pvel
        k1z0 = self.spread*self.displacement
        k2z0 = self.momk*self.displacement
        fxn_str = '\n %s*sech( %s*x + %s*y - %s - %s*t )*exp( -i*( %s*x + %s*y - %s - %s*t ) )\n' % (A0,k1x,k1y,k1z0,w1,k2x,k2y,k2z0,w2)
        return fxn_str

    def evaluate(self,x,y,t):
        ''' 
        This is the function for the initial conditions of a soliton
        with the parameters specified in __init__
        
        '''
        A0 = self.amp
        k1x = self.spread*cos(self.rot*pi/180)
        k1y = self.spread*sin(self.rot*pi/180)
        k2x = self.momk*cos(self.pdir*pi/180)
        k2y = self.momk*sin(self.pdir*pi/180)
        w1 = self.spread*self.gvel
        w2 = self.momk*self.pvel
        k1z0 = self.spread*self.displacement
        k2z0 = self.momk*self.displacement
        evaluation = A0*planewave( k2x*x + k2y*y - k2z0 - w2*t)/cosh( k1x*x + k1y*y - k1z0 - w1*t )
        return evaluation
        
class GaussianPacket: 
    
    def __init__(self, pvel=1, amp=1, std=[1,1], displacement = [0,0],pdir=[0,0]):
        self.amp = amp
        self.std = np.array(std)
        self.pvel = pvel
        self.pdir = np.array(pdir)
        self.mass = pi*self.std[0]*self.std[1]*amp**2
        self.displacement = displacement
        
    def __str__(self):
        ''' 
        This makes printing the soliton nice. 
        
        '''
        A0 = self.amp
        stdx = self.std[0]
        stdy = self.std[1]
        k2x = self.pdir[0]
        k2y = self.pdir[1]
        vp = self.pvel
        z0 = self.displacement
        u0 = self.pdir
        x0 = z0[0]
        y0 = z0[1]
        fxn_str = '\n %s*gaussian(x - %s, y - %s, %s,%s)*planewave( %s*(x-%s) + %s*(y-%s) - %s*t )\n' % (A0,x0,y0,stdx,stdy,k2x,x0,k2y,y0,vp)
        return fxn_str

    def evaluate(self,x,y,t):
        ''' 
        This is the function for the initial conditions of a soliton
        with the parameters specified in __init__
        
        '''
        A0 = self.amp
        stdx = self.std[0]
        stdy = self.std[1]
        k2x = self.pdir[0]
        k2y = self.pdir[1]
        vp = self.pvel
        z0 = self.displacement
        u0 = self.pdir
        x0 = z0[0]
        y0 = z0[1]
        evaluation = A0*planewave( k2x*(x-x0) + k2y*(y-y0) - vp*t )*gaussian(x - x0, y - y0, stdx,stdy)
        return evaluation
        

        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        