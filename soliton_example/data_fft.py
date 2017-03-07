import numpy as np
import os
import matplotlib.pyplot as plt
from radial_profile import radial_profile

NAME = 'soln2d' # what to name the files
#print 'Saving...'
file_abc = NAME + "abc"+ ".npy"

data_abc = np.load(file_abc)


a = np.fft.fft2(data_abc)
A = np.absolute(a)
print a.shape

plt.imshow(A[0])
plt.show()
    
center = (50, 50)
radi = 50

rad = radial_profile(A[5], center)

plt.plot(rad)
plt.show()




