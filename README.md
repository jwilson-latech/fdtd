G-FDTD+ABC PACKAGE
==================

This package implements the G-FDTD (M=1) scheme with absorbing boundary conditions. 

    https://github.com/quantumapoptosi/fdtd

DEPENDENCIES
============

    Python
    NumPy
    PyOpenCL
    MatPlotLib

SETUP AND USE
=============

 1. Ensure that Python 2.7 and dependencies are installed.

 2. Change the 'exact' function in main.py to appropriate initial condition.

 3. At bottom of main.py change main(width,steps) to be the desired 
    spacial resolution and the number of time steps required
    Example: main(100,200) sets up an 100x100 grid and computes 
    200 iterations on the grid.
 
 4. In terminal run main.py with 'python main.py' without quotes.

ABOUT THE CODE
==============

The example server and consumer code provided in this package are
intended to be instructional in the use of this OpenID library.  While
it is not recommended to use the example code in production, the code
should be sufficient to explain the general use of the library.

If you aren't familiar with the Django web framework, you can quickly
start looking at the important code by looking in the 'views' modules:

  djopenid.consumer.views
  djopenid.server.views

Each view is a python callable that responds to an HTTP request.
Regardless of whether you use a framework, your application should
look similar to these example applications.

EXAMPLE
=======
In terminal run main.py with 'python main.py' without quotes.



CONTACT
=======
Please feel free to contact me for questions, comments, or concerns. 

      Joshua Wilson (jwilson@latech.edu)
