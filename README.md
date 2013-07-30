===============================================================================
HYPERPYWS: Hyperbolic Python WENO Solver
===============================================================================

Finite difference WENO solver for hyperbolic problems.  Currently supports
hyperbolic problems in 1D.

===============================================================================
Required libraries/tools/software
===============================================================================

. Working installation of Python including Numpy.

===============================================================================
Installation procedure
===============================================================================

[0] Download the tarball and unzip into a location of your choosing.

This software doesn't require modifying Python paths or environment variables.

===============================================================================
Running a simulation
===============================================================================

[1] Open a terminal and navigate to an application located in the apps folder.

[2] Run the application.  The -h flag provides a list of options that are
available.

For example, here is an example of a run inside the Python interpreter:

#>> cd /path/to/hyperpyws/apps/burgers/test_sine_to_n
#>> run  sin_to_n.py -h
#usage: python sine_to_n.py [-h] [-O {5,7}] [-w {JS,Z,CFD}] [-s X] [-f N] [-v V] [-o [FILE]] mx CFL
#>> run sin_to_n.py 100 0.4

