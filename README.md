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

In order to see the results of your simulation, you need to specify the -f
flag indicating the number of frames you desire.

For example, here is an example of a single run of Burger's equation using 
mx=100 points, and a CFL of 0.8:

~$ cd $HYPERPYWS/apps/burgers/test_sine_to_n
~$ pythong test_sine_to_n.py -h
usage: python sine_to_n.py [-h] [-O {5,7}] [-w {JS,Z,CFD}] [-s X] [-f N] [-v V] [-o [FILE]] mx CFL
...
~$ python sin_to_n.py 100 0.8 -f 10
