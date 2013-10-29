#==============================================================================#
# This file is part of HYPERPYWS: Hyperbolic Python WENO Solver
#
#   *** This software is made available "as is" without any assurance that it
#   *** will work for your purposes.  The software may in fact have defects, so
#   *** use the software at your own risk.
#
# License: GPL, see COPYING for details
#
# Copyright (C) 2013 
#
#    David Seal,  seal@math.msu.edu,  Michigan State University
#    Yaman Guclu, guclu@math.msu.edu, Michigan State University
#
#===============================================================================

try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ Test case definition: shock-tube (SOD) for 1D Euler's equations with 
      outflow boundary conditions.
  """
  # Import external modules
  import numpy as np
  
  # Import modules from library
  from hyperpyws.model_equations.euler  import Euler1D
  from hyperpyws.boundary               import OutflowBC_left, OutflowBC_right
  from hyperpyws.simulation             import TestCase
  
  # Ratio of specific heats
  gamma = 1.4
  
  # Initial conditions
  def q_init (x):

## left half of Woodward and Collela blast wave problem:
#   rho =        np.ones (x.shape)
#   u1  =        np.zeros(x.shape)
#   p   = 1000.0*np.ones (x.shape)
    
#   for i,xi in enumerate(x):
#     if xi > 0.5:
#       p[i] = 0.01

#   eng  = p/(gamma-1.0) + 0.5*rho*u1**2
#   q    = np.empty( 3, dtype=object )
#   q[:] = [ rho, rho*u1, eng ]

    # Ridiculous "Lax" shock-tube problem (ICs due to Harten):
    rho =  0.445  * np.ones(x.shape)
    m   =  0.3111 * np.ones(x.shape)    # u = 0.6991011235955056
    E   =  8.928  * np.ones(x.shape)    # p = 3.5277019280898867
   
    for i,xi in enumerate(x):
      if xi > 0.5:
        rho[i] = 0.5
        m  [i] = 0.                     # u = 0.0
        E  [i] = 1.4275                 # p = 0.5708

    # note that these give the following 
    q    = np.empty( 3, dtype=object )
    q[:] = [ rho, m, E ]
    
   
    return q
  
  # Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    def SetBCs(q,t):
      OutflowBC_left (q,mx,mbc)
      OutflowBC_right(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Euler1D( gamma )
  test.xlims    = [0.0, 1.0]
  test.tend     = 0.16
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = "Run shock-tube (SOD) for 1D Euler's eqs. with outflow BCs."

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
