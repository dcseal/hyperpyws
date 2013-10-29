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
  """ Test case definition: shock-entropy interaction for 1D Euler's equations 
      with outflow boundary conditions.
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
    
    eps   = 0.2
    
    rho = 0*x
    u1  = 0*x
    p   = 0*x
    for i,xi in enumerate(x):
      if xi < -4.0:
        rho[i] = 3.857143
        u1 [i] = 2.629369
        p  [i] = 10.3333
      else:
        rho[i] = 1.0 + eps*np.sin(5.0*xi)
        u1 [i] = 0.
        p  [i] = 1.0
    
    eng  = p/(gamma-1.0) + 0.5*rho*u1**2
    q    = np.empty( 3, dtype=object )
    q[:] = [ rho, rho*u1, eng ]
    
    return q
  
  # Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    def SetBCs(q,t):
      OutflowBC_left (q,mx,mbc)
      OutflowBC_right(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Euler1D (gamma)
  test.xlims    = [-5.0, 5.0]
  test.tend     = 1.8
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = "Run shock-entropy for 1D Euler's eqs. with outflow BCs."

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
