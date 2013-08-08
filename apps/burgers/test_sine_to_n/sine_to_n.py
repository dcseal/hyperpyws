#==============================================================================#
# This file is part of HYPERPYWS: Hyperbolic Python WENO Solver
#
#
#   This software is made available for research and instructional use only.
#   You may copy and use this software without charge for these non-commercial
#   purposes, provided that the copyright notice and associated text is
#   reproduced on all copies.  For all other uses (including distribution of
#   modified versions), please contact the author at the address given below.
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
  """ Test case definition: 1D Burger's equation with periodic boundary
  conditions and smooth initial conditions.  When run to a large enough time,
  this problem develops shocks.
  """
  # Import modules from library
  from hyperpyws.model_equations.burgers    import Burgers1D
  from hyperpyws.boundary                   import PeriodicBCs
  from hyperpyws.simulation  import TestCase
  
  # Initial conditions
  def q_init (x):
    import numpy as np
    return [ 0.5*(1. + np.sin(2.0*np.pi*x) ) ]

  # Spatial derivative of initial conditions 
  # ( this derivative is needed to compute exact solution )
  def q_init_partial_x (x):
    import numpy as np
    return [ np.pi*np.cos(2.0*np.pi*x) ]


  # Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Burgers1D( )
  test.xlims    = [0.0, 1.0]
  test.tend     = 0.25
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = 'Run Buckley-Leverett with square-wave ICs and outflow BCs.'

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
