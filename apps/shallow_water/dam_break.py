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
  """ Test case definition: dam-break for Shallow Water equations with
      outflow boundary conditions.
  """
  # Import external modules
  import numpy as np
  
  # Import modules from library
  from hyperpyws.model_equations.shallow_water import ShallowWater1D
  from hyperpyws.boundary               import OutflowBC_left, OutflowBC_right
  from hyperpyws.simulation             import TestCase
  
  # Specific gravity
  g = 1.0
  
  # Initial conditions
  def q_init (x):
    
    h = 3.0*np.ones (x.shape)
    u =     np.zeros(x.shape)
    
    for i,xi in enumerate(x):
      if xi > 0.5:
        h[i] = 1.0
    
    q    = np.empty( 2, dtype=object )
    q[:] = [ h, h*u ]
    
    return q
  
  # Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    def SetBCs(q,t):
      OutflowBC_left (q,mx,mbc)
      OutflowBC_right(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = ShallowWater1D( g )
  test.xlims    = [0.0, 1.0]
  test.tend     = 0.2
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = "Run dam-break for 1D Shallow Water eqs. with outflow BCs."

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
