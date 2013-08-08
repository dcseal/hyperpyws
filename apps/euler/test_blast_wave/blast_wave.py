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
  """ Test case definition: blast-wave problem for 1D Euler's equations with 
      solid-wall boundary conditions.
  """
  # Import external modules
  import numpy as np
  
  # Import modules from library
  from hyperpyws.model_equations.euler  import Euler1D
  from hyperpyws.simulation             import TestCase
  
  # Ratio of specific heats
  gamma = 1.4
  
  # Initial conditions
  def q_init (x):
    
    rho =     np.ones (x.shape)
    u1  =     np.zeros(x.shape)
    p   = 100*np.ones (x.shape)
    
    for i,xi in enumerate(x):
      if xi < 0.1:
        p[i] = 1000.
      elif 0.1 < xi < 0.9:
        p[i] = 0.01
    
    eng  = p/(gamma-1.0) + 0.5*rho*u1**2
    q    = np.empty( 3, dtype=object )
    q[:] = [ rho, rho*u1, eng ]
    
    return q
  
  # Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    def SetBCs(q,t):
      Euler1D.SolidWallBC_left (q,mx,mbc)  # IMP: These boundary conditions do
      Euler1D.SolidWallBC_right(q,mx,mbc)  # not work properly with TD_RK4, yet!
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Euler1D( gamma )
  test.xlims    = [0.0, 1.0]
  test.tend     =  0.038
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = "Run blast-wave for 1D Euler's eqns. with solid-wall BCs."

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
