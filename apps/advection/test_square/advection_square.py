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
  """ Test case definition: 1D constant advection with non-smooth initial 
      conditions (square wave) and periodic BCs.
  """
  # Import modules from library
  from hyperpyws.model_equations.advection  import Advection1D
  from hyperpyws.boundary                   import PeriodicBCs
  from hyperpyws.simulation                 import TestCase
  
  # Constant velocity
  v = 1.0
  
  # Initial conditions
  def q_init (x):
    q0 = 0.0*x
    for i,xi in enumerate(x):
      if( xi > 0.3 and xi < 0.5 ):
        q0[i] = 1.0
    return [ q0 ]
  
  # Exact solution
  def q_exact (x,t):
    return q_init( (x-v*t) % 1.0 )
  
  # Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Advection1D(v)
  test.xlims    = [0.0, 1.0]
  test.tend     =  1.0
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  test.qexact   = q_exact
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = 'Run 1D advection with square-wave ICs and periodic BCs.'

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
