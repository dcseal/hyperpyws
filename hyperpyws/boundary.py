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

#===============================================================================

def PeriodicBCs (q, mx, mbc):
  """ Impose periodic boundary conditions.
  """
  r = mbc+mx   # Index of first ghost cell on the right
  for qi in q:
    qi[:mbc] = qi[r-mbc:  r]  # Populate ghost cells on the left
    qi[r:  ] = qi[mbc:2*mbc]  # Populate ghost cells on the right

#===============================================================================

def OutflowBC_left  (q, mx, mbc):
  """ Impose outflow boundary condition at the left end.
  """
  for qi in q:
    qi[:mbc] = qi[mbc]

def OutflowBC_right (q, mx, mbc):
  """ Impose outflow boundary condition at the right end.
  """
  r = mbc+mx   # Index of first ghost cell on the right
  for qi in q:
    qi[r:] = qi[r-1]

#===============================================================================
