# coding: utf8

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


import numpy as np
#from numpy import sqrt

from ..flux import Flux1D

#===============================================================================

class ShallowWater1D (Flux1D):
  """ 1D Shallow Water equations.
  """
  _meq = 2  # number of equations in model
  
  #-----------------------------------------------------------------------------
  def __init__ (self, g):
    
    self._g = g
  
  #-----------------------------------------------------------------------------
  @property
  def g(self):
    """ Gravity. """
    
    return self._g
  
  #-----------------------------------------------------------------------------
  def f (self, q):
    """ Flux function f(q). """
   
    # Conservd quantities
    [h, hu] = q

    # Compute primitive variables
    u = hu / h     # --> careful with division by zero!
   
    g = self._g

    # Compute fluxes
    f = np.empty( 2, dtype=object )
    
    f[0] = hu.copy()
    f[1] = h*(u**2) + 0.5*g*(h**2)
    
    return f
  
  #-----------------------------------------------------------------------------
  def J (self, q):
    """ Jacobian matrix of flux function: J(q)[i,j] = ∂f[i]/∂q[j]. """
    
    # Rename conserved quantities
    [h, hu] = q
    
    # Mass-averaged velocity [m/s]
    u = hu/h
        
    # Useful temporary variables
    g   = self._g
    
    # Data structures for matrix J of numpy arrays
    o = np.zeros( q[0].shape )
    e = np.ones ( q[0].shape )
    J = np.empty( (2,2), dtype=object )
    
    # Jacobian matrix
    J[0,:] = [           o,      e ]
    J[1,:] = [ -u**2 + g*h,  2.0*u ]
    
    return J
  
  #-----------------------------------------------------------------------------
  def eig (self, q):
    """ Compute eigenvalues of Jacobian matrix J. """

    # Rename conserved quantities
    [h, hu] = q
    
    # Mass-averaged velocity [m/s]
    u = hu/h
        
    # Useful temporary variables
    sq_gh  = np.sqrt( self._g*h )

    # Eigenvalues
    eig = np.empty( 2, dtype=object )
    eig[:] = [ u - sq_gh, u + sq_gh ]
    
    return eig
  
  #-----------------------------------------------------------------------------
  def R (self, q):
    """ Matrix having the right eigenvectors of J as columns. """

    # Rename conserved quantities
    [h, hu] = q
    
    # Mass-averaged velocity [m/s]
    u = hu/h
        
    # Useful temporary variables
    sq_gh  = np.sqrt( self._g*h )

   
    # Data structure for matrix R of numpy arrays
    e = np.ones ( q[0].shape )
    R = np.empty( (2,2), dtype=object )
    
    # Right eigenvectors of J (along columns)
    R[0,:] = [         e,          e ]
    R[1,:] = [   u-sq_gh,  u + sq_gh ] 
    
    return R
  
  #-----------------------------------------------------------------------------
  def L (self, q):
    """ Matrix having the left eigenvectors of J as rows; L = inv(R). """

    # Rename conserved quantities
    [h, hu] = q
    
    # Mass-averaged velocity [m/s]
    u = hu/h
        
    # Useful temporary variables
    g      = self._g
    sq_gh  = np.sqrt( self._g*h )

    # Data structure for matrix L of numpy arrays
    L = np.empty( (2,2), dtype=object )
    
    # Left eigenvectors of J (along rows)
    L[0,:] = [ 0.5*(sq_gh+u)/sq_gh,   -0.5/sq_gh    ]
    L[1,:] = [ 0.5*(sq_gh-u)/sq_gh,    0.5/sq_gh    ]
    
    return L
  
  #-----------------------------------------------------------------------------
  def MaxWaveSpeed (self, q):
    """ Maximum wave speed in the range of values for q. """
    
    eig  = self.eig(q)
    vmin = min(eig[0])
    vmax = max(eig[1])
    
    return max(abs(vmin),abs(vmax)) 
    
#===============================================================================
