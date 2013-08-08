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

import numpy as np

from ..flux import Flux1D

#===============================================================================

class Advection1D (Flux1D):
  """ 1D constant advection.
  """
  _meq = 1  # number of equations in model
  
  #-----------------------------------------------------------------------------
  def __init__ (self, v):
    
    self._v = v  # Constant velocity
  
  #-----------------------------------------------------------------------------
  def f (self, q):
    
    r    = np.empty( q.shape[0], dtype=object )  # vector
    r[0] = q[0] * self._v
    
    return r
  
  #-----------------------------------------------------------------------------
  def J (self, q):
    
    neq  = q.shape[0]
    dims = q[0].shape
    
    r      = np.empty( (neq,neq), dtype=object )  # matrix
    r[0,0] = np.ones (dims) * self._v
    
    return r
  
  #-----------------------------------------------------------------------------
  def R (self, q):
    
    neq  = q.shape[0]
    dims = q[0].shape
    
    r      = np.empty( (neq,neq), dtype=object )  # matrix
    r[0,0] = np.ones (dims)
    
    return r
  
  #-----------------------------------------------------------------------------
  def L (self, q):
    
    neq  = q.shape[0]
    dims = q[0].shape
    
    r    = np.empty( (neq,neq), dtype=object )  # matrix
    r[0,0] = np.ones (dims)
    
    return r
  
  #-----------------------------------------------------------------------------
  def eig (self, q):
    
    neq  = q.shape[0]
    dims = q[0].shape
    
    r    = np.empty( neq, dtype=object )  # vector
    r[0] = np.ones (dims) * self._v
    
    return r
  
  #-----------------------------------------------------------------------------
  def MaxWaveSpeed (self, q):
    """ Maximum wave speed in the range of values for q. """
    
    return abs(self._v)
    
#===============================================================================  
