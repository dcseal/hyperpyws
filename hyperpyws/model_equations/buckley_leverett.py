# coding: utf8

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

class BuckleyLeverett1D (Flux1D):
  """ 1D Buckley-Leverett two-phase flow through porous media equation.
  """
  _meq = 1  # number of equations in model
  
  #-----------------------------------------------------------------------------
  def __init__ (self, M):
    
    self._M = M
  
  #-----------------------------------------------------------------------------
  @property
  def M (self):
    """ Free parameter in model. """
    
    return self._M
  
  #-----------------------------------------------------------------------------
  def f (self, q):
    """ Flux function f(q). """
    
    u0   = q[0]**2
    u1   = (1.0-q[0])**2
    
    f    = np.empty( 1, dtype=object )  # vector
    f[0] = u0 / (u0 + self._M*u1)
    
    return f
  
  #-----------------------------------------------------------------------------
  def J (self, q):
    """ Jacobian matrix of flux function: J(q)[i,j] = ∂f[i]/∂q[j]. """
    
    u = q[0]
    M = self._M
    
    J      = np.empty( (1,1), dtype=object )  # matrix
    J[0,0] = (2.*M*u*(1.-u)) / (u**2 + M*(1.-u)**2)**2
    
    return J
  
  #-----------------------------------------------------------------------------
  def eig (self, q):
    """ Compute eigenvalues of Jacobian matrix J. """
    
    eig    = np.empty( 1, dtype=object )
    eig[0] = self.J(q)[0,0]
    
    return eig
  
  #-----------------------------------------------------------------------------
  def R (self, q):
    """ Matrix having the right eigenvectors of J as columns. """
    
    R      = np.empty( (1,1), dtype=object )  # matrix
    R[0,0] = np.ones ( q[0].shape )
    
    return R
  
  #-----------------------------------------------------------------------------
  def L (self, q):
    """ Matrix having the left eigenvectors of J as rows; L = inv(R). """
    
    L      = np.empty( (1,1), dtype=object )  # matrix
    L[0,0] = np.ones ( q[0].shape )
    
    return L
  
  #-----------------------------------------------------------------------------
  def MaxWaveSpeed (self, q):
    """ Maximum wave speed in the range of values for q. """
    
    u = np.linspace( min(q[0]), max(q[0]), 100 )
    M = self._M
    eig = (2.*M*u*(1.-u)) / (u**2 + M*(1.-u)**2)**2
    
    return max(abs(eig))
  
#===============================================================================  
