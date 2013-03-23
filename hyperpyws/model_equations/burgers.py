# coding: utf8

import numpy as np

from ..flux import Flux1D

#===============================================================================

class Burgers1D (Flux1D):
  """ 1D Burger's equation.
  """
  _meq = 1  # number of equations in model
  
  #-----------------------------------------------------------------------------
  def __init__ (self):
    pass   
  
  #-----------------------------------------------------------------------------
  def f (self, q):
    """ Flux function f(q). """
    
    f    = np.empty( 1, dtype=object )  # vector
    f[0] = 0.5 * q[0]**2
    
    return f
  
  #-----------------------------------------------------------------------------
  def J (self, q):
    """ Jacobian matrix of flux function: J(q)[i,j] = ∂f[i]/∂q[j]. """
    
    J      = np.empty( (1,1), dtype=object )  # matrix
    J[0,0] = q[0].copy()   # TODO: better way to do this?
    
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
    return max( max(abs(self.eig(q))) )
  
#===============================================================================  
