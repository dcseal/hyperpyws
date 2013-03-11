#coding: utf8

from abc import ABCMeta, abstractmethod

#===============================================================================

class Flux1D(object):
  """
  Abstract base class for any 1D flux: requires implementation of flux function 
  f, its Jacobian matrix J, and eigendecomposition of the latter.
  
  """
  __metaclass__ = ABCMeta
  
  @abstractmethod
  def f(self, q):
    """ Flux function f(q). """
  
  @abstractmethod
  def J(self, q):
    """ Jacobian matrix of flux function: J(q)[i,j] = ∂f[i]/∂q[j]. """
  
  @abstractmethod
  def R(self, q):
    """ Matrix having the right eigenvectors of J as columns. """
  
  @abstractmethod
  def L(self, q):
    """ Matrix having the left eigenvectors of J as rows; L = inv(R). """
  
  @abstractmethod
  def eig(self, q):
    """ Compute eigenvalues of Jacobian matrix J. """
  
  @abstractmethod
  def MaxWaveSpeed (self, q):
    """ Maximum wave speed in the range of values for q. """
  
  @property
  def meq (self):
    """ Return number of equations in model. """
    return self._meq
    
#===============================================================================