#coding: utf8

from abc import ABCMeta, abstractmethod

#===============================================================================

class WenoMetaClass (ABCMeta):
  """ Metaclass that adds the properties 'stencil' and 'mbc' 
      to an abstract base class.
  """
  @property
  def stencil (cls):
    """ Return symmetric stencil for reconstruction of u[i-1/2].
    """
    return cls._stencil
  
  @property
  def mbc (cls):
    """ Return number of ghost-cells required by reconstruction.
    """
    return cls._mbc
  
#===============================================================================

class WenoReconstruction (object):
  """ Abstract base class for any WENO reconstruction algorithm.
  """
  __metaclass__ = WenoMetaClass
  
  def __init__(self):
    raise Exception('Abstract class cannot be instantiated.')
  
  @abstractmethod
  def reconstruct_left  (cls, *u_stencil):
    """ Reconstruct quantity u[i-1/2] using stencil shifted to the left. """
  
  @abstractmethod
  def reconstruct_right (cls, *u_stencil):
    """ Reconstruct quantity u[i-1/2] using stencil shifted to the right. """

#===============================================================================