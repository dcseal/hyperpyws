#coding: utf8

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
