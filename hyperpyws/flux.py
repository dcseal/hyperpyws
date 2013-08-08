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
