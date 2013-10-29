#coding: utf8

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


from ..weno import WenoReconstruction

#===============================================================================
# CLASS: Upwinded Central Finite Difference
#===============================================================================

class CentralFiniteDifference5( WenoReconstruction ):
  """
  Abstract class for performing 5th-order linear conservative reconstruction on 
  a uniform mesh.  Using a 5-point stencil centered about x[i], function u(x) is 
  reconstructed at locations x[i-1/2] and x[i+1/2].  3rd-order recontructions on
  small 3-point stencils follow standard WENO5 algorithm, but these are then 
  linearly combined to provide full 5th-order accuracy of 5-point stencil.
    
  """
  _stencil = [-2,-1,0,+1,+2]  # 5-point stencil
  _mbc     = 3                # required number of ghost-cells
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_right (cls, *u_stencil):
    """ Reconstruct u_{i-1/2}.
    """
    # Use reconstruct_left method, after inverting the order of the stencil:
    return cls.reconstruct_left( *u_stencil[::-1] )
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_left (cls, *u_stencil):
    """ Reconstruct u_{i+1/2}.
    """
    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui
    u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1
    u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2
    
    # Return 5th-order conservative linear reconstruction
    return 0.1*u1 + 0.6*u2 + 0.3*u3
  
#===============================================================================

class CentralFiniteDifference7( WenoReconstruction ):
  """
  Abstract class for performing a 7-th order conservative 
  reconstruction based on the linear weights from a WENO method.
    
  """
  _stencil = [-3,-2,-1,0,+1,+2,+3]  # 7-point stencil
  _mbc     = 4                      # required number of ghost-cells
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_right (cls, *u_stencil):
    """ Reconstruct u_{i-1/2}.
    """
    # Use reconstruct_left method, after inverting the order of the stencil:
    return cls.reconstruct_left( *u_stencil[::-1] )
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_left (cls, *u_stencil):
    """ Reconstruct u_{i+1/2}.
    """
    # lazy, but readable indexing into list that's being passed in:
    uim3, uim2, uim1, ui, uip1, uip2, uip3 = u_stencil
    
    # reconstruction using smaller stencils
    u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui
    u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1
    u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2
    u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3
    
    # Get linear weights
    gamma = [1./35., 12./35., 18./35., 4./35.]
    
    # Return 7th-order conservative reconstruction
    return gamma[0]*u1 + gamma[1]*u2 + gamma[2]*u3 + gamma[3]*u4
 
#===============================================================================
