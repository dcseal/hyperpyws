#coding: utf8

from ..weno import WenoReconstruction

#===============================================================================
# CLASS: Weno5_Linear
#===============================================================================

class Weno5_Linear (WenoReconstruction):
  """
  Abstract class for performing 5th-order WENO conservative reconstruction on a 
  uniform mesh.  Smoothness indicators and non-linear weights follow Jing-Shu's 
  algorithm.  Using a 5-point stencil centered about x[i], function u(x) is 
  reconstructed at locations x[i-1/2] and x[i+1/2].
  
  """
  _stencil = [-2,-1,0,+1,+2]  # 5-point stencil
  _mbc     = 3                # required number of ghost-cells
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_right (cls, *u_stencil):
    """ Reconstruct u_{i-1/2}.
    """
    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = (-1./6.)*uim2 + (5./6.)*uim1 + ( 1./3.)*ui
    u2 = ( 1./3.)*uim1 + (5./6.)*ui   - ( 1./6.)*uip1
    u3 = (11./6.)*ui   - (7./6.)*uip1 + ( 1./3.)*uip2
   
    # Return 5th-order conservative linear reconstruction
    return 0.3*u1 + 0.6*u2 + 0.1*u3
  
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
