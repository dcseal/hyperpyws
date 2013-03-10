#coding: utf8

from ..weno import WenoReconstruction

#===============================================================================
# FACTORY FUNCTION
#===============================================================================

def Weno5 (WenoType='JS'):
  
  if   WenoType == 'JS':  return Weno5_JS
  elif WenoType == 'Z' :  raise  NotImplementedError
  else                 :  raise  ValueError

#===============================================================================
# CLASS: Weno5_JS
#===============================================================================

class Weno5_JS (WenoReconstruction):
  """
  Abstract class for performing 5th-order WENO conservative reconstruction on a 
  uniform mesh.  Smoothness indicators and non-linear weights follow Jing-Shu's 
  algorithm.  Using a 5-point stencil centered about x[i], function u(x) is 
  reconstructed at locations x[i-1/2] and x[i+1/2].
  
  """
  _stencil = [-2,-1,0,+1,+2]  # 5-point stencil
  _eps     = 1.e-12           # regularization parameter (avoids division by 0)
  _mbc     = 3                # required number of ghost-cells
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_right (cls, *u_stencil):
    """ Reconstruct u_{i-1/2}.
    """
    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*3
    beta[0]=(13./12.)*(uim2-2*uim1+ui)**2+0.25*(uim2-4*uim1+3*ui)**2
    beta[1]=(13./12.)*(uim1-2*ui+uip1)**2+0.25*(uim1-uip1)**2
    beta[2]=(13./12.)*(ui-2*uip1+uip2)**2+0.25*(3*ui-4*uip1+uip2)**2
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = (-1./6.)*uim2 + (5./6.)*uim1 + ( 1./3.)*ui
    u2 = ( 1./3.)*uim1 + (5./6.)*ui   - ( 1./6.)*uip1
    u3 = (11./6.)*ui   - (7./6.)*uip1 + ( 1./3.)*uip2
    
    # Get linear weights and regularization parameter
    gamma = [0.3, 0.6, 0.1]
    eps   = cls._eps
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(gamma,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 5th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_left (cls, *u_stencil):
    """ Reconstruct u_{i+1/2}.
    """
    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*3
    beta[0]=(13./12.)*(uim2-2*uim1+ui)**2+0.25*(uim2-4*uim1+3*ui)**2
    beta[1]=(13./12.)*(uim1-2*ui+uip1)**2+0.25*(uim1-uip1)**2
    beta[2]=(13./12.)*(ui-2*uip1+uip2)**2+0.25*(3*ui-4*uip1+uip2)**2
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui
    u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1
    u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2
    
    # Get linear weights and regularization parameter
    gamma = [0.1, 0.6, 0.3]
    eps   = cls._eps
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(gamma,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 5th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3
  
  #-----------------------------------------------------------------------------
  @classmethod
  def get_eps (cls):
    return cls._eps
  
  @classmethod
  def set_eps (cls, eps):
    cls._eps = eps
  
#===============================================================================
