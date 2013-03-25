#coding: utf8

from ..weno import WenoReconstruction

#===============================================================================
# CLASS: Weno5_Z
#
#   See: R. Borges, M. Carmona, B. Costa, and W.S. Don, "An improved weighted 
#        essentially non-oscillatory scheme for hyperbolic conservation laws".
#        J. Comput. Phys. 227 (2008), pp. 3191-3211.  
#        {http://dx.doi.org/10.1016/j.jcp.2007.11.038}.
#===============================================================================

class Weno5_Z( WenoReconstruction ):
  """
  Abstract class for performing 5th-order WENO conservative reconstruction on a 
  uniform mesh.  Smoothness indicators and non-linear weights follow an
  improvement on Jiang-Shu's algorithm.  Using a 5-point stencil centered about 
  x[i], function u(x) is reconstructed at locations x[i-1/2] and x[i+1/2].
  
  """
  _stencil = [-2,-1,0,+1,+2]  # 5-point stencil
  _eps     = 1.e-12           # regularization parameter (avoids division by 0)
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
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*3
    beta[0]=(13./12.)*(uim2-2*uim1+ui)**2+0.25*(uim2-4*uim1+3*ui)**2
    beta[1]=(13./12.)*(uim1-2*ui+uip1)**2+0.25*(uim1-uip1)**2
    beta[2]=(13./12.)*(ui-2*uip1+uip2)**2+0.25*(3*ui-4*uip1+uip2)**2

    # new term not used in JS reconstruction:
    tau5 = abs( beta[0] - beta[2] )

    # (note: these can be wrapped into definition of \tilde{omega}, in order to
    # avoid this extra computation.)
#   beta = [None]*3
#   beta = [ (b+eps) / ( b + tau5 + eps ) for b in beta_JS ]

    # 3rd-order reconstructions using small 3-point stencils
    u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui
    u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1
    u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2
    
    # Get linear weights and regularization parameter
    gamma = [0.1, 0.6, 0.3]
    eps   = cls._eps
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g*(1. + tau5/(b+eps) ) for g,b in zip(gamma,beta) ]  # << NEW for Z
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
