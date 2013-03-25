#coding: utf8

from ..weno import WenoReconstruction

#===============================================================================
# FACTORY FUNCTION
#===============================================================================

def Weno7 (WenoType='JS'):
  
  if   WenoType == 'JS':  return Weno7_JS
  elif WenoType == 'Z' :  raise ValueError
  else                 :  raise  ValueError

#===============================================================================
# CLASS: Weno7_BS
#
#    See: "Monotonicity Preserving Weighted Essentially Non-oscillatory Schemes 
#           with Increasingly High Order of Accuracy", Balsara and Shu, 
#           JCP, 1999.
#
#===============================================================================

class Weno7_BS (WenoReconstruction):
  """
  Abstract class for performing 5th-order WENO conservative reconstruction on a 
  uniform mesh.  Smoothness indicators and non-linear weights follow Jing-Shu's 
  algorithm.  Using two shifted 5-point stencils centered about x[i], 
  function u(x) is reconstructed at locations x[i-1/2] and x[i+1/2].
  
  """
  _stencil = [-3,-2,-1,0,+1,+2,+3]  # 7-point stencil
  _eps     = 1.e-12                 # regularization parameter (avoids division by 0)
  _mbc     = 5                      # required number of ghost-cells
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_right (cls, *u_stencil):
    """ Reconstruct u_{i-1/2} using an extra point to the right for the stencil.
    """
    # lazy, but readable indexing into list that's being passed in:
    uim3, uim2, uim1, ui, uip1, uip2, uip3 = u_stencil

    # we can use the reconstruct_left method, after reorganizing the stencil:
    return cls.reconstruct_left( uip3, uip2, uip1, ui, uim1, uim2, uim3 )

  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_left (cls, *u_stencil):
    """ Reconstruct u_{i+1/2}.
    """

    # lazy, but readable indexing into list that's being passed in:
    uim3, uim2, uim1, ui, uip1, uip2, uip3 = u_stencil
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*4
    beta[0] = \
     uim3*(  547.*uim3 -  3882.*uim2  + 4642.*uim1 - 1854.*ui) + \
     uim2*( 7043.*uim2 - 17246.*uim1  + 7042.*ui) +              \
     uim1*(11003.*uim1 -  9402.*ui  ) + 2107.*ui**2

    beta[1] = \
     uim2*(  267.*uim2 - 1642.*uim1   + 1602.*ui - 494.*uip1) + \
     uim1*( 2843.*uim1 - 5966.*ui     + 1922.*uip1) +        \
     ui*(   3443.*ui   - 2522.*uip1 ) +  547.*uip1**2

    beta[2] = \
     uim1*( 547.*uim1 - 2522.*ui     + 1922.*uip1 - 494.*uip2) + \
     ui  *(3443.*ui   - 5966.*uip1   + 1602.*uip2 )     + \
     uip1*(2843.*uip1 - 1642.*uip2 ) + 267.*uip2**2

    beta[3] = \
      ui*  ( 2107.*ui   -  9402.*uip1   + 7042.*uip2 - 1854.*uip3 ) + \
      uip1*(11003.*uip1 - 17246.*uip2   + 4642.*uip3 )              + \
      uip2*( 7043.*uip2 -  3882.*uip3 ) + 547.*uip3**2 

    # 3rd-order reconstructions using small 3-point stencils
    u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui
    u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1
    u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2
    u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3
    
    # Get optimal weights and regularization parameter
    C     = [1./35., 12./35., 18./35., 4./35.]
    eps   = cls._eps
    
    # Compute optimal weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(C,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 7th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3 + om[3]*u4
 
  #-----------------------------------------------------------------------------
  @classmethod
  def get_eps (cls):
    return cls._eps
  
  @classmethod
  def set_eps (cls, eps):
    cls._eps = eps
  
#===============================================================================
