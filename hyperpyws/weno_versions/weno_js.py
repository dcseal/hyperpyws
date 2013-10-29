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
# CLASS: Weno5_JS
#
#   See: G.-S. Jiang and C.-W. Shu, "Efficient Implementation of Weighted ENO 
#        Schemes". J. Comput. Phys. 126 (1996), pp. 202-228. 
#===============================================================================

class Weno5_JS( WenoReconstruction ):
  """
  Abstract class for performing 5th-order WENO conservative reconstruction on a 
  uniform mesh.  Smoothness indicators and non-linear weights follow Jiang-Shu's 
  algorithm.  Using a 5-point stencil centered about x[i], function u(x) is 
  reconstructed at locations x[i-1/2] and x[i+1/2].
  
  """
  _stencil = [-2,-1,0,+1,+2]  # 5-point stencil
  _eps     = 1.e-12           # regularization parameter (avoids division by 0)
  _mbc     = 3                # required number of ghost-cells
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_right( cls, *u_stencil ):
    """ Reconstruct u_{i-1/2}.
    """
    # Use reconstruct_left method, after inverting the order of the stencil:
    return cls.reconstruct_left( *u_stencil[::-1] )
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_left( cls, *u_stencil ):
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
  def get_eps( cls ):
    return cls._eps
  
  @classmethod
  def set_eps( cls, eps ):
    cls._eps = eps
  
#===============================================================================
# CLASS: Weno7_JS
#
#   See: D.S. Balsara and C.-W. Shu, "Monotonicity Preserving Weighted 
#        Essentially Non-oscillatory Schemes with Increasingly High Order of 
#        Accuracy". J. Comput. Phys. 160 (2000), pp. 405-452.
#===============================================================================

class Weno7_JS( WenoReconstruction ):
  """
  Abstract class for performing 7th-order WENO conservative reconstruction on a 
  uniform mesh.  Smoothness indicators and non-linear weights follow Jiang-Shu's 
  algorithm.  Using a 7-point stencil centered about x[i], function u(x) is 
  reconstructed at locations x[i-1/2] and x[i+1/2].
  
  """
  _stencil = [-3,-2,-1,0,+1,+2,+3]  # 7-point stencil
  _eps     = 1.e-12                 # regularization param. (avoids divide by 0)
  _mbc     = 4                      # required number of ghost-cells
  
  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_right( cls, *u_stencil ):
    """ Reconstruct u_{i-1/2}.
    """
    # Use reconstruct_left method, after inverting the order of the stencil:
    return cls.reconstruct_left( *u_stencil[::-1] )

  #-----------------------------------------------------------------------------
  @classmethod
  def reconstruct_left( cls, *u_stencil ):
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
    
    # Get linear weights and regularization parameter
    C     = [1./35., 12./35., 18./35., 4./35.]
    eps   = cls._eps
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(C,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 7th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3 + om[3]*u4
 
  #-----------------------------------------------------------------------------
  @classmethod
  def get_eps( cls ):
    return cls._eps
  
  @classmethod
  def set_eps( cls, eps ):
    cls._eps = eps
  
#===============================================================================
