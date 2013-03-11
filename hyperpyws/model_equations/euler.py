# coding: utf8

import numpy as np

from ..flux import Flux1D

#===============================================================================

class Euler1D (Flux1D):
  """ 1D Euler equations.
  """
  _meq = 3  # number of equations in model
  
  #-----------------------------------------------------------------------------
  def __init__ (self, gamma):
    
    self._gamma = gamma
  
  #-----------------------------------------------------------------------------
  @property
  def gamma (self):
    """ Heat capacity ratio. """
    
    return self._gamma
  
  #-----------------------------------------------------------------------------
  def f (self, q):
    """ Flux function f(q). """
    
    # Rename conserved quantities
    [rho, mom, eng] = q
    
    # Compute primitive variables
    u1 = mom / rho                          # --> careful with division by zero!
    p  = (self._gamma-1.0) * (eng-0.5*rho*u1**2)
    
    # Compute fluxes
    f = np.empty( 3, dtype=object )
    
    f[0] = rho*u1  # should be momentum, but would pass a reference...
    f[1] = mom*u1 + p
    f[2] = (eng+p)*u1
    
    return f
  
  #-----------------------------------------------------------------------------
  def J (self, q):
    """ Jacobian matrix of flux function: J(q)[i,j] = ∂f[i]/∂q[j]. """
    
    # Rename conserved quantities
    [rho, mom, eng] = q
    
    # Mass-averaged velocity [m/s]
    u1  = mom/rho
        
    # Useful temporary variables
    g   = self._gamma
    gre = g*rho*eng
    r2  = rho**2
    m2  = mom**2
    
    # Data structures for matrix J of numpy arrays
    o = np.zeros( q[0].shape )
    e = np.ones ( q[0].shape )
    J = np.empty( (3,3), dtype=object )
    
    # Jacobian matrix
    J[0,:] = [                      o,                           e,         o ]
    J[1,:] = [      0.5*(g-3.0)*u1**2,                 -u1*(g-3.0), (g-1.0)*e ]
    J[2,:] = [ -u1*(gre-(g-1.)*m2)/r2, .5*(2.*gre-3.*(g-1.)*m2)/r2,      g*u1 ]
    
    return J
  
  #-----------------------------------------------------------------------------
  def eig (self, q):
    """ Compute eigenvalues of Jacobian matrix J. """
    
    # Rename conserved quantities
    [rho, mom, eng] = q
    
    # Mass-averaged velocity [m/s]
    u1 = mom/rho
    
    # Other physical quantities
    p = (self._gamma-1.0) * (eng-0.5*rho*u1**2)  # pressure    [N/m^2] = [J/m^3]
    c = np.sqrt( self._gamma * p/rho )           # speed of sound          [m/s]
    
    # Eigenvalues
    eig = np.empty( 3, dtype=object )
    eig[:] = [ u1-c, u1, u1+c ]
    
    return eig
  
  #-----------------------------------------------------------------------------
  def R (self, q):
    """ Matrix having the right eigenvectors of J as columns. """
    
    # Rename conserved quantities
    [rho, mom, eng] = q
    
    # Mass-averaged velocity [m/s]
    u1 = mom/rho
    
    # Useful temporary variables
    g = self._gamma
    umag2 = u1**2
    
    # Other physical quantities
    p = (g-1.0) * (eng-0.5*rho*umag2)      # pressure          [N/m^2] = [J/m^3]
    c = np.sqrt( g * p / rho )             # speed of sound                [m/s]
    H = (eng+p) / rho                      # total enthalpy per unit mass [J/kg]
    
    # Data structure for matrix R of numpy arrays
    e = np.ones ( q[0].shape )
    R = np.empty( (3,3), dtype=object )
    
    # Right eigenvectors of J (along columns)
    R[0,:] = [      e,          e,       e ]
    R[1,:] = [   u1-c,         u1,    u1+c ]
    R[2,:] = [ H-u1*c,  0.5*umag2,  H+u1*c ]
    
    return R
  
  #-----------------------------------------------------------------------------
  def L (self, q):
    """ Matrix having the left eigenvectors of J as rows; L = inv(R). """
    
    # Rename conserved quantities
    [rho, mom, eng] = q
    
    # Mass-averaged velocity [m/s]
    u1 = mom/rho
    
    # Useful temporary variables
    g = self._gamma
    umag2 = u1**2
    
    # Other physical quantities
    p = (g-1.0) * (eng-0.5*rho*u1**2)      # pressure          [N/m^2] = [J/m^3]
    c = np.sqrt( g * p / rho )             # speed of sound                [m/s]
    M = u1 / c                             # Mach number
    
    # Data structure for matrix L of numpy arrays
    L = np.empty( (3,3), dtype=object )
    
    # Temporary arrays
    t0 = (g-1.)/2.* M**2
    t1 = (g-1.)*M
    t2 = (g-1.)/c**2
    
    # Left eigenvectors of J (along rows)
    L[0,:] = [ 0.5*(t0+M),  -0.5*(t1+1.)/c,  0.5*t2 ]
    L[1,:] = [ 1.0 -t0   ,        t1    /c,     -t2 ]
    L[2,:] = [ 0.5*(t0-M),  -0.5*(t1-1.)/c,  0.5*t2 ]
    
    return L
  
  #-----------------------------------------------------------------------------
  def MaxWaveSpeed (self, q):
    """ Maximum wave speed in the range of values for q. """
    
    eig  = self.eig(q)
    vmin = min(eig[0])
    vmax = max(eig[2])
    
    return max(abs(vmin),abs(vmax)) 
    
  #-----------------------------------------------------------------------------
  @staticmethod
  def SolidWallBC_left (q, mx, mbc):
    """ Left solid wall boundary condition, specific to 1D Euler equations. """
    
    # Impose mirror conditions
    for qi in q:
      qi[:mbc] = qi[2*mbc-1:mbc-1:-1]
    
    # Fix momentum (flip sign of velocity)
    q[1][:mbc] *= -1
  
  #-----------------------------------------------------------------------------
  @staticmethod
  def SolidWallBC_right (q, mx, mbc):
    """ Right solid wall boundary condition, specific to 1D Euler equations. """
    
    # Impose mirror conditions
    r = mbc+mx   # Index of first ghost cell on the right
    for qi in q:
      qi[r:] = qi[r-1:r-mbc-1:-1]

    # Fix momentum (flip sign of velocity)
    q[1][r:] *= -1

#===============================================================================
