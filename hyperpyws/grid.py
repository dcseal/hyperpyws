#coding: utf8

import numpy as np

#===============================================================================

class Grid1D(object):
  """
  Simple object containing information for a 1D uniform grid and storage for
  each component of the solution vector.
  
  Parameters
  ----------
  xlims : list of floats (2)
    Limits of computational domain
  mx : int
    Number of grid cells
  mbc : int
    Number of ghost cells at each boundary
  meq : int
    Number of equations (solving system of meq scalar conservation laws)
    
  """
  def __init__(self, xlims, mx, mbc=3, meq=1):
    
    # Copy input data
    self._xlow  = xlims[0]
    self._xhigh = xlims[1]
    self._mx    = mx
    self._mbc   = mbc
    self._meq   = meq
    
    # Mesh spacing
    self._dx = (xlims[1]-xlims[0])/mx
    
    # Create full domain (with ghost cells at boundaries)
    self._x  = np.linspace( xlims[0]-(mbc-0.5)*self._dx, \
                            xlims[1]+(mbc-0.5)*self._dx, \
                            mx+2*mbc )
    # Interior points
    self._xint = self._x[mbc:(mx+mbc)]
    
    # Solution (1 array for each component of state vector)
    self._q = np.empty( meq, dtype=object )
    for i in range(meq):
      self._q[i] = np.zeros(mx+2*mbc)
  
  #-----------------------------------------------------------------------------
  @property
  def q(self):
    return self._q
  
  @q.setter
  def q(self,qnew):
    assert( len(qnew) == self._meq )
    assert( all([len(qi) == self._mx+2*self._mbc for qi in qnew]) )
    
    for i in range(self._meq):  self._q[i] = qnew[i]
  
  #-----------------------------------------------------------------------------
  @property
  def qint(self):
    a = self._mbc
    b = self._mbc+self._mx
    return [ qi[a:b] for qi in self._q ]
  
#  # Or should we return a numpy array?
#  @property
#  def qint(self):
#    a = self._mbc
#    b = self._mbc+self._mx
#    q_int = np.empty( self._meq, dtype=object )
#    for i,qi in enumerate( self._q ):  q_int[i] = qi[a:b]
#    return q_int
  
#  @q.setter
#  def qint(self,qnew):
#    assert( len(qnew) == self._meq )
#    assert( all([len(qi) == mx for qi in qnew]) )
#    a = self._mbc
#    b = self._mbc+self._mx
#    for qi,qin in zip(self._q,qnew):
#      qi[a:b] = qnew
  
  #-----------------------------------------------------------------------------
  @property
  def x   (self):  return self._x
  
  @property
  def dx  (self):  return self._dx
  
  @property
  def xint(self):  return self._xint
  
  @property
  def mx  (self):  return self._mx
  
  @property
  def mbc (self):  return self._mbc
  
  @property
  def meq (self):  return self._meq
  
#===============================================================================
  