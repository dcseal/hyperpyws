# coding: utf8
"""
This module contains the class 'TimeManager', which manages the simulation time.
TimeManager is like a 'single clock' which every other object synchronizes with.
Of course, time should be advanced in one specific place only.

"""
#
# Author: Yaman Güçlü, January 2013 - Michigan State University
#
# Last revision: 27 Jan 2013
#

__all__ = ['TimeManager']
__docformat__ = 'reStructuredText'

#===============================================================================
# CLASS: Time Manager
#===============================================================================

class TimeManager (object):
  """
  A single clock that gives the simulation time.  One can start the clock from a
  specific time 't0', and can advance it by a specified amount 'dt'.  The clock
  provides the time and the time-step number (which always starts from zero).
  
  Parameters
  ----------
  t0 : float
    Initial time instant.
  
  """
  def __init__(self,t0=0.0):
    self._t  = t0
    self._ts = 0
  
  #-----------------------------------------------------------------------------
  def advance (self,dt):
    """
    Advance time by specified amount, and increase time-step number by 1.
    
    Parameters
    ----------
    dt : float
      Time-step size.
    
    """
    self._t  += dt
    self._ts += 1
  
  #-----------------------------------------------------------------------------
  @property
  def t (self):
    """ Property: get time instant.
    """
    return self._t
  
  #-----------------------------------------------------------------------------
  @property
  def ts(self):
    """ Property: get time-step number.
    """
    return self._ts
