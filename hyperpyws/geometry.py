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
# CLASS: 2D Point
#===============================================================================

class Point( object ):
  
  __slots__ = ['x', 'y']
  
  def __init__( self, x, y ):
    self.x = x
    self.y = y
  
  def __repr__( self ):
    return 'Point( {}, {} )'.format(self.x, self.y)

#===============================================================================
# BASE CLASS: 2D segment
#===============================================================================

class Segment( object ):
  
  __metaclass__ = ABCMeta
  
  @property
  @abstractmethod
  def x( self ):
    pass
  
  @property
  @abstractmethod
  def y( self ):
    pass

#===============================================================================
# FUNCTION: Concatenate segments
#===============================================================================

def Concatenate( *segments ):
  x = []
  y = []
  for s in segments:
    if hasattr(s.x,'__iter__'):
      x.extend( s.x )
      y.extend( s.y )
    else:
      x.append( s.x )
      y.append( s.y )
  
  EPS = 1.e-14
  xn = [float(x[0])]
  yn = [float(y[0])]
  for i in range(1,len(x)):
    if (x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 > EPS**2:
      xn.append( float(x[i]) )
      yn.append( float(y[i]) )
  
  return (xn,yn)

#===============================================================================
# CLASS: Line segment
#===============================================================================

class LineSegment( Segment ):
  
  __slots__ = ['A','B']
  
  def __init__( self, A, B):
    self.A = A
    self.B = B
  
  @property
  def x( self ):
    return [self.A.x, self.B.x]
  
  @property
  def y( self ):
    return [self.A.y, self.B.y]
  
  def __repr__( self ):
    return 'LineSegment([ {}, {} ])'.format(self.A, self.B)

#===============================================================================
# CLASS: Curve segment
#===============================================================================

class CurveSegment( Segment ):
  
  __slots__ = ['_x','_y']
  
  def __init__( self, x, y ):
    self._x = x
    self._y = y
  
  @property
  def x( self ):  return self._x
  
  @property
  def y( self ):  return self._y

#===============================================================================
