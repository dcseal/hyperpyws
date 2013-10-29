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

import matplotlib.pyplot as plt

#===============================================================================

class RealTimeViz (object):
  
  def __init__(self, grid, clock, exact=None, labels=None):
    
    self._grid  = grid
    self._clock = clock
    self._exact = exact
    self._labels= labels
    self._plots = [None]*grid.meq
    
    if labels is None:
      self._labels = [ 'q{:d}'.format(i) for i in range(self._grid.meq) ]
      
  #-----------------------------------------------------------------------------
  def plotq_init (self):
    
    meq  = self._grid.meq
    xint = self._grid.xint
    qint = self._grid.qint
    
    # Create figures, plot numerical solution, and store dictionary
    for m in range(meq):
      
      data = qint[m]
      fig  = plt.figure()
      ax   = fig.add_subplot(1,1,1)
      line,= ax.plot( xint, data, '.' )
      
      ax.grid()
      ax.set_xlabel('x')
      ax.set_ylabel(self._labels[m],rotation='horizontal')
      
      self._reset_ylim( ax, data )
      self._plots[m] = {'fig': fig, 'ax': ax, 'line': line}
    
    # Plot exact solution, if available, and store lines
    if self._exact is not None:
      qex = self._exact( xint, self._clock.t )
      for m in range(meq):
        ax = self._plots[m]['ax']
        self._plots[m]['line_ex'], = ax.plot( xint, qex[m], '-r' )
        ax.legend(['numerical','exact'], loc='upper right')  
    
    # Give a title to all figures, and show them
    title = 't = %.2f' % self._clock.t
    for p in self._plots:
      p[ 'ax'].set_title(title)
      p['fig'].show()
  
  #-----------------------------------------------------------------------------
  def plotq_renew (self):
    
    # Update exact solution, if available
    if self._exact is not None:
      qex = self._exact( self._grid.xint, self._clock.t)
      for p,qi in zip( self._plots, qex):
        p['line_ex'].set_ydata( qi )
    
    # Update numerical solution and title, and redraw figures
    title = 't = %.2f' % self._clock.t
    for p,qi in zip( self._plots, self._grid.qint ):
      p['line'].set_ydata( qi )
      p[  'ax'].set_title(title)
      
      self._reset_ylim( p['ax'], qi )
      
      p[ 'fig'].canvas.draw()
      
  #-----------------------------------------------------------------------------
  @staticmethod
  def _reset_ylim( ax, data ):
    """ Reset limits of y axis if data does not fit in figure 
        or if margins are too wide.
    """
    # Extract y limits from axes and data, and compute theoretical margin
    ylim   = ax.get_ylim()
    ya,yb  = min(data), max(data)
    margin = (yb-ya)*0.05
    
    # If the data limits are identical, do not do anything
    if ya == yb:
      return
    
    # Check if new y limits are needed
    if   ya <= ylim[0] or ya-ylim[0] > 2*margin:  pass
    elif ylim[1] <= yb or ylim[1]-yb > 2*margin:  pass
    else                                       :  return
    
    # Apply margins and set new y limits
    ya -= margin
    yb += margin
    ax.set_ylim( [ya,yb] )

#===============================================================================
