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
      
      fig  = plt.figure()
      ax   = fig.add_subplot(1,1,1)
      line,= ax.plot( xint, qint[m], '.' )
      
      ax.grid()
      ax.set_xlabel('x')
      ax.set_ylabel(self._labels[m],rotation='horizontal')
      
      self._plots[m] = {'fig': fig, 'ax': ax, 'line': line}
    
    # Plot exact solution, if available, and store lines
    if self._exact is not None:
      qex = self._exact( xint, self._clock.t )
      
      print "File 'visualization.py':"
      print len(qex)
      print len(qex[0])
      
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
      p[ 'fig'].canvas.draw()

#===============================================================================
