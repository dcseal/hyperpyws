import time
import numpy as np

from .flux          import Flux1D
from .weno          import WenoReconstruction
from .grid          import Grid1D
from .mol           import MOL
from .timeline      import TimeManager
from .visualization import RealTimeViz

#===============================================================================
# CLASS: test-case
#===============================================================================

class TestCase (object):
  
  __slots__ = ['ModelEqn','xlims','BCs','qinit','qexact','tend']
  
  def __init__(self):
    
    self.ModelEqn = None
    self.xlims    = None
    self.BCs      = None
    self.qinit    = None
    self.qexact   = None
    self.tend     = None
  
  #-----------------------------------------------------------------------------
  def verify (self):
    """ Check that member attributes are of the proper type. """
    
    # Check if all mandatory attributes were set
    mandatory = set(self.__slots__) - set(['qexact'])
    for m in mandatory:
      if getattr(self,m) is None:
        raise ValueError( "Mandatory member '{:s}' not specified".format(m) )
    
    # Check 'ModelEqn'
    if not isinstance( self.ModelEqn, Flux1D ):
      raise TypeError ('ModelEqn must be a Flux1D object')
    
    # Check 'xlims'
    try             :  iter(self.xlims)
    except TypeError:  raise TypeError ('xlims must be list of 2 floats')
    
    if len(self.xlims) != 2:
      raise ValueError('xlims must be list of 2 floats')
    
    if self.xlims[0] > self.xlims[1]:
      raise ValueError('xlims must be ordered')
    
    # Check 'BCs'
    
    #
    # (not sure yet)
    #
    
    # Check 'qinit'
    if not hasattr( self.qinit, '__call__' ):
      raise TypeError ('qinit must be a callable function')
    
    # Check 'qexact'
    if self.qexact is not None:
      if not hasattr( self.qexact, '__call__' ):
        raise TypeError ('qexact must be a callable function')
    
    # Check 'tend'
    if self.tend <= 0.0:
      raise ValueError ('tend must be a positive real number')
    
#===============================================================================
# CLASS: numerical options and objects
#===============================================================================

class Numerics (object):
  
  __slots__ = ['weno','stepper','CFL','mx']
  
  def __init__(self):
    
    self.weno    = None
    self.stepper = None
    self.CFL     = None
    self.mx      = None
  
  #-----------------------------------------------------------------------------
  def verify (self):
    """ Check that member attributes are of the proper type. """
    
    # Check if all attributes were set
    for m in self.__slots__:
      if getattr(self,m) is None:
        raise ValueError( "Mandatory member '{:s}' not specified".format(m) )
    
    # Check 'weno'
    if not issubclass( self.weno, WenoReconstruction ):
      raise TypeError ('weno must be a subclass of WenoReconstruction')
    
    # Check 'integrator'
    if not hasattr( self.stepper, '__call__'):
      raise TypeError ('stepper must be a callable function')
    
    # Check 'CFL'
    if self.CFL <= 0.0:
      raise ValueError ('CFL must be a positive real number')
    
    # Check 'mx'
    if self.mx <= 0:
      raise ValueError ('mx must be a positive integer number')
    
#===============================================================================
# FUNCTION: run simulation
#===============================================================================
from .time_integrators  import *


def RunSimulation ( test, numr, Tout=[0] ):
  
  # Verify input arguments
  assert(isinstance( test, TestCase ));  test.verify()
  assert(isinstance( numr, Numerics ));  numr.verify()
  
  #-----------------------------------------------------------------------------
  # Create objects
  grid   = Grid1D( test.xlims,
                   numr.mx, 
                   numr.weno.mbc, 
                   test.ModelEqn.meq )
  
  solver = MOL   ( grid,
                   test.ModelEqn,
                   numr.weno,
                   test.BCs (numr.mx, numr.weno.mbc) )
  
  clock  = TimeManager()
  viz    = RealTimeViz( grid, clock, test.qexact )
  
  #-----------------------------------------------------------------------------
  # Set initial conditions
  grid.q = test.qinit( grid.x )
  
  # Initialize real-time plots
  viz.plotq_init()
  time.sleep(0.5)
  
  # Time derivatives of state vector  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  if numr.stepper in [TD_RK4, Taylor2]:
    Fc = solver.TimeDerivatives  # multi-deriv. integrator
  else:
    Fc = lambda q,t : solver.TimeDerivatives(q,t)[0]
  
  #-----------------------------------------------------------------------------
  # Real-time plots counter
  pc = 1
  
  # Advance in time
  stop = False
  while clock.t < test.tend and not stop:
    
    # Compute time-step based on CFL number
    v_max = test.ModelEqn.MaxWaveSpeed( grid.q )
    dt = grid.dx / v_max * numr.CFL
    
    # Reduce last time-step if needed
    if clock.t + dt >= test.tend:
      dt   = test.tend - clock.t
      stop = True
    
    # Print information to terminal
    print ('ts = {:3d};  t = {:.3f};  dt = {:.3e}'\
           .format( clock.ts, clock.t, dt ))
      
    # Advance solution
    grid.q = numr.stepper (Fc, grid.q, clock.t, dt)
    
    # Update time and time-step number
    clock.advance( dt )
    
    # Real-time plots
    if clock.t >= Tout[pc]:
      viz.plotq_renew()
      pc += 1
      time.sleep(0.2)
  
  #-----------------------------------------------------------------------------
  # Last plot
  print ('ts = {:3d};  t = {:.3f};  dt = --'.format( clock.ts, clock.t ))
  viz.plotq_renew()
  time.sleep(0.2)

#===============================================================================
