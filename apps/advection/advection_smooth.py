try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ Test case definition: 1D constant advection with smooth initial conditions
      (sine wave) and periodic BCs.
  """
  # Import external modules
  import numpy as np
  
  # Import modules from library
  from hyperpyws.model_equations.advection  import Advection1D
  from hyperpyws.boundary                   import PeriodicBCs
  from hyperpyws.simulation                 import TestCase
  
  # Constant velocity
  v = 1.0
  
  # Initial conditions
  def q_init (x):
    return [ np.sin( 2.0*np.pi*x ) ]
  
  # Exact solution
  def q_exact (x,t):
    return q_init( x-v*t )
  
  # Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Advection1D(v)
  test.xlims    = [0.0, 1.0]
  test.tend     =  1.0
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  test.qexact   = q_exact
    
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = 'Run 1D advection eqn. with smooth ICs and periodic BCs.'

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
