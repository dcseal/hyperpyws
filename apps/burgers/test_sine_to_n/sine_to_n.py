try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ Test case definition: 1D Burger's equation with periodic boundary
  conditions and smooth initial conditions.  When run to a large enough time,
  this problem develops shocks.
  """
  # Import modules from library
  from hyperpyws.model_equations.burgers    import Burgers1D
  from hyperpyws.boundary                   import PeriodicBCs
  from hyperpyws.simulation  import TestCase
  
  # Initial conditions
  def q_init (x):
    import numpy as np
    return [ 0.5*(1. + np.sin(2.0*np.pi*x) ) ]

  # Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Burgers1D( )
  test.xlims    = [0.0, 1.0]
  test.tend     = 0.4
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = 'Run Buckley-Leverett with square-wave ICs and outflow BCs.'

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
