try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ Test case definition: 1D Euler's equations with smooth initial conditions 
      (sine wave) and periodic boundary conditions.
  """
  # Import external modules
  import numpy as np
  
  # Import modules from library
  from hyperpyws.model_equations.euler  import Euler1D
  from hyperpyws.boundary               import PeriodicBCs
  from hyperpyws.simulation             import TestCase
  
  # Ratio of specific heats
  gamma = 1.4
  
  # Initial conditions
  def q_init (x):
    
    rho = 1.0 + 0.2*np.sin(np.pi*x)
    u1  = np.ones(x.shape)
    p   = np.ones(x.shape)
    
    eng  = p/(gamma-1.0) + 0.5*rho*u1**2
    q    = np.empty( 3, dtype=object )
    q[:] = [ rho, rho*u1, eng ]
    
    return q
  
  # Exact solution
  def q_exact (x,t):
    return q_init( x-t )
  
  # Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
    return SetBCs

  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Euler1D( gamma )
  test.xlims    = [0.0, 2.0]
  test.tend     =  2.0
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  test.qexact   = q_exact
  
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = "Run 1D Euler's eqns. with smooth ICs and periodic BCs."

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
