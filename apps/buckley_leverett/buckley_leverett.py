try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ Test case definition: 1D Buckley-Leverett equation with non-smooth 
      initial conditions (square wave) and outflow BCs.
  """
  # Import modules from library
  from hyperpyws.model_equations.buckley_leverett  import BuckleyLeverett1D
  from hyperpyws.boundary    import OutflowBC_left, OutflowBC_right
  from hyperpyws.simulation  import TestCase
  
  # Free parameter
  M = 1./3.
  
  # Initial conditions
  def q_init (x):
    q0 = 0.0*x
  #  q0 = 0.0*x + 0.6
    for i,xi in enumerate(x):
      if -0.5 < xi < 0.0:
        q0[i] = 1.0
    return [ q0 ]
  
  # Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    def SetBCs(q,t):
      OutflowBC_left (q,mx,mbc)
      OutflowBC_right(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = BuckleyLeverett1D( M )
  test.xlims    = [-1.0, 1.0]
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
