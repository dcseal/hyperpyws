import numpy as np

try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ Test case definition: 1D constant advection with smooth initial conditions
      (sine wave) and periodic BCs.
  """
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
# FUNCTION: Parse input arguments
#===============================================================================

def parse_input():
  
  import argparse, sys
  
  parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description = 'Run 1D advection eqn. with smooth ICs and periodic BCs.',
      formatter_class = argparse.RawTextHelpFormatter,
      )
  
  parser.add_argument('mx',
                      type = int,
                      help = 'number of subdivisions along x axis')
  
  parser.add_argument('CFL',
                      type = float,
                      help = 'maximum Courant number in domain')
  
  parser.add_argument('-s','--time_integrator',
                      type    = int,
                      choices = range(8),
                      default =  5,
                      dest    = 'stepper',
                      metavar = 'X',
                      help    = 
  '''choose integrator X for time-stepping:
  0. FE      (forward-Euler)
  1. RK2     (midpoint rule)
  2. RK2-SSP (Heun's method)
  3. RK3
  4. RK3-SSP
  5. RK4     (gold standard)
  6. Taylor2 (Taylor's method, 2nd-order)
  7. TD-RK4  (Two-derivative Runge-Kutta)
(default: 5)''')
  
  parser.add_argument('-f','--frames',
                      type    = int,
                      default = None,
                      metavar = 'N',
                      help    = 'produce N frames as real-time visualization')
  
  return parser.parse_args()

#===============================================================================
# FUNCTION: Main script
#===============================================================================

def main():
  
  # Parse input arguments
  args = parse_input()
  print(args)
  print('')
  
  # Extract time-integrator function
  import hyperpyws.time_integrators as integrators
  stepper_name = integrators.__all__[args.stepper]
  stepper_func = getattr( integrators, stepper_name )
  
  # Extract numerical parameters from input
  from hyperpyws.weno_versions.weno5  import Weno5_JS
  from hyperpyws.simulation           import Numerics, RunSimulation
  
  # Test-Case container
  test_case = DefineTestCase()
  
  # Numerics container
  num_params         = Numerics()
  num_params.weno    = Weno5_JS       # --> should we pass a string?
  num_params.stepper = stepper_func   # --> should we pass a string?
  num_params.CFL     = args.CFL       # CFL parameter
  num_params.mx      = args.mx        # number of mesh cells in domain
  
  # Real-time visualization: time instants for creating an output
  if args.frames is not None:
    Tout = np.linspace( 0.0, test_case.tend, args.frames+1 )
  else:
    Tout = []
  
  # Run simulation: call default library function
  RunSimulation( test_case, num_params, Tout )

  # Keep matplotlib windows open if necessary
  try: __IPYTHON__
  except NameError:
    import sys
    if not sys.flags.interactive:
      from matplotlib.pyplot import show
      show()

#===============================================================================
if __name__ == '__main__':
  #Run as main program
  main()
