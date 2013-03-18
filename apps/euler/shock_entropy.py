import numpy as np

try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ Test case definition: shock-entropy interaction for 1D Euler's equations 
      with outflow boundary conditions.
  """
  # Import modules from library
  from hyperpyws.model_equations.euler  import Euler1D
  from hyperpyws.boundary               import OutflowBC_left, OutflowBC_right
  from hyperpyws.simulation             import TestCase
  
  # Ratio of specific heats
  gamma = 1.4
  
  # Initial conditions
  def q_init (x):
    
    eps   = 0.2
    
    rho = 0*x
    u1  = 0*x
    p   = 0*x
    for i,xi in enumerate(x):
      if xi < -4.0:
        rho[i] = 3.857143
        u1 [i] = 2.629369
        p  [i] = 10.3333
      else:
        rho[i] = 1.0 + eps*np.sin(5.0*xi)
        u1 [i] = 0.
        p  [i] = 1.0
    
    eng  = p/(gamma-1.0) + 0.5*rho*u1**2
    q    = np.empty( 3, dtype=object )
    q[:] = [ rho, rho*u1, eng ]
    
    return q

  # Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    def SetBCs(q,t):
      OutflowBC_left (q,mx,mbc)
      OutflowBC_right(q,mx,mbc)
    return SetBCs

  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Euler1D (gamma)
  test.xlims    = [-5.0, 5.0]
  test.tend     = 1.8
  test.BCs      = CreateBC_func
  test.qinit    = q_init
    
  return test

#===============================================================================
# FUNCTION: Parse input arguments
#===============================================================================

def parse_input():
  
  import argparse, sys
  
  parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description = "Run shock-entropy for 1D Euler's eqs. with outflow BCs",
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
  
  parser.add_argument('-o','--output',
                      nargs   = '?',
                      const   = 'final.dat',
                      default =   None,
                      metavar = 'FILE',
                      help    = 'save final solution to output file' + \
                                ' (default name: final.dat)')
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
  grid = RunSimulation( test_case, num_params, Tout )
  
  # Store final solution to file
  if args.output is not None:
    data = [grid.xint] + list(grid.qint)
    np.savetxt( args.output, np.column_stack(data) )
  
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
