""" Parser and Main functions for running applications scripts interactively, 
    from the command line.  The user must provide a function returning the 
    TestCase object, and a help-message string. 
"""
from __future__ import print_function

#===============================================================================
# FUNCTION: Parse input arguments
#===============================================================================

def parse_input( help_message ):
  
  import argparse, sys
  
  parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description =    help_message,
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
                      choices = range(9),
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
  6. Fehlberg5 (from the 4-5 pair)
  7. Taylor2 (Taylor's method, 2nd-order)
  8. TD-RK4  (Two-derivative Runge-Kutta)
(default: 5)''')
  
  parser.add_argument('-f','--frames',
                      type    = int,
                      default = None,
                      metavar = 'N',
                      help    = 'produce N frames as real-time visualization'+\
                                ' (default: None)')
  
  parser.add_argument('-o','--output',
                      nargs   = '?',
                      const   = 'final.dat',
                      default =  None,
                      metavar = 'FILE',
                      help    = 'save final solution to output file' + \
                                ' (default name: final.dat)')
  return parser.parse_args()

#===============================================================================
# FUNCTION: Main script
#===============================================================================

def main( help_message, DefineTestCase ):
  
  # Parse input arguments
  args = parse_input( help_message )
  print(args)
  print('')
  
  # Test-Case container
  test_case = DefineTestCase()
  
  # Import external modules
  import numpy as np
  
  # Import modules from library
  try               :  import hyperpyws
  except ImportError:  import hyperpyws_path
  
  # Extract time-integrator function
  import hyperpyws.time_integrators as integrators
  stepper_name = integrators.__all__[args.stepper]    # --> dangerous indexing
  stepper_func = getattr( integrators, stepper_name ) 
  
  # Extract numerical parameters from input
  from hyperpyws.weno_versions.weno5  import Weno5_Z
  from hyperpyws.simulation           import Numerics, RunSimulation
  
  num_params         = Numerics()     # Container for numerical parameters
  num_params.weno    = Weno5_Z        # << maybe we want to use Weno5(method) here, ...
  num_params.stepper = stepper_func   # --> should we pass a string?
  num_params.CFL     = args.CFL       # CFL parameter
  num_params.mx      = args.mx        # Number of mesh cells in domain
  
  # Real-time visualization: time instants for creating an output
  if args.frames is not None:
    Tout = np.linspace( 0.0, test_case.tend, args.frames+1 )
  else:
    Tout = []
  
  # Run simulation: call default library function
  grid = RunSimulation( test_case, num_params, Tout )
  
  # Print numerical parameters to file
  if args.output is not None:
    with open( 'num_params.dat', 'w' ) as f:
      print( num_params, file=f )
  
  # Store final solution to file
  if args.output is not None:
    time = test_case.tend
    data = np.column_stack( [grid.xint] + list(grid.qint) )
    fmt  = '%.15e'
    with open( args.output, 'wb' ) as f:
      print( fmt % time, file=f )         # time instant on first row
      np.savetxt( f, data, fmt=fmt )      # data arrays along columns
  
  # Keep matplotlib windows open if necessary
  try: __IPYTHON__
  except NameError:
    import sys
    if not sys.flags.interactive:
      from matplotlib.pyplot import show
      show()

#===============================================================================
