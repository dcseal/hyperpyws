# coding: utf8

from __future__ import print_function

#===============================================================================
# FUNCTION: Parse input arguments
#===============================================================================

def parse_input():
  
  import argparse, sys
  
  parser = argparse.ArgumentParser (
      prog='python '+sys.argv[0],
      description='Run 1D advection simulation with increasing number of cells',
      formatter_class=argparse.RawTextHelpFormatter
      )
  
  parser.add_argument('CFL',
                      type = float,
                      help = 'maximum Courant number in domain')
  
  parser.add_argument('-O','--weno_order',
                      type    = int,
                      choices = [5,7],
                      default =  5,
                      help    = 'order of accuracy for WENO recontruction'+\
                                ' (default: 5)')
  
  parser.add_argument('-v','--weno_version',
                      choices = ['JS','Z','CFD'],
                      default =  'JS',
                      help    = 
  '''choose WENO version:
  JS  = WENO-JS (Jiang-Shu's algorithm)
  Z   = WENO-Z  (Borges-Carmona-Costa-Don's algorithm)
  CFD = central finite difference (uses WENO linear weights)
(default: JS)''')
  
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
  
  parser.add_argument('-r','--range',
                      type    = int,
                      nargs   = 3,
                      default = [50,3200,7],
                      metavar = ('MIN','MAX','NPTS'),
                      help    = 'number of cells (default: [50,3200,7])')
  
  return parser.parse_args()

#===============================================================================
# FUNCTION: Main script
#===============================================================================

def main():
  
  # Parse input arguments
  args = parse_input()
  print(args)
  print('')
  
  #-----------------------------------------------------------------------------
  import numpy as np
  
  # Import test-case definition
  import      advection_smooth
  test_case = advection_smooth.DefineTestCase()
  
  # Import path to library
  try               :  import hyperpyws
  except ImportError:  import hyperpyws_path
  
  # Import modules from library
  from hyperpyws.output_utilities  import TextDB
  from hyperpyws.simulation        import Numerics, RunSimulation
  from hyperpyws.weno_versions     import Weno
  
  # Construct array with number of subdivisions
  logN1   = np.log10(args.range[0])
  logN2   = np.log10(args.range[1])
  npts    = args.range[2]
  Nx_list = np.logspace(logN1,logN2,npts).round().astype(int).tolist()
  
  # Extract time-integrator function
  import hyperpyws.time_integrators as integrators
  stepper_name = integrators.__all__[args.stepper]
  stepper_func = getattr( integrators, stepper_name )
  
  # Extract WENO reconstruction class
  weno_class = Weno( args.weno_order, args.weno_version )
  
  # Create container for numerical parameters
  num_params         = Numerics()
  num_params.weno    = weno_class
  num_params.stepper = stepper_func
  num_params.CFL     = args.CFL
  
  # Create file with final error in solution
  ostream = TextDB('error.dat')
  ostream.SetField('mx',  '5d',   desc='Number of grix points')
  ostream.SetField('L1', '2.15e', desc=     'L1 norm of error')
  ostream.SetField('L2', '2.15e', desc=     'L2 norm of error')
  ostream.SetField('Li', '2.15e', desc=  'L-inf norm of error')
  ostream.open()
  
  #-----------------------------------------------------------------------------
  # Run series of simulations
  for n in range( len(Nx_list) ):

    Nx = Nx_list[n]

    # Assign new number of grid cells
    num_params.mx = Nx
    
    # Run new simulation
    print('Running simulation with {:5d} cells... '.format(Nx), end='')
    grid = RunSimulation( test_case, num_params )
    
    # Compute error in numerical solution, and its norms
    Err = test_case.qexact( grid.xint, test_case.tend )[0] - grid.qint[0]
    L1  = np.sum (abs (Err))   *grid.dx
    L2  = np.sqrt(sum((Err)**2)*grid.dx)
    Li  = np.amax(abs (Err))
    
    # Append data to file
    ostream.write( Nx, L1, L2, Li)

    # live-feed print statements for error analysis:
    if( n > 0 ):

        # compute a ratio of the errors (copied from convergence_analysis.py)
        LogRatio_err = np.log( L2_old / L2      )
        LogRatio_mx  = np.log( dx_old / grid.dx )
        order = LogRatio_err / LogRatio_mx
        print(' done: L2-error = %2.3e; Order = %2.3f' % ( L2, order ) )

    else:
        print(' done: L2-error = %2.3e; Order = x.xxx ' % L2 )

    # save the old errors for the live feed:
    L2_old = L2
    dx_old = grid.dx

  #-----------------------------------------------------------------------------
  
  # Close output file
  ostream.close()

#===============================================================================
if __name__ == '__main__':
  #Run as main program
  main()
