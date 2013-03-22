# coding: utf8

from __future__ import print_function

#===============================================================================
# FUNCTION: Parse input arguments
#===============================================================================

def parse_input():
  
  import argparse, sys
  
  parser = argparse.ArgumentParser (
      prog='python '+sys.argv[0],
      description='Run 1D Euler simulation with increasing number of cells',
      formatter_class=argparse.RawTextHelpFormatter
      )
  
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
  import      euler_smooth
  test_case = euler_smooth.DefineTestCase()
  
  # Import path to library
  try               :  import hyperpyws
  except ImportError:  import hyperpyws_path
  
  # Import modules from library
  from hyperpyws.output_utilities  import TextDB
  from hyperpyws.simulation        import Numerics, RunSimulation
  from hyperpyws.weno_versions     import weno5
  
  # Construct array with number of subdivisions
  logN1   = np.log10(args.range[0])
  logN2   = np.log10(args.range[1])
  npts    = args.range[2]
  Nx_list = np.logspace(logN1,logN2,npts).round().astype(int).tolist()
  
  # Extract time-integrator function
  import hyperpyws.time_integrators as integrators
  stepper_name = integrators.__all__[args.stepper]
  stepper_func = getattr( integrators, stepper_name )
  
  # Create container for numerical parameters
  num_params         = Numerics()
  num_params.weno    = weno5.Weno5_JS
  num_params.stepper = stepper_func
  num_params.CFL     = args.CFL
  
  # Create files with final error in solution
  ostreams = [ TextDB('error_q0.dat'),
               TextDB('error_q1.dat'),
               TextDB('error_q2.dat') ]
  for os in ostreams:
    os.SetField('mx',  '5d', desc='Number of grix points')
    os.SetField('L1', '.3e', desc=     'L1 norm of error')
    os.SetField('L2', '.3e', desc=     'L2 norm of error')
    os.SetField('Li', '.3e', desc=  'L-inf norm of error')
    os.open()
  
  #-----------------------------------------------------------------------------
  # Run series of simulations
  for Nx in Nx_list:
    
    # Assign new number of grid cells
    num_params.mx = Nx
    
    # Run new simulation
    print('Running simulation with {:5d} cells... '.format(Nx), end='')
    grid = RunSimulation( test_case, num_params )
    
    # Compute error in numerical solution, and its norms
    for i,os in enumerate( ostreams ):
      
      Err = test_case.qexact( grid.xint, test_case.tend )[i] - grid.qint[i]
      L1  = np.sum (abs (Err))   *grid.dx
      L2  = np.sqrt(sum((Err)**2)*grid.dx)
      Li  = np.amax(abs (Err))
      
      # Append data to file
      os.write( Nx, L1, L2, Li )
    
    print(' done.')
  #-----------------------------------------------------------------------------
  
  # Close output files
  for os in ostreams:
    os.close()

#===============================================================================
if __name__ == '__main__':
  #Run as main program
  main()
