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
  
  parser.add_argument('input_file',
                      metavar = 'TEST',
                      help    = 'input file containing test-case definition')
  
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
  import os, sys, numpy as np
  
  # Determine directory and name of input file
  file_dir, file_name = os.path.split( args.input_file )
  
  # Move into test-case directory
  origin = os.path.abspath( os.path.curdir )
  os.chdir( file_dir )
  
  # Add directory to import path
  sys.path.insert( 0, os.path.curdir )
  
  # Import input file as module, and create test-case object
  test_case = __import__(file_name.rstrip('.py')).DefineTestCase()
  
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
  
  # Create files with final error in solution
  ostreams = [ TextDB('error_q0.dat'),
               TextDB('error_q1.dat'),
               TextDB('error_q2.dat') ]
  for f in ostreams:
    f.SetField('mx',  '5d', desc='Number of grix points')
    f.SetField('L1', '.3e', desc=     'L1 norm of error')
    f.SetField('L2', '.3e', desc=     'L2 norm of error')
    f.SetField('Li', '.3e', desc=  'L-inf norm of error')
    f.open()
  
  #-----------------------------------------------------------------------------
  # Run series of simulations
  for Nx in Nx_list:
    
    # Assign new number of grid cells
    num_params.mx = Nx
    
    # Run new simulation
    print('Running simulation with {:5d} cells... '.format(Nx), end='')
    grid = RunSimulation( test_case, num_params )
    
    # Compute error in numerical solution, and its norms
    for i,f in enumerate( ostreams ):
      
      Err = test_case.qexact( grid.xint, test_case.tend )[i] - grid.qint[i]
      L1  = np.sum (abs (Err))   *grid.dx
      L2  = np.sqrt(sum((Err)**2)*grid.dx)
      Li  = np.amax(abs (Err))
      
      # Append data to file
      f.write( Nx, L1, L2, Li )
    
    print(' done.')
  #-----------------------------------------------------------------------------
  
  # Close output files
  for f in ostreams:
    f.close()
  
  # Move back to original directory
  os.chdir( origin )

#===============================================================================
if __name__ == '__main__':
  #Run as main program
  main()
