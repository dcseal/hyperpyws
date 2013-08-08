#==============================================================================#
# This file is part of HYPERPYWS: Hyperbolic Python WENO Solver
#
#
#   This software is made available for research and instructional use only.
#   You may copy and use this software without charge for these non-commercial
#   purposes, provided that the copyright notice and associated text is
#   reproduced on all copies.  For all other uses (including distribution of
#   modified versions), please contact the author at the address given below.
# 
#   *** This software is made available "as is" without any assurance that it
#   *** will work for your purposes.  The software may in fact have defects, so
#   *** use the software at your own risk.
#
# License: GPL, see COPYING for details
#
# Copyright (C) 2013 
#
#    David Seal,  seal@math.msu.edu,  Michigan State University
#    Yaman Guclu, guclu@math.msu.edu, Michigan State University
#
#===============================================================================

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
  
  parser.add_argument('-O','--weno_order',
                      type    = int,
                      choices = [5,7],
                      default =  5,
                      help    = 'order of accuracy for WENO recontruction'+\
                                ' (default: 5)')
  
  parser.add_argument('-w','--weno_version',
                      choices = ['JS','Z','CFD'],
                      default =  'Z',
                      help    = 
  '''choose WENO version:
  JS  = WENO-JS (Jiang-Shu's algorithm)
  Z   = WENO-Z  (Borges-Carmona-Costa-Don's algorithm)
  CFD = central finite difference (uses WENO linear weights)
(default: Z)''')
  
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

  parser.add_argument('-v','--verbosity',
                      type    = bool,
                      default = False,
                      metavar = 'V',
                      help    = 'if( V ), turn on verbose output to screen.'
                                ' (default: False)')
  
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
  
  import hyperpyws.time_integrators as integrators
  from   hyperpyws.weno_versions    import Weno
  from   hyperpyws.simulation       import Numerics, RunSimulation
  
  # Extract time-integrator function
  stepper_name = integrators.__all__[args.stepper]    # --> dangerous indexing
  stepper_func = getattr( integrators, stepper_name ) 
  
  # Extract WENO reconstruction class
  weno_class = Weno( args.weno_order, args.weno_version )
  
  # Collect input numerical parameters
  num_params         = Numerics()     # Container for numerical parameters
  num_params.weno    = weno_class     # Weno class used for reconstructions
  num_params.stepper = stepper_func   # --> should we pass a string?
  num_params.CFL     = args.CFL       # CFL parameter
  num_params.mx      = args.mx        # Number of mesh cells in domain
  
  # Real-time visualization: time instants for creating an output
  if args.frames is not None:
    Tout = np.linspace( 0.0, test_case.tend, args.frames+1 )
  else:
    Tout = []
  
  # Run simulation: call default library function
  grid = RunSimulation( test_case, num_params, Tout, args.verbosity )
  
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
