import numpy as np

try               :  import hyperpyws #weno_step
except ImportError:  import hyperpyws_path

#===============================================================================
# Test case: 1D constant advection with square-wave ICs and periodic BCs
#===============================================================================
from hyperpyws.model_equations.advection  import Advection1D
from hyperpyws.boundary                   import PeriodicBCs
from hyperpyws.simulation                 import TestCase

# Constant velocity
v = 1.0

# Initial conditions
def q_init (x):
  q0 = 0.0*x
  for i,xi in enumerate(x):
    if( xi > 0.3 and xi < 0.5 ):
      q0[i] = 1.0
  return [ q0 ]

# Exact solution
def q_exact (x,t):
  return q_init( (x-v*t) % 1.0 )

# Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
def CreateBC_func (mx, mbc):
  SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
  return SetBCs

# Test-Case container
test          = TestCase()
test.ModelEqn = Advection1D(v)
test.xlims    = [0.0, 1.0]
test.tend     =  1.0
test.BCs      = CreateBC_func  # better to simply say 'periodic'...
test.qinit    = q_init
test.qexact   = q_exact

#===============================================================================
# Choose numerical schemes and give them proper parameters
#===============================================================================
from hyperpyws.weno_versions.weno5  import Weno5_JS
from hyperpyws.time_integrators     import rk3_ssp, rk4, TD_RK4
from hyperpyws.simulation           import Numerics

# Numerics container
numr         = Numerics()
numr.weno    = Weno5_JS
numr.stepper = rk4
numr.CFL     = 0.5       # CFL parameter
numr.mx      = 100       # number of mesh cells in domain

#===============================================================================
# Real-time visualization
#===============================================================================

# Number of plots during simulation (alternatively, dT should be specified)
nplots = 20

# Time instants for creating an output
Tout = np.linspace( 0.0, test.tend, nplots+1 )

#===============================================================================
# Run simulation
#===============================================================================
from hyperpyws.simulation import RunSimulation

# Call default library function
RunSimulation( test, numr, Tout )

#===============================================================================
# OUTPUT
#===============================================================================
import matplotlib.pyplot as plt

plt.show()
