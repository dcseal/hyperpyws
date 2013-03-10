import numpy as np

try               :  import hyperpyws #weno_step
except ImportError:  import hyperpyws_path

#===============================================================================
# Test case: 1D Euler's equations with smooth initial cond. and periodic BCs
#===============================================================================
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

# Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
def CreateBC_func (mx, mbc):
  SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
  return SetBCs

# Test-Case container
test          = TestCase()
test.ModelEqn = Euler1D (gamma)
test.xlims    = [0.0, 2.0]
test.tend     =  2.0
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
