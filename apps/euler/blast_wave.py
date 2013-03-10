import numpy as np

try               :  import hyperpyws #weno_step
except ImportError:  import hyperpyws_path

#===============================================================================
# Test case: Blast-wave for 1D Euler's equations with solid-wall BCs
#===============================================================================
from hyperpyws.model_equations.euler  import Euler1D
from hyperpyws.simulation             import TestCase

# Ratio of specific heats
gamma = 1.4

# Initial conditions
def q_init (x):
  
  rho =     np.ones (x.shape)
  u1  =     np.zeros(x.shape)
  p   = 100*np.ones (x.shape)
  
  for i,xi in enumerate(x):
    if xi < 0.1:
      p[i] = 1000.
    elif 0.1 < xi < 0.9:
      p[i] = 0.01
  
  eng  = p/(gamma-1.0) + 0.5*rho*u1**2
  q    = np.empty( 3, dtype=object )
  q[:] = [ rho, rho*u1, eng ]
  
  return q

# Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
def CreateBC_func (mx, mbc):
  def SetBCs(q,t):
    Euler1D.SolidWallBC_left (q,mx,mbc)
    Euler1D.SolidWallBC_right(q,mx,mbc)
  return SetBCs

# Test-Case container
test          = TestCase()
test.ModelEqn = Euler1D (gamma)
test.xlims    = [0.0, 1.0]
test.tend     =  0.038
test.BCs      = CreateBC_func
test.qinit    = q_init

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
numr.CFL     = 0.1      # CFL parameter
numr.mx      = 100       # number of mesh cells in domain

#===============================================================================
# Real-time visualization
#===============================================================================

# Number of plots during simulation (alternatively, dT should be specified)
nplots = 50

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
