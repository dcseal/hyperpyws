import numpy as np

try               :  import hyperpyws #weno_step
except ImportError:  import hyperpyws_path

#===============================================================================
# Test case: Shock-tube (SOD) for 1D Euler's equations with outflow BCs
#===============================================================================
from hyperpyws.model_equations.euler  import Euler1D
from hyperpyws.boundary               import OutflowBC_left, OutflowBC_right
from hyperpyws.simulation             import TestCase

# Ratio of specific heats
gamma = 1.4
eps   = 0.2

# Initial conditions
def q_init (x):

#   if(x<-4.0e0)
#   {
#    37             rho   = 3.857143;
#    38             u1    = 2.629369;
#    39             u2    = 0.0;
#    40             u3    = 0.0;
#    41             press = 10.3333;
#    42         }
#    43         else
#    44         {
#    45             rho   = 1.0 + eps*sin(5.0*x);
#    46             u1    = 0.0;
#    47             u2    = 0.0;
#    48             u3    = 0.0;
#    49             press = 1.0;
#    50         }

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

# Boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
def CreateBC_func (mx, mbc):
  def SetBCs(q,t):
    OutflowBC_left (q,mx,mbc)
    OutflowBC_right(q,mx,mbc)
  return SetBCs

# Test-Case container
test          = TestCase()
test.ModelEqn = Euler1D (gamma)
test.xlims    = [-5.0, 5.0]
test.tend     =  1.8
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
numr.stepper = TD_RK4
numr.CFL     = 0.4       # CFL parameter
numr.mx      = 1600       # number of mesh cells in domain

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
