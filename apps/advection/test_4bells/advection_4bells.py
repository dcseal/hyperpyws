try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FUNCTION: Test-case setup
#===============================================================================

def DefineTestCase ():
  """ 
  Test case definition: 1D constant advection with composite initial 
  conditions and periodic boundary conditions.  The solution contains:
  
    1. smooth but narrow combination of Gaussians;
    2. square wave;
    3. sharp triangle wave;
    4. half ellipse.
  
  See also
  --------
  G.-S. Jiang and C.-W. Shu, "Efficient implementation of Weighted ENO schemes".
  Journal of Computational Physics. 126 (1996), pp. 202-228.
  
  """
  # Import external modules
  import math, numpy as np
  
  # Import modules from library
  from hyperpyws.model_equations.advection  import Advection1D
  from hyperpyws.boundary                   import PeriodicBCs
  from hyperpyws.simulation                 import TestCase
  
  # Constant velocity
  v = 1.0
  
  # Initial conditions
  def q_init (x):
    a     =  0.5
    z     = -0.7
    delta =  0.005
    alpha = 10.0
    beta  = np.log10(2.0)/(36*delta**2)
    
    G = lambda xi, xc : math.exp(-beta*(xi-xc)**2)
    F = lambda xi, xc : math.sqrt(max(1.0-alpha**2*(xi-xc)**2,0.0))
    
    q0 = np.zeros( x.shape )
    
    for i,xi in enumerate( x ):
      if -0.8 <= xi <= -0.6:
        q0[i] = ( G(xi,z-delta) + G(xi,z+delta) + 4.0*G(xi,z) ) / 6.0
      elif -0.4 <= xi <= -0.2:
        q0[i] = 1.0
      elif 0.0 <= xi <= 0.2:
        q0[i] = 1.0-abs(10*(xi-0.1))
      elif 0.4 <= xi <= 0.6:
        q0[i] = ( F(xi,a-delta) + F(xi,a+delta) + 4.0*F(xi,a) ) / 6.0
    
    return [ q0 ]
  
  # Exact solution
  def q_exact (x,t):
    xmin = -1.0
    xmax = +1.0
    return q_init( ((x-v*t)-xmin) % (xmax-xmin) + xmin )
  
  # Periodic boundary conditions  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  (FIX ME)
  def CreateBC_func (mx, mbc):
    SetBCs = lambda q,t : PeriodicBCs(q,mx,mbc)
    return SetBCs
  
  # Test-Case container
  test          = TestCase()
  test.ModelEqn = Advection1D(v)
  test.xlims    = [-1.0, 1.0]
  test.tend     = 8.0
  test.BCs      = CreateBC_func
  test.qinit    = q_init
  test.qexact   = q_exact
    
  return test

#===============================================================================
# SCRIPT: Run as main program
#===============================================================================
help_message = 'Run 1D advection eqn. with 4 different bell-like shapes as ICs.'

if __name__ == '__main__':
  from hyperpyws.interactive import main
  main( help_message, DefineTestCase )
