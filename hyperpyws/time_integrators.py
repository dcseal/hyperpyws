#==============================================================================#
# This file is part of HYPERPYWS: Hyperbolic Python WENO Solver
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

#def f(q,t):
#  qs = ApplyBCs(q,t)
#  qt, qtt = MOL.TimeDerivatives ( qs )
#  return qt, qtt

__all__ = ['fE', 'rk2_midpoint', 'rk2_Heun', 'rk3', 'rk3_ssp', 'rk4',
           'Fehlberg5', 'Taylor2', 'TD_RK3', 'TD_RK4', 'TD_RK5']

#===============================================================================
# Single-derivative Runge-Kutta methods
#===============================================================================

def fE(f,Y,t,dt):
  """ Classical Forward Euler.
  """
  return Y + dt*f(Y,t)

#-------------------------------------------------------------------------------
def rk2_midpoint (f,Y,t,dt):
  """ Classical 2nd-order Runge-Kutta (midpoint rule). 
  """
  k1 = f(Y          , t       )
  k2 = f(Y+0.5*dt*k1, t+0.5*dt)
  return Y + dt*k2

#-------------------------------------------------------------------------------
def rk2_Heun (f,Y,t,dt):
  """ Strong-stability-preserving (SSP) 2nd-order Runge-Kutta (Heun's method).
  """
  k1 = f(Y      , t   )
  k2 = f(Y+dt*k1, t+dt)
  return Y + 0.5*dt*( k1 + k2 )

#-------------------------------------------------------------------------------
def rk3 (f,Y,t,dt):
  """ Classical 3rd-order Runge-Kutta.
  """
  k1 = f(Y                 , t          )
  k2 = f(Y +     0.5*dt*k1 , t + 0.5*dt )
  k3 = f(Y + dt*(-k1+2.*k2), t +     dt )
  return Y + (dt/6.)*( k1 + 4.0*k2 + k3 )

#-------------------------------------------------------------------------------
def rk3_ssp (f,Y,t,dt):
  """ Strong-stability-preserving (SSP) 3rd-order Runge-Kutta.
  """
  Y1 = Y + dt*f(Y,t)
  Y2 =   (3.*Y +    Y1 +    dt*f(Y1,t+dt))*0.25
  return (   Y + 2.*Y2 + 2.*dt*f(Y2,t+dt))/3.
  
#-------------------------------------------------------------------------------
def rk4 (f,Y,t,dt):
  """ Classical 4th-order Runge-Kutta.
  """
  k1 = f(Y          , t       )
  k2 = f(Y+0.5*dt*k1, t+0.5*dt)
  k3 = f(Y+0.5*dt*k2, t+0.5*dt)
  k4 = f(Y+    dt*k3, t+    dt)
  return Y + (dt/6.)*( k1 + 2.0*k2 + 2.0*k3 + k4 )

def Fehlberg5(f, Y, t, dt ):
    """RK-Fehlberg 5 method. """

    # butcher tableau:
    A = [ 
[                                                               ], 
[      0.25                                                     ], 
[     3./32.,       9./32.                                      ], 
[1932./2197., -7200./2197.,  7296./2197.                        ],
[  439./216.,          -8.,   3680./513.,   -845./4104.         ],
[    -8./27.,           2., -3544./2565.,  1859./4104., -11./40.] ]
  
    # time values:
    c = [0., 0.25, 0.375, 12./13., 1., 0.5]

    k  = [ f( Y, t ) ]

    for i in range(1,6):
        ys = Y + dt* ( sum( [a*ki for (a,ki) in zip( A[i], k)] )  )
        ts = t + dt*c[i]
        k.append( f( ys, ts ) )

    # coefficients for the fifth-order scheme:
    b = [16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.]

    # coefficients for the fourth-order scheme:
    # ( These haven't been tested ... )
    # bs = [25./216., 0., 1408./2565., 2197./4104., -0.2, 0.]

    return Y + dt * sum( [bi*ki for (bi,ki) in zip(b,k)] )



#===============================================================================
# Two-derivative methods
#===============================================================================

def Taylor2 (Fc, Y, t, dt):
  """ Two-derivative, 2nd-order explicit Taylor method.
  """
  k1,dk1 = Fc(Y,t)
  
  return Y + dt * ( k1 + (0.5*dt)* dk1 )

#-------------------------------------------------------------------------------
def TD_RK3 (Fc, Y, t, dt):
  """ Two-Derivative, 3rd-order Runge-Kutta method.
  """

  k1,dk1 = Fc(Y, t)
  Ys     = Y + (dt) * ( k1 + (0.5*dt)* dk1 )
  k2,dk2 = Fc(Ys,t+1.0*dt)
  
  return Y + dt * ( (2.*k1+k2)/3. + dt/6. * (dk1) )

#-------------------------------------------------------------------------------
def TD_RK4 (Fc, Y, t, dt):
  """ Two-Derivative, 4th-order Runge-Kutta method.
  """
  k1,dk1 = Fc(Y, t)
  Ys     = Y + (0.5*dt) * ( k1 + (0.25*dt)* dk1 )
  k2,dk2 = Fc(Ys,t+0.5*dt)
  
  return Y + dt * ( k1 + dt/6. * (dk1 + 2.*dk2) )

#-------------------------------------------------------------------------------
def TD_RK5 (Fc, Y, t, dt):
  """ Two-Derivative, 5th-order Runge-Kutta method.
  """

  # Single parameter family defined by c3:
  c3 = 1.0
  #c3 = .8    # << this value should be unstable <<

  # Remaining coefficients are defined by c3
  bs1 = (10.*c3**2-8.*c3+1.)/(12.*c3*(5.*c3-3.))
  bs2 = (25.*(2.*c3-1.)**3)  /(12.*(5.*c3-3.)*(10.*c3**2-10.*c3+3.))
  bs3 = 1./(12.*c3*(10.*c3**2-10.*c3+3.))
  c2  = (5.*c3-3.)/(5.*(2.*c3-1.))

  a32 = (c3*(2.*c3-1.)*(10.*c3**2-10.*c3+3.))/(2.*(5.*c3-3.))
  a31 = 0.5*c3**2 - a32
  a21 = 0.5*c2**2

  # First stage:
  k1,dk1 = Fc(Y, t)

  # Second stage:
  Y2     = Y + c2*dt*k1 + (a21*dt**2)*dk1
  k2,dk2 = Fc(Y2,t+c2*dt)

  # Third stage:
  Y3     = Y + dt*c3*k1 + dt**2*( a31*dk1+a32*dk2 )
  k3,dk3 = Fc(Y3, t + c3*dt)

  # Final update:
  return Y + dt * ( k1 + dt * (bs1*dk1 + bs2*dk2 + bs3*dk3) )

#===============================================================================
# Low-storage SSP methods
#===============================================================================

#-------------------------------------------------------------------------------
#def ssp10_4(F, Y, t, dt):
# """ Fourth-order, ten stage method.

# This is a low-storage, optimal SSP(10,4) method.  

# See: "Highly efficient strong stability-preserving Runge-Kutta methods with 
#        low-storage implementations", Ketcheson, (2008).

# for further details.


# In particular, we follow the low-storage implementation
# suggested in Pseudocode 3.  There's a single "stage" that
# reassigns each register with a linear combination of the
# two registers being used.
# """

# def AdvanceTimeStage( alpha1, q, alpha2, qnew, beta_dt, L ):
#   
#   qnew = alpha1*q    + alpha2*qnew + beta_dt*L
#   t    = alpha1*told + alpha2*tnew + beta_dt

#   return qnew

# def AdvanceTimeStage4( told, qold, q1, q2 ):
#   """See 'Pseudocode 3.' in "Highly efficient strong stability-preserving
#      Runge-Kutta methods with low-storage implementations", (2008)."""

#   # time advance:
#   #   t2 = (told + 9.0*t1)*(1./25.);
#   #   t1 = 15.0*t2 - 5.0*t1;          # i.e., t1 = .6*told+.4*t1

#   q2 = (qold + 9.0*q1)/(25.0)
#   q1 = 15.0*q2 - 5.0*q1

#   return (q1,q2)

# # coefficients defining the method:
# alpha1 = np.ones(10)
# alpha2 = np.zeros(10);          alpha2[9] = 3.0/5.0
# beta   = 1.0/6.0*np.ones(10);   beta  [9] = 0.1

# def AdvanceTimeStage( alpha1, q, alpha2, qnew, beta_dt, L ):

# <<< TODO : finish writing the small section <<<

# for n in range(5):
#   (q1,q2) = AdvanceTimeStage( alpha1[n], q1, alpha2[n], q1, dt*beta[n], f(q1,t) )

# (q1,q2) = AdvanceTimeStage4( dt, Y, q1, q2 )

# for n in range(5,9):
#   (q1,q2) = AdvanceTimeStage( alpha1[n], q1, alpha2[n], q2, dt*beta[n], f(q1,t) )

# (q1,q2) = AdvanceTimeStage( alpha1[n], q1, alpha2[n], q2, dt*beta[n], f(q1,t) )

# return Y + dt*k2
