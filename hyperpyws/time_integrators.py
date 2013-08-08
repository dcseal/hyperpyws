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

#def f(q,t):
#  qs = ApplyBCs(q,t)
#  qt, qtt = MOL.TimeDerivatives ( qs )
#  return qt, qtt

__all__ = ['fE', 'rk2_midpoint', 'rk2_Heun', 'rk3', 'rk3_ssp', 'rk4',
           'Fehlberg5', 'Taylor2', 'TD_RK4']

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
def TD_RK4 (Fc, Y, t, dt):
  """ Two-Derivative, 4th-order Runge-Kutta method.
  """
  k1,dk1 = Fc(Y, t)
  Ys     = Y + (0.5*dt) * ( k1 + (0.25*dt)* dk1 )
  k2,dk2 = Fc(Ys,t+0.5*dt)
  
  return Y + dt * ( k1 + dt/6. * (dk1 + 2.*dk2) )


