#def f(q,t):
#  qs = ApplyBCs(q,t)
#  qt, qtt = MOL.TimeDerivatives ( qs )
#  return qt, qtt

#===============================================================================
# Single-derivative Runge-Kutta methods
#===============================================================================

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
  k2,dk2 = Fc(Ys,t)
  
  return Y + dt * ( k1 + dt/6. * (dk1 + 2.*dk2) )

