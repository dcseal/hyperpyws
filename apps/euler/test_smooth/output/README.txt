PARAMS
======
CFL  = 0.9
weno = Weno5_JS
tend = 2.0

VARIABLES
=========
q0 = rho   (density)
q1 = rho*u (momentum)
q2 = eng   (total energy)

METHODS
=======
M0 = rk3-ssp
M1 = rk4
M2 = TD-RK4
