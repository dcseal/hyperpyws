import numpy             as np
import matplotlib.pyplot as plt


try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# INPUT DATA
#===============================================================================

# gamma, and shortcut names to common terms found here:
gamma = 1.4
gp1 = gamma+1.
gm1 = gamma-1.

#===============================================================================
# Riemann Data (created from riemann_solve_euler.py)
#===============================================================================

# left and right states:
W_l = [1.,    0., 1.0]
W_r = [0.125, 0., 0.1]
rho_l = W_l[0]
rho_r = W_r[0]


# location of Riemann problem:
xs = 0.5 

# intermediate pressure and velocity:
p_star = 3.031301780506468e-01
u_star = 9.274526200489499e-01

# two speeds for the left rarefaction:
sll    = -1.183215956619923e+00
slr    = -7.027281256118334e-02

# intermediate values for the density
rho_sl = 4.263194281784952e-01
rho_sr = 2.655737117053071e-01

# shock speed for right hand shock:
s3 = 1.752155732030178e+00

#===============================================================================
# Data for left Rarefaction:
#===============================================================================
def LeftFan( x, t, W_l ):
    """ Data in the case of a left rarefaction fan.

    User is responsible for knowing correct bounds of x to supply.
    """

    rho_l     = W_l[0]
    u_l       = W_l[1] 
    press_l   = W_l[2]

    c_l = np.sqrt( abs( gamma*press_l / rho_l ) )

    rho   = rho_l   * ( 2./gp1 + (gm1/(c_l*gp1))*( u_l - (x-xs)/t ) )**(  2./gm1)
    u     = 2/gm1   * ( c_l + gm1/2 *u_l + x/t )
    press = press_l * ( 2./gp1 + (gm1/(c_l*gp1))*( u_l - (x-xs)/t ) )**(2*gamma/gm1)

    return [rho, u, press]


time = 0.2
Npts = 20

#===============================================================================
# FLUX FUNCTION
#===============================================================================

#from hyperpyws.model_equations.euler  import Euler1D

#flux = Euler1D( )

#===============================================================================
# ARRAYS
#===============================================================================

xmin = 0.
xmax = 1.

xx = np.linspace( xs+time*sll, xs+time*slr )
rho_xx, JUNK1, JUNK2 = LeftFan( xx, time, W_l )

#left_shock  = -0.5 + dt * flux.eig([ qs_left_problem ])[0]
#fl          = -0.5 + dt * flux.eig([ ql              ])[0]

#right_shock = dt * flux.eig([ qs_right_problem ])[0]
#fr          = dt * flux.eig([ qr               ])[0]

#===============================================================================
# 2D LINE
#===============================================================================

from hyperpyws.geometry import Point, LineSegment, CurveSegment, Concatenate

# Points
A = Point(      xmin,         rho_l    )
B = Point( xs+time*sll,       rho_l    )
C = Point( xs+time*slr,       rho_sl   )
D = Point( xs+time*u_star,    rho_sl   )
E = Point( xs+time*u_star,    rho_sr   )
F = Point( xs+time*s3,        rho_sr   )
G = Point( xs+time*s3,        rho_r    )
H = Point( xmax,              rho_r    )

# Segments
AB =  LineSegment( A, B ) 
BC = CurveSegment(xx, rho_xx )
CD =  LineSegment( C, D )
DE =  LineSegment( D, E )
EF =  LineSegment( E, F )
FG =  LineSegment( F, G )
GH =  LineSegment( G, H )

# Line
x, y = Concatenate( AB, BC, CD, DE, EF, FG, GH )

#===============================================================================
# PLOTS
#===============================================================================

fig = plt.figure()

ax = fig.add_subplot(1,1,1)
ax.plot( x, y, '.-', color='r', linewidth=2, mec='b', mfc='b' )
ax.grid()
ax.set_xlabel('x')
ax.set_ylabel('q',rotation='horizontal')

fig.gca().set_ylim([-0.1,1.1])
fig.show()

#===============================================================================
# SCRIPT
#===============================================================================

# If run as a script from a non-interactive Python session, keep windows open
if __name__=='__main__':
  try: __IPYTHON__ 
  except NameError:
    import sys
    if not sys.flags.interactive:
      from matplotlib.pyplot import show
      show()
