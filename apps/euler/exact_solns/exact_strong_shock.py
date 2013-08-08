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
W_l = [10., 0., 100.]
W_r = [1., 0., 1.]
rho_l = W_l[0]
rho_r = W_r[0]


# location of Riemann problem:
xs = 2.0 

# intermediate pressure and velocity:
p_star = 1.990857789704076e+01
u_star = 3.852457192610188e+00

# two speeds for the left rarefaction:
sll    = -3.741657386773941e+00
slr    = 8.812912443582839e-01

# intermediate values for the density
rho_sl = 3.157289870454345e+00 
rho_sr = 4.649096058491207e+00

# shock speed for right hand shock:
s3 = 4.908186373442732e+00

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


time = 0.4
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
xmax = 5.

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

fig.gca().set_ylim([-0.1,10.1])
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
