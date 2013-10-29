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

from __future__ import print_function   # for saving to file

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

t       = 0.16                 # final time of simulation
Npts    = 2000                 # number of points used for sampling the rarefaction
OUTFILE = 'exact_soln.dat'  # name of output file


#===============================================================================
# Riemann Data (created from riemann_solve_euler.py)
#===============================================================================

# Left and right-hand Riemann states:
W_l = [0.445, 0.6991011235955056, 3.5277019280898867]
W_r = [0.5, 0.0, 0.5709999999999998]

rho_l, u_l, p_l = W_l
rho_r, u_r, p_r = W_r


# location of Riemann problem:
xs = 0.5 

# intermediate pressure and velocity:
p_star = 2.466711993414689e+00
u_star = 1.529035097908221e+00

# two speeds for the left rarefaction:
sll    = -2.632323209684882e+00
slr    = -1.636402440509624e+00

# intermediate values for the density
rho_sl = 3.446505568509241e-01
rho_sr = 1.304261261849049e+00 

# shock speed for right hand shock:
s3 = 2.479618677175032e+00

#===============================================================================
# Values derived from solution of intermediate state:
#===============================================================================

# momentum:
rho_u_l  = rho_l* u_l
rho_u_sl = rho_sl*u_star
rho_u_sr = rho_sr*u_star
rho_u_r  = rho_r* u_r

# energy:
# energy = 0.5*rho*u**2 + press/gm1
Energy_l  = 0.5*rho_u_l*u_l     + p_l/gm1
Energy_sl = 0.5*rho_u_sl*u_star + p_star/gm1
Energy_sr = 0.5*rho_u_sr*u_star + p_star/gm1
Energy_r  = 0.5*rho_u_r *u_r    + p_r   /gm1


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
    u     = 2.0/gp1 * ( c_l + 0.5*gm1*u_l + (x-xs)/t )
    press = press_l * ( 2./gp1 + (gm1/(c_l*gp1))*( u_l - (x-xs)/t ) )**(2*gamma/gm1)

    energy = 0.5*rho*u**2 + press/gm1
    return [rho, rho*u, energy]


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

xx = np.linspace( xs+t*sll, xs+t*slr, Npts )
rho_fan, rho_u_fan, energy_fan = LeftFan( xx, t, W_l )

#===============================================================================
# 2D LINE
#===============================================================================

from hyperpyws.geometry import Point, LineSegment, CurveSegment, Concatenate

# Density:

# Points
A1 = Point(      xmin,      rho_l    )
B1 = Point( xs+t*sll,       rho_l    )
C1 = Point( xs+t*slr,       rho_sl   )
D1 = Point( xs+t*u_star,    rho_sl   )
E1 = Point( xs+t*u_star,    rho_sr   )
F1 = Point( xs+t*s3,        rho_sr   )
G1 = Point( xs+t*s3,        rho_r    )
H1 = Point( xmax,           rho_r    )

# Segments
AB1 =  LineSegment( A1, B1     ) 
BC1 = CurveSegment(xx, rho_fan )
CD1 =  LineSegment( C1, D1     )
DE1 =  LineSegment( D1, E1     )
EF1 =  LineSegment( E1, F1     )
FG1 =  LineSegment( F1, G1     )
GH1 =  LineSegment( G1, H1     )

# Line
x_rho, rho_ex = Concatenate( AB1, BC1, CD1, DE1, EF1, FG1, GH1 )

# Velocity:
A2 = Point(      xmin,      rho_u_l     )
B2 = Point( xs+t*sll,       rho_u_l     )
C2 = Point( xs+t*slr,       rho_u_sl    )
D2 = Point( xs+t*u_star,    rho_u_sl    )
E2 = Point( xs+t*u_star,    rho_u_sr    )
F2 = Point( xs+t*s3,        rho_u_sr    )
G2 = Point( xs+t*s3,        rho_u_r     )
H2 = Point( xmax,           rho_u_r     )

# Segments
AB2 =  LineSegment( A2, B2   ) 
BC2 = CurveSegment(xx, rho_u_fan )
CD2 =  LineSegment( C2, D2   )
DE2 =  LineSegment( D2, E2   )
EF2 =  LineSegment( E2, F2   )
FG2 =  LineSegment( F2, G2   )
GH2 =  LineSegment( G2, H2   )

# Line
x_rho_u, rho_u_ex = Concatenate( AB2, BC2, CD2, DE2, EF2, FG2, GH2 )

# Energy:
A3 = Point(      xmin,      Energy_l     )
B3 = Point( xs+t*sll,       Energy_l     )
C3 = Point( xs+t*slr,       Energy_sl    )
D3 = Point( xs+t*u_star,    Energy_sl    )
E3 = Point( xs+t*u_star,    Energy_sr    )
F3 = Point( xs+t*s3,        Energy_sr    )
G3 = Point( xs+t*s3,        Energy_r     )
H3 = Point( xmax,           Energy_r     )

# Segments
AB3 =  LineSegment( A3, B3   ) 
BC3 = CurveSegment(xx, energy_fan )
CD3 =  LineSegment( C3, D3   )
DE3 =  LineSegment( D3, E3   )
EF3 =  LineSegment( E3, F3   )
FG3 =  LineSegment( F3, G3   )
GH3 =  LineSegment( G3, H3   )

# Line
x_energy, energy_ex = Concatenate( AB3, BC3, CD3, DE3, EF3, FG3, GH3 )


#===============================================================================
# PLOTS
#===============================================================================

# Figure 1 (Density)
fig1 = plt.figure()

ax = fig1.add_subplot(1,1,1)
ax.plot( x_rho, rho_ex, '.-', color='r', linewidth=2, mec='b', mfc='b' )
ax.grid()
ax.set_xlabel('x')
ax.set_ylabel('Density',rotation='horizontal')

fig1.gca().set_ylim([0.2,1.4])
fig1.show()

# Figure 2 (Momentum)
fig2 = plt.figure()

ax = fig2.add_subplot(1,1,1)
ax.plot( x_rho_u, rho_u_ex, '.-', color='r', linewidth=2, mec='b', mfc='b' )
ax.grid()
ax.set_xlabel('x')
ax.set_ylabel('Momentum',rotation='horizontal')

#fig2.gca().set_ylim([-0.5,18.5])
fig2.show()

# Figure 3 (Energy)
fig3 = plt.figure()

ax = fig3.add_subplot(1,1,1)
ax.plot( x_energy, energy_ex, '.-', color='r', linewidth=2, mec='b', mfc='b' )
ax.grid()
ax.set_xlabel('x')
ax.set_ylabel('Energy',rotation='horizontal')

#fig3.gca().set_ylim([0.0,260])
fig3.show()

#===============================================================================
# PRINT TO FILE
#===============================================================================


if OUTFILE is not None:
    data = np.column_stack( [x_energy, rho_ex, rho_u_ex, energy_ex] )
    fmt  = '%.15e'
    with open( OUTFILE, 'wb' ) as f:
        print( fmt % t, file=f )         # time instant on first row
        np.savetxt( f, data, fmt=fmt )   # data arrays along columns



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
