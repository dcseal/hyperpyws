from __future__ import print_function   # for saving to file

import numpy             as np
import matplotlib.pyplot as plt


try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# INPUT DATA
#===============================================================================

# gravity
g       = 1.0
Npts    = 10    # number of points used use for filling in rarefaction fan.
OUTFILE = 'exact_soln.dat'

#===============================================================================
# Riemann Data (created from riemann_solve_sw.py)
#===============================================================================

# Initial data (in characteristic variables):
W_l = [3.0, 0.0   ]
W_r = [1.0, 0.0   ]

hl = W_l[0]; ul = W_l[1]
hr = W_r[0]; ur = W_r[1]

# location of Riemann problem and computational domain:
xs   = 0.5 
xmin = 0.
xmax = 1.

# intermediate height and velocity:
hstar = 1.848576603096757e+00
ustar = 7.448542169801264e-01

#===============================================================================
# Location of feet for rarefaction wave, and location of shock
#===============================================================================

# two speeds for the left rarefaction:
sll    =    ul - np.sqrt( g*hl ); 
slr    = ustar - np.sqrt( g*hstar );

# shock speed for right hand shock:
s2 = (hstar*ustar - hr*ur) / (hstar - hr );

#===============================================================================
# Data for left Rarefaction:
#===============================================================================
def LeftFan( x, t, xs, W ):
    """ Data for values inside the RightFan.

    User is responsible for knowing correct bounds of x to supply.
    """

    h_K = W[0]
    u_K = W[1]

    A = u_K + 2.*np.sqrt( g*h_K )

    h = 1./(9.*g)*( (u_K + 2.*np.sqrt(g*h_K)) - (x-xs)/t )**2
    u = u_K - 2.*( np.sqrt(g*h) - np.sqrt(g*h_K) )
    return h,u

t = 0.2

#===============================================================================
# FLUX FUNCTION
#===============================================================================

#from hyperpyws.model_equations.euler  import Euler1D

#flux = ShallowWater1D( )

#===============================================================================
# ARRAYS
#===============================================================================

# Riemann data:
xx = np.linspace( xs+t*sll, xs+t*slr, Npts )
h_fan, u_fan = LeftFan( xx, t, 0.5, W_l )

#left_shock  = -0.5 + dt * flux.eig([ qs_left_problem ])[0]
#fl          = -0.5 + dt * flux.eig([ ql              ])[0]

#right_shock = dt * flux.eig([ qs_right_problem ])[0]
#fr          = dt * flux.eig([ qr               ])[0]

#===============================================================================
# 2D LINE
#===============================================================================

from hyperpyws.geometry import Point, LineSegment, CurveSegment, Concatenate

#----------------#
# Exact Height   #
#----------------#

# Points
A = Point(       xmin,       hl    )
B = Point( xs + t*sll,       hl    )
C = Point( xs + t*slr,       hstar )
D = Point( xs +  t*s2,       hstar )
E = Point( xs +  t*s2,       hr    )
F = Point(       xmax,       hr    )

# Segments
AB =  LineSegment( A, B     ) 
BC = CurveSegment(xx, h_fan )
CD =  LineSegment( C, D     )
DE =  LineSegment( D, E     )
EF =  LineSegment( E, F     )

# exact height:
x, h_ex = Concatenate( AB, BC, CD, DE, EF )

#----------------#
# Exact Velocity
#----------------#

# Points
A2 = Point(       xmin,       ul    )
B2 = Point( xs + t*sll,       ul    )
C2 = Point( xs + t*slr,       ustar )
D2 = Point( xs +  t*s2,       ustar )
E2 = Point( xs +  t*s2,       ur    )
F2 = Point(       xmax,       ur    )

# Segments
AB2 =  LineSegment( A2, B2     ) 
BC2 = CurveSegment( xx, u_fan  )
CD2 =  LineSegment( C2, D2     )
DE2 =  LineSegment( D2, E2     )
EF2 =  LineSegment( E2, F2     )

# exact velocity
x, u_ex = Concatenate( AB2, BC2, CD2, DE2, EF2 )

#===============================================================================
# PLOTS
#===============================================================================

# Water height
x    = np.array( x    )
h    = np.array( h_ex )
fig1 = plt.figure(1)

ax = fig1.add_subplot(1,1,1)
ax.plot( x, h, '.-', color='r', linewidth=2, mec='b', mfc='b' )
ax.grid()
ax.set_xlabel('x')
ax.set_ylabel('q',rotation='horizontal')

fig1.gca().set_ylim([0.9,3.1])
fig1.show()

# Momentum
hu   = np.array(u_ex)*np.array(h_ex)
fig2 = plt.figure(2)

ax = fig2.add_subplot(1,1,1)
ax.plot( x, hu, '.-', color='r', linewidth=2, mec='b', mfc='b' )
ax.grid()
ax.set_xlabel('x')
ax.set_ylabel('q',rotation='horizontal')

fig2.gca().set_ylim([-0.1,1.6])
fig2.show()

#===============================================================================
# PRINT TO FILE
#===============================================================================

if OUTFILE is not None:
    data = np.column_stack( [x,h,hu] )
    fmt  = '%.15e'
    with open( OUTFILE, 'wb' ) as f:

        t_str = '%.15e' % t
        ### <<< TODO - THIS FIRST LINE CAUSES AN ERROR: <<< ###
        #print( fmt % t, file=f )         # time instant on first row
        print( t_str, file=f )         # time instant on first row
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
