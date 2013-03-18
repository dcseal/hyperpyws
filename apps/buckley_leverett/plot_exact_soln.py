import numpy             as np
import matplotlib.pyplot as plt


try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# INPUT DATA
#===============================================================================

M = 1./3.

qs_left_problem  = 0.1339745962155613
qs_right_problem = 0.5

dt   = 0.4
Npts = 20

#===============================================================================
# FLUX FUNCTION
#===============================================================================

from hyperpyws.model_equations.buckley_leverett  import BuckleyLeverett1D

flux = BuckleyLeverett1D( M )

#===============================================================================
# ARRAYS
#===============================================================================

ql = np.linspace( 0.0,  qs_left_problem, Npts )
qr = np.linspace( 1.0, qs_right_problem, Npts )

left_shock  = -0.5 + dt * flux.eig([ qs_left_problem ])[0]
fl          = -0.5 + dt * flux.eig([ ql              ])[0]

right_shock = dt * flux.eig([ qs_right_problem ])[0]
fr          = dt * flux.eig([ qr               ])[0]

#===============================================================================
# 2D LINE
#===============================================================================

from hyperpyws.geometry import Point, LineSegment, CurveSegment, Concatenate

# Points
A = Point(       -1.0,             0.0 )
B = Point(       -0.5,             0.0 )
C = Point( left_shock, qs_left_problem )
D = Point( left_shock,             1.0 )
E = Point(        0.0,             1.0 )
F = Point(right_shock, qs_right_problem)
G = Point(right_shock,             0.0 )
H = Point(        1.0,             0.0 )

# Segments
AB =  LineSegment( A, B ) 
BC = CurveSegment(fl,ql )
CD =  LineSegment( C, D )
DE =  LineSegment( D, E )
EF = CurveSegment(fr,qr )
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
