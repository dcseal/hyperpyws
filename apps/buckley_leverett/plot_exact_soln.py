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

from __future__ import print_function

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

t    = 0.4
Npts = 1000

OUTFILE = None

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

left_shock  = t * flux.eig([ qs_left_problem  ])[0] - 0.5
fl          = t * flux.eig([ ql               ])[0] - 0.5

right_shock = t * flux.eig([ qs_right_problem ])[0]
fr          = t * flux.eig([ qr               ])[0]

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
# PRINT TO FILE
#===============================================================================

if OUTFILE is not None:
  data = np.column_stack( [x,y] )
  fmt  = '%.15e'
  with open( OUTFILE, 'wb' ) as f:
    print( fmt % t, file=f )         # time instant on first row
    np.savetxt( f, data, fmt=fmt )      # data arrays along columns

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
