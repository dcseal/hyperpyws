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

import math
import numpy             as np
import matplotlib.pyplot as plt

try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# INPUT DATA
#===============================================================================

# parameters used in 4bells example:
a     =  0.5
z     = -0.7
delta =  0.005
alpha = 10.0
beta  = np.log10(2.0)/(36*delta**2)

# We suggest using an odd number of points in order to hit the peaks of each
# curve exactly:
Npts = 2001

# final time (needed to print to file):
t = 8.0
OUTFILE = 'exact.dat'

#===============================================================================
# Initial Conditions (copied from advection_4bells):
#===============================================================================

def q_init( x ):

    G = lambda xi, xc : math.exp(-beta*(xi-xc)**2)
    F = lambda xi, xc : math.sqrt(max(1.0-alpha**2*(xi-xc)**2,0.0))
    
    q0 = np.zeros( x.shape )
    
    for i,xi in enumerate( x ):
      if -0.8 <= xi <= -0.6:
        q0[i] = ( G(xi,z-delta) + G(xi,z+delta) + 4.0*G(xi,z) ) / 6.0
      elif -0.4 <= xi <= -0.2:
        q0[i] = 1.0
      elif 0.0 <= xi <= 0.2:
        q0[i] = 1.0-abs(10*(xi-0.1))
      elif 0.4 <= xi <= 0.6:
        q0[i] = ( F(xi,a-delta) + F(xi,a+delta) + 4.0*F(xi,a) ) / 6.0
    
    return q0

def Ffunc( x ):

    F = lambda xi, xc : math.sqrt(max(1.0-alpha**2*(xi-xc)**2,0.0))
    
    f = np.zeros( x.shape )
    f1 = np.zeros( x.shape )
    f2 = np.zeros( x.shape )
    f3 = np.zeros( x.shape )
    
    for i,xi in enumerate( x ):
      if 0.4 <= xi <= 0.6:
        f[i] = ( F(xi,a-delta) + F(xi,a+delta) + 4.0*F(xi,a) ) / 6.0

        f1[i] = ( F(xi,a-delta) ) / 6.0
        f2[i] = ( 4.0*F(xi,a )  ) / 6.0
        f3[i] = ( F(xi,a+delta) ) / 6.0
    
    return f, f1, f2, f3


def Gfunc( x ):

    G = lambda xi, xc : math.exp(-beta*(xi-xc)**2)
    
    g = np.zeros( x.shape )
    g1 = np.zeros( x.shape )
    g2 = np.zeros( x.shape )
    g3 = np.zeros( x.shape )
    
    for i,xi in enumerate( x ):
      if -0.8 <= xi <= -0.6:
        g[i] = ( G(xi,z-delta) + G(xi,z+delta) + 4.0*G(xi,z) ) / 6.0
        g1[i] = ( G(xi,z-delta) ) / 6.0
        g2[i] = ( G(xi,z+delta) ) / 6.0
        g3[i] = ( 4.0*G(xi,z )  ) / 6.0
    
    return g,g1,g2,g3


#===============================================================================
# ARRAYS
#===============================================================================

x1 = np.linspace( -0.8, -0.6, Npts )
g_list  = Gfunc( x1 )
g = g_list[0]

# This section is a bit more complicated, given that there are three sections:
#
# Over the domain [0.4,0.6], there are exactly three pieces that go into F.
# The reason being the hard switch on the max( *, 0.0 ) term.  The three
# pieces are the following:
#
# F = f1 + f2
xa = np.linspace( 0.4,   0.405, Npts   )

# F = f1 + f2 + f3
xb = np.linspace( 0.405, 0.595, Npts   )

# F =      f2 + f3
xc = np.linspace( 0.595, 0.6,   Npts   )

# Can also use these curves for the two small segments:
#xa = np.linspace( 0.4,   0.405,    3   )
#xc = np.linspace( 0.595, 0.6,      3   )

x2 = np.linspace(  0.4,  0.6, Npts )

# The whole section
f_list = Ffunc( x2 )
f      = f_list[0]

#===============================================================================
# 2D LINE
#===============================================================================

from hyperpyws.geometry import Point, LineSegment, CurveSegment, Concatenate

# Points
A = Point( -1.0, 0.0 )
B = Point( -0.8, 0.0 )
C = Point( -0.6, 0.0 )
D = Point( -0.4, 0.0 )
E = Point( -0.4, 1.0 )
F = Point( -0.2, 1.0 )
G = Point( -0.2, 0.0 )
H = Point(  0.0, 0.0 )
I = Point(  0.1, 1.0 )
J = Point(  0.2, 0.0 )
K = Point(  0.4, 0.0 )
L = Point(  0.6, 0.0 )
M = Point(  1.0, 0.0 )


# Segments
AB =  LineSegment( A, B  ) 
BC = CurveSegment( x1, g )
CD =  LineSegment( C,  D )
DE =  LineSegment( D,  E )
EF =  LineSegment( E,  F )
FG =  LineSegment( F,  G )
GH =  LineSegment( G,  H )
HI =  LineSegment( H,  I )
IJ =  LineSegment( I,  J )
JK =  LineSegment( J,  K )

# This curve is a bit more complicated: it has three segments, dealing with
# when max( fi, 0 ) turns on:
fa  = Ffunc(xa)
fb  = Ffunc(xb)
fc  = Ffunc(xc)

#  ( <<< because of the sqrt, we only get F[0.405] = f1+f2 + 10^{-8} <<< )
#        Do we care about this? (-DS)
#
KLa = CurveSegment ( xa,  fa[1]+fa[2]           )

KLb = CurveSegment ( xb,  fb[1]+fb[2]+fb[3]     )

#  ( <<< because of the sqrt, we only get F[0.595] = f2+f3 + 10^{-8} <<< )
#        Do we care about this? (-DS)
#
KLc = CurveSegment ( xc,        fc[2]+fc[3]     )

KL  = CurveSegment ( x2, f )

LM =  LineSegment  ( L, M )

# Line
x, y = Concatenate( AB, BC, CD, DE, EF, GH, HI, IJ, JK, KLa, KLb, KLc, LM )

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

# unpack the axes immediately 
# (see: http://matplotlib.org/examples/pylab_examples/subplots_demo.html )
fig, (ax1, ax2 ) = plt.subplots(2, 1, sharex=True)
#fig = plt.figure()

ax1.plot( x, y, '.-', color='r', linewidth=2, mec='b', mfc='b' )
ax1.grid()
ax1.set_xlabel('x')
ax1.set_ylabel('q',rotation='horizontal')
ax1.set_ylim([-0.1,1.1])

ax2.plot( x1, g_list[1], '.-', color='k', linewidth=2, mec='b', mfc='b' )
ax2.plot( x1, g_list[2], '.-', color='b', linewidth=2, mec='b', mfc='b' )
ax2.plot( x1, g_list[3], '.-', color='g', linewidth=2, mec='b', mfc='b' )
ax2.plot( x2, f_list[1], '.-', color='k', linewidth=2, mec='b', mfc='b' )
ax2.plot( x2, f_list[2], '.-', color='b', linewidth=2, mec='b', mfc='b' )
ax2.plot( x2, f_list[3], '.-', color='g', linewidth=2, mec='b', mfc='b' )
ax2.grid()
ax2.set_xlabel('x')
ax2.set_ylabel('F,G')
ax2.set_ylim([-0.01,0.71])

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
