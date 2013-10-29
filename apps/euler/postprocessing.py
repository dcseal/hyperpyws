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

"""
Postprocessing script for Euler's equations.  Still lacks analytical solution.
Loads numerical data from .txt file and creates plots of all quantities.
Intended usage from ipython interactive session.

Examples
--------
  In [1]: run postprocessing.py
  
  In [2]: figs
  
  Out[2]: 
  {'H': <matplotlib.figure.Figure at 0x10c6ade10>,
   'M': <matplotlib.figure.Figure at 0x10cc547d0>,
   'c': <matplotlib.figure.Figure at 0x10cc82f10>,
   'eng': <matplotlib.figure.Figure at 0x10b3acfd0>,
   'mom': <matplotlib.figure.Figure at 0x10b38a1d0>,
   'p': <matplotlib.figure.Figure at 0x10c68f610>,
   'rho': <matplotlib.figure.Figure at 0x10ac2a990>,
   'u': <matplotlib.figure.Figure at 0x10c673c50>}
  
  In [3]: figs['H'].savefig( 'enthalpy.pdf', dpi=500 )
  
  In [4]: plt.close('all')
  
"""
import numpy             as np
import matplotlib.pyplot as plt

#===============================================================================
# INPUT DATA
#===============================================================================

# Data file to be loaded
file_name = 'final.dat'

# Specific heats ratio
g = 1.4

#===============================================================================
# NUMERICAL RESULTS
#===============================================================================

# Open file
with open( file_name ) as f:
  # Read time instant
  time = float( f.readline() )
  # Load axis and conserved quantities
  x, rho, mom, eng = np.loadtxt( f, unpack=True )

# Compute primitive variables
u1 = mom/rho                   # Velocity (with sign)
p  = (g-1.0)*(eng-0.5*mom*u1)  # Pressure (>0)
H0 = (eng+p)/rho               # Total enthalpy per unit mass
c  = np.sqrt( g*p/rho )        # Speed of sound (>0)
M  = u1/c                      # Mach number (with sign)

#===============================================================================
# PLOTTING OPTIONS
#===============================================================================

# Parameters for the numerical results
params_NR = {}
params_NR['linestyle'] = 'None'
params_NR['marker'   ] = '.'
params_NR['mec'      ] = 'b'
params_NR['mfc'      ] = 'b'

# Parameters for the exact results
params_EX = {}


#===============================================================================
# PLOTS
#===============================================================================

figs = {}

#def NewFigure( x, y, time, symbol, fullname, params )
#  
#  title  = '{} at t = {}'.format( fullname, time )
#  xlabel = 'x'
#  ylabel = '{}(x)'.format( symbol )
#  
#  fig = plt.figure()
#  ax  = fig.add_subplot(1,1,1)
#  ax.set_title ( title )
#  ax.set_xlabel( xlabel )
#  ax.set_ylabel( ylabel, rotation='horizontal' )
#  ax.grid()

# Density
figs['rho'] = plt.figure()
figs['rho'].gca().plot      ( x, rho, **params_NR )
figs['rho'].gca().set_title ( 'Density at t = {}'.format(time) )
figs['rho'].gca().set_xlabel( 'x' )
figs['rho'].gca().set_ylabel( r'$\rho(x)$', rotation='horizontal' )
figs['rho'].gca().grid()

# Momentum
figs['mom'] = plt.figure()
figs['mom'].gca().plot      ( x, mom, **params_NR )
figs['mom'].gca().set_title ( 'Momentum at t = {}'.format(time) )
figs['mom'].gca().set_xlabel( 'x' )
figs['mom'].gca().set_ylabel( r'$\rho u(x)$', rotation='horizontal' )
figs['mom'].gca().grid()

# Total energy
figs['eng'] = plt.figure()
figs['eng'].gca().plot      ( x, eng, **params_NR )
figs['eng'].gca().set_title ( 'Total energy at t = {}'.format(time) )
figs['eng'].gca().set_xlabel( 'x' )
figs['eng'].gca().set_ylabel( r'$\mathcal{E}(x)$', rotation='horizontal' )
figs['eng'].gca().grid()

# Velocity
figs['u'] = plt.figure()
figs['u'].gca().plot      ( x, u1, **params_NR )
figs['u'].gca().set_title ( 'Velocity at t = {}'.format(time) )
figs['u'].gca().set_xlabel( 'x' )
figs['u'].gca().set_ylabel( 'u(x)', rotation='horizontal' )
figs['u'].gca().grid()

# Pressure
figs['p'] = plt.figure()
figs['p'].gca().plot      ( x, p, **params_NR )
figs['p'].gca().set_title ( 'Pressure at t = {}'.format(time) )
figs['p'].gca().set_xlabel( 'x' )
figs['p'].gca().set_ylabel( 'p(x)', rotation='horizontal' )
figs['p'].gca().grid()

# Enthalpy
figs['H'] = plt.figure()
figs['H'].gca().plot      ( x, H0, **params_NR )
figs['H'].gca().set_title ( 'Total specific enthalpy at t = {}'.format(time) )
figs['H'].gca().set_xlabel( 'x' )
figs['H'].gca().set_ylabel( '$h_0$(x)', rotation='horizontal' )
figs['H'].gca().grid()

# Mach number
figs['M'] = plt.figure()
figs['M'].gca().plot      ( x, M, **params_NR )
figs['M'].gca().set_title ( 'Mach number at t = {}'.format(time) )
figs['M'].gca().set_xlabel( 'x' )
figs['M'].gca().set_ylabel( 'M(x)', rotation='horizontal' )
figs['M'].gca().grid()

# Sound speed
figs['c'] = plt.figure()
figs['c'].gca().plot      ( x, c, **params_NR )
figs['c'].gca().set_title ( 'Speed of sound at t = {}'.format(time) )
figs['c'].gca().set_xlabel( 'x' )
figs['c'].gca().set_ylabel( 'c(x)', rotation='horizontal' )
figs['c'].gca().grid()


for fig in figs.itervalues():
  fig.show()

