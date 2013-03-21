"""
Usage
=====

In [1]: run compare_integrators_buckley_py DG
In [2]: fig.savefig( 'BuckleyLeverett_DG.pdf' )

In [3]: run compare_integrators_buckley_py WENO
In [4]: fig.savefig( 'BuckleyLeverett_WENO.pdf' )

"""
import sys, os, numpy as np
import matplotlib.pyplot as plt

#===============================================================================
# INPUT DATA
#===============================================================================
origin = os.path.abspath( os.path.curdir )

os.chdir( sys.argv[1] )
execfile( 'README.py' )

#===============================================================================
# LOAD DATA FROM FILES
#===============================================================================

exact = { 'file' : 'exact.dat' }

Stepper1_c = { 'file' : 'integrator1_coarse.dat' }
Stepper2_c = { 'file' : 'integrator2_coarse.dat' }
Stepper3_c = { 'file' : 'integrator3_coarse.dat' }

Stepper1_f = { 'file' : 'integrator1_fine.dat' }
Stepper2_f = { 'file' : 'integrator2_fine.dat' }
Stepper3_f = { 'file' : 'integrator3_fine.dat' }

datasets = [ exact, Stepper1_c, Stepper2_c, Stepper3_c,
                    Stepper1_f, Stepper2_f, Stepper3_f ]

for ds in datasets:
  # Open file
  with open( ds['file'] ) as f:
    # Read time instant
    ds['time'] = float( f.readline() )
    # Load x and q(x)
    ds['x'], ds['q'] = np.loadtxt( f, unpack=True )

#===============================================================================
# PLOTTING OPTIONS
#===============================================================================

# Parameters for the exact results
params_EX = {}
params_EX['label'    ] = 'Exact'
params_EX['color'    ] = 'r'
params_EX['linestyle'] = '-'
params_EX['linewidth'] = 2.0

# Integrator 1
params_1 = {}
params_1['label'    ] = Integrators[1]
params_1['linestyle'] = 'None'
params_1['marker'   ] = 'o'
params_1['mec'      ] = 'b'
params_1['mfc'      ] = 'None'

# Integrator 2
params_2 = {}
params_2['label'    ] = Integrators[2]
params_2['linestyle' ] = 'None'
params_2['marker'    ] = '+'
params_2['mec'       ] = 'k'
params_2['mfc'       ] = 'k'
params_2['markersize'] = 8.0

# Integrator 3
params_3 = {}
params_3['label'    ] = Integrators[3]
params_3['color'    ] = 'g'
params_3['linestyle' ] = '--'
params_3['marker'    ] = '.'
params_3['markersize'] = 3.0

#===============================================================================
# PLOTS
#===============================================================================


def AddDefaults( ax ):
  ax.set_xlabel( 'x' )
  ax.set_ylabel( 'q(x)', rotation='horizontal' )
  ax.set_ylim  ([-0.05, 1.05])
  ax.legend()
  ax.grid()

def Zoom( ax ):
  ax.set_xlim([-0.10, 0.10])
  ax.set_ylim([ 0.85, 1.01])

#def Zoom( ax ):
#  ax.set_xlim([-0.12,0.65])
#  ax.set_ylim([ 0.43, 1.03])


t = exact['time']
#(t,x,y) = exact['time'], exact['x'], exact['q']

# Create figure
fig = plt.figure( figsize=[15, 10] )

# [1] Coarse mesh, full view
ax1 = fig.add_subplot(2,2,1)
ax1.set_title( '{} solution (mx={})'.format( Method, mx_coarse ) )
ax1.plot( exact     ['x'], exact     ['q'], **params_EX )
ax1.plot( Stepper1_c['x'], Stepper1_c['q'], **params_1  )
ax1.plot( Stepper2_c['x'], Stepper2_c['q'], **params_2  )
ax1.plot( Stepper3_c['x'], Stepper3_c['q'], **params_3  )
AddDefaults  ( ax1 )

# [2] Coarse mesh, zoom in
ax2 = fig.add_subplot(2,2,2)
ax2.set_title( '{} solution (mx={})'.format( Method, mx_coarse ) )
ax2.plot( exact     ['x'], exact     ['q'], **params_EX )
ax2.plot( Stepper1_c['x'], Stepper1_c['q'], **params_1  )
ax2.plot( Stepper2_c['x'], Stepper2_c['q'], **params_2  )
ax2.plot( Stepper3_c['x'], Stepper3_c['q'], **params_3  )
AddDefaults  ( ax2 )
ax2.set_xlim([-0.10, 0.10])
ax2.set_ylim([ 0.85, 1.01])

# [3] Fine mesh, full view
ax3 = fig.add_subplot(2,2,3)
ax3.set_title( '{} solution (mx={})'.format( Method, mx_fine ) )
ax3.plot( exact     ['x'], exact     ['q'], **params_EX )
ax3.plot( Stepper1_f['x'], Stepper1_f['q'], **params_1  )
ax3.plot( Stepper2_f['x'], Stepper2_f['q'], **params_2  )
ax3.plot( Stepper3_f['x'], Stepper3_f['q'], **params_3  )
AddDefaults  ( ax3 )

# [4] Fine mesh, zoom in
ax4 = fig.add_subplot(2,2,4)
ax4.set_title( '{} solution (mx={})'.format( Method, mx_fine ) )
ax4.plot( exact     ['x'], exact     ['q'], **params_EX )
ax4.plot( Stepper1_f['x'], Stepper1_f['q'], **params_1  )
ax4.plot( Stepper2_f['x'], Stepper2_f['q'], **params_2  )
ax4.plot( Stepper3_f['x'], Stepper3_f['q'], **params_3  )
AddDefaults  ( ax4 )
#ax4.set_xlim([-0.09, 0.03])
#ax4.set_ylim([ 0.98, 1.01])
ax4.set_xlim([-0.07, 0.005])
ax4.set_ylim([ 0.987, 1.003])

# Show figure
plt.tight_layout()
fig.show()

#===============================================================================
os.chdir( origin )
