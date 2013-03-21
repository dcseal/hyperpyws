"""
Postprocessing script for Buckley-Leverett's equations. 
Loads numerical data from output files and create unique figure for comparison.
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

# Name of output figure
OUTFIG = None

#===============================================================================
# LOAD DATA FROM FILES
#===============================================================================

exact  = { 'file' : 'exact.dat'       }
WENO_c = { 'file' : 'weno_coarse.dat' }
WENO_f = { 'file' : 'weno_fine.dat'   }
DG_c   = { 'file' : 'DG_coarse.dat'   }
DG_f   = { 'file' : 'DG_fine.dat'     }

datasets = [ exact, WENO_c, WENO_f, DG_c, DG_f ]

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
params_EX['color'    ] = 'r'
params_EX['linestyle'] = '-'
params_EX['linewidth'] = 2.0

# Parameters for the numerical results on a coarse mesh
params_COARSE = {}
params_COARSE['linestyle'] = 'None'
params_COARSE['marker'   ] = 'o'
params_COARSE['mec'      ] = 'b'
params_COARSE['mfc'      ] = 'None'

# Parameters for the numerical results on a coarse mesh
params_FINE = {}
params_FINE['linestyle' ] = 'None'
params_FINE['marker'    ] = '.'
params_FINE['mec'       ] = 'k'
params_FINE['mfc'       ] = 'k'
params_FINE['markersize'] = 3.0

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
fig = plt.figure()

# [1] WENO, full view
ax1 = fig.add_subplot(2,2,1)
ax1.set_title( 'WENO solution at t = {}'.format( t ) )
ax1.plot     ( exact ['x'], exact ['q'], label='Exact' , **params_EX     )
ax1.plot     ( WENO_c['x'], WENO_c['q'], label='mx=100', **params_COARSE )
ax1.plot     ( WENO_f['x'], WENO_f['q'], label='mx=600', **params_FINE   )
AddDefaults  ( ax1 )

# [2] WENO, zoom in
ax2 = fig.add_subplot(2,2,2)
ax2.set_title( 'Zoom (WENO)' )
ax2.plot     ( exact ['x'], exact ['q'], label='Exact' , **params_EX     )
ax2.plot     ( WENO_c['x'], WENO_c['q'], label='mx=100', **params_COARSE )
ax2.plot     ( WENO_f['x'], WENO_f['q'], label='mx=600', **params_FINE   )
AddDefaults  ( ax2 );  Zoom( ax2 )

# [3] DG, full view
ax3 = fig.add_subplot(2,2,3)
ax3.set_title( 'DG solution at t = {}'.format( t ) )
ax3.plot     ( exact['x'], exact['q'], label='Exact' , **params_EX     )
ax3.plot     ( DG_c ['x'], DG_c ['q'], label='mx=25' , **params_COARSE )
ax3.plot     ( DG_f ['x'], DG_f ['q'], label='mx=100', **params_FINE   )
AddDefaults  ( ax3 )

# [4] DG, full view
ax4 = fig.add_subplot(2,2,4)
ax4.set_title( 'Zoom (DG)' )
ax4.plot     ( exact['x'], exact['q'], label='Exact' , **params_EX     )
ax4.plot     ( DG_c ['x'], DG_c ['q'], label='mx=25' , **params_COARSE )
ax4.plot     ( DG_f ['x'], DG_f ['q'], label='mx=100', **params_FINE   )
AddDefaults  ( ax4 );  Zoom( ax4 )

# Show figure
plt.tight_layout()
fig.show()
