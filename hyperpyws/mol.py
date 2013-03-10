#coding: utf8

import numpy as np

#===============================================================================

class Slice (object):
  """
  Temporary object that permits us to apply slicing to each (array) 
  component of a numpy array of 'object' type.
  
  Parameters
  ----------
  arr : numpy.ndarray
    Multidimensional Numpy array of dtype=object
    
  """
  def __init__(self, arr):
    # Store reference to input array, and create empty array
    self._arr = arr
    self._res = np.empty( arr.shape, dtype=object )
  
  #-----------------------------------------------------------------------------
  def __getitem__(self, *c):
    """
    Pass slice to each array component of the array of objects.
    
    Parameters
    ----------
    c : list of 'slice' objects
      Slices to be applied to each component of the chosen array
    
    Returns
    -------
    res : numpy.ndarray
      Multidimensional array that contains the required slices
    
    """
    # Create flat iterators
    a = self._arr.flat
    r = self._res.flat
    # Cycle over each array component
    for i,a in enumerate(a):
      r[i] = a.__getitem__(*c)
    # Return array with slices
    return self._res

#===============================================================================

class MOL(object):
  """
  Class for 'Method Of Lines' objects, which calculate the time derivative of 
  the solution, and higher time derivatives if possible.
  
  Parameters
  ----------
  grid : Grid1D
    Object containing information about the grid, and the solution arrays
  
  flux : Flux1D
    Object calculating the flux function f(q), its Jacobian matrix J(q), 
    and its eigen-decomposition [R][D][Ri]
  
  weno : WenoReconstruction
    Object performing conservative recontruction of u[i] over stencil
    
  """
  def __init__(self, grid, flux, weno, SetBCs):
    
    # Store references
    self._grid = grid
    self._flux = flux
    self._weno = weno
    self._SetBCs = SetBCs
    
    # Pre-processing
    # (?)
  
  #-----------------------------------------------------------------------------
  def qt (self, q):
    """ Compute first time-derivative of solution vector.
    """
    # Rename variables
    dx  = self._grid.dx
    mx  = self._grid.mx
    mbc = self._grid.mbc
    meq = self._grid.meq    
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 1: Compute Roe averages
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Index with single ghost cell on left   [--> move to preprocessing?]
    Im1 = slice( mbc-1, (mbc+mx-1)+1 )
#    Im1 = range( mbc-1, (mbc+mx-1)+1 )

    # Index with single ghost cell on right  [--> move to preprocessing?]
    I   = slice( mbc  , (mbc+mx  )+1 )
#    I   = range( mbc  , (mbc+mx  )+1 )
    
    # Simple algebraic averages for now
    qs  = np.empty( meq, dtype=object )
    for m in range(meq):
      qs[m] = 0.5*( q[m][Im1] + q[m][I] )
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 2: Compute eigen-decomposition of Jacobian (projection matrices)
    #         and maximum wave speed in domain, using averages
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    R     = self._flux.R (qs)
    L     = self._flux.L (qs)
    alpha = max([ max(abs(a)) for a in self._flux.eig(qs) ])
    
#    print "File 'mol.py':"
#    print alpha # <<<<<<<<<<<<<<< DEBUG
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 3: Project q[i+s] and f[i+s] over local characteristic variables at 
    #         location x[i], for each point x[i+s] in stencil
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # WENO5 reconstructs f[i-1/2] using stencil [-3,-2,-1, 0,+1,+2]
    f  = self._flux.f(q)
    
    # Stencils  [--> move to preprocessing?]
    extended_stencil = [ self._weno.stencil[0]-1 ] + self._weno.stencil
    S  = [ slice( mbc+sh, (mbc+mx+sh)+1 ) for sh in extended_stencil ]
#    S  = [ range( mbc+sh, (mbc+mx+sh)+1 ) for sh in extended_stencil ]
    
    # Shift q and f according to stencil
    qq = [ Slice(q)[Si] for Si in S ]  # for now, slice on 2D array
    ff = [ Slice(f)[Si] for Si in S ]  # for now, slice on 2D array
    
    # Project onto characteristic variables: q -> w, f -> g
    ww = [ np.dot(L,qi) for qi in qq ]
    gg = [ np.dot(L,fi) for fi in ff ]
    
    
#    print "File 'mol.py':"
#    print ff[3] # <<<<<<<<<<<<<<< DEBUG
#    print gg[3] # <<<<<<<<<<<<<<< DEBUG 
#    print ''
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 4: Flux-splitting and WENO reconstruction
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    ghat = np.empty( meq, dtype=object )
    for m in range(meq):
      
      # TODO: can add in an alpha[me] as well if we choose ...
      
      # Flux splitting
      gg_m = [ 0.5*gi[m] - alpha * wi[m] for gi,wi in zip(gg[1 :],ww[1 :]) ]
      gg_p = [ 0.5*gi[m] + alpha * wi[m] for gi,wi in zip(gg[:-1],ww[:-1]) ]
      
      # Weno reconstruction
      gp = self._weno.reconstruct_left  (*gg_p)
      gm = self._weno.reconstruct_right (*gg_m)
      
      # Sum right and left fluxes
      ghat[m] = gp + gm
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 5: Project flux values back onto conserved variables
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    fhat = np.dot(R,ghat)
    
#    print "File 'mol.py':"
#    print fhat # <<<<<<<<<<<<<<< DEBUG
#    print ghat # <<<<<<<<<<<<<<< DEBUG 
#    print ''
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 6: Compute d/dt(q[i]) according to conservative finite-differences
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    c   = -1.0/dx
    q_t =  0.0*q
    for m in range(meq):
      q_t[m][mbc:mbc+mx] = c * ( fhat[m][1:] - fhat[m][:-1] )
    
    return q_t
  
  #-----------------------------------------------------------------------------
  def qtt (self, q, q_t):
    """
    Compute second time-derivative of solution vector, as
    q_{tt} = -[ f(q)_t ]_x = -[ f'(q) q_t ]_x.
    
    Using centered finite differences we get
    -1/12 *( f_t[i-2] - 8 f_t[i-1] + 8 f_t[i+1] - f_t[i+2] ),
    
    where f_t[i] = f'[i] * q_t[i].  q_t[i] is computed by the previous function.
    
    """
    # Q0: q_t as input parameter?
    # Q1: Set boundary on q_t ??
    # Q2: How to get appropriate FD approximation for given order of accuracy?
    
    # Rename variables
    dx  = self._grid.dx
    mx  = self._grid.mx
    mbc = self._grid.mbc
    meq = self._grid.meq
    
    # Stencils  [--> move to preprocessing?]    
    a,b = mbc, mbc+mx
    Im2 = slice( a-2, b-2 )
    Im1 = slice( a-1, b-1 )
    I   = slice( a  , b   )
    Ip1 = slice( a+1, b+1 )
    Ip2 = slice( a+2, b+2 )
#    I   = np.arange(mbc, mbc+mx)
#    Im2 = I-2
#    Im1 = I-1
#    Ip1 = I+1
#    Ip2 = I+2
    
    # Compute time derivative of flux function
    ft = np.dot( self._flux.J(q), q_t )
    
    # Compute 2nd time derivative of state vector by differentiating ft in space
    c   = -1.0/(12.*dx)
    q_tt = 0.0*q_t
    for m in range(meq):
      q_tt[m][I] = c *( ft[m][Im2] - 8.*(ft[m][Im1] - ft[m][Ip1]) - ft[m][Ip2] )
    
    return q_tt
  
  #-----------------------------------------------------------------------------
  def TimeDerivatives (self, q, t):
    """
    Compute 1st and 2nd time-derivatives of solution vector:
      1. Apply boundary conditions to solution (may depend on t);
      2. Compute qt (q);
      3. Compute qtt(q,qt);
      4. Return [ qt, qtt].
    
    Parameters
    ----------
    q : array-like
      Solution vector at time t
    t : float
      Time instant
    
    Returns
    -------
    r : list
      [qt,qtt] - 1st and 2nd time-derivatives of solution vector q
    
    """
    self._SetBCs(q,t)  # Apply boundary conditions 
    q_t  = self.qt (q)       # Compute 1st time derivative
    q_tt = self.qtt(q, q_t)  # Compute 2nd time derivative
    
    return [q_t, q_tt]
  
#===============================================================================
