
#===============================================================================

def PeriodicBCs (q, mx, mbc):
  """ Impose periodic boundary conditions.
  """
  r = mbc+mx   # Index of first ghost cell on the right
  for qi in q:
    qi[:mbc] = qi[r-mbc:  r]  # Populate ghost cells on the left
    qi[r:  ] = qi[mbc:2*mbc]  # Populate ghost cells on the right

#===============================================================================

def OutflowBC_left  (q, mx, mbc):
  """ Impose outflow boundary condition at the left end.
  """
  for qi in q:
    qi[:mbc] = qi[mbc]

def OutflowBC_right (q, mx, mbc):
  """ Impose outflow boundary condition at the right end.
  """
  r = mbc+mx   # Index of first ghost cell on the right
  for qi in q:
    qi[r:] = qi[r-1]

#===============================================================================
