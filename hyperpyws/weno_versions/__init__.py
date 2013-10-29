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

def Weno( order, version ):
  """
  Selector function for the 'weno_versions' module: selects and returns 
  the appropriate class from the correct sub-module.
  
  Parameters
  ----------
  order : int
    Nominal order of reconstruction procedure (only valid for smooth profiles)
  
  version : str
    WENO reconstruction family 
    (different families differ for non-linear weights and smoothness parameters)
  
  Returns
  -------
  cls : type
    Abstract class for performing left and right WENO reconstructions 
    (essentially a container of functions)
  
  """
  # Available WENO order and version
  avail_order   = (5,7) 
  avail_version = ('JS','Z','CFD')
  
  # Check input: order
  if type( order ) is not int:
    raise TypeError( "'order' must be integer" )
  if order not in avail_order:
    raise ValueError( 'Requested order {:d} not available'.format( order ) )
  
  # Check input: version
  if type( version ) is not str:
    raise TypeError( "'version' must be string" )
  if version not in avail_version:
    raise ValueError( 'Requested family {:s} not available'.format( version ) )
  
  # Import correct module
  if   version == 'JS' :  from . import weno_js as mod
  elif version == 'Z'  :  from . import weno_z  as mod
  elif version == 'CFD':  from . import cfd     as mod
  
  # Determine class name
  if   version == 'JS' :  cls_name = 'Weno{:d}_JS'.format( order )
  elif version == 'Z'  :  cls_name = 'Weno{:d}_Z' .format( order )
  elif version == 'CFD':  cls_name = 'CentralFiniteDifference{:d}'.format(order)
  
  # Import class from module
  try:
    cls = getattr( mod, cls_name )
  except AttributeError:
    raise NotImplementedError, "Class '{}' not available".format( cls_name )
  
  # Return class handler
  return cls
