#===============================================================================
# CLASS: Plain text database
#===============================================================================

class TextDB (object):
  """
  Plain text database, with data added line by line according to a user-defined 
  format. Name, description and format of each column is prescribed by using the
  SetField() method. When the new file is open, a header portion is added at the
  beginning of the file.
  
  Parameters
  ----------
  name : str
    Full file name.
  
  """
  def __init__(self,name):
    assert (type(name) is str)
    self._name     = name
    self._fields   = []
    self._ostream  = None
    self._template = None
  
  #-----------------------------------------------------------------------------
  def SetField (self, name, frmt, numtype=float, desc=None):
    """
    Add one field (i.e. one column) to the text database.
    
    Parameters
    ----------
    name : str
      Variable name (useful for debugging and extracting data later on).
    frmt : str
      Format (e.g. '+2.5e') for printing quantity to file.
    numtype : type
      Type of data.
    desc : str
      Description to be included in header section.
    
    """
    assert (type(name) is str)
    assert (type(frmt) is str)
    assert (type(desc) is str)  
    field = {'name': name, 'frmt': frmt, 'type': numtype, 'desc': desc }
    self._fields.append(field)
  
  #-----------------------------------------------------------------------------
  def open (self):
    """
    Open database stream. If opening for the first time, a new plain text file 
    is created with a header section at its beginning. If reopening, the 
    previous file is opened in 'append' mode.
    
    """
    # If opening for the first time, create a new file
    if self._ostream is None:
      self._ostream = open(self._name,'w')
      self._ostream.write (self._header())
      self._template = \
      ' '.join(['{:%s}'% field['frmt'] for field in self._fields]) + '\n'
    # If reopening, append to old file 
    else:
      self._ostream = open(self._ostream.name,'a')
  
  #-----------------------------------------------------------------------------
  def write (self,*data):
    """
    Add a new line to the database, according to the formatting defined in a 
    stored template.  The numerical values are passed as an argument list.
    
    """
    self._ostream.write(self._template.format(*data))
  
  #-----------------------------------------------------------------------------
  def close (self):
    """ Close database file stream.
    """
    self._ostream.close()
  
  #-----------------------------------------------------------------------------
  def _header (self):
    header = '# {:s}\n'.format(self._name)
    for idx,field in enumerate(self._fields):
      l0 = '#'
      l1 = '# Column {:2d}'
      l2 = '# ---------'
      l3 = '# {:s} {:s}'
      l4 = '# {:s}'
      lines = '\n'.join([l0,l1,l2,l3,l4]) + '\n'
      header += lines.format(idx,field['name'],field['type'],field['desc'])
    header += '\n'
    return header
