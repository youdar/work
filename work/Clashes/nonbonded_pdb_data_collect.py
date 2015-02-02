from __future__ import division
from mmtbx.command_line import nonbonded_overlaps
from libtbx.utils import null_out
from libtbx.utils import Sorry
from libtbx import easy_run
import iotbx.pdb
import getpass
import sys
import os

class NBO_results(object):
  """ Container for the test results """
  def __init__(self):
    """
    x1: PDB ID
    x2: Macro molecule overlaps
    x3: Symmetry overlaps
    x4: All overlaps
    x5: Macro molecule overlaps per 1000 atoms
    x6: Symmetry overlaps per 1000 atoms
    x7: All overlaps per 1000 atoms
    x8: year model deposited in PDB
    x9: experiment type
    """
    self.pdb_id = ''
    self.macro_molecule_overlaps = 0
    self.symmetry_overlaps = 0
    self.all_overlaps = 0
    self.macro_molecule_clashscore = 0.0
    self.symmetry_clashscore = 0.0
    self.all_clashscore = 0.0
    self.year = 'xxxx'
    self.experiment_type = 'xxxx'

  def __repr__(self):
    """ return a comma separated sting of data """
    s = ','.join(['{}'] * 9)
    data = [self.pdb_id,self.macro_molecule_overlaps,
            self.symmetry_overlaps,self.all_overlaps,
            round(self.macro_molecule_clashscore,1),
            round(self.symmetry_clashscore,1),round(self.all_clashscore,1),
            self.year,self.experiment_type]
    return s.format(*data)

  def __str__(self):
    """ return table representation of overlaps data """
    s = self.__repr__()
    data = s.split(',')
    header = 'PDB id| macro mol.| symmetry |  all  | year | experiment type'
    l = len(header)
    # test data length is correct
    title = header.split('|')
    data = data[:4] + data[7:]
    assert len(data) == len(title)
    # analyze the header string and create a template for the data
    # find size of each column
    s_pos = [i for i in range(l) if header[i] == '|']
    s_pos.insert(0,0)
    s_pos.append(l-1)
    s_lengths = [s_pos[i+1]-s_pos[i] for i in range(len(s_pos)-1)]
    d_length = [len(d) for d in data]
    # make sure data length is not longer than the title length
    s_lengths = [max(x,y) for x,y in zip(s_lengths,d_length)]
    # create the template
    s_template = ['{:^%d}'%x for x in s_lengths]
    s_template = '|'.join(s_template)
    # combine header and results to a table
    s_out = s_template.format(*title) + '\n'
    l = len(s_out)
    s_out += '-'*l + '\n'
    s_out += s_template.format(*data) + '\n'
    s_out = ('{:^%d}\n'%l).format('non-bonded overlaps') + s_out
    return s_out

def get_pdb_id(file_name):
  '''
  clean a pdb file name, remove path and file extensions

  Args:
    file_name (str): file name that might include path and extension

  Returns:
    file_name (str): 4 letter pdb id
  '''
  file_name = os.path.split(file_name)[-1]
  file_name = file_name.lower()
  file_name = file_name.split('.')[0]
  if len(file_name)>4:
    if 'pdb' in file_name:
      i = file_name.find('pdb')
      file_name = file_name[i+3:i+7]
  return file_name


def set_working_folder():
  """ set working folder path for Youval """
  username = getpass.getuser()
  osType = sys.platform
  if username.lower() == 'youval':
    if osType.startswith('win'):
      dr = r'C:\Phenix\Dev\Work\work\Clashes\wtest'
    else:
      dr = '/net/cci/youval/work/work/Clashes/wtest'
    os.chdir(dr)

def run(file_name,verbose=False,clean_up=True):
  """
  Update file_name to include hydrogens, then get overlaps and
  overlaps per 1000 atoms.

  Args:
    file_name (str): pdb file name including path that should be
      processed using ready_set or reduce
    verbose (bool): when True, return results in a table

  Returns:
    out_str (str):
      if verbose is False: x1,x2,x3,x,4,x5,x6,x7,x8,x9
    x1: PDB ID
    x2: Macro molecule overlaps
    x3: Symmetry overlaps
    x4: All overlaps
    x5: Macro molecule overlaps per 1000 atoms
    x6: Symmetry overlaps per 1000 atoms
    x7: All overlaps per 1000 atoms
    x8: year model deposited in PDB
    x9: experiment type

      if verbose is True, the sting will look like:

        PDB id| macro mol. | symmetry  |  all   | year  | experiment type
        -------------------------------------------------------------------
         1a18 |     29     |     5     |   35   | 1997  |X-RAY DIFFRACTION

    Error Types:
      -1: Other processing issue
      -2: model contains unknown_type_pairs
      -3: multiple models in pdb files
      -4: Bad CRYST1 records, bad crystal symmetry
      -5: File could not be fetched
  """
  # Get file and file name
  file_name,pdb_id,file_to_clean = process_file(file_name)
  nbo_out = NBO_results()
  nbo_out.pdb_id = pdb_id
  #
  if os.path.isfile(file_name):
    # get pdb date and experiment_type
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    nbo_out.year = pdb_inp.extract_header_year()
    nbo_out.experiment_type = pdb_inp.get_experiment_type()
    # get pdb overlap count and clashscore
    option = 'substitute_non_crystallographic_unit_cell_if_necessary=false'
    args = [file_name,'verbose={}'.format(verbose),option]
    try:
      r = nonbonded_overlaps.run(args,out=null_out())
      e = None
    except Sorry as e:
      r = None
    except: # catch *all* exceptions
      e = sys.exc_info()[0]
      r = None
    if e:
      # some error
      if hasattr(e,'message') and (e.message):
        unknown_type_pairs = 'unknown type pairs' in e.message
        multiple_models = 'provide only a single model' in e.message
        bad_cryst1 = 'None valid CRSYT1 records' in e.message
        if unknown_type_pairs:
          macro_molecule = sym = all = -2
        elif multiple_models:
          macro_molecule = sym = all = -3
        elif bad_cryst1:
          macro_molecule = sym = all = -4
        else:
          macro_molecule = sym = all = -1
      else:
        macro_molecule = sym = all = -1
    else:
      # All is good
      macro_molecule = r.result.nb_overlaps_macro_molecule
      sym = r.result.nb_overlaps_due_to_sym_op
      all = r.result.nb_overlaps_all
      #
      nbo_out.macro_molecule_clashscore = r.result.cctbx_clashscore_macro_molecule
      nbo_out.symmetry_clashscore = r.result.cctbx_clashscore_due_to_sym_op
      nbo_out.all_clashscore = r.result.cctbx_clashscore_all
  else:
    # file_name is not a valid pdb file
    macro_molecule = sym = all  = -5
  #
  nbo_out.macro_molecule_overlaps = macro_molecule
  nbo_out.symmetry_overlaps = sym
  nbo_out.all_overlaps = all

  if verbose:
    outstr = nbo_out.__str__()
    '''
    -1: Other processing issue
      -2: model contains unknown_type_pairs
      -3: multiple models in pdb files
      -4: Bad CRYST1 records, bad crystal symmetry
      -5: File could not be fetched
    '''
    error_dict = {
      -1: 'Other processing issue',
      -2: 'Model contains unknown_type_pairs',
      -3: 'multiple models in pdb files',
      -4: 'Bad CRYST1 records, bad crystal symmetry',
      -5: 'File could not be fetched'}
    if macro_molecule < 0:
      # if we have an error, add the error message
      outstr = error_dict[macro_molecule] + '\n' + outstr
  else:
    outstr = nbo_out.__repr__()

  if clean_up:
    # Cleanup files if they where not already in local folder
    if file_to_clean:
      for fn in file_to_clean:
        if os.path.isfile(fn): os.remove(fn)
  return outstr

def process_file(file_name):
  """ Get pdb id and pdb file """
  _,fn = os.path.split(file_name)
  pdb_id = get_pdb_id(fn)
  # Get files in pdb format even when at LBL, to avoid issues with phenix.reduce
  if ('_mirror' in file_name) or (file_name == pdb_id):
    file_name = pdb_id + '.pdb'
  file_to_clean = []
  if not os.path.isfile(file_name):
    # leave file in folder if it was already there
    # cmd = 'phenix.fetch_pdb {} --all'.format(pdb_id)
    cmd = 'phenix.fetch_pdb {}'.format(pdb_id)
    r = easy_run.go(cmd,join_stdout_stderr=False)
    for fn in r.stdout_lines:
      fn = os.path.split(fn)[-1]
      if '.pdb' in fn: file_name = fn
      file_to_clean.append(fn)
      fn = fn.replace('.pdd','_with_h.pdb')
      file_to_clean.append(fn)
  return file_name,pdb_id,file_to_clean

if (__name__ == "__main__"):
  help_str = """
Get clashscores for CCTBX and PROBE after updating the pdb file to include
hydrogens

Examples:
>>>python nonbonded_pdb_data_collect.py --verbose
PDB id| macro mol. | symmetry  |  all   | year  | experiment type
-------------------------------------------------------------------
 1a18 |     29     |     5     |   35   | 1997  |X-RAY DIFFRACTION

>>>python nonbonded_pdb_data_collect.py 1a18.pdb
1a18,29,5,35,13.9,2.4,16.5,1997,X-RAY DIFFRACTION

"""
  args = sys.argv[1:]
  if len(args) == 0:
    print help_str
  else:
    current_folder = os.getcwd()
    file_name = args[0]
    verbose = False
    if len(args) > 1:
      verbose_options = ['verbose=true','--verbose','-v']
      verbose = [x.lower() in verbose_options for x in args[1:]]
      verbose = (verbose.count(True) > 0)
    set_working_folder()
    print run(file_name,verbose=verbose,clean_up=True)
    os.chdir(current_folder)
