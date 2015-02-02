from __future__ import division
import os

def run():
  files_to_read = [
    'list of pbd codes - res 3_5 to 8_5 years 2010 and on.txt',
    # 'MTRIX (gt 0) ncs groups (gt 1) ncs (ne) MTRIX'
    'pdb_157_file_list.txt'
    # 'MTRIX (eq 0) cctbx (gt 0)'
    # 'MTRIX (gt 0) cctbx (gt 0) MTRIX (ne) cctbx'
    # 'MTRIX (gt 0) cctbx (gt 0) MTRIX (eq) ccbtx'
    # 'cctbx (gt 0)'
    # 'cctbx (eq 0)'
    # 'all processedd files'
  ]

  common_files = set()
  for fn in files_to_read:
    new_files = pdb_codes_in_file(fn)
    print '{} files in {}'.format(len(new_files),fn)
    if common_files:
      common_files.intersection_update(new_files)
    else:
      common_files = new_files
    if not common_files: break

  print 'Common files:'
  print '============='
  if len(common_files) < 50:
    print list(common_files)
  print '-------------'
  print_list(common_files)
  print 'Done'

def pdb_codes_in_file(file_name):
  """ read pdb codes from file
  pdb code mast be 4 letter long """
  data = set()
  if os.path.isfile(file_name):
    data = open(file_name,'r').read().splitlines()
    # remove comment lines
    data = [x for x in data if x[0] != '#']
    # split lines, if there is more than one code per line
    data = [x for xi in data for x in xi.split(' ')]
    data = {x for x in data if len(x) == 4}
  return data

def print_list(data):
  """ print list of codes in a table """
  codes_in_row = 8
  n_codes = len(data)
  n_rows = n_codes // codes_in_row
  n_in_last_row = n_codes % codes_in_row
  # Print strings
  out_str = '{} ' * codes_in_row
  last_row_out_str = '{} ' * n_in_last_row
  data = list(data)
  for i in range(n_rows):
    print out_str.format(*data[(i * codes_in_row):((i+1) * codes_in_row)])
  if n_in_last_row > 0:
    print last_row_out_str.format(*data[-n_in_last_row:])
  print '......'


if __name__ == '__main__':
  os.chdir(r'C:\Phenix\Dev\Work\work\NCS\data')
  run()
