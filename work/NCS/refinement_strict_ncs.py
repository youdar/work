from __future__ import division
import phenix.refinement
import phenix.refinement.driver
import sys,os

def run(file_to_refine,mtz_file,number_of_macro_cycles):
  from phenix.refinement import command_line

  args = [file_to_refine,
          mtz_file,
          'strategy=individual_sites',
          'main.number_of_macro_cycles={}'.format(number_of_macro_cycles),
          'output.prefix=refine_strictNCS_output',
          '--overwrite', '--quiet']

  command_line.run(command_name="phenix.refine", args=args)


if __name__=='__main__':
  osType = sys.platform
  if osType.startswith('win'):
    tempdir = (r'C:\Phenix\Dev\Work\work\NCS\junk')
  else:
    tempdir = ('/net/cci/youval/Work/work/NCS/junk')
  os.chdir(tempdir)
  print 'start'
  file_to_refine = 'ncs1_shaken.pdb'
  mtz_file = 'asu0_map.mtz'
  number_of_macro_cycles = 1
  run(file_to_refine,mtz_file,number_of_macro_cycles)
  print 'Done'
