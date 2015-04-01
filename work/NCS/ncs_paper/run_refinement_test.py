from __future__ import division
from libtbx import easy_run
import collect_ncs_files
import shutil
import sys
import os

__author__ = 'Youval'

def run(args):
  """
  Run one of the following refinements:

  -no_ncs
  -cartesian_ncs_restraints
  -torsion_ncs_restraints
  -ncs_constraints_no_operators
  -ncs_constraints_all

  Args:
    args should contain:
      pdb_id (str)
      one of the refinements types
      -again : when this option present, erase previous refinement results

  Examples:
  >>>pyhton run_refinement_test.py 1vcr -no_ncs

  >>>pyhton run_refinement_test.py 1vcr -torsion_ncs_restraints -again
  """
  osType = sys.platform
  assert not ('win' in osType),'Please run on LBL machine'
  delete_existing_results = '-again' in args
  if delete_existing_results:
    args.pop(args.index('-again'))
  assert len(args) == 2, 'Not enough parameters provided'
  # collect files and set working directories
  pdb_id = args[0]
  refine_method = args[1]
  c = collect_ncs_files.ncs_paper_data_collection()
  pdb = os.path.join(c.asu_dir,pdb_id + '.pdb')
  mtz = os.path.join(c.mtz_dir,pdb_id + '.mtz')
  cif = os.path.join(c.cif_dir,pdb_id + '.ligands.cif')
  if not os.path.isfile(cif): cif = ''
  files = '{} {} {} '.format(pdb,mtz,cif)
  refinement_dir_list = collect_ncs_files.get_refinement_folders()
  # refinement string
  s = 'optimize_xyz=true '
  s += 'optimize_adp=true strategy=individual_sites+individual_adp'
  option_list = [
    '-no_ncs',
    '-cartesian_ncs_restraints',
    '-torsion_ncs_restraints',
    '-ncs_constraints_sites_no_operators',
    '-ncs_constraints_sites_operators',
    '-ncs_constraints_adp_operators',
    '-ncs_constraints_all']
  # remove limits on ncs distance
  dist_limit = ' refinement.ncs.excessive_distance_limit=None '
  cmd_option_list = [
    s,
    'main.ncs=True ncs.type=cartesian ' + s,
    'main.ncs=True ' + s,
    'ncs_search=true refine_operators=false strategy=individual_sites',
    'ncs_search=true refine_operators=true strategy=individual_sites',
    'ncs_search=true refine_operators=true strategy=individual_adp',
    'ncs_search=true refine_operators=true strategy=individual_adp+individual_sites']
  i = option_list.index(refine_method)
  cmd_option = cmd_option_list[i]
  out_folder = refinement_dir_list[i]
  #
  cmd = 'phenix.refine {} {} '.format(files,cmd_option) + dist_limit
  # Run refinement
  current_dir = os.getcwd()
  make_new_folder = True
  os.chdir(out_folder)
  if os.path.isdir(pdb_id):
    try:
      # remove if empty
      os.rmdir(pdb_id)
    except OSError:
      if delete_existing_results:
        shutil.rmtree(pdb_id)
      else:
        print 'refinement results already exist for:',pdb_id
        make_new_folder = False
  if make_new_folder:
    os.mkdir(pdb_id)
    os.chdir(pdb_id)
    # r = easy_run.fully_buffered(cmd)
    r = easy_run.go(cmd)
  os.chdir(current_dir)
  print 'Done'

if __name__=='__main__':
  run(sys.argv[1:])
