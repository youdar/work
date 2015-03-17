from __future__ import division
from libtbx import easy_run
import collect_ncs_files
from glob import glob
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

  """
  osType = sys.platform
  assert not ('win' in osType),'Please run on LBL machine'
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
  refinement_dir_list = [
    'refine_no_ncs_dir','refine_cartesian_ncs','refine_torsion_ncs',
    'refine_ncs_con_no_oper','refine_ncs_con_all']
  option_list = [
    '-no_ncs','-cartesian_ncs_restraints','-torsion_ncs_restraints',
    '-ncs_constraints_no_operators','-ncs_constraints_all']
  cmd_option_list = [
    '','main.ncs=True ncs.type=cartesian','main.ncs=True',
    '','']
  i = option_list.index(refine_method)
  cmd_option = cmd_option_list[i]
  out_folder = c.__dict__[refinement_dir_list[i]]
  # refinement string
  s = 'optimize_xyz=true '
  s += 'optimize_adp=true strategy=individual_sites+individual_adp'
  cmd = 'phenix.refine {} {} '.format(files,cmd_option) + s

  # Run refinement
  current_dir = os.getcwd()
  make_new_folder = True
  os.chdir(out_folder)
  if os.path.isdir(pdb_id):
    try:
      os.rmdir(pdb_id)
    except OSError:
      print 'refinement results already exist for:',pdb_id
      make_new_folder = False
  if make_new_folder:
    os.mkdir(pdb_id)
    os.chdir(pdb_id)
    r = easy_run.fully_buffered(cmd)
  os.chdir(current_dir)
  print 'Done'

if __name__=='__main__':
  run(sys.argv[1:])
