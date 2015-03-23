from __future__ import division
import collect_ncs_files
from libtbx.command_line import easy_qsub
from glob import glob
import shutil
import sys
import os

class run_queue_tests(object):

  def __init__(self):
    """
    look for refinements that did not work and redo them

    """
    # set working environment
    paths = "/net/chevy/raid1/youval/Work_chevy/build/setpaths.csh"
    self.phenix_source = paths
    sources = os.environ["workfolder"]
    # command path
    c = collect_ncs_files.ncs_paper_data_collection()
    self.com_path = sources + '/NCS/ncs_paper/run_refinement_test.py'
    # where all queue output will be deposited
    self.where_to_run_dir = sources + '/NCS/ncs_paper/ncs_queue_results'
    self.collect_files_from = c.data_dir
    self.pdb_code = []
    self.pdb_file_with_path = []
    # The commands list is a list that will be sent to the queue for processing
    self.commands = []
    # the number of command in each "chunk" sent to the queue
    self.size_of_chunks = 200


  def get_pdb_files(self):
    """() -> list
    Get a list of strings 'pdb_id  refine type' that failed
    """
     # Run on all PDB - Only on LBL machine
    osType = sys.platform
    msg = 'Please run this only on LBL computer'
    assert not osType.startswith('win'),msg
    self.test_to_run = []
    # the order of refinement_dir_list must be the same as refinement_types
    refinement_dir_list = [
      'refine_no_ncs_dir','refine_cartesian_ncs','refine_torsion_ncs',
      'refine_ncs_con_no_oper','refine_ncs_con_all']
    options = [
      '-no_ncs','-cartesian_ncs_restraints','-torsion_ncs_restraints',
      '-ncs_constraints_no_operators','-ncs_constraints_all']
    c = collect_ncs_files.ncs_paper_data_collection()
    refinement_types = collect_ncs_files.get_refine_test_names()
    for rdn,rt,op in zip(refinement_dir_list,refinement_types,options):
      out_folder = c.__dict__[rdn]
      pdb_id_dirs = glob(os.path.join(out_folder,'*'))
      for pdb_dir in pdb_id_dirs:
        pdb_id = pdb_dir[-4:]
        pdb_info_fn = os.path.join(c.data_dir,'log_' + pdb_id)
        pdb_info = c.read_from_file(pdb_info_fn)
        if pdb_info.refinement_records.has_key(rt):
          if pdb_info.refinement_records[rt].r_free_final == 0:
            # del: the line below will prevent redoing when -no_ncs has no value
            if pdb_info.refinement_records.has_key('no ncs'):
              if pdb_info.refinement_records['no ncs'].r_free_final > 0:
                # no result - rerun refinement
                self.test_to_run.append(pdb_id + ' ' + op)
                # delete existing folder
                print pdb_dir
                # shutil.rmtree(pdb_dir)

    print 'Processing {} files'.format(len(self.test_to_run))
    # for testing
    self.test_to_run = ['2aji -cartesian_ncs_restraints']

  def get_commands(self):
    """
    Build the command list

    get_commands process the command we want to run, with the options we want to use,
    on all the files we want the command to run on

    It produces a list of commands: list of strings,
    containing the command, options and file.
    in the same format that you would use to run is from the command prompt
    """
    for redo_test in self.test_to_run:
      outString = '{} {}'.format(self.com_path,redo_test)
      self.commands.append("python {}".format(outString))

  def send_to_queue(self):
    # Send the job to the queue
    easy_qsub.run(
      phenix_source = self.phenix_source,
      where = self.where_to_run_dir,
      # Optional, when you want all jobs to run on machine_name
      # list of newer hosts: beck, morse, gently, rebus
      # qsub_cmd = 'qsub -q all.q@beck',
      commands = self.commands,
      # set the number of commands to send together to the queue.
      size_of_chunks= self.size_of_chunks)

if __name__ == "__main__":
  queue_job = run_queue_tests()
  queue_job.get_pdb_files()
  queue_job.get_commands()
  # queue_job.send_to_queue()


