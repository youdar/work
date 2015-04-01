from __future__ import division
import collect_ncs_files
from libtbx.command_line import easy_qsub
from glob import glob
import sys
import os

class run_queue_tests(object):

  def __init__(self):
    """
    Run refinements
      '-no_ncs',
      '-cartesian_ncs_restraints',
      '-torsion_ncs_restraints'
      '-ncs_constraints_sites_no_operators',
      '-ncs_constraints_sites_operators',
      '-ncs_constraints_adp_operators',
      '-ncs_constraints_all'
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
    self.size_of_chunks = 4
    self.options = [
      '-no_ncs',
      '-cartesian_ncs_restraints',
      '-torsion_ncs_restraints'
      '-ncs_constraints_sites_no_operators',
      '-ncs_constraints_sites_operators',
      '-ncs_constraints_adp_operators',
      '-ncs_constraints_all']

  def get_pdb_files(self):
    """() -> list
    Get all pdb IDs from LBL PDB mirror index
    """
     # Run on all PDB - Only on LBL machine
    osType = sys.platform
    msg = 'Please run this only on LBL computer'
    assert not osType.startswith('win'),msg
    # set environment
    self.pdb_log_id_list = glob(self.collect_files_from + '/log_*')
    self.pdb_log_id_list = [x[-4:] for x in self.pdb_log_id_list]
    # get list of files not refined yet
    # c = collect_ncs_files.ncs_paper_data_collection()
    # fn = os.path.join(c.ncs_dir,'not_yet_refined.txt')
    # self.pdb_log_id_list = open(fn,'r').read().splitlines()
    print 'Processing {} files'.format(len(self.pdb_log_id_list))
    # for testing
    self.pdb_log_id_list = ['1vcr']

  def get_commands(self):
    """
    Build the command list

    get_commands process the command we want to run, with the options we want to use,
    on all the files we want the command to run on

    It produces a list of commands: list of strings,
    containing the command, options and file.
    in the same format that you would use to run is from the command prompt
    """
    for pdb_id in self.pdb_log_id_list:
      for option in self.options:
        outString = '{} {} {}'.format(self.com_path,pdb_id,option)
        self.commands.append("python {}".format(outString))

  def send_to_queue(self):
    # Send the job to the queue
    easy_qsub.run(
      phenix_source = self.phenix_source,
      where = self.where_to_run_dir,
      # Optional, when you want all jobs to run on machine_name
      # list of newer hosts: beck, morse, gently, rebus
      # qsub_cmd = 'qsub -q all.q@rebus',
      commands = self.commands,
      # set the number of commands to send together to the queue.
      size_of_chunks= self.size_of_chunks)

if __name__ == "__main__":
  queue_job = run_queue_tests()
  queue_job.get_pdb_files()
  queue_job.get_commands()
  queue_job.send_to_queue()


