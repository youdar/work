from __future__ import division
from libtbx.command_line import easy_qsub
import collect_ncs_files
from glob import glob
import sys
import os

class run_queue_tests(object):

  def __init__(self):
    """
    Testing using: /work/NCS/ncs_paper/collect_ncs_files.py

    """
    # set working environment
    paths = "/net/chevy/raid1/youval/Work_chevy/build/setpaths.csh"
    self.phenix_source = paths
    sources = os.environ["workfolder"]
    # command path
    c = collect_ncs_files.ncs_paper_data_collection()
    self.com_path = sources + '/NCS/ncs_paper/get_model_vs_data.py'
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
    Get all pdb IDs from LBL PDB mirror index
    """
     # Run on all PDB - Only on LBL machine
    osType = sys.platform
    msg = 'Please run this only on LBL computer'
    assert not osType.startswith('win'),msg
    # set environment
    self.pdb_log_file_list = glob(self.collect_files_from + '/log_*')
    self.pdb_log_file_list = [x[-4:] for x in self.pdb_log_file_list]
    print 'Processing {} files'.format(len(self.pdb_log_file_list))
    # for testing
    self.pdb_log_file_list = ['2a3x']
    self.pdb_log_file_list = [
      '3hfs', '3ksb', '2p5t', '4h1l', '4gh4', '2ja7', '2hn2', '4nkz', '3km2',
      '4gx2', '3ks2', '1iss', '4h1i', '1n3n', '4gk5', '3hxk', '3g05', '3ksa',
      '4hi2', '1ott', '4kn4', '3nmj', '4hic', '2gw1', '4h37', '4gx5']


  def get_commands(self):
    """
    Build the command list

    get_commands process the command we want to run, with the options we want to use,
    on all the files we want the command to run on

    It produces a list of commands: list of strings,
    containing the command, options and file.
    in the same format that you would use to run is from the command prompt
    """
    for file_name in self.pdb_log_file_list:
      outString = '{0} {1}'.format(self.com_path,file_name)
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
  queue_job.send_to_queue()


