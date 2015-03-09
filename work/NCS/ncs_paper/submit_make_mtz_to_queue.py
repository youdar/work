from __future__ import division
from collect_ncs_files import get_4_letters_pdb_id
from libtbx.command_line import easy_qsub
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
    self.com_path = sources + '/NCS/ncs_paper/get_mtz.py'
    # where all queue output will be deposited
    self.where_to_run_dir = sources + '/NCS/ncs_paper/ncs_queue_results'
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
    pdb_dir = os.environ["PDB_MIRROR_PDB"]
    pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
    self.pdb_code_list = [get_4_letters_pdb_id(x) for x in pdb_files]
    print 'Processing {} files'.format(len(self.pdb_code_list))
    # for testing
    self.pdb_code_list = ['1vcr']

  def get_commands(self):
    """
    Build the command list

    get_commands process the command we want to run, with the options we want to use,
    on all the files we want the command to run on

    It produces a list of commands: list of strings,
    containing the command, options and file.
    in the same format that you would use to run is from the command prompt
    """
    for file_name in self.pdb_code_list:
      outString = '{0} {1}'.format(self.com_path,file_name)
      self.commands.append("python {}".format(outString))

  def send_to_queue(self):
    # Send the job to the queue
    easy_qsub.run(
      phenix_source = self.phenix_source,
      where = self.where_to_run_dir,
      # Optional, when you want all jobs to run on machine_name
      # list of newer hosts: beck, morse, gently, rebus
      qsub_cmd = 'qsub -q all.q@beck',
      commands = self.commands,
      # set the number of commands to send together to the queue.
      size_of_chunks= self.size_of_chunks)

if __name__ == "__main__":
  queue_job = run_queue_tests()
  queue_job.get_pdb_files()
  queue_job.get_commands()
  queue_job.send_to_queue()


