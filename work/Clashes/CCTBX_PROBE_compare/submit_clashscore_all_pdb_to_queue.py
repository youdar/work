from __future__ import division
from compare_clashscores_CCTBX_vs_PROBE import get_pdb_id
from libtbx.command_line import easy_qsub
import os

"""
Perform a clashscore comparison test CCTBX vs PROBE for all PDB files in the
LBL PDB mirror
"""

class run_queue_tests(object):

  def __init__(self, clean_source=False,size_of_chunks=1):
    """ Test using compare_clashscores_CCTBX_vs_PROBE.py """
    # set working environment To "chevy" or "clean" phenix copy
    if clean_source:
      paths = "/net/chevy/raid1/youval/Clean_copy/build/setpaths.csh"
    else:
      paths = "/net/chevy/raid1/youval/Work_chevy/build/setpaths.csh"
    self.phenix_source = paths
    command_path = '/net/cci/youval/work/work/Clashes'
    command = 'compare_CCTBX_vs_PROBE_clashscores.py'
    self.com_path = os.path.join(command_path,command)
    # where all queue output will be deposited
    queue_job = "/net/cci/youval/work/work/Clashes/queue_clash_compare"
    self.where_to_run_dir = queue_job
    self.pdb_dir = os.environ["PDB_MIRROR_PDB"]
    # list of PDB id codes
    self.pdb_code = []
    # file name and path for each pdb structure
    self.pdb_file_with_path = []
    # The commands list is a list that will be sent to the queue for processing
    self.commands = []
    # the number of command in each "chunk" sent to the queue
    self.size_of_chunks = size_of_chunks

  def get_pdb_files(self,use_files_to_run=False):
    """ Get all pdb files form LBL mirror

    Args:
      use_files_to_run (bool): when True get file names fron files_to_run.txt
    """
    fn = '/net/cci/youval/work/work/Clashes/files_to_run.txt'
    if use_files_to_run and os.path.isfile(fn):
      files = open(fn, "r").read().splitlines()
    else:
      files = open(os.path.join(self.pdb_dir, "INDEX"), "r").read().splitlines()
    self.pdb_file_with_path = files

  def get_commands(self):
    """
    get_commands process the command we want to run, with the options we want to use,
    on all the files we want the command to run on

    It produces a list of commands: list of strings,
    containing the command, options and file.
    in the same format that you would use to run is from the command prompt
    """
    # Build the command list
    for fn in self.pdb_file_with_path:
      pdb_id = get_pdb_id(fn)
      # fn = os.path.join(self.pdb_dir,fn)
      # "-c" option leaves only the macro molecule
      outString = '{0} {1} -c >& log_{1}'.format(self.com_path,pdb_id)
      self.commands.append("python {}".format(outString))

  def send_to_queue(self):
    """ Send the job to the queue """""
    easy_qsub.run(
      phenix_source = self.phenix_source,
      where = self.where_to_run_dir,
      # Optional, when you want all jobs to run on machine_name
      # list of newer hosts: beck, morse, gently, rebus
      qsub_cmd = 'qsub -q all.q@rebus',
      commands = self.commands,
      # set the number of commands to send together to the queue.
      size_of_chunks= self.size_of_chunks)


if (__name__ == "__main__"):
  queue_job = run_queue_tests(size_of_chunks=10)
  # queue_job.get_pdb_files(use_files_to_run=True)
  queue_job.get_pdb_files(use_files_to_run=False)
  queue_job.get_commands()
  queue_job.send_to_queue()
  print 'Done'
