# import all needed components
from __future__ import division
from libtbx.command_line import easy_qsub
import os

class run_queue_tests(object):
  def __init__(self):
    # set working environment
    self.phenix_source = "/net/chevy/raid1/user_name/Work_chevy/phenix_build/setpaths.csh"
    self.com_path = '/net/chevy/raid1/user_name/Work_chevy/phenix_sources/location/some_python_command.py'
    # where all queue output will be deposited
    self.where_to_run_dir = "/net/cci-filer2/raid1/home/user_name/Work/queue_job"
    self.pdb_code = []
    self.pdb_file_with_path = []
    # The commands list is a list that will be sent to the queue for processing
    self.commands = []
    # the number of command in each "chunk" sent to the queue
    self.size_of_chunks = 100

  def get_pdb_files(self):
    """() -> list

    Get a list of file names, including path, to process
    Do initial processing that returns the list of files you want to execute
    in using the queue.

    Returns:
    list of PDB file, for all files in in PDB mirror
    """
    pdb_dir = os.environ["PDB_MIRROR_PDB"]
    pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").readlines()
    for pdb_file in pdb_files:
      self.pdb_code.append(pdb_file[-12:-8])
      pdb_file = os.path.join(pdb_dir, pdb_file.strip())
      self.pdb_file_with_path.append(pdb_file)

  def get_commands(self):
    """
    get_commands process the command we want to run, with the options we want to use,
    on all the files we want the command to run on

    It produces a list of commands: list of strings,
    containing the command, options and file.
    in the same format that you would use to run is from the command prompt
    """
    # Set command path to the file/command to be execute.
    value_1 = ''  # Add relevant value
    value_2 = ''  # Add relevant value
    com_options = '--option1={0} --option2={1}'.format(value_1,value_2)

    # Build the command list
    for file_name in self.pdb_file_with_path:
      outString = '{0} {1} {2} {3} &> log_{4}'.format(
        self.com_path,
        file_name,
        com_options,
        self.pdb_code)
      self.commands.append("python {}".format(outString))


  def send_to_queue(self):
    # Send the job to the queue
    easy_qsub.run(
      phenix_source = tself.phenix_source,
      where = self.where_to_run_dir,
  #    qsub_cmd = 'qsub -q all.q@machine_name', # Optional, when you want all jobs to run on machine_name
      commands = self.commands,
      size_of_chunks= self.size_of_chunks)  # set the number of commands to send together to the queue.

if __name__ == "__main__":
  queue_job = run_queue_tests()
  queue_job.get_pdb_files()
  queue_job.get_commands()
  queue_job.send_to_queue()


