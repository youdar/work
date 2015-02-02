from __future__ import division
from libtbx.command_line import easy_qsub
import os

class run_queue_tests(object):

  def __init__(self):
    """
    Testing using: /misc_scripts/strict_ncs/test_ncs_refinement.py
    """
    # set working environment
    paths = "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
    self.phenix_source = paths
    sources = os.environ["sources_chevy"]
    self.com_path = sources + '/misc_scripts/strict_ncs/test_ncs_refinement.py'
    # where all queue output will be deposited
    queue_job = "/net/cci/youval/Work/work/NCS/junk/pdb_test/queue_job"
    self.where_to_run_dir = queue_job
    self.pdb_code = []
    self.pdb_file_with_path = []
    # The commands list is a list that will be sent to the queue for processing
    self.commands = []
    # the number of command in each "chunk" sent to the queue
    self.size_of_chunks = 1

  def get_pdb_files(self):
    """() -> list

    use files from Jan.1.2010 to Nov.18.2014 with resolution of 3.5 - 9.5
    that have NCS relations
    """
    path = '/net/cci-filer2/raid1/home/youval/Work/work/NCS/refinement_test/refinment_subeset'
    file_name = 'jan_1_2010 to nov_18_2014 res 3_5 to 9_5 with NCS'
    fn = os.path.join(path,file_name)
    self.pdb_code_list = open(fn,'r').read().splitlines()
    # for testing
    # self.pdb_code_list = ['4NTT']
    # redo files
    self.pdb_code_list = [
      '1pgl', '1f8v', '1zba', '4aed', '1pgw', '3raa', '2vq0', '4hl8', '2wff',
      '4jgz', '4jgy', '2c4q', '1laj', '1za7', '1b35', '2c50', '3e8k', '3zfe',
      '4iv1', '2ws9', '1lp3', '3bcc', '3vbs', '1tdi', '3oah', '1f2n', '2bs1',
      '1k5m', '3ntt', '1uf2', '2gh8', '2bu1', '1ei7', '3p0s', '2vf9', '1r2j',
      '3vbh', '1tnv', '3vbo', '4ftb', '3s6p', '3vbf', '3dar', '4fts', '2qij',
      '2c4y', '2c4z', '2xpj', '2e0z', '2g34', '2buk', '2g33', '1w39', '2x5i',
      '3fbm', '1js9', '2qqp', '2ztn', '4gh4', '3zfg', '5msf', '4gb3', '1vcr',
      '1llc', '1x9t', '1ny7', '1x9p', '2gtl', '2qzv', '3cji', '1m1c', '1dwn',
      '1ddl', '2wbh', '3ux1', '4gbt', '3vdd', '2bfu', '3r0r', '2vf1', '3nop',
      '1h8t', '1c8n', '3not', '3nou', '1dzl', '1bcc', '6msf', '1ng0', '1vak',
      '3tn9', '3hag', '3n7x', '1wce', '1wcd', '1vb2', '2bq5', '1vb4', '2zzq',
      '1a37', '1a34', '3es5', '2iz9', '2iz8', '4f5x', '3chx', '4g0r', '2w4y',
      '2w4z', '2c51', '3lob', '2fz2', '4iv3', '2fz1', '3ra4', '2wws', '1z7s',
      '1qju', '3ra2', '1qjx', '1qjy', '3zff', '3ra8', '3ra9', '7msf', '1x36',
      '2izw', '1x35', '1vsz', '3kz4', '1ohf', '1ohg', '4bcu', '2wzr', '2izn']

  def get_commands(self):
    """
    get_commands process the command we want to run, with the options we want to use,
    on all the files we want the command to run on

    It produces a list of commands: list of strings,
    containing the command, options and file.
    in the same format that you would use to run is from the command prompt
    """
    # Build the command list
    for file_name in self.pdb_code_list:
      outString = '{0} {1} >& log_{1}'.format(self.com_path,file_name)
      self.commands.append("python {}".format(outString))

  def send_to_queue(self):
    # Send the job to the queue
    easy_qsub.run(
      phenix_source = self.phenix_source,
      where = self.where_to_run_dir,
      # Optional, when you want all jobs to run on machine_name
      # list of newer hosts: beck, morse, gently, rebus
      #qsub_cmd = 'qsub -q all.q@beck',
      commands = self.commands,
      # set the number of commands to send together to the queue.
      size_of_chunks= self.size_of_chunks)

if __name__ == "__main__":
  queue_job = run_queue_tests()
  queue_job.get_pdb_files()
  queue_job.get_commands()
  # queue_job.send_to_queue()


