from __future__ import division
from libtbx.command_line import easy_qsub
import os

class run_queue_tests(object):
  def __init__(self):
    # set working environment
    self.phenix_source = \
      "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
    self.com_path = \
      '/net/cci/youval/Work/work/NCS/test_ncs_refinement.py'
    # where all queue output will be deposited
    self.where_to_run_dir = \
      "/net/cci/youval/Work/work/NCS/junk/pdb_test/queue_job"
    self.pdb_code = []
    self.pdb_file_with_path = []
    # The commands list is a list that will be sent to the queue for processing
    self.commands = []
    # the number of command in each "chunk" sent to the queue
    self.size_of_chunks = 10

  def get_pdb_files(self):
    """() -> list

    the list of 147 file out of 159 that have single NCS with
    MTRIX and CIF files that can be processed

    The list can be found at:
    C:\Phenix\Dev\Work\work\MTRIX\Data\NCS_157_files.csv
    """
    self.pdb_code_list = [
      '3p0s', '2wws', '2xpj', '2bfu', '3n7x', '2iz9', '1dwn', '1ei7','2wzr',
      '2vf1', '1m1c', '1llc', '1dzl', '2vf9', '3ntt', '4ar2','4gmp', '3vdd',
      '3bcc', '3s4g', '3lob', '3qpr', '1tdi', '1ohg', '3e8k', '2qzv', '2e0z',
      '1wcd', '4bcu', '1vcr', '1ng0', '3dar', '4f5x', '4g93', '1bcc', '2izw',
      '1f2n', '1ny7', '3oah', '1vb4', '2gtl', '2g33', '2zzq', '2ws9', '1c8n',
      '2w4z', '1x9t', '3r0r', '4gb3', '1vsz', '2g34', '2c4y', '1z7s', '1ddl',
      '2bq5', '2c4z', '3fbm', '2gh8', '1qjx', '1f8v', '2iz8', '2bs1', '4aqq',
      '1qju', '1x36', '1w39', '1x35', '1pgw', '2wff', '2vq0', '2fz2', '2fz1',
      '1x9p', '3vbu', '3hag', '4gh4', '3chx', '1pgl', '1a37', '1lp3', '3zfe',
      '4fts', '4fsj', '3raa', '2c4q', '1qjy', '4hl8', '3tn9', '3es5', '1js9',
      '4gbt', '4fte', '2x5i', '2izn', '1zba', '1r2j', '1k5m', '2w4y', '2qqp',
      '4jgy', '4ftb', '4ang', '3zfg', '3zff', '3cji', '2c51', '1vak', '1uf2',
      '1ohf', '3ux1', '4g0r', '3s6p', '3ra4', '2bu1', '1laj', '1a34', '7msf',
      '4iv1', '3ra9', '3ra8', '2xbo', '1h8t', '5msf', '4jgz', '3vbo', '2ztn',
      '1b35', '6msf', '4aed', '3vbs', '3vbr', '3vbf', '3ra2', '3kz4', '1za7',
      '1vb2', '1rhq', '4iv3', '3vbh', '3nou', '3not', '3nop', '2xgk', '2wbh',
      '2qij', '2c50', '2buk', '1wce', '1tnv']



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
      outString = '{0} {1} &> log_{1}'.format(self.com_path,file_name)
      self.commands.append("python {}".format(outString))

  def send_to_queue(self):
    # Send the job to the queue
    easy_qsub.run(
    phenix_source = self.phenix_source,
    where = self.where_to_run_dir,
    # Optional, when you want all jobs to run on machine_name
    qsub_cmd = 'qsub -q all.q@beck',
    commands = self.commands,
    size_of_chunks= self.size_of_chunks) # set the number of commands to send together to the queue.

if __name__ == "__main__":
  queue_job = run_queue_tests()
  queue_job.get_pdb_files()
  queue_job.get_commands()
  queue_job.send_to_queue()

__author__ = 'Youval'
