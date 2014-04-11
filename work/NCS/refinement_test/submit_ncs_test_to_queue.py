from __future__ import division
from libtbx.command_line import easy_qsub
import os

class run_queue_tests(object):
  def __init__(self):
    # set working environment
    self.phenix_source = \
      "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
    sources = os.environ["sources_chevy"]
    self.com_path = sources + '/misc_scripts/strict_ncs/test_ncs_refinement.py'
    # where all queue output will be deposited
    self.where_to_run_dir = \
      "/net/cci/youval/Work/work/NCS/junk/pdb_test/queue_job"
    self.pdb_code = []
    self.pdb_file_with_path = []
    # The commands list is a list that will be sent to the queue for processing
    self.commands = []
    # the number of command in each "chunk" sent to the queue
    self.size_of_chunks = 1

  def get_pdb_files(self):
    """() -> list

    the list of 147 file out of 159 that have single NCS with
    MTRIX and CIF files that can be processed

    The list can be found at:
    C:\Phenix\Dev\Work\work\MTRIX\Data\NCS_157_files.csv
    """
    #self.pdb_code_list = [
      #'3p0s', '2wws', '2xpj', '2bfu', '3n7x', '2iz9', '1dwn', '1ei7','2wzr',
      #'2vf1', '1m1c', '1llc', '1dzl', '2vf9', '3ntt', '4ar2','4gmp', '3vdd',
      #'3bcc', '3s4g', '3lob', '3qpr', '1tdi', '1ohg', '3e8k', '2qzv', '2e0z',
      #'1wcd', '4bcu', '1vcr', '1ng0', '3dar', '4f5x', '4g93', '1bcc', '2izw',
      #'1f2n', '1ny7', '3oah', '1vb4', '2gtl', '2g33', '2zzq', '2ws9', '1c8n',
      #'2w4z', '1x9t', '3r0r', '4gb3', '1vsz', '2g34', '2c4y', '1z7s', '1ddl',
      #'2bq5', '2c4z', '3fbm', '2gh8', '1qjx', '1f8v', '2iz8', '2bs1', '4aqq',
      #'1qju', '1x36', '1w39', '1x35', '1pgw', '2wff', '2vq0', '2fz2', '2fz1',
      #'1x9p', '3vbu', '3hag', '4gh4', '3chx', '1pgl', '1a37', '1lp3', '3zfe',
      #'4fts', '4fsj', '3raa', '2c4q', '1qjy', '4hl8', '3tn9', '3es5', '1js9',
      #'4gbt', '4fte', '2x5i', '2izn', '1zba', '1r2j', '1k5m', '2w4y', '2qqp',
      #'4jgy', '4ftb', '4ang', '3zfg', '3zff', '3cji', '2c51', '1vak', '1uf2',
      #'1ohf', '3ux1', '4g0r', '3s6p', '3ra4', '2bu1', '1laj', '1a34', '7msf',
      #'4iv1', '3ra9', '3ra8', '2xbo', '1h8t', '5msf', '4jgz', '3vbo', '2ztn',
      #'1b35', '6msf', '4aed', '3vbs', '3vbr', '3vbf', '3ra2', '3kz4', '1za7',
      #'1vb2', '1rhq', '4iv3', '3vbh', '3nou', '3not', '3nop', '2xgk', '2wbh',
      #'2qij', '2c50', '2buk', '1wce', '1tnv']

    # self.pdb_code_list = [
    #   '3dar', '1vcr', '1r2j', '1a37', '1llc', '1tnv', '1tdi', '1w39', '1ny7',
    #   '1ddl', '1c8n', '2bfu', '4gmp', '3vbr', '3vbu', '3vbo', '4jgy', '3es5',
    #   '3nop', '3not', '3nou', '3bcc', '1bcc', '1z7s', '6msf', '2iz8', '7msf',
    #   '2izn', '2c50', '2c51', '2iz9', '2c4y', '2c4z', '5msf', '2c4q', '2bu1',
    #   '3raa', '3oah', '3ra2', '3ra9', '3ra8', '3ra4', '3qpr', '1ei7', '1a34',
    #   '3chx', '2wbh', '2fz1', '2fz2', '2gh8', '1wcd', '3fbm', '4gb3', '1laj',
    #   '3vbh', '1dzl', '3hag', '4iv3', '1js9', '3n7x', '4gh4', '4jgz', '3tn9',
    #   '4iv1', '1vb2', '1vb4', '1vak', '3s4g', '2buk', '1x36', '4bcu', '1b35',
    #   '2wzr', '1k5m', '2bq5', '1zba', '1pgw', '3vbs', '1x35', '3vbf', '1pgl',
    #   '4fsj', '4fte', '4fts', '2e0z', '4ftb', '2w4y', '2w4z', '2qzv', '3vdd',
    #   '3p0s', '1qjx', '1qjy', '1qju', '3r0r', '2bs1', '2ztn', '1x9t', '2zzq',
    #   '1x9p', '4aqq', '1za7', '4ar2', '2wws', '2xpj', '4hl8', '3ntt', '2vf1',
    #   '3ux1', '2xgk', '2izw', '3cji', '4gbt', '2vq0', '4g93', '2g34', '2qij',
    #   '2g33', '1f2n', '4g0r', '1ng0', '2ws9', '2xbo', '2wff', '1wce', '1dwn',
    #   '2vf9', '3zfe', '3zff', '3zfg', '2x5i', '1h8t', '3lob', '4ang', '2gtl',
    #   '2qqp', '1f8v', '1m1c', '1lp3', '4aed', '3e8k', '1uf2', '1ohg', '1ohf',
    #   '3s6p', '3kz4', '4f5x', '1vsz']

    self.pdb_code_list = [
      '2xpj','1vsz','3e8k','3s6p','4f5x','3dar','1ohg',
      '3kz4','1uf2','3kz4','3zfe','4f5x']

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
      outString = '{0} {1} -quiet &> log_{1}'.format(self.com_path,file_name)
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


