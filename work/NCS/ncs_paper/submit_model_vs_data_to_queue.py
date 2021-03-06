from __future__ import division
import collect_ncs_files
from libtbx.command_line import easy_qsub
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
      '3ksb', '2p5t', '4h1l', '4gh4', '2ja7', '2hn2', '4nkz', '3km2', '4gx2',
      '3ks2', '1iss', '4h1i', '1n3n', '4gk5', '3hxk', '3g05', '3ksa', '4hi2',
      '1ott', '4kn4', '3nmj', '4hic', '2gw1', '4h37', '4gx5']

    self.pdb_log_file_list = [
      '4g1u', '4hj0', '4nry', '4a7p', '3vbs', '4qjh', '3b8c', '2h1j', '2pf4',
      '3pi4', '3bbp', '4u7u', '4l6y', '3vwu', '3n97', '3u60', '1nov', '4od4',
      '4od5', '4lrx', '3u61', '3p2d', '1wce', '4kr7', '2fjh', '2w29', '2ost',
      '2o94', '1f8v', '3l4b', '4u4f', '3wcn', '3asn', '3be0', '3rjr', '4fn9',
      '2fs3', '3fzj', '1tnv', '2r1a', '3oij', '3fm7', '4fqk', '4fnr', '3b61',
      '2xpj', '3tu4', '4fqm', '4x4q', '3u5z', '3rfu', '3hqq', '2xyo', '3nou',
      '4x4r', '4fnu', '4pdw', '2fsy', '3rh7', '3bna', '4u0h', '2vf9', '3v4p',
      '4ryj', '2r0q', '3q4f', '3g76', '1fu2', '3frt', '3uo7', '4hl8', '1uf2',
      '4qsm', '4f5x', '3kz4', '3l73', '3vs7', '3txx', '1ohg', '3t4a', '1gr5',
      '1fub', '3l5j', '4pqp', '4u6g', '4idw', '1m0f', '4ld9', '3ug6', '4aed',
      '4qt0', '2r6f', '4u6u', '4lsx', '4f2m', '3pt6', '3r4d', '4ffz', '2gqq',
      '3l3o', '1lzh', '4puf', '1lp3', '4hkj', '4fb3', '2vf1', '2wws', '2xxa',
      '2w4b', '3gfq', '3gdu', '4fi3', '2frp', '3cz3', '3ze1', '3zjy', '4qvz',
      '2ft1', '3v4v', '2vq0', '4nko', '4gi2', '4hg4', '3uaj', '3zfg']

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


