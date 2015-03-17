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
      '4c5y', '3t5q', '4lyg', '4tkr', '4gg6', '2zom', '4g1u', '3vqk', '4rgf',
      '4gup', '4pa7', '4bsr', '4f5c', '3l9k', '4hj0', '4frt', '4l3b', '4hfc',
      '2h5k', '3f7a', '4nry', '4a7p', '3vbs', '4by0', '4bst', '4fbg', '4qjh',
      '3blv', '3b8c', '2yq4', '3bba', '4mw8', '3ula', '4pbw', '2h1j', '2pf4',
      '4pog', '3ubv', '3tth', '4pa9', '2xdd', '3p56', '3pi4', '4ilb', '3bbp',
      '2nww', '3tru', '2h1l', '1s5l', '2qeq', '3bx3', '4u7u', '4l6y', '4md7',
      '4md8', '3n56', '3uy9', '3r6m', '3vwu', '3n97', '3n57', '3m9s', '4oiz',
      '3rqu', '4md9', '3rbl', '2xwb', '3rwt', '3qkv', '4p2z', '3u60', '1nov',
      '4od4', '4od5', '4lrx', '3muu', '4n1k', '3kwq', '4ai6', '3u61', '4orh',
      '3p2d', '1wce', '4kr7', '4hfd', '4kr8', '4uf6', '4xr7', '2fjh', '3t05',
      '2nwx', '3ukq', '2xz1', '2w29', '4atb', '4mvc', '2ost', '2o94', '4mvd',
      '1f8v', '4aze', '3t07', '4akg', '4aby', '3l4b', '4u4f', '3n5u', '2xk5',
      '3u9u', '3ux4', '3wcn', '3asn', '3be0', '3rjr', '4fn9', '2xny', '3b1k',
      '2fs3', '3fzj', '4tsy', '4mvm', '4r89', '3m9d', '3zef', '3usj', '3usk',
      '1tnv', '3pqy', '3csy', '1xri', '2xeq', '4atv', '2q47', '4av1', '2r1a',
      '3uso', '3rl0', '4x4o', '4mvo', '3oij', '3mf0', '4mvq', '4r8p', '4ri0',
      '4kw6', '3ztl', '2ycb', '3uu2', '3fm7', '4fqk', '3t0t', '3ta9', '4lzo',
      '4fnr', '3b61', '2xpj', '3vr3', '4x5t', '2xka', '3tu4', '4fqm', '3l1r',
      '4asn', '4aso', '4ass', '4pn6', '2qij', '4x4q', '3u5z', '3rfu', '4f7r',
      '4o6o', '4aym', '3fg6', '3jxe', '3hqq', '2xyo', '4e36', '4c31', '4p5h',
      '4mvr', '3nou', '4x4r', '4ddq', '4fnu', '3q16', '4mvu', '4pdw', '4p5o',
      '2fsy', '4eo2', '3q17', '3rh7', '3rcc', '3sku', '2x8c', '2yhm', '2yhn',
      '2xwj', '3waz', '4mvz', '4pr9', '2xkb', '3r3m', '3q41', '4o6p', '4q79',
      '4ubu', '4gko', '4c3o', '3bna', '3t79', '4u0h', '4bql', '4aj5', '2vf9',
      '3zwg', '3u7u', '4mf3', '4au5', '3b54', '3az8', '3q4c', '4m9s', '3ojb',
      '4ddv', '3v4p', '4ryj', '2y6r', '4wp0', '2r0q', '3q4f', '3v2h', '3g93',
      '3sen', '3taf', '4u7b', '3g76', '4f9l', '3zzi', '4cx7', '1fu2', '3frt',
      '3uo7', '4hl8', '2ycw', '1uf2', '4qsm', '4f5x', '4biy', '3kz4', '3o3c',
      '4ddx', '4s0t', '3l73', '3vs7', '3txx', '4u7d', '1ohg', '4q85', '4oia',
      '4kk5', '3m6a', '4m9z', '3t4a', '1gr5', '4fsx', '1fub', '4cxa', '4p2q',
      '4gz8', '3uob', '4gza', '3l5j', '4pqp', '4tzz', '3s4g', '3uux', '3zq4',
      '4mtf', '2qva', '4ut9', '3sqv', '4u6g', '3uqr', '4q58', '4idw', '3mgs',
      '4mtg', '1m0f', '4kk6', '4boz', '4ld9', '3ug6', '4h56', '4ldb', '4fzh',
      '4aed', '4qt0', '3o9u', '2r6f', '4u6u', '4tk2', '2pvs', '4fgn', '4ldi',
      '4p2r', '4txy', '4il4', '4lsx', '3r1a', '3szk', '4tk4', '3lrc', '3r1b',
      '3s24', '4hel', '2raq', '3ria', '4f2m', '3pt6', '3r4d', '4ei5', '4pk6',
      '4ffz', '2gqq', '4h5m', '3l3o', '4kk9', '4i3h', '2wsx', '1lzh', '4kka',
      '4mto', '4puf', '1lp3', '4oyj', '4mtp', '2qjk', '3jua', '4hkj', '4kkb',
      '2yje', '4d5s', '4fb3', '3dru', '4f52', '4f88', '3zxu', '3zrt', '2yjf',
      '4kkc', '4wfg', '4nc9', '2vf1', '4wfh', '2wws', '2xxa', '2w4b', '3lvg',
      '3zla', '4kud', '4gjw', '3gfq', '4x6l', '2w6g', '2w2w', '3gdu', '4fi3',
      '2frp', '3p0h', '4wnw', '3cz3', '3ze1', '3v65', '3wo3', '3zkc', '2xkp',
      '4n0a', '2zp2', '3pbk', '2w6h', '3zqj', '4wyk', '4mw3', '3zjy', '4qvz',
      '2ft1', '3v4v', '3fop', '3p83', '2xte', '4p30', '3uzb', '2w6i', '2vq0',
      '4mhj', '2y4t', '3oj4', '2xro', '4nko', '4q9j', '4gfa', '2px0', '4gi2',
      '2ygc', '3pso', '2pr1', '4q9t', '3m4q', '4hg4', '3h1k', '3uaj', '4x8f',
      '3g6j', '2wy1', '3h0g', '4pa4', '3n6v', '3m05', '4lo7', '2xg8', '4pa6',
      '4c2t', '2yoe', '3g33', '3h1u', '3ny0', '3zfg', '4k95', '2g0b', '3hou',
      '3hfs', '4nkz', '4hi2']

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


