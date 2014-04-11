from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import os

def run():
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").readlines()
  pdb_code_list = [
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

  pdb_code_list = [
    '3dar', '1vcr', '1r2j', '1a37', '1llc', '1tnv', '1tdi', '1w39', '1ny7',
    '1ddl', '1c8n', '2bfu', '4gmp', '3vbr', '3vbu', '3vbo', '4jgy', '3es5',
    '3nop', '3not', '3nou', '3bcc', '1bcc', '1z7s', '6msf', '2iz8', '7msf',
    '2izn', '2c50', '2c51', '2iz9', '2c4y', '2c4z', '5msf', '2c4q', '2bu1',
    '3raa', '3oah', '3ra2', '3ra9', '3ra8', '3ra4', '3qpr', '1ei7', '1a34',
    '3chx', '2wbh', '2fz1', '2fz2', '2gh8', '1wcd', '3fbm', '4gb3', '1laj',
    '3vbh', '1dzl', '3hag', '4iv3', '1js9', '3n7x', '4gh4', '4jgz', '3tn9',
    '4iv1', '1vb2', '1vb4', '1vak', '3s4g', '2buk', '1x36', '4bcu', '1b35',
    '2wzr', '1k5m', '2bq5', '1zba', '1pgw', '3vbs', '1x35', '3vbf', '1pgl',
    '4fsj', '4fte', '4fts', '2e0z', '4ftb', '2w4y', '2w4z', '2qzv', '3vdd',
    '3p0s', '1qjx', '1qjy', '1qju', '3r0r', '2bs1', '2ztn', '1x9t', '2zzq',
    '1x9p', '4aqq', '1za7', '4ar2', '2wws', '2xpj', '4hl8', '3ntt', '2vf1',
    '3ux1', '2xgk', '2izw', '3cji', '4gbt', '2vq0', '4g93', '2g34', '2qij',
    '2g33', '1f2n', '4g0r', '1ng0', '2ws9', '2xbo', '2wff', '1wce', '1dwn',
    '2vf9', '3zfe', '3zff', '3zfg', '2x5i', '1h8t', '3lob', '4ang', '2gtl',
    '2qqp', '1f8v', '1m1c', '1lp3', '4aed', '3e8k', '1uf2', '1ohg', '1ohf',
    '3s6p', '3kz4', '4f5x', '1vsz']



  pdb_length = []
  for fn in pdb_code_list:
    file_path = [x.strip() for x in pdb_files if x[-12:-8]==fn][0]
    print file_path
    full_path = os.path.join(pdb_dir,file_path)
    m = multimer(
      file_name=full_path,
      reconstruction_type='cau',
      error_handle=True,eps=1e-2)
    pdb_length.append([len(m.sites_cart()),fn])
  pdb_length.sort()
  print '='*60
  print pdb_length[:5]
  print '='*60

if __name__=='__main__':
  print 'start'
  run()
  print 'done'










