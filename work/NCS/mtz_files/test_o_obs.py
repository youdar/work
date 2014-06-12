from __future__ import division
from iotbx import mtz
import os

__author__ = 'Youval'


def run():
  # test f_obs relative to Pavel's MTZ collection
  # original location /net/cci-filer2/raid1/share/pdbmtz/mtz_files
  test_path = '/net/cci-filer2/raid1/share/pdbmtz/mtz_files'
  mtz_object = mtz.object(file_name=file_name)
  files_to_test = [
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




def get_structure_factors():
  """
  Get f_obs and r_free_flags
  From cif file if available
  """
  f_obs = None
  r_free_flags = None
  if self.full_path_cif:
    miller_arrays = iotbx.cif.reader(
      file_path=self.full_path_cif).\
      as_miller_arrays(force_symmetry=True)
    # print miller_arrays[0].completeness()
    for ma in miller_arrays:
      if ma.is_xray_amplitude_array():
        # Consider using Bijvoet mates
        ma = ma.average_bijvoet_mates()
        self.f_obs = abs(ma)
        break
      elif not self.f_obs and ma.is_xray_intensity_array():
        # Consider using Bijvoet mates
        ma = ma.average_bijvoet_mates()
        # convert i_obs to f_obs
        self.f_obs = abs(ma.french_wilson(log=null_out()))
  else:
    raise RuntimeError("No cif file.")

  if self.f_obs:
    # self.f_obs.show_summary()
    self.r_free_flags = self.f_obs.generate_r_free_flags()
    # Data completeness: Fraction of unmeasured reflections within the
    # [d_min, d_max] range,where d_min and d_max are highest and lowest
    # resolution of data set correspondingly.
    self.completeness = self.f_obs.array().completeness()
  else:
    raise RuntimeError("Missing amplitude array.")

def process_pdb_and_cif_files(args):
  """
  Fetch pdb file when not working at LBL, if does not exist in working
  directory

  Argument:
  args: (str) 4 letters pdb code or a list where the first variable
        is (str) 4 letters pdb code
  """
  if type(args) == list:
    file_name = args[0]
  else:
    file_name = args
  assert len(file_name) == 4
  assert type(file_name) == str
  # Read pdb file when not working at LBL
  pdb_code = file_name.lower()
  pdb_file_name = pdb_code + '.pdb'
  cif_file_name = pdb_code + '-sf.cif'
  # Get pdb file
  if os.path.isfile(pdb_file_name):
    self.full_path_pdb = os.path.realpath(self.pdb_file_name)
    print 'Using pdb file from local machine (not from MIRROR)'
  else:
    self.full_path_pdb = self.get_full_path(data_type='pdb')
  # Get cif file
  if os.path.isfile(self.cif_file_name):
    self.full_path_cif = os.path.realpath(self.cif_file_name)
    print 'Using cif file from local machine (not from MIRROR)'
  else:
    self.full_path_cif = self.get_full_path(data_type='xray')
  # Collect Fobs and r_free_flags from cif file
  self.get_structure_factors()
  # Process PDB file
  self.process_pdb_file()
