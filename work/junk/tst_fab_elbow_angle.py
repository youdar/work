from __future__ import division
from mmtbx.utils.fab_elbow_angle import fab_elbow_angle
import libtbx.load_env
import os
import iotbx.pdb
from libtbx.test_utils import approx_equal

def exercise():
  """
  Exercise FAB elbow angle calculation
  target angles are taken from Stanfield, et al., JMB 2006
  """
  for it in [("1bbd", 127,114,118,7),
             ("1bbd", 127,107,113,2),
             ("7fab", 132,104,117,10),
             ("1dba", 183,107,113,6),
             ("1plg", 190,112,117,4),
             ("1nl0", 220,107,113,7)]:
    pdb_code, target_angle,limit_light, limit_heavy, tolerance = it
    pdb_file_name = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/mmtbx/fab_elbow_angle/%s.pdb.gz"%pdb_code,
      test=os.path.isfile)
    ph = iotbx.pdb.input(file_name=pdb_file_name).construct_hierarchy()
    elbow_angle = fab_elbow_angle(
      pdb_hierarchy=ph,
      limit_light=limit_light,
      limit_heavy=limit_heavy).fab_elbow_angle
    print "inputs:", it, "output:", elbow_angle
    assert approx_equal(elbow_angle, target_angle, tolerance), \
      [elbow_angle, target_angle]

if __name__ == "__main__":
  exercise()
