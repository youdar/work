from __future__ import division
from libtbx.utils import Sorry
import os,sys
from iotbx import pdb
import cPickle as pickle
from scitbx import matrix


def run(file_name):
    pdb_inp = pdb.input(file_name=file_name)          # read the pdb file data
    try:
        TRASFORM_info = pdb_inp.process_mtrix_records(error_handle=True,eps=1e-2)
    except Sorry as e:
        print os.path.basename(file_name)


def test_matrix(s,eps=1e-3):
    # s is a list of 9 float numbers
    R = matrix.sqr(s)
    t1 = R*R.transpose()
    t2 = matrix.identity(3)
    det_is = R.determinant()
    det_tet = abs(det_is - 1) < eps
    inv_test = abs(t1 - t2) < eps
    return  det_is, det_tet, inv_test


if __name__=='__main__':
    run(sys.argv[1])