from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex

asu = """\
CRYST1   17.101   19.080   20.163  90.00  90.00  90.00 P 1
SCALE1      0.058476  0.000000  0.000000        0.00000
SCALE2      0.000000  0.052411  0.000000        0.00000
SCALE3      0.000000  0.000000  0.049596        0.00000
ATOM      1  N   THR A   1       8.111  11.079  10.645  1.00 20.00           N
ATOM      2  CA  THR A   1       8.000   9.721  10.125  1.00 20.00           C
ATOM      3  C   THR A   1       8.075   8.693  11.249  1.00 20.00           C
ATOM      4  O   THR A   1       8.890   8.817  12.163  1.00 20.00           O
ATOM      5  CB  THR A   1       9.101   9.420   9.092  1.00 20.00           C
ATOM      6  OG1 THR A   1       9.001  10.342   8.000  1.00 20.00           O
ATOM      7  CG2 THR A   1       8.964   8.000   8.565  1.00 20.00           C
TER
ATOM      1  N   THRaa   1       3.097   7.753  15.237  1.00 20.00           N
ATOM      2  CA  THRaa   1       3.613   7.314  13.945  1.00 20.00           C
ATOM      3  C   THRaa   1       4.967   6.628  14.097  1.00 20.00           C
ATOM      4  O   THRaa   1       5.824   7.086  14.852  1.00 20.00           O
ATOM      5  CB  THRaa   1       3.752   8.492  12.963  1.00 20.00           C
ATOM      6  OG1 THRaa   1       2.473   9.106  12.765  1.00 20.00           O
ATOM      7  CG2 THRaa   1       4.291   8.010  11.625  1.00 20.00           C
TER
ATOM      1  N   THRab   1       5.422   0.660  16.493  1.00 20.00           N
ATOM      2  CA  THRab   1       5.208   1.362  15.233  1.00 20.00           C
ATOM      3  C   THRab   1       6.410   2.230  14.875  1.00 20.00           C
ATOM      4  O   THRab   1       6.982   2.901  15.734  1.00 20.00           O
ATOM      5  CB  THRab   1       3.947   2.244  15.284  1.00 20.00           C
ATOM      6  OG1 THRab   1       2.801   1.430  15.559  1.00 20.00           O
ATOM      7  CG2 THRab   1       3.746   2.965  13.960  1.00 20.00           C
TER
"""

def run():
  pdb_inp1 = iotbx.pdb.input(source_info=None, lines=asu)
  ph1 = pdb_inp1.construct_hierarchy()
  sites_cart1 = ph1.atoms().extract_xyz()

  pdb_inp2 = iotbx.pdb.input(file_name='test_asu.pdb')
  ph2 = pdb_inp2.construct_hierarchy()
  sites_cart2 = ph2.atoms().extract_xyz()

  d = flex.sqrt((sites_cart1 - sites_cart2).dot()).as_double().min_max_mean().as_tuple()
  print 'sites_cart from string in file - sites_cart from same pdb, but from file'
  print d
  print '-'*50




if __name__=='__main__':
  run()