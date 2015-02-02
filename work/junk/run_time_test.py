from __future__ import division
from mmtbx.ncs.ncs_search import Score_record
from mmtbx.ncs.ncs_search import get_matching_res_indices
from mmtbx.ncs.ncs_search import initialize_score_matrix
from scitbx import matrix
import cProfile
import pstats
import time

def g(x,y):
  return 0.5*(1 + x + y)*(x + y) + y

def f(r):
  (r11,r12,r13,r21,r22,r23,r31,r32,r33) = (10*r.transpose()).round(0).elems
  return str(int(g(r11,g(r22,r33))))

def run1():
  for i in range(120):
    for j in range(i,120):
      r = matrix.sqr([0.1235,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9])
      print f(r)

def run2():
  for i in range(10000):
    x = 1000
    y = 1014
    if x > y:
      z = x
    else:
      z = y

def run3():
  for i in range(10000):
    x = 1000
    y = 1014
    z = max(x,y)

def run4():
  t = 0
  R = initialize_score_matrix(row=232,col=232,max_score=100)
  n = 10
  for k in xrange(n):
    t0=time.time()
    x = {}
    for i in xrange(1,232):
      for j in xrange(1,232):
        # if j > 220: break
        s = 5
        next_row = 7
        R[i,j] = Score_record(score=s)
        s1,mc = gap_score(R[i,j],40)
        s2,cm = gap_score(R[i,j],40)
        # if R.has_key((0,0)):
        #   s3 = s1*s2
        # if R.has_key((0,1)):
        #   s3 = s2*s2
        # if R.has_key((0,1)):
        #   s3 = s1*s1

        s3 = s1*s2
        s3 = s3*s2
        s3 = s1*s1

        # s2 = s3 = 0
        # if x.has_key(1):
        #   s2 = 0
        # if x.has_key(2):
        #   s3 = 0
    t += time.time()-t0

  print "Time: %6.4f"%(t/n)

def run5():
  td = 0
  tl = 0
  to = 0
  n = 5
  row = 500
  col = 500
  for k in xrange(n):
    t0=time.time()
    use_dict(row,col)
    td += time.time()-t0
    t0=time.time()
    use_list(row,col)
    tl += time.time()-t0
    t0=time.time()
    use_dict_now(row,col)
    to += time.time()-t0
  s = 'dict time: {0:.4f}  list time: {1:.4f}  dict time: {2:.4f}'
  print s.format(td/n,tl/n,to/n)

def run6():
  td = 0
  tl = 0
  n = 10
  for k in xrange(n):
    t0=time.time()
    res_alignment1()
    td += time.time()-t0
    t0=time.time()
    res_alignment2()
    tl += time.time()-t0
  s = 'all dict first: {0:.4f}  build dict during: {1:.4f}'
  print s.format(td/n,tl/n)
  print (td/n - tl/n) * 400

def use_dict(row,col,max_score):
  R = {}
  for i in xrange(row + 1):
    score_row = max_score - i
    if score_row > 0: nr = i + 1
    else: nr = None
    R[i,0] = Score_record(score=score_row,origin=(i-1,0),next_record=nr)
    for j in xrange(1,col + 1):
      if i == 0:
        score_col = max_score - i
        R[i,j] = Score_record(score=score_col,origin=(0,j-1))
      else:
        R[i,j] = Score_record()
  return R

def use_list(row,col):
  R = []
  max_score = 500
  for i in xrange(row + 1):
    score_row = max_score - i
    r  =[Score_record(score=score_row,origin=(i-1,0))]
    for j in xrange(col + 1):
      if i == 0:
        score_col = max_score - i
        r.append(Score_record(score=score_col,origin=(0,j-1)))
      else:
        r.append(Score_record())
    R.append(r)
  return R

def use_dict_now(row,col,max_score):
  R = {}
  for i in xrange(row + 1):
    score = max_score - i
    R[i,0] = Score_record(score=score,origin=(i-1,0))
  for i in xrange(col + 1):
    score = max_score - i
    R[0,i] = Score_record(score=score,origin=(0,i-1))
  return R



def res_alignment1():
  """
  Align two chains hierarchies.

  Penalize misalignment, gaps, when contiguous section is less
  than min_contig_length and when the
  Do not give any points for alignment. (score only change when "bad" things
  happens)

  Args:
    seq_a, seq_b (lists of str): list of characters to compare
    min_fraction_domain (float): min percent of similarity between hierarchies
      similarity define as:
      (number of matching res) / (number of res in longer chain)
    min_contig_length (int): domain < min_contig_length rejected

  Returns:
    aligned_sel_a (list): the indices of the aligned components of seq_a
    aligned_sel_b (list): the indices of the aligned components of seq_b
    similarity (float): actual similarity between hierarchies
  """
  min_contig_length=10
  min_fraction_domain=0.2
  a = len(seq_a)
  b = len(seq_b)
  # Check for the basic cases
  if (a == 0) or (b == 0): return [],[],0
  if seq_a == seq_b: return range(a),range(a),1.0
  min_matches = min(min_contig_length,min(a,b))
  misalign_penalty = 100000
  # limit the number of mis-alignments
  max_mis_align = int((1 - min_fraction_domain) * max(a,b))
  # Starting score according to the required similarity
  score = max_mis_align + 1
  # build score matrix
  R = use_dict(row=a,col=b,max_score=score)
  # populate score matrix
  for j in xrange(1,b + 1):
    for i in xrange(1,a + 1):
      # We want to halt when there is not way to get a good score
      not_aligned = (seq_a[i-1].lower() != seq_b[j-1].lower())
      s1 = R[i-1,j-1].score - misalign_penalty * not_aligned
      s2,mc2 = gap_score(R[i-1,j],min_matches)
      s3,mc3 = gap_score(R[i,j-1],min_matches)
      if (s1 >= s2) and (s1 >= s3):
        s = s1
        i_temp,j_temp = i-1,j-1
        match_count = R[(i_temp,j_temp)].match_count + 1
        consecutive_matches = R[(i_temp,j_temp)].consecutive_matches + 1
      elif (s2 >= s1) and (s2 >= s3):
        i_temp,j_temp = i-1,j
        s = s2
        match_count = mc2
        consecutive_matches = 0
      else:
        i_temp,j_temp = i,j-1
        s = s3
        match_count = mc3
        consecutive_matches = 0
      # if i > (0.9 * a):
      #   aaa = 10
      #   break
      # else:
      #   aaa = -10
      R[i,j].score = s
      R[i,j].origin = (i_temp,j_temp)
      R[i,j].match_count = match_count
      R[i,j].consecutive_matches = consecutive_matches

def get_score(r,min_matches):
  s = r.score - 1
  mc = r.match_count
  if r.consecutive_matches < min_matches:
    s -= r.consecutive_matches
    mc -= r.consecutive_matches
  return s,mc

def res_alignment2(seq_a, seq_b,
                   min_contig_length=10,
                   min_fraction_domain=0.2):
  a = len(seq_a)
  b = len(seq_b)
  # Check for the basic cases
  if (a == 0) or (b == 0): return [],[],0
  if seq_a == seq_b: return range(a),range(a),1.0
  min_matches = min(min_contig_length,min(a,b))
  misalign_penalty = 100000
  # limit the number of mis-alignments
  max_mis_align = int((1 - min_fraction_domain) * max(a,b))
  # Starting score according to the required similarity
  score = max_mis_align + 1
  # build score matrix
  R = initialize_score_matrix(row=a,col=b,max_score=score)
  # populate score matrix
  for j in xrange(1,b + 1):
    for i in xrange(1,a + 1):
      not_aligned = (seq_a[i-1] != seq_b[j-1])
      s1 = R[i-1,j-1].score - misalign_penalty * not_aligned
      #
      s2,mc2 = get_score(R[i-1,j],min_matches)
      # s2 = R[i-1,j].score - 1
      # mc2 = R[i-1,j].match_count
      # if R[i-1,j].consecutive_matches < min_matches:
      #   s2 -= R[i-1,j].consecutive_matches
      #   mc2 -= R[i-1,j].consecutive_matches
      #
      s3,mc3 = get_score(R[i,j-1],min_matches)
      # s3 = R[i,j-1].score - 1
      # mc3 = R[i,j-1].match_count
      # if R[i,j-1].consecutive_matches < min_matches:
      #   s3 -= R[i,j-1].consecutive_matches
      #   mc3 -= R[i,j-1].consecutive_matches
      #
      if (s1 >= s2) and (s1 >= s3):
        s = s1
        i_temp,j_temp = i-1,j-1
        match_count = R[(i_temp,j_temp)].match_count + 1
        consecutive_matches = R[(i_temp,j_temp)].consecutive_matches + 1
      elif (s2 >= s1) and (s2 >= s3):
        i_temp,j_temp = i-1,j
        s = s2
        match_count = mc2
        consecutive_matches = 0
      else:
        i_temp,j_temp = i,j-1
        s = s3
        match_count = mc3
        consecutive_matches = 0
      if s == 0:
        x = 500
        if i > 0.5 * a:
          break
        else:
          x = 500
      R[i,j].score = s
      R[i,j].origin = (i_temp,j_temp)
      R[i,j].match_count = match_count
      R[i,j].consecutive_matches = consecutive_matches
  #
  aligned_sel_a, aligned_sel_b, similarity = get_matching_res_indices(
    R=R,row=a,col=b,min_fraction_domain=min_fraction_domain)
  return aligned_sel_a, aligned_sel_b, similarity
def get_score(r,min_matches):
  s = r.score - 1
  mc = r.match_count
  if r.consecutive_matches < min_matches:
    s -= r.consecutive_matches
    mc -= r.consecutive_matches
  return s,mc

seq_a = ['ARG', 'PHE', 'GLN', 'TYR', 'LEU', 'VAL', 'LYS', 'ASN', 'GLN', 'ASN',
         'LEU', 'HIS', 'ILE', 'ASP', 'TYR', 'LEU', 'ALA', 'LYS', 'LYS', 'LEU',
         'HIS', 'ASP', 'ILE', 'GLU', 'GLU', 'GLU', 'TYR', 'ASN', 'LYS', 'LEU',
         'THR', 'HIS', 'ASP', 'VAL', 'ASP', 'LYS', 'LYS', 'THR', 'ILE', 'ARG',
         'GLN', 'LEU', 'LYS', 'ALA', 'ARG', 'ILE', 'SER', 'ASN', 'LEU', 'GLU',
         'GLU', 'HIS', 'HIS', 'CYS', 'ASP', 'GLU', 'HIS', 'GLU', 'SER', 'GLU',
         'CYS', 'ARG', 'GLY', 'ASP', 'VAL', 'PRO', 'GLU', 'CYS', 'ILE', 'HIS',
         'ASP', 'LEU', 'LEU', 'PHE', 'CYS', 'ASP', 'GLY', 'GLU', 'LYS', 'ASP',
         'CYS', 'ARG', 'ASP', 'GLY', 'SER', 'ASP', 'GLU', 'ASP', 'PRO', 'GLU',
         'THR', 'CYS', 'SER', 'LEU', 'ASN', 'ILE', 'THR', 'HIS', 'VAL', 'GLY',
         'SER', 'SER', 'TYR', 'THR', 'GLY', 'LEU', 'ALA', 'THR', 'TRP', 'THR',
         'SER', 'CYS', 'GLU', 'ASP', 'LEU', 'ASN', 'PRO', 'ASP', 'HIS', 'ALA',
         'ILE', 'VAL', 'THR', 'ILE', 'THR', 'ALA', 'ALA', 'HIS', 'ARG', 'LYS',
         'SER', 'PHE', 'PHE', 'PRO', 'ASN', 'ARG', 'VAL', 'TRP', 'LEU', 'ARG',
         'ALA', 'THR', 'LEU', 'SER', 'TYR', 'GLU', 'LEU', 'ASP', 'GLU', 'HIS',
         'ASP', 'HIS', 'THR', 'VAL', 'SER', 'THR', 'THR', 'GLN', 'LEU', 'ARG',
         'GLY', 'PHE', 'TYR', 'ASN', 'PHE', 'GLY', 'LYS', 'ARG', 'GLU', 'LEU',
         'LEU', 'LEU', 'ALA', 'PRO', 'LEU', 'LYS', 'GLY', 'GLN', 'SER', 'GLU',
         'GLY', 'TYR', 'GLY', 'VAL', 'ILE', 'CYS', 'ASP', 'PHE', 'ASN', 'LEU',
         'GLY', 'ASP', 'ASP', 'ASP', 'HIS', 'ALA', 'ASP', 'CYS', 'LYS', 'ILE',
         'VAL', 'VAL', 'PRO', 'SER', 'SER', 'LEU', 'PHE', 'VAL', 'CYS', 'ALA',
         'HIS', 'PHE', 'ASN', 'ALA', 'GLN', 'ARG', 'TYR']

seq_b = ['LEU', 'ASP', 'PRO', 'ARG', 'LEU', 'GLY', 'ALA', 'ASN', 'ALA', 'PHE',
         'LEU', 'ILE', 'ILE', 'ARG', 'LEU', 'ASP', 'ARG', 'ILE', 'ILE', 'GLU',
         'LYS', 'LEU', 'ARG', 'THR', 'LYS', 'LEU', 'ASP', 'GLU', 'ALA', 'GLU',
         'LYS', 'ILE', 'ASP', 'PRO', 'GLU', 'HIS', 'PHE', 'VAL', 'SER', 'GLU',
         'ILE', 'ASP', 'ALA', 'ARG', 'VAL', 'THR', 'LYS', 'ILE', 'GLU', 'GLY',
         'THR', 'HIS', 'CYS', 'GLU', 'LYS', 'ARG', 'THR', 'PHE', 'GLN', 'CYS',
         'GLY', 'GLY', 'ASN', 'GLU', 'GLN', 'GLU', 'CYS', 'ILE', 'SER', 'ASP',
         'LEU', 'LEU', 'VAL', 'CYS', 'ASP', 'GLY', 'HIS', 'LYS', 'ASP', 'CYS',
         'HIS', 'ASN', 'ALA', 'HIS', 'ASP', 'GLU', 'ASP', 'PRO', 'ASP', 'VAL',
         'CYS', 'ASP', 'THR', 'SER', 'VAL', 'VAL', 'LYS', 'ALA', 'GLY', 'ASN',
         'VAL', 'PHE', 'SER', 'GLY', 'THR', 'SER', 'THR', 'TRP', 'HIS', 'GLY',
         'CYS', 'LEU', 'ALA', 'ARG', 'GLU', 'ASP', 'HIS', 'VAL', 'THR', 'ARG',
         'ILE', 'THR', 'ILE', 'THR', 'ALA', 'SER', 'LYS', 'ARG', 'ARG', 'LYS',
         'PHE', 'PHE', 'THR', 'ALA', 'ARG', 'ILE', 'TRP', 'LEU', 'ARG', 'ALA',
         'LEU', 'VAL', 'GLU', 'SER', 'GLU', 'LEU', 'GLU', 'ARG', 'HIS', 'GLY',
         'GLU', 'ASN', 'VAL', 'THR', 'SER', 'SER', 'PHE', 'ASN', 'ALA', 'LYS',
         'GLY', 'TYR', 'TYR', 'ASN', 'PHE', 'ALA', 'SER', 'ARG', 'ARG', 'LEU',
         'ILE', 'LEU', 'LEU', 'PRO', 'THR', 'ASP', 'ASP', 'HIS', 'ASP', 'ASP',
         'HIS', 'LEU', 'ALA', 'VAL', 'VAL', 'CYS', 'SER', 'PHE', 'ASN', 'ARG',
         'GLY', 'ASP', 'ASN', 'GLU', 'ARG', 'ALA', 'GLU', 'CYS', 'HIS', 'ARG',
         'VAL', 'THR', 'GLU', 'ALA', 'THR', 'LEU', 'HIS', 'GLN', 'CYS', 'ALA',
         'ASP', 'LEU', 'PHE', 'VAL', 'THR', 'LEU', 'GLU', 'GLU', 'HIS', 'ASP']

seq_a = ['ALA', 'GLY', 'TYR', 'ASP', 'ARG', 'HIS', 'ILE', 'THR', 'ILE', 'PHE',
         'SER', 'PRO', 'GLU', 'GLY', 'ARG', 'LEU', 'TYR', 'GLN', 'VAL', 'GLU',
         'TYR', 'ALA', 'PHE', 'LYS', 'ALA', 'THR', 'ASN', 'GLN', 'THR', 'ASN',
         'ILE', 'ASN', 'SER', 'LEU', 'ALA', 'VAL', 'ARG', 'GLY', 'LYS', 'ASP',
         'CYS', 'THR', 'VAL', 'VAL', 'ILE', 'SER', 'GLN', 'LYS', 'LYS', 'VAL',
         'PRO', 'ASP', 'LYS', 'LEU', 'LEU', 'ASP', 'PRO', 'THR', 'THR', 'VAL',
         'SER', 'TYR', 'ILE', 'PHE', 'CYS', 'ILE', 'SER', 'ARG', 'THR', 'ILE',
         'GLY', 'MET', 'VAL', 'VAL', 'ASN', 'GLY', 'PRO', 'ILE', 'PRO', 'ASP',
         'ALA', 'ARG', 'ASN', 'ALA', 'ALA', 'LEU', 'ARG', 'ALA', 'LYS', 'ALA',
         'GLU', 'ALA', 'ALA', 'GLU', 'PHE', 'ARG', 'TYR', 'LYS', 'TYR', 'GLY',
         'TYR', 'ASP', 'MET', 'PRO', 'CYS', 'ASP', 'VAL', 'LEU', 'ALA', 'LYS',
         'ARG', 'MET', 'ALA', 'ASN', 'LEU', 'SER', 'GLN', 'ILE', 'TYR', 'THR',
         'GLN', 'ARG', 'ALA', 'TYR', 'MET', 'ARG', 'PRO', 'LEU', 'GLY', 'VAL',
         'ILE', 'LEU', 'THR', 'PHE', 'VAL', 'SER', 'VAL', 'ASP', 'GLU', 'GLU',
         'LEU', 'GLY', 'PRO', 'SER', 'ILE', 'TYR', 'LYS', 'THR', 'ASP', 'PRO',
         'ALA', 'GLY', 'TYR', 'TYR', 'VAL', 'GLY', 'TYR', 'LYS', 'ALA', 'THR',
         'ALA', 'THR', 'GLY', 'PRO', 'LYS', 'GLN', 'GLN', 'GLU', 'ILE', 'THR',
         'THR', 'ASN', 'LEU', 'GLU', 'ASN', 'HIS', 'PHE', 'LYS', 'LYS', 'SER',
         'LYS', 'ILE', 'ASP', 'HIS', 'ILE', 'ASN', 'GLU', 'GLU', 'SER', 'TRP',
         'GLU', 'LYS', 'VAL', 'VAL', 'GLU', 'PHE', 'ALA', 'ILE', 'THR', 'HIS',
         'MET', 'ILE', 'ASP', 'ALA', 'LEU', 'GLY', 'THR', 'GLU', 'PHE', 'SER',
         'LYS', 'ASN', 'ASP', 'LEU', 'GLU', 'VAL', 'GLY', 'VAL', 'ALA', 'THR',
         'LYS', 'ASP', 'LYS', 'PHE', 'PHE', 'THR', 'LEU', 'SER', 'ALA', 'GLU',
         'ASN', 'ILE', 'GLU', 'GLU', 'ARG', 'LEU', 'VAL', 'ALA', 'ILE', 'ALA',
         'GLU', 'GLN', 'ASP']

seq_b = ['SER', 'ARG', 'ARG', 'TYR', 'ASP', 'SER', 'ARG', 'THR', 'THR', 'ILE',
         'PHE', 'SER', 'PRO', 'GLU', 'GLY', 'ARG', 'LEU', 'TYR', 'GLN', 'VAL',
         'GLU', 'TYR', 'ALA', 'LEU', 'GLU', 'SER', 'ILE', 'SER', 'HIS', 'ALA',
         'GLY', 'THR', 'ALA', 'ILE', 'GLY', 'ILE', 'MET', 'ALA', 'SER', 'ASP',
         'GLY', 'ILE', 'VAL', 'LEU', 'ALA', 'ALA', 'GLU', 'ARG', 'LYS', 'VAL',
         'THR', 'SER', 'THR', 'LEU', 'LEU', 'GLU', 'GLN', 'ASP', 'THR', 'SER',
         'THR', 'GLU', 'LYS', 'LEU', 'TYR', 'LYS', 'LEU', 'ASN', 'ASP', 'LYS',
         'ILE', 'ALA', 'VAL', 'ALA', 'VAL', 'ALA', 'GLY', 'LEU', 'THR', 'ALA',
         'ASP', 'ALA', 'GLU', 'ILE', 'LEU', 'ILE', 'ASN', 'THR', 'ALA', 'ARG',
         'ILE', 'HIS', 'ALA', 'GLN', 'ASN', 'TYR', 'LEU', 'LYS', 'THR', 'TYR',
         'ASN', 'GLU', 'ASP', 'ILE', 'PRO', 'VAL', 'GLU', 'ILE', 'LEU', 'VAL',
         'ARG', 'ARG', 'LEU', 'SER', 'ASP', 'ILE', 'LYS', 'GLN', 'GLY', 'TYR',
         'THR', 'GLN', 'HIS', 'GLY', 'GLY', 'LEU', 'ARG', 'PRO', 'PHE', 'GLY',
         'VAL', 'SER', 'PHE', 'ILE', 'TYR', 'ALA', 'GLY', 'TYR', 'ASP', 'ASP',
         'ARG', 'TYR', 'GLY', 'TYR', 'GLN', 'LEU', 'TYR', 'THR', 'SER', 'ASN',
         'PRO', 'SER', 'GLY', 'ASN', 'TYR', 'THR', 'GLY', 'TRP', 'LYS', 'ALA',
         'ILE', 'SER', 'VAL', 'GLY', 'ALA', 'ASN', 'THR', 'SER', 'ALA', 'ALA',
         'GLN', 'THR', 'LEU', 'LEU', 'GLN', 'MET', 'ASP', 'TYR', 'LYS', 'ASP',
         'ASP', 'MET', 'LYS', 'VAL', 'ASP', 'ASP', 'ALA', 'ILE', 'GLU', 'LEU',
         'ALA', 'LEU', 'LYS', 'THR', 'LEU', 'SER', 'LYS', 'THR', 'THR', 'ASP',
         'SER', 'SER', 'ALA', 'LEU', 'THR', 'TYR', 'ASP', 'ARG', 'LEU', 'GLU',
         'PHE', 'ALA', 'THR', 'ILE', 'ARG', 'LYS', 'GLY', 'ALA', 'ASN', 'ASP',
         'GLY', 'GLU', 'VAL', 'TYR', 'GLN', 'LYS', 'ILE', 'PHE', 'LYS', 'PRO',
         'GLN', 'GLU', 'ILE', 'LYS', 'ASP', 'ILE', 'LEU', 'VAL', 'LYS', 'THR',
         'GLY', 'ILE', 'THR']

if __name__ == "__main__":
  import time
  t0 = time.time()
  for i in range(100):
    res_alignment2(seq_a, seq_b)
  print round(time.time()-t0,4)
  # run6()
  # print '--------------------'
  # print "OK"
  # cProfile.run("res_alignment2(seq_a, seq_b)",filename='cProfile_log')
  # p = pstats.Stats('cProfile_log')
  # p.sort_stats('time').print_stats(15)

