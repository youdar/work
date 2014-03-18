import os
import cPickle as pickle

def run():
  d = r'c:\Phenix\Dev\Work\work\MTRIX\Data'
  os.chdir(d)

  # set file namme to write results to
  file_out_name = 'NCS_157_files.csv'
  # read needed data
  Collect_tested_files = pickle.load(open('Collect_tested_files','r'))
  Collect_tested_fn = {x[0] for x in Collect_tested_files}
  make_collect_tested_files_dict(Collect_tested_files)
  Collect_tested_files_dict = pickle.load(open('Collect_tested_files_dict','r'))
  files_with_good_MTRIX = set(pickle.load(open('files_with_good_MTRIX','r')))
  file_to_year_dict = pickle.load(open('file_to_year_dict','r'))
  # 157 NCS files
  mtrix_not_included_in_pdb = set(open('mtrix_not_included_in_pdb.txt').read().splitlines())

  # Set of files we are working on:
  work_set = get_worklist(Collect_tested_files_dict,mtrix_not_included_in_pdb)

  write_format_2(work_set,Collect_tested_files_dict,file_out_name)


def write_format_1(work_set,Collect_tested_files_dict,file_out_name):

  # records in Collect_tested_files_dict
  # file_name::r_work_expected::r_work_model_pdb::r_work_model::msg
     #x[0]: file_name: 4 charaters PDB file name
     #x[1]: r_work_expected: r_work from pdb file
     #x[2]: r_work_model_pdb: r_work calulated from pdb
     #x[3]: r_work_model:r_work calculated for complete ASU

  # print to a CSV table with the following 6 fields
  # PDB | Date | Processed |                  R-work values
  #     |      |           | Reported in PDB | Calc from PDB | Calc from reconstructured
  f = open(file_out_name,'w')
  for x in work_set:
    # PDB | Date | Processed | Reported in PDB | Calc from PDB | Calc from reconstructured
    l = [x,file_to_year_dict[x]]
    if Collect_tested_files_dict.has_key(x) and Collect_tested_files_dict[x][1]<100:
      l.append('OK')
      l.extend(map(str,Collect_tested_files_dict[x][1:4]))
    else:
      l.append('--')
      l.extend(['--','--','--'])

    l = ','.join(l) + '\n'
    f.write(l)
    l = []
  f.close()

  print 'results were saved in {}'.format(file_out_name)
  print 'Done...'

def write_format_2(work_set,Collect_tested_files_dict,file_out_name):

  # records in Collect_tested_files_dict
  # file_name::r_work_expected::r_work_model_pdb::r_work_model::msg
     #x[0]: file_name: 4 charaters PDB file name
     #x[1]: r_work_expected: r_work from pdb file
     #x[2]: r_work_model_pdb: r_work calulated from pdb
     #x[3]: r_work_model:r_work calculated for complete ASU

  # print to a CSV table with the following 6 fields
  # PDB | Reported in PDB | Calc from PDB | Calc from reconstructured
  f = open(file_out_name,'w')
  for x in work_set:
    # PDB | Reported in PDB | Calc from PDB | Calc from reconstructured
    l = [x]
    if Collect_tested_files_dict.has_key(x) and Collect_tested_files_dict[x][1]<100:
      l.extend(map(str,Collect_tested_files_dict[x][1:4]))
    else:
      l.extend(['--','--','--'])

    l = ','.join(l) + '\n'
    #print l
    f.write(l)
    l = []
  f.close()

  print 'results were saved in {}'.format(file_out_name)
  print 'Done...'

def get_worklist(Collect_tested_files_dict,work_set=[]):

  a = set(['1pgl', '3m8l', '2vq0', '2btv', '3chx', '4aed', '1pgw', '4fte', '4ang', '1f8v',
           '4hl8', '2wff', '4jgz', '4jgy', '2c4q', '1laj', '1za7', '1b35', '2c50', '3e8k',
           '4iv3', '3zfg', '2ws9', '1lp3', '3bcc', '1rhq', '2qij', '1tdi', '3oah', '2c4y',
           '1f2n', '2bs1', '1k5m', '3ntt', '1uf2', '2gh8', '2bu1', '2qzv', '1c8n', '2zah',
           '3p0s', '3nou', '1r2j', '3vbh', '2c51', '3dpr', '3vbo', '1dzl', '4ftb', '3s6p',
           '3vbf', '3dar', '4fts', '3vbs', '3vbr', '3vbu', '2c4z', '2xpj', '2e0z', '2g34',
           '2iz8', '2buk', '2g33', '2w0c', '1w39', '2x5i', '1tnv', '3qpr', '3fbm', '1js9',
           '2qqp', '2ztn', '4gh4', '1qju', '5msf', '4gb3', '1vcr', '1llc', '1x9t', '1ny7',
           '1x9p', '2gtl', '1vak', '3cji', '1m1c', '4aqq', '2xgk', '1dwn', '1ddl', '2wbh',
           '4ar2', '3s4g', '3ux1', '4gbt', '1h8t', '3vdd', '2bfu', '3r0r', '2vf1', '3nop',
           '2bny', '1ei7', '3not', '2vf9', '3hag', '1bcc', '6msf', '1ng0', '3raa', '3tn9',
           '3n7x', '1wce', '1wcd', '1vb2', '2xbo', '1vb4', '2zzq', '1a37', '1a34', '3es5',
           '2iz9', '2bq5', '4f5x', '1zba', '4g93', '4g0r', '2w4y', '2w4z', '2wws', '3lob',
           '2fz2', '3zfe', '2fz1', '3ra4', '1z7s', '4iv1', '3ra2', '1qjx', '1qjy', '3zff',
           '3ra8', '3ra9', '4gmp', '1x33', '7msf', '1x36', '2izw', '1x35', '1vsz', '3kz4',
           '1ohf', '1ohg', '3nap', '4fsj', '4bcu', '2wzr', '2izn'])

  b = set(['1am7', '3m8l', '2btv', '1fui', '2uu7', '1o7j', '1mpr', '4hmg', '3u1n', '1up6', '1up7', '1up4', '3dar',
           '2uwa', '3zll', '1c5e', '2v4j', '3hmg', '4afl', '1gzh', '4b3z', '4ev0', '2vxi', '1dao', '3cf2', '1upm',
           '1vyw', '1upp', '3sod', '5hmg', '3nop', '1ofs', '2bny', '2wtl', '1ugi', '1swo', '3nou', '3dpr', '3n97',
           '1axg', '1ovo', '1hkd', '1w3m', '2rsl', '5cro', '1hh3', '1dzl', '1o8c', '2ux2', '1b7b', '1hcy', '2vkr',
           '2w0c', '2wl6', '1w39', '1b8d', '1tnv', '1hhf', '1hhc', '1hha', '4bl4', '1hhz', '1ei7', '1hhu', '1vcr',
           '1e3d', '1llc', '1av2', '2boy', '1nal', '2yaw', '1e2t', '3q7h', '2xvt', '1qez', '2a5h', '1kmn', '4ejf',
           '1swb', '1swe', '1swd', '1swg', '1swh', '1swk', '1swj', '1swl', '3not', '1swn', '1swq', '1swp', '1swr',
           '4a63', '2mpr', '4fp4', '2bnl', '2x17', '1a0s', '4hkj', '1z6o', '4bik', '1w7q', '1w7r', '1mms', '2w4y',
           '1bxz', '4j26', '2zah', '1wbi', '2xts', '1b26', '1gy8', '1vzm', '1dwn', '1wc7', '1x33', '4fi4', '2hmg',
           '1x35', '1vsz', '1qox', '2bht', '3nap', '2izw', '1o4z'])

  # sort files according to abs(Collect_tested_files_dict[x][1]-Collect_tested_files_dict[x][3])
  results = []
  tmp = []
  #work_set = a
  for x in work_set:
    if Collect_tested_files_dict.has_key(x):
      results.append([abs(Collect_tested_files_dict[x][1]-Collect_tested_files_dict[x][3]),x])
    else:
      tmp.append(['-',x])


  results.sort()
  results.reverse()
  results.extend(tmp)
  results = [x[1] for x in results]
  print '*'*60
  print 'records to be sorted: {}'.format(len(work_set))
  print 'records sent back: {}'.format(len(results))
  print '*'*60
  return results

def make_collect_tested_files_dict(Collect_tested_files):
  tmp = {x[0]:x for x in Collect_tested_files}
  pickle.dump(tmp, open('Collect_tested_files_dict','w'))



if __name__=='__main__':
  run()

  # for a
  #f = open('NCS_pdb_csv.txt','w')
  # for b
  #f = open('file_with_issues_csv.txt','w')
  #n = 15
  #l = []
  #for i,x in enumerate(b):
    #l.append(x)
    #if (i+1)%n == 0:
      #l = ','.join(l)
      #l += '\n'
      #f.write(l)
      #l = []
  #f.close()