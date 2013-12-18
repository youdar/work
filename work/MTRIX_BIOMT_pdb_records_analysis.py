import os
import iotbx.pdb

'''
Collect information on MTRIX and BIOMT records
in the pdb
'''

def full_pdb_file_paths():
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").readlines()
  result = []
  for pdb_file in pdb_files:
    #pdb_code = pdb_file[-12:-8]
    pdb_file = os.path.join(pdb_dir, pdb_file.strip())
    # result is a list of records like: '/net/chevy/raid1/pdb_mirror/pdb/00/pdb100d.ent.gz'
    result.append(pdb_file)
  return result

def run():
  print 'Start testing'
  print '============='
  n_pdb_inp = 0
  n_biomt = 0
  empty_biomt = 0
  n_mtrix = 0
  empty_mtrix = 0
  good_biomt = set()
  good_mtrix = set()
  odd_biomt = set()
  odd_mtrix = set()
  eye_biomt = set()
  eye_mtrix = set()
  # set identity matrix for testing
  eye = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
  # files with BIOMT records that are larger than 1
  f1 = open('biomt_ok.txt','w')
  # files with MTRIX info that is larger than 1
  f2 = open('mtrix_ok.txt','w')
  # files with both MTRIX and BIOMT good info
  f3 = open('biomt_mtrix_ok.txt','w')
  # BIOMT with length 1 if the 1 is not the identity matrix
  f4 = open('biomt_odd.txt','w')
  # MTRIX with length 1 if the 1 is not the identity matrix
  f5 = open('mtrix_odd.txt','w')
  # Goob biomt but no mtrix
  f6 = open('good_bio_bad_mt.txt','w')
  # Good mtrix but no biomt
  f7 = open('bad_bio_good_mt.txt','w')
  f8 = open('eye_biomt.txt','w')
  f9 = open('eye_mtrix.txt','w')
  # start going over all pdb files
  for i_file, file_name in enumerate(full_pdb_file_paths()):
    #print "processing:", i_file, file_name
    try:
      pdb_inp = iotbx.pdb.input(file_name=file_name)
      n_pdb_inp += 1
    except Exception, e:
      #print "FAILED: %s"%str(e)
      #print 'pdb_inp reading failed on {}'.format(file_name)
      pass
    try:
      biomt_obj = pdb_inp.process_BIOMT_records()
      n_biomt += 1
      if biomt_obj==[]:
        empty_biomt += 1
      elif len(biomt_obj) == 1:
        if (biomt_obj[0].values[0] != eye) and (not mtrix_obj[0].coordinates_present):
          odd_biomt.add(file_name)
        else:
          eye_biomt.add(file_name)
      else:
        f1.write(file_name + '\n')
        good_biomt.add(file_name)
    except Exception, e:
      #print "FAILED: %s"%str(e)
      #print 'BIOMT reading failed on {}'.format(file_name)
      pass
    try:
      mtrix_obj = pdb_inp.process_mtrix_records()
      n_mtrix += 1
      if mtrix_obj==[]:
        empty_mtrix += 1
      elif len(mtrix_obj) == 1:
        if (mtrix_obj[0].values[0] != eye) and (not mtrix_obj[0].coordinates_present):
          odd_mtrix.add(file_name)
        else:
          eye_mtrix.add(file_name)
      else:
        f2.write(file_name + '\n')
        good_mtrix.add(file_name)
    except Exception, e:
      #print "FAILED: %s"%str(e)
      #print 'MTRIX reading failed on {}'.format(file_name)
      pass
  tnor = len(full_pdb_file_paths())  # tnr: total number of records
  # Intersect good biomt and mtrix
  both_good = good_biomt & good_mtrix
  f3.writelines([x+'\n' for x in both_good])
  f4.writelines([x+'\n' for x in odd_biomt])
  f5.writelines([x+'\n' for x in odd_mtrix])
  f6.writelines([x+'\n' for x in (good_biomt - good_mtrix)])
  f7.writelines([x+'\n' for x in (good_mtrix - good_biomt)])
  f8.writelines([x+'\n' for x in eye_biomt])
  f9.writelines([x+'\n' for x in eye_mtrix])
  print '========================================================='
  print 'total number of records: {}'.format(tnor)
  print 'number of files that could not be opened:',(tnor - n_pdb_inp)
  print '---------------------------------------------------------'
  print 'number of BIOMT records that can be read (including empty and identity): ', n_biomt
  print 'number of reading issues with biomt: ', (n_pdb_inp - n_biomt)
  print 'number of empty with biomt: ', empty_biomt
  print 'BIOMT with a single record which is not the identitiy and not present: ',len(odd_biomt)
  print 'number of identity BIOMT: ',len(eye_biomt)
  print 'number of good BIOMT data: ', good_biomt
  print '---------------------------------------------------------'
  print 'number of MTRIX records that can be read (including empty and identity): ', n_mtrix
  print 'number of reading issues with mtrix: ', (n_pdb_inp - n_mtrix)
  print 'number of empty with mtrix: ', empty_mtrix
  print 'MTRIX with a single record which is not the identitiy and not present: ',len(odd_mtrix)
  print 'number of identity MTRIX: ',len(eye_mtrix)
  print 'number of goos MTRIX data: ', good_mtrix
  print '---------------------------------------------------------'
  print 'number of good BIOMT and MTRIX: ', len(both_good)
  print 'Good BIOMT and bad MTRIX: ', len(good_biomt - good_mtrix)
  print 'Dab BIOMT and GOOD MTRIX: ', len(good_mtrix - good_biomt)
  print '========================================================='
  f1.close()
  f2.close()
  f3.close()
  f4.close()
  f5.close()
  f6.close()
  f7.close()
  f8.close()
  f9.close()



if __name__ == "__main__":
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  run()
