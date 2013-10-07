import os
import iotbx.pdb
import multiprocessing as mp


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
  n_mtrix = 0
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
    except Exception, e: 
      #print "FAILED: %s"%str(e)   
      #print 'BIOMT reading failed on {}'.format(file_name)
      pass
    try: 
      mtrix_obj = pdb_inp.process_mtrix_records() 
      n_mtrix += 1
    except Exception, e: 
      #print "FAILED: %s"%str(e)   
      #print 'MTRIX reading failed on {}'.format(file_name)
      pass
  tnor = len(full_pdb_file_paths())  # tnr: total number of records
  print '========================================================='
  print 'total number of records: {}'.format(tnor)
  print 'number of files that could not be opened:',(tnor - n_pdb_inp)
  print 'number of issues with biomt: ', (n_pdb_inp - n_biomt)
  print 'number of issues with mtrix: ',(n_pdb_inp - n_mtrix)


def proces_file(f,file_name):
  #print "processing:", i_file, file_name
  try: 
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    n_pdb_inp += 1
  except Exception, e: 
    #f.write('pdb_inp reading failed on {}\n'.format(file_name))
    pass
  try: 
    biomt_obj = pdb_inp.process_BIOMT_records()
    
    n_biomt += 1
  except Exception, e:    
    #f.write('BIOMT reading failed on {}'.format(file_name))
    pass
  try: 
    mtrix_obj = pdb_inp.process_mtrix_records() 
    n_mtrix += 1
  except Exception, e:   
    #f.write('MTRIX reading failed on {}'.format(file_name))
    pass
  

if (__name__ == "__main__"):
  while True:
    try:
      np = int(raw_input('Number of processors available is {}. How many would you like to use? '.format(mp.cpu_count())))
      if np>mp.cpu_count():
        raise Exception
      break
    except ValueError:
      print 'Please enter an integer \n'
    except Exception:
      print 'Please choose smaller number of processors \n'
  # set number of CPUs
  p = mp.Pool(processes=np)
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  f = open('zz_out-parallel.txt','w')
  # Start
  print 'Start testing'
  print '============='
  n_pdb_inp = mp.Value('i',0)
  n_biomt = mp.Value('i',0)
  n_mtrix = mp.Value('i',0)
  
  for i_file, file_name in enumerate(full_pdb_file_paths()):  
    print 'ok'
  #run()
  f.close()
