import cPickle as pickle
import os

os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
print os.getcwd()
structure_factors_files = pickle.load(open('dict_structure_factors_files','r'))
x = '/net/cci/pdb_mirror/structure_factors/'
for k,v in structure_factors_files.items():
    structure_factors_files[k] = x+v


pickle.dump(structure_factors_files,open('dict_structure_factors_files2','w'))
print 'Done'