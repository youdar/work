import os
import cPickle as pickle

directory_path = 'C:\Phenix\Dev\Work\work\Clashes'
os.chdir(directory_path)
file_name = 'pdb_clash_score_and_name'
pdb_clash_score_and_name = pickle.load(open(file_name,'r'))
# pdb_clash_score_and_name = list([score_with_hydrogen,score_without_hydrogen,experment_type,file_name]...)
new_file_name = '{}_data.txt'.format(file_name)


directory_path = 'C:\Phenix\Dev\Work\work\Clashes\Data'
os.chdir(directory_path)
f = open(new_file_name,'w')
for i,x in enumerate(pdb_clash_score_and_name):
    if x[2]=='':
        pdb_clash_score_and_name[i][2] = 'ELECTRON MICROSCOPE'
    outstr = '{3},{0:.1f},{1:.1f},{2}\n'.format(x[0],x[1],x[2],x[3])
    f.write(outstr)
f.close()