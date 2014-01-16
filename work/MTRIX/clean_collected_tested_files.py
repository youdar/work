import os
import cPickle as pickle

def run():
    print os.getcwd()
    data_dir = '/net/cci-filer2/raid1/home/youval/Work/work/MTRIX/Data'
    #data_dir = r'c:\Phenix\Dev\Work\work\MTRIX\Data'
    results = []
    # all 157 files with good MTRIX also have structure factors files
    ## good MTRIX is tested with eps = 0.01
    #files_with_good_MTRIX = set(pickle.load(open(os.path.join(data_dir,'files_with_good_MTRIX'),'r')))
    ## file location dictionary
    #good_MTRIX_pdb_files = pickle.load(open(os.path.join(data_dir,'dict_good_MTRIX_pdb_files'),'r'))
    #structure_factors_files = pickle.load(open(os.path.join(data_dir,'dict_structure_factors_files'),'r'))
    #print 'File names and location dictionaries are loaded...'

    ## look at files that were processed
    #files_with_problems = pickle.load(open(os.path.join(data_dir,'files_with_problems'),"r"))
    #files_with_problems = {x[0] for x in files_with_problems}

    Collect_tested_files = pickle.load(open(os.path.join(data_dir,'Collect_tested_files'),"r"))



    # clean errors
    for i,val in enumerate(Collect_tested_files):
        for j,rec in enumerate(val[4:]):
            if 'determine output label' in rec:
                val[j+4] = 'OK'
        Collect_tested_files[i] = val

    print 'processed {} files'.format(len(Collect_tested_files))
    Collect_tested_files_set = {x[0] for x in Collect_tested_files}
    print 'Number of unique files in Collect_tested_files: {}'.format(len(Collect_tested_files_set))

    while 1:
        file_name_list = [x[0] for x in Collect_tested_files]
        multiple_records = set()
        for fn in Collect_tested_files_set:
            if file_name_list.count(fn) > 1:
                multiple_records.add(fn)
        # stop if no multiple_records
        if not multiple_records:
            break
        else:
            print 'need to remove {} duplicates cycle'.format(len(multiple_records))
        # remove duplicates
        for fn in multiple_records:
            first_rec = file_name_list.index(fn)
            p = Collect_tested_files.pop(first_rec)
            print p
            file_name_list.pop(first_rec)
    # clean repeated OK messege
    for i,fn in enumerate(Collect_tested_files):
        while fn.count('OK') > 1:
            ind = fn.index('OK')
            fn.pop(ind)
            Collect_tested_files[i] = fn


    print 'After clean up'
    print 'Numebr of records in Collected_tested_files: {}'.format(len(Collect_tested_files))
    print 'Number of unique files in Collect_tested_files: {}'.format(len({x[0] for x in Collect_tested_files}))

    # Saving the clean file
    #pickle.dump(Collect_tested_files, open(os.path.join(data_dir,'Collect_tested_files'),'w'))

    print 'Done...'

if __name__=='__main__':
    run()