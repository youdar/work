import os,sys
from libtbx import easy_run
from iotbx.pdb import fetch

def run():
    work_dir = '/net/cci-filer2/raid1/home/youval/Work/work/junk'
    os.chdir(work_dir)
    #
    #file_name = '3zzv'
    #file_name = '4iw4'
    #file_name = '2btv'

    files = set(['3m8l', '2btv', '1vcr', '1llc', '3dpr', '3dar',
                 '2w4y', '2zah', '1tnv', '1dwn', '2w0c', '1w39',
                 '1x33', '2izw', '1x35', '3nop', '2bny', '1ei7',
                 '3not', '3nou', '1dzl', '3lob', '3nap', '1vsz'])


    for fn in files:
        sf_fn = fetch.get_pdb (fn,'xray',mirror='rcsb',log=sys.stdout,format='cif')
        cmd = 'phenix.cif_as_mtz {0} --output-file-name=tmp_file.mtz --merge --remove-systematic-absences --map-to-asu --extend-flags'.format(sf_fn)
        r = easy_run.go(cmd)
        print 'Results for file: {}'.format(fn)
        print '-'*60
        for x in r.stdout_lines:
            print x
        print '='*60


if __name__=='__main__':
    run()
