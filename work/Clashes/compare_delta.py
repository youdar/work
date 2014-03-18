from __future__ import division
import Test_internal_clashscore
import os,sys

'''
Collect clash information from PROBE in resraints_manager
and compare them
'''

def get_files_data():
  '''() -> list,list
  reads files
  RM_clash_results
  PROBE_clash_results
  
  in folder: C:\Phenix\Dev\Work\work\Clashes\junk
  
  Returns:
  RM_clash_dict,PROBE_clash_dict : two dictionaries containing the clash 
                                         information from PROBE and resraints_manager
  clash_in_both,clash_only_MR,clash_only_PROBE : sets of clash keys
  '''
  RM_clash_results = open('RM_clash_results','r').read().splitlines()
  PROBE_clash_results = open('PROBE_clash_results','r').read().splitlines()
  #
  RM_clash_dict = {}
  PROBE_clash_dict = {}
  clash_in_both = set()
  clash_only_MR = set()
  clash_only_PROBE = set()
  clash_MR = set()
  clash_PROBE = set()
  #
  for x in RM_clash_results:
    x = x.split('::')
    RM_clash_dict[x[0]] = [float(x[1]),float(x[2])]
    clash_MR.add(x[0])
  for x in PROBE_clash_results:
    x = x.split('::')
    PROBE_clash_dict[x[0]] = float(x[1])
    clash_PROBE.add(x[0])
  #
  clash_in_both = clash_MR.intersection(clash_PROBE)
  clash_only_MR = clash_MR - clash_PROBE
  clash_only_PROBE = clash_PROBE - clash_MR

  return RM_clash_dict,PROBE_clash_dict,clash_in_both,clash_only_MR,clash_only_PROBE


if __name__=='__main__':
  currentpath = os.getcwd()
  workpath = 'c:\\Phenix\\Dev\\Work\\work\\Clashes\\junk'
  os.chdir(workpath)
  #
  file_name = sys.argv[1]
  file_name = Test_internal_clashscore.get_new_file_name(file_name)
  nb_clashscore,clashscore_probe,time_internal,time_probe = Test_internal_clashscore.call_both_clashscores(file_name)
  output_file_name = Test_internal_clashscore.get_file_name(file_name)
  #
  RM_clash_dict,PROBE_clash_dict,clash_in_both,clash_only_MR,clash_only_PROBE = get_files_data()
  # Print clashes that are only in one of the methods:
  print '\nClash info for: {}'.format(output_file_name)
  print '='*80
  print 'nonbonded_clashscore: {0:.3f}'.format(nb_clashscore[1])
  print 'clashscore          : {0:.3f}'.format(clashscore_probe)
  print '='*80
  print 'Clashes that show up only in PROBE'
  print '-'*80
  for rec in clash_only_PROBE: 
    print '{0:30}{1:^14.3f}'.format(rec,PROBE_clash_dict[rec])
  print '='*80
  print 'Clashes that show up only in restraints_manager'
  print 'Note: those clashes do not include clashes due to symmetry operations'
  print '-'*80
  for rec in clash_only_MR: 
    print '{0:30}{1:^14.3f}'.format(rec,RM_clash_dict[rec][0])
  print '='*80
  #
  print 'Clashes in both'
  outstr =  '{0:30}{1:^14.3f}{2:^14.3f}{3:^14.3f}{4:^14.3f}'
  print '{0:30}{1:^14}{2:^14}{3:^14}{4:^14}'.format('Clash','overlap RM','overlap RPROBE','diff','vdw')
  print '-'*80
  for rec in clash_in_both:
    overlap_RM = RM_clash_dict[rec][0]
    vdw_RM = RM_clash_dict[rec][1]
    overlap_PROBE = PROBE_clash_dict[rec]
    print outstr.format(rec,overlap_RM,overlap_PROBE,overlap_RM-overlap_PROBE,vdw_RM)
  print '='*80
  #
  os.chdir(currentpath)
  print 'Done'