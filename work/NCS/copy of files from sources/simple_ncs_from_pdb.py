"""
  simple_ncs_from_pdb.py
  tct 2006-12-12

**************************************
To use as a method, specify input PDB and any other commands as text arguments:

ncs_search=simple_ncs_from_pdb(args=['input_pdb.pdb','max_rmsd=4.'])

Alternatively, you can specify inputs with a params phil object:
ncs_search=simple_ncs_from_pdb(params=params)

Now ncs_search.ncs_object is an object from "sources/mmtbx/mmtbx/ncs.py"
You can get a summary or text for resolve or phenix.refine with:

ncs_search.ncs_object.display_all()
text=ncs_search.ncs_object.format_all_for_resolve()
text=ncs_search.ncs_object.format_all_for_phenix_refine()
**************************************

Purpose: Figure out the ncs from the chains in this pdb object
and return an "ncs" object from "ncs.py" that represents it. This ncs object
can write out a little file for resolve and a file for phenix.refine
specifying this ncs.

Major assumption: residue numbers are consistent among chains

Approach: Use residue numbers to align the residue names, identify
pairs of chains that can match. Choose groupings of chains that maximize the
smallest number of matching residues between each member of a group and the
first (reference) member of the group. Within a pair of chains, allow some
segments to match and others not. Each pair of segments must have a
length >= min_length and percent identity >=min_percent.  A pair of segments
may not end in a mismatch. An overall pair of chains must have an rmsd
of CA atoms of <= rmsd_max.

If find_invariant_domain is specified then once all chains that can be matched
with the above algorithm are identified, all remaining chains are matched,
allowing the break-up of chains into invariant domains. The invariant
domains each get a separate NCS group.

If residue numbers are not the same for corresponding chains, but
they are simply offset by a constant for each chain, this will be
recognized and the chains will be aligned.

"""
from __future__ import division
import iotbx
import iotbx.phil
from iotbx import pdb
from mmtbx.ncs.ncs_from_pdb import ncs_from_pdb
import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry,null_out
from libtbx import adopt_init_args
from mmtbx.ncs.ncs import ncs
import sys, os, string
from cctbx.array_family import flex
import mmtbx.monomer_library.pdb_interpretation
from phenix.utilities.arg_display_methods import arg_display_methods
from phenix.utilities.composite_params import get_composite_master_params
from phenix.utilities.list_methods import list_methods
from phenix.utilities.headers import headers
from phenix.utilities.is_debug import is_debug
from phenix.utilities import citations
import iotbx.ncs

from phenix.utilities.catenate_equals import catenate_equals

ncs_master_params = iotbx.phil.parse("""
simple_ncs_from_pdb
  .short_caption = Simple NCS from PDB file
{

 pdb_in = None
   .type = path
   .help = 'Input PDB file to be used to identify ncs'
   .short_caption = PDB file
   .style = bold noauto OnChange:load_ncs_pdb_file file_type:pdb no_map
 temp_dir = ""
   .type = path
   .help = "temporary directory (ncs_domain_pdb will be written there)"
   .style = noauto
 min_length= 10
   .help = "minimum number of matching residues in a segment"
   .type = int
    .short_caption = Min. number of matching residues per segment
   .expert_level = 2
 min_fraction_represented = 0.10
   .help = "Minimum fraction of residues represented by NCS to keep."
           "If less...skip ncs entirely"
   .type = float
    .short_caption = Min. fraction represented
   .expert_level = 2
 njump   = 1
   .help = "Take every njumpth residue instead of each 1"
   .type = int
   .short_caption = Number of residues per jump
    .expert_level = 2
 njump_recursion   = 10
   .help = "Take every njump_recursion residue instead of each 1 on recursive call"
   .type = int
   .expert_level = 2
 min_length_recursion = 50
   .help = "minimum number of matching residues in a segment for recursive call"
   .type = int
   .short_caption = Min. length recursion
   .expert_level = 2
 min_percent= 80.
   .help = "min percent identity of matching residues"
   .type = float
   .short_caption = Min. percent identity
   .expert_level = 0
   .style = bold noauto
 max_rmsd = 2.
   .help = "max rmsd of 2 chains. If 0, then only search for domains"
   .type = float
    .short_caption = Max. RMSD
   .expert_level = 0
   .style = bold noauto
 quick = True
   .type = bool
   .help = "If quick is set and all chains match, just look for 1 NCS group"
   .short_caption = Quick search
   .style = noauto
 max_rmsd_user = 3.
   .help=max rmsd of chains suggested by user (i.e., if called from phenix.refine \
         with suggested ncs groups)

   .type = float
   .short_caption = Max. RMSD for user-specified chains
   .expert_level = 2
 maximize_size_of_groups = True
   .type = bool
   .help = '''You can request that the scoring be set up to maximize
       the number of members in NCS groups (maximize_size_of_groups=True)
       or that scoring is set up to maximize the length of the matching
       segments in the NCS group (maximize_size_of_groups=False)'''
 require_equal_start_match = True
   .type = bool
   .help = '''You can require that all matching segments start at the same
           relative residue number for all members of an NCS group,
           trimming the matching region as necessary. This
           is required if residue numbers in different chains are not the
           same, but not otherwise'''
 ncs_domain_pdb_stem  = None
   .type = str
   .help = '''NCS domains will be written to ncs_domain_pdb_stem+"group_"+nn'''
   .style = noauto
 write_ncs_domain_pdb = False
   .type = bool
   .help = '''You can write out PDB files representing NCS domains for
        density modification if you want'''
   .style = noauto
 domain_finding_parameters
    .style = box auto_align noauto
  {
   find_invariant_domains = True
   .type = bool
   .help = "Find the parts of a set of chains that follow NCS"
   initial_rms = 0.5
   .type=float
   .help="Guess of RMS among chains"
   match_radius = 4.0
   .type = float
   .help = '''max allow distance difference between pairs of matching
        atoms of two residues'''
   similarity_threshold = 0.95
   .type=float
   .help="Threshold for similarity between segments"
   smooth_length = 0
   .help = "two segments separated by smooth_length or less get connected"
   .type=int
   min_contig_length = 3
   .help = "segments < min_contig_length rejected"
   .type=int
   min_fraction_domain = 0.2
     .help = "domain must be this fraction of a chain"
     .type = float
   max_rmsd_domain = 2.
     .help = "max rmsd of domains"
     .type = float
 }
 verbose = False
   .type = bool
   .help = "Verbose output"
   .short_caption = Debugging output
    .style = noauto
 raise_sorry = False
   .type = bool
   .help = "Raise sorry if problems"
   .short_caption = raise sorry

 debug = False
   .type = bool
   .help = "Debugging output"
   .short_caption = Debugging output
    .style = noauto
 dry_run = False
   .type = bool
   .help = '''Just read in and check parameter names'''
    .style = noauto
    }
""")
master_params = ncs_master_params

restraint_group_params = iotbx.phil.parse("""
 restraint_group
  .multiple=True
  .optional=True
  .short_caption=Restraint group
  .style = noauto auto_align
{
  reference=None
    .type=atom_selection
    .optional=True
    .short_caption=Reference selection
    .input_size=400
    .style = bold
  selection=None
    .type=atom_selection
    .multiple=True
    .short_caption=Restrained selection
    .input_size=400
    .style = bold
  coordinate_sigma=0.05
    .type=float
  b_factor_weight=10
    .type=float
}
""")

class simple_ncs_from_pdb(arg_display_methods,list_methods,headers):

  def __init__(self, args   = None,
                     params = None,
                     ignore_chains   = [],
                     required_chains = [],
                     exclude_chains = [],
                     ncs_master_params = ncs_master_params,
                     command_name   = "simple_ncs_from_pdb",
                     all_chain_proxies = None,
                     pdb_inp = None,
                     hierarchy = None,
                     suppress_print = False,
                     source_info    = None,
                     pdb_file       = None,
                     groups_only = None,
                     suggested_ncs_groups = None,
                     log = sys.stdout,
                     quiet=False,
                     exclude_h=False,
                     exclude_d=False,
                     write_ncs_domain_pdb=None,
                     ncs_domain_pdb_stem=None,
                     temp_dir=None
               ):
    self.log=log
    self.quiet=quiet
    self.exclude_h=exclude_h
    self.exclude_d=exclude_d
    self.exclude_chains=exclude_chains
    self.required_chains=required_chains
    args=catenate_equals(args).new_args()
    self.process_inputs(args)
    if self.args==[]:
       self.args=None
    args=self.args
    if suggested_ncs_groups is None: # take it from args if not directly given
      suggested_ncs_groups=self.suggested_ncs_groups
    else:
      self.suggested_ncs_groups=suggested_ncs_groups
    allow_recursion=self.allow_recursion
    exact_match_only=self.exact_match_only

    master_params_name='simple_ncs_from_pdb'
    self.Name='simple_ncs_from_pdb'


    if not suppress_print:
      citations.add_citation('phenix','simple_ncs_from_pdb')

    # run a test case
    if args is not None and 'exercise' in args:
      self.exercise()
      return

    args=self.special_cases(args)

    master_params=get_composite_master_params(
         method_list=['simple_ncs_from_pdb'],
         location_list=['phenix.command_line'])

    args=self.get_keyword_table(args,out=self.log)       # set self.keyword_table

    summary,header=self.get_summary_and_header(command_name)

    done,master_params,new_params,changed_params,help=self.get_params(
        command_name,master_params,args,out=self.log)
    if params is None :
      params = new_params
    if done: return
    if not quiet: print >>self.log, header

    if help or (params and params.simple_ncs_from_pdb.verbose):
      print >>self.log, "Values of all params:"
      master_params.format(python_object=params).show(out=log)

    if help or params is None: return

    # Done with standard processing of inputs
    # overwrite with direct inputs, if any:
    if write_ncs_domain_pdb is not None:
      params.simple_ncs_from_pdb.write_ncs_domain_pdb=write_ncs_domain_pdb
    if ncs_domain_pdb_stem is not None:
      params.simple_ncs_from_pdb.ncs_domain_pdb_stem=ncs_domain_pdb_stem
    if temp_dir is not None:
      params.simple_ncs_from_pdb.temp_dir=temp_dir

    # Things that must be defined...

    self.params=params
    if not suppress_print:
      print >>self.log,"Parameters used for simple_ncs_from_pdb:"
      master_params.format(python_object=params).show(out=self.log)
      print >>self.log

    if params.simple_ncs_from_pdb.dry_run:
      print "ARGS: ",args
      return

    # read in the PDB file if needed
    if((all_chain_proxies is None) and (pdb_inp is None and hierarchy is None)):
      if(pdb_file is None):
        if args is not None and args and args[0] and os.path.isfile(args[0]):
          pdb_file=args[0]
        elif params.simple_ncs_from_pdb.pdb_in is not None:
          pdb_file=params.simple_ncs_from_pdb.pdb_in
        else:
          raise Sorry("\nNeed PDB file for simple_ncs_from_pdb"+
             "\n\nPlease make the PDB file the first argument like this: \n"+
             "phenix.simple_ncs_from_pdb mypdb.pdb ...\n")
      if not os.path.isfile(pdb_file):
         raise Sorry("The file "+str(pdb_file)+" is missing?")

      raw_records = flex.std_string()
      raw_records.extend(flex.split_lines(open(pdb_file).read()))
      if pdb_inp is None:
        pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records)
      mon_lib_srv = mmtbx.monomer_library.server.server()
      ener_lib = mmtbx.monomer_library.server.ener_lib()
      processed_pdb= mmtbx.monomer_library.pdb_interpretation.process(
         mon_lib_srv=mon_lib_srv,
         ener_lib=ener_lib,
         params=None,
         raw_records=raw_records,
         strict_conflict_handling=False,
         max_atoms=None,
         log=null_out())

      all_chain_proxies=processed_pdb.all_chain_proxies
      self.source_info=pdb_file
    if hierarchy is None:
      hierarchy = pdb_inp.construct_hierarchy()

    # set input params
    self.verbose=params.simple_ncs_from_pdb.verbose
    if self.verbose: hierarchy.show(out=self.log)
    if source_info:
      self.source_info=source_info
    if not hasattr(self,'source_info') or not self.source_info:
      self.source_info="None"

    ####

    # self.find_invariant_domains=\
    #    params.simple_ncs_from_pdb.domain_finding_parameters.find_invariant_domains
    # self.min_fraction_domain= \
    #    params.simple_ncs_from_pdb.domain_finding_parameters.min_fraction_domain
    # self.max_rmsd_domain=params.simple_ncs_from_pdb.domain_finding_parameters.max_rmsd_domain
    # self.min_percent=params.simple_ncs_from_pdb.min_percent
    # self.quick=params.simple_ncs_from_pdb.quick
    # self.maximize_size_of_groups=params.simple_ncs_from_pdb.maximize_size_of_groups
    # self.require_equal_start_match=params.simple_ncs_from_pdb.require_equal_start_match
    # self.write_ncs_domain_pdb=params.simple_ncs_from_pdb.write_ncs_domain_pdb
    # self.ncs_domain_pdb_stem=params.simple_ncs_from_pdb.ncs_domain_pdb_stem
    # self.all_chain_proxies=all_chain_proxies
    # self.hierarchy=hierarchy
    # self.pdb_inp=pdb_inp
    # self.njump=params.simple_ncs_from_pdb.njump
    # if self.njump<1:
    #   raise Sorry("njump must be >=1")
    # self.min_length=params.simple_ncs_from_pdb.min_length
    # self.min_fraction_represented=params.simple_ncs_from_pdb.min_fraction_represented

    # Use cctbx tool to get ncs object
    find_param = params.simple_ncs_from_pdb.domain_finding_parameters
    inp_param = params.simple_ncs_from_pdb
    ncs_obj = iotbx.ncs.input(
      hierarchy=hierarchy,
      chain_similarity_limit=find_param.similarity_threshold,
      min_contig_length=find_param.min_contig_length,
      min_percent=inp_param.min_percent,
      max_rmsd=inp_param.max_rmsd,
      max_dist_diff=find_param.match_radius,
      use_minimal_master_ncs=True,
      process_similar_chains=True,
      allow_different_size_res=True,
      exclude_misaligned_residues=True,
      check_atom_order=False,
      write_messages=False,
      log=self.log,
      quiet=self.quiet)
    spec_object = ncs_obj.get_ncs_info_as_spec(write=False)

    # # identify chains in the PDB file
    # find_param = self.params.simple_ncs_from_pdb.domain_finding_parameters
    # ncs_process = ncs_from_pdb(
    #   verbose=self.verbose,
    #   log=self.log,
    #   njump=self.njump,
    #   min_length=self.min_length,
    #   min_percent=self.min_percent,
    #   suggested_ncs_groups=self.suggested_ncs_groups,
    #   require_equal_start_match=self.require_equal_start_match,
    #   maximize_size_of_groups=self.maximize_size_of_groups,
    #   required_chains=self.required_chains,
    #   min_fraction_domain=self.min_fraction_domain,
    #   initial_rms=find_param.initial_rms,
    #   match_radius =find_param.match_radius,
    #   similarity_threshold=find_param.similarity_threshold,
    #   min_contig_length=find_param.min_contig_length,
    #   max_rmsd_domain=find_param.max_rmsd_domain,
    #   min_fraction_represented=self.min_fraction_represented)
    #
    # chains,chain_ids,starting_residue_numbers,offset_dict, total_residues= \
    #   ncs_process.get_chain_list(
    #     hierarchy=hierarchy,ignore_chains=self.ignore_chains)

    # self.total_residues = total_residues
    # # so we do not have to pass this around. It does not change ever.
    # self.offset_dict=offset_dict

    # if not suppress_print:
    #   print >>self.log,"Chains in this PDB file: ",chain_ids
    # if self.verbose:
    #  for chain,chain_id,start in zip(chains,chain_ids,starting_residue_numbers):
    #    print >>self.log,"CHAIN: ",chain," ID: ",chain_id," START" ,start
    #  print >>self.log,"OFFSET LIST: "
    #  for id in self.offset_dict.keys():
    #    print >>self.log,id,self.offset_dict[id]


    # # set up temp_dir if needed
    # if params.simple_ncs_from_pdb.temp_dir:
    #   if os.path.isfile(params.simple_ncs_from_pdb.temp_dir):
    #     raise Sorry(
    #      "The directory "+str(params.simple_ncs_from_pdb.temp_dir)+" cannot be created...")
    #   if not os.path.isdir(params.simple_ncs_from_pdb.temp_dir):
    #     os.mkdir(params.simple_ncs_from_pdb.temp_dir)
    #   if not self.quiet:
    #     print >>self.log,"Working in ",params.simple_ncs_from_pdb.temp_dir

    # groups=[]  # [ ['A','B'],['C','D','E]]
    # #              group
    # list_of_residue_range_list=[] # [ [      [       [1,120],[130-250] ]]]
    # #                                 group  member  ranges
    # [all_chains,all_chain_ids,all_starting_residue_numbers]=[chains,chain_ids,
    #   starting_residue_numbers]


    # #========Get suggested NCS groups from phil object ======================
    #
    # #  If called with suggested_ncs_groups phil object, we pull out all those
    # #  chains here...
    # # initialize suggested NCS groups if any
    # ncs_process.suggested_ncs_groups,self.suggested_group_list=\
    #       ncs_process.get_suggested_groups(suggested_ncs_groups,chain_ids)
    #
    # # NOTE: residue_range_list is in original residue numbers, with offsets
    # used_ids=[]
    # if self.suggested_group_list:  # see if we want to pull out these chains:
    #   print >> self.log, "Getting NCS from suggested chains: "
    #   ncs_process.suggested_ncs_groups=[] # not both
    #   print '(0)'
    #   for [group,residue_range_list] in self.suggested_group_list:
    #     rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
    #       ncs_process.get_rmsd(group,hierarchy,chains,chain_ids,
    #          starting_residue_numbers,residue_range_list)
    #     if not rmsd_list or len(rmsd_list)<2:
    #       pass
    #     elif rmsd_list[1]>params.simple_ncs_from_pdb.max_rmsd_user:
    #       print >>self.log,"Warning: requested alignment of ",group,\
    #            ncs_process.add_offsets(residue_range_list,group),\
    #            " \nrejected due \nto rmsd =",rmsd_list[1]," > ",\
    #              params.simple_ncs_from_pdb.max_rmsd_user,\
    #        ". To keep it, set ncs.max_rmsd_user="+str(int(rmsd_list[1]+0.999))
    #     else:
    #       groups.append(group)
    #       list_of_residue_range_list.append(residue_range_list)
    #       print >>self.log,"RMSD for suggested group ",group,\
    #         ncs_process.add_offsets(residue_range_list,group)," is ",rmsd_list[1]
    #   # Now remove all chains in kept groups from list of chains so we only
    #   # look elsewhere
    #   used_ids=ncs_process.add_ids(groups,used_ids)
    #   [chains,chain_ids,starting_residue_numbers]=\
    #     ncs_process.remove_used_chains(
    #       chains,chain_ids,starting_residue_numbers,used_ids)
    #========End of suggested NCS groups from phil object ==============

    #======= Try to get possible NCS groups from direct comparison of ====
    #         all pairs of chains. Use a high njump to go quickly...
    #         This will work if all chains in a group are identical

    # if allow_recursion:
    #   njump_use=params.simple_ncs_from_pdb.njump_recursion
    #   min_length_use=params.simple_ncs_from_pdb.min_length_recursion
    #   if args is not None:
    #    args_use=args
    #   else:
    #    args_use=[]
    #   args_use+=['njump='+str(njump_use),'min_length='+str(min_length_use)]
    #   if ncs_process.suggested_ncs_groups:   # present if NCS groups defined as "ACDE"
    #     args_use.append("suggested_ncs_groups="+str(suggested_ncs_groups))
    #   if self.verbose:
    #     logfile=log
    #   else:
    #     logfile=null_out()
    #   if self.exact_match_only:
    #      args_use.append('exact_match')
    #   args_use.append('no_recursion')
    #
    #   quick_find_groups=simple_ncs_from_pdb (
    #                  args   = args_use,
    #                    # at end of arg list to overwrite..
    #                  pdb_inp = pdb_inp,  # 091608
    #                  params = params,
    #                  ignore_chains   = used_ids,
    #                  all_chain_proxies = all_chain_proxies,
    #                  hierarchy = hierarchy ,
    #                  suppress_print = True,
    #                  source_info    = source_info,
    #                  log=logfile,
    #                  quiet=True,
    #                  groups_only = True)
    #
    #   ncs_process.suggested_ncs_groups=quick_find_groups.sequence_groups
    #   if not suppress_print:
    #     print >>self.log,"GROUPS BASED ON QUICK COMPARISON:",\
    #       quick_find_groups.sequence_groups
    #
    # # ======= End of getting NCS groups from direct comparison =========
    #
    # # ======= Getting groups quickly if this is called by simple_ncs ======
    # #         This does the work for the recursive call above
    # if groups_only:
    #   self.sequence_groups,self.sequence_list_of_residue_range_list=\
    #      ncs_process.find_groups(hierarchy,chains,chain_ids,
    #        starting_residue_numbers,
    #        min_length=params.simple_ncs_from_pdb.min_length,
    #        min_percent=params.simple_ncs_from_pdb.min_percent,
    #        max_rmsd=params.simple_ncs_from_pdb.max_rmsd,
    #        max_rmsd_user=params.simple_ncs_from_pdb.max_rmsd_user,
    #        called_by_self=True,
    #        exact_match_only=exact_match_only)
    #   return
    # ======= End of getting groups quickly if this is called by simple_ncs ======


    # # ====== Get NCS groups using only sequence information and the
    # #        suggested_ncs_groups found above, if any
    # sequence_groups,sequence_list_of_residue_range_list=\
    #      ncs_process.find_groups(hierarchy,chains,chain_ids,
    #        starting_residue_numbers,
    #        min_length=params.simple_ncs_from_pdb.min_length,
    #        min_percent=params.simple_ncs_from_pdb.min_percent,
    #        max_rmsd=999999.,max_rmsd_user=15,called_by_self=True,
    #        exact_match_only=exact_match_only)
    # if self.verbose:
    #   print >>self.log,"SEQUENCE-BASED GROUPS: ",sequence_groups
    #   print >>self.log,sequence_list_of_residue_range_list
    # ====== End of NCS groups using sequence information and suggested groups=======

    # JUST HERE problem sequence_list_of_residue_range_list has offsets in
    # some cases where it should not  021107...

    # # =======Get new groups with domains on all sequence groups=======
    # invariant_groups,invariant_list_of_residue_range_list=\
    #         ncs_process.find_invariant_groups(hierarchy,
    #          sequence_groups,sequence_list_of_residue_range_list,
    #          chains,chain_ids,
    #          starting_residue_numbers)
    # if invariant_groups:
    #     if self.verbose:
    #       print >>self.log,"Invariant groups found:",\
    #           invariant_groups,invariant_list_of_residue_range_list
    #     groups+=invariant_groups
    #     list_of_residue_range_list+=invariant_list_of_residue_range_list

    # =======End of new groups with domains =============

    # [chains,chain_ids,starting_residue_numbers]=[all_chains,  # restore these
    #     all_chain_ids,all_starting_residue_numbers]

    #  ======== Write out results and make an ncs object with them in it===
    # ncs_object=ncs(exclude_h=self.exclude_h,exclude_d=self.exclude_d)
    # count=0
    # for group,residue_range_list in zip(groups,list_of_residue_range_list):
    #   count+=1
    #   # if necessary, add offsets from self.offset_dict to the values in
    #   #  residue_range_list
    #   residue_range_list_with_offsets=ncs_process.add_offsets(residue_range_list,group)
    #   if self.verbose:
    #     print >>self.log,"\nNCS GROUP",count,":",group  ,\
    #         residue_range_list_with_offsets
    #   # group is a list of chain_ids, with the reference one first
    #   # so get rmsd for members of the group from reference
    #   rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
    #     ncs_process.get_rmsd(group,hierarchy,chains,chain_ids,starting_residue_numbers,
    #     residue_range_list)  # NO OFFSET (I know, it's confusing)!
    #
    #   if self.write_ncs_domain_pdb:
    #     ncs_domain_pdb=self.make_ncs_domain_pdb(
    #      stem=self.ncs_domain_pdb_stem,
    #      hierarchy=hierarchy,group_number=count,
    #      group=group,residue_range_list=residue_range_list_with_offsets,
    #      params=params)
    #   else:
    #     ncs_domain_pdb=None
    #   if not rmsd_list:
    #      print >>self.log,"\nNCS GROUP",count,":",group  ,\
    #            residue_range_list_with_offsets
    #      print >>self.log,"No rmsd found...giving up on this group"
    #   else:
    #     chain_residue_id=[group,residue_range_list_with_offsets]
    #     ncs_object.import_ncs_group(ncs_rota_matr=r_list,
    #      rmsd_list=rmsd_list,
    #      residues_in_common_list=residues_in_common_list,
    #      center_orth=center_list,
    #      trans_orth=trans_list,
    #      chain_residue_id=chain_residue_id,
    #      ncs_domain_pdb=ncs_domain_pdb,
    #      source_of_ncs_info=self.source_info)

    ncs_object = spec_object

    # if len(ncs_object.ncs_groups()) >=1:
      # and ncs_process.too_few_residues_represented(ncs_object=ncs_object,
      # total_residues=self.total_residues): # skip entirely
      # print >>self.log,"Skipping NCS. Too few residues represented (< %6.1f percent of total)" %(100.*self.min_fraction_represented)
      # ncs_object=ncs(exclude_h=self.exclude_h,exclude_d=self.exclude_d)

    if len(ncs_object.ncs_groups())<1:
      if not suppress_print:
        if self.source_info:
          print >>self.log,"\nNo NCS found from the chains in ",self.source_info
        else:
          print >>self.log,"\nNo NCS found"

    if not suppress_print:
      ncs_object.display_all(log=self.log)
      f=open("simple_ncs_from_pdb.resolve",'w')
      ncs_object.format_all_for_resolve(log=self.log,out=f)
      f.close()
      f=open("simple_ncs_from_pdb.ncs",'w')
      ncs_object.format_all_for_phenix_refine(log=self.log,out=f)
      f.close()
      f=open("simple_ncs_from_pdb.ncs_spec",'w')
      ncs_object.format_all_for_group_specification(log=self.log,out=f)
      f.close()

    self.ncs_object=ncs_object

    citations.show(out=self.log, source='simple_ncs_from_pdb')

    #  ======== Done with writing out results and making an ncs object ====

  def make_ncs_domain_pdb(self,stem='',hierarchy=None,group_number=None,
         group=[],residue_range_list=[],params=None):
    # write out a file with the name of the group containing all atoms
    # associated with this group
    #if self.exclude_chains is set then EXCLUDE members of exclude_chains

    if hierarchy is None or group_number is None: return None
    if stem is None: stem=""
    file_name=stem+'group_'+str(group_number)+'.pdb'
    full_file_name=os.path.join(params.simple_ncs_from_pdb.temp_dir,file_name)
    f=open(full_file_name,'w')
    start_with_residue=1
    crystal_symmetry=None
    if crystal_symmetry:
      print >> out, pdb.format_cryst1_record(
         crystal_symmetry=crystal_symmetry)
    for model in hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          for residue in conformer.residues():
            resseq_int=residue.resseq_as_int()
            ok=False
            for id,residue_ranges in zip(group,residue_range_list):
                if id in self.exclude_chains: continue
                if chain.id==id:
                  for start,end in residue_ranges:
                    if resseq_int>=start and resseq_int<=end:
                      ok=True
                      break
                if ok: break
            if ok:
              for atom in residue.atoms():
                print >>f, atom.format_atom_record()

    f.close()
    return file_name

  def special_cases(self,args):
    # special cases for input files so user doesn't need to specify:
    new_args=[]
    if not args: return []
    for arg in args:
      if (os.path.isfile(arg)):
        if arg[-3:] in ['pdb','PDB']:
          arg_use='simple_ncs_from_pdb.pdb_in='+arg
        else:
          arg_use=arg
      else:
        arg_use=arg
      new_args.append(arg_use)
    return new_args

  def get_chain_length(self,chains):
    shortest_chain=9999999
    longest_chain=0
    for chain in chains:
      length=len(chain)
      if length<shortest_chain: shortest_chain=length
      if length>longest_chain: longest_chain=length
    return shortest_chain,longest_chain

  def max(self,x,y):
    if x>=y: return x
    return y

  def ncs_object(self):
    # This is an instance of "ncs" as in "$PHENIX/phenix/phenix/autosol/ncs.py"
    if hasattr(self,'ncs_object'):
      return self.ncs_object

  def is_pdb(self,file_name): # XXX REMOVE PDB 2013-08-20
    if not os.path.isfile(file_name): return False
    for line in open(file_name).readlines():
       if line[0:4]=='ATOM' or line[0:6]=='HETATM':
         return True
    return False

  def get_help(self,command_name,master_params,summary):
    print >>self.log,summary
    print >>self.log,"Default parameters:"
    master_params.format(python_object=
          master_params.fetch(sources=[]).extract()).show()

  def raise_missing(self,what):
      raise Sorry("""\
          Missing file name for %(what)s :
          Please add %(what)s=file_name
          to the command line to specify %(what)s .""" % vars())

  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Find ncs among chains in a PDB file "
    header+="\n\n# type phenix.doc for help\n"

    summary= "usage: %s protein.pdb [parameter=value ...]" % command_name
    summary+="\n\nYou can set any parameter by specifying its path: "
    summary+="\nncs.max_rmsd =3 sets rmsd to 3"
    summary+="\nTo test use: %s exercise\n" % command_name
    return summary,header

  def exercise(self):  # run with a few ha sites
    text="""
ATOM      2  CA  GLY A   1      43.603 -11.488  24.325  1.00 35.57      2MLT 113
ATOM      6  CA  ILE A   2      44.200  -8.183  22.475  1.00 27.55      2MLT 117
ATOM     14  CA  GLY A   3      43.999 -10.264  19.329  1.00 21.05      2MLT 125
ATOM     18  CA  ALA A   4      40.378 -11.260  20.106  1.00 21.80      2MLT 129
ATOM     23  CA  VAL A   5      39.355  -7.658  21.083  1.00 19.34      2MLT 134
ATOM     30  CA  LEU A   6      41.062  -6.432  17.957  1.00 17.59      2MLT 141
ATOM     38  CA  LYS A   7      39.079  -8.646  15.636  1.00 22.55      2MLT 149
ATOM     47  CA  VAL A   8      35.792  -7.369  17.211  1.00 20.52      2MLT 158
ATOM     54  CA  LEU A   9      36.899  -3.759  16.725  1.00 16.83      2MLT 165
ATOM     62  CA  THR A  10      37.338  -4.508  13.084  1.00 19.41      2MLT 173
ATOM     69  CA  THR A  11      34.132  -6.405  12.343  1.00 24.14      2MLT 180
ATOM     76  CA  GLY A  12      31.584  -6.595  15.140  1.00 24.17      2MLT 187
ATOM     80  CA  LEU A  13      31.923  -2.919  16.364  1.00 23.24      2MLT 191
ATOM     88  CA  PRO A  14      31.026  -1.278  13.030  1.00 17.52      2MLT 199
ATOM     95  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10      2MLT 206
ATOM    100  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20      2MLT 211
ATOM    108  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41      2MLT 219
ATOM    116  CA  SER A  18      24.825   0.827  12.744  1.00 19.98      2MLT 227
ATOM    122  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72      2MLT 233
ATOM    136  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67      2MLT 247
ATOM    144  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84      2MLT 255
ATOM    153  CA  ARG A  22      18.740   2.273  14.020  1.00 20.38      2MLT 264
ATOM    164  CA  LYS A  23      17.500   1.928  17.550  1.00 22.62      2MLT 275
ATOM    173  CA  ARG A  24      18.059   5.674  18.276  1.00 27.11      2MLT 284
ATOM    184  CA  GLN A  25      15.836   6.730  15.339  1.00 37.50      2MLT 295
ATOM    193  CA  GLN A  26      13.132   4.360  16.583  1.00 46.66      2MLT 304
ATOM    204  CA  GLY B 101      26.196  11.215  25.772  1.00 31.28      2MLT 315
ATOM    208  CA  ILE B 102      25.695   8.457  23.127  1.00 26.61      2MLT 319
ATOM    216  CA  GLY B 103      25.775  10.608  19.984  1.00 20.83      2MLT 327
ATOM    220  CA  ALA B 104      29.261  12.192  20.662  1.00 22.19      2MLT 331
ATOM    225  CA  VAL B 105      30.560   8.605  21.561  1.00 20.43      2MLT 336
ATOM    232  CA  LEU B 106      29.339   7.258  18.306  1.00 15.03      2MLT 343
ATOM    240  CA  LYS B 107      30.789  10.122  16.413  1.00 20.35      2MLT 351
ATOM    249  CA  VAL B 108      34.238   9.637  18.027  1.00 23.38      2MLT 360
ATOM    256  CA  LEU B 109      34.150   5.918  17.183  1.00 20.22      2MLT 367
ATOM    264  CA  THR B 110      33.081   6.344  13.562  1.00 23.85      2MLT 375
ATOM    271  CA  THR B 111      35.833   8.796  12.896  1.00 26.41      2MLT 382
ATOM    278  CA  GLY B 112      38.482   7.602  15.318  1.00 22.53      2MLT 389
ATOM    282  CA  LEU B 113      38.268   3.810  15.278  1.00 23.61      2MLT 393
ATOM    290  CA  PRO B 114      39.945   3.149  11.872  1.00 20.15      2MLT 401
ATOM    297  CA  ALA B 115      43.145   4.994  13.009  1.00 22.26      2MLT 408
ATOM    302  CA  LEU B 116      43.070   3.167  16.352  1.00 18.86      2MLT 413
ATOM    310  CA  ILE B 117      42.952  -0.212  14.601  1.00 14.88      2MLT 421
ATOM    318  CA  SER B 118      45.891   0.689  12.345  1.00 17.71      2MLT 429
ATOM    324  CA  TRP B 119      47.907   1.909  15.276  1.00 17.53      2MLT 435
ATOM    338  CA  ILE B 120      47.296  -1.149  17.359  1.00 14.07      2MLT 449
ATOM    346  CA  LYS B 121      48.314  -3.437  14.512  1.00 20.04      2MLT 457
ATOM    355  CA  ARG B 122      51.517  -1.513  14.068  1.00 25.56      2MLT 466
ATOM    366  CA  LYS B 123      52.367  -1.651  17.742  1.00 23.85      2MLT 477
ATOM    375  CA  ARG B 124      51.670  -5.431  17.781  1.00 28.09      2MLT 486
ATOM    386  CA  GLN B 125      54.153  -6.012  14.932  1.00 40.40      2MLT 497
ATOM    395  CA  GLN B 126      56.818  -4.124  16.883  1.00 45.23      2MLT 506
        """

    expected_result="""refinement.ncs.restraint_group {
reference = chain 'A' and (resseq 1:26 )
selection = chain 'B' and (resseq 101:126 )
}"""

    file_name='temp.pdb'
    f=open(file_name,'w')
    f.write(text)
    f.close()
    args=[file_name]
    args.append("min_length=1")
    args.append("min_percent=10")
    run_simple_ncs_from_pdb=simple_ncs_from_pdb(args=args,log=self.log)

    result=open('simple_ncs_from_pdb.ncs').read()
    if result and \
       string.replace(string.replace(result," ",""),"\n","")== \
       string.replace(string.replace(expected_result," ",""),"\n",""):
      print >>self.log,'OK'
    elif result:
      print >>self.log,"Output does not match. Result: ",result,\
        "\nExpected result: ",expected_result
    else:
      print >>self.log,"No result"

  def process_inputs(self,args_read):
   if not args_read: args_read=[]
   if self.log is None :
    self.log=sys.stdout
   self.args=[]
   sugg=None
   self.suggested_ncs_groups=None
   self.ignore_chains=[]
   self.allow_recursion=True
   self.exact_match_only=False

   for arg in args_read:
    if arg[:5]=='sugg=':
      sugg=arg[5:]
      print >>self.log,"Suggested groups:",sugg
    elif arg[:7]=='ignore=':
      self.ignore_chains=[]
      for char in arg[7:]:
        self.ignore_chains+=char
      print >>self.log,"Ignoring chains:",self.ignore_chains
    elif arg=='exact_match':
      self.exact_match_only=True
    elif arg=='no_recursion':
      self.allow_recursion=False
    else:
      self.args.append(arg)
   self.suggested_ncs_groups=sugg

# FIXME: can't just pass a phil file path for some reason
class launcher (object) :
  def __init__ (self, params, tmp_dir) :
    adopt_init_args(self, locals())

  def __call__ (self) :
    os.chdir(self.tmp_dir)
    ncs = simple_ncs_from_pdb(
      args=None,
      params=self.params,
      quiet=True,
      exclude_h=True,
      exclude_d=True,
      log=sys.stdout)
    return ncs.ncs_object

if (__name__ == "__main__"):
  argument_list=sys.argv[1:]
  if True or is_debug(argument_list).value:
    sys_stdout_sav=sys.stdout
    simple_ncs=simple_ncs_from_pdb(args=sys.argv[1:])
    sys.exit(0)
  try:
    sys_stdout_sav=sys.stdout
    simple_ncs=simple_ncs_from_pdb(args=sys.argv[1:])
  except KeyboardInterrupt:
    pass
  except Exception, e:
    print "\n************************************************"
    print e
    print "\n************************************************"
    if sys.stdout != sys_stdout_sav:
      sys.stdout=sys_stdout_sav # restore output stream so we can print error
      print "\n************************************************"
      print e
      print "\n************************************************"
    from phenix.utilities.is_raise_sorry import is_raise_sorry
    if is_raise_sorry(argument_list).value:
      from libtbx.utils import Sorry
      raise Sorry(e)

