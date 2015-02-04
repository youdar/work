# fixme: fix comments
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
length >= xxxxx and percent identity >=min_percent.  A pair of segments
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
from phenix.utilities.composite_params import get_composite_master_params
from phenix.utilities.arg_display_methods import arg_display_methods
from phenix.utilities.list_methods import list_methods
from phenix.utilities.is_debug import is_debug
from phenix.utilities.headers import headers
from phenix.utilities import citations
from cctbx.array_family import flex
from libtbx import adopt_init_args
from libtbx.utils import Sorry
from libtbx.utils import Usage
import sys, os, string
from iotbx import pdb
import iotbx.phil
import iotbx.ncs
import iotbx

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
 min_percent = 80.
   .help = '''Threshold for similarity between chains, where similarity
   define as: (number of matching res) / (number of res in longer chain)'''
   .type = float
   .short_caption = Min. percent identity
   .expert_level = 0
   .style = bold noauto
 max_rmsd = 2.
   .help = '''limit of rms difference between chains to be considered
      as copies'''
   .type = float
    .short_caption = Max. RMSD
   .expert_level = 0
   .style = bold noauto
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
   match_radius = 4.0
   .type = float
   .help = '''max allow distance difference between pairs of matching
        atoms of two residues'''
   similarity_threshold = 0.95
   .type=float
   .help='''Threshold for similarity between matching chains.
      A smaller value cause more chains to be grouped together and can lower
      the number of common residues'''
   min_contig_length = 3
   .help = "segments < min_contig_length rejected"
   .type=int
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

usage_string = '''
phenix.simple_ncs_from_pdb file.pdb [options ...]

Options:
  pdb_in (str):         Input PDB file to be used to identify ncs
  temp_dir (str):       temp. directory (ncs_domain_pdb will be written there)
  min_percent (80.):    Similarity threshold between chains, where similarity
                        define as:(# of matching res)/(# of res in longer chain)
  max_rmsd (2.0):       limit of rms difference between chains to be considered
  match_radius (4.0):   max allow distance difference between pairs of matching
                        atoms of two residues
  similarity_threshold (0.95):
                        Threshold for similarity between matching chains
                        A smaller value cause more chains to be grouped together
                        and can lower the number of common residues
  min_contig_length (3):
                        reject if segments length < min_contig_length
  ncs_domain_pdb_stem (""):
                        NCS domains written to ncs_domain_pdb_stem+"group_"+nn
  write_ncs_domain_pdb (False):
                        You can write out PDB files representing NCS domains
                        for density modification if you want
  verbose (False):      Verbose output
  raise_sorry (False):  Raise sorry if problems
  debug (False):        Debugging output
  dry_run (False):      Just read in and check parameter names

Example:
  phenix.clashscore xxxx.pdb min_percent=85 min_contig_length=10
'''

class simple_ncs_from_pdb(arg_display_methods,list_methods,headers):

  def __init__(self,
               args   = None,
               params = None,
               ignore_chains   = None,
               exclude_chains = None,
               ncs_master_params = ncs_master_params,
               command_name   = "simple_ncs_from_pdb",
               pdb_inp = None,
               hierarchy = None,
               suppress_print = False,
               source_info    = None,
               pdb_file       = None,
               log = None,
               quiet=False,
               exclude_h=False,
               exclude_d=False,
               write_ncs_domain_pdb=None,
               ncs_domain_pdb_stem=None,
               temp_dir=None
               ):
    """
    Find NCS relation

    Args:
      ignore_chains (list): list of chain IDs to ignore (Write in log)
      exclude_chains (list): list of chain IDs to ignore (do not write in log)
      pdb_inp (object): PDB input object
      hierarchy (object)
      pdb_file (str): pdb file name

    """
    if not log: log = sys.stdout
    if not ignore_chains: self.ignore_chains = []
    if not exclude_chains: self.exclude_chains = []
    if ignore_chains:
      self.exclude_chains.extend(ignore_chains)
    self.log=log
    self.quiet=quiet
    self.exclude_h=exclude_h
    self.exclude_d=exclude_d
    self.process_similar_chains=True
    args=catenate_equals(args).new_args()
    self.process_inputs(args)
    if self.ignore_chains:
      print >>self.log,"Ignoring chains:",self.ignore_chains
    args=self.args
    self.Name='simple_ncs_from_pdb'
    if not suppress_print:
      citations.add_citation('phenix','simple_ncs_from_pdb')
    # run a test case
    if args and 'exercise' in args:
      self.exercise()
      return
    # process parameters
    args=self.special_cases(args)
    master_params=get_composite_master_params(
         method_list=['simple_ncs_from_pdb'],
         location_list=['phenix.command_line'])
    args=self.get_keyword_table(args,out=self.log)   # set self.keyword_table
    summary,header=self.get_summary_and_header(command_name)
    done,master_params,new_params,changed_params,help=self.get_params(
        command_name,master_params,args,out=self.log)
    if params is None :
      params = new_params
    if done: return
    if not quiet:
      print >>self.log, header
    if help or (params and params.simple_ncs_from_pdb.verbose):
      print >>self.log, "Values of all params:"
      master_params.format(python_object=params).show(out=log)
    if help or (params is None): return
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
    #
    if params.simple_ncs_from_pdb.dry_run:
      print >> self.log,"ARGS: ",args
      return
    self.raise_sorry = params.simple_ncs_from_pdb.raise_sorry
    # read in the PDB file if needed
    pdb_str = ''
    if(pdb_inp is None) and (hierarchy is None):
      if(pdb_file is None):
        if args and args[0] and os.path.isfile(args[0]):
          pdb_file=args[0]
        elif params.simple_ncs_from_pdb.pdb_in is not None:
          pdb_file=params.simple_ncs_from_pdb.pdb_in
        else:
          if self.raise_sorry:
            msg = "\nNeed PDB file for simple_ncs_from_pdb\n\n"
            msg += "Please make the PDB file the first argument like this:\n"
            msg += "phenix.simple_ncs_from_pdb mypdb.pdb ...\n"
            raise Sorry(msg)
      if not os.path.isfile(pdb_file):
         raise Sorry("The file "+str(pdb_file)+" is missing?")
      pdb_str = open(pdb_file).read()
      raw_records = flex.std_string()
      # Note the iotbx.pdb.input removes the MTRIX and BIOMT info. So both
      # pdb_inp and the hierarchy created from it, will not contain this info
      raw_records.extend(flex.split_lines(pdb_str))
      pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records)
      self.source_info = pdb_file
    # set input params
    self.verbose=params.simple_ncs_from_pdb.verbose
    if self.verbose:
      if not hierarchy:
        hierarchy = pdb_inp.construct_hierarchy()
        if not pdb_str:
          # create only if not existing, so not to remove MTRIX info
          pdb_str = pdb_inp.as_pdb_string()
      hierarchy.show(out=self.log)
    if source_info:
      self.source_info = source_info
    if (not hasattr(self,'source_info')) or (not self.source_info):
      self.source_info = "None"
    # Use cctbx tool to search for ncs and get ncs object
    find_param = params.simple_ncs_from_pdb.domain_finding_parameters
    inp_param = params.simple_ncs_from_pdb
    # choose input method (the order is important.)
    if pdb_str:
      input_method = {'pdb_string':pdb_str}
    elif hierarchy:
      input_method = {'hierarchy':hierarchy}
    else:
      input_method = {'pdb_inp':pdb_inp}
    ncs_obj = iotbx.ncs.input(
      chain_similarity_limit=find_param.similarity_threshold,
      min_contig_length=find_param.min_contig_length,
      min_percent=inp_param.min_percent,
      max_rmsd=inp_param.max_rmsd,
      max_dist_diff=find_param.match_radius,
      exclude_chains=self.exclude_chains,
      use_minimal_master_ncs=True,
      process_similar_chains=True,
      allow_different_size_res=True,
      exclude_misaligned_residues=True,
      check_atom_order=False,
      write_messages=False,
      log=self.log,
      quiet=self.quiet,
      **input_method)
    spec_object = ncs_obj.get_ncs_info_as_spec(
      write=False,exclude_h=self.exclude_h,exclude_d=self.exclude_d)
    ncs_object = spec_object
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

  def max(self,x,y):
    if x>=y: return x
    return y

  def ncs_object(self):
    # This is an instance of "ncs" as in "$PHENIX/phenix/phenix/autosol/ncs.py"
    if hasattr(self,'ncs_object'):
      return self.ncs_object

  def get_chain_length(self,chains):
    """ return the shortest and longest length of chains, from a list of
    chains """
    shortest_chain=9999999
    longest_chain=0
    for chain in chains:
      length=len(chain)
      if length<shortest_chain: shortest_chain=length
      if length>longest_chain: longest_chain=length
    return shortest_chain,longest_chain

  def is_pdb(self,file_name):
    """ check if file_name is a name of a pdb file that contains atoms """
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
    """ Process exact_match and ignore parameters """
    if not args_read: args_read=[]
    self.args=[]
    for arg in args_read:
      if arg=='exact_match':
        self.process_similar_chains = False
      elif 'ignore=' in arg:
        i = arg.index()
        for char in arg[i+7:]:
          self.ignore_chains += char
          self.exclude_chains += char
      else:
        self.args.append(arg)

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
  arg_list=sys.argv[1:]
  if (not arg_list) or ('--help' in arg_list) or ('-h' in arg_list):
    raise Usage(usage_string)
  if True or is_debug(arg_list).value:
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
    if is_raise_sorry(arg_list).value:
      from libtbx.utils import Sorry
      raise Sorry(e)
