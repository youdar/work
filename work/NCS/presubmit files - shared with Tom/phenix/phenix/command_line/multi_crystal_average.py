"""
multi_crystal_average.py
Applies multi-crystal NCS density modification
"""
from __future__ import division

from phenix.utilities.is_debug import is_debug
from phenix.utilities.composite_params import write_changed_params_to_file
from phenix.command_line.simple_ncs_from_pdb import simple_ncs_from_pdb
from phenix.command_line.get_cc_mtz_mtz import get_cc_mtz_mtz
import iotbx
import iotbx.phil
import libtbx.phil.command_line
from iotbx import pdb
from cctbx.array_family import flex
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import Sorry,null_out
from libtbx import runtime_utils
from phenix.autosol import iotbx_pdb_v0
from phenix.autosol.UserMethods import GeneralMethods
from phenix.autosol.copy_to_temp_dir import copy_to_temp_dir
from phenix.autosol.copy_from_temp_dir import copy_from_temp_dir
from phenix.autosol.trim_file_name import trim_file_name
from phenix.autosol.run_resolve import run_resolve
from phenix.utilities.catenate_equals import catenate_equals
from cStringIO import StringIO
import os,sys

master_params=iotbx.phil.parse("""
    multi{

      crystal
        .multiple = True
        .short_caption = Crystal data
        .caption = At least a PDB file and experimental data are required \
          for each crystal form.  You must also provide map coefficients \
          for at least one crystal, and if you have multiple NCS groups, \
          all crystals must have map coefficients.
        .style = auto_align box caption_width:720
      {

        pdb_file= None
          .type = path
          .short_caption = PDB file
          .help = '''
           PDB files, one for each crystal (One pdb_file=xxx.pdb per crystal)
           These should be in the same order as datafiles and map files.
           They are used to identify the NCS within each crystal and
           between crystals.  You should create these by placing the
           unique set of atoms (the NCS asymmetric unit) in each
           NCS asymmetric unit of each unit cell. Normally you would
           do this by carrying out molecular replacement on each
           crystal with the same search model.
           '''
          .style = bold

        map_coeffs = None
          .type = path
          .short_caption = Map coefficients
          .help = '''
              Mtz files with map coefficients.(One map_coeffs=xxx.mtz per crystal).
              At least one crystal must have map coefficients. Use
              &quot;None&quot; with quotes
              for any crystals that do not have starting maps.
             NOTE: If you have
              multiple NCS groups then you need map coefficients for all
              crystals.
              '''
          .style = bold OnChange:guess_map_coeff_labels

        map_coeffs_labin = None
          .type = str
          .short_caption = Map coefficient labels
          .help = '''
            Optional labin lines for mtz files with map coefficients. (One
          map_coeffs_labin=my_labin per crystal, or none at all)
           They look like map_coeffs_labin=" 'FP=FP PHIB=PHIM FOM=FOMM'"
           Put each set of labin values inside single quotes, and the whole
           list inside double quotes.
           You can leave out a labin statement for a file by putting in &quot;None&quot;
            and the routine will guess the column labels'''

        datafile = None
          .type = path
          .short_caption = Data file
          .help = ''' Mtz files with structure factors
            and optional phases and FOM and optional HL coefficients.
          One datafile for each crystal to be included (One datafile=xxx.mtz
           for each crystal).
          '''
          .style = bold OnChange:guess_data_labels

        datafile_labin = None
          .type = str
          .input_size = 400
          .short_caption = Data labels
          .help = '''
            Optional labin line for mtz file (In same order as mtz file).
            It can contain FP SIGFP [PHIB FOM] [HLA HLB HLC HLD].
           It looks like this:
           datafile_labin='FP=FP SIGFP=SIGFP PHIB=PHIM FOM=FOMM'
           You can leave out a labin statement for a file by putting in
            &quot;None&quot;
            and the routine will guess the column labels
            NOTE: If you supply HL coefficients they will be used in phase
              recombination.  If you supply PHIB or PHIB and FOM and not HL
              coefficients, then HL coefficients will be derived from your
              PHIB and FOM and used in phase recombination. '''

        solvent_content = None
          .type = float(value_min=0,value_max=1)
          .input_size = 64
          .help = '''
            Solvent content (0 to 1, typically 0.5) for each crystal
            (one solvent_content=xxx for each crystal).
            '''

        ha_file = None
          .type = path
          .short_caption = Heavy-atom sites
          .help = '''
           Optional file, normally containing heavy atom sites or other
           coordinates that are present in this structure, but not in the
           others. If supplied, the density near (within 2 A) of atoms in
           this file will not be transferred to other crystals.
           '''

        perfect_map_coeffs = None
          .type = path
          .help = ''' Optional mtz files with perfect map coefficients
            for comparison. '''
          .style = noauto

        perfect_map_coeffs_labin = None
          .type = str
          .help = '''
            Optional labin lines for mtz files with perfect map coefficients.
            For comparison and checking only. Not normally used.
                  '''
          .style = noauto

        update_map = True
          .type = bool
          .help = "Update the map for this crystal on each cycle."
                  "Set to False if this crystal represents a high-resolution"
                  "model or is otherwise expected to be already as good as it"
                  "can be."
          .style = noauto

      }

      averaging {
        cycles = 5
          .type = int
          .help = '''
            Number of cycles of density modification
             '''
          .short_caption = Cycles of density modification
          .input_size = 64
          .style = spinner bold

        resolution = None
          .type = float
          .help = "high-resolution limit for map calculation"
                  "Default is use all data"
          .short_caption = High resolution
          .input_size = 64

        fill = True
          .type = bool
          .help = "Fill in all missing reflections to resolution res_fill."
          .short_caption = Fill in missing reflections
          .style = OnChange:toggle_res_fill

        res_fill = None
          .type = float
          .help = "Resolution for filling in missing data (default = highest"
                  "resolution of any datafile)"
                "Default is use all data"
          .short_caption = Resolution for filling in missing data
          .input_size = 64

        temp_dir = "temp_dir"
          .type = path
          .help = "Optional temporary work directory"
          .style = hidden

        output_dir = ""
          .type = path
          .help = "Output directory where files are to be written"
          .style = hidden

        use_model_mask = False
          .type = bool
          .help = '''
            You can use the PDB files you input to define the solvent
           boundary if you wish. These will partially define the NCS
           asymmetric unit (by limiting it to the non-solvent region) but
           the exact NCS asymmetric unit will always be defined automatically
           (by the overlap of NCS-related density).
           Note that this is different than the command write_ncs_domain_pdb
           which defines individual regions where NCS applies for each
           domain.
          '''
          .short_caption = Define solvent boundary with model

        sharpen = False
          .type = bool
          .help = '''
            You can sharpen the maps or not in the density-modification process.
            (They are unsharpened at the end of the process if so).
            Not normally used, as an anisotropy correction with sharpening
            is normally applied to all the data.
          '''
          .short_caption = Sharpen maps in density modification

        equal_ncs_weight = False
          .type = bool
          .help = '''
            You can fix the NCS weighting to equally weight all copies.
          '''
          .short_caption = Equal weighting for all NCS copies
          .style = OnChange:toggle_weight_ncs

        weight_ncs = None
          .type = float
          .help = '''
           You can set the weighting on NCS symmetry (and cross-crystal averaging)
          '''
          .short_caption = Weight on NCS symmetry
          .input_size = 64

        write_ncs_domain_pdb = False
          .type = bool
          .help = '''
            You can use the input PDB files to define NCS boundaries.
            The atoms in the PDB files will be grouped into domains
            during the analysis of NCS and written out to domain-specific
            PDB files. (If there is only one domain or NCS group then there will
            be only one domain-specific PDB file and it will be the same as
            the starting PDB file.) Then the domain-specific PDB files will
            be used to define the regions over which the corresponding NCS
            operators apply.
            Note that this is different than the command use_model_mask which
            only defines the overall solvent boundary with your model.
              '''
          .short_caption = Define NCS boundaries with model

        mask_cycles = 1
          .type = int
          .help = '''
            Number of mask cycles in each cycle of density modification
          '''
          .short_caption = Mask cycles
          .input_size = 64
          .style = spinner
      }

      aniso
        .short_caption = Anisotropy correction
        .style = auto_align columns:2
      {
        remove_aniso = True
          .type = bool
          .help = "Remove anisotropy from data files before use"
                  "Note: map files are assumed to be already corrected"
          .short_caption = Remove anisotropy from data

        b_iso = None
          .type = float
          .help = "Target overall B value for anisotropy correction. Ignored"
                  "if remove_aniso = False.  If None,"
                  "default is minimum of (max_b_iso, lowest B of datasets,  "
                  "target_b_ratio*resolution)"
          .short_caption = Target overall B value
          .input_size = 64

        max_b_iso = 40.
          .type = float
          .help = "Default maximum overall B value for anisotropy correction. "
                  "Ignored if remove_aniso = False. Ignored if b_iso is set. "
                  "If used, default is minimum of (max_b_iso, "
                  "lowest B of datasets, target_b_ratio*resolution)"
          .short_caption = Maximum overall B value
          .input_size = 64

        target_b_ratio = 10.
          .type = float
          .help = "Default ratio of target B value to resolution for "
                  " anisotropy correction. "
                  "Ignored if remove_aniso = False. Ignored if b_iso is set. "
                  "If used, default is minimum of (max_b_iso, "
                  "lowest B of datasets, target_b_ratio*resolution)"
          .short_caption = Ratio of target B value to resolution
          .input_size = 64

      }
      control {
        verbose = True
          .type = bool
          .help = '''verbose output'''
          .style = hidden

        debug = False
          .type = bool
          .help = '''debugging output'''
          .style = hidden

        raise_sorry = False
          .type = bool
          .help = "Raise sorry if problems"
          .short_caption = raise sorry4
          .style = hidden

        coarse_grid= False
          .type = bool
          .help = '''
            You can set coarse_grid in resolve
          '''

        resolve_size = 12
          .help = "   Size for solve/resolve"
                  "(&quot;&quot;,&quot;_giant&quot;,"
                  "&quot;_huge&quot;,&quot;_extra_huge&quot; or a number "
                  "where 12=giant 18=huge"
          .type = str
          .short_caption = RESOLVE binary size
          .expert_level = 2
          .style = hidden

        resolve_command_list = None
          .help = "   Commands for resolve. One per line in the form:   keyword"
               "value   value can be optional   Examples:   coarse_grid  "
               "resolution 200 2.0   hklin test.mtz   NOTE: for command-line"
               "usage you need to enclose the whole set of commands in double"
               "quotes (&quot;) and each individual command in single quotes (')"
               "like this: "
               "resolve_command_list=&quot;'no_build' 'b_overall 23' &quot; "
          .type = str
          .multiple = True
          .short_caption = RESOLVE commands
          .expert_level = 2


        dry_run = False
          .type = bool
          .help = '''Just read in and check parameter names'''
          .style = hidden
        base_gui_dir = None
          .type = path
          .help = GUI parameter only
          .expert_level = 3
          .short_caption = Output directory
          .style = bold output_dir
      }

 }
include scope phenix.command_line.simple_ncs_from_pdb.ncs_master_params

""", process_includes=True)

class multi_crystal_average(GeneralMethods):
  def __init__(self,args=[],quiet=False,out=sys.stdout):


    command_name='multi'
    from phenix.command_line.multi_crystal_average import master_params
    args=catenate_equals(args).new_args()
    args=self.get_keyword_table(args,out=out)       # set self.keyword_table
    summary,header=self.get_summary_and_header(command_name)
    done,master_params,params,changed_params,help=self.get_params(
        command_name,master_params,args,out=out)
    if done: return
    if not quiet: print >>out, header

    if help or (not quiet) or (params and params.multi.control.verbose):
      print >>out,"Values of all params:"
      master_params.format(python_object=params).show(out=out)

    if help or params is None: return

    debug=params.multi.control.debug

    file_name='multi_crystal_average_updated.eff'
    write_changed_params_to_file(file_name,master_params,changed_params)

    # Set up OutputDir (where files will go) and temp_dir (working directory)
    self.Facts={'OutputDir':''}
    if params.multi.averaging.output_dir:
      self.Facts['OutputDir']=params.multi.averaging.output_dir
      if not os.path.exists(params.multi.averaging.output_dir):
        os.mkdir(params.multi.averaging.output_dir)
    if params.multi.averaging.temp_dir:
      if not os.path.exists(params.multi.averaging.temp_dir):
        os.mkdir(params.multi.averaging.temp_dir)
      self.Facts['temp_dir']=params.multi.averaging.temp_dir
    else:
      self.Facts['temp_dir']=self.create_temp_dir()

    self.remove_stopwizard()

    if params.multi.control.resolve_command_list is None:
      self.Facts['resolve_command_list']=[]
    else:
      self.Facts['resolve_command_list']=params.multi.control.resolve_command_list
    if params.multi.control.coarse_grid:
      self.Facts['resolve_command_list'].append('coarse_grid')
      print >>out,"Using coarse grid in density modification"
    if params.multi.averaging.sharpen:
      self.Facts['resolve_command_list'].append('sharpen')
      print >>out,"Temporarily sharpening maps during density modification"
    if params.multi.averaging.equal_ncs_weight:
      self.Facts['resolve_command_list'].append('equal_ncs_weight')
      print >>out,"Using equal weighting for all NCS and inter-crystal averages"
    if params.multi.averaging.weight_ncs is not None:
      self.Facts['resolve_command_list'].append(
        'weight_ncs '+str(params.multi.averaging.weight_ncs))
      self.Facts['resolve_command_list'].append(
        'weight_ncs_init '+str(params.multi.averaging.weight_ncs))
      print >>out,"Setting ncs weighting to ",params.multi.averaging.weight_ncs

    # Set up a few paths for resolve
    self.Facts['resolve_size']=params.multi.control.resolve_size

    # set the resolution
    self.Facts['resolution']=params.multi.averaging.resolution

    # Make sure all the requested files exist
    for entry in ["map_coeffs",
      "datafile", "perfect_map_coeffs"]:
      for c in params.multi.crystal:
       if getattr(c,entry) is not None:
         self.check_file(getattr(c,entry))

    if params.multi.control.dry_run:
      print >>out,"ARGS: ",args
      return

    # Set up and get labin line for each map or data file or perfect map file

    if params.multi.aniso.remove_aniso and params.multi.aniso.b_iso is None:
      # pre-run datafile conversion to get hires and b_aniso_min
      text_io=StringIO()
      datafile,datafile_labin,hires_datafiles=\
       self.set_up_mtz_files(params,
         require_phib=False,list='datafile',
         program_map_labels= ['FP','SIGFP','PHIB','FOM',
                 'HLA','HLB','HLC','HLD','FreeR_flag'],suffix="_data_",
         remove_aniso=params.multi.aniso.remove_aniso, # correct aniso
         out=text_io,set_b_iso=True) # set b_iso value from datafiles
      # sets self.b_aniso_min
      if self.b_aniso_min is None:
        raise Sorry("Sorry, unable to read any datafiles?")

    # now do it for real
    datafile,datafile_labin,hires_datafiles=\
       self.set_up_mtz_files(params,
         require_phib=False,list='datafile',
         program_map_labels= ['FP','SIGFP','PHIB','FOM',
                 'HLA','HLB','HLC','HLD','FreeR_flag'],suffix="_data_",
         remove_aniso=params.multi.aniso.remove_aniso, # correct aniso
         out=out)
    if hires_datafiles:
      print >>out,"\nHigh resolution for datafiles: %7.2f A" %(hires_datafiles)
    else:
      raise Sorry( "Sorry, need at least 2 datafiles for multi...")

    map_coeffs,map_coeffs_labin,hires_map_coeffs=\
         self.set_up_mtz_files(params,
         require_phib=True,list='map_coeffs',
         program_map_labels= ['FP','PHIB','FOM'],suffix="_map_coeffs_",
         target_map_coeffs=True,out=out) # 2011-10-18 need this to get FWT
    if hires_map_coeffs:
      print >>out,"\nHigh resolution for map coeffs: %7.2f A" %(hires_map_coeffs)
    else:
      raise Sorry( "Sorry, need at least 2 maps for multi...")

    perfect_map_coeffs,perfect_map_coeffs_labin,hires_perfect=\
         self.set_up_mtz_files(params,
         require_phib=True,list='perfect_map_coeffs',
         program_map_labels= ['FP','PHIB','FOM'],suffix="_perfect_",
          target_map_coeffs=True,out=out)
    if hires_perfect:
      print >>out,"\nHigh resolution for perfect data: %7.2f A" %(hires_perfect)

    # Make sure we have some PDB files defining our NCS

    # Make sure we have solvent content and PDB file for each one

    pdb_files=[]
    ha_files=[]
    solvent_content_list=[]
    update_map_list=[]
    for c in params.multi.crystal:
      pdb_files.append(getattr(c,'pdb_file'))
      ha_files.append(getattr(c,'ha_file'))
      solvent_content_list.append(getattr(c,'solvent_content'))
      update_map_list.append(getattr(c,'update_map'))
      if not pdb_files[-1] or not solvent_content_list[-1]:
        raise Sorry("Sorry, please enter PDB files and solvent "+
        "content for each crystal "+\
       "with commands like: solvent_content=0.5 and pdb_file=mypdbfile.pdb")

    # Make sure PDB files exist and get a simple list of their names
    # and copy them to temp_dir
    print >>out,"\nList of PDB files for NCS",pdb_files

    pdb_file_list=self.get_pdb_files(
       pdb_files,self.Facts['temp_dir'])
    print >>out,"\nCopied list of PDB files for NCS",pdb_file_list
    ha_file_list=self.get_pdb_files(
       ha_files,self.Facts['temp_dir'],allow_missing=True)
    print >>out,"\nCopied list of heavy_atom mask files",ha_file_list
    out.flush()

    # Decide if we are going to use these PDB files for solvent masks
    self.use_model_mask=params.multi.averaging.use_model_mask
    out.flush()


    # set the fill resolution if necessary
    if params.multi.averaging.fill and params.multi.averaging.res_fill is None \
      and hires_datafiles is not None:
      params.multi.averaging.res_fill=hires_datafiles
      print >>out,"\nReflections will be filled to resolution of %7.2f A\n" %(
        params.multi.averaging.res_fill)
    self.res_fill=params.multi.averaging.res_fill
    self.fill=params.multi.averaging.fill

    # Get NCS operators for all N orderings of the crystals
    formatted_ncs_list,ncs_objects_list=self.get_ncs_operators_all_orders(
          pdb_file_list,ha_file_list,
          map_coeffs,map_coeffs_labin,
          datafile,datafile_labin,
          solvent_content_list,params,out=out)
    out.flush()


    # Now we are ready to do the work
    print >>out,"\n****** Ready to run multi-crystal NCS density modification****\n"
    out.flush()

    # Get starting CC to perfect map, if any
    if self.contains_something(perfect_map_coeffs):
      print >>out,"\nStarting CC to perfect maps: "
      self.get_cc_values(perfect_map_coeffs,perfect_map_coeffs_labin,
          map_coeffs,map_coeffs_labin)
    out.flush()

    for cycle in xrange(1,params.multi.averaging.cycles+1):
      print >>out,"\nCycle ",cycle
      for i in xrange(len(pdb_file_list)):
        crystal_number=i
        # decide if this particular crystal needs to have map updated
        if not update_map_list[i]:
          print >>out,"Skipping update of map for crystal ",crystal_number+1,datafile[crystal_number]
          continue
        print >>out,"Running density modification on crystal ",crystal_number+1,datafile[crystal_number]
        out.flush()
        # Reorder the datafiles, putting "crystal_number" as the first one:
        local_pdb_file_list,local_ha_file_list, \
          local_map_coeffs,local_map_coeffs_labin,\
          local_datafile,local_datafile_labin, \
          local_ncs_objects_list, \
          local_formatted_ncs_list,local_solvent_content_list =\
            self.get_files_in_order(crystal_number,pdb_file_list,ha_file_list,
              map_coeffs,map_coeffs_labin,
              datafile,datafile_labin,
              ncs_objects_list,
              formatted_ncs_list,solvent_content_list)
        # get the new map coeffs for this datafile
        id='cycle_'+str(cycle)+'_xl_'+str(crystal_number+1)
        new_map_coeff,new_map_coeff_labin= \
          self.run_denmod(crystal_number,id,params.multi.averaging.mask_cycles,
            local_pdb_file_list,local_ha_file_list,
            local_map_coeffs,local_map_coeffs_labin,
            local_datafile,local_datafile_labin,
            local_ncs_objects_list,
            local_formatted_ncs_list,local_solvent_content_list,
            write_ncs_domain_pdb=params.multi.averaging.write_ncs_domain_pdb,
            out=out)
        out.flush()

        # ...and replace the existing ones with the new map coeffs
        map_coeffs[crystal_number]=new_map_coeff
        map_coeffs_labin[crystal_number]=new_map_coeff_labin
        print >>out,"New map for crystal ",crystal_number+1,":",new_map_coeff

      # Get current CC to perfect map, if any
      if self.contains_something(perfect_map_coeffs):
        print >>out,"\nCurrent CC values to perfect maps:"
        self.get_cc_values(perfect_map_coeffs,perfect_map_coeffs_labin,
          map_coeffs,map_coeffs_labin,out=out)

      # Copy out current files
      for new_map_coeffs,dd in zip(map_coeffs,datafile):
        print >>out,"Current best map coeffs for",dd,"are in",new_map_coeffs
        copy_from_temp_dir(self.Facts['temp_dir'],new_map_coeffs,
              new_map_coeffs,OutputDir=self.Facts['OutputDir'])
      out.flush()
    self.map_coeffs = map_coeffs

  def replace_none(self,param_list): # replace all "None" in lists with None
    if type(param_list) != type([1,2,3]): return param_list
    new_list=[]
    for item in param_list:
      if type(item)==type("abc") and item.lower()=='none':
        item=None
      new_list.append(item)
    return new_list

  def contains_something(self,value_list):
    if not value_list: return False
    for x in value_list:
      if x: return True
    return False
  def get_cc_values(self,perfect_map_coeffs,perfect_map_coeffs_labin,
          map_coeffs,map_coeffs_labin,out=sys.stdout):
      for mtz_1,labin_1,mtz_2,labin_2 in zip(
          perfect_map_coeffs,perfect_map_coeffs_labin,
          map_coeffs,map_coeffs_labin):
        if mtz_1 is None or mtz_2 is None:
           continue  # skip anything that has no maps
        args=[]
        args.append("mtz_1="+os.path.join(self.Facts['temp_dir'],mtz_1))
        args.append("labin_1="+labin_1)
        args.append("mtz_2="+os.path.join(self.Facts['temp_dir'],mtz_2))
        args.append("labin_2="+labin_2)
        if self.Facts['resolution']:
          args.append("resolution="+str(self.Facts['resolution']))
        args.append("temp_dir="+self.Facts['temp_dir'])
        args.append("output_dir="+self.Facts['temp_dir'])
        ff=StringIO()
        get_cc=get_cc_mtz_mtz(args,quiet=True,out=ff)
        cc_value=get_cc.final_cc
        print >>out,"Map",mtz_2,"  CC=",self.round(cc_value,2)

  def run_denmod(self,crystal_number,id,mask_cycles,
            local_pdb_file_list,local_ha_file_list,
            local_map_coeffs,local_map_coeffs_labin,
            local_datafile,local_datafile_labin,
            local_ncs_objects_list,
            local_formatted_ncs_list,local_solvent_content_list,
            just_dump_map=False,
            write_ncs_domain_pdb=False,out=sys.stdout):

    # run one cycle of density modification on first file in list and return the
    # updated map coefficients and labin

    # if no starting map coeffs are available, get them first
    if local_map_coeffs[0] is None and not just_dump_map:
       if self.number_of_groups is not None and self.number_of_groups>1 \
            and not write_ncs_domain_pdb:
         raise Sorry("Sorry, if you have multiple NCS groups "+
           "and do not specify\nwrite_ncs_domain_pdb then "+
           "you need to supply starting maps for all crystals")
       local_mask_cycles=1
       local_id=id+'_dump'
       self.run_denmod(crystal_number,local_id,local_mask_cycles,
            local_pdb_file_list,local_ha_file_list,
            local_map_coeffs,local_map_coeffs_labin,
            local_datafile,local_datafile_labin,
            local_ncs_objects_list,
            local_formatted_ncs_list,local_solvent_content_list,
            just_dump_map=True,out=out)
       # now we have a map called ncs_dump.mtz that we can use with
       # pattern_phase to generate phases for this datafile
       hklout=local_id+"_pattern_phase.mtz"
       logfile=local_id+"_pattern_phase.log"
       command_2="cc_map_file dump_ncs.map"
       if self.fill: # fill to full resolution
         command_2+="\nfill \nres_fill %6.2f \nkeep_missing"%(self.res_fill)
       my_run_resolve=run_resolve(
             Facts=self.Facts,temp_dir=self.Facts['temp_dir'],
             resolution=[1000.,self.Facts['resolution']],
             size=self.Facts['resolve_size'],
             hklin=local_datafile[0],
             labin=local_datafile_labin[0],
             solvent_content=local_solvent_content_list[0],
             build="no_build",
             command_1="pattern_phase",
             command_2=command_2,
             ha_file="NONE",
             hklout=hklout,
             logfile=logfile)
       local_map_coeffs[0]=hklout
       local_map_coeffs_labin[0]="FP=FWT PHIB=PHWT"  # 2011-10-19


    # Now ready for density modification using all crystals

    logfile='denmod_'+id+'.log'
    hklout='denmod_'+id+'.mtz'
    command_2=""
    # command_2+="\ninclude_self"
    if self.fill: # fill to full resolution
      command_2+="\nfill \nres_fill %6.2f \nkeep_missing"%(self.res_fill)
    found_maps=False
    total_maps=1
    for file,labin,dd,datafile_labin in zip(
      local_map_coeffs[1:],local_map_coeffs_labin[1:],
      local_datafile[1:],local_datafile_labin[1:]):
      if file is not None:

        command_2+="\ncrystal_map "+file
        command_2+="\nlabin_asis  labin %s" %(labin)
        total_maps+=1
        if not found_maps:
          found_maps=True
      else:  # generate dummy map file with zero FOM
        dummy_map='dummy_map_'+id+'.mtz'
        dummy_log='dummy_map_'+id+'.log'
        labels_list=self.get_labels_from_labin(datafile_labin,
           target_label_list=['FP'],required=['FP'],out=out)
        short_labin="FP="+labels_list.split()[0].rstrip()
        my_run_resolve=run_resolve(
           Facts=self.Facts,temp_dir=self.Facts['temp_dir'],  # copy over
           resolution=[1000.,self.Facts['resolution']],
           size=self.Facts['resolve_size'],
           hklin=dd,
           labin=short_labin,
           hklout=dummy_map,
           mask_cycles=1, minor_cycles=0,build="no_build",
           solvent_content=0.5,
           logfile=dummy_log,
           command_1="keep_missing")

        command_2+="\ncrystal_map "+dummy_map
        command_2+="\nlabin_asis  labin FP=FWT PHIB=PHWT"  # 2011-10-19
        total_maps+=1

    command_1=local_formatted_ncs_list[0]
    if not found_maps:
      print >>out,"NOTE: no information available yet for other crystals."
      command_1=self.remove_ncs_domain_pdb(command_1)
    #command_2+="\nno_free"  # always...try this HERE JUST
    if just_dump_map:
      command_2+="\nno_free \nno_refine_ncs \nforce_ncs "
      command_2+="\nforce_ncs_boundary"
      command_2+="\ndump_ncs \nminor_cycles 1"

    elif found_maps and write_ncs_domain_pdb:
      command_2+="\nforce_ncs_boundary"
      command_2+="\noverlap_min 0.001"
      command_2+="\nfraction_ncs_min 0.001"

    if just_dump_map or self.use_model_mask:
      command_2+="\nmodel "+\
       str(trim_file_name(local_pdb_file_list[0]).trimmed_file)
      command_2+="\nuse_model_mask"

    have_ha=False
    if local_ha_file_list:
      for x in local_ha_file_list:
        if x: have_ha=True
    if have_ha:  # remove heavy atom sites
      file_name=self.create_mask_file(local_ncs_objects_list[0],local_ha_file_list,
         crystal_number=crystal_number)
      command_2+="\nremove_mask %s" %(file_name)

    try:
      my_run_resolve=run_resolve(Facts=self.Facts,
       n_crystal_map_max=total_maps,
       temp_dir=self.Facts['temp_dir'],
       size=self.Facts['resolve_size'],
       hklin=local_datafile[0],
       labin=local_datafile_labin[0],
       hklstart=local_map_coeffs[0],
       labstart=local_map_coeffs_labin[0],
       hklout=hklout,
       mask_cycles=mask_cycles,
       solvent_content=local_solvent_content_list[0],
       resolution=[1000.,self.Facts['resolution']],
       logfile=logfile,
       command_1=command_1,
       command_2=command_2)
    except KeyboardInterrupt: raise
    except Exception,e:
      raise Sorry(str(e)+"\nSorry, failed...error messages above")
    copy_from_temp_dir(self.Facts['temp_dir'],logfile,
              logfile,OutputDir=self.Facts['OutputDir'])

    logfile_output=os.path.join(self.Facts['OutputDir'],logfile)
    if not os.path.isfile(logfile_output):
      raise Sorry("Sorry, unable to run density modification?")

    map_coeff=hklout
    map_coeff_labin="FP=FWT PHIB=PHWT"  # 2011-10-19 map coeffs are FWT PHWT

    return map_coeff,map_coeff_labin

  def create_mask_file(self,ncs_objects,local_ha_file_list,skip_first=False,
    crystal_number=None):
    # add ha file as mask
    all_text=""
    # local_ncs_objects_list has N ncs objects, one for each pair of files,
    # file_1=target file, file_2=other file
    # ncs_objects=local_ncs_objects_list[0] (now passed in directly)
    if skip_first:
      ncs_objects_use=ncs_objects[1:]  # listing ha for files all but target
      local_ha_file_list_use=local_ha_file_list[1:]
    else:
      ncs_objects_use=ncs_objects
      local_ha_file_list_use=local_ha_file_list
    for ncs,ha_file in zip(ncs_objects_use,local_ha_file_list_use):
     if not ha_file or not os.path.isfile(ha_file): continue
     from phenix.command_line.apply_ncs import apply_ncs
     # RUN SEPARATELY FOR EACH CHAIN FROM LOCAL_HA_FILE (ASSUME ATOMS IN
     # ANY ONE CHAIN GO IN SAME NCS AU)
     pdb_input_all=iotbx.pdb.input(file_name=ha_file)
     hierarchy=pdb_input_all.construct_hierarchy()
     chain_list=self.get_chain_list(hierarchy)
     ff=StringIO()
     for c in chain_list:
       pdb_input=self.get_pdb_input_for_chain(hierarchy,c)
       args=['pdb_out=None',
           'max_copies=1',
           'start_copy=1',
           'no_match_copy=1',
        ]
       an=apply_ncs(args,pdb_input=pdb_input,out=ff,ncs_object=ncs)
       output_text=an.output_text
       all_text+=output_text
    lines=all_text.splitlines()
    have_cryst=False
    if crystal_number is not None:
      # this file is going to be for resolve to read in
      file_name='TEMP_HA_'+str(crystal_number+1)+'.pdb'
    else:
      file_name='TEMP_HA.pdb'
    full_file=os.path.join(self.Facts['temp_dir'],file_name)
    f=open(full_file,'w')
    if not os.path.isfile(full_file): # don't print out if it already is there
      print "\nHeavy-atom sites in frame of reference of crystal "+\
      "%s (to \nbe excluded) are in %s" %(crystal_number+1,f.name)
    for line in lines:
      if not line: continue
      if line.startswith("CRYST1"):
        if not have_cryst:
          have_cryst=True
        else:
          continue # already wrote it out
      print >>f,line.rstrip()
    f.close()
    return file_name

  def get_pdb_input_for_chain(self,pdb_hierarchy,c):
    f=StringIO()
    # get pdb_input for just chain c of this  hierarchy
    for model in pdb_hierarchy.models()[:1]:
      for chain in model.chains():
        if chain.id != c:continue
        for conformer in chain.conformers()[:1]:
          for residue in conformer.residues():
            atoms=residue.atoms()
            for atom in atoms:
              print >>f, atom.format_atom_record()
    pdb_text=f.getvalue()
    from cctbx.array_family import flex
    pdb_inp=iotbx.pdb.input(source_info='string',
                  lines=flex.split_lines(pdb_text))
    return pdb_inp

  def get_chain_list(self,hierarchy):
    chain_list=[]
    for model in hierarchy.models():
      for chain in model.chains():
        if not chain.id in chain_list:chain_list.append(chain.id)
    return chain_list

  def remove_ncs_domain_pdb(self,command):
    # just remove line with ncs_domain_pdb
    text=""
    for line in command.split("\n"):
      if not line.lower().find('ncs_domain_pdb')>-1:
         text+=line+"\n"
    return text

  def get_ncs_operators_all_orders(self,
          pdb_file_list,ha_file_list,
          map_coeffs,map_coeffs_labin,
          datafile,datafile_labin,
          solvent_content_list,params,out=sys.stdout):
    # Sequentially define each PDB file as the std one and get NCS operators
    #  relative to that one for all the files
    formatted_ncs_list=[]
    ncs_objects_list=[]
    for crystal_number in xrange(len(pdb_file_list)):
      print >>out,"\n--------------------------------------------------------"
      dummy_formatted_ncs_list=len(pdb_file_list)*[None]
      dummy_ncs_objects_list=len(pdb_file_list)*[None]
      local_pdb_file_list,local_ha_file_list,\
      local_map_coeffs,local_map_coeffs_labin,\
      local_datafile,local_datafile_labin, \
      local_ncs_objects_list, \
      local_formatted_ncs_list,local_solvent_content_list =\
        self.get_files_in_order(crystal_number,pdb_file_list,ha_file_list,
          map_coeffs,map_coeffs_labin,
          datafile,datafile_labin,
          dummy_ncs_objects_list,
          dummy_formatted_ncs_list,solvent_content_list)

      formatted_ncs=self.get_ncs_for_all_crystals(
         local_pdb_file_list,params,out=out)
      #print >>out,"\n\nNCS formatted for resolve using crystal %d (%s) as std : \n" %(crystal_number+1,pdb_file_list[crystal_number])
      #print >>out,formatted_ncs
      ncs_objects_list.append(self.ncs_object_list) # each one has N ncs objects, one for each (first_file, other_file) in list of PDB files in this order.
      formatted_ncs_list.append(formatted_ncs)
    return formatted_ncs_list,ncs_objects_list

  def set_up_mtz_files(self,params,
         require_phib=True,list='map_coeffs',
         program_map_labels= ['FP','PHIB','FOM'],suffix="_",
         target_map_coeffs=None,remove_aniso=None,out=sys.stdout,
         set_b_iso=False):
    # get labels, copy file to standard format in temp_dir, return
    #   name of new file and labels to use.
    # if remove_aniso=True then apply aniso correction
    print >>out,"\n"+list+" files:"
    hires=None
    self.b_aniso_min=None
    map_coeffs=[]
    map_coeffs_labin=[]
    list_labin=list+"_labin"
    count=0
    r=None
    for c in params.multi.crystal:
      coeffs=getattr(c,list)
      labels=getattr(c,list_labin)
      count+=1

      if not coeffs or coeffs.lower()=='none':
        new_coeffs=None
        new_labin=None
      else:
        self.resolution_input_hkl=None
        labin_1=self.set_up_labin_generic(
          mtz_file=coeffs,
          labin_input=labels,
          require_phib=require_phib,name=list,
           program_map_labels=program_map_labels,
           target_map_coeffs=target_map_coeffs)
        if self.resolution_input_hkl and \
           (not hires or self.resolution_input_hkl<hires):
          hires=self.resolution_input_hkl

        # Now read the file and make a new mtz with standard labels...
        copied_coeffs=self.stem(trim_file_name(coeffs).trimmed_file)+\
           suffix+str(count)+".mtz"
        copy_to_temp_dir(self.Facts['temp_dir'],coeffs,copied_coeffs)
        new_coeffs=self.stem(copied_coeffs)+"_ed.mtz"
        try:
          my_run_resolve=run_resolve(
           Facts=self.Facts,temp_dir=self.Facts['temp_dir'],  # copy over
           resolution=[1000.,self.Facts['resolution']],
           size=self.Facts['resolve_size'],
           hklin=copied_coeffs,
           labin=labin_1,
             hklout=new_coeffs,
             mask_cycles=1, minor_cycles=0,build="no_build",
             solvent_content=0.5,
             logfile="copy_"+self.stem(copied_coeffs)+".log",
             command_1="keep_missing")
        except Exception, e:
          raise Sorry(str(e)+\
          "\nPlease check column labels for the file "+str(coeffs)+":"+\
          "\nLabels chosen are: "+labin_1+\
          "\nAvailable labels are: "+str(self.all_input_map_labels))
          #note: self.all_input_map_labels is set with self.set_up_labin_generic
        new_labin=""
        for pair in ["FP=FP","SIGFP=SIGFP","PHIB=PHIM","FOM=FOMM",
          "HLA=HLAM","HLB=HLBM","HLC=HLCM","HLD=HLDM","FreeR_flag=FreeR_flag"]:
          target=pair.split("=")[0]+"="
          if labin_1.find(target)>-1:
            new_labin+=pair+" "
        print >>out,coeffs," LABIN: ",labin_1
        print >>out,"...copied to: ",os.path.join(self.Facts['temp_dir'],new_coeffs),\
         " LABIN: ",new_labin

        if remove_aniso: # NOTE all this is done in temp_dir
          temp_coeffs=new_coeffs
          temp_labin=new_labin
          new_coeffs=temp_coeffs[:-4]+"_aniso.mtz"

          new_labin_list=[]
          extension="_aniso"
          for pair in temp_labin.split():
            spl=pair.split("=")
            if pair in ["FP=FP","SIGFP=SIGFP"]:
              new_labin_list.append(pair+extension)
            else:
              new_labin_list.append(pair)
          new_labin=" ".join(new_labin_list)
          full_temp_coeffs=os.path.join(self.Facts['temp_dir'],temp_coeffs)
          full_new_coeffs=os.path.join(self.Facts['temp_dir'],new_coeffs)
          print >>out,"\nRemoving anisotropy from %s to yield %s" %(
            temp_coeffs,new_coeffs)
          # FP is always FP=FP here...so we don't need to specify it
          args=["mtz_in=%s" %(full_temp_coeffs),
            "mtz_out=%s" %(full_new_coeffs),
            "extension=%s" %(extension)]
          if params.multi.aniso.b_iso:
            args.append("b_iso=%s" %(str(params.multi.aniso.b_iso)))
          from phenix.command_line.remove_aniso import remove_aniso
          r=remove_aniso(args=args,out=out)
          if r.b_aniso_min is not None and \
              (self.b_aniso_min is None or r.b_aniso_min < self.b_aniso_min):
            self.b_aniso_min = r.b_aniso_min
          # new labin line has _aniso appended to FP, SIGFP and is in new_coeffs


      map_coeffs.append(new_coeffs)
      map_coeffs_labin.append(new_labin)

    if r and set_b_iso: #set params b_iso from what we know: Use routine in r
      params.multi.aniso.b_iso=r.get_b_iso(
          max_b_iso=params.multi.aniso.max_b_iso,
          b_iso=params.multi.aniso.b_iso,
          b_aniso_min=self.b_aniso_min,
          target_b_ratio=params.multi.aniso.target_b_ratio,
          hires=hires)
    return map_coeffs,map_coeffs_labin,hires

  def get_files_in_order(self,i,pdb_file_list,ha_file_list,
          map_coeffs,map_coeffs_labin,
          datafile,datafile_labin,
          ncs_objects_list,
          formatted_ncs_list,solvent_content_list):
          # make file i the first

    if i>=len(pdb_file_list): raise AssertionError,"Error in get_files_in_order"

    local_pdb_file_list=[ pdb_file_list[i] ]
    local_ha_file_list=[ ha_file_list[i] ]
    local_map_coeffs=[ map_coeffs[i] ]
    local_map_coeffs_labin=[ map_coeffs_labin[i] ]
    local_datafile=[ datafile[i] ]
    local_datafile_labin=[ datafile_labin[i] ]
    local_ncs_objects_list=[ ncs_objects_list[i] ]
    local_formatted_ncs_list=[ formatted_ncs_list[i] ]
    local_solvent_content_list=[ solvent_content_list[i] ]
    for j in xrange(len(pdb_file_list)):
      if i==j: continue
      local_pdb_file_list.append(pdb_file_list[j])
      local_ha_file_list.append(ha_file_list[j])
      local_map_coeffs.append(map_coeffs[j])
      local_map_coeffs_labin.append(map_coeffs_labin[j])
      local_datafile.append(datafile[j])
      local_datafile_labin.append(datafile_labin[j])
      local_ncs_objects_list.append(ncs_objects_list[j])
      local_formatted_ncs_list.append(formatted_ncs_list[j])
      local_solvent_content_list.append(solvent_content_list[j])
    return  local_pdb_file_list,local_ha_file_list,local_map_coeffs, \
     local_map_coeffs_labin,local_datafile,local_datafile_labin, \
     local_ncs_objects_list, \
     local_formatted_ncs_list,local_solvent_content_list



  def get_ncs_for_all_crystals(self,pdb_file_list,params,
      out=sys.stdout):

    outfile_list=[]
    # first let's see if we can extract just the unique part of the
    # first PDB file

    print >>out,'Finding NCS groups for crystal',pdb_file_list[0],'by itself'
    pdb_inp,hierarchy,all_chain_proxies=self.get_info('',pdb_file_list[0])
    ncs=simple_ncs_from_pdb(pdb_inp=pdb_inp,hierarchy=hierarchy,
       params=params,quiet=True,suppress_print=True)
    number_of_ncs_groups=len(ncs.ncs_object.ncs_groups())
    print >>out,"NCS groups in ",pdb_file_list[0],":",number_of_ncs_groups
    if number_of_ncs_groups==0:
       unique_chain_list=self.get_unique_chains(pdb_file_list[0])
    else:  # select unique chains from this file
      all_chain_list=self.get_unique_chains(pdb_file_list[0])
      # we want all the unique chains in NCS groups, plus anything not
      # in NCS groups that are in all_chain_list
      unique_chain_list=[]
      accounted_for_chain_list=[]
      for ncs_group in ncs.ncs_object.ncs_groups():
        chain_residue_id=ncs_group.chain_residue_id()
        print >>out,"Chains in NCS group: ",chain_residue_id
        [group,residue_id]=chain_residue_id
        unique_chain=group[0]
        accounted_for_chain_list+=group # accounted for by NCS group
        if not unique_chain in unique_chain_list:
          unique_chain_list.append(unique_chain)
      # now add in anything that is not accounted for
      for id in  all_chain_list:
        if not id in accounted_for_chain_list:
          unique_chain_list.append(id)
    # now read file and select just unique chains...
    unique_text,new_unique_chains=self.get_unique_model_with_unique_chains(
         unique_chain_list,pdb_file_list[0],pdb_file_list,out=out)

    f=open('unique_text.pdb','w')
    print >>f, unique_text
    f.close()

    print >>out,"List of unique chains for crystal %s:%s " %(
      pdb_file_list[0],unique_chain_list)
    print >>out,"List of renamed chains: ",new_unique_chains
    # ready with text of edited PDB file for unique chains in pdb_file_list[0]
    i=0
    ncs_group_list_dict={}
    used_file_names=[]
    print >>out,"\nFinding NCS to map each crystal "+\
        "to %s...\n" %(pdb_file_list[0])
    text_io=StringIO()

    ncs_object_list=[]
    for file in pdb_file_list:
      print >>out,"\nFinding NCS to map %s on to %s" %(file,pdb_file_list[0])
      if file in used_file_names:
         raise Sorry ("Sorry, please only enter each file once:"+\
          " duplicate is :"+str(file))
      used_file_names.append(file)
      i+=1
      pdb_inp,hierarchy,all_chain_proxies=self.get_info(unique_text,file)
      if i==1 and params.multi.averaging.write_ncs_domain_pdb:
        write_ncs_domain_pdb=True
        ncs_domain_pdb_stem=trim_file_name(file).trimmed_file+'_'
      else:
        write_ncs_domain_pdb=False # 2012-02-18 fixed was None
        ncs_domain_pdb_stem=None
      ncs=simple_ncs_from_pdb(pdb_inp=pdb_inp,hierarchy=hierarchy,
       params=params,quiet=True,suppress_print=True,
       exclude_chains=new_unique_chains,  # do not write them out
       write_ncs_domain_pdb=write_ncs_domain_pdb,
       ncs_domain_pdb_stem=ncs_domain_pdb_stem,
       temp_dir=self.Facts['temp_dir'])
      n=ncs.ncs_object
      ncs_object_list.append(n) # 2011-11-02
      n.display_all(log=out)
      # format it for resolve now
      out_1=StringIO()
      # crystal_number in ncs file for resolve starts at 0 2011-12-15
      n.format_all_for_resolve(out=out_1,crystal_number=i-1,
         skip_identity=True)


      text=out_1.getvalue()
      text=self.remove_blank_lines(text)
      spl=text.split("new_ncs_group")

      while spl[0]=='':
        if len(spl)<2:
          raise Sorry("Sorry no NCS groups found for the file "+file+
            "?")
        spl=spl[1:]
      if len(spl)>1:
        print >>out,'Model ',file,'has',len(spl),'NCS groups...'
      ncs_group_list_dict[file]=spl
    number_of_groups=None
    for file in pdb_file_list:
      spl=ncs_group_list_dict[file]
      if number_of_groups is None:
        number_of_groups=len(spl)
        print >>out,"Number of NCS groups is ",number_of_groups
      if len(spl)!=number_of_groups:
        line='Sorry, the number of NCS groups for '+str(file)+\
        '('+str(len(spl))+') does not \nmatch previous value of '+\
        str(number_of_groups)+"\nPlease make sure all monomers in"+\
       " the PDB files used for NCS definitions for multi are IDENTICAL"
        raise Sorry(line)
    for i in xrange(number_of_groups):
      print >>text_io, "new_ncs_group"
      for file in pdb_file_list:
        spl=ncs_group_list_dict[file]
        print >> text_io, spl[i]

    self.number_of_groups=number_of_groups
    self.ncs_object_list=ncs_object_list
    return text_io.getvalue()  # has all the NCS groups for this set of crystals
    # in this order

  def remove_blank_lines(self,text):
    spl=text.split("\n")
    new_spl=[]
    for s in spl:
      if s.replace(" ",""):
        new_spl.append(s)
    text="\n".join(new_spl)
    return text

  def get_unique_model_with_unique_chains(self,unique_chains,
        unique_model,pdb_file_list,out=sys.stdout):
    # rename chain on this one to something never used...
    all_chains=self.get_unique_chains(pdb_file_list)
    # decide on some new chain ID's that do not duplicate anything:
    try_list='-+0123456789'
    new_unique_chains=[]
    for i in unique_chains:
      found=False
      for a in try_list:
       if found: break
       for b in try_list:
         if found:break
         t=a+b
         if not t in new_unique_chains+unique_chains+all_chains:
           new_unique_chains.append(t)
           found=True
      if not found:
        raise Sorry("Sorry, unable to find unique chain ID not in "+\
           str(unique_chains)+" or "+str(all_chains))
    print >>out,"Unique chains from",unique_model,'renamed to:',\
       new_unique_chains ,"\n"

    # now replace the chain names in unique_chains with new_unique_chains
    # and make a new file

    unique_text=self.replace_chain_names(unique_model,
       unique_chains,new_unique_chains)

    return unique_text,new_unique_chains

  def replace_chain_names(self,file,unique_chains,new_unique_chains):

    out=StringIO()
    text=open(file).read()
    pdb_input = pdb.input(file_name=file)
    hierarchy = iotbx.pdb.input(file_name=file).construct_hierarchy()
    for model in hierarchy.models():
      chains = model.chains()
      for chain in chains:
        if not chain.id in unique_chains:
          continue  # skip chains we do not want!
        for o,n in zip(unique_chains,new_unique_chains):
          if chain.id==o: new_chain_id=n
        conformers = chain.conformers()
        serial=0
        for conformer in conformers:
          residues = conformer.residues()
          for residue in residues:
            atoms = list(residue.atoms())
            for atom in atoms:
              serial+=1
              print >> out, iotbx_pdb_v0.format_atom_record(
                serial=serial,
                name=atom.name,
                altLoc=conformer.altloc,
                resName=residue.resname,
                chainID=new_chain_id,
                resSeq=residue.resseq,
                iCode=residue.icode,
                site=atom.xyz,
                occupancy=atom.occ,
                tempFactor=atom.b,
                segID=atom.segid,
                element=atom.element,
                charge=atom.charge)

    unique_text=out.getvalue()
    return unique_text

  def get_unique_chains(self,file_list):
     unique_chains=[]
     if type(file_list)!=type([1,2,3]):
       file_list_use=[file_list]
     else:
       file_list_use=file_list
     for file in file_list_use:
       pdb_input = pdb.input(file_name=file)
       hierarchy = iotbx.pdb.input(file_name=file).construct_hierarchy()
       for model in hierarchy.models():
          chains = model.chains()
          for chain in chains:
            if not chain.id in unique_chains:
               unique_chains.append(chain.id)
     return unique_chains


  def get_info(self,unique_text,file):
      raw_records = flex.std_string()
      raw_records.extend(flex.split_lines(unique_text))
      raw_records.extend(flex.split_lines(open(file).read()))
      pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records)
      hierarchy=pdb_inp.construct_hierarchy()
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
      return pdb_inp,hierarchy,all_chain_proxies

  def check_file(self,file):
    if file.lower() != 'none' and not os.path.isfile(file):
      raise Sorry("Sorry, the file "+str(file)+" does not exist?")

  def get_pdb_files(self,file_list,temp_dir,allow_missing=False):
    new_file_list=[]
    count=0
    for file in file_list:
      count+=1
      if not file or not os.path.isfile(file):
        if allow_missing:
          new_file_list.append(None)
        else:
          raise Sorry("Sorry, the file "+str(file)+" does not exist?")
      else:
        new_file=self.stem(trim_file_name(file).trimmed_file)+"_"+str(count)+".pdb"
        copy_to_temp_dir(temp_dir,file,new_file)
        full_new_file=os.path.join(temp_dir,new_file)
        new_file_list.append(full_new_file)
    return new_file_list

  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Apply multi-crystal averaging with statistical density"
    header+="\n# modification "
    header+="\n\n# Type phenix.doc for help"

    summary= "usage: phenix.%s multi.eff " % command_name
    return summary,header

########################################################################
# GUI STUFF
def validate_params (params) :
  params = params.multi
  if (len(params.crystal) < 2) :
    raise Sorry("At least two crystal forms must be provided.")
  n_maps = 0
  for n, xtal in enumerate(params.crystal, start=1) :
    if (xtal.pdb_file is None) :
      raise Sorry("PDB file not defined for crystal #%d." % n)
    elif (not os.path.isfile(xtal.pdb_file)) :
      raise Sorry("The PDB file '%s' does not exist or is not a file." %
        xtal.pdb_file)
    if (xtal.datafile is None) :
      raise Sorry("No data file provided for crystal #%d." % n)
    elif (not os.path.isfile(xtal.datafile)) :
      raise Sorry("The data file '%s' does not exist or is not a file." %
        xtal.datafile)
    if (xtal.map_coeffs is not None) :
      if (not os.path.isfile(xtal.map_coeffs)) :
        raise Sorry(("The map coefficients file '%s' does not exist or is "+
          "not a file.") % xtal.map_coeffs)
      n_maps += 1
    if (xtal.solvent_content is None) or (xtal.solvent_content <= 0) :
      raise Sorry("Solvent content must be greater than 0 and less than 1.")
  if (n_maps == 0) :
    raise Sorry("You must provide map coefficients for at least one crystal "+
      "form.")
  if (params.averaging.cycles is None) :
    raise Sorry("Please specify the number of cycles of density modification "+
      "to run.")
  if (params.control.base_gui_dir is None) :
    raise Sorry("Please specify an output directory.")
  elif (not os.path.isdir(params.control.base_gui_dir)) :
    raise Sorry("The path '%s' does not exist or is not a directory." %
      params.control.base_gui_dir)
  return True

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    if (not os.path.exists(self.output_dir)) :
      os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    mca = multi_crystal_average(args=list(self.args), out=sys.stdout)
    files = getattr(mca, "map_coeffs", None)
    if (files is not None) :
      file_names = []
      for file_name in files :
        file_names.append(os.path.abspath(file_name))
      return file_names
    return None
########################################################################

# now do it
import sys
if __name__=="__main__":
  argument_list=sys.argv[1:]
  if True or is_debug(argument_list).value:

    sys_stdout_sav=sys.stdout
    my_multi_class=multi_crystal_average(args=argument_list)
    sys.exit(0)
  try:
    sys_stdout_sav=sys.stdout
    my_multi_class=multi_crystal_average(args=argument_list)
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
    else:
      # 2011-10-10
      # Print out note on format change for users with old input files
      print "\nNOTE: Format of inputs for pdb_file, datafile, map_coeffs, \nsolvent_content etc have changed (again). Now they are grouped by crystal.  \nSo now you have N crystal groups, each with \npdb_file=pdb_1.pdb solvent_content=0.4 etc.\n"
