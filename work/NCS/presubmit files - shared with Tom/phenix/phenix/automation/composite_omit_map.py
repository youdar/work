
"""
Module for computing composite omit maps using a method similar to the
equivalent script in CNS.  An alternative (and much faster) method is
implemented in mmtbx.maps.composite_omit_map; both of these are used by the
program phenix.composite_omit_map.
"""

from __future__ import division
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out
from libtbx import easy_mp
from libtbx import Auto
import libtbx.phil
import os.path
import time
import sys

omit_map_phil = libtbx.phil.parse("""
  ncs_average = False
    .type = bool
    .short_caption = NCS average composite maps
    .help = NCS average composite maps
  box_cushion_radius = 2.5
    .type = float
    .help = Additional buffer around omitted region within which to set \
      occupancies to zero.
  exclude_bulk_solvent = False
    .type = bool
    .expert_level = 1
    .short_caption = Exclude bulk solvent from omit regions
    .help = Prevent bulk solvent from extending into omitted regions.
  optimize_binning = True
    .type = bool
    .help = Re-group omit regions to reduce the total number of regions while \
      keeping the fraction omitted per region as similar as possible.
  refinement {
    annealing_type = *cartesian torsion
      .type = choice
      .short_caption = Annealing method
    annealing_temperature = Auto
      .type = int
      .help = Defaults to 5000K for Cartesian, 2500K for torsion
    ncs = Auto
      .type = bool
      .short_caption = Use NCS restraints in refinement
      .help = Use NCS restraints in refinement
    harmonic_restraints = True
      .type = bool
      .short_caption = Harmonic restraints on omitted atoms
      .help = Apply harmonic restraints to zero-occupancy atoms to prevent \
        other atoms outside the omitted region from filling in the missing \
        scattering.
    harmonic_restraints_sigma = 0.1
      .type = float
      .short_caption = Standard deviation for harmonic restraints.
    vary_random_seed = False
      .type = bool
      .help = Vary the random seed between runs.
    automatic_linking = False
      .type = bool
      .help = Automatically link covalent ligands (sugars, etc.).
  }
  ncs_averaging {
    include scope mmtbx.map_tools.ncs_averaging_params
  }
  debugging
    .help = Debugging options, not for general use.
    .expert_level = 4
    .style = hidden
  {
    flatten_background = False
      .type = bool
      .expert_level = 4
      .style = hidden
    control_map = False
      .type = bool
      .expert_level = 4
      .style = hidden
  }
""")

class generate_composite_omit_maps (object) :
  def __init__ (self,
      fmodel,
      xray_structure,
      pdb_hierarchy,
      selection,
      params,
      fraction_omit=0.05,
      map_types=("2mFo-DFc",),
      parallel_params=None,
      tmp_file_base="tmp_omit",
      omit_type="refine",
      refine_args_base=(),
      debug=False,
      resolution_factor=0.25,
      processed_pdb_file=None,
      exclude_free_r_reflections=True,
      fill_missing_f_obs=False,
      out=sys.stdout) :
    assert (str(omit_type) in ["refine", "anneal", "None"])
    self.fmodel = fmodel
    self.xray_structure = xray_structure
    self.pdb_hierarchy = pdb_hierarchy
    self.selection = selection
    self.map_types = map_types
    self.params = params
    self.tmp_file_base = tmp_file_base
    self.refine_args_base = refine_args_base
    self.debug = debug
    self.resolution_factor = resolution_factor
    self.omit_type = omit_type
    self.maps_phil_str = ""
    self.exclude_free_r_reflections = exclude_free_r_reflections
    self.fill_missing_f_obs = fill_missing_f_obs
    self.fraction_omit = fraction_omit
    import mmtbx.maps.composite_omit_map
    if (self.omit_type in ["refine", "anneal"]) :
      mtz_data = fmodel.f_obs().as_mtz_dataset(
        column_root_label="F")
      mtz_data.add_miller_array(fmodel.r_free_flags(),
        column_root_label="FreeR_flag")
      if (fmodel.hl_coeffs() is not None) :
        mtz_data.add_miller_array(fmodel.hl_coeffs(),
          column_root_label="HL")
      mtz_data.mtz_object().write("%s_data.mtz" % self.tmp_file_base)
      print >> out, "Wrote data to %s_data.mtz" % self.tmp_file_base
      self.maps_phil_str = ""
      for i_map, map_type in enumerate(map_types) :
        self.maps_phil_str += """
          refinement.electron_density_maps.map_coefficients {
            map_type = %s
            mtz_label_amplitudes = OMIT%d
            mtz_label_phases = PHOMIT%d
            fill_missing_f_obs = %s
            exclude_free_r_reflections = %s
            ncs_average = %s
          }""" % (map_type, i_map, i_map, self.fill_missing_f_obs,
            self.exclude_free_r_reflections, False)
    if (params.exclude_bulk_solvent) :
      fmodel.mask_params.ignore_zero_occupancy_atoms = False
    self.final_map_coeffs = []
    self.single_omit_mode = False
    if (fraction_omit == 1.0) :
      assert (selection is not None) and (selection.count(True) > 0)
      self.single_omit_mode = True
      region = mmtbx.maps.composite_omit_map.omit_regions(
        serial=1,
        selection=selection)
      self.omit_groups = [ region ]
      print >> out, "Only a single omit region will be used."
    else :
      make_sub_header("Generating omit regions", out=out)
      self.omit_groups = mmtbx.maps.composite_omit_map.create_omit_regions(
        xray_structure=xray_structure,
        selection=selection,
        fraction_omit=fraction_omit,
        optimize_binning=params.optimize_binning,
        box_cushion_radius=params.box_cushion_radius,
        log=out)
      for group in self.omit_groups :
        group.show(out=out)
    self.ncs_average = None
    make_sub_header("Calculating individual omit maps", out=out)
    self.fmodel.info().show_targets(out=out)
    t1 = time.time()
    print >> out, ""
    def show_omit_rfactors (result) :
      print >> out, "  omit region %3d: r_work=%6.4f r_free=%6.4f" % \
        (result.serial, result.r_work, result.r_free)
    all_omit_map_coeffs = easy_mp.parallel_map(
      func=self,
      iterable=self.omit_groups,
      params=parallel_params,
      preserve_order=True,
      callback=show_omit_rfactors,
      preserve_exception_message=True)
    t2 = time.time()
    print >> out, "  Total runtime: %.1fs" % (t2 - t1)
    assert len(all_omit_map_coeffs) == len(self.omit_groups)
    self.fmodel.update_xray_structure(self.xray_structure,
      update_f_mask=True,
      update_f_calc=True)
    if (not self.debug) :
      mtz_tmp = self.tmp_file_base + "_data.mtz"
      if os.path.isfile(mtz_tmp) :
        os.remove(mtz_tmp)
    if (params.ncs_average) :
      assert (processed_pdb_file is not None)
      self.processed_pdb_file = processed_pdb_file
      make_sub_header("Setting up NCS for map averaging", out=out)
      self.setup_ncs_averaging(out=out)
    make_sub_header("Assembling omitted regions", out=out)
    for i_map, map_type in enumerate(self.map_types) :
      print >> out, "Compositing %s map..." % map_type
      if (self.single_omit_mode) :
        assert (len(all_omit_map_coeffs) == 1)
        composite_map_coeffs = all_omit_map_coeffs[0].map_coeffs_list[i_map]
      else :
        omit_map_coeffs = []
        for omit_maps in all_omit_map_coeffs :
          omit_map_coeffs.append(omit_maps.map_coeffs_list[i_map])
        background_map_coeffs = fmodel.map_coefficients(
          map_type=map_type,
          fill_missing=self.fill_missing_f_obs,
          exclude_free_r_reflections=self.exclude_free_r_reflections,
          merge_anomalous=True)
        composite_map = mmtbx.maps.composite_omit_map.combine_maps(
          map_arrays=omit_map_coeffs,
          omit_groups=self.omit_groups,
          background_map_coeffs=background_map_coeffs,
          resolution_factor=self.resolution_factor,
          flatten_background=params.debugging.flatten_background,
          control_map=params.debugging.control_map)
        if (debug) :
          from iotbx import ccp4_map
          from scitbx.array_family import flex
          gridding_first = (0,0,0)
          gridding_last = tuple(composite_map.focus())
          ccp4_map.write_ccp4_map(file_name="composite_%s.ccp4" % map_type,
            unit_cell=background_map_coeffs.unit_cell(),
            space_group=background_map_coeffs.space_group(),
            gridding_first=(0,0,0),
            gridding_last=gridding_last,
            map_data=composite_map,
            labels=flex.std_string(["composite %s map" % map_type]))
        composite_map_coeffs = \
          background_map_coeffs.structure_factors_from_map(
            map=composite_map,
            use_sg=True)
      if (self.ncs_average is not None) :
        print >> out, "  running NCS averaging..."
        composite_map_coeffs = self.ncs_average(
          map_coeffs=composite_map_coeffs,
          fmodel=self.fmodel)
      self.final_map_coeffs.append(composite_map_coeffs)

  def __call__ (self, *args) :
    return self.generate_single_omit_maps(*args)

  def generate_single_omit_maps (self, omit_group) :
    omit_map_coeffs = []
    r_work = r_free = None
    if (str(self.omit_type) == "None") :
      fmodel_tmp = self.fmodel.deep_copy()
      xrs_omit = self.xray_structure.deep_copy_scatterers()
      xrs_omit.set_occupancies(0., selection=omit_group.selection)
      fmodel_tmp.update_xray_structure(xrs_omit,
        update_f_mask=False,
        update_f_calc=True)
      r_work = fmodel_tmp.r_work()
      r_free = fmodel_tmp.r_free()
      for map_type in self.map_types :
        map_coeffs = fmodel_tmp.map_coefficients(
          map_type=map_type,
          fill_missing=self.fill_missing_f_obs,
          exclude_free_r_reflections=self.exclude_free_r_reflections,
          merge_anomalous=True)
        omit_map_coeffs.append(map_coeffs)
      if (self.debug) :
        pdb_out = "region_%d.pdb" % omit_group.serial
        atoms = self.pdb_hierarchy.atoms()
        occ = atoms.extract_occ()
        occ_tmp = occ.deep_copy()
        occ_tmp.set_selected(omit_group.selection, 0)
        atoms.set_occ(occ_tmp)
        open(pdb_out, "w").write(self.pdb_hierarchy.as_pdb_string(
          crystal_symmetry=self.xray_structure))
        atoms.set_occ(occ)
    else :
      refine_args = self.refine_args_base + [
        "%s_data.mtz" % self.tmp_file_base,
        "main.ncs=%s" % self.params.refinement.ncs,
        "output.prefix=%s" % self.tmp_file_base,
        "output.serial_format=%d",
        "output.serial=%d" % omit_group.serial,
        "write_def_file=False",
        "write_geo_file=False",
        "--quiet",
        "--overwrite",
      ]
      if (self.params.refinement.automatic_linking) :
        refine_args.append("link_all=True")
      if (self.fmodel.hl_coeffs() is not None) :
        refine_args.append("main.use_experimental_phases=True")
      if (self.params.refinement.vary_random_seed) :
        refine_args.append("main.random_seed=None")
      if (self.omit_type == "anneal") :
        if (self.params.refinement.annealing_type == "torsion") :
          refine_args.append("main.simulated_annealing_torsion=True")
          if isinstance(self.params.refinement.annealing_temperature, int) :
            refine_args.append("tardy.start_temperature_kelvin=%d" %
              self.params.refinement.annealing_temperature)
        else :
          refine_args.append("main.simulated_annealing=True")
          if isinstance(self.params.refinement.annealing_temperature, int) :
            refine_args.append("simulated_annealing.start_temperature=%d" %
              self.params.refinement.annealing_temperature)
      if (self.params.exclude_bulk_solvent) :
        refine_args.append("refinement.mask.ignore_zero_occupancy_atoms=False")
      if (self.params.refinement.harmonic_restraints) :
        refine_args.append("main.harmonic_restraints=True")
        refine_args.append("harmonic_restraints.selection='occupancy = 0'")
        refine_args.append("harmonic_restraints.sigma=%f" %
          self.params.refinement.harmonic_restraints_sigma)
      atoms = self.pdb_hierarchy.atoms()
      occ = atoms.extract_occ()
      occ_tmp = occ.deep_copy()
      occ_tmp.set_selected(omit_group.selection, 0)
      atoms.set_occ(occ_tmp)
      pdb_file = "%s_%d_in.pdb" % (self.tmp_file_base, omit_group.serial)
      open(pdb_file, "w").write(self.pdb_hierarchy.as_pdb_string(
        crystal_symmetry=self.xray_structure))
      refine_args.append(pdb_file)
      atoms.set_occ(occ)
      from phenix.refinement import command_line
      refine_obj = command_line.run(
        command_name="phenix.refine",
        args=refine_args,
        custom_maps_phil_str=self.maps_phil_str,
        out=null_out(),
        replace_stderr=False,
        skip_write_data_mtz_file=True)
      r_work = refine_obj.fmodel.r_work()
      r_free = refine_obj.fmodel.r_free()
      file_prefix = "%s_%d." % (self.tmp_file_base, omit_group.serial)
      mtz_out = file_prefix + "mtz"
      from iotbx import file_reader
      mtz_in = file_reader.any_file(mtz_out)
      for miller_array in mtz_in.file_server.miller_arrays :
        labels = miller_array.info().labels
        if (labels[0].startswith("OMIT")) :
          omit_map_coeffs.append(miller_array)
      if (not self.debug) :
        os.remove(pdb_file)
        # clean up refinement output
        for file_name in os.listdir(os.getcwd()) :
          if file_name.startswith(file_prefix) :
            os.remove(file_name)
    assert len(omit_map_coeffs) == len(self.map_types)
    if (self.debug) :
      import iotbx.mtz
      mtz_out = "region_%d.mtz" % omit_group.serial
      mtz_dataset = omit_map_coeffs[0].as_mtz_dataset(
        column_root_label="FWT",
        label_decorator=iotbx.mtz.ccp4_label_decorator())
      for k, map_coeffs in enumerate(omit_map_coeffs[1:]) :
        mtz_dataset.add_miller_array(map_coeffs,
          column_root_label="OMIT%d" % (k+1))
      mtz_dataset.mtz_object().write(mtz_out)
    import mmtbx.maps.composite_omit_map
    return mmtbx.maps.composite_omit_map.omit_region_results(
      map_coeffs_list=omit_map_coeffs,
      r_work=r_work,
      r_free=r_free,
      serial=omit_group.serial)

  def write_mtz_file (self, file_name) :
    import iotbx.map_tools
    return iotbx.map_tools.write_map_coefficients_generic(
      file_name=file_name,
      map_types=self.map_types,
      map_coeffs=self.final_map_coeffs)

  def setup_ncs_averaging (self, out=sys.stdout) :
    from phenix.command_line import simple_ncs_from_pdb
    from mmtbx import map_tools
    ncs_object = None
    scattering_types = \
      self.xray_structure.scatterers().extract_scattering_types()
    exclude_d = (scattering_types == "D").count(True) > 0
    exclude_h = (scattering_types == "H").count(True) > 0
    # XXX why doesn't this work with just a hierarchy object?
    ncs_search = simple_ncs_from_pdb.simple_ncs_from_pdb(
      params               = None,
      pdb_inp              = self.pdb_hierarchy.as_pdb_input(),
      source_info          = "pdb_hierarchy",
      hierarchy            = self.pdb_hierarchy,
      suppress_print       = False,
      exclude_h            = exclude_h,
      exclude_d            = exclude_d,
      log                  = out)
    self.ncs_average = map_tools.ncs_averager(
      ncs_object=ncs_search.ncs_object,
      params=self.params.ncs_averaging,
      log=out)
