from __future__ import division
import collect_ncs_files

__author__ = 'Youval'


""" update all PDB structure data, by collecting refinements results """
c = collect_ncs_files.ncs_paper_data_collection()
c.collect_refinement_results()
print 'Done...'
