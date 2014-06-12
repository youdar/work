from __future__ import division
from libtbx.phil import parse

# master_phil = parse("""
#   minimization.input {
#     file_name = None
#       .type = path
#     label = None
#       .type = str
#   }
#   """)
#
# user_phil = parse("""
#   minimization.input {
#     file_name = experiment.dat
#   }
#   """)
#
# command_line_phil = parse("minimization.input.label=set2")
#
# working_phil = master_phil.fetch(
#   sources=[user_phil, command_line_phil])
# working_phil.show()
# print '='*30
# working_params = working_phil.extract()
#
# print working_params.minimization.input.file_name
# print working_params.minimization.input.label
# print '='*30
# working_params.minimization.input.label = "set3"
# modified_phil = master_phil.format(python_object=working_params)
# modified_phil.show()
# print '='*30
#
# master_phil = parse("""
#   random_integers = None
#     .type = ints
#   euler_angles = None
#     .type = floats(size=3)
#   unit_cell_parameters = None
#     .type = floats(size_min=1, size_max=6)
#   rotation_part = None
#     .type = floats(size=9, value_min=-1, value_max=1)
#   """)
#
# user_phil = parse("""
#   random_integers = 3 18 5
#   euler_angles = 10 -20 30
#   unit_cell_parameters = 10,20,30
#   rotation_part = "1,0,0;0,-1,0;0,0,-1"
#   """)
#
# working_phil = master_phil.fetch(source=user_phil)
# working_phil.show()
# print
# working_params = working_phil.extract()
# print working_params.random_integers
# print working_params.euler_angles
# print working_params.unit_cell_parameters
# print working_params.rotation_part
# print
# working_phil = master_phil.format(python_object=working_params)
# working_phil.show()

print '='*30
# master_phil = parse("""
#   ncs_refinement {
#     ncs_group
#     .multiple = True
#       {
#         transform
#         .multiple = True
#         {
#           rotations = None
#             .type = floats(size=9, value_min=-1, value_max=1)
#             .multiple = True
#           translations = None
#             .type = floats(size=3)
#             .multiple = True
#           transform_serial_num = None
#             .type = int
#             .multiple = True
#           coordinates_present = None
#             .type = bool
#         }
#         apply_to_all_chains = True
#           .type = bool
#         apply_to_chains_id = None
#           .type = str
#           .multiple = True
#         dont_apply_when_coordinates_present = False
#           .help ="Apply transformations even if coordinates are already present"
#           .type = bool
#       }
#     }
#   """)

# master_phil = parse("""
#   ncs_refinement {
#     dont_apply_when_coordinates_present = True
#       .type = bool
#       .help='''
#       Can apply transformations even if coordinates are already present.
#       Need to provide NCS selection'''
#     ncs_selection = None
#       .type = str
#       .help = '''
#       When applying transforms to coordinates that are present, you need to
#       provide selection string for the NCS portion of the pdb hierarchy'''
#     apply_to_all_chains = True
#       .type = bool
#       .help = '''
#       If any custom ncs_groups are set, apply_to_all_chains will be
#       forced to be False
#       '''
#     ncs_group
#       .multiple = True
#       {
#         transform
#         .multiple = True
#         {
#           rotations = None
#             .type = floats(size=9, value_min=-1, value_max=1)
#           translations = None
#             .type = floats(size=3)
#           coordinates_present = None
#             .type = bool
#         }
#         apply_to_selection = None
#           .help = '''Selection syntax similar to what is used in CNS or PyMOL.
#           multiple selection lines will be concatenates.'''
#           .type = str
#           .multiple = True
#
#       }
#     }
#   """)

# user_phil = parse("""
#   ncs_refinement {
#     ncs_group {
#       transform {
#         rotations = (1.0,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9)
#         translations = (1,2,3)
#         rotations = 0.1,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9
#         translations = (0,2,1)
#         rotations = 1.0,0.2,1.0,0.2,0.5,0.6,0.7,0.8,0.9
#         translations = (-1,3,-2)
#       }
#       apply_to_chains_id = A
#     }
#     ncs_group {
#       transform {
#         rotations = 0.2,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9
#         translations = (0,0,0)
#         rotations = 0.3,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9
#       }
#       apply_to_chains_id = B
#       apply_to_chains_id = C
#     }
#   }
#   """)

master_phil = parse("""
  ncs_group_selection {
    ncs_group
        .multiple = True
        {
        ncs_group_number = None
          .type = int
        ncs_selection = ''
          .type = str
          .help = 'Residue selection string for a single NCS group'
        asu_selection = ''
          .type = str
          .help = 'Residue selection string for NCS copies location in ASU'
          .multiple = True
    }
  }
  """)






# working_phil = master_phil.fetch(source=user_phil).extract()
master_phil.show()

print '='*30

# master_phil = parse("""
#   minimization.input {
#     file_name = None
#       .type = path
#   }
#   minimization.parameters {
#     max_iterations = 10
#       .type = int
#   }
#   """)
#
# user_phil = parse("""
#   minimization.input.file_name = experiment.dat
#   minimization.parameters.max_iterations = 5
#   """)
#
working_params = master_phil.fetch(source=user_phil).extract()
print working_params
# print working_params.minimization.input.file_name
# print working_params.minimization.parameters.max_iterations
