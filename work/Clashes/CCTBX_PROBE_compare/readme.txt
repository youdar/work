To compare CCTBX non-bonded clash scores to PROBE clashscore use the
following tools:
- The test is done with hydrogen atoms

For a single file
-----------------
>>> python compare_clashscores_CCTBX_vs_PROBE.py [-v] [verbose]


All PDB survey
--------------
To perform a PDB wide test, comparing all files in the LBL PDB mirror
use
>>> python submit_clashscore_all_pdb_to_queue.py

Collecting and looking at results
---------------------------------
Use collecting_and_looking_at_data.py to collect and look at the test results
- if test_clean_data and test_data_dict exist,
we are just going to see plots and info about the data in that file.
- if test_clean_data and test_data_dict do not exist, collect info from the
folder containing the tests from the queue and create them
- the file test_data contains all files with data that indicates the type of
failure.
  (-1) : PDB file was not found
  (-2) : proper structure factor file was not available
  (-3) : Some other error occurred

Collecting PDB IDs that were not processed
------------------------------------------
After running submit_clashscore_all_pdb_to_queue.py and
collecting_and_looking_at_data.py some files might not have ran.
Use collect_remaining_files.py to collect the list of those files and create
a text file containing that list.
The test file containing that PDB IDs list is: "files_to_run.txt"
PDB IDs of files that were processed are found in the pickled file "test_data"

Note: This must be done on LBL machine, since we are looking at the index of
PDB files on the LBL PDB mirror, "PDB_MIRROR_PDB"


