Collect non-bonded overlaps and overlaps per 1000 atoms

following tools:
- The test is done with hydrogen atoms

For a single file
-----------------
>>> python nonbonded_pdb_data_collect.py xxxx.pdb [-v] [verbose]

Returns:
  out_str (str):
    if verbose is False: x1,x2,x3,x,4,x5,x6,x7,x8,x9

      x1: PDB ID
      x2: Macro molecule overlaps
      x3: Symmetry overlaps
      x4: All overlaps
      x5: Macro molecule overlaps per 1000 atoms
      x6: Symmetry overlaps per 1000 atoms
      x7: All overlaps per 1000 atoms
      x8: year model deposited in PDB
      x9: experiment type

    if verbose is True, the sting will look like
      PDB id| macro mol. | symmetry  |  all   | year  | experiment type
      -------------------------------------------------------------------
       1a18 |     29     |     5     |   35   | 1997  |X-RAY DIFFRACTION

  Error Types:
    -1: Other processing issue
    -2: model contains unknown_type_pairs
    -3: multiple models in pdb files
    -4: Bad CRYST1 records, bad crystal symmetry
    -5: File could not be fetched

All PDB survey
--------------
To perform a PDB wide test, comparing all files in the LBL PDB mirror
use
>>> python submit_overlap_all_pdb_to_queue.py

Collecting and looking at results
---------------------------------
Use collecting_overlap_data.py to collect and look at the test results

- if "test_clean_data" and "test_clean_dict" exist, plotting is done on
  existing data

- if "test_clean_data" and "test_clean_dict" do not exist, collect info from the
  folder containing the tests, from the queue and create them.

- the file "test_data.txt"  is a comma separated file, containing all data.
- the file "test_clean_data" is a pickled list containing only files that
  could be processed (without error of any kind).
- the file "test_clean_dict" is a pickled dictionary containing all data.

- "test_data.txt" is a comma separated file with the header:
  PDB ID,Macro molecule overlaps,Symmetry overlaps,All overlaps,
  Macro molecule overlaps per 1000 atoms,Symmetry overlaps per 1000 atoms,
  All overlaps per 1000 atoms,year model deposited in PDB,experiment type





