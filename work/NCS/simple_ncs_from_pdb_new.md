Finding NCS in chains from a PDB file with simple\_ncs\_from\_pdb
=================================================================

Author(s)
---------

-  simple\_ncs\_from\_pdb : Tom Terwilliger
-  simple\_ncs\_from\_pdb : Youval Dar
-  Phil command interpreter: Ralf W. Grosse-Kunstleve
-  find\_domain: Peter Zwart

Purpose
-------

The simple\_ncs\_from\_pdb method identifies NCS in the chains in a PDB
file and writes out the NCS operators in forms suitable for
phenix.refine, resolve, and the AutoSol and AutoBuild Wizards.

Usage
-----

How simple\_ncs\_from\_pdb works:
---------------------------------

The basic steps that the simple\_ncs\_from\_pdb carries out are:

-  (1) Identify sets of matching segments in chains in the PDB file by
   sequences. These are potential NCS-related chains.

-  (2) Determine which matching segments have overall rms distance (RMSD) within
   the given tolerance (max\_rmsd, typically 2 A)

-  (3) Remove residues from matching segments if locally spatially
   misaligned (match\_radius, exclude\_misaligned\_residues)

-  (4) Remove atoms from matching resides if not in both residues (for
   example, if a side chain is missing in one of the matching residues)

-  (5) Optionally, can check atoms order residues

-  (6) Group all chains with common matching segments. (keeping minimal common
   segments)

-  (7) Group Master-Copies pairs that share rotation and translation. Can group
   in a way tha minimize the number of chains in the master copy, or that
   minimize the overall numbers of transforms (rotation & translation)


Additional notes on how simple\_ncs\_from\_pdb works:
-----------------------------------------------------

The chains matching is done using dynamic programming alignment of residues
and atoms. The first pass contains restriction on minimal similarity, set by
min/_percent (number of matching residues)/(number of residues in longer chain)

From the matching chains list, we remove chain pairs where the matching segments
exceed the RMSD limit.

Matching segments are scanned for local residue misalignment. Residues where
(max_atom_distance - min_atom_distance) > match/_radius are excluded from
matching segment. This allow local differences in matching chains.

If matching residues have different number of atoms (For example if one
containing the side chain while the other not), only the matching atoms will
be included.

Grouping of chains done in two steps. First, a single master chain to copies
that share matching segments. Second, if two groups of master and copies
share the same rotations and translations, we group the masters and the copies.

When grouping chains together, we limit the similarity level of chains pairs
being grouped using similarity/_threshold. Lower similarity/_threshold will
allow grouping of more chains.

Similarity groping example: Consider four chains (A,B,C,D) where
A has 100 residues, B has 50, C has 50 and D has 50. Consider that A-B have
50 contiguous matching residues and that C and D have 40. Also assume
that A -> B has the same rotation and translation as C -> D.
(A,B) similarity is 50/100 = 0.5 and (C,D) similarity is 40/50 = 0.8.
The 0.5/0.8 = 0.625, if the similarity/_threshold=0.95 we can't group will get
two groups [(A),(B)] , [(C),(D)]. but if the similarity/_threshold=0.6,
we will get a single group [(A,C),(B,D)] where each NCS copy contains two
chains.

The search for master to group is done using spatial proximity, to improve
the probability that chains that are next to each other will be group together.

Alternative confirmation are excluded from matching segments.

The matching is done by the residues name strings, not by the residue numbers,
this allows handling of insertions in PDB file.

The result of the NCS search is combination of NCS related groups and invariant
or non-NCS related regions, to the atom level. In every NCS group all copies
have the same number of atoms and can be reproduced by applying the applying
the appropriate rotation and translation to the master copy.

*Output files from simple\_ncs\_from\_pdb*
------------------------------------------

When running

::

    phenix.simple_ncs_from_pdb 4boz.pdb

The following files will be produced

-  NCS operators written in format for phenix.refine

::

    4boz_simple_ncs_from_pdb.ncs


-  NCS operators written in format for the PHENIX Wizards

::

    4boz_simple_ncs_from_pdb.ncs_spec
    4boz_simple_ncs_from_pdb.resolve


Examples
--------

Standard run of simple\_ncs\_from\_pdb:
---------------------------------------

Running simple\_ncs\_from pdb is easy. For example, you can type:

::

    phenix.simple_ncs_from_pdb 4boz.pdb

Simple\_ncs\_from\_pdb will analyze the chains in 4boz.pdb and identify
any NCS that exists. For this sample run the following output is
produced:

::

    GROUP 1
    Summary of NCS group with 2 operators:
    ID of chain/residue where these apply: 
    [['A', 'D'], [[[147, 150], [152, 211], [213, 275], [280, 305], 
    [307, 308]], [[147, 150], [152, 211], [213, 275], [280, 305], [307, 308]]]]
    RMSD (A) from chain A:  0.0  0.82
    Number of residues matching chain A:[155, 155]
    
    OPERATOR 1
    CENTER:   24.4880  -13.3177  -20.1848
    
    ROTA 1:    1.0000    0.0000    0.0000
    ROTA 2:    0.0000    1.0000    0.0000
    ROTA 3:    0.0000    0.0000    1.0000
    TRANS:     0.0000    0.0000    0.0000
    
    OPERATOR 2
    CENTER:   15.9430   11.8822    0.6609
    
    ROTA 1:    0.7955   -0.5660    0.2164
    ROTA 2:   -0.5511   -0.8242   -0.1300
    ROTA 3:    0.2520   -0.0159   -0.9676
    TRANS:    18.4021    5.3569  -23.3621
    
    
    GROUP 2
    Summary of NCS group with 2 operators:
    ID of chain/residue where these apply: 
    [['B', 'E'], [[[1, 41], [43, 76]], [[1, 41], [43, 76]]]]
    RMSD (A) from chain B:  0.0  0.88
    Number of residues matching chain B:[75, 75]
    
    OPERATOR 1
    CENTER:   46.1427   -9.2567  -26.1789
    
    ROTA 1:    1.0000    0.0000    0.0000
    ROTA 2:    0.0000    1.0000    0.0000
    ROTA 3:    0.0000    0.0000    1.0000
    TRANS:     0.0000    0.0000    0.0000
    
    OPERATOR 2
    CENTER:   28.3065   -4.0803   11.2523
    
    ROTA 1:    0.7550   -0.5940    0.2778
    ROTA 2:   -0.5894   -0.8004   -0.1096
    ROTA 3:    0.2874   -0.0810   -0.9544
    TRANS:    19.2296    5.3926  -23.9141
    
Another way to view the results without creating files is

::

    phenix.simple_ncs_from_pdb 4boz.pdb show_summary=true

    Chains in model:
    ---------------------------------------------------
    A    B    C    D    E    
    . . . . . . . . . . . . . . . . . . . . . . . . . . 
    
    NCS summary:
    ---------------------------------------------------
    Number of NCS groups     :   2
    Group #                  :   1
    Number of copies         :   2
    Chains in master         :   A
    Chains in copies         :   D
    Group #                  :   2
    Number of copies         :   2
    Chains in master         :   B
    Chains in copies         :   E
    . . . . . . . . . . . . . . . . . . . . . . . . . . 
    
    Transforms:
    ---------------------------------------------------
    Group #                  :   1
    Transform #              :   1
    RMSD                     :   0
    ROTA   0    1.0000    0.0000    0.0000
    ROTA   1    0.0000    1.0000    0.0000
    ROTA   2    0.0000    0.0000    1.0000
    TRANS       0.0000    0.0000    0.0000
    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
    Transform #              :   2
    RMSD                     :   0.817
    ROTA   0    0.7955   -0.5511    0.2520
    ROTA   1   -0.5660   -0.8242   -0.0159
    ROTA   2    0.2164   -0.1300   -0.9676
    TRANS      -5.7994   14.4602  -25.8920
    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
    Group #                  :   2
    Transform #              :   1
    RMSD                     :   0
    ROTA   0    1.0000    0.0000    0.0000
    ROTA   1    0.0000    1.0000    0.0000
    ROTA   2    0.0000    0.0000    1.0000
    TRANS       0.0000    0.0000    0.0000
    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
    Transform #              :   2
    RMSD                     :   0.8756
    ROTA   0    0.7550   -0.5894    0.2874
    ROTA   1   -0.5940   -0.8004   -0.0810
    ROTA   2    0.2778   -0.1096   -0.9544
    TRANS      -4.4667   13.8014  -27.5738


There are 5 chains in the PDB file (A,B,C,D,E). In the first group the master
is A the the copy D and in the second group the master is B and the copy is E.
Chain C is not in any group.

::

    RMSD (A) from chain A:  0.0  0.82 
    
show the RMSD of matching atoms between the master and every other copy, 
The list of numbers is the list of matching residues by residue number. 
Note that this is not the exact selection, as it appears in:
4boz\_simple\_ncs\_from\_pdb.ncs. 

ROTA and TRANS are the rotation and translation information. In the 
.ncs\_spec file and the default on-screen representation of the results, the 
rotation and translation are:

::

    Master = ROT x Copy + TRANS
    
While in the summery and in the CCBTX implementation the coomon use is:

::

    Copy = ROT x Master + TRANS


So the rotation/Translation are the inverse of each other in the two formats.

A portion of the contents of the 4boz\_simple\_ncs\_from\_pdb.ncs\_spec file, 
which you can edit if you want and which you can use in the AutoBuild Wizard,
are shown below. NOTE: The ncs operators describe how to map the N'th
ncs-related copy on to the first copy.

::

    Summary of NCS information
    Thu Apr  2 15:44:03 2015
    /net/cci-filer2/raid1/home/...
    
    new_ncs_group
    new_operator
    
    rota_matrix    1.0000    0.0000    0.0000
    rota_matrix    0.0000    1.0000    0.0000
    rota_matrix    0.0000    0.0000    1.0000
    tran_orth     0.0000    0.0000    0.0000
    
    center_orth   24.4880  -13.3177  -20.1848
    CHAIN A
    RMSD 0
    MATCHING 155
      RESSEQ 147:150
      RESSEQ 152:211
      RESSEQ 213:275
      RESSEQ 280:305
      RESSEQ 307:308
    
    new_operator
    
    rota_matrix    0.7955   -0.5660    0.2164
    rota_matrix   -0.5511   -0.8242   -0.1300
    rota_matrix    0.2520   -0.0159   -0.9676
    tran_orth    18.4021    5.3569  -23.3621
    
    center_orth   15.9430   11.8822    0.6609
    CHAIN D
    RMSD 0.817
    MATCHING 155
      RESSEQ 147:150
      RESSEQ 152:211
      RESSEQ 213:275
      RESSEQ 280:305
      RESSEQ 307:308


The file used for refinement is 4boz\_simple\_ncs\_from\_pdb.ncs. This file 
can also be modified if a particular NCS relations need to be reinforced. The
content of that file is the exact selection sting of the atoms in the NCS groups

The content of 1vcr\_simple\_ncs\_from\_pdb.ncs is

::

    refinement.ncs.restraint_group {
      reference = chain A
      selection = chain B
      selection = chain C
      selection = chain D
      selection = chain E
    }
    
If the following option is used

::
    
    phenix.simple_ncs_from_pdb 1vcr.pdb ncs_file_format=constraints

The content of 1vcr\_simple\_ncs\_from\_pdb.ncs will have constraint\_group
instead of restraint\_group

::

    refinement.ncs.constraint_group {
      reference = chain A
      selection = chain B
      selection = chain C
      selection = chain D
      selection = chain E
    }


Possible Problems
-----------------

Specific limitations and problems:
----------------------------------

-  a large min\_percent can prevent matching segments in chains that are of 
   very different size.

-  While the grouping is done using spatial proximity order, there is no 
   guaranty of finding the absolute optimal grouping.  

-  The value of max\_rmsd, the RMSD limit between chains should be set with 
   care. Small value might exclude chains with problems in the model, even 
   though they are NCS related. Large value might cause non-NCS related 
   chains and to be constraint as related.
   
-  The .ncs\_spec files do not support representation of groups with multiple
   chains in master and copies. a group [(A,C),(B,D)] will be represented as 
   [(A),(B)] and [(C),(D)], duplicating the rotation and translation operations.
   
-  The .ncs\_spec files do not support representation atoms, residue name 
   and insertion selection, it only shows the numbers of the matching residues.
   
-  The .ncs\_spec files rotation and translation operation are

:: 

    Master = ROT x Copy + TRANS
    
   and not:
   
::

    Copy = ROT x Master + TRANS
   
   They are the inverse of the rotation and translation that are used in the 
   implementation of the NCS relation.

Literature
----------

Additional information
----------------------


List of all available keywords
------------------------------

{{phil:phenix.command_line.simple_ncs_from_pdb}}
