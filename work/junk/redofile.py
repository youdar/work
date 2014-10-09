from __future__ import division
import os

fn = r'C:\Phenix\Dev\Work\work\junk\origin_shifted.pdb'
data = open(fn,'r').read().splitlines()
data = [x for x in data if x.startswith('ATOM') and (' 217 ' in x) and ('GLY' in x)]


fn = r'C:\Phenix\Dev\Work\work\junk\test.pdb'
data = '\n'.join(data)
open(fn,'w').write(data)

print 'Done'
