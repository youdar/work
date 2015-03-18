from __future__ import division
import collect_ncs_files
import pandas
import os




c = collect_ncs_files.ncs_paper_data_collection()
fn = os.path.join(c.ncs_dir,'ncs_paper_data.csv')
print os.path.isfile(fn)
df = pandas.DataFrame.from_csv(fn)
d








