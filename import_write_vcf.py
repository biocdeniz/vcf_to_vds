#!/bin/env python


import os, sys
import datetime

from hail import *
hc = HailContext()

out = sys.argv[-1]
vcfs = sys.argv[:-1]

# Read the specified .vcf file into memory from the specified path, if `force`
# is set to true, it will attempt to import a gzipp-ed file.
# See: https://www.hail.is/docs/stable/hail.HailContext.html#hail.HailContext.import_vcf
vds = hc.import_vcf(vcfs, force=True)
print("%s Completed reading %s file.\n".format(datetime.datetime.now(), sys.argv[1]))

#table = hc.import_table('s3://avl-hail-dev/laura_stuff/laura_id.txt').key_by('Sample')
#vds = vds.annotate_samples_table(table, root='sa.id')

# Write the converted VDS object to disk. If `parquet_genotypes` is True,
# third party compatible Parque information will be written out as well.
# See: https://www.hail.is/docs/stable/hail.VariantDataset.html#hail.VariantDataset.write
vds.write(out)
print("%s Wrote %s file.\n".format(datetime.datetime.now(), sys.argv[2]))

hc.stop()