#!/usr/bin/env python
import pandas as pd
from Bio import SearchIO

df = pd.DataFrame(columns=['contig_id', 'e-value'])
df = df.set_index('contig_id')
hmmFile = snakemake.input[0] 
with open (hmmFile,'r') as handle: 
    for record in SearchIO.parse(handle, 'hmmer3-text'):
        for hit in (record.hits):
            df.loc[hit.id]=hit.evalue
df.to_csv(snakemake.output[0], sep='\t')
