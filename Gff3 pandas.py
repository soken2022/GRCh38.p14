import numpy
import numpy as np
import gffpandas.gffpandas as gp
import gff3
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import plot
bac = r'A:\PhD\genomes\sequence\bacteri\Streptococcus pyogenes strain NCTC12064.gff3'
human = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
data = gp.read_gff3(human)
filtered1 = data.filter_feature_of_type(['gene'])
filtered2 = data.filter_by_length(500,1000)
filtered3 = data.filter_feature_of_type(['+'])
mRNA=filtered3.df.to_numpy()

for i in mRNA:
    print(i)