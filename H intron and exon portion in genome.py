#Exons = gene - introns
#CDS = gene - introns - UTRs
#therefore also:
#CDS = Exons - UTRs

# for exon simple calculat
#imports

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
import gzip
from Bio import SeqIO
#from memory_profiler import profile

# input information
input_organism_name = 'GRCh38.p14'
input_element_name = 'gene' #study genomic elements
input_element_name_for_saving = '_'.join(input_element_name)
input_start_chromosome_number = 1
input_end_chromosome_number = 24
input_number_of_chromosomes = input_end_chromosome_number - input_start_chromosome_number + 1
input_chromosomes_name = '01'#['01','02','03','04','05','06','07','08','09' ,'10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22']
input_plot_len = 1000000
annotation_file_path_out_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out.gz'
annotation_file_path_out_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
fasta_path_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
fasta_path_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.gz'


from my_functions import open_fasta_file_multi_record_from_refseq
from my_functions import calculat_portion_of_elements
from my_functions import seq_Nucleotide_composition
sequence = open_fasta_file_multi_record_from_refseq(fasta_path_refseq,input_chromosomes_name)

P,C,B = calculat_portion_of_elements(sequence,annotation_file_path_out_refseq,input_chromosomes_name,input_element_name)

seq_Nucleotide_composition(P)

'''A 25.63 T 25.62 C 22.38 G 22.6 N 11658686
a 0.88 t 0.86 c 1.01 g 1.01 n 0
all nucleotide 99.99000000000001'''

'''A 26.29 T 26.3 C 23.11 G 23.33 N 11658686
a 0.22 t 0.18 c 0.29 g 0.28 n 0
all nucleotide 100.0'''
