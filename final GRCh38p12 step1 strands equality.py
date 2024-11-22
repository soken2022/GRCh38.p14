# equality portion between + and - strand

#imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
import gzip
from Bio import SeqIO

# input information
input_organism_name = 'GRCh38.p14'
refseq_Assembly_ID = 'GCF_000001405.40'
input_element_name = ['LINE','SINE','LTR','DNA','RC/Helitron','tRNA','scRNA','Retroposons','RNA','srpRNA'] #study genomic elements
element_name_gff = ['gene','pseudogene']
chromosomes_names = ['8']#['10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22','X','Y']#['10']#['1','2','3','4','5','6','7','8','9' ,'10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22','X','Y']

input_plot_len = 8000
annotation_file_path_out_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out'
annotation_file_path_out_refseq = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_rm.out'
fasta_path_refseq = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
annotation_file_path_GFF = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
fasta_path_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.gz'



from my_functions import  check_and_creat_refseq_file
from functions_for_primats import open_fasta_file_multi_record_UCSC_and_refseq
from functions_for_primats import portion_of_TE_elements_in_geneme_with_UCSC_and_refseq
from functions_for_primats import nucleotide_composition_in_seqenence_simple
from functions_for_primats import nucleotide_composition_in_seqenence_ATCG_2
from functions_for_primats import test
from functions_for_primats import merge_seq_seq_share_with_1234
from functions_for_primats import analysis_composition_of_labled_sequence
from functions_for_primats import only_ATCG_and_atcg

refseq_records_lists = ['NC_000008.11']#['NC_000010.11','NC_000011.10','NC_000012.12','NC_000013.11','NC_000014.9','NC_000015.10','NC_000016.10','NC_000017.11','NC_000018.10','NC_000019.10','NC_000020.11','NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10']# check_and_creat_refseq_file(refseq_Assembly_ID,input_organism_name,chromosomes_names)


#refseq

#analysis for genes

plus_portions = []
comp_portions = []
duble_portions = []
unmask_portions = []
unsequenced = []



for Chr_number,Chr_ID in zip(chromosomes_names,refseq_records_lists):
    Chr = [Chr_number,Chr_ID]
    print(f'start chr{Chr_number,Chr_ID}')
    print('GFF analysis '*10)
    print('step',1)
    Chr_sequence_refseq = open_fasta_file_multi_record_UCSC_and_refseq(fasta_path_refseq,Chr)
    print('step',2)
    plus,comp,both = test(Chr_sequence_refseq,annotation_file_path_GFF,Chr,element_name_gff)
    '''print('nucleotide_composition_in_seqenence in plus')
    nucleotide_composition_in_seqenence_simple(plus)
    print('nucleotide_composition_in_seqenence in comp')
    nucleotide_composition_in_seqenence_simple(comp)
    print('nucleotide_composition_in_seqenence in both')
    nucleotide_composition_in_seqenence_simple(both)'''


    print('GFF finished '+'*'*50)
    #TEs **********************************************
    print('TEs analysis '+ '*'*50)
    print('step',3)
    Plus_strand,Complement_strand,Both = portion_of_TE_elements_in_geneme_with_UCSC_and_refseq(Chr_sequence_refseq,annotation_file_path_out_refseq,Chr,input_element_name)
    '''print('nucleotide_composition_in_seqenence in plus')
    nucleotide_composition_in_seqenence_simple(Plus_strand)
    print('nucleotide_composition_in_seqenence in comp')
    nucleotide_composition_in_seqenence_simple(Complement_strand)
    print('nucleotide_composition_in_seqenence in both')
    nucleotide_composition_in_seqenence_simple(Both)'''


    print('TEs finished '+'*'*50)
    #nucleotide_composition_in_seqenence_atcg(Plus_strand)
    #nucleotide_composition_in_seqenence_atcg(Complement_strand)
    #nucleotide_composition_in_seqenence_atcg(Both)
    print('GFF and TEs analysis'+'*' * 50)
    print('step',4)
    merged_sequence_plus = merge_seq_seq_share_with_1234(plus,Plus_strand)
    print('information about merging of GFF and TEs plus information')
    print('step',5)
    print('GFF+ and TE+ on plus strand ')
    one_masked_portion_type1_p,one_masked_portion_type2_p,Duble_masked_portion_p,unmasked_region_p,unsequenced_region_p = nucleotide_composition_in_seqenence_ATCG_2(merged_sequence_plus)
    print('step',6)
    plus_sequence = only_ATCG_and_atcg(merged_sequence_plus)
    print('step',7)
    print('GFF- and TE- on minis strand')
    merged_sequence_comp = merge_seq_seq_share_with_1234(comp,Complement_strand)
    print('result about merging of GFF and TEs comp information')
    print('step',8)
    one_masked_portion_type1_c,one_masked_portion_type2_c,Duble_masked_portion_c,unmasked_region_c,unsequenced_region_c = nucleotide_composition_in_seqenence_ATCG_2(merged_sequence_comp)
    print('step',9)
    comp_sequence = only_ATCG_and_atcg(merged_sequence_comp)
    print('step',10)
    print('Both (+-)GFF and (+-)TEs on one strand')
    merged_sequence_both = merge_seq_seq_share_with_1234(both,Both)
    print('information about merging of GFF and TEs both strand information')
    print('step',11)
    one_masked_portion_type1_B,one_masked_portion_type2_B,Duble_masked_portion_B,unmasked_region_B,unsequenced_region_B = nucleotide_composition_in_seqenence_ATCG_2(merged_sequence_both)
    print('step',12)
    both_sequecne = only_ATCG_and_atcg(merged_sequence_both)
    print('analysis_composition_of_labled_sequence_final')
    print('step',13)
    P,C,Sh,U = analysis_composition_of_labled_sequence(plus_sequence,comp_sequence,both_sequecne)

    '''plus_portions.append(P)
    comp_portions.append(C)
    duble_portions.append(Sh)
    unmask_portions.append(U)

    print(f'finish chr{Chr_number,Chr_ID}')
    Chr_sequence_refseq = 0
    plus,comp,both = 0,0,0
    Plus_strand,Complement_strand,Both = 0,0,0
    merged_sequence_plus=0
    plus_sequence = 0
    merged_sequence_comp = 0
    comp_sequence = 0
    merged_sequence_both = 0
    both_sequecne = 0
divition_P_to_C = []
for i in range(len(plus_portions)):
    divition_P_to_C.append(plus_portions[i]/comp_portions[i])
import matplotlib.pyplot as plt2
plt2.plot(divition_P_to_C)
plt2.plot(duble_portions)
plt2.plot(unmask_portions)
plt2.show()
print(divition_P_to_C)
print(plus_portions)
print(comp_portions)'''








#UCSC

'''for Chr_number,Chr_ID in zip(chromosomes_names,refseq_records_lists):
    Chr = [Chr_number,Chr_ID]
    chr_sequence_UCSC = open_fasta_file_UCSC_multi_record(fasta_path_UCSC,Chr)
    P,C,B = portion_of_TE_elements_in_geneme_with_UCSC(chr_sequence_UCSC,annotation_file_path_out_UCSC,Chr,input_element_name)
    nucleotide_composition_in_seqenence_atcg(P)
    nucleotide_composition_in_seqenence_atcg(C)
    nucleotide_composition_in_seqenence_atcg(B)'''




'''for Chr in chromosomes_names:
    open_fasta_file_multi_record_from_refseq(fasta_path_refseq,Chr)'''