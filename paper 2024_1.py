

# input information
input_organism_name = 'GRCh38.p14'
refseq_Assembly_ID = 'GCF_000001405.40'

input_plot_len = 8000
annotation_file_path_out_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out'
annotation_file_path_out_refseq = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_rm.out'
fasta_path_refseq = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
annotation_file_path_GFF = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
fasta_path_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.gz'



#mus
fasta = ''





from functions_for_primats import open_fasta_file_multi_record_UCSC_and_refseq
from functions_for_primats import portion_of_TE_elements_in_geneme_with_UCSC_and_refseq

Chr = ['6','NC_000006.12']
chr_sequence_refseq = open_fasta_file_multi_record_UCSC_and_refseq(fasta_path_refseq,Chr)
'''Plus_strand,Complement_strand,Both = portion_of_TE_elements_in_geneme_with_UCSC_and_refseq(chr_sequence_refseq,annotation_file_path_out_refseq,Chr)







print(len(chr_sequence_refseq))
kmer1 = chr_sequence_refseq.count('AATCGC')#kmer
kmer2 = chr_sequence_refseq.count('AATCCC')#sub
kmer3 = chr_sequence_refseq.count('TATCGC')#sub2
kmer4 = chr_sequence_refseq.count('GCGATT')#RC
kmer5 = chr_sequence_refseq.count('CGCTAA')#mirror
kmer6 = chr_sequence_refseq.count('TTAGCG')#C


print('kmer',kmer1)
print('sub',kmer2)
print('sub2',kmer3)
print('RC',kmer4)
print('mirror',kmer5)
print('C',kmer6)


print('plus')


kmer1 = Plus_strand.count('AATCG')#kmer
kmer2 = Plus_strand.count('AATCC')#sub
kmer3 = Plus_strand.count('TATCG')#sub2
kmer4 = Plus_strand.count('CGATT')#RC
kmer5 = Plus_strand.count('GCTAA')#mirror
kmer6 = Plus_strand.count('TTAGC')#C

kmer1c = Complement_strand.count('AATCG')#kmer                                       =
kmer2c = Complement_strand.count('AATCC')#sub
kmer3c = Complement_strand.count('TATCG')#sub2
kmer4c = Complement_strand.count('CGATT')#RC
kmer5c = Complement_strand.count('GCTAA')#mirror
kmer6c = Complement_strand.count('TTAGC')#C



print('kmer',kmer1   ,kmer1c )
print('sub',kmer2    ,kmer2c )
print('sub2',kmer3   ,kmer3c )
print('RC',kmer4     ,kmer4c )
print('mirror',kmer5 ,kmer5c )
print('C',kmer6      ,kmer6c ) '''

from functions_for_primats import d2024_11_18_plus_and_minis_lableing_gff_out

d2024_11_18_plus_and_minis_lableing_gff_out(chr_sequence_refseq,annotation_file_path_GFF,Chr)


