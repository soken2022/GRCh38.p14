

fasta_path_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
fasta_path_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa'

input_chromosomes_name =    ['22']# ['01','02','03','04','05','06','07','08','09' ,'10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22']

import sys

from my_functions import open_fasta_file_UCSC_multi_record
for i in input_chromosomes_name:
    ucsc_seq = list(open_fasta_file_UCSC_multi_record(fasta_path_UCSC,i))

    from my_functions import open_fasta_file_UCSC_multi_record
    fasta_path_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa'
    annotation_file_path_out_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out.gz'
    annotation_file_path_gff_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
    from my_functions import  open_annotation_file_and_prepare_it
    anottations = open_annotation_file_and_prepare_it(annotation_file_path_gff_refseq,i)
    element_name = ['LINE','SINE','LTR','DNA'] #study genomic elements
    anottations.readline()
    anottations.readline()
    anottations.readline()
    anottations.readline()
    anottations.readline()
    anottations.readline()
    anottations.readline()
    anottations.readline()
    anottations.readline()
    anottations.readline()
    from my_functions import target_sequences_lowering
    for annotation in anottations:
        annotation = annotation.split()
        if annotation[2] == 'gene':
            print(annotation)
            ucsc_seq = target_sequences_lowering(ucsc_seq,annotation[3],annotation[4])


    sys.modules[__name__].__dict__.clear()


from  my_functions import seq_Nucleotide_composition
seq_Nucleotide_composition(ucsc_seq)





