

fasta_path_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
fasta_path_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa'

input_chromosomes_name =    ['24']# ['01','02','03','04','05','06','07','08','09' ,'10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22']


from my_functions import open_fasta_file_multi_record_from_refseq
from my_functions import open_fasta_file_UCSC_multi_record
for i in input_chromosomes_name:
    ucsc_seq = open_fasta_file_UCSC_multi_record(fasta_path_UCSC,i)
    refseq  = open_fasta_file_multi_record_from_refseq(fasta_path_refseq,i)
    import sys

    from my_functions import open_fasta_file_multi_record_from_refseq
    from my_functions import open_fasta_file_UCSC_multi_record
    fasta_path_refseq = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
    fasta_path_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa'
    annotation_file_path_out_UCSC = r'A:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out.gz'
    from my_functions import  open_annotation_file_and_prepare_it
    anottations = open_annotation_file_and_prepare_it(annotation_file_path_out_UCSC,i)
    element_name = ['LINE','SINE','LTR','DNA'] #study genomic elements
    anottations.readline()
    anottations.readline()
    anottations.readline()
    for annotation in anottations:
        annotation = annotation.split()
        if (element_name[0] in annotation[10] or element_name[1] in annotation[10] or element_name[2] in annotation[10] or element_name[3] in annotation[10]) and (annotation[4] == f'chr{i}'): # or annotation[4] == f'chr{chromosome_number}_'):
            if annotation[8] == '+':
                sequence_of_element_p = ucsc_seq[int(annotation[5]) - 1:int(annotation[6]) - 1]
                sequence_of_element_refseq = refseq[int(annotation[5]) - 1:int(annotation[6]) - 1]

                print('ucsc',sequence_of_element_p)
                #print(ucsc_seq[int(i[5]):int(i[6])])
                print('refs',sequence_of_element_refseq)

    sys.modules[__name__].__dict__.clear()





