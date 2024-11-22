from gain_gene_sequence import get_gene_sequence
from convert_GeneID_to_Nucleotide_ID import convert_gene_id_to_nucleotide_ids
from Nucleotide_ID_to_Sequence import download_gene_info_from_refseq

def extract_gene_from_gff(gff_file_path):

    gene_info = []

    with open(gff_file_path, 'r') as gff_file:
        for line in gff_file:
            # Skip comment lines
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            len_gene = int(fields[4]) - int(fields[3])
            # Check if the feature is a gene
            if len(fields) > 2 and fields[2] == 'gene':
                #print(type(fields[8]),fields[8])
                #GeneID = fields[8].split(';')[1].split(':')[1].split(',')
                #Gene_Name = fields[8].split(';')[2].split('=')[1]

                attributes = dict(item.split('=') for item in fields[8].split(';'))
                #print('atatribut',attributes)
                alpa = attributes['Dbxref']

                GeneID_dic = dict(i.split(':') for i in alpa.split(','))


                #attributes2 = dict((attributes['Dbxref'].replace(':','=')))
                #print('atatribut2',attributes2)
                #print(attributes.keys())
                GeneID = GeneID_dic['GeneID']
                Gene_Name = fields[8].split(';')[2].split('=')[1]
                #Gene_Name = dict(item.split('=') for item in fields[8].split(';')[2])
                print('AAAAAAAAAAAAAAAAAAAAAAAAAA')
                print('GeneID :',GeneID,'Gene_Name :',Gene_Name)
                print(len_gene)
                print(fields)
                '''if len(GeneID) > 1 :
                    GeneID = [GeneID[0]]
                print('-------------------')
                print(fields)
                print('GeneID',int(GeneID[0]))
                #print('GeneID',int(GeneID[0]),'len',len_gene)
                print('Gene_Name',Gene_Name)
                print(fields)
                #get_gene_sequence_by_id(GeneID)'''
                #get_gene_sequence(Gene_Name)
                print('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')
                Nucleotide_IDs = convert_gene_id_to_nucleotide_ids(GeneID)
                for i in range(len(Nucleotide_IDs)):
                    download_gene_info_from_refseq(Nucleotide_IDs[i])

                attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                #print(attributes.keys())















                ''''# Check if the gene name matches
                if attributes.get('Name') == gene_name:
                    gene_info.append({
                        'seqid': fields[0],
                        'source': fields[1],
                        'feature': fields[2],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'attributes': attributes
                    })'''

    #return gene_info

# Example usage
annotation_GFF_file_path_refseq = r'A:\PhD\genomes\dm6 Drosophila\ncbi_dataset\data\GCF_000001215.4\genomic.gff'
#gene_name_to_extract = 'TP53'

extract_gene_from_gff(annotation_GFF_file_path_refseq)