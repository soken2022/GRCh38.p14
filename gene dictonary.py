
def extract_gene_from_gff(gff_file_path, gene_name):
    """
    Extract information about a specific gene from a GFF file.

    Parameters:
    - gff_file_path (str): Path to the GFF file.
    - gene_name (str): Name of the gene to extract.

    Returns:
    - gene_info (list): List of dictionaries containing information about the gene.
    """
    gene_info = []

    with open(gff_file_path, 'r') as gff_file:
        for line in gff_file:
            # Skip comment lines
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            # Check if the feature is a gene
            if len(fields) > 2 and fields[2] == 'gene':
                attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                print(attributes)

                # Check if the gene name matches
                if attributes.get('Name') == gene_name:
                    gene_info.append({
                        'seqid': fields[0],
                        'source': fields[1],
                        'feature': fields[2],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'attributes': attributes
                    })

    return gene_info

# Example usage
gff_file_path = r'A:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
gene_name_to_extract = 'TP53'

result = extract_gene_from_gff(gff_file_path, gene_name_to_extract)

if result:
    print(f"Information about {gene_name_to_extract}:")
    for entry in result:
        print(entry)
else:
    print(f"Gene '{gene_name_to_extract}' not found in the GFF file.")

print(entry['attributes']['Name'])
print(entry['attributes']['Dbxref'])
