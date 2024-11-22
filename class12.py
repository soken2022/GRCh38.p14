sequence1 = 'AAATACCGTCTTTTGATCAATTGTATCGTCTTGTTCTGTCCTGCGTGGGCGAATATTTTGGGCTCGTGAGCTTCACGGAT'
sequence2 = 'AAATACCGTCTTTTCATCAATTGTATCGTCTTGTTCTGTCCTGCGTGGGCGAATATTTTGGGCTCGTGAGCTTCACGGAT'

start = 10
end = 15

gene1 = sequence1[start:end]
gene2 = sequence2[start:end]
print(gene1)
print(gene2)


e_coli_genome_path = r'D:\PhD\genomes\E coli K-12\ncbi_dataset\data\GCF_000005845.2\GCF_000005845.2_ASM584v2_genomic.fna'

genome_file = open(e_coli_genome_path)
genome_sequence = genome_file.read()
genome_sequence = genome_sequence.replace('\n','')
#print(genome_sequence)
#print(len(genome_sequence))

list_length = []
map_path = r'D:\PhD\genomes\E coli K-12\ncbi_dataset\data\GCF_000005845.2\genomic.gff'
file_mape = open(map_path)
for line in file_mape:
    if '#' not in line:
        line = line.split()
        if line[2] == 'gene':
            start = int(line[3])
            end = int(line[4])
            gene_length = end - start
            list_length.append(gene_length)



print(list_length)
print(max(list_length))
print(min(list_length))

sum_of_list = sum(list_length)

average = sum_of_list/len(list_length)

print('average',average)



