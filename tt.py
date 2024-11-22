'''Max_horof_AT=""
Max_adasd_AT=0
Final_Max_N=0
if A>T:
    Max_horof_AT="A"
    Max_adasd_AT=A
else:
    Max_horof_AT="T"
    Max_adasd_AT=T

Max_horof_CG=""
Max_adasd_CG=0
if C>T:
    Max_horof_CG="C"
    Max_adasd_CG=C
else:
    Max_horof_CG="G"
    Max_adasd_CG=G

if Max_adasd_AT> Max_adasd_CG:
    print(Max_horof_AT)
else:
    print(Max_horof_CG)'''


















from Bio import  SeqIO
import matplotlib.pyplot as plt

def path_to_sequence(path,chromosome_number):
    for chromosome in SeqIO.parse(path,"fasta"):
        print(chromosome.id)
        ch_id=chromosome.id
        if chromosome_number in ch_id:
            sequence=str(chromosome.seq)
            break
    return sequence

def sequence_composition(sequence):
    A = sequence.count('A')
    T = sequence.count('T')
    C = sequence.count('C')
    G = sequence.count('G')
    all = A + T+ C+ G
    if all == 0 :
        print('len was 0')
    else:
        A=round(A/all  *100,2)
        T=round(T/all  *100,2)
        C=round(C/all  *100,2)
        G=round(G/all  *100 ,2)
        #print(all,"A",A,"T",T,"C",C,"G",G)
    return all,A,T,C,G

def gene_composition(sequence,path_gff,chromosome_number):
    gff_file=open(path_gff)
    all_genes_length = []
    all_A_percent = []
    all_T_percent = []
    all_C_percent = []
    all_G_percent = []
    for line in gff_file:
        if line[0]=="#":
            continue
        line=line.split()
        if line[2]=="gene":
            #     continue
            if chromosome_number  in line[0]:
                if line[2]=="gene":
                    # print(line)

                    start=int(line[3])
                    end=int(line[4])
                    gene_length=end-start
                    all_genes_length.append(gene_length)
                    sequence_gene=sequence[start:end]
                    # print(sequence_gene)
                    length,A,T,C,G = sequence_composition(sequence_gene)
                    all_A_percent.append(A)
                    all_T_percent.append(T)
                    all_C_percent.append(C)
                    all_G_percent.append(G)
    return all_genes_length,all_A_percent,all_T_percent,all_C_percent,all_G_percent


def list_to_plot(list):
    #list.sort()
    plt.plot(list)
    #plt.scatter(list)
    plt.show()

def list_to_plot2(L,A,T,C,G):
    #list.sort()
    #plt.plot(L,color = 'black')
    plt.plot(A,color= 'green')
    plt.plot(T,color = 'blue')
    plt.plot(C,color = 'yellow')
    plt.plot(G,color = 'red')

    #plt.scatter(list)
    plt.show()


def sort_legth_with_composition(length,A,T,C,G):
    L = length



    L2 = []
    A2 = []
    T2 = []
    C2 = []
    G2 = []


    for i in range(len(L)):
        max_number = max(L)
        L2.append(max_number)
        index_max_number = L.index(max_number)
        A2.append(A[index_max_number])
        T2.append(T[index_max_number])
        C2.append(C[index_max_number])
        G2.append(G[index_max_number])

        L.pop(index_max_number)
        A.pop(index_max_number)
        T.pop(index_max_number)
        C.pop(index_max_number)
        G.pop(index_max_number)

    print(L2)
    print(A2)
    print(T2)
    print(C2)
    print(G2)

    return L2,A2,T2,C2,G2

