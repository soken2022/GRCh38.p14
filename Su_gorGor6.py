import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gzip
from Bio import SeqIO

# input information

input_organism_name = 'groGro6'
input_element_name = 'SINE'  # ['LINE','SINE','LTR','DNA'] #study genomic elements
input_element_name_for_saving = input_element_name
input_start_chromosome_number = 1
input_end_chromosome_number = 24
input_number_of_chromosomes = input_end_chromosome_number - input_start_chromosome_number + 1
input_chromosomes_name = ['1', '2A', '2B', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
                          '17', '18', '19', '20', '21', '22', '23', '24']
input_plot_len = 2000
annotation_file_path = 'D:\PhD\genomes\gorilla\gorGor6.fa.out.gz'
fasta_path = 'D:\PhD\genomes\gorilla\gorGor6.fa'

CC = input_chromosomes_name
CC = [input('ch_number?')]

sequence = ''
fasta_path = fasta_path
ch_number = '22'
print(ch_number)
for seq_record in SeqIO.parse(fasta_path, "fasta"):
    if seq_record.id == rf'chr{ch_number}':
        print(seq_record.id)
        sequence = str(seq_record.seq)
        break
len_chromosoeme_su  =len(sequence)
print(len_chromosoeme_su)





# class definition

class GroGro6:
    def __init__(self, organism_name, chromosome_sequence, number_of_chromosome, element_name, chr_annotation,len_chromosoeme_su):
        self.organism_name = organism_name
        self.chromosome_sequence = chromosome_sequence
        self.number_of_chromosome = number_of_chromosome
        self.element_name = element_name
        self.chr_annotation = chr_annotation
        self.len_chromosoeme_su = len_chromosoeme_su

    def open_fasta_file_and_prepare_it(self):
        fasta_path = self.chromosome_sequence
        ch_number = self.number_of_chromosome
        fasta_file = open(fasta_path)
        fasta_file.readline()
        sequence = fasta_file.read()
        sequence = sequence.replace('\n', '')
        fasta_file.close()
        print(ch_number)
        len_chromosome = len(sequence)
        return sequence , len_chromosome


    def open_fasta_file_multi_record(self):
        sequence = ''
        fasta_path = self.chromosome_sequence
        ch_number = self.number_of_chromosome
        print(ch_number)
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            if seq_record.id == rf'chr{ch_number}':
                print(seq_record.id)
                sequence = str(seq_record.seq)
                break

        return sequence




    def open_annotation_file_and_prepare_it(self):
        chr_annotation = self.chr_annotation
        ch_number = self.number_of_chromosome
        file_gff = gzip.open(chr_annotation,'rt')
        print('mapping start', ch_number)
        return file_gff

    def counting_of_elements(self, sequence, gff, plot_length):
        counting_of_plus_elements = 0
        counting_of_comp_elements = 0
        sequences_of_plus_elements = ''
        sequences_of_complement_elements = ''
        elements_overlap = list(sequence)
        counting_pluses_in_difference_length = [0] *plot_length
        counting_complements_in_difference_length = [0] * plot_length
        compposintion_in_diffrence_length_plus = [''] * plot_length
        compposintion_in_diffrence_length_comp = [''] * plot_length
        # element_location_d = [0] * len(sequence)
        # Alus_location_c = [0] * len(sequence)

        '''A_portion_in_length = [0] * 4000
        T_portion_in_length = [0] * 4000
        C_portion_in_length = [0] * 4000
        G_portion_in_length = [0] * 4000'''

        chromosome_number = self.number_of_chromosome
        element_name = self.element_name
        if self.number_of_chromosome == 23:
            chromosome_number = 'X'
        if self.number_of_chromosome == 24:
            chromosome_number = 'Y'
        gff.readline()
        gff.readline()
        gff.readline()

        sequence2_for_overlabing = sequence
        for annotation in gff:
            annotation = annotation.split()
            if (element_name in annotation[10] or 'LINE' in annotation[10] or 'LTR' in annotation[10] or 'DNA' in annotation[10]) and (annotation[4] == f'chr{chromosome_number}'): # or annotation[4] == f'chr{chromosome_number}_'):
                if annotation[8] == '+':
                    sequence_of_element_p = sequence[int(annotation[5]) - 1:int(annotation[6]) - 1]

                    sequences_of_plus_elements += sequence_of_element_p
                    elements_overlap[int(annotation[5]) - 1:int(annotation[6]) - 1] = sequence_of_element_p.lower()
                    #element_location_d[int(line[5])] += 10
                    len_element_p = len(sequence_of_element_p)


                    if len_element_p == 0 :
                        #print(annotation)
                        print('start',int(annotation[5]))
                        print('len',len(sequence))


                    if len_element_p > plot_length:
                        print('Not studied Element it is very big',len_element_p)
                    if len_element_p < plot_length:
                        counting_of_plus_elements += 1
                        counting_pluses_in_difference_length[len_element_p] += 1
                        compposintion_in_diffrence_length_plus[len_element_p] += sequence_of_element_p


                    '''A_portion_in_length[len_element] += Eplus.count('A')
                    T_portion_in_length[len_element] += Eplus.count('T')
                    C_portion_in_length[len_element] += Eplus.count('C')
                    G_portion_in_length[len_element] += Eplus.count('G')'''
                elif 'C' in annotation:
                    sequence_of_element_c = sequence[int(annotation[5]) - 1:int(annotation[6]) - 1]
                    sequences_of_complement_elements += sequence_of_element_c
                    elements_overlap[int(annotation[5]) - 1:int(annotation[6]) - 1] = sequence_of_element_c.lower()
                    #Alus_location_c[int(line[5])] += 10
                    len_element_c = len(sequence_of_element_c)
                    if len_element_c > plot_length:
                        print('Not studied Element it is very big',len_element_c)
                    if len_element_c < plot_length:
                        counting_of_comp_elements += 1
                        counting_complements_in_difference_length[len_element_c] += 1
                        compposintion_in_diffrence_length_comp[len_element_c]+=sequence_of_element_c

        print('element counting finished')
        print('  plus    comp')
        print(counting_of_plus_elements,'  ',counting_of_comp_elements)
        plus_elements_portion_in_sequence = round(len(sequences_of_plus_elements)/len(sequence)*100,2)
        complement_elements_portion_in_sequence = round(len(sequences_of_complement_elements)/len(sequence)*100,2)
        print(' portion of elements sequence in + and - strand:\n',plus_elements_portion_in_sequence,'%',complement_elements_portion_in_sequence,'%')
        len_sum_of_pluses = len(sequences_of_plus_elements)
        A_d = round(sequences_of_plus_elements.count('a') / len_sum_of_pluses * 100, 2)
        T_d = round(sequences_of_plus_elements.count('t') / len_sum_of_pluses * 100, 2)
        C_d = round(sequences_of_plus_elements.count('c') / len_sum_of_pluses * 100, 2)
        G_d = round(sequences_of_plus_elements.count('g') / len_sum_of_pluses * 100, 2)
        len_sum_of_complements = len(sequences_of_complement_elements)
        A_c = round(sequences_of_complement_elements.count('a') / len_sum_of_complements * 100, 2)
        T_c = round(sequences_of_complement_elements.count('t') / len_sum_of_complements * 100, 2)
        C_c = round(sequences_of_complement_elements.count('c') / len_sum_of_complements * 100, 2)
        G_c = round(sequences_of_complement_elements.count('g') / len_sum_of_complements * 100, 2)
        print('+composition', '-composition')
        print('A ', A_d, '   ', ' A ', A_c)
        print('T ', T_d, '   ', ' T ', T_c)
        print('C ', C_d, '   ', ' C ', C_c)
        print('G ', G_d, '   ', ' G ', G_c)

        plus_stand_information = [counting_of_plus_elements, counting_pluses_in_difference_length, sequences_of_plus_elements, len(sequence),compposintion_in_diffrence_length_plus]
        comp_strand_information = [counting_of_comp_elements,  counting_complements_in_difference_length, sequences_of_complement_elements, len(sequence),compposintion_in_diffrence_length_comp]


        '''sum_of_all_ch_information_result_list_for_plus = []
        all_ch_p_result_list[0] += counting_of_plus_elements
        all_ch_p_result_list[1] += len(sequences_of_plus_elements)
        all_ch_p_result_list[2] += sequences_of_plus_elements.count('A')
        all_ch_p_result_list[3] += sequences_of_plus_elements.count('T')
        all_ch_p_result_list[4] += sequences_of_plus_elements.count('C')
        all_ch_p_result_list[5] += sequences_of_plus_elements.count('G')
        global all_ch_c_result_list
        all_ch_c_result_list[0] += counting_of_comp_elements
        all_ch_c_result_list[1] += len(sequences_of_complement_elements)
        all_ch_c_result_list[2] += sequences_of_complement_elements.count('A')
        all_ch_c_result_list[3] += sequences_of_complement_elements.count('T')
        all_ch_c_result_list[4] += sequences_of_complement_elements.count('C')
        all_ch_c_result_list[5] += sequences_of_complement_elements.count('G')

        all_ch_Es_in_diff_len_plus
        all_ch_Es_in_diff_len_comp

        all_ch_Es_in_diff_len_plus = [x + y for x, y in zip(all_ch_Es_in_diff_len_plus, counting_pluses_in_difference_length)]
        all_ch_Es_in_diff_len_comp = [x + y for x, y in zip(all_ch_Es_in_diff_len_comp, counting_complements_in_difference_length)]

        #count_over laping'''

        A_o = elements_overlap.count('a')
        T_o = elements_overlap.count('t')
        C_o = elements_overlap.count('c')
        G_o = elements_overlap.count('g')
        length_of_overlaping = A_o + T_o + C_o + G_o
        print('overlaping information',length_of_overlaping,length_of_overlaping/len(sequence)*100,A_o,T_o,C_o,G_o,)
        print('overlap to priyor',length_of_overlaping/(len(sequences_of_plus_elements)+len(sequences_of_complement_elements))*100)
        overlap_information = []
        return plus_stand_information, comp_strand_information, length_of_overlaping

    def scatering_plot_for_frequency_in_differrence_length(self,plus_frequencey,comp_frequencey):
        plt.figure(figsize=(13, 4))
        X1 = np.arange(len(plus_frequencey))
        X2 = np.arange(len(comp_frequencey))

        plt.scatter(X1,plus_frequencey,color= 'green',s = 10)
        plt.scatter(X2,comp_frequencey,color = 'red', s = 10)

        diffrence_length_portion_in_ch_plus = [0]*len(X1)
        for i in range(len(plus_frequencey)):
            diffrence_length_portion_in_ch_plus[i] = plus_frequencey[i] * X1[i]/180

        diffrence_length_portion_in_ch_comp = [0]*len(X1)
        for j in range(len(comp_frequencey)):
            diffrence_length_portion_in_ch_comp[j] = comp_frequencey[j] * X1[j]/180
        plt.scatter(X1,diffrence_length_portion_in_ch_plus,color = 'black', s = 10)
        plt.scatter(X1,diffrence_length_portion_in_ch_comp,color = 'yellow', s = 10)
        Z = self.element_name




        plt.xlabel(Z)
        plt.ylabel('frequency')
        plt.show()

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#objects
#QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQq




for i in CC:
    chromosome_number = i
    Human = GroGro6(input_organism_name, fasta_path, chromosome_number,input_element_name,annotation_file_path,len_chromosoeme_su)
    element_name = Human.element_name
    sequence = Human.open_fasta_file_multi_record()
    mask_map = Human.open_annotation_file_and_prepare_it()
    plus_element_inf, comp_element_inf, overlaping_inf = Human.counting_of_elements(sequence,mask_map,input_plot_len)
    print(plus_element_inf[1])


Human.scatering_plot_for_frequency_in_differrence_length(plus_element_inf[1],comp_element_inf[1])

def scatering_plot_for_Nucleotid_composition_in_differrence_length(plus_element_inf_4,comp_element_inf_4):

    Ap = [round(k.count('a')/(len(k)+0.001)*100,2) for k in plus_element_inf_4]
    Tp = [round(k.count('t')/(len(k)+0.001)*100,2) for k in plus_element_inf_4]
    Cp = [round(k.count('c')/(len(k)+0.001)*100,2) for k in plus_element_inf_4]
    Gp = [round(k.count('g')/(len(k)+0.001)*100,2) for k in plus_element_inf_4]

    Ac = [round(h.count('a')/(len(h)+0.001)*100,2) for h in comp_element_inf_4]
    Tc = [round(h.count('t')/(len(h)+0.001)*100,2) for h in comp_element_inf_4]
    Cc = [round(h.count('c')/(len(h)+0.001)*100,2) for h in comp_element_inf_4]
    Gc = [round(h.count('g')/(len(h)+0.001)*100,2) for h in comp_element_inf_4]
    Ap_Ac = [(Ap[i]+ Ac[i])/2 for i in range(len(Ap))]
    Tp_Tc = [(Tp[i]+ Tc[i])/2 for i in range(len(Ap))]
    Cp_Cc = [(Cp[i]+ Cc[i])/2 for i in range(len(Ap))]
    Gp_Gc = [(Gp[i]+ Gc[i])/2 for i in range(len(Ap))]
    import matplotlib.patches as mpatches

    fig, ax = plt.subplots(figsize=(13, 4))
    X1 = np.arange(len(plus_element_inf_4))
    plt.scatter(X1,Ap,color= 'green',s = 10)
    plt.scatter(X1,Tp,color= 'blue',s = 10)
    plt.scatter(X1,Ac,color= 'yellow',s = 10)
    plt.scatter(X1,Tc,color= 'red',s = 10)

    plt.xlabel('Element length')
    plt.ylabel('Nucleotide composition')
    green_patch = mpatches.Patch(color='green', label='A composition in plus E')
    blue_patch = mpatches.Patch(color='blue', label='T composition in plus E')
    yellow_patch = mpatches.Patch(color='yellow', label='A composition in comp E')
    red_patch = mpatches.Patch(color='red', label='T composition in comp E')
    ax.legend(handles=[green_patch,blue_patch,yellow_patch,red_patch])
    #ax.legend(handles=[blue_patch])

    plt.show()


    fig, ax2 = plt.subplots(figsize=(13, 4))
    X1 = np.arange(len(plus_element_inf_4))
    plt.scatter(X1,Ap_Ac,color= 'green',s = 10)
    plt.scatter(X1,Tp_Tc,color= 'blue',s = 10)
    plt.scatter(X1,Cp_Cc,color= 'yellow',s = 10)
    plt.scatter(X1,Gp_Gc,color= 'red',s = 10)

    plt.ylabel('Nucleotide composition')
    green_patch = mpatches.Patch(color='green', label='Ap_Ac')
    blue_patch = mpatches.Patch(color='blue', label='Tp_Tc')
    yellow_patch = mpatches.Patch(color='yellow', label='Cp_Cc')
    red_patch = mpatches.Patch(color='red', label='Gp_Gc')
    ax2.legend(handles=[green_patch,blue_patch,yellow_patch,red_patch])
    #ax.legend(handles=[blue_patch])

    plt.show()


scatering_plot_for_Nucleotid_composition_in_differrence_length(plus_element_inf[4],comp_element_inf[4])
#plt.plot(C,color = 'yellow')
#plt.plot(G,color = 'red')

