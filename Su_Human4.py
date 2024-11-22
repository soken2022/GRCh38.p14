
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
input_element_name = ['LINE','SINE','LTR','DNA'] #study genomic elements
input_element_name_for_saving = '_'.join(input_element_name)
input_start_chromosome_number = 1
input_end_chromosome_number = 24
input_number_of_chromosomes = input_end_chromosome_number - input_start_chromosome_number + 1
input_chromosomes_name = [22]#['1','2','3','4','5','6','7','8','9' ,'10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22']
input_plot_len = 2000
annotation_file_path_out_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out.gz'
fasta_path_refseq = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
fasta_path_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.gz'


#first informotation prepration
'''def first_prepration(fasta_path,ch_num):
    sequence = ''
    ch_number = str(ch_num)
    print('input ch_number',ch_number)
    if fasta_path[-2:] == 'gz':
        fasta_path = gzip.open(fasta_path,'rt')
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        if seq_record.id == f'chr{ch_number}':
            print('seq_record.id',seq_record.id)
            sequence = str(seq_record.seq)
            break
    len_chromosoeme  =len(sequence)
    print('len_chromosoeme',len_chromosoeme)
    return len_chromosoeme,ch_number
len_chromosoeme,ch_number = first_prepration(fasta_path_UCSC,22)'''



# class definition **************************************************************************
class GRCh38_p14:

    def __init__(self, organism_name, genome_path, annotation_path, number_of_chromosome, len_chromosoeme, element_name):
        self.organism_name = organism_name
        self.chromosome_sequence = genome_path
        self.number_of_chromosome = number_of_chromosome
        self.element_name = element_name
        self.chr_annotation = annotation_path
        self.len_chromosoeme_su = len_chromosoeme

    def open_simple_fasta_file_and_prepare_it(self):
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

    def open_fasta_file_UCSC_multi_record(self):
        sequence = ''
        fasta_path = self.chromosome_sequence
        ch_number = self.number_of_chromosome
        if fasta_path[-2:] == 'gz':
            fasta_path = gzip.open(fasta_path,'rt')
        print('for chromosoeme',ch_number)
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            if seq_record.id == rf'chr{ch_number}':
                print('sequence id is',seq_record.id)
                sequence = str(seq_record.seq)
                fasta_path.close()
                l  = len(sequence)
                print('len chromosome is',l,'bp')

                break

        return sequence

    def open_fasta_file_multi_record_from_refseq(self):
        sequence = ''
        fasta_path = self.chromosome_sequence
        ch_number = self.number_of_chromosome
        print('ch_ref',ch_number)
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            sequence_id = seq_record.id
            if f'NC_0000{ch_number}' in sequence_id :
                with open(r'D:\PhD\genomes\GRCh38.p14\hgdownload.soe.ucsc.edu_goldenPath_hg38_bigZips_p14_hg38.p14.chromAlias.txt') as filex :
                    for line in filex:
                        if sequence_id in line:
                            split_line = line.split()
                            sequence = str(seq_record.seq)

                            break

        return sequence

    def open_annotation_file_and_prepare_it(self):
        chr_annotation = self.chr_annotation
        ch_number = self.number_of_chromosome
        if chr_annotation[-2:] == 'gz':
            file_gff = gzip.open(chr_annotation,'rt')
        else:
            file_gff = open(chr_annotation)
        print('mapping start', ch_number)
        return file_gff

    def element_information_record(self,plus_negative,element_group,sequence_of_element,plot_length):
        # 1 element number 2 sequence 3 length distribution 4 composition in diffrence length


        len_element = len(sequence_of_element)
        if len_element == 0 :
            #print(annotation)
            print('len is 0')
        if len_element > plot_length:
            print('Not studied Element it is very big',len_element)
        if len_element < plot_length:
            #element_information[plus_negative][element_group][0] += 1
            #element_information[plus_negative][element_group][1] += sequence_of_element
            element_information[plus_negative][element_group][2][len_element] += 1
            element_information[plus_negative][element_group][3][len_element] += sequence_of_element
        return element_information

    def counting_of_elements(self, sequence, gff, plot_length):
        chromosome_number = self.number_of_chromosome
        element_name = self.element_name

        # 1 element number 2 sequence 3 length distribution 4 composition in diffrence length
        global element_information
        element_information = [[[[0]*plot_length,['']*plot_length]]*4,[[[0]*plot_length,['']*plot_length]]*4]


        if self.number_of_chromosome == 23:
            chromosome_number = 'X'
        if self.number_of_chromosome == 24:
            chromosome_number = 'Y'
        gff.readline()
        gff.readline()
        gff.readline()

        #sequence2_for_overlabing = sequence
        for annotation in gff:
            annotation = annotation.split()
            sequence_of_element = sequence[int(annotation[5]) - 1:int(annotation[6]) - 1]

            if element_name[0] in annotation[10] and (annotation[4] == f'chr{chromosome_number}'):
                #print('SINE')
                if annotation[8] == '+':
                    Human21.element_information_record(0,0,sequence_of_element,plot_length)
                else:
                    Human21.element_information_record(1,0,sequence_of_element,plot_length)
            if element_name[1] in annotation[10]  and (annotation[4] == f'chr{chromosome_number}'):
                #print('LINE')
                if annotation[8] == '+':
                    Human21.element_information_record(0,1,sequence_of_element,plot_length)
                else:
                    Human21.element_information_record(1,1,sequence_of_element,plot_length)
            if element_name[2] in annotation[10]  and (annotation[4] == f'chr{chromosome_number}'):
                if annotation[8] == '+':
                    Human21.element_information_record(0,2,sequence_of_element,plot_length)
                else:
                    Human21.element_information_record(1,2,sequence_of_element,plot_length)
            if element_name[3] in annotation[10] and (annotation[4] == f'chr{chromosome_number}'): # or annotation[4] == f'chr{chromosome_number}_'):
                if annotation[8] == '+':
                    Human21.element_information_record(0,3,sequence_of_element,plot_length)
                else:
                    Human21.element_information_record(1,3,sequence_of_element,plot_length)
        print(element_information)
        return element_information




        #both making#########################################################
        counting_of_both_elements = counting_of_plus_elements + counting_of_comp_elements
        counting_both_in_difference_length = [sum(i) for i in zip(counting_pluses_in_difference_length,counting_complements_in_difference_length)]
        sequences_of_both_elements = sequences_of_plus_elements + sequences_of_complement_elements
        compposintion_in_diffrence_length_both = [i+j for i,j in zip(compposintion_in_diffrence_length_plus,compposintion_in_diffrence_length_comp)]

        print('element counting finished')
        print('both   plus    comp')
        print(counting_of_both_elements,'  ',counting_of_plus_elements,'  ',counting_of_comp_elements)
        plus_elements_portion_in_sequence = round(len(sequences_of_plus_elements)/len(sequence)*100,2)
        complement_elements_portion_in_sequence = round(len(sequences_of_complement_elements)/len(sequence)*100,2)
        both_elements_portion_in_sequence = round(len(sequences_of_both_elements)/len(sequence)*100,2)
        print(' portion of elements in sequence :\nboth plus comp\n',both_elements_portion_in_sequence,plus_elements_portion_in_sequence,'%',complement_elements_portion_in_sequence,'%')
        len_both_seq = len(sequences_of_both_elements)
        A_b = round(sequences_of_both_elements.count('a') / len_both_seq * 100, 2)
        T_b = round(sequences_of_both_elements.count('t') / len_both_seq * 100, 2)
        C_b = round(sequences_of_both_elements.count('c') / len_both_seq * 100, 2)
        G_b = round(sequences_of_both_elements.count('g') / len_both_seq * 100, 2)
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
        print('both composition','+composition', '-composition')
        print('A ',A_b,'    ','A ', A_d, '   ', ' A ', A_c)
        print('T ',T_b,'    ','T ', T_d, '   ', ' T ', T_c)
        print('C ',C_b,'    ','C ', C_d, '   ', ' C ', C_c)
        print('G ',G_b,'    ','G ', G_d, '   ', ' G ', G_c)
        both_stand_information = [counting_of_both_elements,counting_both_in_difference_length,sequences_of_both_elements,len(sequence),compposintion_in_diffrence_length_both]
        plus_stand_information = [counting_of_plus_elements, counting_pluses_in_difference_length, sequences_of_plus_elements, len(sequence),compposintion_in_diffrence_length_plus]
        comp_strand_information = [counting_of_comp_elements,  counting_complements_in_difference_length, sequences_of_complement_elements, len(sequence),compposintion_in_diffrence_length_comp]
        del sequences_of_plus_elements
        del sequences_of_complement_elements
        del sequences_of_both_elements
        del counting_of_plus_elements
        del counting_of_comp_elements
        del counting_of_both_elements
        del counting_pluses_in_difference_length
        del counting_complements_in_difference_length
        del counting_both_in_difference_length
        del sequence
        del sequence_of_element_p
        del sequence_of_element_c
        gff.close()
        return both_stand_information,plus_stand_information, comp_strand_information, #length_of_overlaping

    def segments_to_composition(self,both_element_inf_4,plus_element_inf_4,comp_element_inf_4):
        Ab = [k.count('a') for k in both_element_inf_4]
        Tb = [k.count('t') for k in both_element_inf_4]
        Cb = [k.count('c') for k in both_element_inf_4]
        Gb = [k.count('g') for k in both_element_inf_4]
        B_group = [Ab,Tb,Cb,Gb]

        Ap = [k.count('a') for k in plus_element_inf_4]
        Tp = [k.count('t') for k in plus_element_inf_4]
        Cp = [k.count('c') for k in plus_element_inf_4]
        Gp = [k.count('g') for k in plus_element_inf_4]
        P_group = [Ap,Tp,Cp,Gp]
        Ac = [h.count('a') for h in comp_element_inf_4]
        Tc = [h.count('t') for h in comp_element_inf_4]
        Cc = [h.count('c') for h in comp_element_inf_4]
        Gc = [h.count('g') for h in comp_element_inf_4]
        C_group = [Ac,Tc,Cc,Gc]
        return B_group,P_group,C_group

    def scatering_plot_for_frequency_in_differrence_length(self,Both_frequencey,plus_frequencey,comp_frequencey):

        X1 = np.arange(len(plus_frequencey))
        X2 = np.arange(len(comp_frequencey))
        '''plt.plot(X1,sum_of_plus_and_comp,color = 'black')
        #plt.show()
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
        #plt.show()'''
        fig,axs = plt.subplots(2)
        axs[0].plot(X1,Both_frequencey, color = 'black')
        axs[1].plot(X1,plus_frequencey,color= 'green')
        axs[1].plot(X2,comp_frequencey,color = 'red')
        plt.show()

    def scatering_plot_for_Nucleotid_composition_in_differrence_length(self,B_group,P_group,C_group):

        B_ATCG = [sum(i) for i in zip(B_group[0],B_group[1],B_group[2],B_group[3])]
        Ab = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[0],start=0)]
        Tb = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[1],start=0)]
        Cb = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[2],start=0)]
        Gb = [round(v/(B_ATCG[i]+0.001)*100,2) for i,v in enumerate(B_group[3],start=0)]


        P_ATCG = [sum(i) for i in zip(P_group[0],P_group[1],P_group[2],P_group[3])]
        Ap = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[0],start=0)]
        Tp = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[1],start=0)]
        Cp = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[2],start=0)]
        Gp = [round(v/(P_ATCG[i]+0.001)*100,2) for i,v in enumerate(P_group[3],start=0)]

        C_ATCG = [sum(i) for i in zip(C_group[0],C_group[1],C_group[2],C_group[3])]
        Ac = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[0],start=0)]
        Tc = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[1],start=0)]
        Cc = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[2],start=0)]
        Gc = [round(v/(C_ATCG[i]+0.001)*100,2) for i,v in enumerate(C_group[3],start=0)]
        Ap_Ac = [(Ap[i]+ Ac[i])/2 for i in range(len(Ap))]
        Tp_Tc = [(Tp[i]+ Tc[i])/2 for i in range(len(Ap))]
        Cp_Cc = [(Cp[i]+ Cc[i])/2 for i in range(len(Ap))]
        Gp_Gc = [(Gp[i]+ Gc[i])/2 for i in range(len(Ap))]
        import matplotlib.patches as mpatches

        plt.figure(figsize=(13, 4))
        X1 = np.arange(len(Ap))
        plt.scatter(X1,Ab,color= 'green',s = 10)
        plt.scatter(X1,Tb,color= 'blue',s = 10)
        plt.scatter(X1,Cb,color= 'yellow',s = 10)
        plt.scatter(X1,Gb,color= 'red',s = 10)
        plt.title('number 3 plot')
        plt.show()


        fig, axs = plt.subplots(3,figsize=(13, 4))
        X1 = np.arange(len(Ap))
        axs[0].scatter(X1,Ap,color= 'green',s = 10)
        axs[0].scatter(X1,Tp,color= 'blue',s = 10)
        axs[0].scatter(X1,Cp,color= 'lime',s = 10)
        axs[0].scatter(X1,Gp,color= 'cyan',s = 10)
        #axs[0].title('plot number 4_P')
        axs[1].scatter(X1,Ac,color= 'yellow',s = 10)
        axs[1].scatter(X1,Tc,color= 'red',s = 10)
        axs[1].scatter(X1,Cc,color= 'deeppink',s = 10)
        axs[1].scatter(X1,Gc,color= 'navy',s = 10)

        axs[2].scatter(X1,Ap,color= 'green',s = 10)
        axs[2].scatter(X1,Tp,color= 'blue',s = 10)
        axs[2].scatter(X1,Cp,color= 'lime',s = 10)
        axs[2].scatter(X1,Gp,color= 'cyan',s = 10)
        #axs[0].title('plot number 4_P')
        axs[2].scatter(X1,Ac,color= 'yellow',s = 10)
        axs[2].scatter(X1,Tc,color= 'red',s = 10)
        axs[2].scatter(X1,Cc,color= 'deeppink',s = 10)
        axs[2].scatter(X1,Gc,color= 'navy',s = 10)



        #axs[1].title('plot number 4_C')
        #plt.xlabel('Elements length')
        #plt.ylabel('Nucleotide composition %')


        '''green_patch = mpatches.Patch(color='green', label='A composition in plus E')
        blue_patch = mpatches.Patch(color='blue', label='T composition in plus E')
        lime_patch = mpatches.Patch(color='lime', label='C composition in plus E')
        cyan_patch = mpatches.Patch(color='cyan', label='G composition in plus E')
        yellow_patch = mpatches.Patch(color='yellow', label='A composition in comp E')
        red_patch = mpatches.Patch(color='red', label='T composition in comp E')
        deeppink_patch = mpatches.Patch(color='deeppink', label='C composition in comp E')
        navy_patch = mpatches.Patch(color='navy', label='G composition in comp E')
        ax.legend(handles=[green_patch,blue_patch,lime_patch,cyan_patch,yellow_patch,red_patch,deeppink_patch,navy_patch])
        #ax.legend(handles=[blue_patch])

        plt.show()'''



        '''X1 = np.arange(len(Ap))
        axs[2].scatter(X1,Ap_Ac,color= 'green',s = 10)
        axs[2].scatter(X1,Tp_Tc,color= 'blue',s = 10)
        axs[2].scatter(X1,Cp_Cc,color= 'yellow',s = 10)
        axs[2].scatter(X1,Gp_Gc,color= 'red',s = 10)'''

        '''plt.ylabel('Nucleotide composition')
        green_patch = mpatches.Patch(color='green', label='Ap_Ac')
        blue_patch = mpatches.Patch(color='blue', label='Tp_Tc')
        yellow_patch = mpatches.Patch(color='yellow', label='Cp_Cc')
        red_patch = mpatches.Patch(color='red', label='Gp_Gc')
        ax2.legend(handles=[green_patch,blue_patch,yellow_patch,red_patch])
        #ax.legend(handles=[blue_patch])'''

        plt.show()


#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#objects
#QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQq
all_Analysis_for_plus = [0] * input_plot_len
all_Analysis_for_comp = [0] * input_plot_len
all_Analysis_for_both = [0] * input_plot_len

all_B_group = [[0]*input_plot_len] * 4
all_P_group = [[0]*input_plot_len] * 4
all_C_group = [[0]*input_plot_len] * 4

for ch in input_chromosomes_name:
    Human21 = GRCh38_p14(input_organism_name,fasta_path_UCSC,annotation_file_path_out_UCSC,ch,20,input_element_name)
    print(Human21.organism_name)
    print(Human21.chr_annotation)
    print(Human21.element_name)
    seq = Human21.open_fasta_file_UCSC_multi_record()
    E_mask = Human21.open_annotation_file_and_prepare_it()
    Human21.counting_of_elements(seq,E_mask,input_plot_len)
    #Analysis_for_both,Analysis_for_plus,Analysis_for_comp = Human21.counting_of_elements(seq,E_mask,input_plot_len)
    all_Analysis_for_plus = [sum(i) for i in zip(all_Analysis_for_plus,Analysis_for_plus[1])]
    all_Analysis_for_comp = [sum(i) for i in zip(all_Analysis_for_comp,Analysis_for_comp[1])]
    all_Analysis_for_both = [sum(i) for i in zip (all_Analysis_for_both,Analysis_for_both[1])]
    B_group,P_group,C_group = Human21.segments_to_composition(Analysis_for_both[4],Analysis_for_plus[4],Analysis_for_comp[4])
    for h in range(4):
        all_B_group[h] = [sum(i) for i in zip(all_B_group[h],B_group[h])]
        all_P_group[h] = [sum(i) for i in zip(all_P_group[h],P_group[h])]
        all_C_group[h] = [sum(i) for i in zip(all_C_group[h],C_group[h])]
    del P_group
    del C_group
    del B_group
    del seq
    del E_mask
    del Analysis_for_plus
    del Analysis_for_comp
    del Analysis_for_both
Human21.scatering_plot_for_frequency_in_differrence_length(all_Analysis_for_both,all_Analysis_for_plus,all_Analysis_for_comp)
Human21.scatering_plot_for_Nucleotid_composition_in_differrence_length(all_B_group,all_P_group,all_C_group)




