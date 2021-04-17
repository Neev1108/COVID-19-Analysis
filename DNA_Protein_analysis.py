#!/usr/bin/env python
# coding: utf-8

# In[104]:


from collections import Counter
import random
import os
import subprocess
import sys
import csv
from Bio import SeqIO


# # Fasta class with similar methods to Biopython, will hold a sequence record

# In[105]:



class fasta:
    
    GC_content = 0
    c = {}
    transcription = ''
    translation = ''
    translation_nospaces = ''
    is_dna = True
    molecular_weight = float(0)
    coronavirus = False
    
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.c = Counter(sequence)
        self.checkDNA_orProtein()
        self.checkCoronavirus()
        
        if self.is_dna == True:
            self.transcription = self.transcribe()
            self.translation = self.translate()
        else:
            self.translation = self.sequence
            self.translation_nospaces = self.sequence
        
    
    def checkCoronavirus(self):
        if 'coronavirus' in self.name:
            self.coronavirus = True
        else:
            self.coronavirus = False
        
        
    def print_translation(self, sequence):
        if len(sequence) > 50:
            return sequence[0:50] + "...."
        else:
            return sequence
    
    def getLength(sequence):
        return len(sequence)
    
    def checkDNA_orProtein(self):
        DNA_set = {'A', 'C', 'G', 'T'}
        sequence_set = set(list(self.sequence))
        difference = sequence_set - DNA_set
        
        if len(difference) != 0:
            self.is_dna = False
        
    def getAccession(self):
        first_index = self.name.index('|')
        return self.name[0:first_index]
        
        
    def toString(self):
        return self.name + " " + self.sequence
    
    def get_GC_content(self):
        G_content = self.c.get('G')
        C_content = self.c.get('C')
        GC_content = (G_content + C_content)/len(self.sequence)
        return "{:.6%}".format(GC_content)
        
    def transcribe(self):
        dna_string = self.sequence
        dna_string = dna_string.replace('T', 'U')
        return dna_string
    
    def printDNACounter(self):
        return 'A = ' + str(self.c.get('A')), 'C = ' + str(self.c.get('C')), 'G = ' + str(self.c.get('G')), 'T = ' + str(self.c.get('T'))
    
    def getName(self):
        first_index = self.name.index('[')
        second_index = self.name.index(']')
        return self.name[first_index+1:second_index]
        
    
    def translate(self):
        rna = self.transcription
        codon_table = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
            "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
            "UAU":"Y", "UAC":"Y", "UAA":"Stop", "UAG":"Stop",
            "UGU":"C", "UGC":"C", "UGA":"Stop", "UGG":"W",
            "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
            "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
            "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
            "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
            "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
            "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
            "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
            "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
            "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
            
        amino_acid_chain = []
        for i in range(0,len(rna),3):
            codon = rna[i:i+3]
            amino_acid = codon_table.get(codon)
            amino_acid_chain.append(amino_acid)
        
        amino_acid_chain_string = '-'.join(str(x) for x in amino_acid_chain)
        self.translation_nospaces = amino_acid_chain_string.split('-')
        return amino_acid_chain_string
    
    def calculateMolecularWeight(self):
        aa = {'A' : 71.03711,'R' : 156.10111,'N' : 114.04293,'D' : 115.02694,'C' : 103.00919,'E' : 129.04259,
              'Q' : 128.05858,'G' : 57.02146,'H' : 137.05891,'I' : 113.08406,'L' : 113.08406,'K' : 128.09496,
              'M' : 131.04049,'F' : 147.06841,'P' : 97.05276,'S' : 87.03203,'T' : 101.04768,'W' : 186.07931,
              'Y' : 163.06333,'V' : 99.06841, 'X' : 0.00000  }
        protein_sequence = self.translation_nospaces
        molecular_weight = float(0)
        
        for amino_acid in protein_sequence:
            aa_weight = aa.get(amino_acid) 
            if aa_weight != None:
                molecular_weight += aa_weight
            else:
                continue
            
        self.molecular_weight = molecular_weight
        return molecular_weight
        
    
    def printDNA_analysis(self):
        print('\n Name:  ', self.getName())
        print('\n \t DNA Analysis: ')
        print('--------------------------')
        print('\t DNA sequence: ', self.sequence)
        print('\t', 'Transcription: ',self.transcription)
        print('\t', 'Number of each base in DNA: ', self.printDNACounter())
        print('\t GC content: ', self.get_GC_content())
    
            
    def printProtein_analysis(self):
        print('\n \t Protein Analysis: ')
        print('--------------------------')
        print('\n \t Amino Acid Chain: ', self.print_translation(self.translation))
        print('\t Number of Amino acids: ', len(self.translation_nospaces))
        print('\t Molecular Weight: ', self.calculateMolecularWeight(), 'Da')
        
    


# # Some Methods you can run on a list of sequences

# In[106]:


def profile_consensus_matrix(fasta_records):
    equal_length = len(fasta_records[0])
    count = 0
    profile_matrix = []
    
    for i in range(4):   
        row = []
        for column in zip(*fasta_records):
            if i == 0:
                row.append(column.count('A'))
                continue
            elif i == 1:
                row.append(column.count('C'))
                continue
            elif i == 2:
                row.append(column.count('G'))
                continue
            else:
                row.append(column.count('T'))
                continue
            
        profile_matrix.append(row)
    
    print("Profile Matrix\n")
    for index,value in enumerate(profile_matrix):
        if index == 0:
            print("A: ", end = " ")
        if index == 1:
            print("C: ", end = " ")
        if index == 2:
            print("G: ", end = " ")
        if index == 3:
            print("T: ", end = " ")
        for c in value:
            print(c, end = " ")
            
        print("\n")
        
    consensus_matrix = []
    for column in zip(*profile_matrix):
        max_num = max(column)
        index_of_max = column.index(max_num)
    
        if index_of_max == 0:
            consensus_matrix.append('A')
        if index_of_max == 1:
            consensus_matrix.append('C')
        if index_of_max == 2:
            consensus_matrix.append('G')
        if index_of_max == 3:
            consensus_matrix.append('T')
            
    print('\nConsensus Matrix \n')
    for x in consensus_matrix:
        print(x, end = " ")
        
        
def findingSharedMotif(fasta_records):
    fasta_records = sorted(fasta_records, key = len)
    shortest_seq = fasta_records[0]
    
    compare_seq = fasta_records[1:]
    motif = ''
    for x in range(len(shortest_seq)):
        for y in range(x, len(shortest_seq)):
            current_seq = shortest_seq[x:y+1]
            found = False
            for seq in compare_seq:
                if current_seq in seq:
                    found = True
                else:
                    found = False
                    break   
            if found == True and len(current_seq) > len(motif):
                motif = current_seq
                
    return motif

       


# # Printing the analysis and then some file IO.

# In[ ]:



f = open('basic_ds.fasta', 'r')
f2 = open('related_ds.fasta', 'r')

def file_into_fasta_records(f):
    fasta_records = []
    for seq_record in SeqIO.parse(f, "fasta"):
        name = seq_record.description
        sequence = str(seq_record.seq)
        fasta_records.append(fasta(name, sequence))
        
    print('\nNumber of records: ', len(fasta_records))
    return fasta_records

def fasta_records_analysis(fasta_records):
    sequence_record = []
    for x in fasta_records:
        sequence_record.append(x.sequence)
    
    print("Shared Motif is: ", findingSharedMotif(sequence_record))
    if x.is_dna == True:
        profile_consensus_matrix(fasta_records)
    for x in fasta_records:
        if x.is_dna == True:
            x.printDNA_analysis()
            x.printProtein_analysis()
        else:
            print('\n Name:  ', x.getName())
            x.printProtein_analysis()
            
def write_fasta(fasta_records, length, file, length_of_first):
    if length_of_first == True:
        for i in range(0,length):
            file.write('\n\n>' + fasta_records[i].name)
            file.write('\n' + fasta_records[i].sequence)
    else:
        for i in range(length, len(fasta_records)):
            file.write('\n\n>' + fasta_records[i].name)
            file.write('\n' + fasta_records[i].sequence)


# In[108]:


fasta_records_basic_ds = file_into_fasta_records(f)
print("\nBasic Dataset \n\n")
fasta_records_analysis(fasta_records_basic_ds)

fasta_records_related_ds = file_into_fasta_records(f2)
print("\nRelated Dataset \n\n")
fasta_records_analysis(fasta_records_related_ds)
                
        


# # Split the datasets into training and test splits. I manually had stratify the data to make sure it was a valid set

# In[110]:


eighty_percent_basic_ds = open('eighty_percent_basic_ds.fasta', 'w')
twenty_percent_basic_ds = open('twenty_percent_basic_ds.fasta', 'w')

try:
    os.mkdir("split_datasets_basic")
except FileExistsError: 
    print("File already exists, will continue")
    
os.rename('eighty_percent_basic_ds.fasta', "split_datasets_basic/eighty_percent_basic_ds.fasta")
os.rename('twenty_percent_basic_ds.fasta', "split_datasets_basic/twenty_percent_basic_ds.fasta")

eighty_percent_length = int(len(fasta_records_basic_ds) * .8)
twenty_percent_length = len(fasta_records_basic_ds) - eighty_percent_length

write_fasta(fasta_records_basic_ds, eighty_percent_length, eighty_percent_basic_ds, True)
write_fasta(fasta_records_basic_ds, len(fasta_records_basic_ds) - twenty_percent_length, twenty_percent_basic_ds, False)


eighty_percent_basic_ds.close()
twenty_percent_basic_ds.close()

eighty_percent_related_ds = open('eighty_percent_related_ds.fasta', 'w')
twenty_percent_related_ds = open('twenty_percent_related_ds.fasta', 'w')

try:
    os.mkdir("split_datasets_related")
except FileExistsError: 
    print("File already exists, will continue")
    
os.rename('eighty_percent_related_ds.fasta', "split_datasets_related/eighty_percent_basic_ds.fasta")
os.rename('twenty_percent_related_ds.fasta', "split_datasets_related/twenty_percent_basic_ds.fasta")

eighty_percent_length = int(len(fasta_records_related_ds) * .8)
twenty_percent_length = len(fasta_records_related_ds) - eighty_percent_length


write_fasta(fasta_records_related_ds, eighty_percent_length, eighty_percent_related_ds, True)
write_fasta(fasta_records_related_ds, len(fasta_records_related_ds) - twenty_percent_length, twenty_percent_related_ds, False)


eighty_percent_related_ds.close()
twenty_percent_related_ds.close()


