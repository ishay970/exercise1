#Ishay Eldar

import sys

#  1.1:
def find_ssr(dna_sec):
    """  The function receives an input of a genetic sequence
         and returns a dictionary of the number of repetitions
         of each segment up to 6 nucleotides long
    :param dna_sec: a given sequence
    :return: dict: A dictionary that contains as a key is all sequences of length 1-6 that repeat at least 3 times,
                    and as a value - number of repetitions
    """
    dict = {}                      #The dictionary that the function will return
    counter = 1                    #number of reapits
    string = dna_sec               #the secuens that we check whether it contains SSR

    for n in range(1, 7) :
        for index in range(len(string) - n):
            counter = 1
            ssr_length = string[index: index+n]                #the SSR that we checking 
            if string[index+n: index + (2*n)] == ssr_length:   
                for k in range(len(string)//n) :
                    if string[index+n + k*n: index+n +(k+1)*n] == ssr_length: #We found an SSR with 3 repeats
                        counter += 1
                    else:
                        break
                if counter >= 3 :                    #We will add the SSR to the dictionary
                    if dict.get(ssr_length,0) == 0 :
                        dict[ssr_length] = counter
                    elif counter > dict.get(ssr_length):
                        dict[ssr_length] = counter
    if not any(dict):                                #We did not find any SSR
        return None
   
    return (dict)

"""for check:
"AAACAAAA" output: 'a' : 4
""ATCAAATCAAATCAAGAGAGAGGGGG" output: 'A': 3, 'G': 5, 'AG': 4, 'GA': 3, 'ATCAA': 3
"""

def dict_to_print(dictionary):
    """
    The function passes the dictionary to a sorted list that adds ";" between the values
    :param dictionary: dict of SSRS and repeats
    :return: linked list
    """
    dict_to_sort = {}
    list_sorted_keys = sorted(dictionary.keys())
    for key in list_sorted_keys:
        dict_to_sort[key] = dictionary[key]
    str_to_print = ';'.join([f"{k},{v}" for k, v in dict_to_sort.items()])
    return str_to_print

def ssr_printer(dna_seq1):
    """The function checks whether the dictionary is empty or the input is empty.
    Otherwise prints the dictionary values as a sorted list
    :param dna_seq1: the input of sequens
    """
    if len(dna_seq1) == 0:
        print("No simple repeats in DNA sequence")
    if  find_ssr(dna_seq1) is None:
        print("No simple repeats in DNA sequence")
    else:
        print(dict_to_print(find_ssr(dna_seq1)))

#---------------------------------------------------------------------------


#1.2:
    """The function receives an input of a DNA sequence,
    and returns the RNA that will be received in the transcription process
    :param dna_seq: a given DNA sequence that we needs to transcribe to RNA
    return: string of RNA sequense
    """
        
def transcribe(dna_seq) :
    dna = dna_seq
    rna_temporary = dna[len(dna)-1::-1].upper().replace('A','U').replace('T', 'A').replace('C','g').replace('G', 'C')
    rna = rna_temporary.replace('g', 'G')
    return rna
    
"""for check:
    dna = "CttGAT"
    print(rna)
    output: AUCAAG"""
"""
    The function prints the RNA sequence obtained from the input DNA sequence
    """
def transcribe_printer(dna_seq2):
    rna_to_print = transcribe(dna_seq2)
    print("RNA sequence:", rna_to_print)
    
#--------------------------------------------------------------------------

#1.3:
    """
    The function accepts an RNA sequence and returns the amino acid sequence
    that will be obtained in the longest possible protein
    :param rna_seq: a given sequense of RNA
    :param readin_frame: The location from which the translation will start 
    :return str_amino_acids_with_tags: The resulting amino acid sequence
                                    in a convenient printable string
    """
def translate(rna_sec, reading_frame):

    str_amino_acids = ""                    #A string of the amino acids that the function will return

    dict_of_amino_acids = {                 #A dictionary of the letters that mark the amino acids
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "", "UAG": "",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

    start_codon = "AUG"
    stop_codons = ["UAA", "UAG", "UGA"]
    
    length_protein = 0          #saves the length of cuurent protein 
    max_length = 0              #saves the length of longest protein
    start_translate = 0         #saves the location of start codon
    end_translate = 0           #saves the location of stop codon or end translate

    #At first we will check what is the longest possible protein sequence :max_length:
    #And then we will translate it to the amino acids
    
    for index in range(reading_frame,len(rna_sec) - (3+reading_frame) ,3) :
                                                  
        if rna_sec[index:index+3] == start_codon :
            length_protein = 1
            for n in range(1,len(rna_sec)//3) :
                if index+(3*(n)) > len(rna_sec) - 3 : #We arrivade at the end of the sequence
                    if length_protein > max_length:
                        max_length = length_protein
                        start_translate = index
                        end_translate = index+3*(max_length)
                    break
                
                elif rna_sec[index+3*n: index+3*(n+1)] in stop_codons: #We arrived at a stop code
                    if length_protein > max_length:
                        max_length = length_protein
                        start_translate = index
                        end_translate = index+3*(max_length)  
                    break
                
                else:
                    length_protein+=1

#Now we will translate the longest sequence into amino acids
                    
    if max_length > 0 :
        for j in range (start_translate, end_translate, 3) :
            codon = rna_sec[j : j+3]
            amino_acid = dict_of_amino_acids[codon]
            str_amino_acids += amino_acid
#Now we will transfer the amino acids from the list
#            to a convenient printable string

    tag = ";"
    str_amino_acids_with_tags = tag.join(str_amino_acids)
    
    return str_amino_acids_with_tags
            
    
def translate_printer(rna_argv, reading_frame_argv):
    """
    The function prints the longest possible amino acid sequence
    from the RNA sequence we received
    :param rna_argv: a given sequense of RNA
    :param reading_frame_argv: The location from which the translation will start 
    """
    translate_to_print = translate(rna_argv, int(reading_frame_argv))
    if translate_to_print == "":
        print("Non-coding RNA")
    else:
        print("Translation: " + translate_to_print) 

#-------------------------------------------------------------------------

if __name__ == '__main__':
    ssr_printer(sys.argv[1])
    transcribe_printer(sys.argv[2])
    translate_printer(sys.argv[3], sys.argv[4])

