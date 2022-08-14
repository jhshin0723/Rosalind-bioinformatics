
def s(dna):
    dna = dna.upper()
    count_A = dna.count('A')
    count_C = dna.count('C')
    count_G = dna.count('G')
    count_T = dna.count('T')
    dict_dna = {}
    dict_dna['A'] =count_A
    dict_dna['C'] =count_C
    dict_dna['G'] =count_G
    dict_dna['T'] =count_T
 
    return dict_dna



def dna2rna(dna):
    rna = ''
    for symbol in dna:
        if symbol == 'A':
            rna = rna + 'A'
        elif symbol == 'T':
            rna = rna + 'U'
        elif symbol == 'G':
            rna = rna + 'G'
        elif symbol == 'C':
            rna = rna + 'C'
            
    return rna

def reverse_complement(dna):
    rev_dna =''
    for char in dna:
        rev_dna = char + rev_dna 
    
    rna = ''
    for symbol in rev_dna:
        if symbol == 'A':
            rna = rna + 'T'
        elif symbol == 'T':
            rna = rna + 'A'
        elif symbol == 'G':
            rna = rna + 'C'
        elif symbol == 'C':
            rna = rna + 'G'
            
    return rna

from math import factorial

def mendels_law(hom, het, rec):
    numAA = hom
    numAa = het
    numaa = rec
    AA_Aa = 0
    AA_aa = 0
    AA_AA = 0
    aa_aa = 0
    Aa_Aa = 0
    Aa_aa = 0
    
    numtot = 0
    numdom = 0
    if(numAA>=2):
      AA_AA = int(factorial(numAA)/(factorial(numAA-2)*factorial(2)))
    if(numaa>=2):
      aa_aa = int(factorial(numaa)/(factorial(numaa-2)*factorial(2)))
    if(numAa>=2):
      Aa_Aa = int(factorial(numAa)/(factorial(numAa-2)*factorial(2)))
    
    AA_Aa = numAA * numAa 
    AA_aa = numAA * numaa 
    Aa_aa = numAa * numaa

    while (AA_AA>0):
        numtot = numtot + 4
        numdom = numdom + 4
        AA_AA = AA_AA - 1
    
    while (aa_aa>0):
        numtot = numtot + 4
        aa_aa = aa_aa - 1
    
    while (Aa_Aa>0):
        numtot = numtot + 4
        numdom = numdom + 3
        Aa_Aa = Aa_Aa - 1
    
    while (AA_Aa>0):
        numtot = numtot + 4
        numdom = numdom + 4
        AA_Aa = AA_Aa - 1
    
    while (AA_aa>0):
        numtot = numtot + 4
        numdom = numdom + 4
        AA_aa = AA_aa - 1
    
    while (Aa_aa>0):
        numtot = numtot + 4
        numdom = numdom + 2
        Aa_aa = Aa_aa - 1

    
    
    Result = float(numdom / numtot)
    
    return Result
   





def fibonacci_rabbits(n,k):
	f1,f2=1,1
	for i in range(n-1):
		f2, f1= f1,f1 + (f2 * k)
	print(f2)
	return f2

def gc_content(dna_list):
	max_gc = 0.0
	index = -1
	for i in range(len(dna_list)):
		count = 0
		for j in range(len(dna_list[i])):
			if dna_list[i][j] == 'G' or dna_list[i][j] == 'C':
				count += 1
	gc_con = ((count*1.0) / len(dna_list[i])) * 100.0
	
	if gc_con > max_gc:
		index = i
		max_gc = gc_con
	return index, round(max_gc,6)




def rna2codon(seq):
  table = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
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
  protein ="" 
  if len(seq)%3 == 0: 
      for i in range(0, len(seq), 3): 
          codon = seq[i:i + 3]
          if table[codon]!="STOP":
            protein+= table[codon] 
          else:
              break
  return protein 






# import re
def locate_substring(dna_snippet, dna):
    # list1 = []
    ret = []
    for i in range(len(dna)):
        # if dna[i:i+len(dna_snippet)] == dna_snippet:
        tmp = dna[i:]
        currIdx = tmp.find(dna_snippet)
        if currIdx >= 0:
            if currIdx + i not in ret:
                ret.append(currIdx + i)
    #     return list1
    # for m in re.finditer(dna_snippet, dna):
    #     ret += [ m.start() ]
        # dna = dna[i:]
    return ret
print(locate_substring("ATAT", "GATATATGCATATACTT"))



def hamming_dist(dna1, dna2):
  diff_char = 0
  
  for i in range(len(dna1)):
    if(dna1[i:i+1]!=dna2[i:i+1]):
      diff_char = diff_char + 1

  return diff_char


def count_dom_phenotype(genotypes):
  numAA_AA = genotypes[0]
  numAA_Aa = genotypes[1]
  numAA_aa = genotypes[2]
  numAa_Aa = genotypes[3]
  numAa_aa = genotypes[4]
  numaa_aa = genotypes[5]
  numdom = 0
  
  while (numAA_AA>0):
    numdom = numdom + 2
    numAA_AA = numAA_AA - 1
    
  while (numaa_aa>0):
    numdom = numdom + 0
    numaa_aa = numaa_aa - 1
    
  while (numAa_Aa>0):
    numdom = numdom + 1.5
    numAa_Aa = numAa_Aa - 1
    
  while (numAA_Aa>0):
    numdom = numdom + 2
    numAA_Aa = numAA_Aa - 1
    
  while (numAA_aa>0):
    numdom = numdom + 2
    numAA_aa = numAA_aa - 1
  
  while (numAa_aa>0):
    numdom = numdom + 2
    numAa_aa = numAa_aa - 1

  return numdom


def source_rna(protein):
  numberapper ={
    'F':2,
    'L':6,
    'I':3,
    'M':1,
    'S':6,
    'T':4,
    'Y':2,
    'N':2,
    'K':2,
    'C':2,
    'W':1,
    'R':6,
    'V':4,
    'A':4,
    'H':2,
    'Q':2,
    'D':2,
    'E':4,
    'G':4,
    '*':3,
  }
  count = 3
  for i in protein:
      count = count*numberapper[i]
  return count

def splice_rna(dna,intron_list):
  for i in intron_list:
    while i in dna:
      dna = dna.replace(i," ")
  rna = dna2rna(dna)
  rna = rna.strip()
  output = rna2codon(rna)
  return output