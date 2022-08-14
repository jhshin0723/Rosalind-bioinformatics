def find_splice(dna_motif, dna):
  #Accepts two strings, returns a list of index positions of the motif characters for the first occurence
  #Once first letter is found, start looking for next motif letter starting at that index

  poslist = [] #output splice positions

  for index,letter in enumerate(dna_motif):
    for index1,letter1 in enumerate(dna): 
      if letter1 == dna_motif[index]:  
        poslist.append(index1) #once matching dna sequence letter is found, add to list
        dna = dna[:index1] + '-' + dna[index1+1:] #replace the found letter with a '-'
        break 
      else:
        dna = dna[:index1] + '-' + dna[index1+1:] #if sequence let doesn't match motif let, replace with '-'
  if len(poslist)!=len(dna_motif): 
    poslist = []
  return poslist
  
def shared_motif(dna_list): 
  #Accepts a list of dna strings, and outputs longest common substring 
  com_substr = '' #common substring
  if len(dna_list) > 1 and len(dna_list[0]) > 0:
    #The for loops here go through all possible starting letters of common substring, and 
    #all lengths of possible substring, within the first dna string
    #Then creates a substring with the starting index and the length 
    for start_i in range(len(dna_list[0])): 
      for lensub in range(len(dna_list[0]) - start_i + 1):
        if lensub > len(com_substr) and all(dna_list[0][start_i:start_i + lensub] in dna for dna in dna_list):
          com_substr = dna_list[0][start_i:start_i + lensub]
          #If the current substring is in all dna strings, and longer than previous com_str
          #make it the com_substr
  return com_substr


def get_edges(dict):
        key_list = dict.keys()
        list=[]
        for i in key_list:
            list.append(i)
        list2=[]
        for i in range(0,len(list)):
            for j in range(i+1,len(list)):

                if(dict[list[i]][:3]==dict[list[j]][-3:] or dict[list[i]][-3:]==dict[list[j]][:3]):
                    list2.append((list[i],list[j]))          
        return list2
    
dict={
    
    "Rosalind_0498":"AAATAAA",
    "Rosalind_2391":"AAATTTT",
    "Rosalind_2323":"TTTTCCC",
    "Rosalind_0442":"AAATCCC",
    "Rosalind_5013":"GGGTGGG"
    
}

adjencey_list=get_edges(dict)
print (adjencey_list)

def assemble_genome(words):
  n = len(words)
  overlaps = [[0 for _ in range(n)] for _ in range(n)]
  for i in range(n):
    for j in range(n):
      if i == j:
        continue
      x, y = words[i], words[j]
      size = len(x)
      for k in range(1, size):
        if y.startswith(x[k:]):
          overlaps[i][j] = size - k
          break

  def helper(i, mask):
    if mask == (1<<n) - 1:
      return words[i]
    ans = '#' * 320
    for j in range(n):
      if mask & (1<<j) == 0:
        k = overlaps[i][j]
        string = helper(j, mask | (1<<j))
        if len(words[i] + string[k:]) < len(ans):
          ans = words[i] + string[k:]
    return ans
  return min([helper(i, 1<<i) for i in range(n)], key=len)

def reverse_complement(dna): #This function is used in rev_palindrome, but from previous milestone
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
def rev_palindrome(dna):
   res = []
   for i in range(len(dna)-4):
       for j in range(i+3,min(len(dna),i+12)):
           s = dna[i:j+1]
           if reverse_complement(dna[i:j+1]) == s:
               res.append((i,j-i+1))

   return res

def random_genome(dna, gc_content): 
  #Accepts dna string and gc_content float list, and returns a list of the logarithm 
  # values of the probabilities that a random string constructed with the GC-content
  # will match the given DNA string exactly. 
  import math #import math
  loglist = [] #output float list
  cg_freq = 0
  for i in gc_content: #for each value in gc_content, calculate frequencies
    cg_freq = (i/2)
    at_freq = (1-i)/2
    result = 1
    for let in dna: 
      #for loop multiplies result by either cg_freq or at_freq 
      # depending on what the letter in dna is
      if let == 'A' or let == 'T':
        result = result * at_freq
      else:
        result = result * cg_freq
    loglist.append(math.log10(result))
  return (loglist)

def perfect_match(rna):
  import math  
  d={"A":0,"C":0,"G":0,"U":0}
  for i in rna:
    d[i]+=1
  if d["A"]==d["U"] and d["C"]==d["G"]:
    return math.factorial(d["A"])*math.factorial(d["C"])
  else:
    return 0
