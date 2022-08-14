def load_file(filename):
  listDNA = []
  datafile = open( filename,'r' )
  data = datafile.read()
  datafile.close()
  data_as_lines = data.split( '\n' )
  for line in data_as_lines:
    listDNA.append(str(line))
  return listDNA



def assemble_genome2(words):
    superstring = words[0]
    print (words)
    i = 0
    while i < len(words)-1:
        if superstring[len(superstring)-8:] == words[i+1][0:8]: # check if last 8 letters match first 8 letter of next string
            superstring = superstring + words[i+1][8:]
            i = i + 1

        elif superstring[0:8] == words[i+1][len(words[i+1])-8:]:
            superstring = words[i+1][0:len(words[i+1])-8] + superstring
            i = i + 1

        else:
            words = words + [words[i+1]]
            i = i + 1
    return superstring
