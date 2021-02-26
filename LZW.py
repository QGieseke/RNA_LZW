import sys


def main(argv):

    if(len(argv) == 0):
        print("Usage: CMD [input file] -[options]")

    fp = open(argv[0], "r")

    full_fasta = fp.read()
    
    full_fasta = full_fasta.split("\n")
    
    #print(full_fasta)
    #remove first line of metadata, join rest w/no newline
    seq = "".join(full_fasta[1:])
    LZW(seq)
#   for seq in range(0, len(full_fasta), 2):
#       print(seq)
#       LZW(full_fasta[seq:seq+2])




def LZW(seq):
    if(len(seq) < 2):
       return 0 
    print(seq)
    #initialize the dictionary to contain all strings of length 1
    dictionary = {}
    dict_val = 0
    for char in "AUCG":
        dictionary[char] = dict_val
        dict_val = dict_val+1

    buff = ""
#    buff = buff + seq[1][0]
    output = ""
    for char in seq:
        buff = buff + char
        #find the longest string W that matches the current input
        try:
            dictionary[buff]
        
        except:
            try: #debugging try catch, not control flow try catch :p
                #emit the dictionary index for W to output and remove W from the input
                output = output + str(dictionary[buff[:-1]]) + " "
                #add W followed by the next symbol in the input to the dictionary
                dictionary[buff] = dict_val         #TODO what do we associate with it (value metric + location of last match?
                dict_val = dict_val + 1
                #remove W from input
                buff = buff[len(buff)-1:]
            except:
                print("BROKE\n\n")
                print(str(dictionary))
                print(buff)
                print(output)
                return
        #GOTO step 2 (loop bottom)

    print(str(dictionary))
    print(output)


if __name__ == "__main__":
    main(sys.argv[1:])
