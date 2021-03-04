import sys
import numpy as np


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


# flip given sequence
def flip_key(key):
    key_flipped = ''
    for char in key:
        if char == 'A':
            key_flipped = key_flipped + 'U'
        elif char == 'U':
            key_flipped = key_flipped + 'A'
        elif char == 'C':
            key_flipped = key_flipped + 'G'
        elif char == 'G':
            key_flipped = key_flipped + 'C'
    return key_flipped


# check if opposing key exists, and how many times
def opposing_keys(dictionary, seq):
    for key in dictionary:
        try:
            print(key + ' exists in dictionary ' + str(len(dictionary[key][1])) + ' times. ' + flip_key(key) + ' exists in dictionary ' + str(len(dictionary[flip_key(key)][1])) + ' times. ')
            print(key + ' exists in sequence ' + str(seq.count(key)) + ' times. ' + flip_key(key) + ' exists in sequence ' + str(seq.count(flip_key(key))) + ' times. ')
        except:
            print(key + ' exists in dictionary ' + str(len(dictionary[key][1])) + ' times. ' + flip_key(key) + ' does not exist in dictionary.')
            print(key + ' exists in sequence ' + str(seq.count(key)) + ' times. ' + flip_key(key) + ' exists in sequence ' + str(seq.count(flip_key(key))) + ' times. ')


def LZW(seq):
    if(len(seq) < 2):
       return 0 
    # print('SEQUENCE: ' + seq)
    #initialize the dictionary to contain all strings of length 1
    dictionary = {}
    dict_val = 0
    for char in "AUCG":
        dictionary[char] = [dict_val, []]
        dict_val = dict_val+1

    buff = ""
#    buff = buff + seq[1][0]
    output = ""
    for index, char in enumerate(seq):
        buff = buff + char
        #find the longest string W that matches the current input
        try:
            dictionary[buff][1].append(index)
            # check if in dict, if in dict add last index?
        
        except:
            try: #debugging try catch, not control flow try catch :p
                #emit the dictionary index for W to output and remove W from the input
                output = output + str(dictionary[buff[:-1]]) + " "
                #add W followed by the next symbol in the input to the dictionary
                dictionary[buff] = [dict_val, [index]]  #TODO what do we associate with it (value metric + location of last match?
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

    print('DICTIONARY: ' + str(dictionary))
    opposing_keys(dictionary, seq)
    # print('OUTPUT: ' + output)


if __name__ == "__main__":
    main(sys.argv[1:])
