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


def test_metrics(dictionary, seq):
    # for keys of size 1-3, 4-7, 8-11, 12+
    dict_13 = {}
    dict_47 = {}
    dict_811 = {}
    dict_12plus = {}
    dict_notable = {}

    for key in dictionary:
        if len(key) >= 1 and len(key) <= 3:
            dict_13[key] = dictionary[key]
        if len(key) >= 4 and len(key) <= 7:
            dict_47[key] = dictionary[key]
        if len(key) >= 8 and len(key) <= 11:
            dict_811[key] = dictionary[key]
        if len(key) >= 12:
            dict_12plus[key] = dictionary[key]
    
    per_vals = []
    for key in dict_13:
        dict_num = len(dict_13[key][1])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)

    print(str(len(dict_13.keys())) + ' keys of len 1 - 3 found instances with an avg accuracy of ' + str(sum(per_vals)/len(dict_13.keys())) + '% accuracy')

    per_vals = []
    for key in dict_47:
        dict_num = len(dict_47[key][1])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)
    print(str(len(dict_47.keys())) + ' keys of len 4 - 7 found instances with an avg accuracy of ' + str(sum(per_vals)/len(dict_47.keys())) + '% accuracy')

    per_vals = []
    for key in dict_811:
        dict_num = len(dict_811[key][1])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)
    print(str(len(dict_811.keys())) + ' keys of len 8 - 11 found instances with an avg accuracy of ' + str(sum(per_vals)/len(dict_811.keys())) + '% accuracy')

    per_vals = []
    for key in dict_12plus:
        dict_num = len(dict_12plus[key][1])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)
    print(str(len(dict_12plus.keys())) + ' keys of len 12+ found instances with an avg of ' + str(sum(per_vals)/len(dict_12plus.keys())) + '%% accuracy')

    print('The following keys had an accuracy >100%')
    for key in dict_notable.keys():
        print(key + ' exists in dictionary, seq ' + str(dict_notable[key][0]) + ', ' + str(dict_notable[key][1]) + ' times respectively')

    max_len = -1
    max_list = []
    for key in dictionary.keys(): 
        if len(key) > max_len: 
            max_len = len(key) 
            max_list = [key]
        elif len(key) == max_len:
            max_list.append(key)
    
    print('the longest subsequences found were ' + str(max_list))


# check if opposing key exists, and how many times
def opposing_keys_check(dictionary, seq):
    for key in dictionary:
        try:
            print(key + ' exists in dictionary ' + str(len(dictionary[key][1])) + ' times. ' + flip_key(key) + ' exists in dictionary ' + str(len(dictionary[flip_key(key)][1])) + ' times. ')
            # comment out the below line for faster output
            print(key + ' exists in sequence ' + str(seq.count(key)) + ' times. ' + flip_key(key) + ' exists in sequence ' + str(seq.count(flip_key(key))) + ' times. ')
        except:
            print(key + ' exists in dictionary ' + str(len(dictionary[key][1])) + ' times. ' + flip_key(key) + ' does not exist in dictionary.')
            # comment out the below line for faster output
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

    for i in range(10):
        buff = ""
    #  buff = buff + seq[1][0]
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

    for key in dictionary:
        temp_list = []
        for i in dictionary[key][1]:
            if i not in temp_list:
                temp_list.append(i)
        dictionary[key][1] = temp_list

    print('DICTIONARY: ' + str(dictionary))
    test_metrics(dictionary, seq)
    # opposing_keys_check(dictionary, seq)
    # print('OUTPUT: ' + output)


if __name__ == "__main__":
    main(sys.argv[1:])
