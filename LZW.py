import sys
import numpy as np


struct_list = []


def main(argv):

    # if(len(argv) == 0):
    #     print("Usage: CMD [input file] -[options]")
# 
    # fp = open(argv[0], "r")
# 
    # full_fasta = fp.read()
    # full_fasta = full_fasta.split("\n")
    # 
    # #remove first line of metadata, join rest w/no newline
    # seq = "".join(full_fasta[1:])
    seq = 'UUGUCACUGGACGAAGUGAAUGGGUCAAAUGGGCUUGUCUAAGUUCCGACCCAGGAAAGUCCCCGGGCUAUUCUCCGCAGAGAUGCGGCCUUGUCAACGAGAGAGUCAUACGGUGGAGUCGAUCCAAUCAGAGCUGAGACGUGUGUAUGGAGCUCACGUGGUCCCUCUGCCGAUGAUAUCAUGGCCGUGAUAACACAGUUAUUCACCUGGUUCUUAAUUGGACUAACCGAACGGGGCUUACGUUCCAGUGACA'
    dictionary = LZW(seq)
    dot_bracket = ['-' for i in seq]
    # print(str(dot_bracket))

    all_dot_brackets(0,dot_bracket,dictionary,5)

    global struct_list
    # print(str(struct_list))
    string_list = []
    for struct in struct_list:
        dot_bracket = ''
        for ele in struct:
            dot_bracket+=ele
        string_list.append(dot_bracket)
    print(str(string_list))
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
        elif char == 'T':  # remove from final version
            key_flipped = key_flipped + 'A'
    return key_flipped[::-1]


def expand_dict(dictionary,threshold):
    for x in range(5):
        dict_copy = dictionary.copy()
        for key in dictionary:  # iterate through each key in dict
            if len(key) >= threshold:  # if length of key is over the threshold, begin to break down
                if key[1:] in dictionary:  # if subkey removing first ele of key is already in dictionary, append the indeces of key location to list
                    dict_copy[key[1:]] = list(dict.fromkeys(dictionary[key[1:]] + dictionary[key]))
                else:  # else, create dictionary entry for subkey which contains curr keys indeces
                    dict_copy[key[1:]] = list(dict.fromkeys(dictionary[key]))

                to_append = [i - 1 for i in dictionary[key]]  # alter indeces for shortened key

                if key[:-1] in dictionary:  # if subkey removing last ele of key is already in dictionary, append the altered indexes of key location to list
                    dict_copy[key[:-1]] = list(dict.fromkeys(dictionary[key[:-1]] + to_append))
                else:  # else, create dictionary entry for subkey which contains altered indeces
                    dict_copy[key[:-1]] = list(dict.fromkeys(to_append))
        dictionary = dict_copy
    dict_copy = dictionary.copy()
    for key in dictionary:
        if len(key) < threshold:
            del dict_copy[key]
    return dict_copy


def flip_check(key,dictionary):
    try:
        dictionary[flip_key(key)]
        return True
    except:
        return False


def all_dot_brackets(i,dot_bracket,dictionary,threshold):
    # if '-' not in dot_bracket:
    if i >= len(dot_bracket) - 1 or '-' not in dot_bracket:
        print('FINISHED A DB')
        global struct_list
        struct_list.append(dot_bracket)
        return

    for key in dictionary:
        print('BEGINNING SEARCH FOR ' + key + ' and ' + flip_key(key) + ' which exists: ' + str(flip_check(key,dictionary)))
        if len(key) >= threshold and flip_check(key,dictionary):
            print(key + ' is ' + str(len(key)) + ' long')
            for index in dictionary[key]:
                location_start = index - len(key)
                dash_count = dot_bracket[index - len(key):index].count('-')
                len_segment = len(dot_bracket[index - len(key):index])
                if location_start == i and dash_count == len_segment:
                    for match_index in dictionary[flip_key(key)]:
                        if dot_bracket[match_index - len(flip_key(key)):match_index].count('-') == len(dot_bracket[match_index - len(flip_key(key)):match_index]):
                            db_copy = dot_bracket.copy()
                            if index < match_index:
                                for x in range(index - len(key), index): db_copy[x] = '('
                                for x in range(match_index - len(key), match_index): db_copy[x] = ')'
                            else:
                                for x in range(index - len(key), index): db_copy[x] = ')'
                                for x in range(match_index - len(key), match_index): db_copy[x] = '('
                            print(str(db_copy))
                            all_dot_brackets(index + 1,db_copy,dictionary,threshold)
    db_copy = dot_bracket.copy()
    # print(str(db_copy))
    db_copy[i] = '.'
    all_dot_brackets(i + 1,db_copy,dictionary,threshold)


def LZW(seq):
    if(len(seq) < 2):
       return 0 
    # print('SEQUENCE: ' + seq)
    #initialize the dictionary to contain all strings of length 1
    dictionary = {}
    # dict_val = 0

    for char in "AUCG":
        dictionary[char] = []
        # dict_val = dict_val+1

    for i in range(5):
        buff = ""
    #  buff = buff + seq[1][0]
        output = []
        for index, char in enumerate(seq):
            buff = buff + char
            #find the longest string W that matches the current input
            try:
                dictionary[buff].append(index)
                # check if in dict, if in dict add last index?

            except:
                try: #debugging try catch, not control flow try catch :p
                    #emit the dictionary index for W to output and remove W from the input
                    output.append(buff[:-1])
                    # output.append(dictionary[buff[:-1]])
                    #add W followed by the next symbol in the input to the dictionary
                    dictionary[buff] = [index]  #TODO what do we associate with it (value metric + location of last match?
                    # dict_val = dict_val + 1
                    #remove W from input
                    buff = buff[len(buff)-1:]
                except:
                    print("BROKE\n\n")
                    print(str(dictionary))
                    print(buff)
                    print(str(output))
                    return
            #GOTO step 2 (loop bottom)

    for key in dictionary:  # remove repeat instances of indexes
        temp_list = []
        for i in dictionary[key]:
            if i not in temp_list:
                temp_list.append(i)
        dictionary[key] = temp_list

    # print('DICTIONARY: ' + str(dictionary))
    # print('POTENTIAL STEMS: ' + str(stem_match(dictionary,8,seq)))
    dictionary = expand_dict(dictionary,5)
    return dictionary
    # opposing_keys_check(dictionary, seq)


if __name__ == "__main__":
    main(sys.argv[1:])
