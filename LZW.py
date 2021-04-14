import sys
import numpy as np


struct_list = []


def main(argv):

    if(len(argv) == 0):
        print("Usage: CMD [input file] -[options]")

    fp = open(argv[0], "r")

    full_fasta = fp.read()
    full_fasta = full_fasta.split("\n")
    
    #remove first line of metadata, join rest w/no newline
    seq = "".join(full_fasta[1:])
    dictionary = LZW(seq)
    dot_bracket = ['-' for i in seq]
    print(str(dot_bracket))

    all_dot_brackets(0,dot_bracket,dictionary,8)

    global struct_list
    print(str(struct_list))
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


def test_metrics(dictionary, seq):
    print('current dictionary size: ' + str(len(dictionary)))
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
    
    # per_vals = []
    # for key in dict_13:
    #     dict_num = len(dict_13[key][1])
    #     seq_num = seq.count(key)
    #     if dict_num > seq_num:
    #         dict_notable[key] = [dict_num,seq_num]
    #     else:
    #         per_vals.append((dict_num/seq_num)*100)

    # print(str(len(dict_13.keys())) + ' keys of len 1 - 3 found instances with an avg accuracy of ' + str(sum(per_vals)/len(dict_13.keys())) + '% accuracy')

    per_vals = []
    for key in dict_47:
        dict_num = len(dict_47[key])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)
    print(str(len(dict_47.keys())) + ' keys of len 4 - 7 found instances with an avg of ' + str(sum(per_vals)/len(dict_47.keys())) + '% accuracy')

    per_vals = []
    for key in dict_811:
        dict_num = len(dict_811[key])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)
    print(str(len(dict_811.keys())) + ' keys of len 8 - 11 found instances with an avg of ' + str(sum(per_vals)/len(dict_811.keys())) + '% accuracy')

    per_vals = []
    for key in dict_12plus:
        dict_num = len(dict_12plus[key])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)
    try:
        print(str(len(dict_12plus.keys())) + ' keys of len 12+ found instances with an avg of ' + str(sum(per_vals)/len(dict_12plus.keys())) + '%% accuracy')
    except Exception:
        pass

    # print('The following keys had an accuracy >100%')
    # for key in dict_notable.keys():
    #     print(key + ' exists in dictionary, seq ' + str(dict_notable[key][0]) + ', ' + str(dict_notable[key][1]) + ' times respectively')

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


def expand_dict(dictionary,threshold):
    for x in range(5):
        dict_copy = dictionary.copy()
        for key in dictionary:  # iterate through each key in dict
            if len(key) >= threshold:  # if length of key is over the threshold, begin to break down
                if key[1:] in dictionary:  # if subkey removing first ele of key is already in dictionary, append the indeces of key location to list
                    dict_copy[key[1:]] = list(dict.fromkeys(dictionary[key[1:]] + dictionary[key]))
                else:  # else, create dictionary entry for subkey which contains curr keys indeces
                    dict_copy[key[1:]] = list(dict.fromkeys(dictionary[key]))

                to_append = [i - 1 for i in dictionary[key[:-1]]]  # alter indeces for shortened key

                if key[:-1] in dictionary:  # if subkey removing last ele of key is already in dictionary, append the altered indexes of key location to list
                    dict_copy[key[:-1]] = list(dict.fromkeys(dictionary[key[:-1]] + to_append))
                else:  # else, create dictionary entry for subkey which contains altered indeces
                    dict_copy[key[:-1]] = list(dict.fromkeys(to_append))
        dictionary = dict_copy
    return dictionary


def all_dot_brackets(i,dot_bracket,dictionary,threshold):
    # if '-' not in dot_bracket:
    if i == len(dot_bracket) - 1:
        global struct_list
        struct_list.append(dot_bracket)
        return

    for key in dictionary:
        try: 
            if len(key) >= threshold and dictionary[flip_key(key)]:
                for index in dictionary[key]:
                    if index - len(key) -1 == i and dot_bracket[index - len(key) -1:index].count('-') == len(dot_bracket[index - len(key) -1:index]):
                        for match_index in dictionary[flip_key(key)]:
                            if dot_bracket[match_index - len(key) -1:match_index].count('-') == len(dot_bracket[match_index - len(key) -1:match_index]):
                                db_copy = dot_bracket.copy()
                                if index < match_index:
                                    for x in range(index - len(key) -1, index): db_copy[x] = '('
                                    for x in range(match_index - len(key) -1, match_index): db_copy[x] = ')'
                                else:
                                    for x in range(index - len(key) -1, index): db_copy[x] = ')'
                                    for x in range(match_index - len(key) -1, match_index): db_copy[x] = '('
                                # print(str(db_copy))
                                all_dot_brackets(index + 1,db_copy,dictionary,threshold)
        except:
            continue
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
    dict_val = 0

    if 'T' in seq:
        for char in "ATCG":
            dictionary[char] = [dict_val, []]
            dict_val = dict_val+1
    else:
        for char in "AUCG":
            dictionary[char] = [dict_val, []]
            dict_val = dict_val+1

    for i in range(20):
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
                    dict_val = dict_val + 1
                    #remove W from input
                    buff = buff[len(buff)-1:]
                except:
                    print("BROKE\n\n")
                    print(str(dictionary))
                    print(buff)
                    print(str(output))
                    return
            #GOTO step 2 (loop bottom)

        # for key in dictionary:  # only in loop for testing
        #     temp_list = []
        #     for i in dictionary[key][1]:
        #         if i not in temp_list:
        #             temp_list.append(i)
        #     dictionary[key][1] = temp_list

        # print('Pass ' + str(k))
        # test_metrics(dictionary, seq)

    # for key in dictionary:  # remove repeat instances of indexes
    #     temp_list = []
    #     for i in dictionary[key]:
    #         if i not in temp_list:
    #             temp_list.append(i)
    #     dictionary[key][1] = temp_list

    # print('DICTIONARY: ' + str(dictionary))
    # print('POTENTIAL STEMS: ' + str(stem_match(dictionary,8,seq)))
    dictionary = expand_dict(dictionary,5)
    return dictionary
    # opposing_keys_check(dictionary, seq)


if __name__ == "__main__":
    main(sys.argv[1:])
