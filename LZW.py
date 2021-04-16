import sys
import numpy as np
from collections import deque


MIN_MATCH_LEN = 2       #heuristic, no basis

def main(argv):

    if(len(argv) == 0):
        print("Usage: CMD [input file] -[options]")

    fp = open(argv[0], "r")

    full_fasta = fp.read()
    full_fasta = full_fasta.split("\n")
    
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
        elif char == 'T':  # remove from final version
            key_flipped = key_flipped + 'A'
    return key_flipped


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
    print(str(len(dict_47.keys())) + ' keys of len 4 - 7 found instances with an avg of ' + str(sum(per_vals)/len(dict_47.keys())) + '% accuracy')

    per_vals = []
    for key in dict_811:
        dict_num = len(dict_811[key][1])
        seq_num = seq.count(key)
        if dict_num > seq_num:
            dict_notable[key] = [dict_num,seq_num]
        else:
            per_vals.append((dict_num/seq_num)*100)
    print(str(len(dict_811.keys())) + ' keys of len 8 - 11 found instances with an avg of ' + str(sum(per_vals)/len(dict_811.keys())) + '% accuracy')

    per_vals = []
    for key in dict_12plus:
        dict_num = len(dict_12plus[key][1])
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


def stem_match(dictionary, threshold, seq):
    stem_keys = []
    for key in dictionary:
        if len(key) >= threshold:
            stem_keys.append([key,'FALSE'])

    print(str(len(stem_keys)) + ' keys over length ' + str(threshold))

    for ele in stem_keys:
        if seq.count(flip_key(ele[0])) > 0 or seq.count(flip_key(ele[0])[::-1]) > 0:
            ele[1] = 'SHOULD BE TRUE'
        if flip_key(ele[0]) in dictionary:
            ele[1] = 'TRUE1'
        elif flip_key(ele[0])[::-1] in dictionary:
            ele[1] = 'TRUE2'

    return stem_keys

#returns 0 on no match 
#NO INLINE MATCH 1 on inverted match            AUCG -> UAGC
#1 on fliped+inverted match     AUGC -> CGAU
def check_match(seq, index, dictionary):
    try:
        if(index not in dictionary[flip_key(seq[::-1])][1]):
            dictionary[flip_key(seq[::-1])][1].append(index)
        
    except:
        return 0

def LZW(seq):
    if(len(seq) < 2):
       return 0 
    # print('SEQUENCE: ' + seq)
    #initialize the dictionary to contain all strings of length 1
    dictionary = {}
    dict_val = 0

    #Dictionary is [key] -> [dict entry #, match list, appearence list]
    if 'T' in seq:
        for char in "ATCG":
            dictionary[char] = [dict_val, [], []]
            dict_val = dict_val+1
    else:
        for char in "AUCG":
            dictionary[char] = [dict_val, [], []]
            dict_val = dict_val+1

    output = list(seq) #Split output into a char array, for keeping track of match indecies 
    for i in range(50):
        buff = ""
    #  buff = buff + seq[1][0]
        for index, char in enumerate(seq):
            buff = buff + char

            #find the longest string W that matches the current input
            if(buff in dictionary):
                if (len(buff) > MIN_MATCH_LEN):
                    if(index not in dictionary[buff][2]):
                        dictionary[buff][2].append(index)
                    # check if in dict, if in dict add to appearence 
                    check_match(buff, index,dictionary)    #If the inverse is in the dictionary,add the current index to its match list

            else:
                try: #debugging try catch, not control flow try catch :p
                    #emit the dictionary index for W to output and remove W from the input
                    if(len(buff) > MIN_MATCH_LEN):
                        check_match(buff, index, dictionary)

                    #add W followed by the next symbol in the input to the dictionary
                    dictionary[buff] = [dict_val,[], [index]]  #first entry in dictionary = first occurance of sequence
                    dict_val = dict_val + 1
                    #remove W from input
                    buff = buff[len(buff)-1:]
                except Exception as E:
                    print("BROKE\n\n")
                    print(str(dictionary))
                    print(buff)
                    print(str(output))
                    raise(E)
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

    for key in dictionary:  # remove repeat instances of indexes
        temp_list = []
        for i in dictionary[key][1]:
            if i not in temp_list:
                temp_list.append(i)
        dictionary[key][1] = temp_list

   # print('DICTIONARY: ' + str(dictionary))
    small_dict = {}
    for key in dictionary:
        if(len(dictionary[key][1]) > 0):
            #print(key, dictionary[key])
            small_dict[key] = dictionary[key]   #remove all bad entries, makes the next step faster


    print_dot_brace(output, small_dict)
    #print('POTENTIAL STEMS: ' + str(stem_match(dictionary,8,seq)))
    #print(str(output));
    # test_metrics(dictionary, seq)
    # opposing_keys_check(dictionary, seq)

#[idx, match, appear]
def print_dot_brace(output, dictionary):
    
    #TODO this finds 1 best match for the given sequence, even when the sequence might have 8 locations and 8 matches, meaning that it could have up to 8 good stems, fixing that would be some nasty optimization maybe
    #orrr, add each dictionary tuples, where you pick the best match tuple, remove those, find the next best match tuple, until one list is empty, then later go through the tuples.. Ill do that next.
    to_remove = []
    for key in dictionary:
        match_idx = -1
        loc_idx = len(output) + 1
        for match in dictionary[key][1]:
            for loc in dictionary[key][2]:
                if(abs(match-loc) < len(key)):     #if the found match overlaps itself, UGAUAGCA does this. (inverted palindrome?)
                    continue
                if(abs(match-loc) < abs(match_idx - loc_idx)):    #if the distance between the new pair is less
                    match_idx = match
                    loc_idx = loc           #replace previous match with closer match (assumming closer structures of same length are more likely to form)
        if(match_idx == -1):
            to_remove.append(key)
            continue;   #no match found, probably since only matches were overlaps, remove and continue
        dictionary[key][1] = match_idx
        dictionary[key][2] = loc_idx    #replace dictionary with only the best matches

    for key in to_remove:
        del dictionary[key]

    sorted_dict = sorted(dictionary.items(), key=lambda k:len(k[0]), reverse=True)  #actually an array of [dict keys, [vals]]
    #print(sorted_dict)
    for key in sorted_dict:
        #print(len(key[0]))
        match_clr = check_output(output, key[1][1], len(key[0]))
        app_clr = check_output(output, key[1][2], len(key[0]))  #this is checking the match "backwards" but since we just wanna see if its clear, its ok
        if(match_clr == -1 and app_clr == -1):
            #print("test", key[1][1], key[1][2])
            set_output(output, max(key[1][1], key[1][2]), len(key[0]), (max(key[1][1], key[1][2]), 1))
            set_output(output, min(key[1][1], key[1][2]), len(key[0]), (max(key[1][1], key[1][2]), 0))
        else:
            z = 0
            #there is a stem overlapping what we want to fill in, which is gonna be smaller or equal since we are iterating in order of length
            #TODO make this maybe check same length stems to see if another one is close?
        #print('|'.join([str(k) for k in output]))

    #At this point, the output array contains matched sets of characters for each stem, in essence being dot-bracket with a variation of bracket for each stem, and emptyness for the dots. The following code is my attempt at an elegant algorithm for simplifing down to a minimum number of bracket styles while maintaining proper matching
    #the numbers for the brackets are the lowest (left-most) index of the right matched pair so '..((..))..' would be represented at this point in the code as '..66..66..' (but with whatever random bases were there as the dots) 
    #print('|'.join([str(k) for k in output]))
    prev_val = 0
    stack = []
    brace_style = {}
    for idx, val in enumerate(output):
        if(type(val) == str):
            continue            #skip gaps entirely
        if(val == prev_val):
            continue
        prev_val = val  #glides over continuous streaks of the same character
        if(val not in brace_style):
            brace_style[val] = 0

        if(peek(stack, val)):  #top of stack (through transparent) matches, pop
            while(stack[-1][1] == 0):
                stack=stack[:-1]        #pop as many trasnparent entries as possible
            stack = stack[:-1]

        elif(idx ==  val):      #curr idx is a right match, but not on top of stack, pseudoknot found
            for i in stack[::-1]: #iterate backwards through the stack, finding the largest bracket style
                max_brace = 0 
                if(i == val):
                    break   #we found the matching brace
                    max_brace = max(max_brace, brace_style[i])
            brace_style[val] = max_brace+1
            #push transparent 
            stack.append((val, 0))
            #make previous appearence transparent as well
            for (i, check_val) in enumerate(stack[::-1]):
                if (check_val == val):
                    stack[-(i+1)] = (val, 0)    #TODO early terminate this for performance?

        else:                   #non-right, non-matched thing, push opaque to stack
            stack.append((val, 1))


    print(brace_style)
    #now that we have brace styles, go through and 1 for 1 swap
    for idx, val in enumerate(output):
        if(type(val) == str):
            output[idx] = '.'
        elif val[1] == 1:
            output[idx] = brace_close(brace_style[val])
        else:
            output[idx] = brace_open(brace_style[val])

    print(''.join(output))

#peek through transparent values
def peek(stack, val):
    for j in stack[::-1]:
        if j[0] == val:
            return 1
        if j[1] == 1:
            return 0
    return 0

def brace_open(i):
    braces = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return braces[i] #if we need more than this, send help

def brace_close(i):
    braces = ")]}>abcdefghijklmnopqrstuvwxyz"
    return braces[i]

def set_output(output, start, length, val):
    for i in range(length):
        output[start-length+i+1] = val

def check_output(output, start, length):
    for i in range(start -length, start):
       # print(len(output), start, length, output[i])
        if (not type(output[i]) == str):
            return i
    return -1   #This is bad since -1 is the "good" case, but

if __name__ == "__main__":
    main(sys.argv[1:])
