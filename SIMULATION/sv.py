


# insert 
def find_discontinuous_keys(k, original_dict, tolerance):
    previous_discontinuous_dict = {}
    next_discontinuous_dict = {}

    keys = list(original_dict.keys())
    values = list(original_dict.values())

    for i in range(len(keys)):
        key = keys[i]
        value = values[i]

        if i > 0 and key != keys[i - 1] + 1:
            if abs(value[0] - values[i - 1][0]) <= tolerance and\
                abs(value[0] - values[i - 1][0])>=k/2:
                next_discontinuous_dict[key] = value

        if i < len(keys) - 1 and key != keys[i + 1] - 1:
            if abs(value[0] - values[i + 1][0]) <= tolerance and\
                abs(value[0] - values[i + 1][0])>=k/2:
                previous_discontinuous_dict[key] = value

    return previous_discontinuous_dict, next_discontinuous_dict


def compare_and_check_equal(data, start_i, start_j, count):
    equal_count = 0
    total_comparisons = 0

    for i in range(start_i, start_i + count + 1):
        j = start_j + (i - start_i)  # Calculate corresponding j
        total_comparisons += 1
        if data[i] == data[j]:
            equal_count += 1

    if equal_count / total_comparisons > 0.8:
        return True
    else:
        return False

def insert_dicts(chr_number,k,data,previous_dict, next_dict):
    items_previous = list(previous_dict.items())
    items_next = list(next_dict.items())

    results = {}

    for i in range(len(items_previous)):
        key_previous, value_previous = items_previous[i]
        key_next, value_next = items_next[i]
       
        first_letters = ""  # Initialize an empty string to store the first letters

       
        first_letters = ''.join(data[key][0][0] for key in range(key_previous + k-1, key_next))

        result = 'insert'+str(value_previous[0]+k-1)
        results[result] = (chr_number,value_previous[0]+k-1, first_letters)

    return results

# repeat
def find_repeat_keys(k, original_dict, tolerance):
    previous_discontinuous_dict = {}
    next_discontinuous_dict = {}

    keys = list(original_dict.keys())
    values = list(original_dict.values())

    for i in range(len(keys)):
        key = keys[i]
        value = values[i]

        if i > 0 and key != keys[i - 1] + 1:
            if abs(value[0] - values[i - 1][0]) <= tolerance and\
                abs(value[0] - values[i - 1][0])>0:
                next_discontinuous_dict[key] = value

        if i < len(keys) - 1 and key != keys[i + 1] - 1:
            if abs(value[0] - values[i + 1][0]) <= tolerance and\
                abs(value[0] - values[i + 1][0])>0:
                previous_discontinuous_dict[key] = value

    return previous_discontinuous_dict, next_discontinuous_dict


def compare_repeat(data, start_i, start_j, count):
    equal_count = 0
    total_comparisons = 0

    for i in range(start_i, start_i + count + 1):
        j = start_j + (i - start_i)  # Calculate corresponding j
        total_comparisons += 1
        if data[i] == data[j]:
            equal_count += 1

    if equal_count / total_comparisons > 0.8:
        return True
    else:
        return False

def find_repeats(dna_sequence, similarity_threshold):
    length = len(dna_sequence)
    for size in range(1, length // 2 + 1):
        repeat = dna_sequence[:size]
        repeat_count = dna_sequence.count(repeat)
        repeat_length = size * repeat_count
        similarity = repeat_length / length
        if similarity >= similarity_threshold:
            return repeat_count
    return 1


def repeat_dicts(chr_number,k,data,previous_dict, next_dict):
    items_previous = list(previous_dict.items())
    items_next = list(next_dict.items())

    results = {}

    for i in range(len(items_previous)):
        key_previous, value_previous = items_previous[i]
        key_next, value_next = items_next[i]
        # count = key_next - key_previous - 1
        # start_i = key_previous - count + 1
        # start_j = key_previous + 1
        first_letters = ""  # Initialize an empty string to store the first letters

        # for key in range(key_previous+k, key_next+k-1):
        #     first_sequence = data[key][0]  # Get the first sequence from the tuple
        #     first_letter = first_sequence[0]  # Get the first letter
        #     first_letters += first_letter  # Add the first letter to the string
        first_letters = ''.join(data[key][0][0] for key in range(key_previous+1,  key_next))

        dup_count =find_repeats(first_letters,0.8)
        result = 'repeat'+str(int(value_previous[0]+len(first_letters)/dup_count))
        # results[result] = (chr_number,int(value_previous[0]+len(first_letters)/dup_count), first_letters, dup_count)
        results[result] = (chr_number,int(value_previous[0]+len(first_letters)/dup_count), data[key_previous+1][0][0], dup_count)

    return results

# delete
def find_delete_keys(k,original_dict):
    previous_discontinuous_dict = {}
    next_discontinuous_dict = {}

    keys = list(original_dict.keys())
    values = list(original_dict.values())

    for i in range(len(keys)):
        key = keys[i]
        value = values[i]

        if i > 0 and key == keys[i - 1] + k:
            # if value[0] == values[i - 1][0] +k:
            next_discontinuous_dict[key] = value

        if i < len(keys) - 1 and key == keys[i + 1] - k:
            # if value[0] == values[i + 1][0] -k:
            previous_discontinuous_dict[key] = value

    return previous_discontinuous_dict, next_discontinuous_dict

def dele_dicts(chr_number,k,data,previous_dict, next_dict):
    items_previous = list(previous_dict.items())
    items_next = list(next_dict.items())

    results = {}

    for i in range(len(items_previous)):
        key_previous, value_previous = items_previous[i]
        key_next, value_next = items_next[i]
        # count = key_next - key_previous - 1
        # start_i = key_previous - count + 1
        # start_j = key_previous + 1
        first_letters = ""  # Initialize an empty string to store the first letters

        
        first_sequence = data[key_previous][0]  # Get the first sequence from the tuple
        first_letter = first_sequence[-1]  # Get the first letter
             # Add the first letter to the string
        # first_letters = (data[key][0][-1] for key in range(key_previous))

       
        result = 'del'+str(value_previous[0]+k-1)
        count= value_next[0]- value_previous[0]-k
        results[result] = (chr_number,value_previous[0]+k-1, first_letter, count)

    return results

# invert
def find_invert_keys(k,original_dict):
    previous_discontinuous_dict = {}
    next_discontinuous_dict = {}

    keys = list(original_dict.keys())
    values = list(original_dict.values())

    for i in range(len(keys)):
        key = keys[i]
        value = values[i]

        if i > 0 and key -(keys[i - 1] + 1)>=150 :
            if value[0] != values[i - 1][0] +1:
                next_discontinuous_dict[key] = value

        if i < len(keys) - 1 and (keys[i + 1] - 1)-key>=150 :
            if value[0] != values[i + 1][0] -1:
                previous_discontinuous_dict[key] = value

    return previous_discontinuous_dict, next_discontinuous_dict

def undefined_dicts(k,data,previous_dict, next_dict):
    items_previous = list(previous_dict.items())
    items_next = list(next_dict.items())

    

    for i in range(len(items_previous)):
        key_previous, value_previous = items_previous[i]
        key_next, value_next = items_next[i]
       
        first_letters = ""  # Initialize an empty string to store the first letters
        first_letters = None
        seq =[]
       
        first_letters = ''.join(data[key][0][0] for key in range(key_previous+k, key_next))

        seq.append(first_letters)

    return seq
    
def find_inverted_sequences(k, readdictt, previous_dict, next_dict, chromosomes, chr_number):
    """
    Given some parameters, this function attempts to identify and return inverted sequences.
    
    Args:
    - k (int): Presumably the k-mer length.
    - readdictt (dict): Dictionary related to sequence reads.
    - previous_dict (dict): Dictionary representing previous sequences or nodes.
    - next_dict (dict): Dictionary representing next sequences or nodes.
    - chromosomes (dict): Dictionary containing chromosome data.
    - chr_number (int): The chromosome number to consider.
    
    Returns:
    - dict: Dictionary containing the inverted sequences found.
    """
    
    undefined_seq = undefined_dicts(k, readdictt, previous_dict, next_dict)
    invert_dict = {}
    
    items_previous = list(previous_dict.items())
    items_next = list(next_dict.items())
    
    for i in range(len(items_previous)):
        key_previous, value_previous = items_previous[i]
        key_next, value_next = items_next[i]

        reversed_sequence = undefined_seq[i][::-1]
        mat = chromosomes[list(chromosomes.keys())[chr_number-1]].find_sequence_location(reversed_sequence, lennode=k-1)[0]
        
        if len(mat) >= (len(reversed_sequence)-k+2)*0.8:
            result = 'invert'+str(value_previous[0]+k)
            count = value_next[0]-1-(value_previous[0]+k)+1
            invert_dict[result] = (chr_number,value_previous[0]+k, value_next[0]-1,count)
    
    return invert_dict

# tra
def find_tra_keys(k,original_dict):
    previous_discontinuous_dict = {}
    next_discontinuous_dict = {}

    keys = list(original_dict.keys())
    values = list(original_dict.values())

    for i in range(len(keys)):
        key = keys[i]
        value = values[i]

        if i > 0 and key -(keys[i - 1] + 1)>1.5*k :
            # if value[0] == values[i - 1][0] +k:
            next_discontinuous_dict[key] = value

        if i < len(keys) - 1 and (keys[i + 1] - 1)-key>1.5*k :
            # if value[0] == values[i + 1][0] -k:
            previous_discontinuous_dict[key] = value

    return previous_discontinuous_dict, next_discontinuous_dict


