import random

# 定义碱基
bases = ['A', 'T', 'G', 'C']
# 生成随机序列
random_sequence = ''.join(random.choice(bases) for _ in range(50))
print(random_sequence)
def insert(chromosome,Result = '',result='' ):
    bases = ['A', 'T', 'G', 'C']
    random_sequence = ''.join(random.choice(bases) for _ in range(150))
    random_numbers = [i for i in range(20000000, 20010000, random.randint(1900, 2000))]
    for index in random_numbers:
            oneinsert =str(index)+'insert'
            chromosome =chromosome
            chromosome_start = index
            randomseed =1
            
            list = [oneinsert,chromosome,index
                    ,randomseed]
            converted_list = [str(item) for item in list]
            result = ','.join(converted_list)
            Result = Result + ":"+result
            print(index)
    return  Result[1:]


def tandem_repeat(chromosome,Result = '',result=''):
    random_numbers = [i for i in range(20000000, 20010000, random.randint(1900, 3000))]
    for index in random_numbers:
            onedelete = str(index)+'tandem_repeat'
            chromosome = chromosome
            chromosome_start = index
            chromosome_end = index+random.randint(35, 55)
            randomseed = 1
            m = random.randint(1, 3)
            list= [onedelete,chromosome,chromosome_start
                    ,chromosome_end,m,randomseed]
            converted_list = [str(item) for item in list]
            result = ','.join(converted_list)
            Result = Result + ":"+result
    return  Result[1:]

command_part1 = 'python2 -B Wangshj.py -out /home/cloudam/my_folder_graph -depth 100 '
random_sequence = str(random_sequence)
tandem_repeat =  '-tandem_repeat'+' '+tandem_repeat(chromosome='chr1') 
insert_function_call = '-insert ' + '' + str(insert(chromosome='chr1'))
# tandem_repeat_call = 
bed = '-bed_file mybed.bed'

# 使用空格进行字符串连接
full_command = command_part1 +  ' ' + tandem_repeat
print(full_command)


