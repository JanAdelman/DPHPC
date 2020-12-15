import os
import random
import string
#from termcolor import colored
import subprocess
import pandas as pd

#cores = [4,8,16,32,64]
cores = [4]
#word_len = [10**5, 10*6, 10**7, 10**8, 10**9]
word_len = [10**5, 10**4]
#WORD_LEN = 10**6
alphabet = ["A","C","T","G"]
k = 32
perf_data = list()


#alphabet = string.ascii_uppercase
for core in cores:
    for length in word_len:
        for repetition in range(0,1):

            input_string = ''.join(random.choices(alphabet, k=length))
            """
            kmers = dict()
            for i in range(len(input_string)-k+1):
                kmers[i] = input_string[i:i+k]

            expected_output = list({k: v for k, v in sorted(kmers.items(), key=lambda item: item[1])}.keys())
            expected_output = ','.join(map(str, expected_output)) + ','

            with open("expected.txt", "w") as file:
                file.write(expected_output)
            """

            with open("input.txt", "w") as file:
                file.write(input_string)


            #print("MPI Running")
            process = subprocess.Popen('mpiexec -n ' + str(core) + ' --oversubscribe'+' ./main', shell=True, stdout=subprocess.PIPE)
            process.wait()

            timed = process.communicate()[0]

            #print("MPI done!")
            """
            output = str()
            with open("result.txt", "r") as file:
                output = file.read()
            """

            perf_data.append([core, length,float(timed.decode().split('*')[1])])

            #print('Dictionary: ', perf_data)

            #print(output == expected_output)

df = pd.DataFrame(perf_data, columns = ['nr_cores', 'str_length', 'time'])
df.to_csv("performance.csv")
