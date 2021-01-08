import os
import random
import string
#from termcolor import colored
import subprocess
import pandas as pd

cores = [4,6]
#cores = [4]
#word_len = [10**5, 10*6, 10**7, 10**8, 10**9]
word_len = [10**4, 10**5]
#WORD_LEN = 10**6
alphabet = ["A","C","T","G"]
k = 32
perf_data = list()


#alphabet = string.ascii_uppercase
for core in cores:
    for length in word_len:
        for repetition in range(3):

            input_string = ''.join(random.choices(alphabet, k=length))

            with open("data/input.txt", "w") as file:
                file.write(input_string)

            process = subprocess.Popen('mpiexec -n ' + str(core) + ' ./main', shell=True, stdout=subprocess.PIPE)
            process.wait()

            timed = process.communicate()[0]

            perf_data.append([core, length,float(timed.decode().split('*')[1])])



        df = pd.DataFrame(perf_data, columns = ['nr_cores', 'str_length', 'time'])
        df.to_csv("performance.csv")
