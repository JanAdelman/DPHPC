import os
import random
import string
import subprocess
import pandas as pd

cores = [4,6]
word_len = [10**5, 10*6, 10**7, 10**8, 10**9]
alphabet = ["A","C","T","G"]
k = 32
perf_data = list()

#alphabet = string.ascii_uppercase
for core in cores:
    for length in word_len:
        for repetition in range(2):

            input_string = ''.join(random.choices(alphabet, k=length))

            with open("data/input.txt", "w") as file:
                file.write(input_string)

            process = subprocess.Popen('mpiexec -n ' + str(core) + ' ./main', shell=True, stdout=subprocess.PIPE)
            process.wait()

            timed = process.communicate()[0]
            overall_time = float(timed.decode().split('*')[1])
            min_time = float(timed.decode().split('*')[3])
            max_time = float(timed.decode().split('*')[5])
            average_time = float(timed.decode().split('*')[7])
   
            perf_data.append([core, length, overall_time, min_time, max_time, average_time])

        df = pd.DataFrame(perf_data, columns = ['nr_cores', 'str_length', 'overall_time', 'min_time', 'max_time', 'average_time'])
        df.to_csv("performance.csv")
