import os
import random
import string
#from termcolor import colored
import subprocess

WORD_LEN = 10**6
 
#alphabet = string.ascii_uppercase
alphabet = ["A","C","T","G"]
input_string = ''.join(random.choices(alphabet, k=WORD_LEN))
#input_string="ATATAACCATAAATATAC"
WORD_LEN = len(input_string)

k = 32
terminator="$"*(k-1)
input_string=input_string+terminator


kmers = dict()
for i in range(len(input_string)-k+1):
    kmers[i] = input_string[i:]

print("Python Sorting")
expected_output = list({k: v for k, v in sorted(kmers.items(), key=lambda item: item[1])}.keys())
expected_output = ','.join(map(str, expected_output)) + ','

with open("expected.txt", "w") as file:
    file.write(expected_output)

print("Python Sorted")



with open("input.txt", "w") as file:
    file.write(input_string)

print("MPI Running")
process = subprocess.Popen('make compile-run', shell=True, stdout=subprocess.PIPE)
process.wait()
print("MPI done!")

output = str()
with open("result.txt", "r") as file:
    output = file.read()



print("-" * 20)
print("MPI: " + output[:20] + " : " + str(len(output)))
print("EXP: " + expected_output[:20] +  " : "  + str(len(expected_output)))
print(output == expected_output)


