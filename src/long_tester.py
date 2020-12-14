import os
import random
import string
from termcolor import colored
import subprocess

WORD_LEN = 10**3
 
#alphabet = string.ascii_uppercase
alphabet = ["A","C","T","G"]
input_string = ''.join(random.choices(alphabet, k=WORD_LEN))

k = 32

kmers = dict()
for i in range(len(input_string)-k+1):
    kmers[i] = input_string[i:i+k]

print("Python Sorting")
expected_output = list({k: v for k, v in sorted(kmers.items(), key=lambda item: item[1])}.keys())
expected_output = ','.join(map(str, expected_output)) + ','

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

if (output == expected_output):
    color = "green"
else:
    color = "red"

print("-" * 20)
print("MPI: " + output[:20])
print("EXP: " + expected_output[:20])
print(colored(output == expected_output, color))

