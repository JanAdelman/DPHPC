import os
import random
import string
#from termcolor import colored
import subprocess

WORD_LEN = 10**8
 
#alphabet = string.ascii_uppercase
alphabet = ["A","C","T","G"]
input_string = ''.join(random.choices(alphabet, k=WORD_LEN))

k = 32
terminator="$"*(k-1)
input_string=input_string+terminator

with open("input.txt", "w") as file:
    file.write(input_string)
