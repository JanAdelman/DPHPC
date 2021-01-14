import os
import random
import string
#from termcolor import colored
import subprocess

WORD_LEN = 10**6
 
#alphabet = string.ascii_uppercase
alphabet = ["A","C","T","G"]
input_string = ''.join(random.choices(alphabet, k=WORD_LEN))

with open("./data/input.txt", "w") as file:
    file.write(input_string)
