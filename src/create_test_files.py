import os
import random
import string
import subprocess

k = 32

for pow in [6,7,8,9]:
    
    WORD_LEN = 10**pow

    print("Generating Sequence...")
    #alphabet = string.ascii_uppercase
    alphabet = ["A","C","T","G"]
    input_string = ''.join(random.choices(alphabet, k=WORD_LEN))
    print("...Generated Sequence")

    print("Creating Kmers...")
    kmers = dict()
    for i in range(len(input_string)-k+1):
        kmers[i] = input_string[i:i+k]
    print("..Created Kmers")

    print("Sorting...")
    expected_output = list({k: v for k, v in sorted(kmers.items(), key=lambda item: item[1])}.keys())
    expected_output = ','.join(map(str, expected_output)) + ','
    print("...Sorted")


    print("Writing Files...")
    with open("expected" + str(pow) + ".txt", "w") as file:
        file.write(expected_output)

    with open("input" + str(pow) + ".txt", "w") as file:
        file.write(input_string)
    print("...done!")

