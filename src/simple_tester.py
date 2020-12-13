import os
import random
import string
from termcolor import colored

NUM_TRIES = 10
MIN_LEN_STR = 1000000
MAX_LEN_STR = 1000000
MIN_K = 32
MAX_K = 32

for i in range(NUM_TRIES):
    wordlen = random.randint(MIN_LEN_STR, MAX_LEN_STR)
    #alphabet = string.ascii_uppercase
    alphabet = ["A","C","T","G"]
    input_string = ''.join(random.choices(alphabet, k=wordlen))
    k = random.randint(MIN_K, min(MAX_K, wordlen))

    kmers = dict()
    for i in range(len(input_string)-k+1):
        kmers[i] = input_string[i:i+k]

    expected_output = list({k: v for k, v in sorted(kmers.items(), key=lambda item: item[1])}.keys())
    expected_output = ','.join(map(str, expected_output)) + ','


    stream = os.popen('mpiexec -n 4 ./main ' + input_string + " " + str(k))
    output = stream.read()

    if (output == expected_output):
        color = "green"
    else:
        color = "red"

    print("-" * 20)
    print("K = " + str(k) + " / LEN = " + str(wordlen))
    #print(input_string)
    print("MPI: " + output[:20])
    print("EXP: " + expected_output[:20])
    print(colored(output == expected_output, color))

