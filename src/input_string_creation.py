import os
import random
import string
import csv

string_length = [10**5, 10**6, 10**7, 10**8, 10**9]
Len = [10_5, 10_6, 10_7, 10_8, 10_9]
alphabet = ["A","C","T","G"]

for i in range(0, len(string_length)):
	
	input = string_length[i]
	name = Len[i]
	input_string = ''.join(random.choices(alphabet, k=input))

	with open("./data/input_{}.txt".format(Len[i]), "w") as file:
		file.write(input_string)


    