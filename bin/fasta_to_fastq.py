#!/usr/bin/python

import sys

FILE=sys.argv[1]

header = False
with open("{}".format(FILE), "r") as f :
    for line in f:
        if line.startswith(">") :
            if header : 

                print(f"@{header}\n{sequence}\n+\n{qualscore}")

                header = line.strip().replace(">", "")
            else : 
                header = line.strip().replace(">", "")
        else : 
            sequence = line.strip()
            qualscore = "I"*len(sequence)

f.close()