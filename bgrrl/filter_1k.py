import sys

from ktio.ktio import readFasta

if len(sys.argv) > 2:
    minlen = int(sys.argv[2])
else:
    minlen = 1000

# >1 length=630452 depth=0.94x
for _id,_seq in readFasta(sys.argv[1]):
    if len(_seq) >= minlen:
        print(_id, _seq, sep="\n")
