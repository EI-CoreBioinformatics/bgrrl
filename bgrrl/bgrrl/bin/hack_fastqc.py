import re
import os
import sys

"""
This module was written to counteract mis-classified Phred-encoding by FastQC on low-readcount samples.
"""

def hackFastQCData(_in, out=sys.stdout, allow_editing=True):
    hack_encoding = False
    mode = ""
    for line in _in:
        if line.startswith("Encoding"):
            line = line.strip().split("\t")
            encoding = line[1]
            # print(encoding)
            if encoding.startswith("Illumina 1.5"):
                # print("hack_encoding=", allow_editing)
                hack_encoding = allow_editing       
                if hack_encoding:
                    line = [line[0], "Illumina 1.9"]
            line = "\t".join(line)
        elif line.startswith(">>Per sequence quality scores") and hack_encoding:
            # print(line)
            mode = "qscore"
            # print("mode=", mode)
        elif line.startswith(">>Per base sequence quality") and hack_encoding:
            # print(line)
            mode = "qseq"
            # print("mode=", mode)
        elif line.startswith(">>END_MODULE"):
            if mode == "qscore" or mode == "qseq":
                mode = ""
        elif mode == "qscore" and not line.startswith("#"):
            # print("WAS:", line)
            line = line.strip().split("\t") 
            line[0] = "{:d}".format(int(line[0]) + 31)
            line = "\t".join(line)
            # print("IS:", line)
        elif mode == "qseq" and not line.startswith("#"):
            # print("WAS:", line)
            line = line.strip().split("\t")
            line[1] = str(float(line[1]) + 31)
            line = "\t".join(line)
            # print("IS:", line)

        print(line.strip(), file=out)

def run(fn):

    with open(fn) as _in, open(fn + ".bak", "w") as _out:
        hackFastQCData(_in, _out, allow_editing=False)
    with open(fn + ".bak") as _in, open(fn, "w") as _out:
        hackFastQCData(_in, _out, allow_editing=True)

    return None 

def main(args):
 
    for cdir, dirs, files in os.walk(args[0]):
        for f in filter(lambda x:x == "fastqc_data.txt", files):
            print(cdir)
            run(os.path.join(cdir, f))
    


if __name__ == "__main__": 
    main(sys.argv[1:])
