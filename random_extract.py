# Random extract the sequences by given number. By Zegang Wei. Email: david_nwpu@163.com

import os,sys,getopt,random
import numpy as np
if __name__ == "__main__":
    usage = """

usage: python random_extract.py -i <sequence file> -n <subset number> -o <out file>

Example: python random_extract.py -i seq.fa -n 1000 -o seq_random.fa

--input/-i	(required) The sequence file.
--out/-o	(required) The output file.
--num/-n	(optional) subset number [default: all sequences].
--help/-h	Help

"""

    opts,arg=getopt.getopt(sys.argv[1:],"i:n:h:o:",['input=','num=','help=', 'out='])
    parameters=[a[0] for a in opts]
    if '-h' in parameters or '--help' in parameters:
        print(usage)
        sys.exit(1)
    if len(parameters)==0:
        print(usage)
        sys.exit(1)
    if '-i' not in parameters and '--input' not in parameters:
        print("***Error, a input file is requred.***\n")
        print(usage)
        sys.exit(1)
    if '-o' not in parameters and '--out' not in parameters:
        print("***Error, the out file is requred.***\n")
        print(usage)
        sys.exit(1)
    if '-n' not in parameters and '--num' not in parameters:
        print("***Error, subset numbe is requred.***\n")
        print(usage)
        sys.exit(1)
    num_subset = -1
    for i,a in opts:
        if i in ("--input","-i"):
            if not os.path.isfile(a):
                print("***%s is not found.***"%(a))
                print(usage)
                sys.exit(1)
            seq_file = a
        if i in ("--out","-o"):
            out_file = a
        if i in ("--num","-n"):
            num_subset = int(a)
    ff = open(seq_file)
    #ll = open(annoted_file)
    line = ff.readline()
    #labels = ll.readlines()
    #ll.close()
    num = -1
    heads = []
    seqs = []
    nn = 0
    seq_is = ''
    while line:
        if line.startswith('>'):
            line = line.strip()
            num += 1
            heads.append(line)
            if num > 0:
                seqs.append(seq_is)
                seq_is = ''
        else:
            line = line.strip()
            seq_is += line
        line = ff.readline()
    ff.close()
    seqs.append(seq_is)
    print('\n  Total sequences: %d' % len(seqs))
    if len(seqs) != len(heads):
        print("\nSequences number not equal header number! Program quit!")
        sys.exit(1)
    if num_subset == -1:
        num_subset = num

    random_list = np.arange(num_subset)
    random.shuffle(random_list)
    mm = -1
    out = open(out_file, 'w')
    for i in random_list:
        mm += 1
        out.write(heads[i] + "\n")
        out.write(seqs[i] + "\n")
        if mm + 1 == num_subset:
            break
    out.close()
    #print('\n  Total sequences: %d\n' % num)
    print("\n  Random sequence selected finished.\n")
