# calculate the ARI value from the clustering results. By Ze-Gang Wei.
# You just need to provide the right labels and clustered labels.

# Notice: the label you clustered is the index of the sequence, So you also need to provide the sequence file

import os,sys,getopt
from sklearn import metrics

if __name__ == "__main__":
        usage = """

usage: python ARI_compute_index.py -s <sequence file> -i <clustered labels file>  -r <right labels file> 

Example: python ARI_compute_index.py -s seq.fa -i labels_clustered.txt -r labels_right.txt
----------------------- All file path needs absolute path.
--input/-i      (required) The labels file clustered by one method.
--ref/-r        (required) The right labels file with two columns: 1) sequence header and 2) genus.
--seq/-s        (required) The sequence file.
--help/-h       Help

"""
        opts,arg=getopt.getopt(sys.argv[1:],"i:r:s:h:",['input=', 'ref=', 'seq=', 'help='],)
	parameters=[a[0] for a in opts]
	if '-h' in parameters or '--help' in parameters:
		print usage
                sys.exit(1)
        if len(parameters)==0:
                print usage
                sys.exit(1)
        if '-i' not in parameters and '--input' not in parameters:
                print "***Error, a input file is requred.***\n"
                print usage
                sys.exit(1)
	if '-r' not in parameters and '--ref' not in parameters:
                print "***Error, the annotated input file is requred.***\n"
                print usage
                sys.exit(1)
	if '-s' not in parameters and '--seq' not in parameters:
                print "***Error, the sequence file is requred.***\n"
                print usage
                sys.exit(1)


        for i,a in opts:
                 if i in ("--input","-i"):
                        if not os.path.isfile(a):
                                print "***%s is not found.***"%(a)
                                print usage
                                sys.exit(1)
                        your_label_file = a

                 if i in ("--ref","-r"):
                        if not os.path.isfile(a):
                                print "***%s is not found.***"%(a)
                                print usage
                                sys.exit(1)
                        label_file = a
	         if i in ("--seq","-s"):
			if not os.path.isfile(a):
				print "***%s is not found.***"%(a)
				print usage
				sys.exit(1)
			sequence_file = a

        seqs = []
        ss = open(sequence_file)
        line = ss.readline().strip()
        while line:
            if line.startswith('>'):
                seqs.append(line[1:])
            line = ss.readline().strip()
        ss.close()
	label_file
	label = {}
	ff = open(label_file)
	line = ff.readline().strip()
	while line:
		seq_header = line.split()[0]
		seq_genus  = line.split()[1]
		label[seq_header] = seq_genus
		line = ff.readline().strip()
	ff.close()
	file_name = os.path.basename(your_label_file)
	oo = open(your_label_file)
	row = oo.readline().strip()
	label_predicted = {}
	while row:
		seq_header = int(row.split()[0])
                seq_genus  = row.split()[-1]
                label_predicted[seqs[seq_header - 1]] = seq_genus
		row = oo.readline().strip()
	oo.close()

	label_right = []
	label_you = []
	for i,j in label.items():
	    label_right.append(j)
            label_you.append(label_predicted[i])
        ARI = metrics.adjusted_rand_score(label_right, label_you)
        print '\n\n ARI is: ', ARI
