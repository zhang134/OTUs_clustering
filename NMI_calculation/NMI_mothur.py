# calculate the NMI value from the clustering results of mothur. By Ze-Gang Wei.

import os,sys,getopt

if __name__ == "__main__":
        usage = """

usage: python NMI_mothur.py --input <otu file>  -r <reference NMI label file>

Example: python NMI_mothur.py -i mix_merge_clean_otus.txt -r tax_assignments_NMI_label.txt
----------------------- All file path needs absolute path.
--input/-i      (required) The OTU clustered output file by mothur program.
--ref/-r         (required) The annotated file with two columns: 1) sequence header and 2) genus.
--help/-h       Help

"""
	opts,arg=getopt.getopt(sys.argv[1:],"i:r:h:",['input=','ref=','help='],)
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

        for i,a in opts:
                 if i in ("--input","-i"):
                        if not os.path.isfile(a):
                                print "***%s is not found.***"%(a)
                                print usage
                                sys.exit(1)
                        otu_file = a

                 if i in ("--ref","-r"):
                        if not os.path.isfile(a):
                                print "***%s is not found.***"%(a)
                                print usage
                                sys.exit(1)
                        label_file = a


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
	otu_num = 0
	cluster = {}
	oo = open(otu_file)
	row = oo.readline().strip()
	while row:
		if not row.startswith('>'):
			row = oo.readline().strip()
			continue
		otu_num += 1 
		row_list = row.split()
		otu_id = row_list[1]
		seq = row_list[0]
		seq = seq[1:]
		cluster[seq] = otu_id
		row = oo.readline().strip()
	oo.close()

	nmi_file = otu_file.split('.')[0:-1]
	nmi_file = '_'.join(nmi_file) + "_nmi.txt"
	f = open(nmi_file, 'w')
	for i,j in cluster.items():
		if label.get(i, "False") == 'False':
			continue
		genus = label[i]
		f.write(i + '\t' + genus + '\t' + j + '\n')
	f.close()
	print "\nOTU number by mothur: ", len(set(cluster.values()))
	print "True OTU number is  : ", len(set(label.values()))
	#cmd = "python /home/disk1/weizegang/mywork/my_script/NMI_compute/NMI_calculation.py -i " + nmi_file
        cmd = "python2 /home/disk2/weizegang/mywork/my_script/NMI_compute/NMI_calculation.py -i " + nmi_file
	os.system(cmd)
	print "\noutput file:"
	print nmi_file
		





