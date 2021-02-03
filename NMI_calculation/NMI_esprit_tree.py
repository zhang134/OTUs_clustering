# calculate the NMI value from the clustering results of uclust. By Ze-Gang Wei.

import os,sys,getopt

if __name__ == "__main__":
        usage = """

usage: python NMI_esprit_tree.py --input <otu file>  -r <reference NMI label file> -s <sequence file>

Example: python NMI_esprit_tree.py -i esprit_out..Clusters -r tax_assignments_NMI_label.txt -s sequence.fa
----------------------- All file path needs absolute path.
--input/-i      (required) The OTU clustered output file by ESPRIT-Tree (*.Clusters).
--ref/-r        (required) The annotated file with two columns: 1) sequence header and 2) genus.
--seq/-s	(required) The sequence file.
--help/-h       Help

"""
	opts,arg=getopt.getopt(sys.argv[1:],"i:r:s:h:",['input=','ref=','seq=','help='],)
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
                        otu_file = a

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
	
	sequence = []
	qq = open(sequence_file)
	ss = qq.readline()
	while ss:
		if ss.startswith('>'):
			ss = ss.strip()
			ss = ss[1:]
			sequence.append(ss)
		ss = qq.readline()
	qq.close()
	esprittree = {}
	#otu_num = 0
	oo = open(otu_file)
	row = oo.readline().strip()
	cutoffss = []
	while row:
		#otu_num = 0 
		cluster = {}
		row_list = row.split('|')
		cutoff = row_list[0].strip()
		cutoffss.append(cutoff)
		for x in range(1,len(row_list)):
			otu_id = x
			otu_is = row_list[x]
			for index in otu_is.split():
				seq = sequence[int(index)]
				cluster[seq] = str(otu_id)
		esprittree[cutoff] = cluster
		row = oo.readline().strip()
	oo.close()

#	nmi_file = otu_file.split('.')[0:-1]
	nmi_all = {}
	otu_all = {}
	for cutoff in cutoffss:
		group = esprittree[cutoff]
		nmi_file = otu_file.split('.')[0:-1]
		path = otu_file.split('.')[0:-1]
		nmi_file = '_'.join(nmi_file) + "_" + cutoff + "_nmi.txt"
		f = open(nmi_file, 'w')
		for i,j in group.items():
			if label.get(i, "False") == 'False':
				continue
			genus = label[i]
			f.write(i + '\t' + genus + '\t' + j + '\n')
		f.close()
		print cutoff
		otu_all[cutoff] = len(set(group.values()))
		#print "---------------------------------------------------------------------"
		#print "\nOTU number by ESPRIT-Tree: ", len(set(group.values()))
		#print "True OTU number is  : ", len(set(label.values()))
		cmd = "python ~/mywork/my_script/NMI_compute/NMI_calculation.py -i " + nmi_file + '>' +  '_'.join(path) + '_NMI.log'
		os.system(cmd)
		#print "output file:"
		#print nmi_file
		hh = open('_'.join(path) + '_NMI.log')
		#print '_'.join(path) + '_NMI.log'
		all_lines = hh.readlines()
		hh.close()
		nmi_all[cutoff] = all_lines[1].split()[2]
	print '\n\nNMI number from 0.01 to 0.10:'
	haha = 0
	for oo in cutoffss:
		haha += 1
		if haha%2 == 1 and haha < 21:
	#print 'OTU number: from 0.01 to 0.15':
			print nmi_all[oo],
	haha = 0
	print '\n\nOTU value from 0.01 to 0.10:'
	for oo in cutoffss:
                #print 'OTU number: from 0.01 to 0.15':
		haha += 1
		if haha%2 == 1 and haha < 21:
			print otu_all[oo],
