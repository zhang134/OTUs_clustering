# calculate the NMI value from the clustering results of CD-HIT for all cufoffs. By Ze-Gang Wei.

import os,sys,getopt

if __name__ == "__main__":
        usage = """

usage: python NMI_cdhit_all_cutoff.py --input <otu directory>  -r <reference NMI label file> -s <sequence file>

Example: python NMI_cdhit_all_cutoff.py -i ~/cdhit/ -r tax_assignments_NMI_label.txt -s sequence.fa
----------------------- All file path needs absolute path.
--input/-i      (required) The OTU clustered output file directory by CD-HIT, need absolute path.
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
                        if not os.path.exists(a):
                                print "***%s directory not found.***"%(a)
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

	all_file_names = os.listdir(otu_file)
	key_words = 'cdhit'
	cdhit = {}
	cutoffss = []
	cutoff_all = ['cdhit99','cdhit98','cdhit97','cdhit96','cdhit95','cdhit94','cdhit93','cdhit92','cdhit91','cdhit90',\
		      'cdhit89','cdhit88','cdhit87','cdhit86','cdhit85']
	cutoff_files = {}
	for cc in cutoff_all:
		for ff in all_file_names:
			if cc in ff:
				cutoff = cc[5:]
				ff_is = otu_file + ff
				oo = open(ff_is)
				row = oo.readline().strip()
				otu_id = 0
				cluster = {}
				cutoffss.append(cutoff)
				while row:
					if row.startswith('>Cluster'):
						otu_id += 1 
					row_list = row.split('...')[0]
					header = row_list.split('>')[1]
					cluster[header] = str(otu_id)
					row = oo.readline().strip()
				cdhit[cutoff] = cluster
				cutoff_files[cutoff] = ff_is
				oo.close()
				break
	nmi_all = {}
	otu_all = {}
	for cutoff in cutoffss:
		group = cdhit[cutoff]
		file_name = cutoff_files[cutoff]
		nmi_file = file_name.split('.')[0:-1]
		#print 'nmi_file = ', nmi_file
		#path = otu_file
		#print 'path = ', path
		nmi_file = "_".join(nmi_file) + '_' + cutoff + "_nmi.txt"
		print 'cutoff = ',cutoff
		f = open(nmi_file, 'w')
		for i,j in group.items():
			if label.get(i, "False") == 'False':
				continue
			genus = label[i]
			f.write(i + '\t' + genus + '\t' + j + '\n')
		f.close()
		#print cutoff
		otu_all[cutoff] = len(set(group.values()))
		#print "---------------------------------------------------------------------"
		#print "\nOTU number by ESPRIT-Tree: ", len(set(group.values()))
		#print "True OTU number is  : ", len(set(label.values()))
		cmd = "python ~/mywork/my_script/NMI_compute/NMI_calculation.py -i " + nmi_file + '>' +  "_".join(nmi_file.split('.')[0:-1]) + '_NMI.log'
		os.system(cmd)
		#print "output file:"
		#print nmi_file
		hh = open("_".join(nmi_file.split('.')[0:-1]) + '_NMI.log')
		#print '_'.join(path) + '_NMI.log'
		all_lines = hh.readlines()
		hh.close()
		nmi_all[cutoff] = all_lines[1].split()[2]
	print '\n\nNMI number from 0.01 to 0.15:'
	for oo in cutoffss:
		#print 'OTU number: from 0.01 to 0.15':
		print nmi_all[oo],
	print '\n\nOTU value from 0.01 to 0.15:'
	for oo in cutoffss:
                #print 'OTU number: from 0.01 to 0.15':
                print otu_all[oo],
