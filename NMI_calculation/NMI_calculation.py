import math
import sys
import os
import getopt

def readList(listfile):
    goldDict={}
    testDict={}
    readsNum=0
    for line in open(listfile,'r'):
        readsNum+=1
        lineList=line.strip().split()
        if len(lineList) != 3:
            print "***There is an error in input file, at %d row.\nthree column that separated by tab is requred.***"
            sys.exit(1)
        seqname=lineList[0]
        specie=lineList[1]
        cluster=lineList[2]
        #insert to gold dictionary
        goldList=goldDict.get(specie,[])
        goldList.append(seqname)
        goldDict[specie]=goldList
        #insert test dictionary
        testList=testDict.get(cluster,[])
        testList.append(seqname)
        testDict[cluster]=testList
    return goldDict,testDict,readsNum

def NMIcalculater(listfile):
    goldDict,testDict,readsNum=readList(listfile)
    I_value=0
    Hs_value=0
    Hc_value=0

    N=float(readsNum) #number of all reads

    #H(C)
    for  cluster in testDict:
        clusterSeq=testDict[cluster]
        Cj=float(len(clusterSeq))
        sub_Hc=-((Cj/N)*math.log(Cj/N))
        Hc_value+=sub_Hc
    #H(S)
    for specie in goldDict:
        specieSeq=goldDict[specie]
        Si=float(len(specieSeq))
        sub_Hs=-Si/N*math.log(Si/N)
        Hs_value+=sub_Hs  #H(S)
    #I(S,C)
    for specie in goldDict:
        specieSeq=goldDict[specie]
        Si=float(len(specieSeq))
        for cluster in testDict:
            clusterSeq=testDict[cluster]
            intersection=list(set(specieSeq).intersection(set(clusterSeq)))
            if len(intersection)==0:
                continue
            Aij=float(len(intersection))
            Cj=float(len(clusterSeq))
            subI_log=math.log((Aij/N)/((Si*Cj)/(N*N)))
            sub_I=(Aij/N)*subI_log
            I_value+=sub_I
    NMI_value=(2*I_value)/(Hs_value+Hc_value)
    return NMI_value
if __name__=="__main__":
    usage="""usage: python NMI_caculation.py --input <inputfile>

 --input/-i	the list file (three column separated by tab) for all sequence names. The first column is the sequence names, the second column is the species, and third column is the cluster names.
 --help/-h	help

    """
    opts,arg=getopt.getopt(sys.argv[1:],"i:h:",['input=','help='],)
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
    #
    for i,a in opts:
        if i in ("--input","-i"):
            if not os.path.isfile(a):
                print "***%s is not found.***"%(a)
                print usage
                sys.exit(1)
            inputfile=a
    nmi=NMIcalculater(inputfile)
    print "\nNMI value: %f\n"%(nmi)
