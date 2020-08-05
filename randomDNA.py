import argparse
import os 
import sys
import random
my_parser = argparse.ArgumentParser(description='local alignment params')
my_parser.add_argument('-num', metavar='num_of_seq',
                       type=int,
                       help='number of random sequences to generate')
my_parser.add_argument('-len', metavar='length_of_seq', type=int,
                       help='length of random generated seq')
args = my_parser.parse_args()
#information parsed in by user
seqL=args.len
seqC= args.num
#used to calculate final nucleotide frequencies 
totalCount= seqC*seqL
def printLine(a, l, randfile):
    i=0
    randfile.write(">\n")
    while i< len(a[l]):
        randfile.write(a[l][i])
        i+=1
    randfile.write("\n")

def randDNAGenerator():
    randfile= open("randSeq.txt", "w")
    #count of each nucleotide
    Acount=0
    Gcount=0
    Tcount=0
    Ccount=0
    #counts outer index of list
    currSeqeunce=0
    #counts inner index of list
    currIndex=0
    list=[['']*seqL]*seqC
    while currSeqeunce<seqC:
        #choose random 0-4 num, determines Nucletotide
        rand=random.randint(0,3)
        #append random nucleotide to the list 
        if rand==0:
            Acount+=1
            list[currSeqeunce][currIndex]='A'
        elif rand==1:
            Ccount+=1
            list[currSeqeunce][currIndex]='C'

        elif rand==2:
            Gcount+=1
            list[currSeqeunce][currIndex]='G'

        else:  
            Tcount+=1
            list[currSeqeunce][currIndex]='T'

        currIndex+=1
        if currIndex>=seqL:
            printLine(list, currSeqeunce, randfile)
            currSeqeunce+=1
            currIndex=0
    #calculate freq of each nucleotide and output
    Afreq= float(Acount)/totalCount
    Cfreq= float(Ccount)/totalCount
    Gfreq= float(Gcount)/totalCount
    Tfreq= float(Tcount)/totalCount
    randfile.write(">frequency of A's: "+str(Afreq)+"\n")
    randfile.write(">frequency of C's: "+str(Cfreq)+"\n")
    randfile.write(">frequency of G's: "+str(Gfreq)+"\n")
    randfile.write(">frequency of T's: "+str(Tfreq)+"\n")
    randfile.close()
randDNAGenerator()
