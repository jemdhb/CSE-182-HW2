import argparse
import os 
import sys
import numpy as np
import matplotlib.pyplot as plt 
#arg parser
my_parser = argparse.ArgumentParser(description='local alignment params')
my_parser.add_argument('seq_file', metavar='seq_file',
                       type=str,
                       help='the path to the file w/ the seqs')
args = my_parser.parse_args()
seqFile=open(args.seq_file, 'r')
localQ1= open("localQ1.txt", "w")
averages= open("average.txt", "w")
#reads files with more than 2 sequences 
def readFilePairs(file, pairnum):
    seqFile= open(file, "r")
    list=['']*pairnum*2
    with seqFile as infile:
        i=0
        j=0
        for line in infile:
            if line.startswith(">"):
                continue
            else:
                #alternate which string to read
                #based on the modulus
                if i%2==1:
                    x=line[0: len(line)-1]
                    list[i]+=x
                    i+=1
                else:
                    y=line[0: len(line)-1]
                    list[i]+= y
                    i+=1
                #done for then done read non sequence data at 
                # the bottom of a file
                if i>pairnum*2: break

    return list
#modified local for plotting purposes
def locAL(x,y, match, mismatch, indel):
    #used for scores
    A= np.zeros((len(x) + 1, len(y) + 1), np.int)
    #used for directionality 
    strings = np.empty((len(x) + 1, len(y) + 1), dtype='string')
    #best score and location of the best score
    best = 0
    optloc = [0,0]
    # fill in A with scores and strings with directions 
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            score= A[i,j]
            A[i][j], strings[i,j] = max(
            (A[i][j-1] + indel, "i"), 
            (A[i-1][j] + indel, "d"),
            (A[i-1][j-1] + (match if x[i-1] == y[j-1] else mismatch),"m"), (0,"x"))
    # return the opt score and the best location
    best=np.max(A)
    i,j = np.unravel_index(A.argmax(), A.shape)
    optloc[0]=i
    optloc[1]=j
    #will have aligned sequences
    seq1=""
    seq2=""
    #local alignment 
    while (optloc[0]>0 and optloc[1]>0 ):
        #if match
        xx= optloc[0]
        yy=optloc[1]
        #exit if zero score
        if A[xx,yy]==0:
            break
        #diagonal
        if strings[xx][yy]=="m":
            seq1=x[xx-1]+seq1
            seq2=y[yy-1]+seq2
            optloc[0] -= 1
            optloc[1] -= 1
        #vertical movement
        elif strings[xx][yy]=="i":
            seq1="-"+seq1
            seq2=y[yy-1]+seq2
            optloc[1] -= 1
        #horizontal movement 
        elif strings[xx][yy] == "d":
            seq1=x[xx-1]+seq1
            seq2='-'+seq2
            optloc[0] -= 1
    #outputs file in a readable manner
    ii=60
    jj=60
    while ii<len(seq1) and jj<len(seq2):
        localQ1.write(seq1[ii-60:ii]+"\n")
        localQ1.write(seq2[jj-60:jj]+"\n\n")
        ii+=60
        jj+=60
    localQ1.write(seq1[ii-60:len(seq1)]+"\n")
    localQ1.write(seq2[jj-60:len(seq2)]+"\n")
    localQ1.write("high score of: "+str(best)+" with length: "+ str(len(seq1))+"\n")
    return best, len(seq1)
#used to plot average alignment lengths 
def trackavg(matchvals, mismatchvals):
    #read previously calculated pairs
    l=readFilePairs("randseq.txt",52)
    #while loop tracking variables
    k=-1
    m=0
    average=0
    lengths=[]*len(matchvals)
    #run outer until local has been run on all pairs 
    #with all match values 
    while k< (len(matchvals)-1):
        lengths.append(average)
        k+=1
        print("k: "+ str(k))
        num=[0]*len(matchvals)
        avgL=0
        i=0
        j=0
        #run local until all pairs have been read 
        while i+1<len(l)/2 and j < len(l)/2:
            p=locAL(l[i],l[i+1],matchvals[k],mismatchvals[k],mismatchvals[k])
            i+=2
            #grab average value from local results 
            avgL+=p[1]
            j+=1
            #calc avg
            average= float(avgL)/j
        
    # heights of bar--> avg sequence length 
    height = lengths
  
    # plotting bar chart 
    plt.bar(mismatchvals, height,
        width = 0.8, color ='green') 
  
    # x-axis name
    plt.xlabel('Mismatch = Indel') 
    # y-axis name
    plt.ylabel('Average alignment length') 
    # plot title 
    plt.title('Average length with different Scoring values') 
  
    # function to show the plot 
    plt.show() 
#plots 3a
matchvals=[8,8,8,8,8,8,8,8]
mismatchvals=[0,-.25,-.33,-.5,-1,-10,-20,-30]
trackavg(matchvals, mismatchvals)
#plots 3b
matchvals=[8,8,8,8,8]
mismatchvals=[-10,-15,-20,-25,-30]
trackavg(matchvals, mismatchvals)
