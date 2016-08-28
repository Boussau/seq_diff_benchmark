
#! /usr/bin/env python
#computes pairwise distances between sequences in a fasta file

import sys


# Reading a fasta file
def readFasta (file):
    try:
        f=open(file, 'r')
    except IOError:
        print ("Unknown file: ",file)
        sys.exit()
    seqs = list()
    nams = list()
    seq=""
    for l in f:
        if l.startswith('>'):
            if seq != "" and seqname != "":
                seqs.append(seq)
                nams.append(seqname)
            seq=""
            liste=l.split(":")
            if (liste[1].strip() != "0"):
                seqname = liste[0].replace(">","")
            else :
                seqname = ""
        else:
            seq=seq+l.strip()
    # last sequence
    if seq != "" and seqname != "":
        seqs.append(seq)
        nams.append(seqname)
    print ("Number of non-0 haplotypes: "+str(len(seqs)))
    f.close()
    return (nams, seqs )


# Computing the pairwise differences between two sequences
# Here zip is much faster than the loop programmed by hand.
def computePairwiseDifferences( char* x, char* y ):
   cdef int dif = 0
   cdef char xk
   cdef char yk
   for k in range(len(x) ):
      xk = x[k]
      yk = y[k]
      if xk != yk and xk != 'N' and yk != 'N':
         dif += 1
   return dif
    #dif = sum(1 for xk,yk in zip(x, y) if xk != yk and xk != "N" and yk != "N")


def computeAndOutputPairwiseDifferences(out, nams, seqs):
    try:
        fout=open(out, 'w')
    except IOError:
        print ("Unknown file: ",out)
        sys.exit()

    fout.write("from\tto\tdist\n")
    for seqi in range(len(seqs)-1):
        namei = nams[seqi]
        for seqj in range(seqi+1,len(seqs)):
            dif = computePairwiseDifferences( seqs[seqi], seqs[seqj] )
            fout.write(namei+"\t"+nams[seqj]+"\t"+str(dif)+"\n")
    #    if ( (seqi % 10) == 0):
    #        print ("seqi: "+str(seqi))
    fout.close()


if __name__ == '__main__':

  argv=sys.argv[1:]
  file=argv[0]
  out=argv[1]

  (nams, seqs) = readFasta(file)
  computeAndOutputPairwiseDifferences(out, nams, seqs)
