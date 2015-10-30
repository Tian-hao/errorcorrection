#!/usr/bin/env python
#For alignment
import os
import sys
import glob
import operator
import numpy
from multiprocessing import Pool, Process
from collections import Counter
from string import atof
from itertools import imap
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def hashinref(infile):
  refhash = {}
  handle = open(infile, "rU")
  for record in SeqIO.parse(handle, "fasta"):
    refhash[str(record.id)] = str(record.seq)
  handle.close()
  return refhash

def callindel(seqaln,refaln):
  I   = []
  D   = []
  pos = 0
  trimseq = ''
  for n in range(len(refaln)):
    if refaln[n] != '-': pos += 1; trimseq += seqaln[n]
    if refaln[n] == '-' and seqaln[n] != '-': I.append(pos) #Calling Insertion
    if refaln[n] != '-' and seqaln[n] == '-': D.append(pos) #Calling deletion
  I = Counter(I)
  D = Counter(D)
  Iprofile = []
  Dprofile = []
  [Iprofile.append(str(pos)+'I'+str(I[pos])) for pos in sorted(I.keys())]
  [Dprofile.append(str(pos)+'D'+str(D[pos])) for pos in sorted(D.keys())]
  if len(Iprofile) == 0: Iprofile = '-'
  else: Iprofile = "-".join(Iprofile)
  if len(Dprofile) == 0: Dprofile = '-'
  else: Dprofile = "-".join(Dprofile)
  return trimseq,Iprofile,Dprofile
    

def alignment(refseq,infile,outfile):
  print 'processing file: %s' % infile
  infile  = open(infile,'r')
  outfile = open(outfile,'w')
  totalcount = 0
  goodcount  = 0
  for line in infile.xreadlines():
    line   = line.rstrip().rsplit("\t")
    seq    = line[0]
    count  = int(line[1])
    algns  = pairwise2.align.globalms(seq,refseq,1,-1,-1,-0.5)
    dist   = hamming(seq,refseq)
    seqaln = algns[0][0]
    refaln = algns[0][1]
    trimseq, Iprofile, Dprofile = callindel(seqaln,refaln)
    totalcount += count
    if '-' not in seqaln: goodcount += count
    outfile.write("\t".join(map(str,[trimseq,count,dist,Iprofile,Dprofile]))+"\n")
  infile.close()
  outfile.close()

def multi_run_wrapper(args):
   return alignment(*args)

def main():
  reffile  = "../Fasta/SeqInfo.fa"
  refhash  = hashinref(reffile)
  refseq  = refhash['ProteinG_WT']
  infiles  = sorted(glob.glob('../errorcorrect/*.dna'))
  outfiles = [infile.replace('errorcorrect/','alignment/') for infile in infiles]
  Inputs   = [(refseq,infiles[n],outfiles[n]) for n in range(len(infiles))]
  pool   = Pool(processes=1)
  pool.map(multi_run_wrapper,Inputs)
  pool.close()
  pool.join()

if __name__ == "__main__":
    main()

