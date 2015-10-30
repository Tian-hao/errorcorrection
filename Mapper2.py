#!/usr/bin/python
import os
import sys
import glob
from collections import Counter

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def errorcorrect(H):
  EFreads = []
  [EFreads.append(H[bc][0]) for bc in H.keys() if len(set(H[bc])) == 1 and len(H[bc]) >= 3]
  return EFreads

def hashingread(infile):
  H = {}
  infile = open(infile,'r')
  countread = 0
  for line in infile.xreadlines():
    countread += 1
    if countread%5000000==0: print 'hasing %d reads' % countread
    line = line.rstrip().rsplit("\t")
    bc   = line[0]
    read = line[1]
    if H.has_key(bc): H[bc].append(read)
    else: H[bc] = [read]
  infile.close()
  return H

def outputing(H,outfile):
  outfile  = open(outfile,'w')
  for read in H.keys():
    outfile.write(read+"\t"+str(H[read])+"\n")
  outfile.close()

def readin(infile):
  EFreads = []
  infile = open(infile,'r')
  countread = 0
  for line in infile.xreadlines():
    countread += 1
    if countread%5000000==0: print 'hasing %d reads' % countread
    line = line.rstrip().rsplit("\t")
    bc   = line[0]
    read = line[1]
    EFreads.append(read)
  infile.close()
  return EFreads

def main():
  filenames = sorted(glob.glob('./paired/*'))
  for infile in filenames:
    print 'working on %s' % infile
    readhash = {}
    EFreads  = []
    readhash = hashingread(infile)
    EFreads  = errorcorrect(readhash)
    readhash = Counter(EFreads); del EFreads
    DNAfile  = infile.replace('paired/','errorcorrect/')+'Tag.dna'
    outputing(readhash,DNAfile)

    readhash = {}
    EFreads  = []
    EFreads  = readin(infile)
    readhash = Counter(EFreads); del EFreads
    DNAfile  = infile.replace('paired/','errorcorrect/')+'.dna'
    outputing(readhash,DNAfile)

if __name__ == "__main__":
  main()

