#!/usr/bin/env python

import os
import sys
import glob
import operator
from string import atof
from itertools import imap
from Bio import SeqIO

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

def errorcorrect(H):
  EFreads = {}
  for bc in H.keys():
    if len(set(H[bc])) == 1 and len(H[bc]) >= 3:
      EFreads[bc] = len(H[bc])
  return EFreads

def main():
  outfile  = open('ec34.cov','w')
  filenames = sorted(glob.glob('../paired/*'))
  for infile in filenames:
    print 'working on %s' % infile
    readhash = {}
    EFreads  = {}
    readhash = hashingread(infile)
    EFreads  = errorcorrect(readhash)
    cov = 0
    for bc in EFreads.keys():
      cov += EFreads[bc]
    outfile.write(infile+"\t"+str(cov)+"\n")
  outfile.close()

if __name__ == "__main__":
    main()
