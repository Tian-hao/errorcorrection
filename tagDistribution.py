#!/usr/bin/env python

import os
import sys
import glob
from collections import Counter

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
  outfile = open(outfile,'w')
  for tag in H.keys():
    outfile.write(tag+"\t"+str(H[tag])+"\n")
  outfile.close()

def noerrorcorrect(H):
  EFreads = {}
  for bc in H.keys():
    EFreads[bc] = len(H[bc])
  return EFreads

def errorcorrect(H):
  EFreads = {}
  for bc in H.keys():
    if len(H[bc]) >= 3 and len(set(H[bc])) == 1:
      EFreads[bc] = len(H[bc])
  return EFreads

filenames = sorted(glob.glob('../paired/*'))
for infile in filenames:
  print 'working on %s' % infile
  readhash = {}
  count = {}
  readhash = hashingread(infile)
  noEFreads = noerrorcorrect(readhash)
  EFreads = errorcorrect(readhash)
  DNAfile  = infile.replace('paired/','tagDistribution/')+'.cov'
  outputing(noEFreads,DNAfile)
  DNAfile = infile.replace('paired/','tagDistribution/Tag')+'.cov'
  outputing(EFreads,DNAfile)
