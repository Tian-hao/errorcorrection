#!/usr/bin/env python
import os
import sys
import glob
import string
import operator
from string import atof
from itertools import imap
from Bio import SeqIO

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

def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def hamming(str1, str2):
  assert len(str1) == len(str2)
  return sum(imap(operator.ne, str1, str2))

def callmut(refseq):
  pass

def readrecord(R1file,R2file,refseq,bchash,outfiles,outfiles_forward,outfiles_reverse):
  R2handle = SeqIO.parse(R2file,"fastq")
  countread = 0
  for record_R1 in SeqIO.parse(R1file,"fastq"):
    countread += 1
    if countread%1000000 == 0: print 'finish processing %d reads' % countread
    record_R2 = R2handle.next()
    ID_R1   = record_R1.id
    seq_R1  = str(record_R1.seq)
    ID_R2   = record_R2.id
    seq_R2  = str(record_R2.seq)
    assert(ID_R1==ID_R2)
    MID_R1     = seq_R1[4:7]
    MID_R2     = seq_R2[4:7]
    barcode_R1 = seq_R1[0:4]+seq_R1[7:11]
    barcode_R2 = seq_R2[0:4]+seq_R2[7:11]
    '''
    primer_R1  = seq_R1[11:31]
    primer_R2  = seq_R2[11:31]
    ROI_R1     = seq_R1[31:40]+seq_R1[76:79]
    ROI_R2     = seq_R2[31:34]+seq_R2[70:79]
    Qual_R1    = record_R1.letter_annotations["phred_quality"]
    Qual_R2    = record_R2.letter_annotations["phred_quality"]
    Qual_R1    = Qual_R1[31:40]+Qual_R1[76:79]
    Qual_R2    = Qual_R2[31:34]+Qual_R2[70:79]
    '''
    ROI_R1     = seq_R1[11:99]
    ROI_R2     = seq_R2[11:99]
    Qual_R1    = record_R1.letter_annotations["phred_quality"]
    Qual_R2    = record_R2.letter_annotations["phred_quality"]
    Qual_R1    = Qual_R1[11:99]
    Qual_R2    = Qual_R2[11:99]
    BC = barcode_R1+barcode_R2
    if min(Qual_R1) >= 0 and MID_R1==MID_R2 and MID_R1 in bchash.keys():
      outfiles_forward[MID_R1].write(BC+"\t"+ROI_R1+"\n")
    if min(Qual_R2) >= 0 and MID_R1==MID_R2 and MID_R2 in bchash.keys():
      outfiles_reverse[MID_R2].write(BC+"\t"+rc(ROI_R2)+"\n")
    if min(Qual_R1) >= 0 and min(Qual_R2) >= 0:
      if MID_R1==MID_R2 and MID_R1 in bchash.keys() and ROI_R1==rc(ROI_R2):
        outfiles[MID_R1].write(BC+"\t"+ROI_R1+"\n")
   # if min(Qual_R1) < 30 or min(Qual_R2) < 30: continue
   # if MID_R1==MID_R2 and MID_R1 in bchash.keys() and ROI_R1==rc(ROI_R2):
   # if min(Qual_R2) > 30 and MID_R2 in bchash.keys():
   #   outfiles[MID_R2].write(BC+"\t"+rc(ROI_R2)+"\n")

def hashinref(infile):
  refhash = {}
  handle = open(infile, "rU")
  for record in SeqIO.parse(handle, "fasta"):
    refhash[str(record.id)] = str(record.seq)
  handle.close()
  return refhash

def hashinbarcode(infile):
  bchash = {}
  infile = open(infile,'r')
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    bchash[line[1]] = line[0]
  return bchash

def main():
  refhash = hashinref("../Fasta/SeqInfo.fa") 
  bchash  = hashinbarcode('../Fasta/Barcode.fa')
  refseq  = refhash['ProteinG_AO'].replace(refhash['ForwardPrimer'],'').replace(rc(refhash['ReversePrimer']),'')
  filenames = sorted(glob.glob('../fastq/*R1*.fastq'))
  outfiles = {}
  outfiles_forward = {}
  outfiles_reverse = {}
  for bc in bchash.keys():
    outfiles[bc] = open('../paired/'+bchash[bc],'w')
    outfiles_forward[bc] = open('../paired/'+bchash[bc]+'Forward','w')
    outfiles_reverse[bc] = open('../paired/'+bchash[bc]+'Reverse','w')
  countfile = 0
  for R1file in filenames:
    countfile += 1
    print 'working on %d PE file' % countfile
    R2file = R1file.replace('_R1_','_R2_')
    readrecord(R1file,R2file,refseq,bchash,outfiles,outfiles_forward,outfiles_reverse)
  for bc in bchash.keys(): outfiles[bc].close()  

if __name__ == "__main__":
  main()

