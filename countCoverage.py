#!/usr/bin/env python
import os
import sys
import string

total = 0
infile = open(sys.argv[1],'r')
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  count = int(line[1])
  total += count
print sys.argv[1]
print total
