#!/usr/bin/python

import csv, sys

gene_file  = sys.argv[1] if len(sys.argv)>1 else "GeneProductSet.txt"
TUSet_file = sys.argv[2] if len(sys.argv)>2 else "TUSet.txt"

blattner = dict()
with file(gene_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    if row[0].startswith('#'): continue
    if row[2]=="":
      # if gene has no blattner id we use the gene name as id
      name = row[1]
    else:
      name = row[2]
    blattner[row[1]] = name
    
with file(TUSet_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    if row[0].startswith('#'): continue
    for i in row[3].split(","):
      print blattner[i]+"\t"+row[2]


