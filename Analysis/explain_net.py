#!/usr/bin/python

# here we build the initial putative network and compare it to the gold-standard
# reference, to determine why there are missing arcs in the prediction

import csv, sys
from math import log
from getopt import getopt,GetoptError

try:
  opts, argv = getopt(sys.argv[1:], "lds", ["log","dump","sum"])
except GetoptError as err:
  # print help information and exit:
  print(err) # will print something like "option -a not recognized"
  sys.exit(2)

dump = False
log_transform = False
sum_values = False
for o, a in opts:
  if o == "-d":
    dump = True
  elif o == "-l":
    log_transform = True
  elif o == "-s":
    sum_values = True
  else:
    assert False, "unhandled option"

fimo_file  = argv[0] # "prodoric-NC_000913.fimo.txt"
blast_file = argv[1] # "TF-prodoric.blastp.txt"
tf_file    = argv[2] # "TF-motif-prodoric.txt"
net_file   = argv[3] # 

def str2log(s):
  x = float(s)
  if x == 0.0:
    ans = -500
  else:
    ans = int(log(x)/log(2))
  return ans

# in the FIMO file lower scores are better. We keep the minimum score for each arc
FIMO = dict() # keys are motif_id, values are dictionaries with gene_id as key and score as value
with file(fimo_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    if row[0].startswith('#'): continue
    # score can be p-value or q-value
    motif_id = row[0]
    gene_id  = row[1]
    p = str2log(row[2]) if log_transform else float(row[2])
    if FIMO.has_key(motif_id): # if we saw the motif_id before...
      old = FIMO[motif_id].get(gene_id, p+1) # in the upstream of the same gene...
      if p<old:                            # or upstream a new gene...
        FIMO[motif_id][gene_id] = p          # we keep the best score for the pair
    else:
      FIMO[motif_id] = {gene_id:p} # in the motif is new we create a dictionary

tf_motif = dict() # keys are protein_id, values are motif_id
with file(tf_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    if len(row)<2: continue
    if row[0].startswith('#'): continue
    tf_motif[row[1]] = row[0]

inf = float("inf")
# we read BLASTP output in tabular format
BLAST = dict() # keys are gene_id, values are (motif_id, score) if the gene codes a TF
with file(blast_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    # ignore comment lines
    if row[0].startswith('#'): continue
    # score can be E-value or bit-score
    p = str2log(row[10]) if log_transform else float(row[10]) # OR -float(row[11])
    gene_id = row[0] # gene id
    prot_id = row[1] # TF protein ID
    oldmotif, oldp = BLAST.get(gene_id, (None,inf))   # update only if score is better
    if p<oldp: BLAST[gene_id] = (prot_id, p)

with file(net_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    # ignore comment lines
    if row[0].startswith('#'): continue
    gene_src = row[0]
    gene_tgt = row[1]
    if gene_src=="" or gene_tgt=="":
      print "NO", gene_src, gene_tgt
      continue
    if BLAST.has_key(gene_src):
      (prot_id, p) = BLAST[gene_src]
      if tf_motif.has_key(prot_id):
        motif_id = tf_motif[prot_id]
        lst = FIMO.get(motif_id,{})
        if lst.has_key(gene_tgt):
          print "OK", gene_src, gene_tgt, prot_id, p, motif_id, lst[gene_tgt]
        else:
          print "BS", gene_src, gene_tgt, prot_id, p, motif_id
      else:
        print "MT", gene_src, gene_tgt, prot_id, p
    else:
      print "TF", gene_src, gene_tgt

