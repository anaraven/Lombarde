#!/usr/bin/python

import csv, sys
from math import log
from getopt import getopt,GetoptError

try:
  opts, argv = getopt(sys.argv[1:], "lds1", ["log","dump","sum","single"])
except GetoptError as err:
  # print help information and exit:
  print(err) # will print something like "option -a not recognized"
  sys.exit(2)

dump = False
log_transform = False
sum_values = False
single = False
for o, a in opts:
  if o == "-d":
    dump = True
  elif o == "-l":
    log_transform = True
  elif o == "-s":
    sum_values = True
  elif o == "-1":
    single = True
  else:
    assert False, "unhandled option"

fimo_file  = argv[0] # "prodoric-NC_000913.fimo.txt"
blast_file = argv[1] # "TF-prodoric.blastp.txt"
tf_file    = argv[2] # "TF-motif-prodoric.txt"

def str2log(s):
  x = float(s)
  if x == 0.0:
    ans = -8000
  else:
    ans = int(10*log(x)/log(2))
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

print >>sys.stderr, "FIMO", len(FIMO)
if dump:
  for motif_id, lst in FIMO.iteritems():
    for gene_id, p in lst.iteritems():
      print "FIMO", motif_id, gene_id, p

tf_motif = dict() # keys are protein_id, values are motif_id
with file(tf_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    if len(row)<2: continue
    if row[0].startswith('#'): continue
    tf_motif[row[1]] = row[0]

# we read BLASTP output in tabular format
BLAST = dict() # keys are gene_id, values are (motif_id, score) if the gene codes a TF
with file(blast_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    # ignore comment lines
    if row[0].startswith('#'): continue
    # score can be E-value or bit-score
    p = str2log(row[10]) if log_transform else float(row[10]) # OR -float(row[11])
    gene_id = row[0]                    # gene id
    motif_id = tf_motif.get(row[1],None) # motif associated to the TF
    if motif_id is None: continue        # ignore genes without known motif 
    if BLAST.has_key(gene_id):             # if we saw this gene before
      oldmotif_id, oldp = BLAST[gene_id]   # update only if score is better
      if p<oldp:
        BLAST[gene_id] = (motif_id, p)
    else:
      BLAST[gene_id] = (motif_id, p)

      #oldmotif_id, oldp = BLAST.get(gene_id, (None,float(inf))   # update only if score is better
      #if p<oldp: BLAST[gene_id] = (motif_id, p)

print >>sys.stderr, "BLAST", len(BLAST)
if dump:
  for gene_id, mp in BLAST.iteritems():
    motif_id, p = mp
    print "BLAST", gene_id, motif_id, p

for gene_id,xp in BLAST.iteritems():
    motif_id, p = xp
    targets = FIMO.get(motif_id, {})
    #print "#", gene_id, motif_id, p, len(targets)
    for tgt_gene,V in targets.iteritems():
      if single:
        print "%s\t%s\t%G" % (gene_id, tgt_gene, 1-(1-p)*(1-V))
      elif sum_values:
        print "%s\t%s\t%G" % (gene_id, tgt_gene, p+V)
      else:
        print "%s\t%s\t%G\t%G" % (gene_id, tgt_gene, p, V)
