#!/usr/bin/python

import csv, igraph, sys

op_file = sys.argv[2] # operon_list.txt
op = dict()
with file(op_file, "r") as fp:
  for row in csv.reader(fp, delimiter="\t"):
    if row[0].startswith('#'): continue
    op[row[0]] = row[1]

net_file = sys.argv[1] # gold_std_genes.txt
net = igraph.Graph.Read_Ncol(net_file, weights=True, names=True, directed=True)
print >>sys.stderr,net.vcount(),"vertices,",net.ecount(),"edges read"
net.delete_vertices([v.index for v in net.vs if not op.has_key(v["name"])])
print >>sys.stderr,net.vcount(),"vertices,",net.ecount(),"edges with assigned class"
equ = [op[v["name"]] for v in net.vs]
#unq = {x:1 for x in equ}
#nid = {y:x for x,y in enumerate(unq.keys())}
#oid = {x:y for x,y in enumerate(unq.keys())}
unq = dict([(x,1) for x in equ])
nid = dict([(y,x) for x,y in enumerate(unq.keys())])
oid = dict([(x,y) for x,y in enumerate(unq.keys())])
equiv = [nid[v] for v in equ]

net.contract_vertices(equiv)
print >>sys.stderr,net.vcount(),"vertices,",net.ecount(),"edges after contraction"
net.vs["name"] = [oid[v.index] for v in net.vs]
net.simplify(combine_edges="min")
print >>sys.stderr,net.vcount(),"vertices,",net.ecount(),"edges after reduction"
net.write_ncol(sys.stdout)

