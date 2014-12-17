# Lombarde
*Andrés Aravena*, Dec 17, 2014

Combines a putative transcriptional regulatory network (such as those predicted using MEME) and co-expression data (such as microarray experiments results) to produce a realistic transcriptional regulatory network.

*Lombarde* models transcriptional regulatory networks in a scheme that integrates putative transcriptional regulatory networks with co-expression data to determine the simplest and most confident sub-network that explains the observed co-expressions.
              
The output will be a subgraph of the putative transcriptional regulatory network that satisfies the *Lombarde* criteria: each pair of co-expressed vertices should share a common regulator (either direct or via a regulation cascade), among all the common regulators select the most confident ones.

## Two tools
The *Lombarde* model is directly implemented on the `lombarde.R` script, whose inputs are graphs and whose outputs are a new graph and optionally a log on how this new graph is produced. This is the fists tool provided.

These inputs graphs are derived from experimental data that is preprocessed by standard tools like `MEME/FIMO`, `BLAST` and `MRNET`, and by ad-hoc scripts such as `build_fimo_blast_net.py`, `discretize-weight.R` and `contract.R`.

For a first approach to this suite of tools we also provide a “all-in-one” tool, named `lombarde-full.sh`. This is essentially a wrapper to all the ad-hoc scripts so all preprocessing is carried on automatically.

## Input
Currently most of the files represent graphs in NCOL format as parsed by *igraph* library on *R* (<http://igraph.org/r/doc/read.graph.html>).

The basic `lombarde.R` tool is then invoked as:

		lombarde.R -o output.ncol -a output.log coexp.ncol putativeTRN.ncol

where `coexp.ncol` represent the set of co-expressed elements, one pair per line, and `putativeTRN.ncol` represents the putative transcriptional network built based on the output of BLAST and MEME/FIMO.

This last input file can be built doing:

		build_fimo_blast_net.py fimo.txt blastp.txt coupling.txt > putativeTRN.txt