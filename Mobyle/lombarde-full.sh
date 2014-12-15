#!/bin/bash

# wrapper to encapusulate all stages in one script.
# all parameters are given in the command line.

# This wrapper is intended to be called from Mobyle interface. A more flexible
# approach is to use each sub-program separately.

# where to find the R and python programs
H=/genoma/homes/mobyle
CODE=$H/softs/lombarde/current
export LD_LIBRARY_PATH=$H/softs/R/3.1.1/lib64/R/lib/

# the 1st parameter is the output filename
OUT=$1

# second parameter is Binding Site prediction file. FIMO output format
FIMO=$2

# 3rd parameter is Transcriptin factor prediction file. BLAST -m 8 format
BLAST=$3

# 4th parameter is the TF-motif relationship provided by the TF database
COUPLING=$4

# 5th parameter. Pairs of co-expressed genes. One pair per line
COEXP=$5

# 6th param. Number of discretization levels for the Binding Sites p-values
LVL2=$6

# 7th param. Operon membership.
OPERONS=$7

# hidden parameter. Number of discretization levels for Transcription Factor E-values.
LVL1=3
# shift

# echo out=$OUT
# echo FIMP=$FIMO
# echo BLAST=$BLAST
# echo COUPLING=$COUPLING
# echo LEVELS=$LVL2
# echo COEXP=$COEXP
# echo OPERONS=$OPERONS

# preprocessing: cleaning of FIMO output to keep only BS in the same strand as the regulated gene
awk -F"\t" -v OFS="\t" '/#/ {next} $3>$4 {next} {print $1, $2, $6}' $FIMO > fimo-clean.txt

# combination of FIMO, BLAST and COUPLING to determine a putative gene-to-gene regulation network.
python $CODE/build_fimo_blast_net.py fimo-clean.txt $BLAST $COUPLING > putativeTRN-genes.txt 2> /dev/null

# each line of `putativeTRN-genes.txt` has an arc per line: vertex-from, vertex-to, blast E-value and FIMO p-value

# now the E-values and p-values are discretized and converted into a single discrete value for each arc
Rscript --vanilla $CODE/discretize-weight.R -i putativeTRN-genes.txt --n1 $LVL1 --n2 $LVL2 -o putativeTRN-genes.ncol
# the new file is on NCOL format. See <http://igraph.org/r/doc/read.graph.html>

# now the putative gene-to-gene regulation network is contracted to become a operon-to-operon one
Rscript --vanilla $CODE/contract.R putativeTRN-genes.ncol $OPERONS putativeTRN-operons.ncol > /dev/null 2> /dev/null

# the same contraction is applyed to the pairs of co-expressed genes
Rscript --vanilla $CODE/contract.R $COEXP $OPERONS coexp-operons.ncol > /dev/null 2> /dev/null

# finally the LOMBARDE method is applied
Rscript --vanilla $CODE/lombarde.R -o $OUT coexp-operons.ncol putativeTRN-operons.ncol 2> /dev/null

# postprocessing: clean all intermendiate files
rm -f fimo-clean.txt putativeTRN-genes.ncol putativeTRN-genes.txt putativeTRN-operons.ncol coexp-operons.ncol
