# combination of FIMO, BLAST and COUPLING to determine a putative gene-to-gene regulation network.
./build_fimo_blast_net.py Data/prodoric-NC_000913.fimo.clean.txt Data/TF-prodoric.blastp.clean.txt Data/TF-motif-prodoric.clean.txt > putativeTRN-genes.txt

# each line of `putativeTRN-genes.txt` has an arc per line: vertex-from, vertex-to, blast E-value and FIMO p-value

# now the E-values and p-values are discretized and converted into a single discrete value for each arc
./discretize-weight.R -i putativeTRN-genes.txt --n1 3 --n2 7 -o putativeTRN-genes.ncol
# the new file is on NCOL format. See <http://igraph.org/r/doc/read.graph.html>

# now the putative gene-to-gene regulation network is contracted to become a operon-to-operon one
./contract.R putativeTRN-genes.ncol Data/operon_list.clean.txt putativeTRN-operons.ncol

# the same contraction is applyed to the pairs of co-expressed genes
./contract.R Data/coexp-genes.ncol Data/operon_list.clean.txt coexp-operons.ncol

# finally the LOMBARDE method is applied
./lombarde.R -o output.ncol -a lombarde.out coexp-operons.ncol putativeTRN-operons.ncol

### optional step, requieres Potassco suite

# gringo optimize.lp | unclasp
