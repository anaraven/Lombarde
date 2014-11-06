#!/usr/bin/awk -f

FILENAME==ARGV[1] {gs[$1,$2]=1; next}

/explanation/ {n_coexp++}

/arcInVshape/ {
	arc[$3,$4]=1
}

/explanation/ {
	n_arcs=0
	n_valid=0
	for(x in arc) {
		n_arcs++
		n_valid += gs[x]
	}
	print n_coexp, n_arcs, n_valid
}

