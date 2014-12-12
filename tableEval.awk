#!/usr/bin/awk -f

# Usage:
# ./tableEval.awk Input/RDB8.1/gold-std.ncol Out1/EvalBase/*/RegulonDB/MEME*/9/*/gl.out | tr / '\t' | cut -f 3-7,9-

BEGIN {
	if(!MIN) MIN = 1
	if(!MAX) MAX = 1000000
	if(N)    MIN = MAX = N
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	    "out1/eval/Coexpr/PWM/Network/Levels/Base/gl.out",
	    "n_coexp", "n_arcs", "n_valid", "n_vertex", "avg_degree", "n_explanation", "n_vshape"
}

function print_line() {
    if(n_coexp>=MIN && n_coexp<=MAX) {
	printf "%s\t%d\t%d\t%d\t%d\t%5.1f\t%d\t%d\n", FILENAME, n_coexp, n_arcs, n_valid,
		n_vertex, n_arcs/n_vertex, n_explanation, n_vshape
    }
}

FILENAME==ARGV[1] {
	gs[$1,$2] = 1
       	next
}

FILENAME!=last_file {
	print_line()
	n_arcs   = 0
	n_valid  = 0
	n_coexp  = 0
	n_vertex = 0
	n_explanation = 0
	n_vshape = 0
	delete arc
	delete vertex
}

/arcInVshape/ && !arc[$3,$4] {
	n_arcs++
	arc[$3,$4] = 1
	if(gs[$3,$4]) n_valid++
}

/arcInVshape/ && !vertex[$3] {
	n_vertex++
	vertex[$3] = 1
}

/arcInVshape/ && !vertex[$4] {
	n_vertex++
	vertex[$4] = 1
}

/vshape/ {
	if(n_coexp != $5)
	    print_line()
	n_coexp = $5
	n_vshape++
	last_file=FILENAME
}

/explanation/ {
	n_explanation++
}

END {
	print_line()
}

# TODO: FIX: last `explanation` is not counted properly
