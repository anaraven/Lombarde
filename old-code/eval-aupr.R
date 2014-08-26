library(minet)
names <- c("prodoric-corr-pearson-3.Rdata",
"prodoric-corr-pearson-5.Rdata",
"prodoric-mi-pearson-aracne-3.Rdata",
"prodoric-mi-pearson-aracne-5.Rdata",
"prodoric-mi-pearson-c3-3.Rdata",
"prodoric-mi-pearson-c3-5.Rdata",
"prodoric-mi-pearson-clr-3.Rdata",
"prodoric-mi-pearson-clr-5.Rdata",
"prodoric-mi-pearson-mrnet-3.Rdata",
"prodoric-mi-pearson-mrnet-5.Rdata",
"rdb-corr-pearson-3.Rdata",
"rdb-mi-pearson-aracne-3.Rdata",
"rdb-mi-pearson-c3-3.Rdata",
"rdb-mi-pearson-clr-3.Rdata",
"rdb-mi-pearson-mrnet-3.Rdata",
"rdb-corr-pearson-5.Rdata",
"rdb-mi-pearson-aracne-5.Rdata",
"rdb-mi-pearson-c3-5.Rdata",
"rdb-mi-pearson-clr-5.Rdata",
"rdb-mi-pearson-mrnet-5.Rdata")

for (fname in names) {
	load(fname)
	cat(fname, auc.pr(table.found(found)), auc.pr(table.found(found.restr)),"\n",sep="\t")
}
