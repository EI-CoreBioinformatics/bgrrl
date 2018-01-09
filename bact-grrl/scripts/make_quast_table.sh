#!/bin/bash

echo $'Assembly\n# contigs (>= 0 bp)\n# contigs (>= 1000 bp)\n# contigs (>= 5000 bp)\n# contigs (>= 10000 bp)\n' \
     $'# contigs (>= 25000 bp)\n# contigs (>= 50000 bp)\nTotal length (>= 0 bp)\nTotal length (>= 1000 bp)\n' \
     $'Total length (>= 5000 bp)\nTotal length (>= 10000 bp)\nTotal length (>= 25000 bp)\nTotal length (>= 50000 bp)\n' \
     $'# contigs\nLargest contig\nTotal length\nGC (%)\nN50\nN75\nL50\nL75\n# Ns per 100 kbp' > quast.tsv

for f in $(find . -name 'report.tsv' | sort); do 
 cp quast.tsv quast.tsv.tmp
 paste quast.tsv.tmp <(cut -f 2 $f) > quast.tsv
done

awk 'BEGIN { FS=OFS="\t" }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}' quast.tsv > quast.tsv.tmp

mv quast.tsv.tmp quast.tsv
