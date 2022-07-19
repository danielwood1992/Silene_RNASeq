fasta=$1;
grep '>' $fasta | cut -f1,2 -d' ' | sed 's/^.//' | sed 's/len=//' | sed 's/ /\t/' | sed 's/_i.*\t/\t/' > $fasta.lengths;
