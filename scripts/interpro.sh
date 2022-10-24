for i in *
do
/data/zeynep/interproscan-5.47-82.0/interproscan.sh -cpu 30 -d /data/zeynep/genomics/interpro/new_sp -f tsv,json,gff3 -goterms -i $i -iprlookup --pathways > $i\_.log &
done
