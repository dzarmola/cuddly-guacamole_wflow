inp=$(echo $1 | sed 's/\..*//g')
inf=$(echo $2 | sed 's/\.//g')
mcxload -abc $inp.out --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $inp.mci -write-tab $inp.tab
mcl $inp.mci -I $2
mcxdump -icl out.$inp.mci.I$inf -tabr $inp.tab -o dump.$inp.mci.I$inf
cp dump.$inp.mci.I$inf $3 #clustering_results_"$inp"_"$2"
