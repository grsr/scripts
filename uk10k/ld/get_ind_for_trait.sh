for CHR in {1..22} X
do 
    perl find_independent_snps.pl $1 $CHR > ~/ls113/uk10k/ioanna/pruned_02r2_01maf/${1}_${CHR}_ind.txt
done
