for TRAIT in BMI Weight Height Hip Waist WHR HipBMIadj WaistBMIadj WHRBMIadj TRFM TFM TLM 
do
    export TRAIT
#    intersectBed -a ../pruned_02r2_01maf/${TRAIT}.bed -b e74_coding_genes_50kb_flanks_normal_chrs.bed -wb > traits/${TRAIT}/snps_genes.bed
#    intersectBed -a ../pruned_02r2_01maf/${TRAIT}_novel.bed -b e74_coding_genes_50kb_flanks_normal_chrs.bed -wb > traits/${TRAIT}/novel_snps_genes.bed
#    python correct_pvals.py traits/${TRAIT}/snps_genes.bed traits/${TRAIT}/genes_pvals.csv
#    python correct_pvals.py traits/${TRAIT}/novel_snps_genes.bed traits/${TRAIT}/novel_genes_pvals.csv
#    python run_geneset.py traits/${TRAIT}/genes_pvals.csv msigdb/kegg_20_200_genes.txt traits/${TRAIT}/kegg
#    python run_geneset.py traits/${TRAIT}/novel_genes_pvals.csv msigdb/kegg_20_200_genes.txt traits/${TRAIT}/novel_kegg
#    mkdir traits/${TRAIT}/kegg_results
#    mkdir traits/${TRAIT}/novel_kegg_results
#    gsub -q normal -a 158 -m 7 "perl gsea_wrapper.pl traits/${TRAIT}/kegg/ traits/${TRAIT}/kegg_results/ 10000"
#    gsub -q normal -a 158 -m 7 "perl gsea_wrapper.pl traits/${TRAIT}/novel_kegg/ traits/${TRAIT}/novel_kegg_results/ 10000"
#    python compute_pw_stats.py traits/${TRAIT}/kegg_results/ traits/${TRAIT}/kegg_summary.csv
#    python compute_pw_stats.py traits/${TRAIT}/novel_kegg_results/ traits/${TRAIT}/novel_kegg_summary.csv
    perl -F, -lane 'next if /fwer/; print $ENV{TRAIT}.",".$_ if $F[2] < 0.01' < traits/${TRAIT}/novel_kegg_summary.csv
done

