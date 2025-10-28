mkdir -p ref; perl enrich2gmt.pl /public/Database/pipRef/Gallus_gallus/Ensembl_release110/ref4SingleCell/Ensembl_release110/annot/ref ref 
perl ready_gsea.pl all.xls ./ref ./
perl run_GSEA.pl  -run_go 0 -run_ko 1 -run_reactome 1 -run_do 0 -run_sge yes -compare group_diff.list -group group.list -enrich_gmt ref -outdir ./ -exp exp.xls 
