# this folder will contain results from scripts for comparison of aSNP frequency
# such files can be obtained using plink:

plink --vcf ${pop}_chr${chr}.vcf --maf 0.05 --r2 --ld-window-kb 200 --ld-window 99999 \
--ld-window-r2 0.8 --out /${pop}_chr${chr}_maf0.05_LD
