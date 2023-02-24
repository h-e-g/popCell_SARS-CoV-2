source(sprintf('%s/single_cell/resources/template_scripts/querySNPs.R',EIP))
SNP_info=getMap(annotate=TRUE)
setnames(SNP_info,"ID","rsID")
setnames(SNP_info,"chromposID_hg38","ID")

freqs=fread(sprintf("%s/users/Javier/results/freq/CHS_CEU_YRI_allCHR_popCellSNVs_nomono_hg38.frq.strat.gz",EIP))

freqs_wide=merge(freqs,SNP_info[,.(ID,REF,ALT,ANCESTRAL)],by.x='SNP',by.y='ID')
freqs_wide[,DAF_or_MAF:=ifelse(A2==ANCESTRAL,MAF,ifelse(A1==ANCESTRAL,1-MAF,MAF))]
freqs_wide[,ALT_DERANC:=ifelse(A2==ANCESTRAL,'DER',ifelse(A1==ANCESTRAL,"ANC","DER"))]
freqs_wide=dcast(unique(freqs_wide),CHR+SNP+A1+A2+ANCESTRAL+ALT_DERANC~CLST,value.var=c('DAF_or_MAF'))
freqs_wide[,DAF_or_MAF_CEU:=CEU]
freqs_wide[,DAF_or_MAF_CHS:=CHS]
freqs_wide[,DAF_or_MAF_YRI:=YRI]
SNP_info=merge(SNP_info,freqs_wide[,.(SNP,A1,A2,ALT_DERANC,DAF_or_MAF_CEU,DAF_or_MAF_CHS,DAF_or_MAF_YRI)],by.x=c('ID','REF'),by.y=c('SNP','A2'),all.x=TRUE)

fst=fread(sprintf("%s/users/Javier/results/pbs/FST_CEU_CHS_YRI_allCHR_popCellSNVs_nomono_hg38.pbs.txt.gz",EIP))

natsel_res_all=lapply(c("CHS","CEU","YRI"),function(pop){
  pbs_res=fread(sprintf("%s/users/Javier/results/annotation/%s_ALLCHR_popCellSNVs_nomono_hg38_PBSwindows.txt.gz",EIP,pop))
  # pbs_res=pbs_res[ID%in%unique(eqtl_file$ID),]
  pbs_res[,POP:=pop]
  print(sprintf("%s.",pop))
  return(pbs_res)
}) %>% rbindlist()

natsel_res_dcast=dcast(natsel_res_all,CHROM+POS+ID~POP,value.var=c('PBS','P_PBS'))

SNP_info=merge(SNP_info,freqs_wide[,.(SNP,A1,A2,ALT_DERANC,DAF_or_MAF_CEU,DAF_or_MAF_CHS,DAF_or_MAF_YRI)],by.x=c('ID','REF'),by.y=c('SNP','A2'),all.x=TRUE)

SNP_info=merge(SNP_info,natsel_res_dcast,by='ID',all.x=TRUE)

SNP_info_basics=SNP_info[,.(rsID, posID=ID, CHROM, POS_B38=POS, POS_B37=core_pos, REF, ALT, ANCESTRAL, DR2, CLASS, ALT_freq_AFB=AFB, ALT_freq_ASH=ASH, ALT_freq_EUB=EUB, max_MAF, DistGene, NearestGene)]
fwrite(SNP_info_basics,file=sprintf('%s/users/Maxime/popCell_SARS-CoV-2/7__genotyping_data/data/SNP_info_basics.tsv.gz',EIP),sep='\t')

SNP_info_covid=SNP_info[,.(rsID, posID=ID, covid_A2_pval, covid_A2_beta, covid_A2_sebeta, covid_B2_pval, covid_B2_beta, covid_B2_sebeta, covid_C2_pval, covid_C2_beta, covid_C2_sebeta)]
fwrite(SNP_info_covid,file=sprintf('%s/users/Maxime/popCell_SARS-CoV-2/7__genotyping_data/data/SNP_info_covid.tsv.gz',EIP),sep='\t')

# add LD scores
ld_score=list()
for (pop in c('CEU','CHS','YRI')){
  ld_score[[pop]]=fread(sprintf("%s/users/Javier/results/ld_scores/%s_ALLCHR_popCellSNVs_nomono_maf0.05_hg38_ldscores.txt.gz",EIP,pop))
}
ld_score=rbindlist(ld_score,idcol='pop')
ld_score_pop=dcast(ld_score,CHROM+POS+ID~pop,value.var="LDscore")

SNP_info=merge(SNP_info,ld_score_pop[,.(ID,LD_score_YRI=YRI,LD_score_CEU=CEU,LD_score_CHS=CHS)],by='ID',all.x=TRUE)
SNP_info=merge(SNP_info,fst,by.x='posID',by.y='SNP',all.x=TRUE)

SNP_info_popgen=SNP_info[,.(rsID, posID=ID, present_1kG_HiCov_b38=!is.na(A1),DAF_or_MAF_YRI, DAF_or_MAF_CEU, DAF_or_MAF_CHS, PBS_YRI, PBS_CEU, PBS_CHS, P_PBS_YRI, P_PBS_CEU, P_PBS_CHS,LD_score_YRI,LD_score_CEU,LD_score_CHS)]

fwrite(SNP_info_popgen,file=sprintf('%s/users/Maxime/popCell_SARS-CoV-2/7__genotyping_data/data/SNP_info_popgen.tsv.gz',EIP),sep='\t')

#### add annotation aSNPs


introgressed_snps <- fread(sprintf("%s/users/Javier/results/asnps/Final_lenient_Adaptive_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt",EIP))
archaic_origin_snps = fread(sprintf("%s/users/Javier/results/asnps/Final_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt",EIP))
archaic_origin_snps=dcast(archaic_origin_snps,CHROM+POS+ID+ANC+DER+ASNP+ORIGIN~Population,value.var=c('ASNP_FREQ','METHOD'))
SNP_info_introgressed=dcast(introgressed_snps[HAP_ORIGIN!='ANY',],CHROM+POS+REF+ALT+ID~Population,value.var=c('HAP_ORIGIN','Q99','NUM_aSNPs'))

SNP_info=merge(SNP_info,SNP_info_introgressed[,.(posID=ID,HAP_ORIGIN_CEU, HAP_ORIGIN_CHS, Q99_CEU, Q99_CHS, NUM_aSNPs_CEU, NUM_aSNPs_CHS)],by='posID',all.x=TRUE)

SNP_info[is.na(HAP_ORIGIN_CEU),HAP_ORIGIN_CEU:='']
SNP_info[is.na(HAP_ORIGIN_CHS),HAP_ORIGIN_CHS:='']
SNP_info[,Q99_CEU:=ifelse(is.na(Q99_CEU) | Q99_CEU==FALSE,'','yes')]
SNP_info[,Q99_CHS:=ifelse(is.na(Q99_CHS) | Q99_CHS==FALSE,'','yes')]
SNP_info[is.na(NUM_aSNPs_CEU),NUM_aSNPs_CEU:=0]
SNP_info[is.na(NUM_aSNPs_CHS),NUM_aSNPs_CHS:=0]
SNP_info=merge(SNP_info,archaic_origin_snps[,.(posID=ID,ANC, DER, ASNP,  ORIGIN, ASNP_FREQ_CEU, ASNP_FREQ_CHS, METHOD_CEU, METHOD_CHS)],by='posID',all.x=TRUE)
SNP_info[is.na(METHOD_CEU),METHOD_CEU:='']
SNP_info[is.na(METHOD_CHS),METHOD_CHS:='']

SNP_info_archaics=SNP_info[,.(rsID,posID,ANC, DER, ASNP, ORIGIN, ASNP_FREQ_CEU, ASNP_FREQ_CHS, METHOD_CEU, METHOD_CHS, HAP_ORIGIN_CEU, HAP_ORIGIN_CHS, Q99_CEU, Q99_CHS, NUM_aSNPs_CEU, NUM_aSNPs_CHS)]
LD_pruned_CEU=fread(sprintf("%s/users/Javier/VCF/prune/tagSNPs_CEU_maf0.05pct_LD0.8_allCHR.prune.in",EIP),header=F)$V1
LD_pruned_CHS=fread(sprintf("%s/users/Javier/VCF/prune/tagSNPs_CHS_maf0.05pct_LD0.8_allCHR.prune.in",EIP),header=F)$V1
SNP_info_archaics[,tagSNP_CEU:=posID%in%LD_pruned_CEU]
SNP_info_archaics[,tagSNP_CHS:=posID%in%LD_pruned_CHS]
fwrite(SNP_info_archaics,file=sprintf('%s/users/Maxime/popCell_SARS-CoV-2/7__genotyping_data/data/SNP_info_archaics.tsv.gz',EIP),sep='\t')

CHS_LD=fread("/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/CHS_LDblocks_GOOD_aSNPs_popCellSNVs_nomono_hg38_biSNPs_REF2AA_phased.txt")
CEU_LD=fread("/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/CEU_LDblocks_GOOD_aSNPs_popCellSNVs_nomono_hg38_biSNPs_REF2AA_phased.txt")

# SNP_info_ARCHAICS=SNP_info[,.(aSNP_yes_no, introgressed_yes_no, Archaic_freq_CEU, Archaic_freq_CHS, tagSNP_YRI, tagSNP_CEU, tagSNP_CHS)]
