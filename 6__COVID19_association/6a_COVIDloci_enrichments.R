
# SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
# EQTL_SCRIPT_DIR="09_eQTLmapping"
# EQTL_SCRIPT_DIR=""
# sbatch --array=1-132 --parsable --mem=30G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_gwas_enrich_%A_%a.log -J gwas_enrich ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/11c_COVIDloci_enrichments.R --nsamp 10000 --set
# sbatch --array=1-256 --parsable --mem=30G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_gwas_enrich_%A_%a.log -J gwas_enrich ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/11c_COVIDloci_enrichments.R --nsamp 10000 --set
# sbatch --qos=geh -p geh --array=2,8,14,88,89,93,94,99,107,108,109,226,227 --parsable --mem=80G -o ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_%A_%a.log -J gwas_enrich ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/11c_COVIDloci_enrichments.R --nsamp 10 --set

# sbatch --qos=geh -p geh --array=107-109,190-197,225-227,255-256 --parsable --mem=120G -o ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_%A_%a.log -J gwas_enrich ${SCRIPT_DIR}/00_Rscript.sh 11c_COVIDloci_enrichments.R --nsamp 10000 --set
# job 26129377
# TODO: look for Errors (TIME LINIT | memory) in logs
# tail ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129377_*.log
# ls ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129377_*.log | grep error
# ls ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129377_*.log | grep LIMIT
# ls ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129377_*.log | grep memory

# sbatch --qos=geh -p geh --array=1-22,87-91,110-131,175-179,198-224,228-254%26 --parsable --mem=120G -o ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_%A_%a.log -J gwas_enrich ${SCRIPT_DIR}/00_Rscript.sh 11c_COVIDloci_enrichments.R --nsamp 10000 --set
#job 26129095

#TODO: look for Errors (TIME LINIT | memory)
# ls ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129095_*.log | wc
# ls ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129095_*.log | grep error
# ls ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129095_*.log | grep LIMIT
# ls ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_gwas_enrich_26129095_*.log | grep memory

options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(mashr))
suppressMessages(library(tictoc))
suppressMessages(library(readr))
suppressMessages(library(rtracklayer))

EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
FIGURE_DIR = sprintf("%s/single_cell/project/pop_eQTL/figures/",EVO_IMMUNO_POP_ZEUS)
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)
eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
OUT_DIR = eQTL_DIR
SCRATCH="/pasteur/appa/scratch/mrotival/"

theme_set(theme_bw())
theme_update(
  text=element_text(family="serif",size=12),
  panel.grid=element_blank(),legend.position="bottom",
  strip.background=element_rect(fill="#012158"),strip.text=element_text(color="white")
)
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EVO_IMMUNO_POP_ZEUS))
source(sprintf("%s/single_cell/resources/template_scripts/GOSeq.R",EVO_IMMUNO_POP_ZEUS))
source(sprintf("%s/single_cell/resources/template_scripts/querySNPs.R",EVO_IMMUNO_POP_ZEUS))


NSAMPLE=10000
SET_NUM=1

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--set' | cmd[i]=='-s' ){SET_NUM = cmd[i+1]} # ID of the set to test
  if (cmd[i]=='--nsamp' | cmd[i]=='-n' ){NSAMPLE = cmd[i+1]} # ID of the set to test
}

plot_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/Enrichments_matchedDist',EVO_IMMUNO_POP_ZEUS)
dir.create(plot_DIR)

CIS_DIST=1e5

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
#dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",EVO_IMMUNO_POP_ZEUS))

###############################################
########### SNP informations genome wide ######
###############################################
tic('loding SNP info')
SNP_info=getMap(annotate=TRUE)
toc()

# COVID_A2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
# COVID_A2[,SNP:=gsub(':','_',SNP)]
# mm=match(SNP_info$posID,COVID_A2$SNP)
# SNP_info[,covid_A2_pval:=COVID_A2[mm,all_inv_var_meta_p]]
# SNP_info[,covid_A2_beta:=COVID_A2[mm,all_inv_var_meta_beta]]
# SNP_info[,covid_A2_sebeta:=COVID_A2[mm,all_inv_var_meta_sebeta]]
#
# COVID_B2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
# COVID_B2[,SNP:=gsub(':','_',SNP)]
# mm=match(SNP_info$posID,COVID_B2$SNP)
# SNP_info[,covid_B2_pval:=COVID_B2[mm,all_inv_var_meta_p]]
# SNP_info[,covid_B2_beta:=COVID_B2[mm,all_inv_var_meta_beta]]
# SNP_info[,covid_B2_sebeta:=COVID_B2[mm,all_inv_var_meta_sebeta]]
#
# COVID_C2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
# COVID_C2[,SNP:=gsub(':','_',SNP)]
# mm=match(SNP_info$posID,COVID_C2$SNP)
# SNP_info[,covid_C2_pval:=COVID_C2[mm,all_inv_var_meta_p]]
# SNP_info[,covid_C2_beta:=COVID_C2[mm,all_inv_var_meta_beta]]
# SNP_info[,covid_C2_sebeta:=COVID_C2[mm,all_inv_var_meta_sebeta]]
#
# ### Define distance to nearest gene
# #SNP_info[,DistGene:=XXXX]
# ## get gene list
# gtf <- rtracklayer::import(sprintf("%s/single_cell/resources/references/RNA/human_iav_sars/genes/genes.gtf",EVO_IMMUNO_POP_ZEUS))
# # gtf <- rtracklayer::import(sprintf("%s/single_cell/resources/references/RNA/human_iav_sars/genes/genes.gtf",EVO_IMMUNO_POP_ZEUS))
# Feature_annot=as.data.table(gtf)[type=='gene',.(gene_id,gene_name,seqnames, start, end, strand,gene_type)]
# Feature_annot[is.na(gene_name),gene_name:=gene_id]
# Feature_annot[,seqnames:=as.character(seqnames)]
# # for IAV_M and IAV_NS  we have two lines with the same genename (2 transcripts)
# # we only keep the longest transcript
# Feature_annot=Feature_annot[!(gene_id=="IAV_M" & end==784) & !(gene_id=="IAV_NS" & end==719),]
#
# ## get autosomes only
# gl <- as_tibble(as.data.frame(Feature_annot))
# unique(gl$gene_type)
# gl <- na.omit(gl)
# gl <- gl[gl$seqnames %in% paste0("chr",seq(1:22)),]
# library(purrr)
# gl$CHROM <- as.numeric(map(strsplit(gl$seqnames, "chr"), ~.x[2]) %>% unlist())
#
# ## get all sites
# m_maf_gr <- GRanges(seqnames = SNP_info$CHROM,
#         ranges = IRanges(start = SNP_info$POS, end = SNP_info$POS))
#
# gl_gr <- GRanges(seqnames = gl$seqnames,
#                     ranges = IRanges(start = gl$start, end = gl$end))
# dist_gr <- distanceToNearest(m_maf_gr, gl_gr)
# dist_df <- as.data.table(dist_gr)
# rm(dist_gr)
#
# all(dist_df$queryHits==seq(1:nrow(SNP_info)))
# SNP_info$DistGene <- dist_df$distance
# SNP_info$NearestGene <- gl$gene_name[dist_df$subjectHits]
# fwrite(SNP_info,file=sprintf("%s/popCell_data/02_ImputedGenotypes/Imputed/b38/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct_annotated.covid.balancing.iHS.FST.Age.DistGene.ANCESTRAL.CHS.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')

###############################################
###############################################
########### resampling analyses ###############
###############################################
###############################################

########### resampling functions
#### step 1: resample SNPs and extract relevant informations
# resample_matched_on_MAF=function(snpSet,NSAMP=1000){
#   eQTL_AF=SNP_info[match(snpSet,ID),AF]
#   SNP_info[,is_eQTL:=ID%chin%snpSet]
#
#   #resamp_eQTL=replicate(NSAMP,list(SNP_info[max_MAF>.05 & ! is.na(covid_A2_pval*covid_B2_pval*covid_C2_pval),.(resampID=sample(ID,sum(is_eQTL),replace=FALSE)),by=cut(AF,seq(0,1,by=.05))]))
#   resamp_eQTL=replicate(NSAMP,list(SNP_info[max_MAF>.05 & ! is.na(covid_A2_pval*covid_B2_pval*covid_C2_pval),.(resampID=sample(ID,sum(is_eQTL),replace=FALSE)),by=.(cut(AF,seq(0,1,by=.05)),cut(DistGene/1000,c(0,1,5,10,20,50,100,Inf)))]))
#   resamp_eQTL=rbindlist(resamp_eQTL,idcol='resamp')
#   resamp_eQTL=rbind(data.table(resamp=0,cut=cut(eQTL_AF,seq(0,1,by=.05)),resampID=snpSet),resamp_eQTL)
#
#   resamp_eQTL[,covid_A2_pval:=SNP_info[match(resampID,ID),covid_A2_pval]]
#   resamp_eQTL[,covid_B2_pval:=SNP_info[match(resampID,ID),covid_B2_pval]]
#   resamp_eQTL[,covid_C2_pval:=SNP_info[match(resampID,ID),covid_C2_pval]]
#
#   resamp_eQTL
# }

cat('OK')
resample_matched_on_MAF_DIST=function(snpSet,NSAMP=1000){
  SNP_metadata=SNP_info[max_MAF>.05 & !is.na(covid_A2_pval*covid_B2_pval*covid_C2_pval) & DistGene<1e5,]
  SNP_metadata[,is_eQTL:=ID%chin%snpSet]
  SNP_metadata[,BIN:=paste(cut(pmin(AF_global,1-AF_global),seq(0,1,by=.01)),cut(DistGene/1000,c(-1,0,1,5,10,20,50,100,Inf)))]
  resamp_eQTL=list()
  resamp_eQTL[[1]]=SNP_metadata[is_eQTL==TRUE,.(BIN,covid_A2_pval,covid_B2_pval,covid_C2_pval)]

    for (i in 1:NSAMP){
      ## print
      tic(paste(i,"/",NSAMP))
      resamp_eQTL[[i+1]] <- SNP_metadata[,.SD[sample(1:.N,sum(is_eQTL),replace=F),.(covid_A2_pval,covid_B2_pval,covid_C2_pval)],by=BIN]
      toc()
    }
    resamp_eQTL=rbindlist(resamp_eQTL,idcol='resamp')
    resamp_eQTL[,resamp:=resamp-1]
  # resamp_eQTL[,covid_A2_pval:=SNP_info[match(resampID,ID),covid_A2_pval]]
  # resamp_eQTL[,covid_B2_pval:=SNP_info[match(resampID,ID),covid_B2_pval]]
  # resamp_eQTL[,covid_C2_pval:=SNP_info[match(resampID,ID),covid_C2_pval]]

  resamp_eQTL
}
cat('OK')

#### step 2: compute Fold enricgments and Pvalues from a set of resampled SNPs.
generate_result_table=function(resamp_eQTL){
  # we add 1 pseudo count to the number of significant pvalues to avoid numerical issues
  resamp_eQTL_PctSignif=as.data.table(resamp_eQTL)[,.(
              Pct_A2_5pct=mean(c(1,covid_A2_pval<0.05),na.rm=T),
              Pct_A2_1pct=mean(c(1,covid_A2_pval<0.01),na.rm=T),
              Pct_A2_01pct=mean(c(1,covid_A2_pval<0.001),na.rm=T),
              Pct_A2_001pct=mean(c(1,covid_A2_pval<0.0001),na.rm=T),
              Pct_B2_5pct=mean(c(1,covid_B2_pval<0.05),na.rm=T),
              Pct_B2_1pct=mean(c(1,covid_B2_pval<0.01),na.rm=T),
              Pct_B2_01pct=mean(c(1,covid_B2_pval<0.001),na.rm=T),
              Pct_B2_001pct=mean(c(1,covid_B2_pval<0.0001),na.rm=T),
              Pct_C2_5pct=mean(c(1,covid_C2_pval<0.05),na.rm=T),
              Pct_C2_1pct=mean(c(1,covid_C2_pval<0.01),na.rm=T),
              Pct_C2_01pct=mean(c(1,covid_C2_pval<0.001),na.rm=T),
              Pct_C2_001pct=mean(c(1,covid_C2_pval<0.0001),na.rm=T)),by=resamp]

cat('OK')
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_A2_5pct>=resamp_eQTL_PctSignif[resamp==0, Pct_A2_5pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_A2_1pct>=resamp_eQTL_PctSignif[resamp==0, Pct_A2_1pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_A2_01pct>=resamp_eQTL_PctSignif[resamp==0, Pct_A2_01pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_A2_001pct>=resamp_eQTL_PctSignif[resamp==0, Pct_A2_001pct])]
  #
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_B2_5pct>=resamp_eQTL_PctSignif[resamp==0, Pct_B2_5pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_B2_1pct>=resamp_eQTL_PctSignif[resamp==0, Pct_B2_1pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_B2_01pct>=resamp_eQTL_PctSignif[resamp==0, Pct_B2_01pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_B2_001pct>=resamp_eQTL_PctSignif[resamp==0, Pct_B2_001pct])]
  #
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_C2_5pct>=resamp_eQTL_PctSignif[resamp==0, Pct_C2_5pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_C2_1pct>=resamp_eQTL_PctSignif[resamp==0, Pct_C2_1pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_C2_01pct>=resamp_eQTL_PctSignif[resamp==0, Pct_C2_01pct])]
  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_C2_001pct>=resamp_eQTL_PctSignif[resamp==0, Pct_C2_001pct])]

  resamp_eQTL_PctSignif_long=melt(resamp_eQTL_PctSignif,id.vars='resamp')
  resamp_eQTL_PctSignif_long[,FE:=value[resamp==0]/value,by=variable]
  resamp_eQTL_PctSignif_long[,threshold:=gsub('Pct_([A-C]2)_(.*)pct','\\2',variable)]
  resamp_eQTL_PctSignif_long[,threshold:=factor(paste0('0.0',threshold),c('0.05','0.01','0.001','0.0001'))]
  resamp_eQTL_PctSignif_long[,GWAS_type:=factor(gsub('Pct_([A-C]2)_(.*)pct','\\1',variable),c('C2','B2','A2'))]
  resamp_eQTL_PctSignif_long[,GWAS_type_full:=case_when(GWAS_type=='C2'~'reported',GWAS_type=='B2'~'hospitalized',GWAS_type=='A2'~'critical')]
  resamp_eQTL_PctSignif_long[,GWAS_type_full:=factor(GWAS_type_full,c('reported','hospitalized','critical'))]

  # pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_eQTL_PctSignif.pdf",EVO_IMMUNO_POP_ZEUS))
  # p <- ggplot(resamp_eQTL_PctSignif_long[resamp!=0,])+geom_violin(aes(x=threshold,y=FE,fill=GWAS_type_full),scale='width') + theme_yann()
  # p <- p + geom_hline(yintercept=1,col='grey') + xlab('p-value threshold')+ ylab('Fold Enrichment of eQTLs in significant associations') + ylim(c(0,10))
  # print(p)
  # dev.off()

  # resamp_eQTL_PctSignif[resamp!=0,mean(Pct_A2_5pct>=resamp_eQTL_PctSignif[resamp==0, Pct_A2_5pct])]


result_table=resamp_eQTL_PctSignif_long[resamp!=0,.(FE_q025=quantile(FE,.025,na.rm=T),
                                              FE=mean(FE,na.rm=T),
                                              FE_q975=quantile(FE,.975,na.rm=T),
                                              pvalue=pnorm(0, mean(log(FE),na.rm=T),
                                                  sd(log(FE),na.rm=T),lower=T)/2),by=.(threshold,GWAS_type_full,GWAS_type)]

  for( i in 1:result_table[,.N]){
    colname=paste0('Pct_',result_table[i,GWAS_type],'_',gsub('0.0','',result_table[i,threshold]),'pct')
    observed=resamp_eQTL_PctSignif[resamp==0,get(colname)]
    result_table[i,P_resamp:=resamp_eQTL_PctSignif[,mean(get(colname)>=observed)]]
  }
  result_table[1:.N]
}

#### step 3: make an enrichment plot from an OR table


plot_enrichment=function(Enrich_table,DIR,pname){
  # Estim_PctSignif=resamp_eQTL_PctSignif_long[resamp!=0,][,.(FE=mean(FE),lower=mean(FE)-2*sd(FE),upper=mean(FE)+2*sd(FE)),by=.(threshold,GWAS_type_full)]
  Estim_PctSignif=Enrich_table[is.finite(FE),][,.(FE=mean(FE),lower=quantile(FE,.025),upper=quantile(FE,.975)),by=.(threshold,GWAS_type_full)]
  color_COVID=c(critical="#DD5555", hospitalized="#E68686", reported="#F1BABA")

  pdf(sprintf("%s/%s",DIR,pname),width=4,height=6)
  p <- ggplot(Enrich_table)+geom_point(aes(x=threshold,y=FE,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5))
  p <- p + geom_errorbar(aes(x=threshold,y=FE,ymin=FE_q025,ymax=FE_q975,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5)) + theme_yann()
  p <- p + geom_hline(yintercept=1,col='grey') + xlab('p-value threshold')+ ylab('Fold Enrichment of eQTLs in significant associations') + ylim(c(0,15))
  p <- p + scale_color_manual(values=color_COVID)
  print(p)
  dev.off()
}

snpSET=snpSets[num==SET_NUM,snps]
set.NAME=snpSets[num==SET_NUM,make.names(unique(set))]


# eQTL snps
snpSET=snpSET[!is.na(SNP_info[match(snpSET,ID),covid_A2_pval*covid_B2_pval*covid_C2_pval])]
cat('OK')
resamp_eQTL=resample_matched_on_MAF_DIST(snpSET,NSAMP=NSAMPLE)
cat('OK')
fwrite(resamp_eQTL,file=sprintf("%s/resamp_%sx_%s_220409.tsv.gz",plot_DIR,NSAMPLE,set.NAME),sep='\t')
resamp_eQTL_Enrichment=generate_result_table(resamp_eQTL)
fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/resamp_COVID19_Enrichment_%sx_%s_220409.tsv.gz",plot_DIR,NSAMPLE,set.NAME),sep='\t')
plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,sprintf('resamp_COVID19_Enrichment_%sx_%s_220409_plot.pdf',NSAMPLE,set.NAME))

q('no')
#
#
# # eQTL snps
# eQTL_snps=eQTL_Signif_both[,unique(snps)]
# eQTL_snps=eQTL_snps[!is.na(SNP_info[match(eQTL_snps,ID),covid_A2_pval*covid_B2_pval*covid_C2_pval])]
# resamp_eQTL=resample_matched_on_MAF(eQTL_snps,NSAMP=NSAMP)
# fwrite(resamp_eQTL,file=sprintf("%s/resamp_1000x_eQTL_FINEMAPPED_220409.tsv.gz",plot_DIR),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_eQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/resamp_COVID19_Enrichment_1000x_eQTL_FINEMAPPED_220409.tsv.gz",plot_DIR),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_eQTL_FINEMAPPED_220409_plot.pdf')
#
# # reQTL snps
# reQTL_snps=reQTL_Signif_both[,unique(snps)]
# resamp_reQTL=resample_matched_on_MAF(reQTL_snps,NSAMP=NSAMP)
# fwrite(resamp_reQTL,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_1000x_reQTL_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_reQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_COVID19_Enrichment_1000x_reQTL_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_reQTL_FINEMAPPED_220409_plot.pdf')
# #
# # COV_reQTL_snps=reQTL_lineage_compare[lfsr_COV<0.01 & lfsr_IAV>0.01,unique(snps)]
# # IAV_reQTL_snps=reQTL_lineage_compare[lfsr_COV>0.01 & lfsr_IAV<0.01,unique(snps)]
# # shared_reQTL_snps=reQTL_lineage_compare[lfsr_COV<0.01 & lfsr_IAV<0.01,unique(snps)]
#
# COV_reQTL_snps=reQTL_lineage_compare[pvalue_COV<0.01 & pvalue_IAV>0.01,unique(snps)]
# IAV_reQTL_snps=reQTL_lineage_compare[pvalue_COV>0.01 & pvalue_IAV<0.01,unique(snps)]
# shared_reQTL_snps=reQTL_lineage_compare[pvalue_COV<0.01 & pvalue_IAV<0.01,unique(snps)]
#
# ###### COV_reQTL
# resamp_reQTL=resample_matched_on_MAF(COV_reQTL_snps,NSAMP=NSAMP)
# fwrite(resamp_reQTL,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_1000x_COV_reQTL_FINEMAPPED_220409_praw.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_reQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_COVID19_Enrichment_1000x_COV_reQTL_FINEMAPPED_220409_praw.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_COV_reQTL_FINEMAPPED_220409_praw_plot.pdf')
#
# ###### IAV_reQTL
# resamp_reQTL=resample_matched_on_MAF(IAV_reQTL_snps,NSAMP=NSAMP)
# fwrite(resamp_reQTL,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_1000x_IAV_reQTL_FINEMAPPED_220409_praw.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_reQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_COVID19_Enrichment_1000x_IAV_reQTL_FINEMAPPED_220409_praw.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_IAV_reQTL_FINEMAPPED_220409_praw_plot.pdf')
#
# ###### shared_reQTL
# resamp_reQTL=resample_matched_on_MAF(shared_reQTL_snps,NSAMP=NSAMP)
# fwrite(resamp_reQTL,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_1000x_shared_reQTL_FINEMAPPED_220409_praw.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_reQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_COVID19_Enrichment_1000x_shared_reQTL_FINEMAPPED_220409_praw.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_shared_reQTL_FINEMAPPED_220409_praw_plot.pdf')
#
# COV_reQTL_snps_Diff=reQTL_lineage_compare[abs(t_diff_COV_IAV)>2 & abs(beta_COV)>abs(beta_IAV),unique(snps)]
# IAV_reQTL_snps_Diff=reQTL_lineage_compare[abs(t_diff_COV_IAV)>2  & abs(beta_IAV)>abs(beta_COV),unique(snps)]
# shared_reQTL_snps_Diff=reQTL_lineage_compare[abs(t_diff)< 2,unique(snps)]
#
# ###### COV_reQTL
# resamp_reQTL=resample_matched_on_MAF(COV_reQTL_snps_Diff,NSAMP=NSAMP)
# fwrite(resamp_reQTL,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_1000x_COV_reQTL_Diff_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_reQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_COVID19_Enrichment_1000x_COV_reQTL_Diff_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_COV_reQTL_Diff_FINEMAPPED_220409_plot.pdf')
#
# ###### IAV_reQTL
# resamp_reQTL=resample_matched_on_MAF(IAV_reQTL_snps_Diff,NSAMP=NSAMP)
# fwrite(resamp_reQTL,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_1000x_IAV_reQTL_Diff_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_reQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_COVID19_Enrichment_1000x_IAV_reQTL_Diff_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_IAV_reQTL_Diff_FINEMAPPED_220409_plot.pdf')
#
# ###### shared_reQTL
# resamp_reQTL=resample_matched_on_MAF(shared_reQTL_snps_Diff,NSAMP=NSAMP)
# fwrite(resamp_reQTL,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_1000x_shared_reQTL_Diff_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# resamp_eQTL_Enrichment=generate_result_table(resamp_reQTL)
# fwrite(resamp_eQTL_Enrichment,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/resamp_COVID19_Enrichment_1000x_shared_reQTL_Diff_FINEMAPPED_220409.tsv.gz",EVO_IMMUNO_POP_ZEUS),sep='\t')
# plot_enrichment(resamp_eQTL_Enrichment,plot_DIR,'resamp_COVID19_Enrichment_1000x_shared_reQTL_Diff_FINEMAPPED_220409_plot.pdf')
