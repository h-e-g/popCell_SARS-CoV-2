# run
# SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
# EQTL_SCRIPT_DIR="09_eQTLmapping"
# sbatch --parsable --mem=30G -o ${SCRIPT_DIR}/${EQTL_SCRIPT_DIR}/log/logfile_selectForColoc_%A.log -J sel_coloc ${SCRIPT_DIR}/00_Rscript.sh ${EQTL_SCRIPT_DIR}/11c_COVIDloci_select.R


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

RUN_NAME="/lineage_condition_logFC___CellPropLineage_SVs_220409"
CIS_DIST=1e5

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",EVO_IMMUNO_POP_ZEUS))


###############################################
########### SNP informations genome wide ######
###############################################

SNP_info=getMap(annotate=TRUE)


###############################################
######         ADD COVID SUMSTATS        ######
###############################################

COVID_A2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
COVID_A2[,SNP:=gsub(':','_',SNP)]
mm=match(SNP_info$posID,COVID_A2$SNP)
SNP_info[,covid_A2_pval:=COVID_A2[mm,all_inv_var_meta_p]]
SNP_info[,covid_A2_beta:=COVID_A2[mm,all_inv_var_meta_beta]]
SNP_info[,covid_A2_sebeta:=COVID_A2[mm,all_inv_var_meta_sebeta]]

COVID_B2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
COVID_B2[,SNP:=gsub(':','_',SNP)]
mm=match(SNP_info$posID,COVID_B2$SNP)
SNP_info[,covid_B2_pval:=COVID_B2[mm,all_inv_var_meta_p]]
SNP_info[,covid_B2_beta:=COVID_B2[mm,all_inv_var_meta_beta]]
SNP_info[,covid_B2_sebeta:=COVID_B2[mm,all_inv_var_meta_sebeta]]

COVID_C2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
COVID_C2[,SNP:=gsub(':','_',SNP)]
mm=match(SNP_info$posID,COVID_C2$SNP)
SNP_info[,covid_C2_pval:=COVID_C2[mm,all_inv_var_meta_p]]
SNP_info[,covid_C2_beta:=COVID_C2[mm,all_inv_var_meta_beta]]
SNP_info[,covid_C2_sebeta:=COVID_C2[mm,all_inv_var_meta_sebeta]]

###############################################
######    define COVID associated loci   ######
###############################################

Putative_COVID_susceptibility_snps=SNP_info[covid_C2_pval<1e-5 & max_MAF>.05,head(.SD[order(covid_C2_pval),],1),by=.(CHROM,round(POS/1e6,0))]
Putative_COVID_susceptibility_snps[,trait:='reported']
Putative_COVID_hospitalization_snps=SNP_info[covid_B2_pval<1e-5 & max_MAF>.05,head(.SD[order(covid_B2_pval),],1),by=.(CHROM,round(POS/1e6,0))]
Putative_COVID_hospitalization_snps[,trait:='hospitalized']
Putative_COVID_ICU_snps=SNP_info[covid_A2_pval<1e-5 & max_MAF>.05,head(.SD[order(covid_A2_pval),],1),by=.(CHROM,round(POS/1e6,0))]
Putative_COVID_ICU_snps[,trait:='critical']

Putative_COVID_snps=rbindlist(list(Putative_COVID_ICU_snps,Putative_COVID_hospitalization_snps,Putative_COVID_susceptibility_snps))
Putative_COVID_snps[,position:=as.numeric(gsub('([0-9]+)_([0-9]+)','\\2',core_posID))]
setnames(Putative_COVID_snps,'chrom','CHROM_NUM')
Putative_COVID_snps=Putative_COVID_snps[order(CHROM,position),]
fwrite(Putative_COVID_snps,file=sprintf('%s/single_cell/resources/references/hgCOVID19/Putative_COVID_snps.tsv.gz',EVO_IMMUNO_POP_ZEUS),sep='\t')

Covid_associated_susceptibility_SNP=SNP_info[covid_C2_pval<1e-5 & max_MAF>.05,]
Covid_associated_susceptibility_SNP[,trait:='reported']
Covid_associated_hospitalization_SNP=SNP_info[covid_B2_pval<1e-5 & max_MAF>.05,]
Covid_associated_hospitalization_SNP[,trait:='hospitalized']
Covid_associated_ICU_SNP=SNP_info[covid_A2_pval<1e-5 & max_MAF>.05,]
Covid_associated_ICU_SNP[,trait:='critical']

Covid_associated_SNP=rbindlist(list(Covid_associated_susceptibility_SNP,Covid_associated_hospitalization_SNP,Covid_associated_ICU_SNP))
setnames(Covid_associated_SNP,'chrom','CHROM_NUM')

fwrite(Covid_associated_SNP,file=sprintf('%s/single_cell/resources/references/hgCOVID19/Covid_associated_SNP_p1e5_MAF5pct.tsv.gz',EVO_IMMUNO_POP_ZEUS),sep='\t')

################################################################################
######  overlap COVID associated loci with eQTL/reQTL with 100kb window   ######
################################################################################


library(GenomicRanges)
GR_covid_susceptibility=makeGRangesFromDataFrame(Covid_associated_SNP[trait=='reported',], seqnames='CHROM', start.field='POS', end.field='POS', strand='*')
GR_covid_susceptibility=union(GR_covid_susceptibility,flank(GR_covid_susceptibility,1e5,both=TRUE))
GR_covid_susceptibility$TRAIT='reported'

GR_covid_hospitalization=makeGRangesFromDataFrame(Covid_associated_SNP[trait=='hospitalized',], seqnames='CHROM', start.field='POS', end.field='POS', strand='*')
GR_covid_hospitalization=union(GR_covid_hospitalization,flank(GR_covid_hospitalization,1e5,both=TRUE))
GR_covid_hospitalization$TRAIT='hospitalized'

GR_covid_ICU=makeGRangesFromDataFrame(Covid_associated_SNP[trait=='critical',], seqnames='CHROM', start.field='POS', end.field='POS', strand='*')
GR_covid_ICU=union(GR_covid_ICU,flank(GR_covid_ICU,1e5,both=TRUE))
GR_covid_ICU$TRAIT='critical'



GR_reQTL=makeGRangesFromDataFrame(SNP_info[ID%chin%unique(reQTL_Signif_both$snps),.(ID,CHROM,POS)], seqnames='CHROM', start.field='POS', end.field='POS', strand='*',keep.extra=TRUE)
GR_eQTL=makeGRangesFromDataFrame(SNP_info[ID%chin%unique(eQTL_Signif_both$snps),.(ID,CHROM,POS)], seqnames='CHROM', start.field='POS', end.field='POS', strand='*',keep.extra=TRUE)

GWAS_overlaps=list()
oo=findOverlaps(GR_eQTL,GR_covid_susceptibility)
GWAS_overlaps[['eQTL_reported']]=data.table(CHR=as.character(seqnames(GR_eQTL))[queryHits(oo)],snps=GR_eQTL$ID[queryHits(oo)],Trait='reported',type='eQTL')
oo=findOverlaps(GR_eQTL,GR_covid_hospitalization)
GWAS_overlaps[['eQTL_hospitalized']]=data.table(CHR=as.character(seqnames(GR_eQTL))[queryHits(oo)],snps=GR_eQTL$ID[queryHits(oo)],Trait='hospitalized',type='eQTL')
oo=findOverlaps(GR_eQTL,GR_covid_ICU)
GWAS_overlaps[['eQTL_critical']]=data.table(CHR=as.character(seqnames(GR_eQTL))[queryHits(oo)],snps=GR_eQTL$ID[queryHits(oo)],Trait='critical',type='eQTL')

oo=findOverlaps(GR_reQTL,GR_covid_susceptibility)
GWAS_overlaps[['reQTL_reported']]=data.table(CHR=as.character(seqnames(GR_reQTL))[queryHits(oo)],snps=GR_reQTL$ID[queryHits(oo)],Trait='reported',type='reQTL')
oo=findOverlaps(GR_reQTL,GR_covid_hospitalization)
GWAS_overlaps[['reQTL_hospitalized']]=data.table(CHR=as.character(seqnames(GR_reQTL))[queryHits(oo)],snps=GR_reQTL$ID[queryHits(oo)],Trait='hospitalized',type='reQTL')
oo=findOverlaps(GR_reQTL,GR_covid_ICU)
GWAS_overlaps[['reQTL_critical']]=data.table(CHR=as.character(seqnames(GR_reQTL))[queryHits(oo)],snps=GR_reQTL$ID[queryHits(oo)],Trait='critical',type='reQTL')

GWAS_overlaps=rbindlist(GWAS_overlaps)
GWAS_overlaps[,CHR:=gsub('chr','',CHR)]

GWAS_overlaps_reQTL=merge(reQTL_Signif_both[,.(snps,gene,gene_name,celltype,state)],GWAS_overlaps[type=='reQTL',])
GWAS_overlaps_eQTL=merge(eQTL_Signif_both[,.(snps,gene,gene_name,celltype,state)],GWAS_overlaps[type=='eQTL',])
GWAS_overlaps_final=rbind(GWAS_overlaps_reQTL,GWAS_overlaps_eQTL)
GWAS_overlaps_final[,RUN:=case_when(type=='eQTL' & celltype%in% lineage_5 ~ RUN_EQTL_LINEAGE,
                                    type=='eQTL' & celltype%in% celltype_22 ~ RUN_EQTL_CELLTYPE,
                                    type=='reQTL' & celltype%in% lineage_5 ~ RUN_REQTL_LINEAGE,
                                    type=='reQTL' & celltype%in% celltype_22 ~ RUN_REQTL_CELLTYPE)]

#### all overlaps
# 1 line per snps, gene, celltype, condition, QTL type
fwrite(GWAS_overlaps_final,file=sprintf('%s/overlap_COVID_1e5_within_100kb_of_eQTL_and_reQTL_peak.txt',OUT_DIR),sep='\t')


# 1 line per snps, gene, celltype, QTL type
fwrite(unique(GWAS_overlaps_final[,.(snps,gene,gene_name,CHR,type,RUN)]),file=sprintf('%s/overlap_COVID_1e5_within_100kb_of_eQTL_and_reQTL_peak_no_Condition.txt',OUT_DIR),sep='\t')

# 1 line per snps, gene
fwrite(unique(GWAS_overlaps_final[,.(snps,gene,gene_name,CHR)]),file=sprintf('%s/overlap_COVID_1e5_within_100kb_of_eQTL_and_reQTL_peak_snp_gene_only.txt',OUT_DIR),sep='\t')

################################################################################
######  overlap COVID associated loci with eQTL/reQTL with 100kb window   ######
################################################################################

#
# LD_FULL=list()
# for(CHR in 11:22){
#   for(POP in c('CEU','YRI','CHS')){
#   cat(CHR,POP,sep='\t')
#   LD_CHR_POP=fread(sprintf('%s/single_cell/resources/references/1kg_ldmap/%s_chr%s_biSNPs_nodups_nomono_hg38.ld.gz',EVO_IMMUNO_POP_ZEUS,POP,CHR))
#   LD_CHR_POP[,core_posID_A:=paste(CHR_A,BP_A,sep='_')]
#   LD_CHR_POP[,core_posID_B:=paste(CHR_B,BP_B,sep='_')]
#   LD_CHR_POP=rbindlist(list(SNP_info[CHROM==paste0("chr",CHR),.(core_posID=core_posID,core_posID_linked=core_posID,r2=1,pop=POP,CHR=CHR)],
#                       LD_CHR_POP[,.(core_posID=core_posID_A,core_posID_linked=core_posID_B,r2=R2,pop=POP,CHR=CHR_A)],
#                       LD_CHR_POP[,.(core_posID=core_posID_B,core_posID_linked=core_posID_A,r2=R2,pop=POP,CHR=CHR_A)]))
#   LD_FULL[[paste(CHR,POP)]]=LD_CHR_POP
#   }
# }
# LD_FULL=rbindlist(LD_FULL)
#
# Putative_COVID_snps_LD80=merge(Putative_COVID_snps,LD_FULL,by='core_posID')
# fwrite(Putative_COVID_snps_LD80,file=sprintf('%s/single_cell/resources/references/hgCOVID19/Putative_COVID_snps_LD80.tsv.gz',EVO_IMMUNO_POP_ZEUS),sep='\t')

###### merge the set of SNPs extracted (in LD with a 10e-5 COVID hit) with 95% CI of eQTLs and reQTLs

# eQTL_CI_lineage[,core_posID:=SNP_info[match(eQTL_CI_lineage$snps,SNP_info$ID),core_posID]]
# reQTL_CI_lineage[,core_posID:=SNP_info[match(reQTL_CI_lineage$snps,SNP_info$ID),core_posID]]
#
# overlaps=list(merge(Putative_COVID_snps_LD80,eQTL_CI_lineage,by.x='core_posID_linked',by.y='core_posID'),
#             merge(Putative_COVID_snps_LD80,eQTL_CI_lineage,by.x='core_posID_linked',by.y='core_posID'))
# overlaps=rbindlist(overlaps)
# overlaps[,pop:=factor(pop,c('YRI','CEU','CHS'))]
#
# overlap_unique=overlaps[order(Symbol,trait,type,-pip_top_component,pop),head(.SD,1),by=.(Symbol,trait,type)]
# setnames(overlap_unique,c('core_posID','core_posID_linked'),c('COVID_peakSNP','eQTL_peakSNP'))

################################################################################################################
############## sets of genes that are worth testing for colocalization with COVID19 ############################
################################################################################################################

# fwrite(overlap_unique,file=sprintf('%s/%s/dist_%s/overlap_COVID_1e5_with_eQTL_and_reQTL_95CI.txt',OUT_DIR,RUN_NAME,CIS_DIST_TEXT),sep='\t')
# overlap_unique=fread(sprintf('%s/%s/dist_%s/overlap_COVID_1e5_with_eQTL_and_reQTL_95CI.txt',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))

################################################################################################################
############## sets of genes that are worth testing for colocalization with COVID19 ############################
################################################################################################################

# library(ggrastr)
# overlap_unique_gene=unique(overlap_unique[,.(Symbol,gene,trait,type,CHROM,posID,COVID_peakSNP,eQTL_peakSNP,covid_A2_pval,covid_B2_pval,covid_C2_pval,cellstate,snps,pvalue)])
# WINDOW_SIZE=2e5
# RUN_NAME=RUN_EQTL_LINEAGE
#
# FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID',EVO_IMMUNO_POP_ZEUS)
# for (i in 1:overlap_unique_gene[,.N]){
#   CHR=overlap_unique_gene[i,CHROM]
#   position=overlap_unique_gene[i,as.numeric(gsub('[0-9]+_([0-9]+)','\\1',eQTL_peakSNP))]
#   GENE=overlap_unique_gene[i,gene]
#   SYMBOL=overlap_unique_gene[i,Symbol]
#   TRAIT=overlap_unique_gene[i,trait]
#   myCT_S=overlap_unique_gene[i,cellstate]
#   SNP=overlap_unique_gene[i,eQTL_peakSNP]
#   myCELLTYPE__STATE=overlap_unique_gene[i,cellstate]
# #   COVID_stats=SNP_info[CHROM==CHR & POS> position-WINDOW_SIZE & POS< position+WINDOW_SIZE,.(ID,CHROM,POS, REF, ALT,covid_A2_pval, covid_B2_pval, covid_C2_pval, covid_A2_beta,  covid_B2_beta,covid_C2_beta,covid_A2_sebeta,  covid_B2_sebeta,covid_C2_sebeta)]
# #   COVID_stats=melt(COVID_stats,id.vars=c('ID','CHROM','POS','REF','ALT'))
# #   COVID_stats[,stat:=gsub('covid_([ABC]2)_(beta|pval|sebeta)','\\2',variable)]
# #   COVID_stats[,GWAS_type:=gsub('covid_([ABC]2)_(beta|pval|sebeta)','\\1',variable)]
# #   COVID_stats=dcast(COVID_stats, ID+CHROM+POS+REF+ALT+GWAS_type~stat)
# #   COVID_stats[,GWAS_type_full:=case_when(GWAS_type=='A2'~'critical',
# #                                         GWAS_type=='B2'~'hospitalized',
# #                                         GWAS_type=='C2'~'reported',
# #                                         TRUE~'NA')]
# #
# #   p1 <- ggplot(COVID_stats,aes(x=POS,y=-log10(pval),col=GWAS_type_full,fill=GWAS_type_full,shape=ifelse(beta>0,'pos','neg')))
# #   p1 <- p1 + rasterize(geom_point(alpha=.7),dpi=400) + facet_grid(GWAS_type_full~1) + scale_fill_manual(values=color_COVID)
# #   p1 <- p1 + scale_color_manual(values=color_COVID) + scale_shape_manual(values=c("pos"=24,"neg"=25))
# #   p1 <- p1 + geom_vline(col='grey',xintercept=position)
# #
# #   RUN_NAME='lineage_condition___CellPropLineage_SVs_220409'
# #   cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')
# #
# #   eQTL_stats=list()
# #   for (CELLTYPE__STATE in cellstates){
# #     cat(CELLTYPE__STATE,'')
# #     cmd=sprintf('gunzip -c %s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/%s/%s/eQTL_ALL_%s_assoc.txt.gz | grep %s',EVO_IMMUNO_POP_ZEUS,RUN_NAME,CELLTYPE__STATE,CHR,GENE)
# #     eQTL_stats[[CELLTYPE__STATE]]=fread(cmd)
# #   }
# # eQTL_stats=rbindlist(eQTL_stats,idcol='cellstate')
# # setnames(eQTL_stats,c('V1','V2','V3','V4','V5','V6','V7','V8','V9'),c('snps','gene','statistic','pval','beta','celltype','R2','sebeta','CisDist'),skip_absent=TRUE)
# # #setnames(eQTL_stats,c('pvalue','se'),c('pval','sebeta'),skip_absent=TRUE)
# #
# # eQTL_stats[,state:=gsub('(.*)__(.*)','\\2',cellstate)]
# # mm=match(eQTL_stats$snps,SNP_info$ID)
# #   eQTL_stats[,POS:=SNP_info[mm,POS]]
# #   eQTL_stats=eQTL_stats[POS>position-WINDOW_SIZE & POS<position+WINDOW_SIZE,]
# #   p2 <- ggplot(eQTL_stats,aes(x=POS,y=-log10(pval),col=state,fill=state,shape=ifelse(beta>0,'pos','neg')))+facet_grid(celltype~1)
# #   p2 <- p2 + rasterize(geom_point(alpha=.7),dpi=400)+ scale_fill_manual(values=color_conditions)
# #   p2 <- p2 + scale_color_manual(values=color_conditions) + scale_shape_manual(values=c('pos'=24,'neg'=25))
# #   p2 <- p2 + geom_vline(col='grey',xintercept=position)
#
# dir.create(sprintf('%s/%s_%s_%s/',FIGURE_DIR,SYMBOL,TRAIT,SNP))
# #fwrite(COVID_stats,file=sprintf('%s/%s_%s_%s/%s_%s_%s_COVID.tsv.gz',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP),sep='\t')
# COVID_stats=fread(sprintf('%s/%s_%s_%s/%s_%s_%s_COVID.tsv.gz',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#   # saveRDS(p1,file=sprintf('%s/%s_%s_%s/%s_%s_%s_COVID_plot.RDS',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#   # pdf(sprintf('%s/%s_%s_%s/%s_%s_%s_COVID.pdf',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#   # print(p1)
#   # dev.off()
#
# #fwrite(eQTL_stats,file=sprintf('%s/%s_%s_%s/%s_%s_%s_eQTL.tsv.gz',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP),sep='\t')
# eQTL_stats=fread(sprintf('%s/%s_%s_%s/%s_%s_%s_eQTL.tsv.gz',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#
#   # saveRDS(p2,file=sprintf('%s/%s_%s_%s/%s_%s_%s_eQTL_plot.RDS',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#   # pdf(sprintf('%s/%s_%s_%s/%s_%s_%s_eQTL.pdf',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#   # print(p2)
#   # dev.off()
#
#   snp_set_eQTL=eQTL_stats[cellstate==myCT_S,snps]
#   snp_set_COVID=COVID_stats[GWAS_type_full==TRAIT & !is.na(beta),ID]
#   snp_set=intersect(snp_set_COVID,snp_set_eQTL)
#
#   MinCell_perCOND=500
#   keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
#
#   Genotypes=queryRange(CHR,position-WINDOW_SIZE,position+WINDOW_SIZE)
#   Genotypes=melt(Genotypes,id.vars=c("ID"),variable.name='IID')[IID%in%keptIID & ID%chin%snp_set,]
#   setkey(Genotypes,ID,IID)
#   Genotypes[,value:=as.numeric(value)]
#
#   Geno_mat=dcast(Genotypes,IID~ID)
#   Geno_res=apply(as.matrix(Geno_mat[,-'IID']),2,function(x){lm(x~substr(Geno_mat$IID,1,3))$res})
#   LDMAP=cor(Geno_res)
#
#   eQTL_stats_locus=eQTL_stats[cellstate==myCT_S & snps%chin%snp_set,][order(POS),]
#   COVID_stats_locus=COVID_stats[GWAS_type_full==TRAIT & ID%chin%snp_set,][order(POS),]
#   compare_beta=merge(COVID_stats,eQTL_stats,by.x=c('ID','POS'),by.y=c('snps','POS'),suffix=c('.covid','.eQTL'),allow.cartesian=TRUE)
#
#   snp_set=eQTL_stats_locus[,snps]
#
#   data_eQTL=list(beta=setNames(eQTL_stats_locus[,beta],eQTL_stats_locus[,snps]),
#                   varbeta=setNames(eQTL_stats_locus[,sebeta^2],eQTL_stats_locus[,snps]),
#                   snp=eQTL_stats_locus[,snps],
#                   position=eQTL_stats_locus[,POS],
#                   LD=LDMAP[snp_set,snp_set],
#                   type='quant',
#                   sdY=1)
#
#   data_COVID=list(beta=setNames(COVID_stats_locus[,beta],COVID_stats_locus[,ID]),
#                   varbeta=setNames(COVID_stats_locus[,sebeta^2],COVID_stats_locus[,ID]),
#                   snp=COVID_stats_locus[,ID],
#                   position=COVID_stats_locus[,POS],
#                   LD=LDMAP[snp_set,snp_set],
#                   type='quant',
#                   sdY=1)
#   #### For a chosen eQTL/condition: show a direct comparison of betas, and a plot colored by R2 with the peaks eQTL SNP
#   # S_eQTL=try(runsusie(data_eQTL,repeat_until_convergence=FALSE))
#   # S_COVID=try(runsusie(data_eQTL,repeat_until_convergence=FALSE))
#   # if(class(S_eQTL)=='try-error' | class(S_COVID)=='try-error'){
#   #   coloc.res=list(summary=data.table(nsnps=length(snp_set), hit1=NA, hit2=NA, PP.H0.abf=1, PP.H1.abf=0, PP.H2.abf=0, PP.H3.abf=0,PP.H4.abf=0,idx1=NA, idx2=NA,converged=FALSE))
#   #   }else{
#   #     coloc.res=coloc.susie(S_eQTL,S_COVID)
#   #     if(is.null(coloc.res$summary)){
#   #       coloc.res=list(summary=data.table(nsnps=length(snp_set), hit1=NA, hit2=NA, PP.H0.abf=1, PP.H1.abf=0, PP.H2.abf=0, PP.H3.abf=0,PP.H4.abf=0,idx1=NA, idx2=NA,converged=TRUE))
#   #       }else{
#   #         coloc.res$summary$converged=TRUE
#   #       }
#   #     }
#   # fwrite(coloc.res$summary,file=sprintf('%s/%s_%s_%s/%s_%s_%s_%s_coloc.tsv',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP,myCT_S),sep='\t')
#   #### run another colocalization technique
#   coloc.res=coloc.signals(data_eQTL,data_COVID,p12=1e-5)
#   fwrite(coloc.res$summary,file=sprintf('%s/%s_%s_%s/%s_%s_%s_%s_coloc_signal.tsv',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP,myCT_S),sep='\t')
#   #p3 <- ggplot(compare_beta[pval.covid<1e-3 | pval.eQTL<1e-3,],aes(x=beta.eQTL,y=beta.covid,col=state,fill=state))+facet_grid(celltype~GWAS_type_full)
#   # p3 <- ggplot(compare_beta[GWAS_type_full==TRAIT,],aes(x=beta.eQTL/sebeta.eQTL,y=beta.covid/sebeta.covid,col=state,fill=state))+facet_grid(celltype~state)
#   # p3 <- p3 + rasterize(geom_point(alpha=.7),dpi=400)+ scale_fill_manual(values=color_conditions)
#   # p3 <- p3 + scale_color_manual(values=color_conditions) + scale_shape_manual(values=c('pos'=24,'neg'=25))
#   # p3 <- p3 + geom_smooth(method='lm',formula=y~0+x,col='black',fill='black') + geom_hline(col='grey',yintercept=0)+ geom_vline(col='grey',xintercept=0)
#   # p3 <- p3 + geom_hline(col='lightgrey',yintercept=3,linetype="dashed")+ geom_vline(col='grey',xintercept=3,linetype="dashed")
#   # p3 <- p3 + geom_hline(col='lightgrey',yintercept=-3,linetype="dashed")+ geom_vline(col='grey',xintercept=-3,linetype="dashed")
#   #   saveRDS(p3,file=sprintf('%s/%s_%s_%s/%s_%s_%s_COVID_eQTL_betacompare_plot.RDS',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#   #   pdf(sprintf('%s/%s_%s_%s/%s_%s_%s_COVID_eQTL_betacompare_plot.pdf',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
#   #   print(p3)
#   #   dev.off()
# }
#
# ############# aggregate results
# tested_loci=dir(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID",EVO_IMMUNO_POP_ZEUS))
# coloc.res=list()
# coloc.signal.res=list()
# for (i in tested_loci){
#   file=dir(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID/%s",EVO_IMMUNO_POP_ZEUS,i),pattern='coloc.tsv')
#   file.signal=dir(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID/%s",EVO_IMMUNO_POP_ZEUS,i),pattern='coloc_signal.tsv')
#   for (j in file){
#     coloc.res[[paste(i,j,sep='/')]]=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID/%s/%s",EVO_IMMUNO_POP_ZEUS,i,j))
#   }
#   for(j in file.signal){
#     coloc.signal.res[[paste(i,j,sep='/')]]=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID/%s/%s",EVO_IMMUNO_POP_ZEUS,i,j))
#   }
# }
# coloc.res=rbindlist(coloc.res,idcol='locus',fill=TRUE)
# coloc.res[,.N,by=.(is.na(hit1),converged,PP.H4.abf>.5)]
#
# # is.na converged PP.H4.abf  N
# # 1:  TRUE      TRUE     FALSE 24
# # 2:  TRUE     FALSE     FALSE 22
# # 3: FALSE      TRUE     FALSE 19
# # 4: FALSE      TRUE      TRUE  5
#
# coloc.signal.res=rbindlist(coloc.signal.res,idcol='locus',fill=TRUE)
# coloc.signal.res[,.N,keyby=.(PP.H4.abf>.2, PP.H4.abf>.5, PP.H4.abf>.8)]
# fwrite(coloc.signal.res,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID/colocalization_signals.txt",EVO_IMMUNO_POP_ZEUS),sep='\t')

##  $ beta    : Named num [1:50] 0.288 0.3 0.334 0.444 0.494 ...
##   ..- attr(*, "names")= chr [1:50] "s1" "s2" "s3" "s4" ...
##  $ varbeta : Named num [1:50] 0.00681 0.0105 0.00733 0.00591 0.01514 ...
##   ..- attr(*, "names")= chr [1:50] "s1" "s2" "s3" "s4" ...
##  $ snp     : chr [1:50] "s1" "s2" "s3" "s4" ...
##  $ position: int [1:50] 1 2 3 4 5 6 7 8 9 10 ...
##  $ type    : chr "quant"
##  $ sdY     : num 1.12
