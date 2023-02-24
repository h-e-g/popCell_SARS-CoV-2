################################################################################
################################################################################
# File name: load_eQTLs.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: load eQTL data
# AIM: create (r)eQTL related objects for easy manipulation
# same script is used for reQTLs
# Effector script
#################################################################
#################################################################
# The following objects are created:

### 1 - RUN_EQTL_LINEAGE, RUN_REQTL_LINEAGE, RUN_EQTL_CELLTYPE, RUN_REQTL_CELLTYPE
# define the RUN folder where to retrieve information on the various eQTL mapping levels

### 2 -  Feature_annot
# annotation of genes used in the data set (chr, position, etc...)

### 3 -  lineage_5, celltype_22, celltype_2lineage
# list of 5 lineages/22 celltypes used and their correspondance.

### 4 - eQTL_Signif_both / reQTL_Signif_both
# list of genome wide significant eQTL/reQTLs at lineage & cell type level.
# celltype and lineage eQTLs have been merged so eQTLs with overlapping 95 CI are represented by the same SNP & a single SNP is provided for each independent gene x SNP association.
# one line per gene x SNP x celltype/lineage x condition

#### 5 - eQTL_Stats_lineage, eQTL_Stats_celltype, reQTL_Stats_lineage, reQTL_Stats_celltype
# list of all summary statistics for significant eQTL/reQTLs at lineage & cell type level.
# all SNP/gene pairs present in QTL_Signif_both / reQTL_Signif_both are considered
# one line per ,gene x SNP x celltype/lineage x condition : all celltypes/lineages & conditions are present for each eQTL regardless of significance.
# only lines with pvalue <0.01 should be considered considered as eQTLs in that tissue (of lfsr<.01 if using mash to call eQTLs)

#### 6 - eQTL_celltype_compare, eQTL_lineage_compare, reQTL_celltype_compare, reQTL_lineage_compare
# same as 5, but differnt conditions are split into distinct columns top allow caharcterization of condition specific eQTLs/reQTLs
# one line per ,gene x SNP x celltype/lineage : all celltypes/lineages  are present for each eQTL regardless of significance.
# only lines with pvalue <0.01 for at least one condition be considered considered as eQTLs in that tissue (of lfsr<.01 if using mash to call eQTLs)
# t_diff is a Zvalue that can be used to assess significance of differences in effect size between conditions.

#### 7 - eQTL_CI_lineage, eQTL_CI_celltype, reQTL_CI_lineage, reQTL_CI_celltype
# 95% CI for all eQTLs and reQTL that pass a lbf>3 threshold in each celltype/lineage (regardless of FDR)

#### 8 - snpSets
# list of eQTLs/reQTLs to consider for enrichment analyses (celltype/condition specific)

#### 9 - eQTL_genotypes
# file containing eQTLs genotypes across 222 indiivduals for all 13839 eQTL/reQTL snps

#################################################################
#################################################################
#################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
EQTL_DIR = "3__eQTL_mapping/SumStats"

RUN_ID='RUN_1'
CIS_DIST_TEXT='100kb'
RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_RUN1"
RUN_REQTL_LINEAGE="lineage_condition_logFC__CellPropLineage_SVs_RUN1"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_RUN1"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__CellPropLineage_SVs_RUN1"

##### compare reQTL effect sizes
celltype_2lineage=fread('1__transcriptome_processing/data/lineage_celltype.tsv')
lineage_5=unique(celltype_2lineage$lineage)
celltype_22=unique(celltype_2lineage$celltype)

###############################################################
########### extract eQTLs (lineage or celltype)  ##############
###############################################################
tic('extracting eQTL information')
eQTL_Signif_both=fread(sprintf('%s/%s_eQTL/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR.txt.gz',EQTL_DIR,RUN_ID))

###############################################################
########### extract lineage level information on eQTLs  #######
###############################################################

eQTL_Stats_lineage=fread(sprintf('%s/%s_eQTL/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_lineage_allCHR_FDR.txt.gz',EQTL_DIR,RUN_ID))
setnames(eQTL_Stats_lineage,'gene_id','gene')
eQTL_Stats_lineage=merge(eQTL_Stats_lineage[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2)],
                              eQTL_Signif_both[celltype%chin%lineage_5,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)

Number_lineage_where_signif=eQTL_Stats_lineage[,.(Number_lineage_where_signif_praw1=length(unique(celltype[pvalue<0.01]))),by=.(gene,snps)]
eQTL_Stats_lineage=merge(eQTL_Stats_lineage,Number_lineage_where_signif,all.x=TRUE)
eQTL_Stats_lineage[,cellstate:=paste(celltype,state,sep='__')]
rm(Number_lineage_where_signif)


eQTL_lineage_compare=dcast(eQTL_Stats_lineage,gene+gene_name+snps+celltype~state,value.var=list('beta','pvalue','se'))
eQTL_lineage_compare[,t_diff_COV_IAV:=(beta_COV-beta_IAV)/sqrt(se_COV^2+se_IAV^2)]
eQTL_lineage_compare[,p_diff_COV_IAV:=2*pnorm(abs(t_diff_COV_IAV),lower=F)]
eQTL_lineage_compare[,fdr_diff_COV_IAV:=p.adjust(p_diff_COV_IAV,'fdr')]

eQTL_lineage_compare[,t_diff_COV_NS:=(beta_COV-beta_NS)/sqrt(se_COV^2+se_NS^2)]
eQTL_lineage_compare[,p_diff_COV_NS:=2*pnorm(abs(t_diff_COV_NS),lower=F)]
eQTL_lineage_compare[,fdr_diff_COV_NS:=p.adjust(p_diff_COV_NS,'fdr')]

eQTL_lineage_compare[,t_diff_IAV_NS:=(beta_IAV-beta_NS)/sqrt(se_IAV^2+se_NS^2)]
eQTL_lineage_compare[,p_diff_IAV_NS:=2*pnorm(abs(t_diff_IAV_NS),lower=F)]
eQTL_lineage_compare[,fdr_diff_IAV_NS:=p.adjust(p_diff_IAV_NS,'fdr')]

###############################################################
########### extract celltype level information on eQTLs #######
###############################################################

eQTL_Stats_celltype=fread(sprintf('%s/%s_eQTL/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_celltype_allCHR_FDR.txt.gz',EQTL_DIR,RUN_ID))
setnames(eQTL_Stats_celltype,'gene_id','gene')
eQTL_Stats_celltype=merge(eQTL_Stats_celltype[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2)],
                              eQTL_Signif_both[celltype%chin%celltype_22,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)

Number_celltype_where_signif=eQTL_Stats_celltype[,.(Number_celltype_where_signif_praw1=length(unique(celltype[pvalue<0.01])))
                                                                ),by=.(gene,snps)]
eQTL_Stats_celltype=merge(eQTL_Stats_celltype,Number_celltype_where_signif,all.x=TRUE)
eQTL_Stats_celltype[,cellstate:=paste(celltype,state,sep='__')]
rm(Number_celltype_where_signif_padj)


eQTL_celltype_compare=dcast(eQTL_Stats_celltype,gene+gene_name+snps+celltype~state,value.var=list('beta','pvalue','se'))
eQTL_celltype_compare[,t_diff_COV_IAV:=(beta_COV-beta_IAV)/sqrt(se_COV^2+se_IAV^2)]
eQTL_celltype_compare[,p_diff_COV_IAV:=2*pnorm(abs(t_diff_COV_IAV),lower=F)]
eQTL_celltype_compare[,fdr_diff_COV_IAV:=p.adjust(p_diff_COV_IAV,'fdr')]

eQTL_celltype_compare[,t_diff_COV_NS:=(beta_COV-beta_NS)/sqrt(se_COV^2+se_NS^2)]
eQTL_celltype_compare[,p_diff_COV_NS:=2*pnorm(abs(t_diff_COV_NS),lower=F)]
eQTL_celltype_compare[,fdr_diff_COV_NS:=p.adjust(p_diff_COV_NS,'fdr')]

eQTL_celltype_compare[,t_diff_IAV_NS:=(beta_IAV-beta_NS)/sqrt(se_IAV^2+se_NS^2)]
eQTL_celltype_compare[,p_diff_IAV_NS:=2*pnorm(abs(t_diff_IAV_NS),lower=F)]
eQTL_celltype_compare[,fdr_diff_IAV_NS:=p.adjust(p_diff_IAV_NS,'fdr')]
###############################################################
toc()

#########################################################
########### extract reQTLs (lineage only) ###############
#########################################################
tic('extracting reQTL information')
reQTL_Signif_both=fread(sprintf('%s/%s_reQTL/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR.txt.gz',EQTL_DIR,RUN_ID))

###############################################################
########### extract lineage level information on reQTLs #######
###############################################################

# extract stats for both cell type and lineage eQTLs
reQTL_Stats_lineage=fread(sprintf('%s/%s_reQTL/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_lineage_allCHR_FDR.txt.gz',EQTL_DIR,RUN_ID))
setnames(reQTL_Stats_lineage,'gene_id','gene')
reQTL_Stats_lineage=merge(reQTL_Stats_lineage[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2)],
                              reQTL_Signif_both[celltype%chin%lineage_5,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
#
Number_lineage_where_signif=reQTL_Stats_lineage[,.(Number_lineage_where_signif_praw1=length(unique(celltype[pvalue<0.01]))),by=.(gene,snps)]
reQTL_Stats_lineage=merge(reQTL_Stats_lineage,Number_lineage_where_signif,all.x=TRUE)
reQTL_Stats_lineage[,cellstate:=paste(celltype,state,sep='__')]

reQTL_lineage_compare=dcast(reQTL_Stats_lineage,gene+gene_name+snps+celltype~state,value.var=list('beta','pvalue','se'))
reQTL_lineage_compare[,t_diff_COV_IAV:=(beta_COV-beta_IAV)/sqrt(se_COV^2+se_IAV^2)]
reQTL_lineage_compare[,p_diff_COV_IAV:=2*pnorm(abs(t_diff_COV_IAV),lower=F)]
reQTL_lineage_compare[,fdr_diff_COV_IAV:=p.adjust(p_diff_COV_IAV,'fdr')]


###############################################################
########### extract celltype level information on reQTLs ######
###############################################################
reQTL_Stats_celltype=fread(sprintf('%s/%s_reQTL/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_celltype_allCHR_FDR.txt.gz',EQTL_DIR,RUN_ID))
setnames(reQTL_Stats_celltype,'gene_id','gene')
reQTL_Stats_celltype=merge(reQTL_Stats_celltype[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2)],
                              reQTL_Signif_both[celltype%chin%celltype_22,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)

Number_celltype_where_signif=reQTL_Stats_celltype[,.(Number_celltype_where_signif_praw1=length(unique(celltype[pvalue<0.01]))),by=.(gene,snps)]
reQTL_Stats_celltype=merge(reQTL_Stats_celltype,Number_celltype_where_signif,all.x=TRUE)
reQTL_Stats_celltype[,cellstate:=paste(celltype,state,sep='__')]

reQTL_celltype_compare=dcast(reQTL_Stats_celltype,gene+gene_name+snps+celltype~state,value.var=list('beta','pvalue','se'))
reQTL_celltype_compare[,t_diff_COV_IAV:=(beta_COV-beta_IAV)/sqrt(se_COV^2+se_IAV^2)]
reQTL_celltype_compare[,p_diff_COV_IAV:=2*pnorm(abs(t_diff_COV_IAV),lower=F)]
reQTL_celltype_compare[,fdr_diff_COV_IAV:=p.adjust(p_diff_COV_IAV,'fdr')]

#############################################################################################################################################
########################## define the set of all significant eQTLs & reQTLs for comparison of beta between COV and IAV ######################
#############################################################################################################################################

Signif_reQTL_compare=merge(reQTL_Signif_both,rbind(reQTL_celltype_compare,reQTL_lineage_compare),all.x=TRUE,by=c('snps','gene','gene_name','celltype'))
Signif_reQTL_compare[,lineage:=ifelse(celltype %chin% celltype_22,celltype_2lineage[match(Signif_reQTL_compare$celltype,celltype),lineage],celltype)]

Signif_eQTL_compare=merge(eQTL_Signif_both,rbind(eQTL_celltype_compare,eQTL_lineage_compare),all.x=TRUE,by=c('snps','gene','gene_name','celltype'))
Signif_eQTL_compare[,lineage:=ifelse(celltype %chin% celltype_22,celltype_2lineage[match(Signif_eQTL_compare$celltype,celltype),lineage],celltype)]

################################################
########### extract eQTL 95% CIs ##############
################################################
tic('extracting 95% credible intervals')
# eQTL lineage
RUN_NAME=RUN_EQTL_LINEAGE
eQTL_CI_lineage=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',EQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
eQTL_CI_lineage=merge(eQTL_CI_lineage, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
eQTL_CI_lineage[,type:='eQTL']
eQTL_CI_lineage[,state:=gsub('.*__(.*)','\\2',cellstate)]

# eQTL celltype
RUN_NAME=RUN_EQTL_CELLTYPE
eQTL_CI_celltype=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',EQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
eQTL_CI_celltype=merge(eQTL_CI_celltype, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
eQTL_CI_celltype[,type:='eQTL']
eQTL_CI_celltype[,state:=gsub('.*__(.*)','\\2',cellstate)]

# reQTL lineage
RUN_NAME=RUN_REQTL_LINEAGE
reQTL_CI_lineage=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',EQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
reQTL_CI_lineage=merge(reQTL_CI_lineage, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
reQTL_CI_lineage[,type:='reQTL']
reQTL_CI_lineage[,state:=gsub('.*__(.*)','\\2',cellstate)]

# reQTL celltype
RUN_NAME=RUN_REQTL_CELLTYPE
reQTL_CI_celltype=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',EQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
reQTL_CI_celltype=merge(eQTL_CI_celltype, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
reQTL_CI_celltype[,type:='reQTL']
reQTL_CI_celltype[,state:=gsub('.*__(.*)','\\2',cellstate)]
toc()

#############@ write eQTL genotype file for downstream analyses
tic('extracting eQTL genotypes')
all_eQTL_snps=unique(c(eQTL_Signif_both$snps,reQTL_Signif_both$snps))
toc()

cols_to_keep=c("gene","snps","celltype","state","gene_name","pvalue","beta","se","R2")
eQTL_stats=rbind(eQTL_Stats_celltype[,mget(cols_to_keep)][,type:='eQTL'],
                eQTL_Stats_lineage[,mget(cols_to_keep)][,type:='eQTL'],
                reQTL_Stats_celltype[,mget(cols_to_keep)][,type:='reQTL'],
                reQTL_Stats_lineage[,mget(cols_to_keep)][,type:='reQTL'])

#############@ define snp sets to consider for enrichment analyses
tic('defining eQTL snps sets')
snpSets_state=rbind(
  eQTL_Stats_lineage[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',state,sep='_'))],
  eQTL_Stats_celltype[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',state,sep='_'))],
  reQTL_Stats_lineage[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',state,sep='_'))],
  reQTL_Stats_celltype[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',state,sep='_'))])

snpSets_celltype=rbind(
      eQTL_Stats_celltype[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',celltype,sep='_'))],
      eQTL_Stats_lineage[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',celltype,sep='_'))],
      reQTL_Stats_celltype[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',celltype,sep='_'))],
      reQTL_Stats_lineage[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',celltype,sep='_'))])

snpSets_cellstate=rbind(
      eQTL_Stats_celltype[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',cellstate,sep='_'))],
      eQTL_Stats_lineage[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',cellstate,sep='_'))],
      reQTL_Stats_celltype[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',cellstate,sep='_'))],
      reQTL_Stats_lineage[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',cellstate,sep='_'))])

snpSets_reQTL=rbind(reQTL_celltype_compare[pvalue_COV<.01 & pvalue_IAV>0.01,.(snps=unique(snps),set='reQTL_COV_specific')],
                    reQTL_lineage_compare[pvalue_COV<.01 & pvalue_IAV>0.01,.(snps=unique(snps),set='reQTL_COV_specific')],
                    reQTL_celltype_compare[pvalue_COV>.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_IAV_specific')],
                    reQTL_lineage_compare[pvalue_COV>.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_IAV_specific')],
                    reQTL_celltype_compare[pvalue_COV<.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_shared')],
                    reQTL_lineage_compare[pvalue_COV<.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_shared')],
snpSets_reQTL=snpSets_reQTL[,.(set,snps)]

snpSets=unique(rbindlist(list(snpSets_state,snpSets_celltype,snpSets_cellstate,snpSets_reQTL)))
snpSets[,num:=cumsum(!duplicated(set))]
snpSets[,set:=make.names(set)]

allowed_celltypes=paste(c(lineage_5,celltype_22),collapse='|')
regex=sprintf('^(r?eQTL)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\2',set)]
snpSets[,state:=gsub(regex,'\\3',set)]
snpSets[,specificity:=gsub(regex,'\\4',set)]
toc()

fwrite(snpSets,file=sprintf('%s/All_eQTL_snpsSets.txt.gz',EQTL_DIR),sep='\t')
