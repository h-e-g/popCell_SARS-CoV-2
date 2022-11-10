# ZZ_load_eQTLs
options(stringsAsFactors=FALSE, max.print=9999, width=200, datatable.fread.input.cmd.message=FALSE)
EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'

.libPaths(sprintf("%s/single_cell/resources/R_libs/4.1.0",EVO_IMMUNO_POP_ZEUS))

suppressMessages(library(data.table))
suppressMessages(library(tictoc))
suppressMessages(library(rtracklayer))

eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)

#################################################################
#################################################################
##### AIM: create eQTL related objects for eacy manipulation ####
#################################################################
#################################################################

# The following objects are created:

### 1 - RUN_EQTL_LINEAGE, RUN_REQTL_LINEAGE, RUN_EQTL_CELLTYPE, RUN_REQTL_CELLTYPE
# define the RUN folder where to retrived information on the various eQTL mapping levels

### 2 -  Feature_annot
# annotation of genes used in the data set (chr, position, etc...)

### 3 -  lineage_5, celltype_22, celltype_2lineage
# list of 5 lineages/22 celltypes used and their correspondance.

### 4 - eQTL_Signif_both / reQTL_Signif_both
# list of geneome wide significant eQTL/reQTLs at lineage & cell type level.
# celltype and lineage eQTLs have been merged so eQTLs with overlapping 95 CI are represented by the same SNP & a single SNP is provided for each independent gene x SNP association.
# one line per ,gene x SNP x celltype/lineage x condition

#### 5 - eQTL_Stats_lineage_mash, eQTL_Stats_celltype_mash, reQTL_Stats_lineage_mash, reQTL_Stats_celltype_mash
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
# list of eQTLs/reQTLs to consider for enrichment analyses (celltype/condition specfic)
# 132 sets in total
# use   snpSets[!grepl('__',set),]   to avoid considering all 89 combinations of celltype/lineage x condition (43 sets left)

#### 9 - eQTL_genotypes
# file containing eQTLs genotypes across 222 indiivduals for all 13839 eQTL/reQTL snps

#################################################################
#################################################################
#################################################################


CIS_DIST_TEXT='100kb'
RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_220409"
RUN_REQTL_LINEAGE="lineage_condition_logFC__logFC__CellPropLineage_SVs_220409"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_220409"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__logFC__CellPropLineage_SVs_220409"


tic('get gene_name for eQTL annotation')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

gtf <- rtracklayer::import(sprintf("%s/single_cell/resources/references/RNA/human_iav_sars/genes/genes.gtf",EVO_IMMUNO_POP_ZEUS))
Feature_annot=as.data.table(gtf)[type=='gene',.(gene_id,gene_name,seqnames, start, end, strand,gene_type)]
Feature_annot[is.na(gene_name),gene_name:=gene_id]
Feature_annot[,seqnames:=as.character(seqnames)]
# for IAV_M and IAV_NS  we have two lines with the same genename (2 transcripts)
# we only keep the longest transcript
Feature_annot=Feature_annot[!(gene_id=="IAV_M" & end==784) & !(gene_id=="IAV_NS" & end==719),]
# we match the remaining features to feature_toUse
Feature_annot=Feature_annot[match(feature_toUse,gene_id),]
toc()

##### compare reQTL effect sizes
celltype_2lineage=data.table(celltype=c("B.M.K","B.M.L","B.N.K","B.N.L","Plasmablast","MONO.CD14","MONO.CD14.INFECTED", "MONO.CD16", "cDC","pDC","NK.CD56brt","NK.CD56dim","NK.M.LIKE","T.CD4.E","T.CD4.N","T.Reg","T.CD8.CM.EM","T.CD8.EMRA","T.CD8.N","ILC","MAIT","T.gd"),
lineage=c("B","B","B","B","B","MONO","MONO", "MONO", "MONO","MONO","NK","NK","NK","T.CD4","T.CD4","T.CD4","T.CD8","T.CD8","T.CD8","T.CD8","T.CD8","T.CD8"))

lineage_5=unique(celltype_2lineage$lineage)
celltype_22=unique(celltype_2lineage$celltype)

###############################################################
########### extracte eQTLs (lineage or celltype) ##############
###############################################################
tic('extracting eQTL information')

# RUN_NAME=RUN_EQTL_LINEAGE
# eQTL_Signif_lineage=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_1pctFDR.txt.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT))
# RUN_NAME=RUN_EQTL_CELLTYPE
# eQTL_Signif_celltype=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_1pctFDR.txt.gz',eQTL_DIR,RUN_EQTL_CELLTYPE,CIS_DIST_TEXT))

RUN_NAME=RUN_EQTL_LINEAGE
eQTL_Signif_both=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
eQTL_Signif_both[,gene_name:=Feature_annot[match(gene,gene_id),gene_name]]
###############################################################
########### extracte lineage level information on eQTLs #######
###############################################################

eQTL_Stats_lineage_mash=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_mashr_and_FDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
setnames(eQTL_Stats_lineage_mash,'gene_id','gene')
eQTL_Stats_lineage_mash=merge(eQTL_Stats_lineage_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],
                              eQTL_Signif_both[celltype%chin%lineage_5,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
#
#eQTL_fineMapped_compare=dcast(eQTL_fineMapped,gene_id+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
eQTL_Stats_lineage_mash[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]
Number_lineage_where_signif_padj=eQTL_Stats_lineage_mash[,.(Number_lineage_where_signif_padj=length(unique(celltype[p.adj<0.01])),
                                                                Number_lineage_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
                                                                Number_lineage_where_signif_lfsr1=length(unique(celltype[lfsr<0.01]))
                                                                ),by=.(gene,snps)]
eQTL_Stats_lineage_mash=merge(eQTL_Stats_lineage_mash,Number_lineage_where_signif_padj,all.x=TRUE)
eQTL_Stats_lineage_mash[,cellstate:=paste(celltype,state,sep='__')]
rm(Number_lineage_where_signif_padj)


eQTL_lineage_compare=dcast(eQTL_Stats_lineage_mash,gene+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
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

eQTL_Stats_celltype_mash=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStatscelltype_mashr_and_FDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
setnames(eQTL_Stats_celltype_mash,'gene_id','gene')
eQTL_Stats_celltype_mash=merge(eQTL_Stats_celltype_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],
                              eQTL_Signif_both[celltype%chin%celltype_22,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
                              eQTL_Stats_celltype_mash[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]
Number_celltype_where_signif_padj=eQTL_Stats_celltype_mash[,.(Number_celltype_where_signif_padj=length(unique(celltype[p.adj<0.01])),
                                                                Number_celltype_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
                                                                Number_celltype_where_signif_lfsr1=length(unique(celltype[lfsr<0.01]))
                                                                ),by=.(gene,snps)]
eQTL_Stats_celltype_mash=merge(eQTL_Stats_celltype_mash,Number_celltype_where_signif_padj,all.x=TRUE)
eQTL_Stats_celltype_mash[,cellstate:=paste(celltype,state,sep='__')]
rm(Number_celltype_where_signif_padj)


eQTL_celltype_compare=dcast(eQTL_Stats_celltype_mash,gene+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
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

# reQTL_Signif_lineage=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_1pctFDR.txt.gz',eQTL_DIR,RUN_REQTL_LINEAGE,CIS_DIST_TEXT))
# reQTL_Signif_celltype=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_1pctFDR.txt.gz',eQTL_DIR,RUN_REQTL_CELLTYPE,CIS_DIST_TEXT))

RUN_NAME=RUN_REQTL_LINEAGE
reQTL_Signif_both=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
reQTL_Signif_both[,gene_name:=Feature_annot[match(gene,gene_id),gene_name]]

###############################################################
########### extract lineage level information on reQTLs #######
###############################################################

# TODO: extract stats for both cell type and lineage eQTLs
reQTL_Stats_lineage_mash=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_mashr_and_FDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
setnames(reQTL_Stats_lineage_mash,'gene_id','gene')
reQTL_Stats_lineage_mash=merge(reQTL_Stats_lineage_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],
                              reQTL_Signif_both[celltype%chin%lineage_5,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
#
#eQTL_fineMapped_compare=dcast(eQTL_fineMapped,gene_id+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
reQTL_Stats_lineage_mash[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]
Number_lineage_where_signif_padj=reQTL_Stats_lineage_mash[,.(Number_lineage_where_signif_padj=length(unique(celltype[p.adj<0.01])),
                                                                Number_lineage_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
                                                                Number_lineage_where_signif_lfsr1=length(unique(celltype[lfsr<0.01]))
                                                                ),by=.(gene,snps)]
reQTL_Stats_lineage_mash=merge(reQTL_Stats_lineage_mash,Number_lineage_where_signif_padj,all.x=TRUE)
reQTL_Stats_lineage_mash[,cellstate:=paste(celltype,state,sep='__')]

reQTL_lineage_compare=dcast(reQTL_Stats_lineage_mash,gene+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
reQTL_lineage_compare[,t_diff_COV_IAV:=(beta_COV-beta_IAV)/sqrt(se_COV^2+se_IAV^2)]
reQTL_lineage_compare[,p_diff_COV_IAV:=2*pnorm(abs(t_diff_COV_IAV),lower=F)]
reQTL_lineage_compare[,fdr_diff_COV_IAV:=p.adjust(p_diff_COV_IAV,'fdr')]


###############################################################
########### extract celltype level information on reQTLs ######
###############################################################
reQTL_Stats_celltype_mash=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStatscelltype_mashr_and_FDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
setnames(reQTL_Stats_celltype_mash,'gene_id','gene')
reQTL_Stats_celltype_mash=merge(reQTL_Stats_celltype_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],
                              reQTL_Signif_both[celltype%chin%celltype_22,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
#
#eQTL_fineMapped_compare=dcast(eQTL_fineMapped,gene_id+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
reQTL_Stats_celltype_mash[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]
Number_celltype_where_signif_padj=reQTL_Stats_celltype_mash[,.(Number_celltype_where_signif_padj=length(unique(celltype[p.adj<0.01])),
                                                                Number_celltype_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
                                                                Number_celltype_where_signif_lfsr1=length(unique(celltype[lfsr<0.01]))
                                                                ),by=.(gene,snps)]
reQTL_Stats_celltype_mash=merge(reQTL_Stats_celltype_mash,Number_celltype_where_signif_padj,all.x=TRUE)
reQTL_Stats_celltype_mash[,cellstate:=paste(celltype,state,sep='__')]

reQTL_celltype_compare=dcast(reQTL_Stats_celltype_mash,gene+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
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

toc()

################################################
########### extract eQTL 95% CIs ##############
################################################
tic('extracting 95% credible intervals')
# eQTL lineage
RUN_NAME=RUN_EQTL_LINEAGE
eQTL_CI_lineage=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
eQTL_CI_lineage=merge(eQTL_CI_lineage, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
eQTL_CI_lineage[,type:='eQTL']
eQTL_CI_lineage[,state:=gsub('.*__(.*)','\\2',cellstate)]

# eQTL celltype
RUN_NAME=RUN_EQTL_CELLTYPE
eQTL_CI_celltype=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
eQTL_CI_celltype=merge(eQTL_CI_celltype, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
eQTL_CI_celltype[,type:='eQTL']
eQTL_CI_celltype[,state:=gsub('.*__(.*)','\\2',cellstate)]

# reQTL lineage
RUN_NAME=RUN_REQTL_LINEAGE
reQTL_CI_lineage=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
reQTL_CI_lineage=merge(reQTL_CI_lineage, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
reQTL_CI_lineage[,type:='reQTL']
reQTL_CI_lineage[,state:=gsub('.*__(.*)','\\2',cellstate)]

# reQTL celltype
RUN_NAME=RUN_REQTL_CELLTYPE
reQTL_CI_celltype=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
reQTL_CI_celltype=merge(eQTL_CI_celltype, Feature_annot[,.(gene=gene_id, Symbol=gene_name)],by='gene');
reQTL_CI_celltype[,type:='reQTL']
reQTL_CI_celltype[,state:=gsub('.*__(.*)','\\2',cellstate)]
toc()

#############@ define snp sets to consider for enrichment analyses
tic('defining eQTL snps sets')
snpSets=rbind(eQTL_Stats_celltype_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',celltype,sep='_'))],
      eQTL_Stats_celltype_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',cellstate,sep='_'))],
      eQTL_Stats_lineage_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',celltype,sep='_'))],
      eQTL_Stats_lineage_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',cellstate,sep='_'))],
      eQTL_Stats_lineage_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',state,sep='_'))],
      eQTL_Stats_celltype_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('eQTL',state,sep='_'))],
      reQTL_Stats_celltype_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',celltype,sep='_'))],
      reQTL_Stats_celltype_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',cellstate,sep='_'))],
      reQTL_Stats_lineage_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',celltype,sep='_'))],
      reQTL_Stats_lineage_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',cellstate,sep='_'))],
      reQTL_Stats_lineage_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',state,sep='_'))],
      reQTL_Stats_celltype_mash[pvalue<.01,.(snps=unique(snps)),by=.(set=paste('reQTL',state,sep='_'))])

snpSets_reQTL=rbind(reQTL_celltype_compare[pvalue_COV<.01 & pvalue_IAV>0.01,.(snps=unique(snps),set='reQTL_COV_specific')],
                    reQTL_lineage_compare[pvalue_COV<.01 & pvalue_IAV>0.01,.(snps=unique(snps),set='reQTL_COV_specific')],
                    reQTL_celltype_compare[pvalue_COV>.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_IAV_specific')],
                    reQTL_lineage_compare[pvalue_COV>.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_IAV_specific')],
                    reQTL_celltype_compare[pvalue_COV<.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_shared')],
                    reQTL_lineage_compare[pvalue_COV<.01 & pvalue_IAV<0.01,.(snps=unique(snps),set='reQTL_shared')],
                    reQTL_celltype_compare[pvalue_COV<.01 & abs(t_diff_COV_IAV)>2 & abs(beta_COV)>abs(beta_IAV),.(snps=unique(snps),set='reQTL_COV_stronger')],
                    reQTL_lineage_compare[pvalue_COV<.01 & abs(t_diff_COV_IAV)>2 & abs(beta_COV)>abs(beta_IAV),.(snps=unique(snps),set='reQTL_COV_stronger')],
                    reQTL_celltype_compare[pvalue_IAV<.01 & abs(t_diff_COV_IAV)>2 & abs(beta_IAV)>abs(beta_COV),.(snps=unique(snps),set='reQTL_IAV_stronger')],
                    reQTL_lineage_compare[pvalue_IAV<.01 & abs(t_diff_COV_IAV)>2 & abs(beta_IAV)>abs(beta_COV),.(snps=unique(snps),set='reQTL_IAV_stronger')],
                    reQTL_celltype_compare[(pvalue_IAV<.01 | pvalue_COV<.01) & abs(t_diff_COV_IAV)<=2 ,.(snps=unique(snps),set='reQTL_same_strength')],
                    reQTL_lineage_compare[(pvalue_IAV<.01 | pvalue_COV<.01) & abs(t_diff_COV_IAV)<=2 ,.(snps=unique(snps),set='reQTL_same_strength')])
snpSets_reQTL=snpSets_reQTL[,.(set,snps)]

snpSets_Signif=rbind(eQTL_Signif_both[,.(snps=unique(snps)),by=.(set=paste('eQTL_GWsignif',celltype,sep='_'))],
#                      eQTL_Signif_both[,.(snps=unique(snps)),by=.(set=paste('eQTL_GWsignif',cellstate,sep='_'))],
                      eQTL_Signif_both[,.(snps=unique(snps)),by=.(set=paste('eQTL_GWsignif',state,sep='_'))],
                      reQTL_Signif_both[,.(snps=unique(snps)),by=.(set=paste('reQTL_GWsignif',celltype,sep='_'))],
#                      reQTL_Signif_both[,.(snps=unique(snps)),by=.(set=paste('reQTL_GWsignif',cellstate,sep='_'))],
                      reQTL_Signif_both[,.(snps=unique(snps)),by=.(set=paste('reQTL_GWsignif',state,sep='_'))])






snpSets=unique(rbindlist(list(snpSets,snpSets_reQTL,snpSets_Signif)))
snpSets[,num:=cumsum(!duplicated(set))]
snpSets[,set:=make.names(set)]
# reassign number of eQTLs that are detected only at lineage level
snpSets[,num:=snpSets[!duplicated(set),.(num,newset=set)][match(set,newset),num]]
# fwrite(snpSets,file=sprintf('%s/%s/dist_%s/All_eQTL_snpsSets.txt.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT),sep='\t')
# snpSets=fread(sprintf('%s/%s/dist_%s/All_eQTL_snpsSets.txt.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT))
allowed_celltypes=paste(c(lineage_5,celltype_22),collapse='|')
regex=sprintf('^(r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\3',set)]
snpSets[,state:=gsub(regex,'\\4',set)]
snpSets[,specificity:=gsub(regex,'\\5',set)]
toc()

###### create SNP sets with fixed number of eQTL
# snpSets_top1kSNPs=eQTL_Stats_celltype_mash[order(pvalue),head(snps[!duplicated(snps)],1000),by=.(set=paste('eQTL.Top1000',celltype,sep='_'))]
# snpSets_top500SNPs=eQTL_Stats_celltype_mash[order(pvalue),head(snps[!duplicated(snps)],500),by=.(set=paste('eQTL.Top500',celltype,sep='_'))]
# eQTL_Stats_celltype_mash[,beta_lower:=sign(beta)*pmax(abs(beta)-2*se,0)]
# snpSets_top1kSNPs_beta=eQTL_Stats_celltype_mash[order(-abs(beta_lower)),head(snps[!duplicated(snps)],500),by=.(set=paste('eQTL.Top1000.beta',celltype,sep='_'))]
# snpSets_top500SNPs_beta=eQTL_Stats_celltype_mash[order(-abs(beta_lower)),head(snps[!duplicated(snps)],500),by=.(set=paste('eQTL.Top500.beta',celltype,sep='_'))]
#
# snpSets_topSNPs=rbind(snpSets_top1kSNPs,snpSets_top500SNPs,snpSets_top1kSNPs_beta,snpSets_top500SNPs_beta)
# snpSets_topSNPs[,num:=256+cumsum(!duplicated(set))]
# setnames(snpSets_topSNPs,'V1','snps')
# snpSets_topSNPs=rbind(snpSets[,.(set,snps,num)],snpSets_topSNPs[,.(set,snps,num)])
# fwrite(snpSets_topSNPs,sprintf('%s/users/Javier/data/snp_sets/snpSets_topSNPs_nov2022.txt.gz',EIP),sep='\t')

#############@ write eQTL genotype fiole for downstream anlayses
tic('extracting eQTL genotypes')
all_eQTL_snps=unique(c(reQTL_Signif_both$snps,eQTL_Signif_both$snps))

# source(sprintf("%s/single_cell/resources/template_scripts/querySNPs.R",EVO_IMMUNO_POP_ZEUS))
#
# SNP_info=getMap(annotate=TRUE)
#
# tic('loading eQTL genotypes')
# eQTL_genotypes=getSNP(all_eQTL_snps,vector=FALSE, Map=SNP_info)
# eQTL_genotypes[,IID:=as.character(IID)]
# toc()
# MinCell_perCOND=500
# keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
# eQTL_genotypes=eQTL_genotypes[ID%chin%all_eQTL_snps & IID%chin%keptIID,]
# fwrite(eQTL_genotypes,file=sprintf('%s/%s/dist_%s/All_eQTL_and_reQTL_genotypes.tsv.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT),sep='\t')
eQTL_genotypes=fread(sprintf('%s/All_eQTL_and_reQTL_genotypes.tsv.gz',eQTL_DIR))
toc()

cols_to_keep=c("gene","snps","celltype","state","gene_name","pvalue","beta","se","R2")
eQTL_stats=rbind(eQTL_Stats_celltype_mash[,mget(cols_to_keep)][,type:='eQTL'],
                eQTL_Stats_lineage_mash[,mget(cols_to_keep)][,type:='eQTL'],
                reQTL_Stats_celltype_mash[,mget(cols_to_keep)][,type:='reQTL'],
                reQTL_Stats_lineage_mash[,mget(cols_to_keep)][,type:='reQTL'])
eQTL_stats[,beta_low:=sign(beta)*pmax(0,abs(beta)-2*se)]
eQTL_stats[,beta_high:=sign(beta)*pmax(0,abs(beta)+2*se)]
