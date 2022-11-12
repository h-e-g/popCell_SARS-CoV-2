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

source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",EVO_IMMUNO_POP_ZEUS))


##### How many eQTLs/reQTLs ?

reQTL_Signif_both[,length(unique(snps))]
reQTL_Signif_both[,length(unique(gene))]

##### what percentage of eQTLS detected in IAV are replicated upon SARS-COV-2 stimulation
IAV_reQTL=merge(reQTL_Signif_both[state=='IAV',],rbind(reQTL_celltype_compare,reQTL_lineage_compare),all.x=TRUE,by=c('snps','gene','gene_name','celltype'))
IAV_reQTL[,lineage:=ifelse(celltype %chin% celltype_22,celltype_2lineage[match(IAV_reQTL$celltype,celltype),lineage],celltype)]


IAV_reQTL[celltype %chin% lineage_5 ,any(pvalue_COV<.01),by=.(celltype,snps,gene_name)][,mean(V1),by=ifelse(celltype=='MONO','MONO','LYMPHOID')]
#     ifelse        V1
# 1: LYMPHOID 0.8198582
# 2:     MONO 0.6258993

COV_reQTL=merge(reQTL_Signif_both[state=='COV',],rbind(reQTL_celltype_compare,reQTL_lineage_compare),all.x=TRUE,by=c('snps','gene','gene_name','celltype'))
COV_reQTL[,lineage:=ifelse(celltype %chin% celltype_22,celltype_2lineage[match(COV_reQTL$celltype,celltype),lineage],celltype)]


COV_reQTL[celltype %chin% lineage_5 ,any(pvalue_IAV<.01),by=.(celltype,snps,gene_name)][,mean(V1),by=ifelse(celltype=='MONO','MONO','LYMPHOID')]
# ifelse        V1
# 1: LYMPHOID 0.9385027
# 2:     MONO 0.5652174

IAV_reQTL[celltype %chin% lineage_5 ,any(pvalue_COV<.01,na.rm=T),by=.(lineage,snps,gene_name)][,mean(V1,na.rm=T),by=lineage]
# celltype        V1
# 1:        B 0.8303571
# 2:    T.CD4 0.8155340
# 3:    T.CD8 0.9005848
# 4:     MONO 0.6258993
# 5:       NK 0.6991150

###### What fraction of reQTLs differ in strength between IAV and COV
Signif_reQTL_compare[,any(p_diff_COV_IAV<.01,na.rm=T),by=.(lineage,snps,gene_name)][,mean(V1,na.rm=T),by=ifelse(lineage=='MONO','MONO','LYMPHOID')]
# 1: LYMPHOID 0.1168521
# 2:     MONO 0.4983974

Signif_reQTL_compare[,any(p_diff_COV_IAV<.01,na.rm=T),by=.(lineage,snps,gene_name)][,mean(V1),by=lineage]
# lineage         V1
# 1:       B 0.13170732
# 2:   T.CD4 0.09839357
# 3:    MONO 0.49839744
# 4:   T.CD8 0.09482759
# 5:      NK 0.18357488

# focusing on lineage eQTL
Signif_reQTL_compare[celltype %chin% lineage_5 ,any(p_diff_COV_IAV<.01,na.rm=T),by=.(lineage,snps,gene_name)][,mean(V1,na.rm=T),by=ifelse(lineage=='MONO','MONO','LYMPHOID')]
# ifelse         V1
# 1: LYMPHOID 0.07701711
# 2:     MONO 0.48940678


Signif_reQTL_compare[celltype %chin% lineage_5 ,any(p_diff_COV_IAV<.01),by=.(lineage,snps,gene_name)][,mean(V1),by=lineage]
# lineage         V1
# 1:       B 0.06611570
# 2:   T.CD4 0.07651715
# 3:   T.CD8 0.04864865
# 4:    MONO 0.48940678
# 5:      NK 0.12781955

Signif_reQTL_compare[celltype %chin% lineage_5 & p_diff_COV_IAV<.01 & abs(beta_COV)>abs(beta_IAV),][order(p_diff_COV_IAV),]


FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/',EVO_IMMUNO_POP_ZEUS)

#########################################################

##### compare reQTL effect sizes
 # QTL_fineMapped=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_01pctFDR_mashr_and_FDR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
 # reQTL_Signif_lineage_compare=dcast(reQTL_Stats_lineage_mash,gene+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
 # reQTL_Signif_lineage_compare[,t_diff:=(beta_COV-beta_IAV)/sqrt(se_COV^2+se_IAV^2)]
 # reQTL_Signif_lineage_compare[,p_diff:=2*pnorm(abs(t_diff),lower=F)]
 # reQTL_Signif_lineage_compare[,fdr_diff:=p.adjust(p_diff,'fdr')]

# library(ggrastr)
# library(ggrepel)
# p <- ggplot(Signif_reQTL_compare[sample(1:.N),],aes(x=beta_COV,y=beta_IAV,col=lineage))
# p <- p + rasterize(geom_point(),dpi=400)+theme_yann()+xlim(c(-1.5,1.5))+ylim(c(-1.5,1.5))
# p <- p + scale_color_manual(values=lineage_color)
# p <- p + geom_label_repel(data=Signif_reQTL_compare[(p_diff_COV_IAV<.01 & abs(beta_COV-beta_IAV)>.7) | abs(beta_COV)>1 | abs(beta_IAV)>1,],aes(x=beta_COV,y=beta_IAV,col=lineage,label = gene_name))
# pdf(sprintf('%s/Fig3D_reQTL_Signif_lineage_compared_220409_withLabels.pdf',FIGURE_DIR))
# print(p)
# dev.off()

###########@ Figure 3D final -- compare reQTL effect sizes
library(ggrastr)
library(ggrepel)
  Fig3D_data=Signif_reQTL_compare[celltype%in%lineage_5,][sample(1:.N),]
  Fig3D_data=Fig3D_data[!duplicated(paste(snps,gene,lineage)),]
  p <- ggplot(Fig3D_data,aes(x=beta_COV,y=beta_IAV,col=lineage))+xlab(expression(beta[COV]))+ylab(expression(beta[IAV]))
  p <- p + rasterize(geom_point(alpha=.5),dpi=400)+xlim(c(-1.5,1.5))+ylim(c(-1.5,1.5))#+theme_yann()
  p <- p + scale_color_manual(values=lineage_color)
  p <- p + geom_hline(yintercept=0,col='lightgrey',linetype=2)+ geom_vline(xintercept=0,col='lightgrey',linetype=2)
#  p <- p + geom_label_repel(data=Fig3D_data[fdr_diff_COV_IAV<.01 & abs(beta_COV-beta_IAV)>1,],aes(x=beta_COV,y=beta_IAV,col=lineage,label = gene_name))
  pdf(sprintf('%s/Fig3D_reQTL_Signif_lineage_compared_220409_lineage_only_withLabels.pdf',FIGURE_DIR),width=5,height=5)
  print(p)
  dev.off()

Fig3D=p+theme(text = element_text(size = textSize))+guides(color=guide_legend(nrow=2,overide.aes=list(alpha=1), title.theme =element_blank()))
fwrite(Fig3D_data,file=sprintf('%s/Final/Fig3Ddata_reQTL_Signif_lineage_compared.tsv',FIGURE_DIR),sep='\t')
saveRDS(Fig3D,file=sprintf('%s/Final/Fig3D_reQTL_Signif_lineage_compared.RDS',FIGURE_DIR))
pdf(sprintf('%s/Final/Fig3D_reQTL_Signif_lineage_compared.pdf',FIGURE_DIR),width=7.2*.28,height=6.7*.4)
print(Fig3D)
dev.off()

######### doing per cell type is too noisy
##### lineage + celltypes
# library(ggrastr)
# p <- ggplot(Signif_reQTL_compare[sample(1:.N),],aes(x=beta_COV,y=beta_IAV,col=lineage))
# p <- p + rasterize(geom_point(),dpi=400)+theme_yann()+xlim(c(-1.5,1.5))+ylim(c(-1.5,1.5))
# p <- p + scale_color_manual(values=lineage_color)
# pdf(sprintf('%s/Fig3D_reQTL_Signif_lineage_compared_220409.pdf',FIGURE_DIR))
# print(p)
# dev.off()

##### lineage + celltypes
# p <- ggplot(Signif_reQTL_compare[celltype%in%celltype_22,][sample(1:.N),],aes(x=beta_COV,y=beta_IAV,col=celltype))
# p <- p + rasterize(geom_point(),dpi=400)+theme_yann()+xlim(c(-1.5,1.5))+ylim(c(-1.5,1.5))
# p <- p + scale_color_manual(values=celltype_color)
# p <- p + geom_label_repel(data=Signif_reQTL_compare[celltype%in%celltype_22,][(p_diff_COV_IAV<.01 & abs(beta_COV-beta_IAV)>.7) | abs(beta_COV)>1 | abs(beta_IAV)>1,],aes(x=beta_COV,y=beta_IAV,col=lineage,label = gene_name))
# pdf(sprintf('%s/Fig3D_reQTL_Signif_lineage_compared_220409_celltype_only_withLabels.pdf',FIGURE_DIR))
# print(p)
# dev.off()


# Signif_reQTL_compare[celltype%in%lineage_5,.(number_different=sum(p_diff_COV_IAV<0.01),keyby=.(celltype,sign_where_strongest=ifelse(abs(beta_COV)>abs(beta_IAV),sign(beta_COV),sign(beta_IAV)),stronger_in_COV=ifelse(abs(beta_COV)>abs(beta_IAV),TRUE,FALSE))]
Fig3E_data=Signif_reQTL_compare[celltype%in%lineage_5,][!duplicated(paste(snps,gene_name,lineage)),]
Fig3E_data=Fig3E_data[,.(number_different_praw1=sum(p_diff_COV_IAV<0.01),number_different_fdr1=sum(fdr_diff_COV_IAV<0.01)),keyby=.(lineage=celltype,stronger_in=ifelse(abs(beta_COV)>abs(beta_IAV),'COV','IAV'))]
Fig3E_data

sort(reQTL_Signif_lineage_compare[fdr_diff<0.01 & abs(beta_COV)>abs(beta_IAV),gene_name])
 # "AC119673.2" "CXCL16"     "D2HGDH"     "DSTYK"      "DTX2"       "HLA-DQA1"   "MMP1"       "MS4A7"      "POMZP3"     "PPP1R17"    "PRKCH"      "SERPINE2"   "VNN1"       "ZFHX3"

fwrite(Fig3E_data,file=sprintf('%s/Final/Fig3Edata_number_virus_specific_reQTLs_bylineage.tsv',FIGURE_DIR),sep='\t')

 p <- ggplot(Fig3E_data,aes(x=lineage,y=number_different_praw1,alpha=gsub('COV','Sars-CoV-2',stronger_in),fill=lineage))
 p <- p + geom_bar(position='stack',stat='Identity')+ theme_yann() + scale_alpha_manual(values=c('Sars-CoV-2'=1,IAV=.5))
 p <- p + scale_fill_manual(values=lineage_color) + ylab('Number of virus-dependent reQTL')
 p <- p + guides(fill= "none") + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1)) # theme(legend.title = element_text('reQTL stronger in response to'))
 pdf(sprintf('%s/Fig3E_number_virus_specific_reQTLs_bylineage_praw1.pdf',FIGURE_DIR),height=5,width=3)
 print(p)
 dev.off()

Fig3E=p+theme(text = element_text(size = textSize),axis.title.x=element_blank())+theme(legend.key.height = unit(.1,"cm"), legend.key.width = unit(.2,"cm"))+ guides(alpha = guide_legend(nrow = 2))
saveRDS(Fig3E,file=sprintf('%s/Final/Fig3E_number_virus_specific_reQTLs_bylineage_praw1.RDS',FIGURE_DIR))
pdf(sprintf('%s/Final/Fig3E_number_virus_specific_reQTLs_bylineage_praw1.pdf',FIGURE_DIR),width=7.2*.15,height=6.7*.4)
print(Fig3E)
dev.off()

 # Fig3E_data=Signif_reQTL_compare[celltype%in%celltype_22,][!duplicated(paste(snps,gene_name,lineage)),]
 # Fig3E_data=Fig3E_data[,.(number_different_praw1=sum(p_diff_COV_IAV<0.01),number_different_fdr1=sum(fdr_diff_COV_IAV<0.01)),keyby=.(lineage=celltype,stronger_in=ifelse(abs(beta_COV)>abs(beta_IAV),'COV','IAV'))]
 # Fig3E_data
 # p <- ggplot(Fig3E_data,aes(x=lineage,y=number_different_praw1,alpha=gsub('COV','Sars-CoV-2',stronger_in),fill=lineage))
 # p <- p + geom_bar(position='stack',stat='Identity')+ theme_yann() + scale_alpha_manual(values=c('Sars-CoV-2'=1,IAV=.5))
 # p <- p + scale_fill_manual(values=celltype_color) + ylab('Number of virus-dependent reQTL')
 # p <- p + guides(fill= "none") + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1)) # theme(legend.title = element_text('reQTL stronger in response to'))
 # pdf(sprintf('%s/Fig3E_number_virus_specific_reQTLs_bycelltype_praw1.pdf',FIGURE_DIR),height=5,width=3)
 # print(p)
 # dev.off()

 # p <- ggplot(Fig3E_data,aes(x=lineage,y=number_different_fdr1,alpha=gsub('COV','Sars-CoV-2',stronger_in),fill=lineage))
 # p <- p + geom_bar(position='stack',stat='Identity')+ theme_yann() + scale_alpha_manual(values=c('Sars-CoV-2'=1,IAV=.5))
 # p <- p + scale_fill_manual(values=color_cellTypes_6level[names(color_cellTypes_6level)!='DC']) + ylab('Number of virus-dependent reQTL')
 # p <- p + guides(fill= "none") + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1)) # theme(legend.title = element_text('reQTL stronger in response to'))
 # pdf(sprintf('%s/number_virus_specific_reQTLs_bylineage_fdr1.pdf',FIGURE_DIR),height=5,width=3)
 # print(p)
 # dev.off()




 ################################################################
 #########   Fig S5B  example of a pDC eQTL       ###############
 ################################################################


 getExpr=function(gene_name,resolution=c('lineage','celltype'),metric=c('logCPM','logFC')){
   resolution=match.arg(resolution)
   metric=match.arg(metric)
   gene_effect=fread(sprintf('gunzip -c %s/2_population_differences/BatchAdjusted_%s_125libs__per_%s_condition_annotated.tsv.gz | grep -e "%s\\|Symbol"',DATA_DIR,metric,resolution,gene_name))
   MinCell_perCOND=500
   keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
   gene_effect[IID%chin%keptIID,]
   }

 get_eQTL=function(rs_id,gene_name,...){
   myGene=getExpr(gene_name,...)
   mySNP=getSNP(rs_id,vector=FALSE)
   DT=merge(myGene,mySNP,by='IID')
 }
 myGene=getExpr('MIR155HG',resolution='celltype',metric='logCPM')
 mySNP=getSNP('rs114273142',vector=FALSE)


 # rs76415310 in HLA-DRB1 pDC COV
 # rs115005201 in HLA-DRB1 pDC IAV

 eQTL_MMP1=get_eQTL('rs534191','MMP1',resolution='lineage',metric='logCPM')
 eQTL_MMP1=eQTL_MMP1[Symbol=='MMP1',]
 eQTL_MMP1[,state:=factor(state,c('NS','COV','IAV'))]
 nIQR=3
 eQTL_MMP1[,logCPM_IQR:=pmax(pmin(logCPM,median(logCPM)+nIQR*IQR(logCPM)),median(logCPM)-nIQR*IQR(logCPM)),by=.(state,celltype)]
 eQTL_MMP1[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]
 irnt=function(x){qnorm(rank(x)/(length(x)+1),mean(x),sd(x))}
 eQTL_MMP1[,logCPM_irnt:=irnt(logCPM),by=.(state,celltype)]



 reQTL_MMP1=get_eQTL('rs534191','MMP1',resolution='lineage',metric='logFC')
 reQTL_MMP1=reQTL_MMP1[Symbol=='MMP1',]

 reQTL_MMP1[,state:=factor(state,c('COV','IAV'))]
 nIQR=3
 reQTL_MMP1[,logFC_IQR:=pmax(pmin(logFC,median(logFC)+nIQR*IQR(logFC)),median(logFC)-nIQR*IQR(logFC)),by=.(state,celltype)]
 reQTL_MMP1[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]
 irnt=function(x){qnorm(rank(x)/(length(x)+1),mean(x),sd(x))}
 reQTL_MMP1[,logFC_irnt:=irnt(logFC),by=.(state,celltype)]



 eQTL_CXCL16=get_eQTL('rs9901226','CXCL16',resolution='lineage',metric='logCPM')
 eQTL_CXCL16[,state:=factor(state,c('NS','COV','IAV'))]
 nIQR=3
 eQTL_CXCL16[,logCPM_IQR:=pmax(pmin(logCPM,median(logCPM)+nIQR*IQR(logCPM)),median(logCPM)-nIQR*IQR(logCPM)),by=.(state,celltype)]
 eQTL_CXCL16[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]
 irnt=function(x){qnorm(rank(x)/(length(x)+1),mean(x),sd(x))}
 eQTL_CXCL16[,logCPM_irnt:=irnt(logCPM),by=.(state,celltype)]

  reQTL_CXCL16=get_eQTL('rs9901226','CXCL16',resolution='lineage',metric='logFC')
  reQTL_CXCL16[,state:=factor(state,c('COV','IAV'))]
  nIQR=3
  reQTL_CXCL16[,logFC_IQR:=pmax(pmin(logFC,median(logFC)+nIQR*IQR(logFC)),median(logFC)-nIQR*IQR(logFC)),by=.(state,celltype)]
  reQTL_CXCL16[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]
  irnt=function(x){qnorm(rank(x)/(length(x)+1),mean(x),sd(x))}
  reQTL_CXCL16[,logFC_irnt:=irnt(logFC),by=.(state,celltype)]

 p <- ggplot(eQTL_MMP1[Symbol=='MMP1' & celltype%in%c('MONO'),],aes(x=Number_of_ALT_alelle,y=logCPM_IQR,fill=state))
 p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=TRUE)
 p <- p + theme_yann() + scale_fill_manual(values=color_conditions)# + scale_color_manual(values=color_populations)
 p <- p + facet_grid(state~celltype) +ylab('logCPM')

 ####------------------------------------------------------------------------------------####
 ####------------------------------------------------------------------------------------####
 Fig3F=p+theme(legend.position='none',text = element_text(size = textSize))
 Fig3F_data=eQTL_MMP1
 fwrite(Fig3F_data,file=sprintf('%s/Final/Fig3F_data_eQTL_MMP1_MONO.tsv',FIGURE_DIR))
 saveRDS(Fig3F,file=sprintf('%s/Final/Fig3F_data_eQTL_MMP1_MONO.RDS',FIGURE_DIR))
 pdf(sprintf('%s/Final/Fig3F_eQTL_MMP1_MONO.pdf',FIGURE_DIR),width=7.2*.2,height=6.7*.4)
 print(Fig3F)
 dev.off()



 p <- ggplot(eQTL_MMP1[Symbol=='MMP1' & celltype%in%c('MONO'),],aes(x=Number_of_ALT_alelle,y=logCPM_IQR,fill=state))
 p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=TRUE)
 p <- p + theme_yann() + scale_fill_manual(values=color_conditions)# + scale_color_manual(values=color_populations)
 p <- p + facet_grid(celltype~state) +ylab('logCPM')

 ####------------------------------------------------------------------------------------####
 ####------------------------------------------------------------------------------------####
 Fig3F=p+theme(legend.position='none',text = element_text(size = textSize))
 Fig3F_data=eQTL_MMP1
 fwrite(Fig3F_data,file=sprintf('%s/Final/Fig3F_data_eQTL_MMP1_MONO.tsv',FIGURE_DIR))
 saveRDS(Fig3F,file=sprintf('%s/Final/Fig3F_data_eQTL_MMP1_MONO.RDS',FIGURE_DIR))
 pdf(sprintf('%s/Final/Fig3F_eQTL_MMP1_MONO.pdf',FIGURE_DIR),width=7.2*.2,height=6.7*.4)
 print(Fig3F)
 dev.off()

 ################################################################################
 ######### read enrichments statistis, create a table and make a plot ###########
 ################################################################################

 ########################### TODO: run enrichment with new sets

 all_Enrich_results=list()

counter=0
EnrichFiles=dir(sprintf("%s/Enrichments_matchedDist/",FIGURE_DIR),pattern='resamp_COVID19_Enrichment_10000x_.*_220409.tsv.gz')
 for (file in EnrichFiles){
     counter=counter+1
     all_Enrich_results[[file]]=try(fread(sprintf("%s/Enrichments_matchedDist/%s",FIGURE_DIR,file)))
     if('try-error'%in%class(all_Enrich_results[[file]])){
       all_Enrich_results[[file]]=NULL
     }
     cat(counter,'')
   }

 all_Enrich_results=rbindlist(all_Enrich_results,idcol='ID')
 all_Enrich_results[,.N,by=ID]

########## fucntion to plot enrichments
 plot_enrichment=function(Enrich_table,DIR,pname,facet=FALSE,angle=0,ymax=15,width=7,...){
   color_COVID=c(reported="#F1BABA", hospitalized="#E68686",critical="#DD5555")
   Enrich_table[,GWAS_type_full:=factor(GWAS_type_full,names(color_COVID))]
   pdf(sprintf("%s/%s",DIR,pname),width=width,...)
   p <- ggplot(Enrich_table)+geom_point(aes(x=eQTL_group,y=FE,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5))
   p <- p + geom_errorbar(aes(x=eQTL_group,y=FE,ymin=FE_q025,ymax=FE_q975,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5)) + theme_yann(rotate.x=45)
   p <- p + geom_hline(yintercept=1,col='grey') + xlab('p-value threshold')+ ylab('Fold enrichment in GWAS loci') + coord_cartesian(ylim=c(0,ymax))
   p <- p + scale_color_manual(values=color_COVID)
   if(facet){
     p <- p + facet_grid(rows=vars(eQTLtype))
     p <- p + theme(axis.text.x=element_text(angle=angle,hjust=1,vjust=1))
   }
   print(p)
   dev.off()
   return(p)
 }

########################### TODO: create enrichment plot
#
# GW_signif_enrich=all_Enrich_results[threshold==0.0001,][grepl('eQTL_GWsignif',ID),]
# GW_signif_enrich[,ID:=gsub('resamp_COVID19_Enrichment_10000x_(.*)_220409.tsv.gz','\\1',ID)]
# GW_signif_enrich[,eQTLtype:=gsub('(.*)_GWsignif_(.*)','\\1',ID)]
# GW_signif_enrich[,eQTL_group:=gsub('(.*)_GWsignif_(.*)','\\2',ID)]
# GW_signif_enrich_lineage=GW_signif_enrich[eQTL_group%in%lineage_5,]
# GW_signif_enrich_celltype=GW_signif_enrich[eQTL_group%in%celltype_22,]
#
# plot_enrichment(GW_signif_enrich_lineage,FIGURE_DIR,"Fig3X_GW_signif_enrich_lineage.pdf",facet=TRUE)
# plot_enrichment(GW_signif_enrich_celltype,FIGURE_DIR,"Fig3X_GW_signif_enrich_celltype.pdf",facet=TRUE,width=10,angle=45)

standard_enrich=all_Enrich_results[threshold==0.0001 & !grepl('__',ID) & !grepl('eQTL\\.',ID) & !grepl('_GWsignif_',ID),]
standard_enrich[,ID:=gsub('resamp_COVID19_Enrichment_10000x_(.*)_220409.tsv.gz','\\1',ID)]
standard_enrich[,eQTLtype:=gsub('(.*)_(.*)','\\1',ID)]
standard_enrich[,eQTL_group:=gsub('(.*)_(.*)','\\2',ID)]
standard_enrich_lineage=standard_enrich[eQTL_group%in%lineage_5,]
# standard_enrich_celltype=standard_enrich[eQTL_group%in%celltype_22,]
standard_enrich_cond=standard_enrich[grepl('eQTL_(NS|IAV|COV|shared|same)',ID),]

# p <- plot_enrichment(standard_enrich_lineage,FIGURE_DIR,"Fig3X_enrich_lineage.pdf",facet=TRUE)
# p <- plot_enrichment(standard_enrich_celltype,FIGURE_DIR,"Fig3X_enrich_celltype.pdf",facet=TRUE,width=10,angle=45)

# standard_enrich_cond1=standard_enrich_cond[ID%in%c('eQTL_NS','eQTL_COV','eQTL_IAV'),]
# standard_enrich_cond1[,eQTL_group:=factor(ID,c('eQTL_NS','eQTL_COV','eQTL_IAV'))]
standard_enrich_cond2=standard_enrich_cond[ID%in%c('eQTL_NS','reQTL_COV','reQTL_IAV'),]
standard_enrich_cond2[,eQTL_group:=factor(ID,c('eQTL_NS','reQTL_COV','reQTL_IAV'))]
# standard_enrich_cond3=standard_enrich_cond[ID%in%c('eQTL_NS','reQTL_same_strength','reQTL_COV_stronger','reQTL_IAV_stronger'),]
# standard_enrich_cond3[,eQTL_group:=factor(ID,c('eQTL_NS','reQTL_same_strength','reQTL_COV_stronger','reQTL_IAV_stronger'))]
# standard_enrich_cond4=standard_enrich_cond[ID%in%c('eQTL_NS','reQTL_shared','reQTL_COV_specific','reQTL_IAV_specific'),]
# standard_enrich_cond4[,eQTL_group:=factor(ID,c('eQTL_NS','reQTL_shared','reQTL_COV_specific','reQTL_IAV_specific'))]

# p <- plot_enrichment(standard_enrich_cond1,FIGURE_DIR,"Fig3X_eQTL_enrich_condition1.pdf",facet=FALSE,angle=45)
# p <- plot_enrichment(standard_enrich_cond2,FIGURE_DIR,"Fig3X_eQTL_enrich_condition2.pdf",facet=FALSE,angle=45)
# p <- plot_enrichment(standard_enrich_cond3,FIGURE_DIR,"Fig3X_eQTL_enrich_condition3.pdf",facet=FALSE,angle=45)
# p <- plot_enrichment(standard_enrich_cond4,FIGURE_DIR,"Fig3X_eQTL_enrich_condition4.pdf",facet=FALSE,angle=45)
#
# Signif_reQTL_compare[celltype %chin% lineage_5 ,.(any(p_diff_COV_IAV<.01),any(p_diff_COV_IAV<.01 & abs(beta_COV)>abs(beta_IAV)),any(p_diff_COV_IAV<.01 & abs(beta_COV)<abs(beta_IAV))),by=.(lineage,snps,gene_name)][,.(.N,sum(V1),sum(V2),sum(V3)),by=lineage]
# Signif_reQTL_compare[,.(any(p_diff_COV_IAV<.01),any(p_diff_COV_IAV<.01 & abs(beta_COV)>abs(beta_IAV)),any(p_diff_COV_IAV<.01 & abs(beta_COV)<abs(beta_IAV))),by=.(lineage,snps,gene_name)][,.(.N,sum(V1,na.rm=T),sum(V2,na.rm=T),sum(V3,na.rm=T)),by=lineage]


#
# Fig3G <- plot_enrichment(standard_enrich_cond2,FIGURE_DIR,"Final/Fig3G_eQTL_enrich_by_condition.pdf",facet=FALSE,angle=45,ymax=10,width=3,height=5)
# Fig3G=Fig3G+guides(col=guide_legend(ncol=1,byrow = TRUE))+theme(text = element_text(size = textSize),axis.title.x=element_blank(),legend.key.height = unit(.01, 'mm'))
# saveRDS(Fig3G,file=sprintf('%s/Final/Fig3G_enrichment_COVID_loci.RDS',FIGURE_DIR))
# fwrite(standard_enrich_cond2,file=sprintf('%s/Final/Fig3Gdata_enrichment_COVID_loci.tsv',FIGURE_DIR),sep='\t')
# pdf(Fig3G,file=sprintf('%s/Final/Fig3G_enrichment_COVID_loci.pdf',FIGURE_DIR),width=7.2*.2,height=6.7*.4)
# print(Fig3G)
# dev.off()
# 

####### update Figure V9
textSize=8;
FIGURE_DIR_V9="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/paper_draft/V9/figureMaterial/Fig3"
Fig3G <- plot_enrichment(standard_enrich_cond2,FIGURE_DIR_V9,"/Fig3G_eQTL_enrich_by_condition.pdf",facet=FALSE,angle=45,ymax=10,width=3,height=5)
Fig3G=Fig3G+guides(col=guide_legend(ncol=1,byrow = TRUE))+theme(text = element_text(size = textSize),axis.title.x=element_blank(),legend.key.height = unit(.01, 'mm'))+ xlab('')
saveRDS(Fig3G,file=sprintf('%s/Fig3G_enrichment_COVID_loci.RDS',FIGURE_DIR_V9))
fwrite(standard_enrich_cond2,file=sprintf('%s/Fig3Gdata_enrichment_COVID_loci.tsv',FIGURE_DIR_V9),sep='\t')
pdf(Fig3G,file=sprintf('%s/Fig3G_enrichment_COVID_loci.pdf',FIGURE_DIR_V9),width=7.2*.2,height=6.7*.4)
print(Fig3G)
dev.off()

FigS6C <- plot_enrichment(standard_enrich_lineage,FIGURE_DIR_V9,"FigS6C_enrichment_COVID_loci_by_lineage.pdf",facet=TRUE,angle=45)
FigS6C=FigS6C+guides(col=guide_legend(ncol=1))+theme(text = element_text(size = textSize))
saveRDS(FigS6C,file=sprintf('%s/FigS6C_enrichment_COVID_loci_by_lineage.RDS',FIGURE_DIR_V9))
fwrite(standard_enrich_lineage,file=sprintf('%s/FigS6Cdata_enrichment_COVID_loci_by_lineage.tsv',FIGURE_DIR_V9),sep='\t')
pdf(FigS6C,file=sprintf('%s/FigS6C_enrichment_COVID_loci_by_lineage.pdf',FIGURE_DIR_V9),width=7.2*.4,height=6.7*.5)
print(FigS6C)
dev.off()


SNP_info=get_SNPinfo(annotate=TRUE)
Signif_reQTL_compare[gene_name=='CXCL16',]
SNP_info[ID%chin%Signif_reQTL_compare[gene_name=='CXCL16',snps],]
Signif_reQTL_compare[gene_name=='MMP1',]
SNP_info[ID%chin%Signif_reQTL_compare[gene_name=='MMP1',snps],]

RDS_files=dir(sprintf("%s/Final/",FIGURE_DIR),pattern='.RDS')


# Fig3A_legend=get_legend(Fig3A)
# Fig3B_legend=get_legend(Fig3B)
# Fig3C_legend=get_legend(Fig3C)
# Fig3D_legend=get_legend(Fig3D)
# Fig3E_legend=get_legend(Fig3E)
# Fig3F_legend=get_legend(Fig3F)
# Fig3G_legend=get_legend(Fig3G)

Fig3H1=readRDS(sprintf('%s/colocalization_COVID/locusZoom/IRF1_rs10066378/locusZoom_IRF1_rs10066378_covid_hospitalized.RDS',FIGURE_DIR))
Fig3H2=readRDS(sprintf('%s/colocalization_COVID/locusZoom/IRF1_rs10066378/locusZoom_IRF1_rs10066378_eQTL_NK.CD56dim__NS.RDS',FIGURE_DIR))

library(gridExtra)
pname=sprintf("%s/Final/Fig3.pdf",FIGURE_DIR)

  textSize=8;
  pdf(pname,width=7.2,height=6.7)
    grid.arrange(
      grobs=list(
        ggplotGrob(Fig3A+theme(legend.position="none",text = element_text(size = textSize),strip.background=element_rect(size=1))),
  #      Fig3A_legend,
        ggplotGrob(Fig3B+theme(axis.title.x=element_blank(),legend.position="top",text = element_text(size = textSize),legend.key.height = unit(.1,"cm"), legend.key.width = unit(.2,"cm"))),
  #     Fig3B_legend,Fig3B+theme()
        #grid.grabExpr(draw(Fig3C,heatmap_legend_side = "bottom")+theme(text = element_text(size = textSize))),
        grid.rect(gp=gpar(col="white")),
        #Fig3C_legend,
        ggplotGrob(Fig3D+theme(text = element_text(size = textSize),legend.key.height = unit(.05, 'mm'))),#+theme(legend.position="bottom",text = element_text(size = textSize)) + guides(shape = guide_legend(override.aes = list(shape = 16)))),
        # Fig3D_legend,
        #grid.rect(gp=gpar(col="white")),
        ggplotGrob(Fig3E+theme(text = element_text(size = textSize),axis.title.x=element_blank(),legend.key.height = unit(.01, 'mm'))),#+theme(legend.position="none",text = element_text(size = textSize))),
        # Fig3E_legend,
        #grid.rect(gp=gpar(col="white")),
        ggplotGrob(Fig3F+theme(strip.background=element_rect(size=.5),plot.margin = unit(c(0,0.1,0.1,0.1), "cm"))),
        ggplotGrob(Fig3G),#+theme(legend.position="none",text = element_text(size = textSize))),
        # Fig3F_legend,
        rbind(ggplotGrob(Fig3H1+theme(axis.title.x=element_blank())), ggplotGrob(Fig3H2)),
        # grid.rect(gp=gpar(col="white")),grid.rect(gp=gpar(col="white")),grid.rect(gp=gpar(col="white")),
        grid.rect(gp=gpar(col="white"))
      ),
      layout_matrix=rbind(
        c(1,2,2,2,3,3),
        c(1,2,2,2,3,3),
        c(6,6,7,9,10,10),
        #c(6,6,6,7,8,8)
        c(8,8,8,9,10,10),
        c(8,8,8,11,10,10)
      ), heights=c(3,.6,3,.3,1.1), widths=c(.22,.06,.15,.21,.16,.2)
    )
  dev.off()

pdf(sprintf("%s/Final/Fig3C.pdf",FIGURE_DIR))
print(Fig3C+theme(legend.position="none"))
dev.off()


FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3',EVO_IMMUNO_POP_ZEUS)
Fig3H1=readRDS(sprintf('%s/colocalization_COVID/locusZoom/IRF1_rs10066378/locusZoom_IRF1_rs10066378_covid_hospitalized.RDS',FIGURE_DIR))
Fig3H2=readRDS(sprintf('%s/colocalization_COVID/locusZoom/IRF1_rs10066378/locusZoom_IRF1_rs10066378_eQTL_NK.CD56dim__NS.RDS',FIGURE_DIR))
pname=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V11/newFigure_6/Fig6b_locusZoom_IRF1_rs10066378.pdf',EVO_IMMUNO_POP_ZEUS)

textSize=8;
pdf(pname,width=7.2*.38,height=6.7*4.4/8)
  grid.arrange(
    grobs=list(
        rbind(ggplotGrob(Fig3H1+theme_yann(lpos='none')+theme(axis.title.x=element_blank())), ggplotGrob(Fig3H2+theme_yann(lpos='none')))
      ),  layout_matrix=matrix(1,1,1),
      heights=c(1), widths=c(1)
    )
  dev.off()

  Fig3H1=readRDS(sprintf('%s/colocalization_COVID/locusZoom/IFNAR2_rs2834158/locusZoom_IFNAR2_rs2834158_covid_hospitalized.RDS',FIGURE_DIR))
  Fig3H2=readRDS(sprintf('%s/colocalization_COVID/locusZoom/IFNAR2_rs2834158/locusZoom_IFNAR2_rs2834158_eQTL_T.CD4__NS.RDS',FIGURE_DIR))
pname=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V11/newFigure_6/Fig6c_locusZoom_IFNAR2_rs2834158.pdf',EVO_IMMUNO_POP_ZEUS)
textSize=8;
pdf(pname,width=7.2*.38,height=6.7*4.4/8)
  grid.arrange(
    grobs=list(
        rbind(ggplotGrob(Fig3H1+theme_yann(lpos='none')+theme(axis.title.x=element_blank())), ggplotGrob(Fig3H2+theme_yann(lpos='none')))
      ),  layout_matrix=matrix(1,1,1),
      heights=c(1), widths=c(1)
    )
  dev.off()


########################### TODO: create enrichment plot

WOS_genes=fread(sprintf('%s/single_cell/resources/references/PATHWAYS/science.Randolph_2021_table_s10.tsv',EVO_IMMUNO_POP_ZEUS))
WOS_genes_melt=melt(WOS_genes,measure=patterns("lfsr$", "beta$"),value.name=c('lfsr', 'beta'))
WOS_genes_melt[,celltype:=c('MONO','T.CD8','T.CD4','NK','B')[variable]]
WOS_genes_melt[lfsr<.01 & beta>0,sort(unique(genes))]
severeCOVID_associated_genes=WOS_genes_melt[lfsr<.01 & beta>0,sort(unique(genes))]
mildCOVID_associated_genes=WOS_genes_melt[lfsr<.01 & beta<0,sort(unique(genes))]
background_genes=WOS_genes_melt[,sort(unique(genes))]
WOS_genes_melt

WOS_rEQTL=merge(WOS_genes_melt,Signif_reQTL_compare,by.x=c('genes','celltype'),by.y=c('gene_name','lineage'),suffix=c('.WOS','.reQTL'))
WOS_eQTL=merge(WOS_genes_melt,Signif_eQTL_compare,by.x=c('genes','celltype'),by.y=c('gene_name','lineage'),suffix=c('.WOS','.eQTL'))

WOS_rEQTL[lfsr<.01 & p_diff_COV_IAV<.01 & abs(beta_COV)>abs(beta_IAV) & !duplicated(genes),]
WOS_eQTL[lfsr<.01 & p_diff_COV_IAV<.01 & abs(beta_COV)>abs(beta_IAV) & !duplicated(genes),]


IAV_reQTL[celltype %chin% lineage_5 ,any(pvalue_COV<.01),by=.(celltype,snps,gene_name)][,mean(V1),by=ifelse(celltype=='MONO','MONO','LYMPHOID')]

IAV_reQTL[celltype %chin% lineage_5 ,any(pvalue_COV<.01),by=.(celltype,snps,gene_name)][,mean(V1),by=celltype]


IAV_reQTL[celltype %chin% lineage_5 ,any(pvalue_COV<.01),by=.(celltype,snps,gene_name)][,mean(V1),by=celltype=='MONO']

IAV_reQTL[celltype %chin% lineage_5 & lineage!='MONO',any(pvalue_COV<.01),by=.(celltype,snps,gene_name)][,mean(V1),by=celltype]

rbind(reQTL_celltype_compare,reQTL_lineage_compare)[celltype!='MONO.CD14.INFECTED',mean(2*pnorm(abs(t_diff_COV_IAV),lower=F)<0.01)]

rbind(reQTL_celltype_compare,reQTL_lineage_compare)[celltype!='MONO.CD14.INFECTED',mean(2*pnorm(abs(t_diff_COV_IAV),lower=F)<0.01,na.rm=T)]
# 0.022

# celltype_2lineage=data.table(celltype=c("B.M.K","B.M.L","B.N.K","B.N.L","Plasmablast","MONO.CD14","MONO.CD14.INFECTED", "MONO.CD16", "cDC","pDC","NK.CD56brt","NK.CD56dim","NK.M.LIKE","T.CD4.E","T.CD4.N","T.Reg","T.CD8.CM.EM","T.CD8.EMRA","T.CD8.N","ILC","MAIT","T.gd"),
# lineage=c("B","B","B","B","B","MONO","MONO", "MONO", "MONO","MONO","NK","NK","NK","T.CD4","T.CD4","T.CD4","T.CD8","T.CD8","T.CD8","T.CD8","T.CD8","T.CD8"))
#
# lineage_5=unique(celltype_2lineage$lineage)
# celltype_22=unique(celltype_2lineage$celltype)
#
#
# ###############################################################
# ########### extracte eQTLS (lineage or celltype) ##############
# ###############################################################
# RUN_NAME="lineage_condition___CellPropLineage_SVs_220409"
# # eQTL_Signif_lineage=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_1pctFDR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
# # eQTL_Signif_celltype=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/celltype_condition___CellPropLineage_SVs_220409/dist_100kb/independent_eQTLs_allcond_1pctFDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
#
# eQTL_Signif_both=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
#
# ###############################################################
# ########### extracte lineage level information on eQTLs #######
# ###############################################################
#
# eQTL_Stats_lineage_mash=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_mashr_and_FDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
# setnames(eQTL_Stats_lineage_mash,'gene_id','gene')
# eQTL_Stats_lineage_mash=merge(eQTL_Stats_lineage_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],
#                               eQTL_Signif_both[celltype%chin%lineage_5,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
# #
# #eQTL_fineMapped_compare=dcast(eQTL_fineMapped,gene_id+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
# eQTL_Stats_lineage_mash[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]
# Number_lineage_where_signif_padj=eQTL_Stats_lineage_mash[,.(Number_lineage_where_signif_padj=length(unique(celltype[p.adj<0.01])),
#                                                                 Number_lineage_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
#                                                                 Number_lineage_where_signif_lfsr1=length(unique(celltype[lfsr<0.01]))
#                                                                 ),by=.(gene,snps)]
# eQTL_Stats_lineage_mash=merge(eQTL_Stats_lineage_mash,Number_lineage_where_signif_padj,all.x=TRUE)
#
# ###############################################################
# ########### extract celltype level information on eQTLs #######
# ###############################################################
#
# eQTL_Stats_celltype_mash=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStatscelltype_mashr_and_FDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
# setnames(eQTL_Stats_celltype_mash,'gene_id','gene')
# eQTL_Stats_celltype_mash=merge(eQTL_Stats_celltype_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],
#                               eQTL_Signif_both[celltype%chin%celltype_22,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
#                               eQTL_Stats_celltype_mash[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]
# Number_celltype_where_signif_padj=eQTL_Stats_celltype_mash[,.(Number_celltype_where_signif_padj=length(unique(celltype[p.adj<0.01])),
#                                                                 Number_celltype_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
#                                                                 Number_celltype_where_signif_lfsr1=length(unique(celltype[lfsr<0.01]))
#                                                                 ),by=.(gene,snps)]
# eQTL_Stats_celltype_mash=merge(eQTL_Stats_celltype_mash,Number_celltype_where_signif_padj,all.x=TRUE)
#
# #########################################################
# ########### extract reQTLs (lineage only) ###############
# #########################################################
# RUN_NAME="lineage_condition_logFC___CellPropLineage_SVs_220409"
#
# # TODO: replace with Signif both
# reQTL_Signif_lineage=fread(sprintf('%s/%s/dist_%s/independent_eQTLs_allcond_1pctFDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
# reQTL_Signif_lineage[,state:=gsub('(.*)__(.*)','\\2',cellstate)]
# reQTL_Signif_lineage[,celltype:=gsub('(.*)__(.*)','\\1',cellstate)]
#
# ###############################################################
# ########### extract lineage level information on reQTLs #######
# ###############################################################
#
# reQTL_Stats_lineage_mash=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_allcond_1pctFDR_allStats_mashr_and_FDR.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
# setnames(reQTL_Stats_lineage_mash,'gene_id','gene')
# reQTL_Stats_lineage_mash=merge(reQTL_Stats_lineage_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],
#                               reQTL_Signif_lineage[celltype%chin%lineage_5,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],all.x=TRUE)
# #
# #eQTL_fineMapped_compare=dcast(eQTL_fineMapped,gene_id+gene_name+snps+celltype~state,value.var=list('beta','beta_mash','lfsr','pvalue','se'))
# reQTL_Stats_lineage_mash[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]
# Number_lineage_where_signif_padj=reQTL_Stats_lineage_mash[,.(Number_celltype_where_signif_padj=length(unique(celltype[p.adj<0.01])),
#                                                                 Number_celltype_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
#                                                                 Number_celltype_where_signif_lfsr1=length(unique(celltype[lfsr<0.01]))
#                                                                 ),by=.(gene,snps)]
# reQTL_Stats_lineage_mash=merge(reQTL_Stats_lineage_mash,Number_lineage_where_signif_padj,all.x=TRUE)
#
# ###############################################################
# ########### extract celltype level information on reQTLs ######
# ###############################################################
# # TODO
#
#
# ################################################
# ########### extracte eQTL 95% CIs ##############
# ################################################
#
# # eQTL lineage
# RUN_NAME="lineage_condition___CellPropLineage_SVs_220409"
# eQTL_CI_lineage=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
#
# # eQTL celltype
# RUN_NAME="celltype_condition___CellPropLineage_SVs_220409"
# eQTL_CI_celltype=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
#
# # reQTL lineage
# RUN_NAME="lineage_condition_logFC___CellPropLineage_SVs_220409"
# reQTL_CI_lineage=fread(sprintf('%s/%s/dist_%s/SusieR_95pct_CredibleSets_lbfover3_eQTLs_allCond.txt.gz',eQTL_DIR,RUN_NAME,CIS_DIST_TEXT))
#
# # reQTL celltype
# # TODO


#########################################
######## GO enrichment of reQTL #########
#########################################

RUN_NAME="lineage_condition_logFC___CellPropLineage_SVs_220409"
QTL_fineMapped=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_01pctFDR_mashr_and_FDR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))

resGO=GOSeq(unique(QTL_fineMapped$gene_id),feature_toUse)
resGO[1:10,1:10]
# Neutrophil degranulation/vesicle

resGO=GOSeq(unique(QTL_fineMapped_compare[lfsr_COV<.01 & lfsr_IAV<.01,gene_id]),feature_toUse)
resGO[1:10,1:10]
# Neutrophil degranulation/vesicle

resGO=GOSeq(unique(QTL_fineMapped_compare[t_diff>2,gene_id]),feature_toUse)
resGO[1:10,1:10]
# HLA

resGO=GOSeq(unique(QTL_fineMapped_compare[t_diff< -2,gene_id]),feature_toUse)
resGO[1:10,1:10]
# nothing

#########################################
######## Number of eQTL per gene ########
#########################################

RUN_NAME="lineage_condition___CellPropLineage_SVs_220409"
QTL_fineMapped=fread(sprintf('%s/%s/dist_%s/FineMapping/independent_eQTLs_01pctFDR_mashr_and_FDR.txt.gz',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
Nb_eQTL_perGene=QTL_fineMapped[,.(counts=length(unique(snps))),by=gene_id]

Nb_eQTL_perGene[,.N,keyby=counts]

p <- ggplot(Nb_eQTL_perGene[,.N,keyby=counts])+theme_yann()+geom_bar(aes(x=factor(counts,1:11),y=N),stat='identity')+ylim(c(0,3000))
p <- p + geom_text(aes(x=counts,y=N,label=N), vjust=-1) +

pdf(sprintf('%s/count_eQTL_per_gene_FineMapped_perCT_220409.pdf',FIGURE_DIR),height=4)
print(p)
dev.off()
