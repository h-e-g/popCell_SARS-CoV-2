################################################################################
################################################################################
# File name: Fig1.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Figure 1
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="LIBRARY"
source(sprintf("./1a__quality_control__lib.R",LIB_DIR))
library(ggrastr)

# declare shortcuts
MISC_DIR="MISC"
source(sprintf("%s/shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("%s/misc_functions.R",MISC_DIR))

# declare useful functions and variables for plotting
source(sprintf("%s/set_colors.R",MISC_DIR))
source(sprintf("%s/misc_plots.R",MISC_DIR))

# read-in library ID
args <- commandArgs(TRUE)
LIB=args[1]

################################################################################
# Fig. 1a

# Scheme

################################################################################
# Fig. 1b

# load meta data created in 1b1__celltype_identification.R
meta_clean=fread(sprintf("%s/data/sce_clean_metadata.tsv",AGGR_QC_DIR))

fig1b_data=meta_clean[,.(Barcode_lib,umap1,umap2,COND,lineage,celltype)]

fig1b_plot=ggplot(fig1b_data,aes(x=umap1,y=umap2))+
  theme_plot() +
  rasterize(
    geom_point(
      aes(col=celltype,fill=celltype),
      size=0.01,alpha=0.5),
    dpi=480
  )+
  xlab('UMAP1')+ylab('UMAP2') +
  scale_color_manual(values=celltype_color,guide="none")+
  scale_fill_manual(
    values=celltype_color,
    breaks=celltype_order,
    labels=celltype_label,
    guide=guide_legend(override.aes=list(shape=21,size=2,color="black",alpha=1),nrow=23,byrow=F)
  )+
  theme(legend.position="right",legend.spacing.y=unit(-5,"pt"))

################################################################################
# Fig. 1c

fig1c_data <- fig1b_data
fig1c_data[,COND:=factor(COND,condition_order)]

# set level number for contour plots
nlevel=2

fig1c_plot=ggplot(fig1c_data,aes(umap1,umap2))+
  theme_plot() +
  scale_fill_distiller('Blues',direction = 1) +
  stat_density2d(data=fig1c_data[COND=="IAV",.(umap1,umap2,COND="IAV")],mapping=aes(fill = stat(nlevel)), geom = "polygon", n = 100, bins = 10, fill=condition_color["IAV"],size=0.1) +
  stat_density2d(data=fig1c_data[COND=="NS",.(umap1,umap2,COND="NS")],mapping=aes(fill = stat(nlevel)), geom = "polygon", n = 100, bins = 10, fill=condition_color["NS"],size=0.1) +
  stat_density2d(data=fig1c_data[COND=="COV",.(umap1,umap2,COND="COV")],mapping=aes(fill = stat(nlevel)), geom = "polygon", n = 100, bins = 10, fill=condition_color["COV"],size=0.1) +
  stat_density2d(data=fig1c_data[,.(umap1,umap2)], fill='lightgrey', alpha=0, geom = "polygon", n = 100, bins = 10, colour="#000000AA",size=0.1) +
  xlab('UMAP1')+ylab('UMAP2')+facet_grid(cols=vars(factor(COND,condition_order)))

################################################################################
# Fig. 1d

# load differential expression results computed in 1e1__differential_expression__condition.R
DE_COND=fread("1__transcriptome_processing/data/DE_condition_by_lineage.tsv.gz")

REF_DIR="1__transcriptome_processing/data"

# load hallmark interferon-stimulated pathway genes
hallmark_ifng_response_genes=fread(sprintf("%s/hallmark_ifng_response.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(DE_COND[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
hallmark_ifna_response_genes=fread(sprintf("%s/hallmark_ifna_response.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(DE_COND[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
hallmark_ifn_response_genes=unique(c(hallmark_ifng_response_genes,hallmark_ifna_response_genes))
names(hallmark_ifn_response_genes)=unique(DE_COND[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_ifn_response_genes]

# load hallmark inflammatory pathway genes
hallmark_inflammatory_response_genes=fread(sprintf("%s/hallmark_inflammatory_response.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(DE_COND[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
names(hallmark_inflammatory_response_genes)=unique(DE_COND[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_inflammatory_response_genes]

hallmark_inflammatory_exclusive=which(!hallmark_inflammatory_response_genes%in%hallmark_ifn_response_genes)%>%hallmark_inflammatory_response_genes[.]
hallmark_ifn_exclusive=which(!hallmark_ifn_response_genes%in%hallmark_inflammatory_response_genes)%>%hallmark_ifn_response_genes[.]

# annotate genes
features=unique(DE_COND[,.(Symbol,ID)])
features$pathway=case_when(
  features$ID%in%names(hallmark_inflammatory_response_genes)~"Inflammatory",
  features$ID%in%names(hallmark_ifn_response_genes)~"IFN-stimulated",
  T~"None"
)

fig1d_data=merge(x=DE_COND,y=features[,.(ID,pathway)],by="ID")
fig1d_data=fig1d_data[order(factor(pathway,c("None","IFN-stimulated","Inflammatory"))),]
fig1d_data[,celltype:=factor(celltype,lineage_order)]
pathway_colors=c("IFN-stimulated"="blue","Inflammatory"="orange","None"="gray")

fig1d_plot=ggplot(fig1d_data)+
  geom_hline(yintercept=0, linetype="dashed",size=0.1)+
  geom_vline(xintercept=0, linetype="dashed",size=0.1)+
  geom_abline(intercept=0,slope=1,size=0.1)+
  rasterise(
    geom_point(
      data=fig1d_data[(FDR_COV>=0.01|abs(logFC_COV)<=0.5)&(FDR_IAV>=0.01|abs(logFC_IAV)<=0.5 ),],
      mapping=aes(logFC_IAV,logFC_COV),alpha=0.2,color="lightgrey",fill="lightgrey",size=0.5),
    dpi=480
  )+
  rasterise(
    geom_point(
      data=fig1d_data[((abs(logFC_COV)>0.5&FDR_COV<0.01)|(abs(logFC_IAV)>0.5&FDR_IAV<0.01)),],
      mapping=aes(logFC_IAV,logFC_COV,color=pathway,fill=pathway),alpha=0.4,size=0.5),
    dpi=480
  )+
  xlab(expression("log"[2]*"FC(IAV/NS)"))+
  ylab(expression("log"[2]*"FC(COV/NS)"))+
  lims(x=c(-5,8.5),y=c(-5,8.5))+
  scale_fill_manual(name="Pathway",values=pathway_colors,breaks=c("Inflammatory","IFN-stimulated"))+
  facet_wrap(~celltype,ncol=5,labeller=labeller(celltype=c("MONO"="MONO","B"="B","T.CD4"="CD4+ T","T.CD8"="CD8+ T","NK"="NK")))+
  scale_color_manual(name="Pathway",values=pathway_colors,guide="none")+
  guides(fill=guide_legend(override.aes=list(pch=21,color="black",size=2,alpha=1),nrow=2))

################################################################################
# Fig. 1e

# load batch-adjusted counts computed in 1c2__pseudobulk_batch_correction.R
Expression=fread("1__transcriptome_processing/data/adjusted_pseudobulk_133libs__per_lineage_condition_IID.tsv.gz")

feature_data=Expression[,.(ID,Symbol)]%>%unique

IFNA1=feature_data[grepl("IFNA[0-9]",Symbol),setNames(ID,Symbol)]
type_I_IFNs=c(feature_data$Symbol[grepl("^IFN[ABKZ][0-9]",feature_data$Symbol)],"IFNE")%>%{setNames(feature_data$ID,feature_data$Symbol)[.]} %>% {.[!is.na(.)]}

# compute CPMs from logCPMs for IFNA genes
Expression_IFNA=Expression[ID%chin%IFNA1,.(ID,Symbol,celltype,IID,state,CPM=2^(logCPM)-1)]
# sum across all IFNA genes
Expression_IFNA=Expression_IFNA[,.(CPM=sum(CPM)),by=.(IID,celltype,condition=state)]
# count number of UMI in millions per cell type in each sample
celltype_condition_IID_counts=meta_clean[,.(N=.N,total_UMI_in_million=sum(sum)/1e6),by=.(IID,celltype,condition)]
Expression_IFNA=merge(Expression_IFNA,celltype_condition_IID_counts,by=c('IID',"celltype","condition"))

# compute number of counts in each sample and cell type
Expression_IFNA[,counts:=CPM*total_UMI_in_million]

# average across individuals
Expression_IFNA_avg=Expression_IFNA[,.(meanIFNAcount_perInd=mean(counts),meanIFNA_perCell=mean(CPM)),by=.(celltype,condition)]
Expression_IFNA_avg[,.(round(meanIFNAcount_perInd,1),round(meanIFNA_perCell,1)),by=.(celltype,condition)]

fig1e_data=Expression_IFNA_avg

fig1e_plot=ggplot(Expression_IFNA_avg[condition!="NS"])+
  theme_plot(rotate.x=90)+
  geom_segment(aes(x=0,xend=meanIFNAcount_perInd,y=celltype,yend=celltype,color=celltype))+
  scale_x_continuous(trans='sqrt',limits=c(-1,4500), breaks=c(10,250,1000,2500),labels=c('10','250','1000','2500'))+
  geom_point(aes(x=meanIFNAcount_perInd,y=celltype,fill=celltype,size=meanIFNA_perCell),pch=21,alpha=0.5)+
  geom_point(aes(x=meanIFNAcount_perInd,y=celltype,color=celltype,size=meanIFNA_perCell),pch=10,fill=NA)+
  scale_color_manual(aesthetics=c("color","fill"),guide="none",values=celltype_color) +
  xlab(expression("Mean number of IFN-"*alpha*" transcript per sample")) + ylab("Cell type") +
  #facet_grid(cols=vars(COND_ifn),labeller=labeller(COND_ifn=COND_ifn_labels))+guides(size=guide_legend(name=""))+
  facet_grid(cols=vars(condition))+
  guides(size=guide_legend(name=""))+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(angle=45,vjust=1))

################################################################################
# Fig. 1f

# load ISG scores computed in 1e3__scores.R
ISG_scores=fread(sprintf('1__transcriptome_processing/data/scores__by__IID_condition.tsv',DATA_DIR))
ISG_scores[,POP:=substr(IID,1,3)]

# load sample meta data
MORTALITY_FILE="../../0__barcode_processing/data/per_library_mortality.txt"
CellMortality=fread(MORTALITY_FILE)
setnames(CellMortality,c('MEAN','MORTALITY'),c('CellCount','Mortality'))
CellMortality=CellMortality[!duplicated(IID),.(IID,CellCount,Mortality)]
mod=CellMortality[,lm(Mortality~CellCount)]
CellMortality[is.na(Mortality),Mortality:=round(predict(mod,newdata=data.frame(CellCount=CellCount)))]

# adjust ISG score
ISG_scores=merge(ISG_scores,CellMortality,by=c("IID"))
ISG_scores[,hallmark_isg_score_adj:=adjust_on_x(hallmark_isg_score,Mortality,model.matrix(~1+POP)[,-1])]

# reshape data
fig1g_data=dcast(ISG_scores,IID+POP~condition,value.var="hallmark_isg_score_adj")

fig1g_plot=ggplot(fig1g_data,aes(x=IAV,y=COV,col=POP,fill=POP)) + theme_plot() + xlim(c(0,3)) + ylim(c(0,3)) + geom_abline(Intercept=0,slope=1,col='lightgrey',linetype=2) +
  geom_point(alpha=.5) + scale_color_manual(values=color_populations,guide="none") +
  scale_fill_manual(values=color_populations,breaks=c("AFB","EUB","ASH"),name="") +
  geom_smooth(method='lm',col='black',fill="grey",show.legend=FALSE)+
  ylab('ISG activity (COV)')+
  xlab("ISG activity (IAV)") +
  guides(fill=guide_legend(override.aes=list(pch=21,alpha=1,size=2,color="black"),ncol=1,byrow=TRUE))
