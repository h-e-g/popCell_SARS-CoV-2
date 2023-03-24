################################################################################
################################################################################
# File name: FigS6.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Figure S6
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="LIBRARY"
source(sprintf("./00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="MISC"
source(sprintf("%s/shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("%s/misc_functions.R",MISC_DIR))

# declare useful functions and variables for plotting
source(sprintf("%s/misc_plots.R",MISC_DIR))
source(sprintf("%s/load_eQTLs.R",MISC_DIR))
source(sprintf("%s/querySNPs.R",MISC_DIR))



################################################################################
# Fig. S6a : Number of eQTL/eGene per resolution

celltype_only_eQTLs=setdiff(eQTL_Signif_both[celltype%chin%celltype_order,unique(paste(snps,gene))],eQTL_Signif_both[celltype%chin%lineage_order,unique(paste(snps,gene))])
lineage_only_eQTLs=setdiff(eQTL_Signif_both[celltype%chin%lineage_order,unique(paste(snps,gene))],eQTL_Signif_both[celltype%chin%celltype_order,unique(paste(snps,gene))])
celltype_lineage_shared_eQTLs=intersect(eQTL_Signif_both[celltype%chin%lineage_order,unique(paste(snps,gene))],eQTL_Signif_both[celltype%chin%celltype_order,unique(paste(snps,gene))])

celltype_only_eGene=setdiff(eQTL_Signif_both[celltype%chin%celltype_order,unique(gene)],eQTL_Signif_both[celltype%chin%lineage_order,unique(gene)])
lineage_only_eGene=setdiff(eQTL_Signif_both[celltype%chin%lineage_order,unique(gene)],eQTL_Signif_both[celltype%chin%celltype_order,unique(gene)])
celltype_lineage_shared_eGene=intersect(eQTL_Signif_both[celltype%chin%lineage_order,unique(gene)],eQTL_Signif_both[celltype%chin%celltype_order,unique(gene)])


FigS6A_data=eQTL_Signif_both[!duplicated(paste(snps,gene)),.(snps,gene,celltype_only_eQTLs=paste(snps,gene)%in%celltype_only_eQTLs,
                                        lineage_only_eQTLs=paste(snps,gene)%in%lineage_only_eQTLs,
                                        celltype_lineage_shared_eQTLs=paste(snps,gene)%in%celltype_lineage_shared_eQTLs,
                                        celltype_only_eGene=gene%in%celltype_only_eGene,
                                        celltype_lineage_shared_eGene=gene%in%celltype_lineage_shared_eGene,
                                        lineage_only_eGene=gene%in%lineage_only_eGene)]
fwrite(FigS6A_data,file=sprintf("%s/FigS6A_data.tsv.gz",SOURCE_DATA),sep='\t')

library(ComplexUpset)
base_annotations=list(
    'Intersection size'=intersection_size(
        text=list(
            hjust=-0.05 ,
            vjust=0.5,
            angle=90
        )))

FigS6A <- upset(upsetData, intersect = c('celltype_only_eQTLs', 'lineage_only_eQTLs', 'celltype_lineage_shared_eQTLs', 'celltype_only_eGene', 'celltype_lineage_shared_eGene', 'lineage_only_eGene'), base_annotations=base_annotations)+theme_plot()


####------------------------------------------------------------------------------------####
FigS6A_data=upsetData
fwrite(FigS6A_data,file=sprintf('%s/Final/FigS6Adata_Nb_eQTL_by_method.tsv',FIGURE_DIR),sep='\t')
####------------------------------------------------------------------------------------####
pdf(sprintf('%s/Final/FigS6A_Nb_eQTL_by_method.pdf',FIGURE_DIR),height=5, width=5)
print(FigS6A)
dev.off()
####------------------------------------------------------------------------------------####


################################################################################
# Fig. S6b :  example of a pDC eQTL

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

eQTL_MIR155HG=get_eQTL('rs114273142','MIR155HG',resolution='celltype',metric='logCPM')
eQTL_MIR155HG[,state:=factor(state,c('NS','COV','IAV'))]
nIQR=3
eQTL_MIR155HG[,logCPM_IQR:=pmax(pmin(logCPM,median(logCPM)+nIQR*IQR(logCPM)),median(logCPM)-nIQR*IQR(logCPM)),by=.(state,celltype)]
eQTL_MIR155HG[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]

p <- ggplot(eQTL_MIR155HG[celltype%in%c('pDC',"MONO.CD14"),],aes(x=Number_of_ALT_alelle,y=logCPM_IQR,fill=celltype))
p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=TRUE)
p <- p + theme_plot() + scale_fill_manual(values=celltype_color[c('pDC','MONO.CD14')])# + scale_color_manual(values=color_populations)
p <- p + facet_grid(celltype~state,scales="free_y")
FigS6B=p+theme(legend.position='none',text = element_text(size = textSize))

####------------------------------------------------------------------------------------####
FigS6B_data=eQTL_MIR155HG
fwrite(FigS6B_data,file=sprintf('%s/FigS6Bdata_eQTL_MIR155HG_pDC.tsv',SOURCE_DATA_DIR))
####------------------------------------------------------------------------------------####
pdf(sprintf('%s/Final/FigS6B_eQTL_MIR155HG_pDC.pdf',FIGURE_DIR),width=7.2*.3,height=6.7*.4)
print(FigS6B)
dev.off()
####------------------------------------------------------------------------------------####
