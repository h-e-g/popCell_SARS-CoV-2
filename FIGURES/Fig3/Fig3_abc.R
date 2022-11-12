options(stringsAsFactors=FALSE, max.print=9999, width=200, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(mashr))
suppressMessages(library(tictoc))
suppressMessages(library(readr))
suppressMessages(library(rtracklayer))
suppressMessages(require(susieR))
suppressMessages(library(ComplexHeatmap))


EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3',EVO_IMMUNO_POP_ZEUS)
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)
eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
SCRATCH="/pasteur/appa/scratch/mrotival/"

theme_set(theme_bw())
theme_update(
  text=element_text(family="serif",size=12),
  panel.grid=element_blank(),legend.position="bottom",
  strip.background=element_rect(fill="#012158"),strip.text=element_text(color="white")
)
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EVO_IMMUNO_POP_ZEUS))

source(sprintf("%s/single_cell/resources/template_scripts/querySNPs.R",EVO_IMMUNO_POP_ZEUS))
source(sprintf("%s/single_cell/resources/template_scripts/GOSeq.R",EVO_IMMUNO_POP_ZEUS))

source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",EVO_IMMUNO_POP_ZEUS))

# eQTL_Stats_lineage_mash Stats_withmash=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/FineMapping/independent_eQTLs_allcond_1pctFDR_allStats_mashr_and_FDR.txt.gz',EVO_IMMUNO_POP_ZEUS))

# eQTL_lineage_signif=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/independent_eQTLs_allcond_1pctFDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# eQTL_lineage_signif[,gene_name:=eQTL_lineage[match(eQTL_lineage_signif$gene,eQTL_lineage$gene_id),gene_name]]
#
#
# eQTL_Signif_lineage=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/independent_eQTLs_allcond_1pctFDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# eQTL_Signif_lineage[,length(unique(snps))]
# eQTL_Signif_lineage[,length(unique(gene))]
#
# eQTL_Signif_celltype=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/celltype_condition___CellPropLineage_SVs_220409/dist_100kb/independent_eQTLs_allcond_1pctFDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# eQTL_Signif_celltype[,length(unique(snps))]
# eQTL_Signif_celltype[,length(unique(gene))]
# celltype_22=unique(eQTL_Signif_celltype$celltype)
# lineage_5=unique(eQTL_Signif_lineage$celltype)
#
#
# eQTL_celltype=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/celltype_condition___CellPropLineage_SVs_220409/dist_100kb/FineMapping/independent_eQTLs_allcond_1pctFDR_allStats_mashr_and_FDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# #eQTL_celltype=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/celltype_condition___CellPropLineage_SVs_220409/dist_100kb/FineMapping/independent_eQTLs_allcond_1pctFDR_allStats_allCHR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# eQTL_Stats_lineage=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_allCHR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# eQTL_Stats_celltype=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStatscelltype_allCHR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# eQTL_Stats_celltype_mash=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStatscelltype_mashr_and_FDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# setnames(eQTL_Stats_celltype_mash,'gene_id','gene')
# # eQTL_Stats_lineage_mash=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/FineMapping/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStats_mashr_and_FDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# # setnames(eQTL_Stats_lineage_mash,'gene_id','gene')
# #


# eQTL_Signif_both=fread(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition___CellPropLineage_SVs_220409/dist_100kb/independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR.txt.gz',EVO_IMMUNO_POP_ZEUS))
# eQTL_Signif_both[,length(unique(snps))]
# eQTL_Signif_both[,length(unique(gene))]
# eQTL_Signif_both[,gene_name:=eQTL_Stats_celltype_mash$gene_name[match(eQTL_Signif_both$gene,eQTL_Stats_celltype_mash$'gene')]]
# eQTL_Signif_both[celltype%in%lineage_5,.N,by=.(gene,celltype)][,.N,keyby=.(celltype,pmin(N,5)]

eQTL_Signif_both[celltype%chin%celltype_22,length(unique(snps))]
eQTL_Signif_both[celltype%chin%celltype_22,length(unique(gene))]
eQTL_Signif_both[celltype%chin%lineage_5,length(unique(snps))]
eQTL_Signif_both[celltype%chin%lineage_5,length(unique(gene))]

######
celltype_only_eQTLs=setdiff(eQTL_Signif_both[celltype%chin%celltype_22,unique(paste(snps,gene))],eQTL_Signif_both[celltype%chin%lineage_5,unique(paste(snps,gene))])
lineage_only_eQTLs=setdiff(eQTL_Signif_both[celltype%chin%lineage_5,unique(paste(snps,gene))],eQTL_Signif_both[celltype%chin%celltype_22,unique(paste(snps,gene))])
celltype_lineage_shared_eQTLs=intersect(eQTL_Signif_both[celltype%chin%lineage_5,unique(paste(snps,gene))],eQTL_Signif_both[celltype%chin%celltype_22,unique(paste(snps,gene))])
length(celltype_only_eQTLs) # 3681
length(lineage_only_eQTLs) # 2991
length(celltype_lineage_shared_eQTLs) # 6354

celltype_only_eGene=setdiff(eQTL_Signif_both[celltype%chin%celltype_22,unique(gene)],eQTL_Signif_both[celltype%chin%lineage_5,unique(gene)])
lineage_only_eGene=setdiff(eQTL_Signif_both[celltype%chin%lineage_5,unique(gene)],eQTL_Signif_both[celltype%chin%celltype_22,unique(gene)])
celltype_lineage_shared_eGene=intersect(eQTL_Signif_both[celltype%chin%lineage_5,unique(gene)],eQTL_Signif_both[celltype%chin%celltype_22,unique(gene)])
length(celltype_only_eGene) # 777
length(lineage_only_eGene) # 1057
length(celltype_lineage_shared_eGene) # 4141

celltype_only_snps=setdiff(eQTL_Signif_both[celltype%chin%celltype_22,unique(snps)],eQTL_Signif_both[celltype%chin%lineage_5,unique(snps)])
lineage_only_snps=setdiff(eQTL_Signif_both[celltype%chin%lineage_5,unique(snps)],eQTL_Signif_both[celltype%chin%celltype_22,unique(snps)])
celltype_lineage_shared_snps=intersect(eQTL_Signif_both[celltype%chin%lineage_5,unique(snps)],eQTL_Signif_both[celltype%chin%celltype_22,unique(snps)])
length(celltype_only_snps) # 3603
length(lineage_only_snps) # 2919
length(celltype_lineage_shared_snps) # 6231


###################################################################################
#########   Fig 3A_00 : Number of eQTL_per Gene in each lineage        ############
###################################################################################

eQTL_perGene_lineage=eQTL_Signif_both[celltype%chin%lineage_5,.(N_eQTL_perGene=length(unique(snps))),by=.(celltype,gene)][,.N,keyby=.(celltype,N_eQTL_perGene)]
eQTL_perGene_lineage[,N_eQTL_perGene_maxed:=factor(ifelse(N_eQTL_perGene<4,as.character(N_eQTL_perGene),'4+'),rev(c('1','2','3','4+')))]
fwrite(eQTL_perGene_lineage,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/Fig3A_eQTL_perGene_lineage_data.tsv",EVO_IMMUNO_POP_ZEUS),sep='\t')

eQTL_perGene_lineage=eQTL_perGene_lineage[,.(Ngene=sum(N)),keyby=.(celltype,N_eQTL_perGene_maxed)]
p <- ggplot(eQTL_perGene_lineage,aes(x=N_eQTL_perGene_maxed, y=Ngene, fill=celltype))+geom_bar(stat='Identity')+facet_grid(row=vars(celltype))+ theme_yann()
p <- p + scale_fill_manual(values=lineage_color) + xlab('Number of eQTL per gene') + ylab('Number of genes') + guides(fill= "none") + coord_flip()
# pdf(sprintf('%s/Fig3A_00_eQTL_perGene_lineage.pdf',FIGURE_DIR),height=6, width=2.5)
# print(p)
# dev.off()

####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####
Fig3A=p
Fig3A_data=eQTL_perGene_lineage
fwrite(Fig3A_data,file=sprintf('%s/Final/Fig3data_eQTL_perGene_lineage.tsv',FIGURE_DIR))
saveRDS(Fig3A,file=sprintf('%s/Final/Fig3A_eQTL_perGene_lineage.RDS',FIGURE_DIR))
pdf(sprintf('%s/Final/Fig3A_eQTL_perGene_lineage.pdf',FIGURE_DIR),width=7.2*.22,height=6.7*.4)
print(Fig3A)
dev.off()
####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####


eQTL_perGene_lineage=eQTL_Signif_both[celltype%chin%lineage_5,.(N_eQTL_perGene=length(unique(snps))),by=.(celltype,gene)][,.N,keyby=.(celltype,N_eQTL_perGene)]

###################################################################################
#########   Fig 3A_0 : overlap between lineage and celltype eQTLs     #############
###################################################################################

eQTL_Signif_both[,eQTL_type:=case_when(paste(snps,gene)%in%celltype_only_eQTLs~'celltype_only',
                                      paste(snps,gene)%in%lineage_only_eQTLs~'lineage_only',
                                      paste(snps,gene)%in%celltype_lineage_shared_eQTLs~'shared')]

eQTL_Signif_both[,eGene_type:=case_when(gene%in%celltype_only_eGene~'celltype_only',
                                        gene%in%lineage_only_eGene~'lineage_only',
                                        gene%in%celltype_lineage_shared_eGene~'shared')]

eQTL_Signif_both[!duplicated(paste(snps,gene)),.(.N,length(unique(gene))),by=.(eGene_type,eQTL_type)]
eQTL_count=eQTL_Signif_both[!duplicated(paste(snps,gene)),.(.N,Ngene=length(unique(gene))),by=.(eGene_type,eQTL_type)][c(1,3,2,5,4),]
eQTL_count
fwrite(eQTL_count,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/Fig3A_eQTLcounts_data.tsv",EVO_IMMUNO_POP_ZEUS),sep='\t')

# eGene_type     eQTL_type    N Ngene
# 1:  lineage_only  lineage_only 1255  1057
# 2:        shared  lineage_only 1736  1373
# 3:        shared        shared 6354  3822
# 4:        shared celltype_only 2835  1907
# 5: celltype_only celltype_only  846   777

# eQTL_count=data.table(x=1:5,N=c(1255,1736,6354,2835,846),Ngene=c(1057,1373,3822,1907,777))

upsetData=eQTL_Signif_both[!duplicated(paste(snps,gene)),.(snps,gene,celltype_only_eQTLs=paste(snps,gene)%in%celltype_only_eQTLs,
                                        lineage_only_eQTLs=paste(snps,gene)%in%lineage_only_eQTLs,
                                        celltype_lineage_shared_eQTLs=paste(snps,gene)%in%celltype_lineage_shared_eQTLs,
                                        celltype_only_eGene=gene%in%celltype_only_eGene,
                                        celltype_lineage_shared_eGene=gene%in%celltype_lineage_shared_eGene,
                                        lineage_only_eGene=gene%in%lineage_only_eGene)]
fwrite(upsetData,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/upsetData.tsv",EVO_IMMUNO_POP_ZEUS),sep='\t')

library(ComplexUpset)
base_annotations=list(
    'Intersection size'=intersection_size(
        text=list(
            hjust=-0.05 ,
            vjust=0.5,
            angle=90
        )))

pdf(sprintf('%s/Fig3A_0Nb_eQTL_by_method.pdf',FIGURE_DIR),height=5.5, width=5)
p <- upset(upsetData, intersect =c('celltype_only_eQTLs','lineage_only_eQTLs','celltype_lineage_shared_eQTLs','celltype_only_eGene','celltype_lineage_shared_eGene','lineage_only_eGene'),
            base_annotations=base_annotations)+theme_yann()
print(p)
dev.off()


####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####
FigS6A_data=upsetData
fwrite(FigS6A_data,file=sprintf('%s/Final/FigS6Adata_Nb_eQTL_by_method.tsv',FIGURE_DIR),sep='\t')

FigS6A=p
saveRDS(FigS6A,file=sprintf('%s/Final/FigS6A_Nb_eQTL_by_method.RDS',FIGURE_DIR))
pdf(sprintf('%s/Final/FigS6A_Nb_eQTL_by_method.pdf',FIGURE_DIR),height=5, width=5)
print(FigS6A)
dev.off()

####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####


Nb_celltype_pereQTL=eQTL_Signif_both[celltype%chin%celltype_22,][,.(N_celltype=length(unique(celltype))),by=.(snps,gene)]
distrib_Nb_celltype_pereQTL=Nb_celltype_pereQTL[,.N,keyby=N_celltype]
distrib_Nb_celltype_pereQTL[,Pct:=100*N/sum(N)]
distrib_Nb_celltype_pereQTL[1:4,]
#    N_celltype    N       Pct
# 1:          1 5876 58.555057
# 2:          2 1652 16.462382
# 3:          3  895  8.918784
# 4:          4  511  5.092177

# eQTL_Signif_celltype_annot=merge(eQTL_Signif_both[celltype%chin%celltype_22,.(gene,snps,celltype,state,top_component,pip_top_component,lbf_top_component)],
#                                   eQTL_Stats_celltype_mash[,.(gene,snps,celltype,state,gene_name,pvalue,beta,se,R2,lfsr,beta_mash,se_mash)],all.y=TRUE)
# eQTL_Signif_celltype_annot[,p.adj:=p.adjust(pvalue,'bonferroni'),by=.(gene,snps)]

# Number_celltype_where_signif_padj=eQTL_Signif_celltype_annot[,.(Number_celltype_where_signif_padj=length(unique(celltype[p.adj<0.01])),
#                                                                 Number_celltype_where_signif_praw1=length(unique(celltype[pvalue<0.01])),
#                                                                 Number_celltype_where_signif_praw01=length(unique(celltype[pvalue<0.001])),
#                                                                 Number_celltype_where_signif_lfsr1=length(unique(celltype[lfsr<0.01])),
#                                                                 Number_celltype_where_signif_lfsr10=length(unique(celltype[lfsr<0.1]))
#                                                                 ),by=.(gene,snps)]
# eQTL_Signif_celltype_annot=merge(eQTL_Signif_celltype_annot,Number_celltype_where_signif_padj,all.x=TRUE)
Nb_celltype_pereQTL_praw1=eQTL_Stats_celltype_mash[,.(Nb_celltype_pereQTL_praw1=length(unique(celltype[pvalue<0.01]))),by=.(snps,gene)]

Nb_celltype_pereQTL_praw1

Nb_lineage_pereQTL_praw1=eQTL_Stats_lineage_mash[,.(Nb_celltype_pereQTL_praw1=length(unique(celltype[pvalue<0.01]))),by=.(snps,gene)]


celltype_specifc_eQTLs=Nb_celltype_pereQTL[N_celltype==1,]
celltype_specifc_eQTLs_annot=merge(unique(eQTL_Signif_both[celltype%chin%celltype_22,.(snps,gene,celltype)]), celltype_specifc_eQTLs, by=c('snps','gene'))
celltype_specifc_eQTLs_annot=merge(celltype_specifc_eQTLs_annot,Nb_celltype_pereQTL_praw1,by=c('snps','gene'))

nrow(celltype_specifc_eQTLs_annot[Nb_celltype_pereQTL_praw1==1,])
# 812
nrow(celltype_specifc_eQTLs_annot[Nb_celltype_pereQTL_praw1==1 & celltype=='MONO.CD14',])
#322
celltype_specifc_eQTLs_annot[Nb_celltype_pereQTL_praw1==1 & celltype%chin%c('MONO.CD14','MONO.CD16','MONO.CD14.INFECTED','cDC','pDC')][,.N,by=celltype][,sum(N)]
# 365

# In which condition are cell stype specific eQTLs found
celltype_specifc_eQTLs_annot_full=merge(unique(eQTL_Signif_both[celltype%chin%celltype_22,]), celltype_specifc_eQTLs, by=c('snps','gene'))
celltype_specifc_eQTLs_annot_full=merge(celltype_specifc_eQTLs_annot_full,Nb_celltype_pereQTL_praw1,by=c('snps','gene'))
celltype_specifc_eQTLs_annot_full[celltype=='MONO.CD14' & Nb_celltype_pereQTL_praw1==1,.N,by=state]
# state   N
# 1:    NS 171
# 2:   COV 161
# 3:   IAV  22

# celltype_2lineage=data.table(celltype=c("B.M.K","B.M.L","B.N.K","B.N.L","Plasmablast","MONO.CD14","MONO.CD14.INFECTED", "MONO.CD16", "cDC","pDC","NK.CD56brt","NK.CD56dim","NK.M.LIKE","T.CD4.E","T.CD4.N","T.Reg","T.CD8.CM.EM","T.CD8.EMRA","T.CD8.N","ILC","MAIT","T.gd"),
  # lineage=c("B","B","B","B","B","MONO","MONO", "MONO", "MONO","MONO","NK","NK","NK","T.CD4","T.CD4","T.CD4","T.CD8","T.CD8","T.CD8","T.CD8","T.CD8","T.CD8"))

eQTL_Stats_celltype_mash=merge(eQTL_Stats_celltype_mash,celltype_2lineage,by='celltype')

Nb_of_eQTL_per_celltype=eQTL_Stats_celltype_mash[celltype%chin%celltype_22,.(N_eQTL_GWsignif=length(unique(snps[!is.na(top_component)])),
                                                                              N_eQTL_GWsignif_CTspecific_praw1=length(unique(snps[!is.na(top_component) & Number_celltype_where_signif_praw1==1])),
#                                                                              N_eQTL_GWsignif_CTspecific_praw01=length(unique(snps[!is.na(top_component) & Number_celltype_where_signif_praw01==1])),
                                                                              N_eQTL_GWsignif_CTspecific_padj=length(unique(snps[!is.na(top_component) & Number_celltype_where_signif_padj==1])),
                                                                              N_eQTL_GWsignif_CTspecific_lfsr1=length(unique(snps[!is.na(top_component) & Number_celltype_where_signif_lfsr1==1])),
#                                                                              N_eQTL_GWsignif_CTspecific_lfsr10=length(unique(snps[!is.na(top_component) & Number_celltype_where_signif_lfsr10==1])),
                                                                              N_eQTL_signif_praw1=length(unique(snps[pvalue<0.01])),
#                                                                              N_eQTL_signif_praw01=length(unique(snps[pvalue<0.001])),
                                                                              N_eQTL_signif_padj=length(unique(snps[p.adj<0.01])),
                                                                              N_eQTL_signif_lfsr1=length(unique(snps[lfsr<0.01]))
#                                                                              N_eQTL_signif_lfsr10=length(unique(snps[lfsr<0.1]))
),by=.(celltype,lineage)]
Nb_of_eQTL_per_celltype=Nb_of_eQTL_per_celltype[order(lineage,-N_eQTL_GWsignif)]
fwrite(Nb_of_eQTL_per_celltype,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/Fig3A_Nb_of_eQTL_per_celltype__data.tsv",EVO_IMMUNO_POP_ZEUS),sep='\t')



################################################################
#########   Fig 3A  eQTL count per cell type     ###############
################################################################

Nb_of_eQTL_per_celltype[,N_eQTL_GWsignif_shared_praw1:=N_eQTL_GWsignif-N_eQTL_GWsignif_CTspecific_praw1]
Nb_of_eQTL_per_celltype[,N_eQTL_signif_praw1_notGW:=N_eQTL_signif_praw1-N_eQTL_GWsignif]
Nb_of_eQTL_per_celltype[,sqrtN_eQTL_GWsignif_CTspecific_praw1:=sqrt(N_eQTL_GWsignif_CTspecific_praw1)]
Nb_of_eQTL_per_celltype[,sqrtN_eQTL_GWsignif_shared_praw1:=sqrt(N_eQTL_GWsignif_shared_praw1)-sqrt(N_eQTL_GWsignif_CTspecific_praw1)]
Nb_of_eQTL_per_celltype[,sqrtN_eQTL_signif_praw1_notGW:=sqrt(N_eQTL_signif_praw1)-sqrt(N_eQTL_GWsignif)]


Nb_of_eQTL_per_celltype_long=melt(Nb_of_eQTL_per_celltype[,.(lineage,celltype,sqrtN_eQTL_GWsignif_CTspecific_praw1,sqrtN_eQTL_GWsignif_shared_praw1,sqrtN_eQTL_signif_praw1_notGW)],id.vars=c('celltype','lineage'))
Nb_of_eQTL_per_celltype_long[,celltype:=factor(celltype,Nb_of_eQTL_per_celltype$celltype)]
Nb_of_eQTL_per_celltype_long[,variable_clear:=case_when(variable=='sqrtN_eQTL_GWsignif_CTspecific_praw1'~"FDR<1%, cell-type-specific",
                                                        variable=='sqrtN_eQTL_GWsignif_shared_praw1'~"FDR<1%, shared",
                                                        variable=='sqrtN_eQTL_signif_praw1_notGW'~"p<0.01, shared",
                                                        TRUE~'debug this')]
Nb_of_eQTL_per_celltype_long[,variable_clear:=factor(variable_clear,c("p<0.01, shared","FDR<1%, shared","FDR<1%, cell-type-specific"))]


library(ggpattern)
#ticks_pos=c(10,50,100,250,500,1000,2500,5000)
ticks_pos=c(10,100,250,500,1000,2500,5000)
p <- ggplot(Nb_of_eQTL_per_celltype_long,aes(x=celltype, fill = celltype, y= value, pattern = variable_clear,alpha=variable_clear)) + scale_y_continuous(breaks=sqrt(ticks_pos),labels=ticks_pos)
p <- p + geom_bar_pattern(position='stack',stat='Identity', color='black', aes(pattern_angle=variable_clear,pattern_color=variable_clear), pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6,width=.7) #color = "black", pattern_fill = "black", pattern_angle = 45,
p <- p + scale_fill_manual(values=celltype_color) + scale_pattern_manual(values = c("FDR<1%, cell-type-specific" = "stripe", "FDR<1%, shared" = "none","p<0.01, shared"="stripe"))
p <- p + ylab('Number of eQTL') + guides(fill= "none") + theme_yann() + theme(legend.title = element_blank()) + theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
p <- p + scale_pattern_angle_manual(values = c("FDR<1%, cell-type-specific" = -45, "FDR<1%, shared" = 45, "p<0.01, shared"=45))
p <- p + scale_pattern_colour_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" = NA,"p<0.01, shared"="black"))
p <- p + scale_pattern_fill_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" =NA,"p<0.01, shared"="black"))
FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3',EVO_IMMUNO_POP_ZEUS)
pdf(sprintf('%s/Fig3A_Nb_eQTL_by_celltype.pdf',FIGURE_DIR),height=6.7*.4, width=7.2*.4)
print(p)
dev.off()


####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####
Fig3B_data=Nb_of_eQTL_per_celltype_long
Fig3B=p+theme(axis.title.x=element_blank(),legend.position="top",text = element_text(size = 8))
fwrite(Fig3B_data,file=sprintf('%s/Final/Fig3data_Nb_eQTL_by_celltype.tsv',FIGURE_DIR),sep='\t')
saveRDS(Fig3B,file=sprintf('%s/Final/Fig3B_Nb_eQTL_by_celltype.RDS',FIGURE_DIR))

pdf(sprintf('%s/Final/Fig3B_Nb_eQTL_by_celltype.pdf',FIGURE_DIR),width=7.2*.4,height=6.7*.4)
print(Fig3B+theme(legend.key.height = unit(.1,"cm"), legend.key.width = unit(.2,"cm")))
dev.off()

####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####


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

eQTL_MIR155HG=get_eQTL('rs114273142','MIR155HG',resolution='celltype',metric='logCPM')
eQTL_MIR155HG[,state:=factor(state,c('NS','COV','IAV'))]
winsorize=0
eQTL_MIR155HG[,logCPM_winsor:=pmax(logCPM,pmin(logCPM,quantile(logCPM,1-winsorize/2)),quantile(logCPM,winsorize/2)),by=.(state,celltype)]
nIQR=3
eQTL_MIR155HG[,logCPM_IQR:=pmax(pmin(logCPM,median(logCPM)+nIQR*IQR(logCPM)),median(logCPM)-nIQR*IQR(logCPM)),by=.(state,celltype)]

eQTL_MIR155HG[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]
irnt=function(x){qnorm(rank(x)/(length(x)+1),mean(x),sd(x))}
eQTL_MIR155HG[,logCPM_irnt:=irnt(logCPM),by=.(state,celltype)]

p <- ggplot(eQTL_MIR155HG[celltype%in%c('pDC',"MONO.CD14"),],aes(x=Number_of_ALT_alelle,y=logCPM_IQR,fill=celltype))
p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=TRUE)
p <- p + theme_yann() + scale_fill_manual(values=celltype_color[c('pDC','MONO.CD14')])# + scale_color_manual(values=color_populations)
p <- p + facet_grid(celltype~state,scales="free_y")

####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####
FigS6B=p+theme(legend.position='none',text = element_text(size = textSize))
FigS6B_data=eQTL_MIR155HG
fwrite(FigS6B_data,file=sprintf('%s/Final/FigS6Bdata_eQTL_MIR155HG_pDC.tsv',FIGURE_DIR))
saveRDS(FigS6B,file=sprintf('%s/Final/FigS6B_eQTL_MIR155HG_pDC.RDS',FIGURE_DIR))
pdf(sprintf('%s/Final/FigS6B_eQTL_MIR155HG_pDC.pdf',FIGURE_DIR),width=7.2*.3,height=6.7*.4)
print(FigS6B)
dev.off()
####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####



#####################################################################################
#########   Fig 3A  eQTL count per cell type - mash based lfsr<.01    ###############
#####################################################################################

Nb_of_eQTL_per_celltype[,N_eQTL_GWsignif_shared_lfsr1:=N_eQTL_GWsignif-N_eQTL_GWsignif_CTspecific_lfsr1]
Nb_of_eQTL_per_celltype[,N_eQTL_signif_lfsr1_notGW:=N_eQTL_signif_lfsr1-N_eQTL_GWsignif]
Nb_of_eQTL_per_celltype[,sqrtN_eQTL_GWsignif_CTspecific_lfsr1:=sqrt(N_eQTL_GWsignif_CTspecific_lfsr1)]
Nb_of_eQTL_per_celltype[,sqrtN_eQTL_GWsignif_shared_lfsr1:=sqrt(N_eQTL_GWsignif_shared_lfsr1)-sqrt(N_eQTL_GWsignif_CTspecific_lfsr1)]
Nb_of_eQTL_per_celltype[,sqrtN_eQTL_signif_lfsr1_notGW:=sqrt(N_eQTL_signif_lfsr1)-sqrt(N_eQTL_GWsignif)]


Nb_of_eQTL_per_celltype_long=melt(Nb_of_eQTL_per_celltype[,.(lineage,celltype,sqrtN_eQTL_GWsignif_CTspecific_lfsr1,sqrtN_eQTL_GWsignif_shared_lfsr1,sqrtN_eQTL_signif_lfsr1_notGW)],id.vars=c('celltype','lineage'))
Nb_of_eQTL_per_celltype_long[,celltype:=factor(celltype,Nb_of_eQTL_per_celltype$celltype)]
Nb_of_eQTL_per_celltype_long[,variable_clear:=case_when(variable=='sqrtN_eQTL_GWsignif_CTspecific_lfsr1'~"FDR<1%, cell-type-specific",
                                                        variable=='sqrtN_eQTL_GWsignif_shared_lfsr1'~"FDR<1%, shared",
                                                        variable=='sqrtN_eQTL_signif_lfsr1_notGW'~"p<0.01, shared",
                                                        TRUE~'debug this')]
Nb_of_eQTL_per_celltype_long[,variable_clear:=factor(variable_clear,c("p<0.01, shared","FDR<1%, shared","FDR<1%, cell-type-specific"))]


library(ggpattern)
ticks_pos=c(10,50,100,250,500,1000,2500,5000)
p <- ggplot(Nb_of_eQTL_per_celltype_long,aes(x=celltype, fill = celltype, y= value, pattern = variable_clear,alpha=variable_clear)) + scale_y_continuous(breaks=sqrt(ticks_pos),labels=ticks_pos)
p <- p + geom_bar_pattern(position='stack',stat='Identity', color='black', aes(pattern_angle=variable_clear,pattern_color=variable_clear), pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6,width=.7) #color = "black", pattern_fill = "black", pattern_angle = 45,
p <- p + scale_fill_manual(values=celltype_color) + scale_pattern_manual(values = c("FDR<1%, cell-type-specific" = "stripe", "FDR<1%, shared" = "none","p<0.01, shared"="stripe"))
p <- p + ylab('Number of eQTL') + guides(fill= "none") + theme_yann() + theme(legend.title = element_blank()) + theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
p <- p + scale_pattern_angle_manual(values = c("FDR<1%, cell-type-specific" = -45, "FDR<1%, shared" = 45, "p<0.01, shared"=45))
p <- p + scale_pattern_colour_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" = NA,"p<0.01, shared"="black"))
p <- p + scale_pattern_fill_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" =NA,"p<0.01, shared"="black"))
FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3',EVO_IMMUNO_POP_ZEUS)
pdf(sprintf('%s/Fig3Abis_Nb_eQTL_by_celltype_lfsr1.pdf',FIGURE_DIR),height=5.5, width=6)
print(p)
dev.off()


##########################################################################
#########   Fig 3B sharing of eQTLs between cell type      ###############
##########################################################################

# Beta_cor=cor(beta_mat[,-(1:4)],use='p')[celltype_2lineage$celltype,celltype_2lineage$celltype]
#
# library(ggcorrplot)
# pdf(sprintf('%s/3Bv1_Corrplot_pearson_eQTLbeta_Signif.pdf',FIGURE_DIR),height=9,width=9)
# p <- ggcorrplot(Beta_cor, method = "circle", hc.order = FALSE)
# print(p)
# dev.off()
#
# Beta_cor_NS=cor(beta_mat[state=='NS',-(1:4)],use='p')[celltype_2lineage$celltype,celltype_2lineage$celltype]
# Beta_cor_IAV=cor(beta_mat[state=='IAV',-(1:4)],use='p')[celltype_2lineage$celltype,celltype_2lineage$celltype]
# Beta_cor_COV=cor(beta_mat[state=='COV',-(1:4)],use='p')[celltype_2lineage$celltype,celltype_2lineage$celltype]


get_pairwise_sharing_raw=function (beta_mat, FDR_mat, factor = 0.5, FDR_thresh = 0.01,cor=FALSE,method='s'){
    R = ncol(beta_mat)
    S = matrix(NA, nrow = R, ncol = R)
    A = matrix(NA, nrow = R, ncol = R)
    for (i in 1:R) {
        for (j in i:R) {
            if(!all(is.na(beta_mat[, i]) | is.na(beta_mat[, j])) ){
              a = FDR_mat[,i]<FDR_thresh | FDR_mat[,j]<FDR_thresh
              if(cor){
              S[i, j] = cor(beta_mat[a, i],beta_mat[a, j],method=method)
              A[i, j]=sum(a)
               if(cor.test(beta_mat[a, i],beta_mat[a, j],method=method)$p.value>0.01){
                  S[i, j] = 0
                }
              }else{
              ratio = beta_mat[a, i]/beta_mat[a, j]
              S[i, j] = mean(ratio > factor & ratio < (1/factor))
              }
            }
        }
    }
    S[lower.tri(S, diag = FALSE)] = t(S)[lower.tri(S, diag = FALSE)]
    colnames(S) = row.names(S) = colnames(beta_mat)
    return(S)
}

celltype_2lineage_ordered=data.table(celltype=c("B.M.K","B.M.L","B.N.K","B.N.L","Plasmablast","MONO.CD14","MONO.CD14.INFECTED", "MONO.CD16", "cDC","pDC","NK.CD56brt","NK.CD56dim","NK.M.LIKE","T.CD8.EMRA","MAIT","T.CD8.CM.EM","T.CD8.N","T.CD4.E","T.CD4.N","T.Reg","T.gd","ILC"),
lineage=c("B","B","B","B","B","MONO","MONO", "MONO", "MONO","MONO","NK","NK","NK","T.CD8","T.CD8","T.CD8","T.CD8","T.CD4","T.CD4","T.CD4","T.CD8","T.CD8"))
celltype_2lineage_noINF=celltype_2lineage_ordered[celltype!='MONO.CD14.INFECTED',]

beta_mat=dcast(eQTL_Stats_celltype_mash,snps+gene+gene_name+state~celltype,value.var='beta')
beta_mat_NS=as.matrix(beta_mat[state=='NS',-c('snps','gene','gene_name','state')])
beta_mat_COV=as.matrix(beta_mat[state=='COV',-c('snps','gene','gene_name','state')])
beta_mat_IAV=as.matrix(beta_mat[state=='IAV',-c('snps','gene','gene_name','state')])

pval_mat=dcast(eQTL_Stats_celltype_mash,snps+gene+gene_name+state~celltype,value.var='pvalue')
pval_mat_NS=as.matrix(pval_mat[state=='NS',-c('snps','gene','gene_name','state')])
pval_mat_COV=as.matrix(pval_mat[state=='COV',-c('snps','gene','gene_name','state')])
pval_mat_IAV=as.matrix(pval_mat[state=='IAV',-c('snps','gene','gene_name','state')])

#sharing_raw_NS=get_pairwise_sharing_raw(beta_mat_NS,pval_mat_NS)
sharing_cor_NS=get_pairwise_sharing_raw(beta_mat_NS,pval_mat_NS,cor=TRUE,FDR_thresh=0.01)
sharing_cor_IAV=get_pairwise_sharing_raw(beta_mat_IAV,pval_mat_IAV,cor=TRUE,FDR_thresh=0.01)
sharing_cor_COV=get_pairwise_sharing_raw(beta_mat_COV,pval_mat_COV,cor=TRUE,FDR_thresh=0.01)

ann_colors = list(
    celltype = celltype_color,
    lineage = color_cellTypes_6level[lineage_5])

library(ComplexHeatmap)
pdf(sprintf('%s/Fig3Bv6_eQTLsharing_raw_NS.pdf',FIGURE_DIR),height=5,width=5)
pheatmap(sharing_cor_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype] ,cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
  annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
dev.off()





################# TODO : add sharing of log FC on the lower diagonal: to show that logFC are more celltype-specific than eQTLs
# pdf(sprintf('%s/3Bv6_eQTLsharing_raw_NS_clustered.pdf',FIGURE_DIR),height=5,width=5,)
# pheatmap(sharing_cor_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=TRUE,cluster_rows=TRUE,fontsize=8,
#   annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
# dev.off()


#### if we switch to spearman corr
COR_comp=as.data.table(reshape2::melt(sharing_cor_NS,varnames=c('celltype_1','celltype_2'),value.name='correlation'))
COR_comp=merge(COR_comp,celltype_2lineage,by.x='celltype_1',by.y='celltype')
COR_comp=merge(COR_comp,celltype_2lineage,by.x='celltype_2',by.y='celltype',suffix=c('_1','_2'))
COR_comp=COR_comp[celltype_2!=celltype_1,]
COR_comp[,comb:=paste(sort(c(celltype_1,celltype_2)),collapse='_'),by=.(celltype_1,celltype_2)]
COR_comp=COR_comp[!duplicated(comb),]
COR_comp[,mean(correlation,na.rm=T),by=.(lineage_1==lineage_2)]


# valuers are for pearson correlation . p=0.02 with spearman -> consider meta lineage instead
#  lineage_1        V1
# 1:      TRUE 0.5901654
# 2:     FALSE 0.4614593
wilcox.test(COR_comp[lineage_1==lineage_2,correlation],COR_comp[lineage_1!=lineage_2,correlation])
# W = 4353, p-value = 0.0005951


COR_comp[,meta_lineage_1:=ifelse(lineage_1%in%c('T.CD4','T.CD8','NK'),'TNK',lineage_1)]
COR_comp[,meta_lineage_2:=ifelse(lineage_2%in%c('T.CD4','T.CD8','NK'),'TNK',lineage_2)]
COR_comp[,mean(correlation,na.rm=T),by=.(meta_lineage_1==meta_lineage_2)]
wilcox.test(COR_comp[meta_lineage_1 ==meta_lineage_2,correlation],COR_comp[meta_lineage_1!=meta_lineage_2,correlation])
# meta_lineage_1        V1
# 1:           TRUE 0.6003208
# 2:          FALSE 0.4737219
# > wilcox.test(COR_comp[meta_lineage_1 ==meta_lineage_2,correlation],COR_comp[meta_lineage_1!=meta_lineage_2,correlation])
#
#      Wilcoxon rank sum test with continuity correction
#
# data:  COR_comp[meta_lineage_1 == meta_lineage_2, correlation] and COR_comp[meta_lineage_1 != meta_lineage_2, correlation]
# W = 7190, p-value = 6.203e-06
# alternative hypothesis: true location shift is not equal to 0


######### check sharing of reQTL
beta_mat=dcast(reQTL_Stats_celltype_mash,snps+gene+gene_name+state~celltype,value.var='beta')
beta_mat_COV=as.matrix(beta_mat[state=='COV',-c('snps','gene','gene_name','state')])
beta_mat_IAV=as.matrix(beta_mat[state=='IAV',-c('snps','gene','gene_name','state')])

pval_mat=dcast(reQTL_Stats_celltype_mash,snps+gene+gene_name+state~celltype,value.var='pvalue')
pval_mat_COV=as.matrix(pval_mat[state=='COV',-c('snps','gene','gene_name','state')])
pval_mat_IAV=as.matrix(pval_mat[state=='IAV',-c('snps','gene','gene_name','state')])

#sharing_raw_NS=get_pairwise_sharing_raw(beta_mat_NS,pval_mat_NS)
sharing_cor_COV=get_pairwise_sharing_raw(beta_mat_COV,pval_mat_COV,cor=TRUE,FDR_thresh=0.01)
sharing_cor_IAV=get_pairwise_sharing_raw(beta_mat_IAV,pval_mat_IAV,cor=TRUE,FDR_thresh=0.01)

ann_colors = list(
    celltype = celltype_color,
    lineage = lineage_color)

pdf(sprintf('%s/Fig3Bv6_lower_reQTLsharing_raw_COV.pdf',FIGURE_DIR),height=5,width=5)
pheatmap(sharing_cor_COV[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
  annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
dev.off()

sharing_cor_NSCOV=sharing_cor_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype]
for (i in 1:(nrow(celltype_2lineage_noINF)-1)){
  for (j in (i+1):nrow(celltype_2lineage_noINF)){
    sharing_cor_NSCOV[i,j]=sharing_cor_COV[celltype_2lineage_noINF$celltype[i],celltype_2lineage_noINF$celltype[j]]
  }
}

# compare correlation eQTL VS reQTL
COR_comp2=as.data.table(reshape2::melt(sharing_cor_NSCOV,varnames=c('celltype_1','celltype_2'),value.name='correlation'))
COR_comp2=merge(COR_comp2,cbind(celltype_2lineage,N=1:nrow(celltype_2lineage)),by.x='celltype_1',by.y='celltype')
COR_comp2=merge(COR_comp2,cbind(celltype_2lineage,N=1:nrow(celltype_2lineage)),by.x='celltype_2',by.y='celltype',suffix=c('_1','_2'))
COR_comp2=COR_comp2[celltype_2!=celltype_1,]
COR_comp[,mean(correlation,na.rm=T),by=.(lineage_1==lineage_2)]
COR_comp2[,COND:=ifelse(N_1<N_2,'COV','NS')]
# COND        V1
# 1:   NS 0.5082113
# 2:  COV 0.3654303
COR_comp2[,wilcox.test(correlation[COND=='COV'],correlation[COND=='NS'])$p.value]
Fig3C_data=data.table(sharing_cor_NSCOV,celltype_2lineage_noINF,celltype_color=celltype_color[celltype_2lineage_noINF$celltype],lineage_color=lineage_color[celltype_2lineage_noINF$lineage])

celltype_2lineage_noINF=Fig3C_data[,.(celltype,lineage)]
sharing=as.matrix(Fig3C_data[,mget(celltype_2lineage_noINF$celltype)])
rownames(sharing)=celltype_2lineage_noINF$celltype
ann_colors=list(celltype=Fig3C_data[,setNames(celltype_color,celltype)],lineage=Fig3C_data[!duplicated(lineage),setNames(lineage_color,lineage)])

pdf(sprintf('%s/Fig3Bv7_eQTL_NS_reQTL_COV_sharing.pdf',FIGURE_DIR),height=4,width=4)
Fig3C <- ComplexHeatmap::pheatmap(sharing,cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
  annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE,legend=FALSE)
dev.off()

pdf(sprintf('%s/Fig3Bv7_eQTL_NS_reQTL_COV_sharing_with_legend.pdf',FIGURE_DIR),height=4,width=4)
Fig3C <- ComplexHeatmap::pheatmap(sharing,cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
  annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE,legend=TRUE)
dev.off()

library(circlize)
library(RColorBrewer)
pdf(sprintf('%s/Fig3Bv7_eQTL_NS_reQTL_COV_sharing_complexHeatmaps.pdf',FIGURE_DIR),height=4,width=4)
Fig3C <- ComplexHeatmap::Heatmap(sharing,cluster_columns=FALSE,cluster_rows=FALSE,column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col='darkgrey'),
                                col = circlize::colorRamp2(seq(0,1,l=7), rev(brewer.pal(n = 7, name = "RdYlBu"))),
                                top_annotation = HeatmapAnnotation(df = celltype_2lineage_noINF,col=ann_colors,show_legend=FALSE,simple_anno_size = unit(.3, "cm"),gp = gpar(col='grey'),row_names_gp = gpar(fontsize = 8)),
                                left_annotation = rowAnnotation(df = celltype_2lineage_noINF,col=ann_colors,show_legend=FALSE,simple_anno_size = unit(.3, "cm"),gp = gpar(col='grey'),column_names_gp = gpar(fontsize = 8)),
                                show_heatmap_legend=FALSE)
draw(Fig3C)
dev.off()

####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####
fwrite(Fig3C_data,file=sprintf('%s/Final/Fig3Cdata_eQTLsharing_raw_NSCOV.tsv',FIGURE_DIR),sep='\t')
saveRDS(Fig3C,file=sprintf('%s/Final/Fig3C_eQTLsharing_raw_NSCOV.RDS',FIGURE_DIR))

pdf(sprintf('%s/Final/Fig3C_eQTLsharing_raw_NSCOV.pdf',FIGURE_DIR),height=6.7*.4,width=7.2*.38)
Fig3C <- ComplexHeatmap::Heatmap(sharing,cluster_columns=FALSE,cluster_rows=FALSE,column_names_gp = gpar(fontsize = 7),row_names_gp = gpar(fontsize = 7),rect_gp = gpar(col='darkgrey'),
                                col = circlize::colorRamp2(seq(0,1,l=7), rev(brewer.pal(n = 7, name = "RdYlBu"))),
                                top_annotation = HeatmapAnnotation(df = celltype_2lineage_noINF,col=ann_colors,show_legend=FALSE,simple_anno_size = unit(.2, "cm"),gp = gpar(col='grey')),
                                left_annotation = rowAnnotation(df = celltype_2lineage_noINF,col=ann_colors,show_legend=FALSE,simple_anno_size = unit(.2, "cm"),gp = gpar(col='grey')),
                                show_heatmap_legend=FALSE)

draw(Fig3C)
dev.off()

####------------------------------------------------------------------------------------####
####------------------------------------------------------------------------------------####


#################################################################################
#########   Fig 3B sharing of eQTLs between cell type - mash      ###############
#################################################################################

m2=readRDS(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/celltype_condition___CellPropLineage_SVs_220409/dist_100kb/mash_estimated_params_observed_independent_eQTLs_allcond_1pctFDR_allStats.RDS',EVO_IMMUNO_POP_ZEUS))

pi_sharing=get_pairwise_sharing(m2, lfsr_thresh = 0.01)
pi_sharing[grep('__NS',rownames(pi_sharing)),grep('__NS',colnames(pi_sharing))]

pi_mash_NS=pi_sharing[grep('__NS',rownames(pi_sharing)),grep('__NS',colnames(pi_sharing))]
colnames(pi_mash_NS)=gsub('(.*)__.*','\\1',colnames(pi_mash_NS))
rownames(pi_mash_NS)=gsub('(.*)__.*','\\1',rownames(pi_mash_NS))

pdf(sprintf('%s/3Bv5_eQTLsharing_mash_NS.pdf',FIGURE_DIR),height=5,width=5)
pheatmap(pi_mash_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
  annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
dev.off()

m2_reQTL=readRDS(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/lineage_condition_logFC__logFC__CellPropLineage_SVs_220409/dist_100kb/mash_estimated_params_observed_independent_eQTLs_allcond_celltype_and_lineageLevel_1pctFDR_allStatscelltype.RDS',EVO_IMMUNO_POP_ZEUS))

pi_sharing=get_pairwise_sharing(m2_reQTL, lfsr_thresh = 0.01)
pi_mash_COV=pi_sharing[grep('__COV',rownames(pi_sharing)),grep('__COV',colnames(pi_sharing))]
colnames(pi_mash_COV)=gsub('(.*)__.*','\\1',colnames(pi_mash_COV))
rownames(pi_mash_COV)=gsub('(.*)__.*','\\1',rownames(pi_mash_COV))

pdf(sprintf('%s/3Bv5_lower_reQTLsharing_mash_COV.pdf',FIGURE_DIR),height=5,width=5)
pheatmap(pi_mash_COV[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
  annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
dev.off()

########################################################################################
#########   Fig 3C : comparison of response eQTLs between conditions     ###############
########################################################################################

# see reQTL_and_ISG_v2.R



#    p=0.0005951
#
# pdf(sprintf('%s/3Bv6_eQTLsharing_raw_IAV.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(sharing_cor_IAV[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
#   annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
# dev.off()
#
# pdf(sprintf('%s/3Bv6_eQTLsharing_raw_COV.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(sharing_cor_COV[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
#   annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
# dev.off()
#
# pdf(sprintf('%s/3Bv5_eQTLsharing_mash_IAV.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(pi_mash_IAV[celltype_2lineage$celltype,celltype_2lineage$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8)
# dev.off()
#
# pdf(sprintf('%s/3Bv5_eQTLsharing_mash_COV.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(pi_mash_COV[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8)
# dev.off()
#
#
# COR_comp=as.data.table(reshape2::melt(pi_mash_NS,varnames=c('celltype_1','celltype_2'),value.name='correlation'))
# COR_comp=merge(COR_comp,celltype_2lineage,by.x='celltype_1',by.y='celltype')
# COR_comp=merge(COR_comp,celltype_2lineage,by.x='celltype_2',by.y='celltype',suffix=c('_1','_2'))
# COR_comp=COR_comp[celltype_2!=celltype_1,]
# COR_comp[,comb:=paste(sort(c(celltype_1,celltype_2)),collapse='_'),by=.(celltype_1,celltype_2)]
# COR_comp=COR_comp[!duplicated(comb),]
# wilcox.test(COR_comp[lineage_1==lineage_2,correlation],COR_comp[lineage_1!=lineage_2,correlation])

#
# celltype_2lineage=data.table(celltype=c("B.M.K","B.M.L","B.N.K","B.N.L","Plasmablast","MONO.CD14","MONO.CD14.INFECTED", "MONO.CD16", "cDC","pDC","NK.CD56brt","NK.CD56dim","NK.M.LIKE","T.CD8.EMRA","MAIT","T.CD8.CM.EM","T.CD8.N","T.CD4.E","T.CD4.N","T.Reg","T.gd","ILC"),
# lineage=c("B","B","B","B","B","MONO","MONO", "MONO", "MONO","MONO","NK","NK","NK","T.CD8","T.CD8","T.CD8","T.CD8","T.CD4","T.CD4","T.CD4","T.CD8","T.CD8"))
#
#
# pdf(sprintf('%s/3Bv2_Corrplot_pearson_eQTLbeta_NS_Signif.pdf',FIGURE_DIR),height=9,width=9)
# p <- ggcorrplot(Beta_cor_NS, method = "circle", hc.order = FALSE)
# print(p)
# dev.off()
# pdf(sprintf('%s/3Bv2_Corrplot_pearson_eQTLbeta_NS_Signif_pheatmap.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(Beta_cor_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8)
# dev.off()
#
# pdf(sprintf('%s/3Bv3_Corrplot_pearson_eQTLbeta_COV_Signif.pdf',FIGURE_DIR),height=9,width=9)
# p <- ggcorrplot(Beta_cor_COV, method = "circle", hc.order = FALSE)
# print(p)
# dev.off()
#
# pdf(sprintf('%s/3Bv3_Corrplot_pearson_eQTLbeta_IAV_Signif.pdf',FIGURE_DIR),height=9,width=9)
# p <- ggcorrplot(Beta_cor_IAV, method = "circle", hc.order = FALSE)
# print(p)
# dev.off()
#
# pval_mat=dcast(eQTL_Signif_celltype_annot,snps+gene+gene_name+state~celltype,value.var='pvalue')
# m2=readRDS(sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/celltype_condition___CellPropLineage_SVs_220409/dist_100kb/mash_estimated_params_observed_independent_eQTLs_allcond_1pctFDR_allStats.RDS',EVO_IMMUNO_POP_ZEUS))
#
# pi_sharing=get_pairwise_sharing(m2, lfsr_thresh = 0.01)
# pi_sharing[grep('__NS',rownames(pi_sharing)),grep('__NS',colnames(pi_sharing))]
#
# pval_mat_2=merge(eQTL_Signif_celltype_annot[!is.na(top_component),.(snps,gene,gene_name,state,celltype)],eQTL_Signif_celltype_annot[,.(snps,gene,gene_name,state,celltype,pvalue)],by=c('snps','gene','gene_name','state'),suffix=c('.discovery','.replication'))
#
# pi1_estim=pval_mat_2[,.(pi1_p01=mean(pvalue<0.01),pi1_50=pmin(pmax(1-mean(pvalue>0.5)/0.5,0),1)),by=.(celltype.discovery,celltype.replication,state)]
#
# pi1_50_NS=dcast(pi1_estim[state=='NS',],celltype.discovery~celltype.replication,value.var='pi1_50')
# ct=pi1_50_NS[,celltype.discovery]
# pi1_50_NS=as.matrix(pi1_50_NS[,-'celltype.discovery'])
# rownames(pi1_50_NS)=ct
# pdf(sprintf('%s/3Bv4_eQTLsharing_pi1_50_NS.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(pi1_50_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8)
# dev.off()
#
#
# pi1_50_COV=dcast(pi1_estim[state=='COV',],celltype.discovery~celltype.replication,value.var='pi1_50')
# pi1_50_COV=dcast(pi1_estim[state=='IAV',],celltype.discovery~celltype.replication,value.var='pi1_50')
#
# suppressMessages(library(ComplexHeatmap))
# pdf(sprintf('%s/3Bv4_eQTLsharing_pi1_50_NS.pdf',FIGURE_DIR),height=9,width=9)
# pheatmap(pi1_50_NS[,-1],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8)
# dev.off()
#
# pi1_p01_NS=dcast(pi1_estim[state=='NS',],celltype.discovery~celltype.replication,value.var='pi1_p01')
# pi1_p01_COV=dcast(pi1_estim[state=='COV',],celltype.discovery~celltype.replication,value.var='pi1_p01')
# pi1_p01_COV=dcast(pi1_estim[state=='IAV',],celltype.discovery~celltype.replication,value.var='pi1_p01')
#
# pi_mash_NS=pi_sharing[grep('__NS',rownames(pi_sharing)),grep('__NS',colnames(pi_sharing))]
# colnames(pi_mash_NS)=gsub('(.*)__.*','\\1',colnames(pi_mash_NS))
# rownames(pi_mash_NS)=gsub('(.*)__.*','\\1',rownames(pi_mash_NS))
#
# pi_mash_COV=pi_sharing[grep('__COV',rownames(pi_sharing)),grep('__COV',colnames(pi_sharing))]
# colnames(pi_mash_COV)=gsub('(.*)__.*','\\1',colnames(pi_mash_COV))
# rownames(pi_mash_COV)=gsub('(.*)__.*','\\1',rownames(pi_mash_COV))
#
# pi_mash_IAV=pi_sharing[grep('__IAV',rownames(pi_sharing)),grep('__IAV',colnames(pi_sharing))]
# colnames(pi_mash_IAV)=gsub('(.*)__.*','\\1',colnames(pi_mash_IAV))
# rownames(pi_mash_IAV)=gsub('(.*)__.*','\\1',rownames(pi_mash_IAV))
#
# pdf(sprintf('%s/3Bv5_eQTLsharing_mash_NS.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(pi_mash_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
#   annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
# dev.off()
#
# pdf(sprintf('%s/3Bv5_eQTLsharing_mash_NS_clustered.pdf',FIGURE_DIR),height=5,width=5,)
# pheatmap(pi_mash_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=TRUE,cluster_rows=TRUE,fontsize=8,
#   annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE)
# dev.off()
#
# pdf(sprintf('%s/3Bv5_eQTLsharing_mash_IAV.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(pi_mash_IAV[celltype_2lineage$celltype,celltype_2lineage$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8)
# dev.off()
#
# pdf(sprintf('%s/3Bv5_eQTLsharing_mash_COV.pdf',FIGURE_DIR),height=5,width=5)
# pheatmap(pi_mash_COV[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype],cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8)
# dev.off()


# Nb_of_eQTL_per_celltype_long=melt(Nb_of_eQTL_per_celltype[,.(lineage,celltype,N_eQTL_GWsignif_CTspecific_praw1,N_eQTL_GWsignif_shared_praw1,N_eQTL_signif_praw1_notGW)],id.vars=c('celltype','lineage'))
# Nb_of_eQTL_per_celltype_long[,celltype:=factor(celltype,Nb_of_eQTL_per_celltype$celltype)]
# Nb_of_eQTL_per_celltype_long[,variable_clear:=case_when(variable=='N_eQTL_GWsignif_CTspecific_praw1'~"FDR<1%, cell-type-specific",
#                                                         variable=='N_eQTL_GWsignif_shared_praw1'~"FDR<1%, shared",
#                                                         variable=='N_eQTL_signif_praw1_notGW'~"p<0.01, shared",
#                                                         TRUE~'debug this')]
# Nb_of_eQTL_per_celltype_long[,variable_clear:=factor(variable_clear,c("p<0.01, shared","FDR<1%, shared","FDR<1%, cell-type-specific"))]
#
#
# library(ggpattern)
# p <- ggplot(Nb_of_eQTL_per_celltype_long,aes(x=celltype, fill = celltype, y= value, pattern = variable_clear,alpha=variable_clear)) + scale_y_continuous(trans='sqrt')#+ scale_y_continuous(breaks=sqrt(10,50,100,500,1000,5000,10000))
# p <- p + geom_bar_pattern(position='dodge',stat='Identity', color = "black", pattern_fill = "black", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6)
# p <- p + scale_fill_manual(values=celltype_color) + scale_pattern_manual(values = c("FDR<1%, cell-type-specific" = "stripe", "FDR<1%, shared" = "none", "p<0.01, shared"="stripe"))
# p <- p + scale_pattern_angle_manual(values = c("FDR<1%, cell-type-specific" = 0, "FDR<1%, shared" = 45, "p<0.01, shared"=45))
# p <- p + scale_pattern_colour_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" = "circle","p<0.01, shared"="black"))
# p <- p + scale_pattern_fill_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" = NA,"p<0.01, shared"="black"))
# p <- p + theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1)) + ylab('Number of eQTL')#+scale_alpha_discrete()
# FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3',EVO_IMMUNO_POP_ZEUS)
# pdf(sprintf('%s/Fig3A_Nb_eQTL_by_celltype.pdf',FIGURE_DIR))
# print(p)
# dev.off()


# ggplot(data = df, aes(x = Class, fill = StudyTime, pattern = Nerd)) +
#   geom_bar_pattern(position = position_dodge(preserve = "single"),
#                    color = "black",
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.1,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 0.6) +
#   scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(4)) +
#   scale_pattern_manual(values = c(Nerd = "stripe", NotNerd = "none")) +
#   labs(x = "Class", y = "Number of Students", pattern = "Nerd?") +
#   guides(pattern = guide_legend(override.aes = list(fill = "white")),
#          fill = guide_legend(override.aes = list(pattern = "none")))


# Pct_of_all_eQTL_under_H1_celltype=eQTL_Signif_celltype_annot[celltype%chin%celltype_22, .(pi1=1-mean(pvalue>.5)/0.5),by=celltype]
