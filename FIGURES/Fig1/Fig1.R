library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

EVO_IMMUNO_POP_ZEUS="/pasteur/zeus/projets/p02/evo_immuno_pop"
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EVO_IMMUNO_POP_ZEUS))
DATA_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data"
FIG_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/paper_draft/V2/figureMaterial"
FIGURE_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/paper_draft/V2/figureMaterial/Fig1"

################################################################################
# Fig. 1A

# Scheme

################################################################################
# Fig. 1B

meta_clean=fread(sprintf("%s/2_population_differences/sce_clean_metadata.tsv",DATA_DIR))

fig1b_data=meta_clean[,.(Barcode_lib,umap1,umap2,COND,lineage,celltype)]

celltype_color=c(
"T.CD4.N"="#169FD8","T.CD4.E"="#005274","T.Reg"="#03C2C0","T.gd"="#ffce73",
"T.CD8.N"="#00BB54","T.CD8.CM.EM"="#0EAD20","T.CD8.EMRA"="#048c41","MAIT"="#004A21",
"NK.M.LIKE"="#63C425","NK.CD56dim"="#B7DC2A","NK.CD56brt"="#7D9E00","ILC"="#3F4F00",
"B.N.K"="#9281DE","B.N.L"="#a681de","B.M.K"="#6045d9","B.M.L"="#8045d9",
"Plasmablast"="#6801A4","MONO.CD14"="#B7637E","MONO.CD16"="#8F568A","pDC"="#8E6B00",
"cDC"="#E7b315","MONO.CD14.INFECTED"="#B9364B")

celltype_order=c("MONO.CD14","MONO.CD16","MONO.CD14.INFECTED","cDC","pDC","B.N.K","B.N.L","B.M.K","B.M.L","Plasmablast","T.CD4.N","T.CD4.E","T.Reg","T.gd","T.CD8.N","T.CD8.CM.EM","T.CD8.EMRA","ILC","MAIT","NK.CD56dim","NK.CD56brt","NK.M.LIKE")
celltype_label=c("MONO CD14+","MONO CD16+","MONO IAV+","cDC","pDC","B N k","B N l","B M k","B M l","Plasmablast","T CD4+ N","T CD4+ E","T Reg","T gd","T CD8+ N","T CD8+ CM/EM","T CD8+ EMRA","ILC","MAIT","NK CD56dim","NK CD56brt","NK mem")

fig1b_plot=ggplot(fig1b_data,aes(x=umap1,y=umap2))+
  theme_yann() +
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
  theme(text=element_text(size=7),legend.position="right",legend.spacing.y=unit(-5,"pt"))

pname=sprintf("%s/Fig1/Fig1b.pdf",FIG_DIR)
pdf(pname,width=6,height=5)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig1b_plot+theme(legend.position="none",text=element_text(size=10))),
      get_legend(fig1b_plot)#,
#      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
#      c(3,1,3),
      c(1,2)
    ), widths=c(7,3)#, widths=c(2,6,2)
  )
dev.off()

################################################################################
# Fig. 1C
cond_order=c("NS","COV","IAV")

fig1c_data=fig1b_data
write_tsv(fig1c_data,file=sprintf("%s/Fig1C_data__ISG_activity_byCOND.tsv",FIGURE_DIR))

nlevel=2

fig1c_plot=ggplot(fig1c_data%>%mutate(COND,factor(COND,cond_order)),aes(x=umap1,y=umap2))+
  theme_yann() +
  scale_fill_distiller('Blues',direction = 1) +
  stat_density2d(data=fig1c_data[COND=="IAV",.(umap1,umap2,COND="IAV")],mapping=aes(fill = stat(nlevel)), geom = "polygon", n = 100, bins = 10, fill="#6699CC",size=0.1) +
  stat_density2d(data=fig1c_data[COND=="NS",.(umap1,umap2,COND="NS")],mapping=aes(fill = stat(nlevel)), geom = "polygon", n = 100, bins = 10, fill="#888888",size=0.1) +
  stat_density2d(data=fig1c_data[COND=="COV",.(umap1,umap2,COND="COV")],mapping=aes(fill = stat(nlevel)), geom = "polygon", n = 100, bins = 10, fill="#DD5555",size=0.1) +
  stat_density2d(data=fig1c_data[,.(umap1,umap2)], fill='lightgrey', alpha=0, geom = "polygon", n = 100, bins = 10, colour="#000000AA",size=0.1) +
  xlab('UMAP1')+ylab('UMAP2')+facet_grid(cols=vars(factor(COND,cond_order)))+
  theme(text=element_text(size=7),panel.spacing=unit(0,"pt"))

################################################################################
################################################################################
#################  overall ISG and inflammatory score   ########################
################################################################################
################################################################################

  EXPRESSION_FILE="BatchAdjusted_logCPM_125libs__per_condition_annotated.tsv.gz"

  Expression=fread(sprintf("%s/single_cell/project/pop_eQTL/data/2_population_differences/%s",EVO_IMMUNO_POP_ZEUS,EXPRESSION_FILE))
  Expression_matrix=dcast(Expression[IID!="ASH071"],ID~IID+state,value.var="logCPM")%>%as.matrix(rownames="ID")
  seurat=CreateSeuratObject(counts=Expression_matrix)
  seurat$IID=str_split(colnames(seurat),"_",simplify=T)[,1]
  seurat$condition=str_split(colnames(seurat),"_",simplify=T)[,2]

  REF_DIR=sprintf("%s/single_cell/resources/references/PATHWAYS",EVO_IMMUNO_POP_ZEUS)
  feature_data=Expression[,.(ID,Symbol)]%>%unique

  hallmark_ifng_response_genes=fread(sprintf("%s/HALLMARK_INTERFERON_GAMMA_RESPONSE.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(feature_data[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
  names(hallmark_ifng_response_genes)=unique(feature_data[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_ifng_response_genes]
  hallmark_ifna_response_genes=fread(sprintf("%s/HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(feature_data[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
  names(hallmark_ifna_response_genes)=unique(feature_data[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_ifna_response_genes]
  hallmark_ifn_response_genes=unique(c(hallmark_ifng_response_genes,hallmark_ifna_response_genes))
  names(hallmark_ifn_response_genes)=unique(feature_data[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_ifn_response_genes]

  hallmark_inflammatory_response_genes=fread(sprintf("%s/HALLMARK_INFLAMMATORY_RESPONSE.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(feature_data[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
  names(hallmark_inflammatory_response_genes)=unique(feature_data[,.(Symbol,ID)])[,setNames(ID,Symbol)][hallmark_inflammatory_response_genes]

  hallmark_inflammatory_exclusive=which(!hallmark_inflammatory_response_genes%in%hallmark_ifn_response_genes)%>%hallmark_inflammatory_response_genes[.]
  hallmark_ifn_exclusive=which(!hallmark_ifn_response_genes%in%hallmark_inflammatory_response_genes)%>%hallmark_ifn_response_genes[.]

  go_ifng_response_genes=fread(sprintf("%s/GO_0034341_INTERFERON_GAMMA_RESPONSE.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(feature_data[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
  names(go_ifng_response_genes)=unique(feature_data[,.(Symbol,ID)])[,setNames(ID,Symbol)][go_ifng_response_genes]
  go_ifna_response_genes=fread(sprintf("%s/GO_0034340_INTERFERON_ALPHA_RESPONSE.txt",REF_DIR),col.names="V1")[,setNames(V1,V1)][unique(feature_data[,.(Symbol,ID)])[,setNames(Symbol,ID)]]%>%{.[!is.na(.)]}
  names(go_ifna_response_genes)=unique(feature_data[,.(Symbol,ID)])[,setNames(ID,Symbol)][go_ifna_response_genes]
  go_ifn_response_genes=unique(c(go_ifng_response_genes,go_ifna_response_genes))
  names(go_ifn_response_genes)=unique(feature_data[,.(Symbol,ID)])[,setNames(ID,Symbol)][go_ifn_response_genes]
  type_I_IFNs=c(feature_data$Symbol[grepl("^IFN[ABKZ][0-9]",feature_data$Symbol)],"IFNE")%>%{setNames(feature_data$ID,feature_data$Symbol)[.]} %>% {.[!is.na(.)]}
  type_II_IFNs=c("IFNG") %>% {setNames(feature_data$ID,feature_data$Symbol)[.]} %>% {.[!is.na(.)]}
  IFNA1=feature_data[grepl("IFNA[0-9]",Symbol),setNames(ID,Symbol)]

  # Create Seurat object to run AddModuleScores on
  seurat=AddModuleScore(seurat,
  features=list(
  hallmark_ifna_response_score=names(hallmark_ifna_response_genes),
  hallmark_ifng_response_score=names(hallmark_ifng_response_genes),
  hallmark_isg_score=names(hallmark_ifn_exclusive),
  hallmark_inflammatory_score=names(hallmark_inflammatory_exclusive),
  go_ifna_response_score=names(go_ifna_response_genes),
  go_ifng_response_score=names(go_ifng_response_genes),
  go_isg_score=names(go_ifn_response_genes),
  type_1_ifn_score=type_I_IFNs,
  type_2_ifn_score=type_II_IFNs,
  ifna_score=IFNA1
  ),name=c("hallmark_ifna_response_score","hallmark_ifng_response_score","hallmark_isg_score","hallmark_inflammatory_score","go_ifna_response_score","go_ifng_response_score","go_isg_score","type_1_ifn_score","type_2_ifn_score","ifna_score"))
  colnames(seurat@meta.data) %<>% str_replace("[0-9]+$","")

  as.data.table(seurat@meta.data)[,-c("orig.ident","nCount_RNA","nFeature_RNA")]%>%
  fwrite(.,"/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/1_dataset_description/scores__by__IID_condition.tsv",sep="\t")
  # scores_IID_cond=fread("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/1_dataset_description/scores__by__IID_condition.tsv")

################# load Cell Mortality ###########################

  CellMortality=fread(sprintf('%s/single_cell/project/pop_eQTL/Cell_count_library_mortality.txt',EVO_IMMUNO_POP_ZEUS))
  setnames(CellMortality,c('MEAN','MORTALITY'),c('CellCount','Mortality'))
  CellMortality=CellMortality[!duplicated(IID),.(IID,CellCount,Mortality)]
  mod=CellMortality[,lm(Mortality~CellCount)]
  CellMortality[is.na(Mortality),Mortality:=round(predict(mod,newdata=data.frame(CellCount=CellCount)))]

################################################################################
# Fig. 1D

#fig1d_data=fread(sprintf("%s/Fig1D_data__ISG_score_byCOND.tsv",FIGURE_DIR))
fig1d_data=fread(sprintf("%s/1_dataset_description/scores__by__IID_condition.tsv",DATA_DIR))
fig1d_data[,condition:=factor(condition,cond_order)]

fig1d_data=merge(fig1d_data,CellMortality,by=c("IID"))
fig1d_data[,POP:=substr(IID,1,3)]

adjust_on_x=function(y,matx,matz){
  # y : la variable à ajuster (un vecteur)
  # mat_x: un vecteur une matrix ou un data.frame contenant les variables sur lesquelles on veut ajuster : (eg les cell types)
  # mat_z:  un vecteur une matrix ou un data.frame contenant les variables à inclure dans le modèle sur lesquelles on ne veut pas ajuster : (eg la population)
  if(is.data.frame(matx)){matx=as.matrix(matx)}
  if(is.vector(matx)){matx=matrix(matx,ncol=1)}
  if(is.data.frame(matz)){matz=as.matrix(matz)}
  if(is.vector(matz)){matz=matrix(matz,ncol=1)}

    y_z=lm(y~.,data=as.data.frame(matz))$res
    matx_z=apply(matx,2,function(x){lm(x~.,data=as.data.frame(matz))$res})
    beta_x_z=lm(y_z~.,data=as.data.frame(matx_z))$coef[-1]
    y_adj=y-matx%*%matrix(beta_x_z,ncol=1)
    as.vector(y_adj)
}
fig1d_data[,hallmark_isg_score_adj:=adjust_on_x(hallmark_isg_score,Mortality,model.matrix(~1+POP)[,-1])]

rankTransform = function(x){
    percentile=rank(x,ties.method='random',na.last = NA)/(length(x)+1)
    mean_level=mean(x,na.rm=TRUE)
    sd_level=sd(x,na.rm=TRUE)
    qnorm(percentile,mean_level,sd_level)
    }

# increased variance in response to SARS-COV-2 is driven by a handful of outliers...
bartlett.test(hallmark_isg_score_adj~condition,data=fig1d_data[IID!="ASH071" & condition!='NS',]) # p< 1E-9

bartlett.test(rankTransform(hallmark_isg_score_adj)~condition,data=fig1d_data[IID!="ASH071" & condition!='NS',]) # p = 0.33
fligner.test(hallmark_isg_score_adj~condition,data=fig1d_data[IID!="ASH071" & condition!='NS',]) #  p=0.03


#fig1d_data[,hallmark_isg_score_adj:=lm(hallmark_isg_score~Mortality)$res+mean(hallmark_isg_score),by=condition]

fig1d_plot=ggplot(fig1d_data[IID!="ASH071",]) + theme_yann() + ylim(c(0,3)) +
  geom_violin(aes(x=condition,y=hallmark_isg_score_adj,fill=condition,col=condition),alpha=0.5,scale='width') +
  geom_boxplot(aes(x=condition,y=hallmark_isg_score_adj,col=condition),notch=TRUE,fill='white',alpha=0.5) +
  scale_fill_manual(values=cond_color)+
  scale_color_manual(guide="none",values=color_conditions)+
  theme(text=element_text(size=7))+ylab('Mean ISG activity per individual') + xlab('')

  saveRDS(fig1d_plot,file=sprintf('%s/Fig1/fig1d_plot.RDS',FIG_DIR))
  pname=sprintf("%s/Fig1/Fig1d.pdf",FIG_DIR)

  pdf(pname,width=2,height=4)
    grid.arrange(
      grobs=list(
        ggplotGrob(fig1d_plot+theme(legend.position="none",text=element_text(size=10))),
  #      get_legend(fig1b_plot)#,
        grid.rect(gp=gpar(col="white"))
      ),
      layout_matrix=rbind(
        c(1,1),
        c(2,2)
      ), heights=c(7,3)
    )
  dev.off()


################################################################################
# Fig. 1E
SIMOA=fread(sprintf('%s/single_cell/project/pop_eQTL/Covariates/Simoa_IFNs_clean.txt',EVO_IMMUNO_POP_ZEUS))
SIMOA[,logIFNa17:=log10(IFNa17)]
SIMOA[,logIFNb:=log10(IFNb*1000)]
SIMOA[,logIFNg:=log10(IFNg)]

IID_metadata=fread(sprintf('%s/popCell_data/00_CRF/scrnaseq_popbased_metadata_full_long.tsv',EVO_IMMUNO_POP_ZEUS))
IID_metadata=unique(IID_metadata[,.(IID,Age,Gender)])

IID_metadata=merge(IID_metadata,CellMortality,by='IID')
IID_metadata[,POP:=substr(IID,1,3)]
SIMOA=merge(SIMOA,IID_metadata,by='IID')

ISG_scores=fread(sprintf('%s/1_dataset_description/scores__by__IID_lineage_condition.tsv',DATA_DIR))
setnames(ISG_scores,"condition","COND")
SIMOA=merge(ISG_scores,SIMOA,by=c('IID','COND'))

# load cell cell propoortions and covariates (age, sex, pop, mortality)

cellProps=fread(sprintf('%s/2_population_differences/Cellcounts/Pct_celltype_by_IID_condition.tsv.gz',DATA_DIR))
current_celltypes=c( "B.M.K","B.M.L", "B.N.K","B.N.L", "ILC","MAIT","MONO.CD14","MONO.CD16", "NK.CD56brt","NK.CD56dim","NK.M.LIKE","Plasmablast","T.CD4.E","T.CD4.N", "T.CD8.CM.EM","T.CD8.EMRA","T.CD8.N","T.Reg","T.gd","cDC","pDC","MONO.CD14.INFECTED")

cellProps[,POP:=substr(IID,1,3)]
setnames(cellProps,'state','COND')

CellProp_and_Covariates=list()
CellProp_and_Covariates[['NS']]=fread(sprintf('%s/single_cell/project/pop_eQTL/data/2_population_differences/Covariates/lineage_condition__CellPropCelltype/Covariates__T.CD4_NS.tsv.gz',EVO_IMMUNO_POP_ZEUS))
CellProp_and_Covariates[['NS']][,MONO.CD14.INFECTED:=0]
CellProp_and_Covariates[['COV']]=fread(sprintf('%s/single_cell/project/pop_eQTL/data/2_population_differences/Covariates/lineage_condition__CellPropCelltype/Covariates__T.CD4_COV.tsv.gz',EVO_IMMUNO_POP_ZEUS))
CellProp_and_Covariates[['COV']][,MONO.CD14.INFECTED:=0]
CellProp_and_Covariates[['IAV']]=fread(sprintf('%s/single_cell/project/pop_eQTL/data/2_population_differences/Covariates/lineage_condition__CellPropCelltype/Covariates__T.CD4_IAV.tsv.gz',EVO_IMMUNO_POP_ZEUS))
CellProp_and_Covariates=rbindlist(CellProp_and_Covariates,idcol='COND', use.names=TRUE)

# add to SIMOA data
SIMOA=merge(CellProp_and_Covariates[,-c("Age","Mortality")],SIMOA,by=c('IID','COND'))

test_CAR=function(DT,responseVar,covariates,lambda=NULL,covariates_noPenalty=NULL,remove){
    require(care)
      data_to_analyse=DT[,mget(c(covariates,covariates_noPenalty))]
      complete_rows=apply(!is.na(data_to_analyse),1,all)
      Y_var=DT[complete_rows==TRUE,get(responseVar)]
      Y_var=matrix(Y_var,ncol=length(responseVar))
  #   Y_var=apply(Y_var,2,rank_transform)
      temp = model.matrix(~0+.,data=data_to_analyse)
      if(is.null(lambda)){
        CAR <- carscore(as.matrix(temp), Y_var)
      }else{
        CAR <- carscore(as.matrix(temp), Y_var,lambda=lambda)
      }
      list(cov=names(CAR),CAR=CAR)
    }

cov_noPenalty=c("POP","Age","Gender","Mortality")

Explained_VAR=SIMOA[celltype%in%c('MONO','T.CD4','T.CD8','NK','B'),test_CAR(.SD,"hallmark_isg_score",c('logIFNa17','logIFNg','logIFNb'),covariates_noPenalty=cov_noPenalty),by=.(celltype,COND)]
Explained_VAR[,is_IFN:=grepl('logIFN',cov)]
Explained_VAR[,VAR:=CAR^2]
fig1e_data=Explained_VAR[is_IFN==TRUE,]

fig1e_data$cov=str_replace(fig1e_data$cov,"^log","")

fwrite(fig1e_data,"/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/paper_draft/V2/figureMaterial/Fig1/Fig1E__data_SiMOA.tsv",sep="\t")

ifn_color=c('IFNa17'='#227c9d','IFNg'='#ffcb77','IFNb'='#17c3b2')
ifn_order=c('IFNa17','IFNb','IFNg')
lineage_order=c("MONO","B","T.CD4","T.CD8","NK")

fig1e_data$celltype=factor(fig1e_data$celltype,lineage_order)
fig1e_data$cov=factor(fig1e_data$cov,rev(ifn_order))
fig1e_data$COND=factor(fig1e_data$COND,cond_order)

fig1e_plot=ggplot(fig1e_data,aes(x=celltype,y=100*VAR,fill=cov))+
geom_bar(stat='Identity')+
scale_fill_manual(values=ifn_color)+
facet_grid(~factor(COND,c('NS','COV','IAV')))+theme_yann()+
scale_x_discrete(breaks=lineage_order,labels=c("MONO","B","CD4+ T","CD8+ T","NK"))+
ylab('Percentage of ISG variance explained')+xlab("Cell types")+
guides(fill=guide_legend(override.aes=list(pch=15,alpha=1,size=0.5),ncol=1,byrow=TRUE)) +
theme(panel.spacing=unit(0,"mm"),text=element_text(size=7),axis.text.x=element_text(angle=45,hjust=1))

saveRDS(fig1e_plot,file=sprintf('%s/Fig1/fig1e_plot.RDS',FIG_DIR))
pname=sprintf("%s/Fig1/Fig1e.pdf",FIG_DIR)

pdf(pname,width=4,height=5)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig1e_plot+theme(legend.position="none",text=element_text(size=10))),
#      get_legend(fig1b_plot)#,
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(1,1),
      c(2,2)
    ), heights=c(7,3)
  )
dev.off()

# pdf(sprintf('%s/SIMOA_ISG_score_varexplained_new.pdf',FIGURE_DIR),width=12)
# print(p)
# dev.off()

################################################################################
# Fig. 1F
EXPRESSION_FILE="BatchAdjusted_logCPM_125libs__per_celltype_condition_annotated.tsv.gz"
Expression=fread(sprintf("%s/single_cell/project/pop_eQTL/data/2_population_differences/%s",EVO_IMMUNO_POP_ZEUS,EXPRESSION_FILE))[IID!='ASH071',]

feature_data=Expression[,.(ID,Symbol)]%>%unique

IFNA1=feature_data[grepl("IFNA[0-9]",Symbol),setNames(ID,Symbol)]
type_I_IFNs=c(feature_data$Symbol[grepl("^IFN[ABKZ][0-9]",feature_data$Symbol)],"IFNE")%>%{setNames(feature_data$ID,feature_data$Symbol)[.]} %>% {.[!is.na(.)]}

# compute CPMs from logCPMs for IFN-a genes
Expression_IFNA=Expression[ID%in%IFNA1,.(ID,Symbol,celltype,IID,state,CPM=2^(logCPM)-1)]
# sum across all IFNA genes
Expression_IFNA=Expression_IFNA[,.(CPM=sum(CPM)),by=.(IID,celltype,condition=state)]
# count number of UMI in million for each cell type in each sample
celltype_condition_IID_counts=meta_clean[,.(N=.N,total_UMI_in_million=sum(sum)/1e6),by=.(IID,celltype,condition)]
Expression_IFNA=merge(Expression_IFNA,celltype_condition_IID_counts,by=c('IID',"celltype","condition"))

# compute number of counts in each sample and cell type
Expression_IFNA[,counts:=CPM*total_UMI_in_million]

Expression_IFNA[celltype=='pDC',summary(lm(log2(CPM+1)~factor(condition,c('NS','COV','IAV'))))]
# Call:
# lm(formula = log2(CPM + 1) ~ factor(condition, c("NS", "COV", "IAV")))
# Coefficients:
#                                             Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                   0.2246     0.3001   0.748    0.455
# factor(condition, c("NS", "COV", "IAV"))COV   6.3652     0.4457  14.282   <2e-16 ***
# factor(condition, c("NS", "COV", "IAV"))IAV  12.1698     0.4174  29.158   <2e-16 ***

Expression_IFNA[celltype=='pDC',summary(lm(log2(CPM+1)~factor(condition,c('NS','COV','IAV'))))]
FoldChanges_pDC=mutate(Expression_IFNA[celltype=='pDC',],logCPM=log2(1+pmax(CPM,0))) %>% dcast(IID~condition,value.var='logCPM') %>% mutate(logFC_COV=COV-NS,logFC_IAV=IAV-NS)
mean(FoldChanges_pDC$logFC_COV,na.rm=T)
mean(FoldChanges_pDC$logFC_IAV,na.rm=T)
wilcox.test(FoldChanges_pDC$logFC_COV,FoldChanges_pDC$logFC_IAV)
# average across individuals

Expression_IFNA_avg=Expression_IFNA[,.(meanIFNAcount_perInd=mean(counts),meanIFNA_perCell=mean(CPM)),by=.(celltype,condition)]
Expression_IFNA_avg[,.(round(meanIFNAcount_perInd,1),round(meanIFNA_perCell,1)),by=.(celltype,condition)]

fig1f_data=Expression_IFNA_avg
fwrite(fig1f_data,"/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/paper_draft/V2/figureMaterial/Fig1/Fig1F__data_IFNA.tsv",sep="\t")

fig1f_plot=ggplot(Expression_IFNA_avg[condition!="NS"])+
  theme_yann(rotate.x=90)+
  geom_segment(aes(x=0, xend=meanIFNAcount_perInd,y=celltype,yend=celltype,color=celltype))+
  scale_x_continuous(trans='sqrt',limits=c(-1,4500), breaks=c(10,250,1000,2500),labels=c('10','250','1000','2500'))+
#  scale_x_continuous(trans='sqrt') +
  geom_point(aes(x=meanIFNAcount_perInd,y=celltype,fill=celltype,size=meanIFNA_perCell),pch=21,alpha=0.5)+
  geom_point(aes(x=meanIFNAcount_perInd,y=celltype,color=celltype,size=meanIFNA_perCell),pch=10,fill=NA)+
  scale_color_manual(aesthetics=c("color","fill"),guide="none",values=celltype_color) +
  scale_size_continuous(breaks=c(50,100,150))+
  xlab(expression("Mean number of IFN-"*alpha*" transcript per sample")) + ylab("Cell type") +
  #facet_grid(cols=vars(COND_ifn),labeller=labeller(COND_ifn=COND_ifn_labels))+guides(size=guide_legend(name=""))+
  facet_grid(cols=vars(condition))+
  guides(size=guide_legend(name=""))+
  theme(text=element_text(size=7),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing=unit(0,"pt"),axis.text.x=element_text(angle=45,vjust=1))

saveRDS(fig1f_plot,file=sprintf('%s/Fig1/fig1f_plot.RDS',FIG_DIR))
pname=sprintf("%s/Fig1/Fig1f.pdf",FIG_DIR)

pdf(pname,width=4,height=5)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig1f_plot+theme(legend.position="none",text=element_text(size=10))),
#      get_legend(fig1b_plot)#,
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(1,1),
      c(2,2)
    ), heights=c(7,3)
  )
dev.off()

################################################################################
# Fig. 1G

ISG_scores=fread(sprintf('%s/1_dataset_description/scores__by__IID_condition.tsv',DATA_DIR))
ISG_scores[,POP:=substr(IID,1,3)]

ISG_scores=merge(ISG_scores,CellMortality,by=c("IID"))
#ISG_scores[,hallmark_isg_score_adj:=lm(hallmark_isg_score~Mortality)$res+mean(hallmark_isg_score),by=condition]
ISG_scores[,hallmark_isg_score_adj:=adjust_on_x(hallmark_isg_score,Mortality,model.matrix(~1+POP)[,-1])]

fig1g_data=dcast(ISG_scores,IID+POP~condition,value.var="hallmark_isg_score_adj")
fig1g_data[,.(COR=cor(COV,IAV),P=cor.test(COV,IAV)$p.value)]
# COR            P
# 1: 0.5953208 1.115494e-22

#fig1g_data=dcast(ISG_scores,IID+POP~condition,value.var="hallmark_isg_score")

# fig1g_data=fread(sprintf("%s/Fig1F_data__ISG_score_vs_Celltype.tsv",FIGURE_DIR))
# fig1g_data=melt(fig1f_data,id.vars=c("IID", "B.M", "B.N", "MAIT", "Myeloid.CD14", "Myeloid.CD16","NK.CD56brt", "NK.CD56dim", "Plasmablast", "Progenitors", "T.CD4.E","T.CD4.N", "T.CD8.CM.EM", "T.CD8.N", "T.CD8.TE", "T.Reg", "T.gd","cDC", "mDC", "pDC", "Pct_IAV", "Myeloid.CD14_nInf", "CellCount","Mortality", "POP",'Gender','Age'))
# CellPct_and_scores_ind_withQC[,COND:=gsub('^(.*)_(IAV|COV|NS)$','\\2',variable)]
# CellPct_and_scores_ind_withQC[,variable:=gsub('^(.*)_(IAV|COV|NS)$','\\1',variable)]
#
# CellPct_and_scores_ind_cond_withQC=merge(CellPct_and_scores,scores_pbmc,by='IID')

#
# CellPct_and_scores_select=CellPct_and_scores_ind_withQC[COND!='NS',.(IID, ISG_score, COND, T.CD8.N, Myeloid.CD14, pDC, Mortality, POP)]
# CellPct_and_scores_select=melt(CellPct_and_scores_select,id.vars= c('IID','ISG_score','POP','Mortality','COND'))
# CellPct_and_scores_select[,ISG_score_adj:=lm(ISG_score~Mortality+POP)$res,by=.(COND,variable)]
# CellPct_and_scores_select[,Percentage_adj:=lm(value~Mortality+POP)$res,by=.(COND,variable)]

fig1g_plot=ggplot(fig1g_data,aes(x=IAV,y=COV,col=POP,fill=POP)) + theme_yann() + xlim(c(0,3)) + ylim(c(0,3)) + geom_abline(Intercept=0,slope=1,col='lightgrey',linetype=2) +
      geom_point(alpha=.5) + scale_color_manual(values=color_populations,guide="none") +
      scale_fill_manual(values=color_populations,breaks=c("AFB","EUB","ASH"),name="") +
      geom_smooth(method='lm',col='black',fill="grey",show.legend=FALSE) + ylab('ISG activity (COV)') + xlab("ISG activity (IAV)") +
      guides(fill=guide_legend(override.aes=list(pch=21,alpha=1,size=2,color="black"),ncol=1,byrow=TRUE)) +
      theme(text=element_text(size=7))
saveRDS(fig1g_plot,file=sprintf('%s/Fig1/fig1g_plot.RDS',FIG_DIR))
pname=sprintf("%s/Fig1/fig1g.pdf",FIG_DIR)

pdf(pname,width=4,height=5)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig1g_plot+theme(legend.position="none",text=element_text(size=10))),
#      get_legend(fig1b_plot)#,
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(1,1),
      c(2,2)
    ), heights=c(7,3)
  )
dev.off()


#### check IFNa adjusted correlation
SIMOA=fread(sprintf('%s/single_cell/project/pop_eQTL/Covariates/Simoa_IFNs_clean.txt',EVO_IMMUNO_POP_ZEUS))
SIMOA[,logIFNa17:=log10(IFNa17)]
SIMOA[,logIFNb:=log10(IFNb*1000)]
SIMOA[,logIFNg:=log10(IFNg)]
setnames(SIMOA,'COND','condition',skip_absent=TRUE)
ISG_scores_SIMOA=merge(ISG_scores,SIMOA,by=c('IID','condition'))
ISG_scores_SIMOA[,hallmark_isg_score_IFNadj:=adjust_on_x(hallmark_isg_score,cbind(Mortality,logIFNa17),model.matrix(~1+POP)[,-1])]

ISG_scores_IFNadj=dcast(ISG_scores_SIMOA,IID+POP~condition,value.var="hallmark_isg_score_IFNadj")

ISG_scores_IFNadj[,.(COR=cor(COV,IAV),P=cor.test(COV,IAV)$p.value)]
# COR            P
# 1: 0.5325759 4.477188e-17


################################################################################

pname=sprintf("%s/Fig1_final_legend.pdf",FIGURE_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      grid.rect(gp=gpar(col="white")),
      get_legend(fig1b_plot)
    ),
    layout_matrix=rbind(
      c(1,1,1,1,1,1,2),
      c(1,1,1,1,1,1,2),
      c(1,1,1,1,1,1,2),
      c(1,1,1,1,1,1,2)
    ), heights=c(4,2,3,1), widths=c(2,2,1,1,1,1,2)
  )
dev.off()



# legend.position="none")),
#       ggplotGrob(fig1d_plot+theme(legend.position="none")),
#       ggplotGrob(fig1e_plot+theme(legend.position="none",axis.title.x=element_blank())),
#       ggplotGrob(fig1f_plot+theme(legend.position="none")),
#       ggplotGrob(fig1g_plot+theme(legend.position=c(0.8,0.3)))
#     ),
#     layout_matrix=rbind(
#       c(1,1,2,2,2,2,1),
#       c(1,1,3,3,3,3,1),
#       c(4,5,5,6,6,7,7),
#       c(1,1,1,1,1,1,1)
#     ), heights=c(4,2,3,1), widths=c(2,2,1,1,1,1,2)
#   )
# dev.off()
# Fig1/Fig1.R [+]
pname=sprintf("%s/Fig1_final.pdf",FIGURE_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      grid.rect(gp=gpar(col="white")),
      ggplotGrob(fig1b_plot+theme(legend.position="none")),
      ggplotGrob(fig1c_plot+theme(legend.position="none")),
      ggplotGrob(fig1d_plot+theme(legend.position="none")),
      ggplotGrob(fig1e_plot+theme(legend.position="none",axis.title.x=element_blank())),
      ggplotGrob(fig1f_plot+theme(legend.position="none")),
      ggplotGrob(fig1g_plot+theme(legend.position=c(0.2,0.7)))
    ),
    layout_matrix=rbind(
      c(1,1,2,2,2,2,1),
      c(1,1,3,3,3,3,1),
      c(4,5,5,6,6,7,7),
      c(1,1,1,1,1,1,1)
    ), heights=c(4,2,3,1), widths=c(2,2,1,1,1,1,2)
  )
dev.off()


pname=sprintf("%s/Fig1_final.pdf",FIGURE_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      grid.rect(gp=gpar(col="white")),
      grid.rect(gp=gpar(col="white")),
      grid.rect(gp=gpar(col="white")),
      # ggplotGrob(fig1b_plot+theme(legend.position="none")),
      # ggplotGrob(fig1c_plot+theme(legend.position="none")),
      ggplotGrob(fig1d_plot+theme(legend.position="none")),
      ggplotGrob(fig1e_plot+theme(legend.position=c(0.15,0.8),legend.spacing.y=unit(-1,'mm'),axis.title.x=element_blank(),legend.key.size=unit(4,"mm"))),
      ggplotGrob(fig1f_plot+theme(legend.position="none")),
      ggplotGrob(fig1g_plot+theme(legend.position=c(0.15,0.75),legend.spacing.y=unit(-2,'mm')))
    ),
    layout_matrix=rbind(
      c(1,1,2,2,2,2,1),
      c(1,1,3,3,3,3,1),
      c(4,5,5,6,6,7,7),
      c(1,1,1,1,1,1,1)
    ), heights=c(4,2,3,1), widths=c(2,2,1,1,1,1,2)
  )
dev.off()
