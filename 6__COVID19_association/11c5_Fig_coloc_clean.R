
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
 eQTL_IRF1=get_eQTL('rs10774671','OAS1',resolution='celltype',metric='logCPM')

 nIQR=3
  eQTL_IRF1=eQTL_IRF1[Symbol=='OAS1',]
  eQTL_IRF1[,state:=factor(state,c('NS','COV','IAV'))]
  nIQR=3
  eQTL_IRF1[,logCPM_IQR:=pmax(pmin(logCPM,median(logCPM)+nIQR*IQR(logCPM)),median(logCPM)-nIQR*IQR(logCPM)),by=.(state,celltype)]
  eQTL_IRF1[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]
  irnt=function(x){qnorm(rank(x)/(length(x)+1),mean(x),sd(x))}
  eQTL_IRF1[,logCPM_irnt:=irnt(logCPM),by=.(state,celltype)]


p <- ggplot(eQTL_IRF1[Symbol=='OAS1' & celltype%in%c('T.CD8.EMRA','T.CD8.N','NK.CD56dim','NK.CD56dim','T.gd'),],aes(x=Number_of_ALT_alelle,y=logCPM_IQR,fill=state))
p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=TRUE)
p <- p + theme_yann() + scale_fill_manual(values=color_conditions)# + scale_color_manual(values=color_populations)
p <- p + facet_grid(state~celltype) +ylab('logCPM')
p <- p+theme_yann()

  ####------------------------------------------------------------------------------------####
  ####------------------------------------------------------------------------------------####
  FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V11/newFigure_6/',EIP)
  pdf(sprintf('%s/Fig6C_eQTL_IRF1.pdf',FIGURE_DIR),width=7.2*1,height=6.7*1)
  print(p)
  dev.off()

  Fig3F_data=eQTL_MMP1
  fwrite(Fig3F_data,file=sprintf('%s/Final/Fig3F_data_eQTL_MMP1_MONO.tsv',FIGURE_DIR))
  saveRDS(Fig3F,file=sprintf('%s/Final/Fig3F_data_eQTL_MMP1_MONO.RDS',FIGURE_DIR))
