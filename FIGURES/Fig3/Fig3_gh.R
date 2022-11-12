################################################################################
# header

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"
FIGURE="4"
DAT_DIR=sprintf("%s/single_cell/project/pop_eQTL/data",EIP)
PAP_DIR=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7",EIP)
FIG_DIR=sprintf("%s/testsFigure",PAP_DIR)
FIGURE_DIR=sprintf("%s/FigureMaterial/Figure%s",PAP_DIR,FIGURE)
POP_DIR=sprintf("%s/2_population_differences/popDE",DAT_DIR)

source(sprintf("%s/template_scripts/processing_pipeline/00_set_colors.R",RES_DIR))
source(sprintf("%s/template_scripts/querySNPs.R",RES_DIR))
source(sprintf("%s/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",RES_DIR))

celltype_color=c(
"T.CD4.N"="#169FD8","T.CD4.E"="#005274","T.Reg"="#03C2C0","T.gd"="#ffce73",
"T.CD8.N"="#00BB54","T.CD8.CM.EM"="#0EAD20","T.CD8.EMRA"="#048c41","MAIT"="#004A21",
"NK.M.LIKE"="#63C425","NK.CD56dim"="#B7DC2A","NK.CD56brt"="#7D9E00","ILC"="#3F4F00",
"B.N.K"="#9281DE","B.N.L"="#a681de","B.M.K"="#6045d9","B.M.L"="#8045d9",
"Plasmablast"="#6801A4","MONO.CD14"="#B7637E","MONO.CD16"="#8F568A","pDC"="#8E6B00",
"cDC"="#E7b315","MONO.CD14.INFECTED"="#B9364B")

lineage_order=c("MONO","B","T.CD4","T.CD8","NK")

lineage_label=c("MYELOID","B","CD4+ T","CD8+ T","NK")

lineage_celltype=fread(sprintf("%s/1_dataset_description/lineages_celltype.tsv",DAT_DIR))

################################################################################
ANALYSE="lineage_condition_AFBEUB__lineage_condition__220409"
popDE_noCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDE_noCellProp[,type:="expression"]
ANALYSE="lineage_condition_AFBEUB__logFC_lineage_condition_logFC__220409"
popDR_noCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDR_noCellProp[,type:="response"]

# popDEGs tested within each lineage, accounting for finer cell proportions within the lineage considered
# (NB: MASH run on all lineages together)
ANALYSE="lineage_condition_AFBEUB__lineage_condition__CellPropLineage_220409"
popDE_lineageCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDE_lineageCellProp[,type:="expression"]
ANALYSE="lineage_condition_AFBEUB__logFC_lineage_condition_logFC__CellPropLineage_220409"
popDR_lineageCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDR_lineageCellProp[,type:="response"]

################################################################################
# Fig. 4a

eqtls<-eQTL_Signif_both[celltype%in%lineages,]
setnames(eqtls,"gene_name","Symbol")

reqtls<-reQTL_Signif_both[celltype%in%lineages,]
setnames(reqtls,"gene_name","Symbol")

popdiff_losses=fread(sprintf("%s/2_population_differences/cell_composition_adjustment_beta_changes.tsv",DAT_DIR))

# fraction of all genes with an eQTL
eqtls_genome_wide=eqtls[,length(unique(Symbol))/12672*100,by=celltype]
setnames(eqtls_genome_wide,"V1","frac_eqtl")
eqtls_genome_wide$adj="gw"
eqtls_genome_wide$type="DE"

# fraction of all genes with an reQTL
reqtls_genome_wide=reqtls[,length(unique(Symbol))/12672*100,by=celltype]
setnames(reqtls_genome_wide,"V1","frac_eqtl")
reqtls_genome_wide$adj="gw"
reqtls_genome_wide$type="DR"

# fraction of popDEGs (raw) with an eQTL
fig4a_data=lapply(c("DE","DR"),function(z){
  lapply(c("no","lineage"),function(y){
    popdegs=get(sprintf("pop%s_%sCellProp",z,y))
    lapply(lineage_order,function(x){
      if (y=="lineage"){
      pdegs_in_lineage=popdegs[FDR<0.01&abs(beta)>0.2&celltype==x,unique(Symbol)]
      signif_losses=popdiff_losses[type==z&celltype==x&adjustment_phase=="no_lineage"&difference_signif=="Different"&difference_type=="difference_lost",unique(Symbol)]
      pdegs_in_lineage=pdegs_in_lineage[pdegs_in_lineage%nin%signif_losses]
      } else {
      pdegs_in_lineage=popdegs[FDR<0.01&abs(beta)>0.2&celltype==x,unique(Symbol)]
      }
      pdegs_in_lineage_w_eqtl=eqtls[celltype==x&Symbol%in%pdegs_in_lineage,length(unique(Symbol))]
      n=pdegs_in_lineage_w_eqtl/length(pdegs_in_lineage)*100
      data.table(celltype=x,frac_eqtl=n,adj=y,type=z)
    })%>%rbindlist()
  })%>%rbindlist()
})%>%rbindlist()
fig4a_data=rbind(fig4a_data,eqtls_genome_wide,reqtls_genome_wide)
fig4a_data$celltype=factor(fig4a_data$celltype,rev(lineage_order))
fig4a_data$adj=factor(fig4a_data$adj,rev(c("gw","no","lineage")))

fwrite(fig4a_data,sprintf("%s/Fig4a_data.tsv",FIGURE_DIR),sep="\t")

fig4a_plot=ggplot(fig4a_data[type=="DE"],aes(frac_eqtl,celltype,fill=celltype,alpha=adj,group=adj))+
  geom_vline(data=fig4a_data[type=="DE",mean(frac_eqtl),by=adj],mapping=aes(xintercept=V1,alpha=adj),linetype="dashed",size=0.2)+
  geom_col(position="dodge")+
  scale_x_continuous(limits=c(0,60),breaks=c(0,15,30,45,60),position="bottom")+
  scale_fill_manual(values=lineage_color,guide="none")+
  scale_y_discrete(labels=rev(c("MONO","B","CD4+ T","CD8+ T","NK")),position="left")+
  scale_alpha_manual(values=c("no"=0.6,"lineage"=1,"gw"=0.2),breaks=c("gw","no","lineage"),labels=c("Genome-wide","Raw","Intra-lineage"),
		     guide=guide_legend(override.aes=list(shape=21,color="black")))+
  xlab("Genes with an eQTL (%)")+
  theme_yann()+
  theme(text=element_text(size=5),axis.title.y=element_text(color="white"))+
  theme(panel.spacing=unit(0,"pt"),legend.position="bottom")

p=ggplot(fig4a_data,aes(celltype,frac_eqtl,fill=celltype,alpha=adj,group=adj))+
  geom_point()+
  scale_fill_discrete(guide="none")+
  scale_alpha_manual(values=c("no"=0.6,"lineage"=1,"gw"=0.2),breaks=c("gw","no","lineage"),
		     labels=c("Genome-wide        ","Raw             ","Adjusted"),guide=guide_legend(override.aes=list(shape=21,color=NA,fill="black",size=2),nrow=1,title.position="top",title.hjust=1),name="Cell proportion adjustment")+
  theme_yann()+
  theme(legend.position="bottom",text=element_text(size=5),legend.title=element_text())
fig4a_legend=get_legend(p+theme(text=element_text(size=5),legend.spacing.x=unit(-1,'mm')))

################################################################################
# Fig. 4b

run_id="220617"
popdiff=T
exp_resp="DE"
comparison="EUB_AFB"
MED_DIR=sprintf("%s/2_population_differences/mediation_analysis",DAT_DIR)

#if(exp_resp=="DE"){
#  popdegs<-popDE_noCellProp
#}else{
#  popdegs<-popDR_noCellProp
#}
#popdegs<-popdegs[FDR<0.01&abs(beta)>0.2,]
#
#bin_order=c("VL","L","M","H","VH")
#breaks=popdegs[,quantile(abs(beta),c(0,0.2,0.4,0.6,0.8,1))]
#popdegs[,bin:=cut(abs(beta),breaks,include.lowest=T,right=F,labels=bin_order)]
#
#if(exp_resp=="DE"){
#  eQTL_file=fread(sprintf("%s/2_population_differences/best_eQTL_per_gene_and_comparison_MAF005.tsv",DAT_DIR))
#}else{
#  eQTL_file=fread(sprintf("%s/2_population_differences/best_reQTL_per_gene_and_comparison_MAF005.tsv",DAT_DIR))
#}
#
#mediation_genetics<-fread(sprintf("%s/genetics/mediation_analysis_genetics_%s%s%s_%s.tsv",MED_DIR,comparison,ifelse(popdiff==T,"_popD","_"),gsub("^D","",exp_resp),run_id))
#mediation_genetics[,mediator_type:='genetics']
#
#mediation_celltype<-fread(sprintf("%s/celltype_frequencies/mediation_analysis_celltype_%s%s%s_%s.tsv",MED_DIR,comparison,ifelse(popdiff==T,"_popD","_"),gsub("^D","",exp_resp),run_id),sep="\t")
#celltype_mediators<-c("MONO.CD16","B.M.K","T.CD4.E","T.CD8.EMRA","NK.M.LIKE")
#mediation_celltype_best=mediation_celltype[mediator%in%celltype_mediators,]
#mediation_celltype_best[,mediator_type:='celltype']
#mediation_celltype_best[,type:=exp_resp]
#
#mediation<-rbind(mediation_genetics,mediation_celltype_best[,mget(colnames(mediation_genetics))])
#mediation[,celltype:=factor(celltype,lineage_order)]
#mediation[,state:=factor(state,cond_order)]
#
#mediation$betapop=popdegs[,setNames(beta,sprintf("%s%s%s",celltype,state,Symbol))][mediation[,sprintf("%s%s%s",celltype,state,Symbol)]]
#mediation$bin=popdegs[,setNames(bin,sprintf("%s%s%s",celltype,state,Symbol))][mediation[,sprintf("%s%s%s",celltype,state,Symbol)]]
#
#mediation$betaeqtl=eQTL_file[comparison=="EUB_AFB",setNames(beta,sprintf("%s%s%s",celltype,state,snps))][mediation[,sprintf("%s%s%s",celltype,state,mediator)]]
#mediation$peqtl=eQTL_file[comparison=="EUB_AFB",setNames(pvalue,sprintf("%s%s%s",celltype,state,snps))][mediation[,sprintf("%s%s%s",celltype,state,mediator)]]
#
#if(exp_resp=="DR"){
#  eqtl_signif_info=merge(eQTL_Signif_both,SNP_info[,.(snps=ID,AFB,EUB,ASH)],by="snps",all.x=T)
#}else{
#  eqtl_signif_info=merge(reQTL_Signif_both,SNP_info[,.(snps=ID,AFB,EUB,ASH)],by="snps",all.x=T)
#}
#
## compute MAF
#eqtl_signif_info[,MAF_AFB:=pmin(AFB,1-AFB)]
#eqtl_signif_info[,MAF_EUB:=pmin(EUB,1-EUB)]
#eqtl_signif_info[,MAF_ASH:=pmin(ASH,1-ASH)]
#
#mediation$eGene=eqtl_signif_info[celltype%in%lineage_order&(MAF_AFB>0.05|MAF_EUB>0.05),setNames(rep(T,nrow(eQTL_Signif_both)),sprintf("%s%s%s",celltype,state,gene_name))][mediation[,sprintf("%s%s%s",celltype,state,Symbol)]]
#mediation[,eGene:=ifelse(is.na(eGene),F,eGene)]
#
#mediation[,FDR_global:=p.adjust(pval,"fdr")]
#
#if(exp_resp=="DE"){
#  fwrite(mediation,file=sprintf("%s/2_population_differences/mediation_analysis/mediation_w_betapop_popDE.tsv",DAT_DIR))
#}else{
#  fwrite(mediation,file=sprintf("%s/2_population_differences/mediation_analysis/mediation_w_betapop_popDR.tsv",DAT_DIR))
#}

mediation_DE=fread(sprintf("%s/2_population_differences/mediation_analysis/mediation_w_betapop_popDE.tsv",DAT_DIR))
mediation_DR=fread(sprintf("%s/2_population_differences/mediation_analysis/mediation_w_betapop_popDR.tsv",DAT_DIR))

mediation=rbind(mediation_DE,mediation_DR)

fwrite(mediation[,.(comp,type,celltype,state,Symbol,ID,betapop,mediator_type,mediator,frac_var,pval,FDR=FDR_global,eGene,betaeqtl,peqtl)],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V8/suppTables/TableS6/TableS6_mediation_analyses.tsv',EIP),sep='\t')

#color_mediator=c(color_cellTypes_24level[celltype_mediators],'Genetics'=grey(.3))
color_mediator=lineage_color[lineage_celltype[,setNames(lineage,celltype)][celltype_mediators]]
attr(color_mediator,"names")=celltype_mediators
color_mediator=c(color_mediator,'Genetics'=grey(.3))


fig4b_data_All<-mediation[,.(N=mean(FDR_global<.01), Npraw1=mean(pval<.01),
			     frac_var_mean=mean(frac_var_std),frac_var_meanpraw1=mean(frac_var_std[pval<.01])),by=.(celltype,state,mediator_type,type,mediator=ifelse(grepl('rs|ss|esv',mediator),'Genetics',mediator))]
fig4b_data_All[,state:=factor(state,c('NS','COV','IAV'))]
fig4b_data_All[,celltype:=factor(celltype,lineage_order)]
fig4b_data_All[,set:=ifelse(type=='DE','popDE','popDR')]
fig4b_data_All[,set:='popDE/DR']


fig4b_data_eGene <- mediation[eGene==TRUE,.(N=mean(FDR_global<.01), Npraw1=mean(pval<.01),frac_var_mean=mean(frac_var_std),frac_var_meanpraw1=mean(frac_var_std[pval<.01])),by=.(celltype,state,mediator_type,type,mediator=ifelse(grepl('rs|ss|esv',mediator),'Genetics',mediator))]
fig4b_data_eGene[,state:=factor(state,c('NS','COV','IAV'))]
fig4b_data_eGene[,celltype:=factor(celltype,lineage_order)]
fig4b_data_eGene[,set:=ifelse(type=='DE','eQTL-popDE','reQTL-popDR')]
fig4b_data_eGene[,set:='eQTL-popDE/DR']

fig4b_data=rbind(fig4b_data_All,fig4b_data_eGene)
fig4b_data[,set:=factor(set,c('popDE/DR','eQTL-popDE/DR'))]

fwrite(fig4b_data,sprintf("%s/Fig4b_data.tsv",FIGURE_DIR),sep="\t")

fig4b_data[,celltype:=factor(celltype,rev(lineage_order))]

fig4b_plot<-ggplot(fig4b_data[type=="DE",],aes(frac_var_mean*100,celltype, size=N*100))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0),alpha=0.6)+
  scale_fill_manual(values=color_mediator,guide='none')+
  scale_size_continuous(name="Significantly mediated popDEGs (%)",breaks = (c(10, 50, 100, 200, 500, 1000)/10),labels=paste0(c(10, 50, 100, 200, 500, 1000)/10,'%'),
		  guide=guide_legend(nrow=1,title.position="top",title.hjust=0.5))+
  facet_grid(cols=vars(state),rows=vars(set),labeller=labeller(set=c("popDE/DR"="All popDEGs","eQTL-popDE/DR"="eQTL popDEGs")))+
  scale_x_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  scale_y_discrete(labels=rev(c("MONO","B","CD4+ T","CD8+ T","NK")),breaks=rev(c("MONO","B","T.CD4","T.CD8","NK")),position="left")+
  xlab("Average popDE explained (%)")+
  theme_yann()+
  theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),
				axis.text.x=element_text(angle=0,hjust=0.5),legend.title=element_text(size=5),axis.title.y=element_blank())#+
#  theme(legend.spacing.x= unit(1, 'mm'),legend.spacing.y=unit(2, 'mm'),legend.text=element_text(size=6),
#	legend.box = "vertical", legend.margin=margin(),legend.key.size = unit(1, 'mm'),axis.title.y=element_blank())

fig4b_legend<-cowplot::get_legend(fig4b_plot)

################################################################################

pname=sprintf("%s/Fig4.pdf",FIGURE_DIR)
pdf(pname,width=6.7,height=7.2)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig4a_plot+theme(legend.position="none",text=element_text(size=7))),
      fig4a_legend,
      ggplotGrob(fig4b_plot+theme(legend.position="none",text=element_text(size=7))),
      fig4b_legend,
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(1,3,3),
      c(2,4,4),
      c(5,5,5),
      c(5,5,5)
    ), heights=c(3.5,0.5,3.5,2.5), widths=c(4,2,4)
  )
dev.off()










###### popDE only
fig4d_plot<-ggplot(fig4d_data[type=='DE',],aes(celltype,frac_var_mean*100,size=N*100))+

  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0),alpha=0.5)+
  #scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 200, 500, 1000)/10, labels=paste(c(10, 50, 100, 200, 500, 1000)/10,'%'), max_size = 4)+

  # scale_alpha_manual(values=setNames(c(0.75,0.25),c(T,F)))+
#  facet_grid(cols=vars(state))+
  #facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  facet_grid(cols=vars(state),rows=vars(set))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDE/DR explained (%)")+
  theme_yann()+
  guides(fill=guide_legend(nrow=2,order=1,override.aes = list(size = 2)),alpha='none')+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())+
  theme(legend.spacing.x= unit(1, 'mm'),legend.spacing.y=unit(2, 'mm'),legend.text=element_text(size=6),legend.key = element_blank(),legend.box = "vertical", legend.margin=margin(),legend.key.size = unit(1, 'mm'),)


dir.create(sprintf('%s/Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_overview_DE__All_and_eGenes.pdf',FIG_DIR),width=.4*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()


###### invert size and y
fig4d_plot<-ggplot(fig4d_data[(type=='DE' & state=='NS') | type=='DR',],aes(celltype,y=N*100,size=frac_var_mean*100))+

  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0),alpha=0.5)+
  #scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 200, 500, 1000)/10, labels=paste(c(10, 50, 100, 200, 500, 1000)/10,'%'), max_size = 4)+

  # scale_alpha_manual(values=setNames(c(0.75,0.25),c(T,F)))+
#  facet_grid(cols=vars(state))+
  #facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  facet_grid(cols=vars(state),rows=vars(set))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDE/DR explained (%)")+
  theme_yann()+
  guides(fill=guide_legend(nrow=2,order=1,override.aes = list(size = 2)),alpha='none')+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())+
  theme(legend.spacing.x= unit(1, 'mm'),legend.spacing.y=unit(2, 'mm'),legend.text=element_text(size=6),legend.key = element_blank(),legend.box = "vertical", legend.margin=margin(),legend.key.size = unit(1, 'mm'),)


dir.create(sprintf('%s/Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_overview_DER_eGene_PctGene_as_y__All_and_eGenes.pdf',FIG_DIR),width=.4*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()



###### invert size and y popDE only
fig4d_plot<-ggplot(fig4d_data[type=='DE',],aes(celltype,y=N*100,size=frac_var_mean*100))+

  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0),alpha=0.5)+
  #scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 200, 500, 1000)/10, labels=paste(c(10, 50, 100, 200, 500, 1000)/10,'%'), max_size = 4)+

  # scale_alpha_manual(values=setNames(c(0.75,0.25),c(T,F)))+
#  facet_grid(cols=vars(state))+
  #facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  facet_grid(cols=vars(state),rows=vars(set))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDE/DR explained (%)")+
  theme_yann()+
  guides(fill=guide_legend(nrow=2,order=1,override.aes = list(size = 2)),alpha='none')+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())+
  theme(legend.spacing.x= unit(1, 'mm'),legend.spacing.y=unit(2, 'mm'),legend.text=element_text(size=6),legend.key = element_blank(),legend.box = "vertical", legend.margin=margin(),legend.key.size = unit(1, 'mm'),)


dir.create(sprintf('%s/Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_overview_DE_PctGene_as_y__All_and_eGenes.pdf',FIG_DIR),width=.4*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()



fig4dbis_plot<-ggplot(fig4d_data[(type=='DR' & state=="COV"),],aes(celltype,frac_var_mean*100,size=N,alpha=.5))+

  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0))+
  #scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 500, 1000)/1000)+

  # scale_alpha_manual(values=setNames(c(0.75,0.25),c(T,F)))+
#  facet_grid(cols=vars(state))+
  #facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  facet_grid(cols=vars(set))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDR explained (%)")+
  theme_yann()+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())

dir.create(sprintf('%s/Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_overview_DE_NS__All_and_eGenes.pdf',FIG_DIR),width=.4*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()


#no eGene, no bin
fig4d_data <- mediation[,.(N=mean(FDR_global<.01), Npraw1=mean(pval<.01),frac_var_mean=mean(frac_var_std),frac_var_meanpraw1=mean(frac_var_std[pval<.01])),by=.(celltype,state,mediator_type,type,eGene,mediator=ifelse(grepl('rs|ss|esv',mediator),'genetics',mediator))]
fig4d_data[,state:=factor(state,c('NS','COV','IAV'))]
fig4d_data[,celltype:=factor(celltype,lineage_order)]
fig4d_plot<-ggplot(fig4d_data[(type=='DE' & state=="NS") | type=='DR',],aes(celltype,frac_var_mean*100,size=N,alpha=.5))+

  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0))+
  #scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 500, 1000)/1000)+

  # scale_alpha_manual(values=setNames(c(0.75,0.25),c(T,F)))+
#  facet_grid(cols=vars(state))+
  #facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  facet_grid(rows=vars(state),cols=vars(eGene))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDE/R explained (%)")+
  theme_yann()+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())

dir.create(sprintf('%s/Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_overview_DER_spliteGene.pdf',FIG_DIR),width=.4*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()



fig4d_plot<-
#  ggplot(fig4d_data[(type=='DE' & state=="NS") | type=='DR',],aes(paste(celltype,bin),frac_var_mean*100,size=N,alpha=bin))+
  ggplot(fig4d_data[eGene==TRUE,][(type=='DE' & state=="NS") | type=='DR',],aes(paste(celltype),frac_var_mean*100,size=N,alpha=.5))+

  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0))+
  #scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 500, 1000))+
  # scale_alpha_manual(values=setNames(c(0.75,0.25),c(T,F)))+
#  facet_grid(cols=vars(state))+
  #facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  facet_grid(cols=vars(state))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDE explained (%)")+
  theme_yann()+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())

dir.create(sprintf('%s/Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_overview_DER_eGene.pdf',FIG_DIR),width=.4*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()

mediation[,mediator_short:=ifelse(grepl('rs|ss|esv',mediator),'genetics',mediator)]
mediation[,state:=factor(state,c('NS','COV','IAV'))]
fig4d_plot<-
#  ggplot(fig4d_data[(type=='DE' & state=="NS") | type=='DR',],aes(paste(celltype,bin),frac_var_mean*100,size=N,alpha=bin))+
  ggplot(mediation[(type=='DE' & state=="NS") | type=='DR',],aes(paste(celltype,mediator_short),frac_var_std*100,alpha=.5))+
  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_violin(aes(fill=mediator_short),position=position_dodge(width=0),scale='width')+
  geom_boxplot(notch=T,fill='white',alpha=0)+
  #scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 500, 1000))+
  # scale_alpha_manual(values=setNames(c(0.75,0.25),c(T,F)))+
#  facet_grid(cols=vars(state))+
  #facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  facet_grid(cols=vars(state))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDE explained (%)")+
  theme_yann()+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())

dir.create(sprintf('%s/Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_violin_DER.pdf',FIG_DIR),width=.7*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()



################################################################################
############################ don't split eQTL and others #######################
################################################################################

res1<-mediation[FDR_global<0.01,.N,by=.(celltype,state,mediator_type,type,mediator=ifelse(grepl('rs',mediator),'genetics',mediator))]
res2<-mediation[,mean(frac_var_std),by=.(celltype,state,mediator_type,type,mediator=ifelse(grepl('rs',mediator),'genetics',mediator))]
setnames(res2,"V1","frac_var_mean")


fig4d_data<-merge(res1,res2,by=c("celltype","state","mediator_type","type","mediator"))
# fig4d_data[grep('rs',mediator),mediator:='genetics']

fig4d_plot<-
  ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N))+
  #ggplot(fig4d_data[state=="NS"],aes(celltype,frac_var_mean*100,size=N,group=eGene_mod,alpha=eGene_mod))+
  geom_point(color="black",pch=21,position=position_dodge(width=0),alpha=1)+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0),alpha=0.5)+
  scale_alpha_manual(values=setNames(c(0.15,0.3,0.45,0.6,0.75),bin_order))+
  scale_fill_manual(values=color_mediator)+
  scale_size_area(breaks = c(10, 50, 100, 500, 1000))+
#  facet_grid(cols=vars(state))+
  # facet_grid(cols=vars(eGene),labeller=labeller(eGene=setNames(c("eQTL-popDEG","popDEG"),c(TRUE,FALSE))))+
  scale_y_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  ylab("Avg. popDE explained (%)")+
  theme_yann()+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())

  dir.create(sprintf('%s//Fig4/',FIG_DIR))
pdf(sprintf('%s/Fig4/Fig4d_mediationPct_overview_DE_nosplit.pdf',FIG_DIR),width=.3*6.7,height=.35*7.2)
print(fig4d_plot)
dev.off()

################################################################################
# Fig. 4e

fig4e_data<-mediation[FDR_global<0.01]
setnames(fig4e_data,"celltype","lineage")
fig4e_data=dcast(fig4e_data,lineage+state+Symbol+ID+type~mediator_type,value.var="frac_var_std",fill=0)
#fig4e_data[,ratio_celltype:=ifelse(genetics==0,celltype,celltype/pmin((celltype+genetics),1))]
#fig4e_data[,ratio_genetics:=ifelse(celltype==0,genetics,genetics/pmin((celltype+genetics),1))]
fig4e_data[,ratio_celltype:=ifelse(genetics==1,0,celltype)]
fig4e_data[,ratio_genetics:=ifelse(celltype==1,0,genetics)]
fig4e_data[,idx:=sprintf("%s_%s_%s_%s_%s_%s",lineage,state,Symbol,ID,celltype,genetics)]

ratio_corrected=fig4e_data[ratio_celltype+ratio_genetics>1,.(idx,ratio_celltype=ratio_celltype/(ratio_celltype+ratio_genetics),ratio_genetics=ratio_genetics/(ratio_celltype+ratio_genetics))]

ratio_corrected_celltype=ratio_corrected[,setNames(ratio_celltype,idx)]
ratio_corrected_genetics=ratio_corrected[,setNames(ratio_genetics,idx)]

fig4e_data[,ratio_celltype_corrected:=ifelse(idx%in%names(ratio_corrected_celltype),ratio_corrected_celltype[idx],ratio_celltype)]
fig4e_data[,ratio_genetics_corrected:=ifelse(idx%in%names(ratio_corrected_genetics),ratio_corrected_genetics[idx],ratio_genetics)]

fig4e_data[,unexplained:=1-(ratio_genetics_corrected+ratio_celltype_corrected)]
#fig4e_data[,unexplained_melt:=unexplained]
fig4e_data_melt=melt(fig4e_data[,.(lineage, state, Symbol, ID, type,ratio_genetics,ratio_celltype, celltype=ratio_celltype_corrected,genetics=ratio_genetics_corrected,unexplained)], measure.vars=c("celltype","genetics","unexplained"),variable.name="mediator_type",value.name="ratio")


fig4e_data_melt=melt(fig4e_data,measure.vars=c("ratio_celltype_corrected","ratio_genetics_corrected","unexplained_melt"),variable.name="mediator_type",value.name="ratio")
# fig4e_data_melt[,type:=gsub("^ratio_","",mediator_type)]
# fig4e_data_melt[,type:=gsub("_corrected$","",mediator_type)]
fig4e_data_melt[,cellstate:=sprintf("%s_%s",lineage,state)]
# fig4e_data_melt[,cellstate_Symbol:=sprintf("%s_%s",cellstate,Symbol)]
fig4e_data_melt_sub=fig4e_data_melt[order(ratio_genetics+ratio_celltype),]
#fig4e_data_melt_sub=fig4e_data_melt[cellstate%in%c("NK_NS","NK_COV","NK_IAV")][order(-unexplained),]
fig4e_data_melt_sub=fig4e_data_melt[cellstate%in%c("T.CD4_NS","NK_COV")]#[order(-unexplained),]
# fig4e_data_melt_sub[,mediator_type:=case_when(
#   mediator_type=="celltype"&lineage=="MONO"~"MONO.CD16",
#   mediator_type=="celltype"&lineage=="B"~"B.M.K",
#   mediator_type=="celltype"&lineage=="T.CD4"~"T.CD4.E",
#   mediator_type=="celltype"&lineage=="T.CD8"~"T.CD8.EMRA",
#   mediator_type=="celltype"&lineage=="NK"~"NK.M.LIKE",
#   T~mediator_type)]
# fig4e_data_melt_sub[,Symbol:=factor(Symbol,unique(Symbol))]
# fig4e_data_melt_sub[,cellstate_Symbol:=factor(cellstate_Symbol,unique(cellstate_Symbol))]

fig4e_plot=ggplot(fig4e_data_melt_sub[type=='DE',])+
  geom_col(aes(ratio*100,cellstate_Symbol,fill=type.1),position="stack",color=NA)+
#  facet_grid(cols=vars(lineage),rows=vars(state),scales="free")+
#  geom_col(aes(ratio*100,Symbol,fill=type),position="stack",color=NA)+
  facet_grid(rows=vars(lineage),scales="free")+
  ylab('popDEGs')+
  xlab('popDE explained (%)')+
  scale_fill_manual(values=c(celltype_color,"genetics"=gray(.3),"unexplained"="lightgray"),breaks=c("T.CD4.E","genetics","unexplained"),labels=c("CD4+ T E","eQTL","Unexplained"),guide=guide_legend(nrow=1,override.aes=list(size=0.2,pch=21,color="black")))+
  theme_yann()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing=unit(0,"mm"))

  pdf(sprintf('%s/Fig4/Fig4e_mediationPct_overview_DE.pdf',FIG_DIR),width=.3*6.7,height=.35*7.2)
  print(fig4e_plot)
  dev.off()
