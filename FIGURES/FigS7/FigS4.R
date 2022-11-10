EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"
DAT_DIR=sprintf("%s/single_cell/project/pop_eQTL/data",EIP)
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EIP))
PAP_DIR=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V2",EIP)
FIG_DIR=sprintf("%s/figureMaterial",PAP_DIR)

celltype_color=c(
"T.CD4.N"="#169FD8","T.CD4.E"="#005274","T.Reg"="#03C2C0","T.gd"="#ffce73",
"T.CD8.N"="#00BB54","T.CD8.CM.EM"="#0EAD20","T.CD8.EMRA"="#048c41","MAIT"="#004A21",
"NK.M.LIKE"="#63C425","NK.CD56dim"="#B7DC2A","NK.CD56brt"="#7D9E00","ILC"="#3F4F00",
"B.N.K"="#9281DE","B.N.L"="#a681de","B.M.K"="#6045d9","B.M.L"="#8045d9",
"Plasmablast"="#6801A4","MONO.CD14"="#B7637E","MONO.CD16"="#8F568A","pDC"="#8E6B00",
"cDC"="#E7b315","MONO.CD14.INFECTED"="#B9364B")

lineage_order=c("MONO","B","T.CD4","T.CD8","NK")

lineage_label=c("MONO","B","CD4+ T","CD8+ T","NK")
attr(lineage_label,"names")=c("MONO","B","T.CD4","T.CD8","NK")

celltype_order=c("MONO.CD14","MONO.CD16","MONO.CD14.INFECTED","cDC","pDC","B.N.K","B.N.L","B.M.K","B.M.L","Plasmablast","T.CD4.N","T.CD4.E","T.Reg","T.gd","T.CD8.N","T.CD8.CM.EM","T.CD8.EMRA","ILC","MAIT","NK.CD56dim","NK.CD56brt","NK.M.LIKE")
celltype_label=c("MONO CD14+","MONO CD16+","MONO IAV+","cDC","pDC","B N k","B N l","B M k","B M l","Plasmablast","T CD4+ N","T CD4+ E","T Reg","T gd","T CD8+ N","T CD8+ CM/EM","T CD8+ EMRA","ILC","MAIT","NK CD56dim","NK CD56brt","NK mem")
celltype_label_exp=expression("CD14+ monocyte","CD16+ monocyte","IAV+ monocyte","cDC","pDC",kappa*"-LC naive B",lambda*"-LC naive B",kappa*"-LC memory B",lambda*"-LC memory B","Plasmablast","CD4+ naive T","CD4+ effector T","Regulatory T",gamma*delta*" T","CD8+ naive T","CD8+ CM/EM T","CD8+ EMRA T","ILC","MAIT","CD56"["dim"]*" NK","CD56"["brt"]*" NK","Memory-like NK")
names(celltype_label_exp)=celltype_order

################################################################################
# Fig. S5a
meta=fread(sprintf("%s/2_population_differences/sce_clean_metadata.tsv",DAT_DIR))
meta$POP=substr(meta$IID,1,3)

lineage_order=c("MONO","B","T.CD4","T.CD8","NK")

celltype_focus=c("MONO.CD16","B.M.K","T.CD4.E","T.CD8.EMRA","NK.M.LIKE")

celltype_alpha_fill=setNames(ifelse(names(celltype_color)%in%c("MONO.CD16","B.M.K","T.CD4.E","T.CD8.EMRA","NK.M.LIKE"),0.9,0.5),names(celltype_color))
celltype_alpha_color=setNames(ifelse(names(celltype_color)%in%c("MONO.CD16","B.M.K","T.CD4.E"),1,0),names(celltype_color))

celltype_freq=lapply(unique(meta$lineage),function(X){
  compartment_sizes=meta[POP!="ASH"&lineage==X,.N,keyby=.(IID,lineage)]
  setnames(compartment_sizes,"N","TOT")
  iid_sizes=meta[POP!="ASH",.N,keyby=.(IID)]
  setnames(iid_sizes,"N","TOT_IID")
  celltype_sizes=meta[POP!="ASH"&lineage==X,.N,keyby=.(IID,celltype)]
  celltype_freq=merge(celltype_sizes,compartment_sizes,by="IID")
  celltype_freq=merge(celltype_freq,iid_sizes,by="IID")
  celltype_freq[,FREQ_CT:=N/TOT*100]
  celltype_freq[,FREQ_ID:=N/TOT_IID*100]
  celltype_freq[,POP:=substr(IID,1,3)]
  celltype_freq[,lineage:=X]
  return(celltype_freq)
})%>%{do.call("rbind",.)}

celltype_freq$lineage=factor(celltype_freq$lineage,lineage_order)

figs5a_data=celltype_freq[celltype%in%celltype_focus,]
figs5a_data[,celltype:=factor(celltype,celltype_focus)]

figs5a_plot=ggplot(figs5a_data)+
  geom_boxplot(aes(celltype,FREQ_CT,fill=celltype,group=POP),color=NA,outlier.shape=NA,show.legend=F,alpha=0.5)+
  geom_boxplot(aes(celltype,FREQ_CT,color=celltype,group=POP),fill=NA,outlier.shape=NA,show.legend=F,lwd=0.2)+
  geom_point(aes(celltype,FREQ_CT,fill=celltype),alpha=0)+
  ylab("Fraction of immune lineage")+
  xlab("Cell type")+
  scale_y_continuous(breaks=c(0,100),labels=c(0,1))+
scale_x_discrete(labels=c("CD16+ monocyte","Memory B","CD4+ effector T","CD8+ EMRA T","Memory-like NK"))+
  scale_fill_manual(values=celltype_color,guide="none")+
  scale_color_manual(values=celltype_color,guide="none")+
  facet_grid(cols=vars(lineage),labeller=labeller(lineage=c("MONO"="MONO","B"="B","T.CD4"="CD4+ T","T.CD8"="CD8+ T","NK"="NK")),scales="free_x")+
  theme_yann()+
  theme(text=element_text(size=6),panel.spacing=unit(0,"pt"))
 
celltypes_show=celltype_order%in%c("T.CD4.E","MONO.CD16","B.M.K","NK.M.LIKE","T.CD8.EMRA")
figs5a_legend=ggplot(celltype_freq)+
  geom_point(aes(POP,FREQ_CT,fill=celltype))+
  scale_fill_manual(values=celltype_color[celltype_order][celltypes_show],
                    breaks=celltype_order[celltypes_show],
                    labels=celltype_label_exp[celltype_order][celltypes_show],
                    name="",
                    guide=guide_legend(override.aes=list(shape=21,size=2,color="black"),ncol=1,byrow=T))+
  theme(legend.position="bottom")+
  theme(text=element_text(size=5),panel.spacing=unit(0,"pt"),axis.title=element_blank(),legend.spacing.y=unit(-2,'mm'),legend.spacing.x=unit(-0,'mm'))
figs5a_legend=get_legend(figs5a_legend)

pname=sprintf("%s/FigS5/Figs5b_boxplot.pdf",FIG_DIR)
pdf(pname,width=6,height=5)
  grid.arrange(
    grobs=list(
      ggplotGrob(figs5a_plot+theme(legend.position="none",text=element_text(size=7))),
      figs5a_legend#,
#      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
#      c(3,1,3),
      c(1,2)
    ), widths=c(7,3)#, widths=c(2,6,2)
  )
dev.off()

################################################################################
# Fig. S5b
library(ggbeeswarm)

cmv_status=fread(sprintf("%s/single_cell/project/pop_eQTL/Covariates/CMV_all.txt",EIP))
setnames(cmv_status,c("CMV (y/n)","Standard Units"),c("CMV","standard_units"))
cmv_status=cmv_status[type=="sample",.(IID=IID_w_ASH,CMV,standard_units)]
cmv_status[,standard_units:=as.numeric(gsub(",",".",standard_units))]
cmv_status$CMV=ifelse(cmv_status$CMV=="y","CMV","NO_CMV")
cmv_status[,POP:=factor(substr(IID,1,3),c("AFB","EUB","ASH"))]

figs5b_plot=ggplot(cmv_status,aes(POP,standard_units,color=POP,fill=POP))+
  geom_hline(yintercept=11,size=0.1)+
  #geom_point(alpha=0.5,position=position_jitter(height=0),size=0.8)+
  geom_quasirandom(alpha=0.5,size=0.8)+
  theme_yann()+ylab("Standard units")+
  scale_color_manual(aesthetics=c("fill","color"),values=color_populations,guide="none")+
  theme(text=element_text(size=7),axis.title.x=element_text(color="white"))

################################################################################
# Fig. S5c
pct_celltype_sample=fread(sprintf("%s/2_population_differences/Cellcounts/Pct_celltype_by_IID_condition.tsv.gz",DAT_DIR))
celltypes=c("NK.M.LIKE","T.CD8.EMRA")
pct_celltype_sample=pct_celltype_sample[,mget(c("IID","state",celltypes))]

pct_celltype_sample$CMV=cmv_status[,setNames(CMV,IID)][pct_celltype_sample$IID]
pct_celltype_sample=pct_celltype_sample[!is.na(CMV),]
pct_celltype_sample[,POP:=substr(IID,1,3)]
pct_celltype_sample[,CMV_POP:=sprintf("%s_%s",CMV,POP)]
figs5c_data=melt(pct_celltype_sample,measure.vars=celltypes,value.name="fraction",variable.name="celltype")
figs5c_data[,POP:=factor(POP,c("AFB","EUB","ASH"))]
figs5c_data[,CMV_POP:=factor(CMV_POP,c("CMV_AFB","NO_CMV_AFB","CMV_EUB","NO_CMV_EUB","CMV_ASH","NO_CMV_ASH"))]

figs5c_plot=ggplot(figs5c_data[state=="NS"&CMV_POP!="NO_CMV_AFB",],aes(CMV_POP,fraction))+
  geom_violin(aes(fill=celltype),color=NA,alpha=0.5,scale="width")+
  geom_violin(aes(color=celltype),fill=NA,size=0.1,scale="width",draw_quantiles=0.5)+
  ylab("Fraction of all cells")+
  facet_grid(cols=vars(POP),rows=vars(celltype),scales="free_x",space="free_x",labeller=labeller(celltype=c("NK.M.LIKE"="Memory-like NK","T.CD8.EMRA"="CD8+ EMRA T")))+
  scale_color_manual(aesthetics=c("color","fill"),values=celltype_color,guide="none")+
  theme_yann()+theme(text=element_text(size=7),axis.title.x=element_text(color="white"),panel.spacing=unit(0,"mm"))

################################################################################
# Fig. S5d
exp_resp="DE"

memory_markers=fread(sprintf("%s/suppTables/TableS3C_DE_memory_vs_nonmemory_subsets.tsv",PAP_DIR))
memory_markers[,reg:=ifelse(beta.tt>0,"upreg","dnreg")]

qrs=memory_markers[,unique(query)]

memory_markers[,celltype:=case_when(
  query=="MONO.CD16"~"MONO",
  query=="B.M.K"~"B",
  query=="T.CD4.E"~"T.CD4",
  query=="T.CD8.EMRA"~"T.CD8",
  query=="NK.M.LIKE"~"NK"
)]

bgs=c("MONO.CD14","B.N.K","T.CD4.N","T.CD8.N","NK.CD56dim")

ANALYSE=sprintf("lineage_condition_AFBEUB__%slineage_condition%s__220409",ifelse(exp_resp=="DR","logFC_",""),ifelse(exp_resp=="DR","_logFC",""))
popde_raw=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popde_raw[,type:="expression"]
popde_raw=popde_raw[,.(Symbol,state,betapop.raw=beta,FDRpop.raw=FDR,celltype)]
ANALYSE="lineage_condition_AFBEUB__lineage_condition__CellPropLineage_220409"
ANALYSE=sprintf("lineage_condition_AFBEUB__%slineage_condition%s__CellPropLineage_220409",ifelse(exp_resp=="DR","logFC_",""),ifelse(exp_resp=="DR","_logFC",""))
popde_adj=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popde_adj[,type:="expression"]
popde_adj=popde_adj[,.(Symbol,state,betapop.adj=beta,FDRpop.adj=FDR,celltype)]

popdiff_changes=fread(sprintf("%s/2_population_differences/cell_composition_adjustment_beta_changes.tsv",DAT_DIR))
popdiff_losses=popdiff_changes[difference_signif=="Different"&adjustment_phase=="no_lineage"&type==exp_resp,.(celltype,state,Symbol,group="loss")]
memory_markers_sub=merge(memory_markers[type==exp_resp&background%in%bgs,],popdiff_losses,by=c("state","Symbol","celltype"),all.x=T)
memory_markers_sub[,group:=ifelse(is.na(group),"no_loss",group)]

memory_markers_sub=merge(memory_markers_sub,popde_raw,by=c("state","Symbol","celltype"),all.x=T)
memory_markers_sub=merge(memory_markers_sub,popde_adj,by=c("state","Symbol","celltype"),all.x=T)

figs5d_data=memory_markers_sub
figs5d_data[,state:=factor(state,c("NS","COV","IAV"))]
figs5d_data[,celltype:=factor(celltype,lineage_order)]

figs5d_plot=ggplot(figs5d_data[FDRpop.raw<0.01&abs(betapop.raw)>0.2&FDR.wx<0.01],aes(abs(beta.tt),abs(betapop.raw-betapop.adj),color=state,fill=state))+
  geom_point(size=0.1,alpha=0.4)+
  geom_smooth(formula=y~x,method="lm",size=0.2)+
  scale_fill_manual(values=cond_color,guide="none")+
  scale_color_manual(values=cond_color,guide="none")+
  xlab(expression("|log"[2]*"FC(memory/non-memory subsets)|"))+
  ylab(expression("|log"[2]*"FC"["r"]*"(EUB/AFB)-log"[2]*"FC"["a"]*"(EUB/AFB)|"))+
  facet_grid(cols=vars(celltype),rows=vars(state),labeller=labeller(celltype=lineage_label))+
  theme_yann()+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"))

figs5d_legend=ggplot(figs5d_data[FDRpop.raw<0.01&abs(betapop.raw)>0.2&FDR.wx<0.01],aes(abs(beta.tt),abs(betapop.raw-betapop.adj),color=state,fill=state))+
  geom_point(size=0.1,alpha=0.4)+
  scale_color_manual(values=cond_color,guide="none")+
  theme_yann()+
  scale_fill_manual(values=cond_color,guide=guide_legend(nrow=1,title="Condition",title.position="left",title.vjust=0.5,override.aes=list(color="black",size=1,pch=21,alpha=1)))



figs5d_legend=cowplot::get_legend(figs5d_legend+theme(text=element_text(size=7)))

################################################################################
# other option
#nbins=25
#increment=figs5d_data[,max(beta.tt)/nbins]
#breaks=figs5d_data[,seq(min(beta.tt)-increment,max(beta.tt)+increment,increment)]
#figs5d_data$bin=cut(figs5d_data$beta.tt,breaks,include.lowest=T,right=F)
#
#figs5dbis_data=melt(figs5d_data,measure.vars=c("betapop.raw","betapop.adj"),value.name="betapop",variable.name="adj")
#figs5dbis_data=figs5dbis_data[!is.na(betapop),mean(betapop),by=.(bin,celltype,state,adj)]
#
#figs5d_plot=ggplot(figs5dbis_data,aes(bin,V1,group=adj,fill=adj))+
#  geom_col(position="dodge")+
#  scale_fill_manual(values=c("betapop.raw"="blue","betapop.adj"="red"),guide="none")+
#  facet_grid(cols=vars(celltype),rows=vars(state),labeller=labeller(celltype=lineage_label))+
#  theme_yann()+theme(text=element_text(size=7),panel.spacing=unit(0,"mm"))
################################################################################

pname=sprintf("%s/FigureS5/Figs5d.pdf",FIG_DIR)
pdf(pname,width=6,height=5)
  grid.arrange(
    grobs=list(
      ggplotGrob(figs5a_plot+theme(legend.position="none",text=element_text(size=7))),
      ggplotGrob(figs5b_plot+theme(legend.position="none",text=element_text(size=7))),
      ggplotGrob(figs5c_plot+theme(legend.position="none",text=element_text(size=7))),
      ggplotGrob(figs5d_plot+theme(legend.position="none",text=element_text(size=7)))
    ),
    layout_matrix=rbind(
      c(1,2,3),
      c(4,4,4)
    ), heights=c(4.5,5.5), widths=c(3.5,2,3.5)
  )
dev.off()
