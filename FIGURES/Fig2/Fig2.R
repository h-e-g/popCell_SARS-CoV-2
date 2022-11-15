################################################################################
################################################################################
# File name: Fig2.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Figure 2
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="LIBRARY"
source(sprintf("./2a__popDEGs_popDRGs__lib.R",LIB_DIR))

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
# Fig. 2a

# load meta data created in 1b1__celltype_identification.R
meta_clean=fread(sprintf("%s/data/sce_clean_metadata.tsv",AGGR_QC_DIR))

meta_clean[,POP:=substr(IID,1,3)]
setnames(meta,"condition","state")

lineage_celltype=fread("1__transcriptome_processing/data/lineage_celltype.tsv")

lineages=lineage_order
attr(lineages,"names")=lineages

fig2a_data=lapply(lineages,function(X){
  if (X=="MONO") {
    celltypes=lineage_celltype[lineage==X&!grepl("INFECTED",celltype),celltype]
  } else {
    celltypes=lineage_celltype[lineage==X,celltype]
  }
  data=meta[POP!="ASH"&state=="NS"&lineage==X,]
  lineage_sizes=data[,.N,keyby=.(IID,lineage)]
  setnames(lineage_sizes,"N","TOT")
  celltype_sizes=data[,.N,keyby=.(IID,celltype=factor(celltype,celltypes))]
  celltype_sizes=dcast(celltype_sizes,IID~celltype,value.var="N",fill=0)
  celltype_sizes=melt(celltype_sizes,measure.vars=celltypes,value.name="N",variable.name="celltype")
  celltype_freq=merge(celltype_sizes,lineage_sizes,by="IID",all.x=T)
  celltype_freq[,FREQ:=N/TOT]
  celltype_freq[,POP:=substr(IID,1,3)]
  celltype_freq_mean=celltype_freq[,mean(FREQ),keyby=.(celltype,POP)]
  celltype_freq_mean[,lineage:=X]
  return(celltype_freq_mean)
})%>%rbindlist()

fig2a_data[,lineage:=factor(lineage,lineage_order)]
fig2a_data[,celltype:=factor(celltype,celltype_order)]

fig2a_plot=ggplot(fig2a_data)+
  geom_col(aes(POP,V1,fill=celltype),alpha=0.7,color=NA)+
  geom_col(aes(POP,V1,color=celltype),fill=NA,size=0.01)+
  theme_yann()+
  ylab("Fraction of immune lineage")+
  xlab("Population")+
  scale_fill_manual(values=celltype_color,guide="none")+
  scale_color_manual(values=celltype_color,guide="none")+
  scale_y_continuous(breaks=c(0,1),labels=c(0,1))+
  facet_grid(cols=vars(lineage),labeller=labeller(lineage=c("MONO"="Myeloid","B"="B","T.CD4"="CD4+ T","T.CD8"="CD8+ T","NK"="NK")))+
  theme(text=element_text(size=7),panel.spacing=unit(0,"pt"))


fig2a_legend=ggplot(fig2a_data)+
  geom_point(aes(POP,V1,fill=celltype))+
  theme_yann()+
  scale_fill_manual(values=celltype_color,
                    breaks=celltype_order,
                    labels=celltype_label_exp,
                    guide=guide_legend(override.aes=list(shape=21,size=2,color="black"),nrow=2,byrow=T))+
  theme(legend.position="right")+
  theme(text=element_text(size=7),panel.spacing=unit(0,"pt"),axis.title=element_blank())
fig2a_legend=get_legend(fig2a_legend)

################################################################################
# Fig. 2b

fig2b_data=fread(sprintf("%s/data/Fig2A_data_popDEGs_adjusted_betaraw_corrected.tsv",DAT_POPDIFF_DIR))
setnames(fig2b_data,c("value"),c("count"))
fig2b_data[,cellprop_adjust:=case_when(adjustment_phase=="n_popD"~"none",adjustment_phase=="n_popD_no_lineage"~"intralineage",adjustment_phase=="n_popD_lineage_all"~"crosslineage")]
fig2b_data$cellprop_adjust=factor(fig2b_data$cellprop_adjust,rev(c("crosslineage","intralineage","none")))

fig2b_data$state=factor(fig2b_data$state,condition_order)
fig2b_data$celltype=case_when(fig2b_data$celltype=="Monocyte"~"MYELOID",T~fig2b_data$celltype)
fig2b_data$celltype=factor(fig2b_data$celltype,rev(c("MYELOID","B","T.CD4","T.CD8","NK")))

fig2b_data[,cellprop_adjust:=as.character(cellprop_adjust)]
fig2b_data=fig2b_data[cellprop_adjust!="crosslineage",]
fig2b_data$cellprop_adjust=factor(fig2b_data$cellprop_adjust,rev(c("intralineage","none")))
fig2b_data=dcast(fig2b_data,celltype+state~cellprop_adjust,value.var="count")%>%
  .[,.(celltype,state,none=none-intralineage,intralineage)]%>%
  melt(measure.vars=c("none","intralineage"),variable.name="cellprop_adjust",value.name="count")

fig2b_plot=ggplot(fig2b_data)+
  geom_col(aes(count,celltype,fill=celltype,alpha=cellprop_adjust,group=cellprop_adjust),position="stack")+
  facet_grid(rows=vars(state),labeller=labeller(state=c("NS"="NS (logCPM)","COV"="COV (logFC)","IAV"="IAV (logFC)")))+
  theme_plot()+
  scale_fill_manual(values=celltype.6level_colors[-2],guide="none")+
  scale_y_discrete(labels=rev(c("Myeloid","B","CD4+ T","CD8+ T","NK")))+
  scale_alpha_manual(values=c("none"=0.6,"intralineage"=1),breaks=c("none","intralineage"),labels=c("Raw","Adjusted"),guide=guide_legend(override.aes=list(shape=21,color="black")))+
  xlab("Number of expression differences")+
  theme(text=element_text(size=7),axis.title.y=element_blank())+
  theme(panel.spacing=unit(0,"pt"),legend.position="bottom")

p=ggplot(fig2b_data)+
  geom_point(aes(celltype,count,fill=celltype,alpha=cellprop_adjust))+
  scale_fill_discrete(guide="none")+
  scale_alpha_manual(
    values=c("none"=0.6,"intralineage"=1),
    breaks=c("none","intralineage"),
    labels=c("Raw","Adjusted"),
    guide=guide_legend(override.aes=list(shape=21,color=NA,fill="black",size=2),nrow=1)
  )+
  theme_plot()+
  theme(legend.position="bottom")
fig2b_legend=get_legend(p)

################################################################################
# Fig. 2c

# define default parameters
NLIBS=125
CELLTYPE='lineage' # celltype variable to use. Will be used for naming of output files
STATE='condition'

# load batch-adjusted counts computed in 1c2__pseudobulk_batch_correction.R
Expr=fread(sprintf("1__transcriptome_processing/data/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz",NLIBS,CELLTYPE,STATE))

celltypes_to_plot=c("MONO","T.CD4")
genes_to_plot=c("GBP7","CCL23")

fig2c_data=Expr[POP!="ASH"&celltype%chin%celltypes_to_plot&Symbol%chin%genes_to_plot,]
fig2c_data[,state:=factor(state,condition_order)]
fig2c_data[,POP:=factor(POP,c("AFB","EUB"))]
fig2c_data[,celltype:=factor(celltype,c("MONO","T.CD4"))]
fig2c_data[,Symbol:=factor(Symbol,c("GBP7","CCL23"))]

fig2c_plot=ggplot(fig2c_data,aes(state,logCPM))+
  geom_boxplot(aes(fill=POP),alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  geom_boxplot(aes(color=POP),fill=NA,outlier.shape=NA,size=0.1,notch=T)+
  facet_grid(rows=vars(Symbol),cols=vars(celltype),labeller=labeller(celltype=c("MONO"="MYELOID","T.CD4"="CD4+ T")))+
  scale_color_manual(name="Population",aesthetics=c("color","fill"),values=color_populations,breaks=c("AFB","EUB"))+
  ylab("Gene expression (logCPM)")+
  theme_plot() +
  theme(axis.title.x=element_text(color="white"),legend.position="none")

################################################################################
# Fig. 2d

DAT_DIR=sprintf("%s/single_cell/project/pop_eQTL/data/2_population_differences",EIP)
POP_DIR=sprintf("%s/popDE",DAT_DIR)

fig2e_data=fread(sprintf("%s/NK.M.LIKE_markers.tsv",DAT_DIR))
fig2e_data[,reg:=ifelse(beta.tt>0,"upreg","dnreg")]

exp_resp="DE"

ANALYSE=sprintf("lineage_condition_AFBEUB__%slineage_condition%s__220409",ifelse(exp_resp=="DR","logFC_",""),ifelse(exp_resp=="DR","_logFC",""))
popde_raw=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popde_raw[,type:="expression"]
popde_raw=popde_raw[celltype=="NK",.(Symbol,state,betapop.raw=beta,FDRpop.raw=FDR)]
ANALYSE="lineage_condition_AFBEUB__lineage_condition__CellPropLineage_220409"
ANALYSE=sprintf("lineage_condition_AFBEUB__%slineage_condition%s__CellPropLineage_220409",ifelse(exp_resp=="DR","logFC_",""),ifelse(exp_resp=="DR","_logFC",""))
popde_adj=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popde_adj[,type:="expression"]
popde_adj=popde_adj[celltype=="NK",.(Symbol,state,betapop.adj=beta,FDRpop.adj=FDR)]

popdiff_losses=fread(sprintf("%s/cell_composition_adjustment_beta_changes.tsv",DAT_DIR))
popdiff_losses=popdiff_losses[difference_signif=="Different",]
popdiff_losses_nk=popdiff_losses[adjustment_phase=="no_lineage"&type==exp_resp&celltype=="NK",.(state,Symbol,group="loss")]
fig2e_data=merge(fig2e_data,popdiff_losses_nk,by=c("state","Symbol"),all.x=T)
fig2e_data[,group:=ifelse(is.na(group),"no_loss",group)]

fig2e_data=merge(fig2e_data,popde_raw,by=c("state","Symbol"),all.x=T)
fig2e_data=merge(fig2e_data,popde_adj,by=c("state","Symbol"),all.x=T)
fig2e_data[,facet:="NK cells"]

bg="NK.CD56dim"
st="COV"

fig2e_data[state==st&type==exp_resp&background==bg,cor.test(beta.tt,(betapop.adj-betapop.raw))]

genes_to_plot_r=c("CCL3","CCL4","IL2RA","IL2RB","IL18RAP","FCER1G","AREG","CD69","TMIGD2")
genes_to_plot_l=c("CADM1","TPRG1","LAG3")

fig2e_data_sub=fig2e_data[background==bg&type==exp_resp&FDRpop.raw<0.01&FDR.tt<0.01,]

library(ggnewscale)

fig2e_plot=ggplot()+
  scale_color_gradient(low=color_populations["AFB"],high=color_populations["EUB"],guide="none")+
  geom_rug(data=data.table(x=seq(-3,3,length.out=1000),y=seq(-3,3,length.out=1000)),mapping=aes(x,y,color=x),sides="b")+
  geom_hline(yintercept=0,size=0.1,color="gray",linetype="dashed")+
  geom_vline(xintercept=0,size=0.1,color="gray",linetype="dashed")+
  geom_abline(slope=1,size=0.1)+
  new_scale_color()+
  geom_point(data=fig2e_data_sub[state!="NS"&group=="no_loss"&abs(betapop.raw)>0.2&FDRpop.adj<0.01&abs(betapop.adj)>0.2],mapping=aes(betapop.raw,betapop.adj),color="gray",alpha=0.2)+
  geom_point(data=fig2e_data_sub[order(-beta.tt)][state==st&group=="loss"],mapping=aes(betapop.raw,betapop.adj,color=state,fill=state,shape=reg,alpha=scale(abs(beta.tt)),size=(abs(beta.tt))))+
  geom_text_repel(data=fig2e_data_sub[order(-beta.tt)][Symbol%in%genes_to_plot_r&state==st&group=="loss"],mapping=aes(betapop.raw,betapop.adj,color=state,label=Symbol),segment.alpha=0.5,segment.size=0.1,max.overlaps=Inf,size=1.5,xlim=c(2.5,3),ylim=c(-2,2),direction="y",seed=2014)+
  geom_text_repel(data=fig2e_data_sub[order(-beta.tt)][Symbol%in%genes_to_plot_l&state==st&group=="loss"],mapping=aes(betapop.raw,betapop.adj,color=state,label=Symbol),segment.alpha=0.5,segment.size=0.1,max.overlaps=Inf,size=1.5,xlim=c(0,-1),ylim=c(-1.5,-2.5),direction="y",seed=2014)+
  scale_color_manual(values=color_conditions,guide="none")+
  scale_fill_manual(values=color_conditions,guide="none")+
  scale_alpha_continuous(guide="none")+
  scale_shape_manual(values=c("upreg"=24,"dnreg"=25),breaks=c("upreg","dnreg"),labels=c("Upregulated","Downregulated"),guide=guide_legend(ncol=1))+
  scale_size_continuous(range=c(0,1.5),guide="none")+
  scale_x_continuous(limits=c(-2.85,2.85))+
  scale_y_continuous(limits=c(-2.85,2.85))+
  new_scale_color()+
#  scale_color_viridis_c()+
#  geom_rug(data=fig2e_data[order(-abs(beta.tt))][background==bg&type==exp_resp&FDRpop.raw<0.01&abs(betapop.raw)>0.2&FDR.tt<0.01&state==st&group=="loss"],mapping=aes(betapop.raw,color=abs(beta.tt)),sides="t",size=0.1)+
  xlab(expression("Raw population effect size (log"[2]*"FC"["r"]*")"))+
  ylab(expression("Adj. population effect size (log"[2]*"FC"["a"]*")"))+
  theme_yann()+theme(text=element_text(size=7))

pname=sprintf("%s/Fig2/Fig2d_scatter_%s_%s.pdf",FIG_DIR,exp_resp,st)
pdf(pname,width=5,height=5)
fig2e_plot+theme(text=element_text(size=10))
dev.off()

################################################################################
# Fig. 2E

exp_resp="DR"
ct_lin="lineage"

lineages=lineage_order
celltypes=celltype_order

if (exp_resp=="DE"&ct_lin=="lineage") {
  states=c("NS","COV","IAV")
  celltype_state=str_c(rep(lineages,each=3),rep(states,5),sep="_")
} else if (exp_resp=="DR"&ct_lin=="lineage") {
  states=c("COV","IAV")
  celltype_state=str_c(rep(lineages,each=2),rep(states,5),sep="_")
} else if (exp_resp=="DE"&ct_lin=="celltype") {
  states=c("NS","COV","IAV")
  celltype_state=str_c(rep(celltypes,each=3),rep(states,22),sep="_")
} else if (exp_resp=="DR"&ct_lin=="celltype") {
  states=c("COV","IAV")
  celltype_state=str_c(rep(celltypes,each=2),rep(states,22),sep="_")
}
attr(celltype_state,"names")=celltype_state

run_id='220409'

# popDEGs tested within each lineage
ANALYSE=sprintf("%s_condition_AFBEUB__%s_condition__%s",ct_lin,ct_lin,run_id)
popDE_noCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDE_noCellProp[,type:="expression"]
ANALYSE=sprintf("%s_condition_AFBEUB__logFC_%s_condition_logFC__%s",ct_lin,ct_lin,run_id)
popDR_noCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDR_noCellProp[,type:="response"]

# popDEGs tested within each lineage, accounting for finer cell proportions within the lineage considered
# (NB: MASH run on all lineages together)
ANALYSE=sprintf("%s_condition_AFBEUB__%s_condition__CellPropLineage_%s",ct_lin,ct_lin,run_id)
popDE_lineageCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDE_lineageCellProp[,type:="expression"]
ANALYSE=sprintf("%s_condition_AFBEUB__logFC_%s_condition_logFC__CellPropLineage_%s",ct_lin,ct_lin,run_id)
popDR_lineageCellProp=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
popDR_lineageCellProp[,type:="response"]

adjust=c("lineage","no")
expresp=c("DE","DR")

for (ad in adjust) {
  for (ex in expresp) {
    gsea_list=lapply(celltype_state,function(x){
      c=str_split(x,"_",simplify=T)[,1]
      s=str_split(x,"_",simplify=T)[,2]
      popde=get(sprintf("pop%s_%sCellProp",ex,ad))
      popde[,FDR:=ifelse(is.na(FDR),1,FDR)]
      popde[,effect_size:=beta]
      popde[,significance:=FDR]
      popde=popde[celltype==c&state==s]
      popde[order(-effect_size),setNames(effect_size,Symbol)]
    })
    assign(sprintf("gsea_list_%s_%s",ad,ex),gsea_list)
  }
}

hallmarks=gmtPathways("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/references/PATHWAYS/Human_GOBP_AllPathways_no_GO_iea_December_01_2021_symbol__only_GO.gmt")

pthwy="POSITIVE REGULATION OF CELL MIGRATION%GOBP%GO:0030335_NK"
targetPATHWAYS=list()
targetPATHWAYS[['NK']]=data.table(pathway=c('POSITIVE REGULATION OF CELL MIGRATION%GOBP%GO:0030335'),celltype=rep('NK',1),state=rep('COV',1))
  targetPATHWAYS=rbindlist(targetPATHWAYS)
pathway_list=paste(targetPATHWAYS$pathway,targetPATHWAYS$celltype,sep='_')

exp_resp="DR"

gsea_list_all_DE <- gsea_list_lineage_DE
gsea_list_all_DR <- gsea_list_lineage_DR

for (pthwy in pathway_list){
  cat(pthwy,exp_resp,'\n')
  ct=str_split(pthwy,"_",simplify=T)[,2]
  dir.create(sprintf("%s/Fig2/",FIG_DIR))
  pname=sprintf("%s/Fig2/gsea_%s_%s_%s.pdf",FIG_DIR,exp_resp,str_split(pthwy,"%",simplify=T)[,1],ct)
      plottitle=str_split(pthwy,"%",simplify=T)[,1]
      plottitle=paste0(substr(plottitle,1,1),tolower(substr(plottitle,2,nchar(plottitle))))
      plottitle=sprintf("%s (%s)",plottitle,ct)
      pathway_no_COV=hallmarks[[str_split(pthwy,"_",simplify=T)[,1]]]
      pathway_all_COV=hallmarks[[str_split(pthwy,"_",simplify=T)[,1]]]
      pathway_no_IAV=hallmarks[[str_split(pthwy,"_",simplify=T)[,1]]]
      pathway_all_IAV=hallmarks[[str_split(pthwy,"_",simplify=T)[,1]]]
  for (st in c("COV","IAV")) {
  ct_st=sprintf("%s_%s",ct,st)
    for (adj in c("no","all")){
      pathway_null_COV=hallmarks[[str_split(pthwy,"_",simplify=T)[,1]]]
      pathway_null_IAV=hallmarks[[str_split(pthwy,"_",simplify=T)[,1]]]
        assign(sprintf("stats_%s_%s",adj,st),get(sprintf("gsea_list_%s_%s",adj,exp_resp))[[ct_st]]*(-1))
        gseaParam=1
        ticksSize=0.7
        assign(sprintf("rnk_%s_%s",adj,st),rank(-get(sprintf("stats_%s_%s",adj,st))))
        assign(sprintf("ord_%s_%s",adj,st),order(get(sprintf("rnk_%s_%s",adj,st))))
        assign(sprintf("statsAdj_%s_%s",adj,st),get(sprintf("stats_%s_%s",adj,st))[get(sprintf("ord_%s_%s",adj,st))])
        assign(sprintf("statsAdj_%s_%s",adj,st),sign(get(sprintf("statsAdj_%s_%s",adj,st)))*(abs(get(sprintf("statsAdj_%s_%s",adj,st)))^gseaParam))
        assign(sprintf("statsAdj_%s_%s",adj,st),get(sprintf("statsAdj_%s_%s",adj,st))/max(abs(get(sprintf("statsAdj_%s_%s",adj,st)))))
        assign(sprintf("pathway_%s_%s",adj,st),unname(as.vector(na.omit(match(get(sprintf("pathway_%s_%s",adj,st)), names(get(sprintf("statsAdj_%s_%s",adj,st))))))))
        assign(sprintf("pathway_%s_%s",adj,st),sort(get(sprintf("pathway_%s_%s",adj,st))))
        assign(sprintf("gseaRes_%s_%s",adj,st),calcGseaStat(get(sprintf("statsAdj_%s_%s",adj,st)), selectedStats = get(sprintf("pathway_%s_%s",adj,st)),returnAllExtremes = TRUE))
        assign(sprintf("bottoms_%s_%s",adj,st),get(sprintf("gseaRes_%s_%s",adj,st))$bottoms)
        assign(sprintf("tops_%s_%s",adj,st),get(sprintf("gseaRes_%s_%s",adj,st))$tops)
        assign(sprintf("n_%s_%s",adj,st),length(get(sprintf("statsAdj_%s_%s",adj,st))))
        assign(sprintf("xs_%s_%s",adj,st),as.vector(rbind(get(sprintf("pathway_%s_%s",adj,st))-1,get(sprintf("pathway_%s_%s",adj,st)))))
        assign(sprintf("ys_%s_%s",adj,st),as.vector(rbind(get(sprintf("bottoms_%s_%s",adj,st)),get(sprintf("tops_%s_%s",adj,st)))))
        assign(sprintf("toPlot_%s_%s",adj,st),data.frame(x = c(0, get(sprintf("xs_%s_%s",adj,st)), get(sprintf("n_%s_%s",adj,st)) + 1), y = c(0, get(sprintf("ys_%s_%s",adj,st)), 0)))
        toPlot_iter=lapply(1:100,function(i){
        set.seed(i)
        assign(sprintf("stats_null_%s",st),setNames(get(sprintf("stats_%s_%s",adj,st)),sample(1:length(get(sprintf("stats_%s_%s",adj,st))),length(get(sprintf("stats_%s_%s",adj,st))),F)%>%{names(get(sprintf("stats_%s_%s",adj,st)))[.]}))
        assign(sprintf("rnk_null_%s",st),rank(-get(sprintf("stats_null_%s",st))))
        assign(sprintf("ord_null_%s",st),order(get(sprintf("rnk_null_%s",st))))
        assign(sprintf("statsAdj_null_%s",st),get(sprintf("stats_null_%s",st))[get(sprintf("ord_null_%s",st))])
        assign(sprintf("statsAdj_null_%s",st),sign(get(sprintf("statsAdj_null_%s",st)))*(abs(get(sprintf("statsAdj_null_%s",st)))^gseaParam))
        assign(sprintf("statsAdj_null_%s",st),get(sprintf("statsAdj_null_%s",st))/max(abs(get(sprintf("statsAdj_null_%s",st)))))
        assign(sprintf("pathway_null_%s",st),as.vector(na.omit(match(get(sprintf("pathway_null_%s",st)), names(get(sprintf("statsAdj_null_%s",st)))))))
        assign(sprintf("pathway_null_%s",st),sort(get(sprintf("pathway_null_%s",st))))
        assign(sprintf("gseaRes_null_%s",st),calcGseaStat(get(sprintf("statsAdj_null_%s",st)), selectedStats = get(sprintf("pathway_null_%s",st)),returnAllExtremes = TRUE))
        assign(sprintf("bottoms_null_%s",st),get(sprintf("gseaRes_null_%s",st))$bottoms)
        assign(sprintf("tops_null_%s",st),get(sprintf("gseaRes_null_%s",st))$tops)
        assign(sprintf("n_null_%s",st),length(get(sprintf("statsAdj_null_%s",st))))
        assign(sprintf("xs_null_%s",st),as.vector(rbind(get(sprintf("pathway_null_%s",st))-1,get(sprintf("pathway_null_%s",st)))))
        assign(sprintf("ys_null_%s",st),as.vector(rbind(get(sprintf("bottoms_null_%s",st)),get(sprintf("tops_null_%s",st)))))
        assign(sprintf("toPlot_null_%s",st),data.frame(x = c(0, get(sprintf("xs_null_%s",st)), get(sprintf("n_null_%s",st)) + 1), y = c(0, get(sprintf("ys_null_%s",st)), 0),iter=i))
        return(get(sprintf("toPlot_null_%s",st)))
        })%>%rbindlist()
  #      nbins=50
  #      increment=round(12672/nbins)
  #      breaks=seq(0,12672+increment,increment)
  #      toPlot_iter$bin=cut(toPlot_iter$x,breaks,include.lowest=T,right=F)
  #      toPlot_iter_IC=toPlot_iter[,quantile(y,c(0.025,0.975)),by=bin]
  #      toPlot_iter_IC$borne=rep(c("min","max"),nrow(toPlot_iter_IC)/2)
  #      toPlot_iter_IC=dcast(toPlot_iter_IC,bin~borne,value.var="V1")
  #      toPlot_iter_IC$x=toPlot_iter[,round(median(x)),by=bin]$V1
        assign(sprintf("toPlot_null_iter_%s",st),toPlot_iter)
    }
  }
      diff <- (max(tops_no_COV) - min(bottoms_no_COV))/16

  xmax=max(toPlot_null_iter_COV$x)
  toPlot_null_iter_COV_CI95=toPlot_null_iter_COV[,.(x=1:xmax,y=approx(x=x,y=y,xout=1:xmax)$y),by=iter][,.(lowerCI95=quantile(y,0.025),upperCI95=quantile(y,0.975)),by=x]
  #toPlot_null_iter_IAV_CI95=toPlot_null_iter_IAV[,.(x=1:xmax,y=approx(x=x,y=y,xout=1:xmax)$y),by=iter][,.(lowerCI95=quantile(y,0.025),upperCI95=quantile(y,0.975)),by=x]


  x = y = NULL
  g <- ggplot()+
  geom_ribbon(data=toPlot_null_iter_COV_CI95,mapping=aes(x=x,ymax=-lowerCI95,ymin=-upperCI95),alpha=0.25)+
  # geom_line(data=toPlot_null_iter_COV,mapping=aes(x = x, y = y,group=iter),color = color_conditions["COV"],alpha=0.25,size=0.1)+
  # geom_line(data=toPlot_null_iter_IAV,mapping=aes(x = x, y = y,group=iter),color = color_conditions["IAV"],alpha=0.25,size=0.1)+
  geom_line(data=toPlot_all_COV,mapping=aes(x = x, y = -y),color = color_conditions["COV"],alpha=0.75,linetype="dashed",size=0.3)+
  geom_line(data=toPlot_all_IAV,mapping=aes(x = x, y = -y),color = color_conditions["IAV"],alpha=0.75,linetype="dashed",size=0.3)+
  geom_line(data=toPlot_no_COV,mapping=aes(x = x, y = -y),color = color_conditions["COV"],size=0.3) +
  geom_line(data=toPlot_no_IAV,mapping=aes(x = x, y = -y),color = color_conditions["IAV"],size=0.3) + theme_yann() +
  geom_segment(data = data.frame(x = pathway_all_COV), mapping = aes(x = x,y = 0, xend = x, yend = diff/2), size = ticksSize,color=color_conditions["COV"],alpha=0.25)+
  geom_segment(data = data.frame(x = pathway_all_IAV), mapping = aes(x = x,y = -diff/2, xend = x, yend = 0), size = ticksSize,color=color_conditions["IAV"],alpha=0.25)+
  geom_segment(data = data.frame(x = pathway_no_COV), mapping = aes(x = x,y = 0, xend = x, yend = diff/2), size = ticksSize,color=color_conditions["COV"])+
  geom_segment(data = data.frame(x = pathway_no_IAV), mapping = aes(x = x,y = -diff/2, xend = x, yend = 0), size = ticksSize,color=color_conditions["IAV"])+
  geom_hline(yintercept = 0,colour = "black")+
  labs(x = expression(beta["AFB"]*">"*beta["EUB"]*"        "*beta["AFB"]*"<"*beta["EUB"]), y = "ES")+
  theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  ggtitle(plottitle)

  fig2f_plot <- g

  pdf(pname,width=6,height=5)
    print(g)
  dev.off()
}

################################################################################
# Fig. 2F

celltype_freq=lapply(lineage_order,function(X){
  if (X=="MONO") {
    celltypes=lineage_celltype[lineage==X&!grepl("INFECTED",celltype),celltype]
  } else {
    celltypes=lineage_celltype[lineage==X,celltype]
  }
  data=meta[POP!="ASH"&state=="NS",]
  lineage_sizes=data[lineage==X,.N,keyby=.(IID,lineage)]
  setnames(lineage_sizes,"N","TOT")
  iid_sizes=data[,.N,keyby=.(IID)]
  setnames(iid_sizes,"N","TOT_IID")
  celltype_sizes=data[lineage==X,.N,keyby=.(IID,celltype)]
  celltype_sizes=dcast(celltype_sizes,IID~celltype,value.var="N",fill=0)
  celltype_sizes=melt(celltype_sizes,measure.vars=celltypes,value.name="N",variable.name="celltype")
  celltype_freq=merge(celltype_sizes,lineage_sizes,by="IID")
  celltype_freq=merge(celltype_freq,iid_sizes,by="IID")
  celltype_freq[,FREQ_CT:=N/TOT]
  celltype_freq[,FREQ_ID:=N/TOT_IID]
  celltype_freq[,POP:=substr(IID,1,3)]
  celltype_freq[,lineage:=X]
  return(celltype_freq)
})%>%rbindlist()

celltype_freq$lineage=factor(celltype_freq$lineage,c("MONO","B","T.CD4","T.CD8","NK"))

cmv_status=fread(sprintf("%s/data/cmv_status.txt",DAT_POPDIFF_DIR))
setnames(cmv_status,c("CMV (y/n)","Standard Units"),c("CMV","standard_units"))
cmv_status=cmv_status[type=="sample",.(IID=IID_w_ASH,CMV)]
cmv_status$CMV=ifelse(cmv_status$CMV=="y","CMV","NO_CMV")

fig2f_data=merge(celltype_freq,cmv_status,by="IID",all.x=T)
fig2f_data=fig2f_data[!is.na(CMV),]
fig2f_data[,POP_CMV:=sprintf("%s_%s",POP,CMV)]

fig2f_data$celltype=factor(fig2f_data$celltype,celltype_order)

fig2f_data$CMV_int=ifelse(fig2f_data$CMV=="CMV",1,0)

fig2f_plot=ggplot(fig2f_data[POP_CMV!="AFB_NO_CMV"&celltype%in%c("T.CD8.EMRA","NK.M.LIKE"),])+
# option 1
  #geom_boxplot(aes(POP_CMV,FREQ_CT,fill=celltype),alpha=0.5,color=NA,outlier.shape=NA,notch=F)+
  #geom_boxplot(aes(POP_CMV,FREQ_CT,color=celltype),fill=NA,outlier.shape=NA,size=0.1,notch=F)+
# option 2
  #geom_boxplot(aes(POP_CMV,FREQ_CT,fill=celltype),alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  #geom_boxplot(aes(POP_CMV,FREQ_CT,color=celltype),fill=NA,outlier.shape=NA,size=0.1,notch=T)+
# option 3
  geom_violin(aes(POP_CMV,FREQ_CT,fill=celltype),alpha=0.5,color=NA,scale="width")+
  geom_violin(aes(POP_CMV,FREQ_CT,color=celltype),fill=NA,size=0.1,scale="width")+
  geom_boxplot(aes(POP_CMV,FREQ_CT),fill="white",color=NA,outlier.shape=NA,show.legend=F,alpha=0.5,notch=T)+
  geom_boxplot(aes(POP_CMV,FREQ_CT),fill=NA,outlier.shape=NA,show.legend=F,size=0.1,notch=T)+
#
  #geom_violin(aes(POP_CMV,FREQ_CT,fill=celltype),scale="width",alpha=0.7,color=NA)+
  #geom_violin(aes(POP_CMV,FREQ_CT,color=celltype),scale="width",draw_quantiles=0.5,fill=NA)+
  facet_grid(cols=vars(celltype),labeller=labeller(celltype=c("T.CD8.EMRA"="CD8+ EMRA T","NK.M.LIKE"="Memory-like NK")))+
  ylab("Fraction of immune lineage")+
  xlab("CMV serostatus")+
  scale_fill_manual(values=celltype_color,guide="none")+
  scale_y_continuous(limits=c(0,1),breaks=c(0,1),labels=c(0,1))+
  scale_color_manual(values=celltype_color,guide="none")+
  theme_plot()

  #pn="fig2f_freq" # option 1
  #pn="fig2f_freq_notch" # option 2
  pn="fig2f_freq_violin" # option 3
  pname=sprintf("%s/Fig2/%s.pdf",FIG_DIR,pn)
  pdf(pname,width=2,height=2)
  print(fig2f_plot)
  dev.off()

pname=sprintf("%s/Fig2/Fig2f_boxplot.pdf",FIG_DIR)
pdf(pname,width=6,height=5)
fig2f_plot+theme(text=element_text(size=10))
dev.off()
