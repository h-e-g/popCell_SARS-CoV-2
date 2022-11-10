################################################################################
################################################################################
# File name: FigS7.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Extended Data Figure 7
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="LIBRARY"
source(sprintf("./1a__quality_control__lib.R",LIB_DIR))

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

theme_set(theme_bw())
theme_update(
  text=element_text(family="sans",size=7),
  panel.grid=element_blank(),legend.position="bottom",
  strip.background=element_rect(fill="#012158"),strip.text=element_text(color="white")
)

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"
#EIP="/home/yaaquino/evo_immuno_pop"
DAT_DIR=sprintf("%s/single_cell/project/pop_eQTL/data",EIP)
CLUES_DIR=sprintf("%s/4_natural_selection/clues",DAT_DIR)
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EIP))
FIG_DIR=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V9/figureMaterial",EIP)

celltype_color=c(
"T.CD4.N"="#169FD8","T.CD4.E"="#005274","T.Reg"="#03C2C0","T.gd"="#ffce73",
"T.CD8.N"="#00BB54","T.CD8.CM.EM"="#0EAD20","T.CD8.EMRA"="#048c41","MAIT"="#004A21",
"NK.M.LIKE"="#63C425","NK.CD56dim"="#B7DC2A","NK.CD56brt"="#7D9E00","ILC"="#3F4F00",
"B.N.K"="#9281DE","B.N.L"="#a681de","B.M.K"="#6045d9","B.M.L"="#8045d9",
"Plasmablast"="#6801A4","MONO.CD14"="#B7637E","MONO.CD16"="#8F568A","pDC"="#8E6B00",
"cDC"="#E7b315","MONO.CD14.INFECTED"="#B9364B")

lineage_order=c("MONO","B","T.CD4","T.CD8","NK")

celltype_order=c("MONO.CD14","MONO.CD16","MONO.CD14.INFECTED","cDC","pDC","B.N.K","B.N.L","B.M.K","B.M.L","Plasmablast","T.CD4.N","T.CD4.E","T.Reg","T.gd","T.CD8.N","T.CD8.CM.EM","T.CD8.EMRA","ILC","MAIT","NK.CD56dim","NK.CD56brt","NK.M.LIKE")
celltype_label=c("MONO CD14+","MONO CD16+","MONO IAV+","cDC","pDC","B N k","B N l","B M k","B M l","Plasmablast","T CD4+ N","T CD4+ E","T Reg","T gd","T CD8+ N","T CD8+ CM/EM","T CD8+ EMRA","ILC","MAIT","NK CD56dim","NK CD56brt","NK mem")
celltype_label_exp=expression("CD14+ monocyte","CD16+ monocyte","IAV+ monocyte","cDC","pDC",kappa*"-LC naive B",lambda*"-LC naive B",kappa*"-LC memory B",lambda*"-LC memory B","Plasmablast","CD4+ naive T","CD4+ effector T","Regulatory T",gamma*delta*" T","CD8+ naive T","CD8+ CM/EM T","CD8+ EMRA T","ILC","MAIT","CD56"["dim"]*" NK","CD56"["brt"]*" NK","Memory-like NK")
names(celltype_label_exp)=celltype_order

################################################################################

eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EIP)
CIS_DIST_TEXT='100kb'
RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_220409"
RUN_REQTL_LINEAGE="lineage_condition_logFC__logFC__CellPropLineage_SVs_220409"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_220409"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__logFC__CellPropLineage_SVs_220409"

snpSets=fread(sprintf('%s/%s/dist_%s/All_eQTL_snpsSets.txt.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT))
allowed_celltypes=paste(c(lineage_order,celltype_order),collapse='|')
regex=sprintf('^(r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\3',set)]
snpSets[,state:=gsub(regex,'\\4',set)]
snpSets[,specificity:=gsub(regex,'\\5',set)]
snpSets[,type:=ifelse(specificity!='','reQTL_breakdown',type)]

sets_of_interest=sprintf("%s_%s",rep(c("eQTL","reQTL"),each=2),rep(c("COV","COV_specific"),2))

################################################################################
# Fig. S8a
clues=fread(sprintf("%s/clues_trajectories.tsv.gz",CLUES_DIR))

selection_date=clues[abs(z_smooth)>3,.(start=max(epoch),end=min(epoch)),by=.(rsID,gene_name,type,pop)]
selection_date$PBS=clues[,setNames(PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$P_PBS=clues[,setNames(P_PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$z_max <- clues[!is.na(z_smooth),max(abs(z_smooth)),by=.(rsID,type,pop)][,setNames(V1,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date[,set:=ifelse(rsID%in%snpSets[set=="reQTL_COV",unique(snps)],'reQTLCOV','other')]
selection_date[,selected:=ifelse(P_PBS<0.05,"PBSSig",'PBSNSig')]

figs8ab_data=selection_date[(set=="reQTLCOV")&type=="reQTL",]
figs8ab_data$idx <- figs8ab_data[order(pop,-start),setNames(1:nrow(figs8ab_data),sprintf("%s%s%s",rsID,type,pop))][figs8ab_data[,sprintf("%s%s%s",rsID,type,pop)]]

mypop="YRI"
breaklims=figs8ab_data[pop==mypop,c(min(idx),max(idx))]
figs8a_plot <- ggplot(figs8ab_data[order(PBS),][pop==mypop,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=0,ymax=800),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#008000",alpha=0.2)+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=selected),size=0.25)+
scale_color_manual(values=c("PBSSig"="#008000","PBSNSig"="#888888"),breaks="PBSSig",labels="Top 5% PBS",guide=guide_legend(title=""))+
scale_y_continuous(breaks=breaklims,labels=c("-3","3"))+
facet_grid(rows=vars(pop))+
coord_cartesian(ylim=breaklims)+
theme(legend.position=c(0.15,0.32),axis.ticks.y=element_line(color=NA),axis.text.y=element_text(color=NA))+
ylab("SNPs")+xlab("Generations before present")


################################################################################
# Fig. S8b

mypop="CEU"
breaklims=figs8ab_data[pop==mypop,c(min(idx),max(idx))]
figs8b_plot <- ggplot(figs8ab_data[order(PBS),][pop==mypop,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-100,ymax=800),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#eba206",alpha=0.2)+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=selected),size=0.25)+
scale_color_manual(values=c("PBSSig"="#eba206","PBSNSig"="#888888"),breaks="PBSSig",labels="Top 5% PBS",guide=guide_legend(title=""))+
scale_y_continuous(breaks=breaklims,labels=c("-3","3"))+
facet_grid(rows=vars(pop))+
coord_cartesian(ylim=breaklims)+
theme(legend.position=c(0.15,0.32),axis.ticks.y=element_line(color=NA),axis.text.y=element_text(color=NA))+
ylab("SNPs")+xlab("Generations before present")

################################################################################
# Fig. S8c

figs8ab_data[,snp_of_interest:=ifelse(rsID%in%int_snps,'yes','no')]

mypop="CHS"
breaklims=figs8ab_data[pop==mypop,c(min(idx),max(idx))]
figs8c_plot <- ggplot(figs8ab_data[order(PBS),][pop==mypop,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-100,ymax=800),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=snp_of_interest),size=0.25)+
facet_grid(rows=vars(pop))+
geom_segment(data=data.table(x=c(721,1203),y=c(411,320),yend=c(600,600)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(951,1634),y=c(411,320),yend=c(600,600)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
scale_color_manual(values=c("yes"="#71458d","no"="#888888"),breaks="yes",labels="SNP of interest",guide=guide_legend(title=""))+
scale_y_continuous(breaks=breaklims,labels=c("-3","3"))+
coord_cartesian(ylim=breaklims)+
theme(legend.position=c(0.15,0.32),axis.ticks.y=element_line(color=NA),axis.text.y=element_text(color=NA))+
ylab("SNPs")+xlab("Generations before present")

################################################################################
# Fig. S8d

int_snps=c("rs4806787","rs1028396")
rand_snps=clues[pop=="CHS"&rsID%nin%int_snps,sample(unique(rsID),100)]

figs8d_data=clues[pop=="CHS"&rsID%in%c(rand_snps,int_snps),]
figs8d_data[,pop:=factor(pop,rev(c('YRI','CEU','CHS')))]
figs8d_data[,random:=ifelse(rsID%chin%int_snps,F,T)]

figs8d_plot=ggplot()+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
  geom_line(data=figs8d_data[random==T,],mapping=aes(epoch,allele_pp_max_smooth,group=rsID,color=pop),alpha=0.2,color="gray")+
  geom_line(data=figs8d_data[random==F,],mapping=aes(epoch,allele_pp_max_smooth,group=rsID,color=pop))+
  scale_color_manual(values=color_populations_1kg)+
  scale_y_continuous(breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"))+
  xlab("Generations from present")+ylab("Frequency (f)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(pop))+
  theme(panel.spacing=unit(0,'pt'))

################################################################################
# Fig. S8e

figs8e_plot=ggplot()+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
  geom_line(data=figs8d_data[random==T,],mapping=aes(epoch,dotallele_norm_smooth*1000,group=rsID,color=pop),alpha=0.2,color="gray")+
  geom_line(data=figs8d_data[random==F,],mapping=aes(epoch,dotallele_norm_smooth*1000,group=rsID,color=pop))+
  scale_color_manual(values=color_populations_1kg)+
coord_cartesian(ylim=c(-3.5,4.5))+
  xlab("Generations from present")+ylab("d/dt(f)*1000")+
  facet_grid(rows=vars(pop))+
  theme(panel.spacing=unit(0,'pt'))

################################################################################
# Fig. S8f

figs8f_plot=ggplot()+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-20,ymax=20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
  geom_line(data=figs8d_data[random==T,],mapping=aes(epoch,z_smooth,color=pop,group=rsID),alpha=0.2,color="gray")+
  geom_line(data=figs8d_data[random==F,],mapping=aes(epoch,z_smooth,color=pop,group=rsID))+
  geom_hline(yintercept=3,linetype='dashed',size=0.1)+
  geom_hline(yintercept=-3,linetype='dashed',size=0.1)+
geom_segment(data=data.table(x=c(721,1203),y=c(3,3),yend=c(-20,-20)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(951,1634),y=c(3,3),yend=c(-20,-20)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
  scale_color_manual(values=color_populations_1kg)+
  scale_y_continuous(breaks=c(-3,0,3))+
coord_cartesian(ylim=c(-5,10))+
  xlab("Generations from present")+ylab("Z-score")+
  facet_grid(rows=vars(pop))+
  theme(panel.spacing=unit(0,'pt'))

################################################################################
pname=sprintf("%s/FigS8/FigS8.pdf",FIG_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      ggplotGrob(figs8a_plot+theme(legend.position="none")),
      ggplotGrob(figs8b_plot+theme(legend.position="none")),
      ggplotGrob(figs8c_plot+theme(legend.position="none")),
      ggplotGrob(figs8d_plot+theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())),
      ggplotGrob(figs8e_plot+theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())),
      ggplotGrob(figs8f_plot+theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())),
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(4,4,1,1),
      c(5,5,1,1),
      c(6,6,1,1),
      c(3,3,2,2),
      c(3,3,2,2),
      c(3,3,2,2),
      c(7,7,7,7)
    ), heights=c(1.5,1.5,1.5,1.5,1.5,1.5,1), widths=c(1,3.5,1.7,2.8)
  )
dev.off()

################################################################################
# Fig. S7b
meta=fread(sprintf("%s/data/metadata.tsv",DAT_DES_DIR))
meta[,POP:=substr(IID,1,3)]
setnames(meta,"condition","state")

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

figs7b_data=celltype_freq[celltype%in%celltype_focus,]
figs7b_data[,celltype:=factor(celltype,celltype_focus)]

figs7b_plot=ggplot(figs7b_data)+
# option 1
  #geom_boxplot(aes(celltype,FREQ_CT,fill=celltype,group=POP),color=NA,outlier.shape=NA,show.legend=F,alpha=0.5,notch=F)+
  #geom_boxplot(aes(celltype,FREQ_CT,color=celltype,group=POP),fill=NA,outlier.shape=NA,show.legend=F,size=0.1,notch=F)+
# option 2
  #geom_boxplot(aes(celltype,FREQ_CT,fill=celltype,group=POP),color=NA,outlier.shape=NA,show.legend=F,alpha=0.5,notch=T)+
  #geom_boxplot(aes(celltype,FREQ_CT,color=celltype,group=POP),fill=NA,outlier.shape=NA,show.legend=F,size=0.1,notch=T)+
# option 3
  geom_violin(aes(celltype,FREQ_CT,fill=celltype,group=POP),alpha=0.5,color=NA,scale="width")+
  geom_violin(aes(celltype,FREQ_CT,color=celltype,group=POP),fill=NA,size=0.1,scale="width")+
  geom_boxplot(aes(celltype,FREQ_CT,group=interaction(celltype,POP)),fill="white",color=NA,outlier.shape=NA,show.legend=F,alpha=0.5,notch=T)+
  geom_boxplot(aes(celltype,FREQ_CT,group=interaction(celltype,POP)),fill=NA,outlier.shape=NA,show.legend=F,size=0.1,notch=T)+
  geom_point(aes(celltype,FREQ_CT,fill=celltype),alpha=0)+
  ylab("Fraction of immune lineage")+
  xlab("Cell type")+
  scale_y_continuous(breaks=c(0,100),labels=c(0,1))+
scale_x_discrete(labels=c("CD16+ monocyte","Memory B","CD4+ effector T","CD8+ EMRA T","Memory-like NK"))+
  scale_fill_manual(values=celltype_color,guide="none")+
  scale_color_manual(values=celltype_color,guide="none")+
  facet_grid(cols=vars(lineage),labeller=labeller(lineage=c("MONO"="MYELOID","B"="B","T.CD4"="CD4+ T","T.CD8"="CD8+ T","NK"="NK")),scales="free_x")+
  theme_plot()+
  theme(text=element_text(size=6),panel.spacing=unit(0,"pt"))

  #pn="figs7b_freq" # option 1
  #pn="figs7b_freq_notch" # option 2
  pn="figs7b_freq_violin" # option 3
  pname=sprintf("%s/FigS7/%s.pdf",FIG_DIR,pn)
  pdf(pname,width=4,height=2)
  print(figs7b_plot)
  dev.off()
