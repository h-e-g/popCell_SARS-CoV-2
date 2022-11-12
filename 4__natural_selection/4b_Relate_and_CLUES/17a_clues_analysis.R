.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(dynamicTreeCut))
suppressMessages(library(BiocNeighbors))
suppressMessages(library(CelliD))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(ggrastr))
suppressMessages(library(scales))
suppressMessages(library(data.table))
suppressMessages(library(DropletUtils))
suppressMessages(library(SoupX))
suppressMessages(library(batchelor))
suppressMessages(library(kBET))
suppressMessages(library(HGC))
suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

theme_set(theme_bw())
theme_update(
  text=element_text(family="serif",size=12),
  panel.grid=element_blank(),legend.position="bottom",
  strip.background=element_rect(fill="#012158"),strip.text=element_text(color="white")
)

`%nin%`=Negate(`%in%`)

source("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R")

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"
DAT_DIR=sprintf("%s/single_cell/project/pop_eQTL/data",EIP)
PAP_DIR=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V9",EIP)
FIG_DIR=sprintf("%s/figureMaterial",PAP_DIR)
RES_DIR=sprintf("%s/single_cell/resources",EIP)

source(sprintf("%s/template_scripts/processing_pipeline/00_set_colors.R",RES_DIR))
source(sprintf("%s/template_scripts/querySNPs.R",RES_DIR))
source(sprintf("%s/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",RES_DIR))

CLUES_DIR=sprintf("%s/4_natural_selection/clues",DAT_DIR)
p<- "CHS"
clues <- lapply(c("YRI","CEU","CHS"),function(p){
  data <- lapply(c("eQTL","reQTL"),function(expresp){
  data=fread(sprintf("%s/%s_CutOff2000/SNPs/%s/allele_freq_max__all_SNPs.tsv.gz",CLUES_DIR,p,expresp))
  colnames(data)[1:5]=c("epoch","allele_pp_max","rsID","lowerCI","upperCI")
  data[,allele_pp_max_smooth:=predict(loess(allele_pp_max~epoch,span=0.1),epoch),by=rsID]
#  data[,allele_pp_max_smooth:=ksmooth(y=allele_pp_max,x=epoch,kernel="normal",bandwidth=0.1),by=rsID]
  data[order(rsID,-epoch),dotallele:=c(diff(allele_pp_max),0),by=.(rsID)]
  data[,dotallele_norm:=dotallele/sqrt(allele_pp_max*(1-allele_pp_max))]
  data[,z:=dotallele_norm/sd(dotallele_norm),by=epoch]
  data[order(rsID,-epoch),dotallele_smooth:=c(diff(allele_pp_max_smooth),0),by=.(rsID)]
  data[,dotallele_norm_smooth:=dotallele_smooth/sqrt(allele_pp_max_smooth*(1-allele_pp_max_smooth))]
  data[,z_smooth:=dotallele_norm_smooth/sd(dotallele_norm_smooth),by=epoch]
  data[,type:=expresp]
  data[,pop:=p]
  data$gene_name=get(sprintf("%s_Signif_both",expresp))[,setNames(gene_name,snps)][data$rsID]
  return(data)
  })%>%rbindlist()
})%>%rbindlist()

#clues_chs_eqtl=fread(sprintf("%s/4_natural_selection/clues/CHS_CutOff2000/SNPs/eQTL/allele_freq_max__all_SNPs.tsv.gz",DAT_DIR),col.names=c("epoch","allele_pp_max","rsID","lowerCI","upperCI"))
#clues_chs_eqtl[,type:="eQTL"]
#clues_chs_eqtl[,pop:="CHS"]
#clues_chs_eqtl$gene_name=eQTL_Signif_both[,setNames(gene_name,snps)][clues_chs_eqtl$rsID]
#
#clues_chs_reqtl=fread(sprintf("%s/4_natural_selection/clues/CHS_CutOff2000/SNPs/reQTL/allele_freq_max__all_SNPs.tsv.gz",DAT_DIR),col.names=c("epoch","allele_pp_max","rsID","lowerCI","upperCI"))
#clues_chs_reqtl[,type:="reQTL"]
#clues_chs_reqtl[,pop:="CHS"]
#clues_chs_reqtl$gene_name=reQTL_Signif_both[,setNames(gene_name,snps)][clues_chs_reqtl$rsID]
#
#clues_chs <- rbind(clues_chs_eqtl,clues_chs_reqtl)

color_populations_1kg <- color_populations
attr(color_populations_1kg,"names")=c("YRI","CEU","CHS")

rsid_types=clues[,unique(sprintf("%s_%s",rsID,type))]

i <- 1
for (rsid_type in rsid_types[8147:11120]) {
  rsid=gsub("_.+$","",rsid_type)
  expresp=gsub("^.+_","",rsid_type)
  print(sprintf("%s out of 11120.",i))
  n_genes <- length(clues[rsID==rsid,unique(gene_name)])
  if (n_genes>1){
    for (g in clues[rsID==rsid,unique(gene_name)]){
      pname <- sprintf("%s/testsFigure/CLUES/trajectories/%s_%s_%s_smooth.png",PAP_DIR,g,rsid,expresp)
      p <- ggplot(clues[rsID==rsid&type==expresp&gene_name==g,])+
      geom_vline(xintercept=c(770,970),linetype="dashed")+
      geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
      geom_point(aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.2)+
      scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
      ggtitle(sprintf("%s (%s, %s, %s)",g,rsid,expresp,clues[rsID==rsid,unique(pop)]))+
      xlab("Generations from present")+ylab("Expected allele frequency")+
      scale_y_continuous(limits=c(0,1))+theme(legend.position="none")
    png(pname,res=480,width=8,height=6,units='in')
    print(p)
    dev.off()
    i <- i+1
    }
  } else{
  pname <- sprintf("%s/testsFigure/CLUES/trajectories/%s_%s_%s_smooth.png",PAP_DIR,clues[rsID==rsid,unique(gene_name)],rsid,expresp)
    p <- ggplot(clues[rsID==rsid&type==expresp,])+
    geom_vline(xintercept=c(770,970),linetype="dashed")+
    geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
    geom_point(aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.2)+
    scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
    ggtitle(sprintf("%s (%s, %s, %s)",clues[rsID==rsid,unique(gene_name)],rsid,expresp,clues[rsID==rsid,unique(pop)]))+
    xlab("Generations from present")+ylab("Expected allele frequency")+
    scale_y_continuous(limits=c(0,1))+theme(legend.position="none")
  png(pname,res=480,width=8,height=6,units='in')
  print(p)
  dev.off()
  i <- i+1
  }
}

rsid="rs1028396"
expresp="reQTL"
  pname <- sprintf("%s/testsFigure/CLUES/trajectories/%s_%s_%s_zscore.png",PAP_DIR,clues[rsID==rsid,unique(gene_name)],rsid,expresp)
    p <- ggplot(clues[rsID==rsid&type==expresp,])+
    geom_vline(xintercept=c(770,970),linetype="dashed")+
#    geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
    geom_point(aes(epoch,z_smooth,color=pop),alpha=1,size=0.2)+
    scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
    ggtitle(sprintf("%s (%s, %s, %s)",clues[rsID==rsid,unique(gene_name)],rsid,expresp,clues[rsID==rsid,unique(pop)]))+
    xlab("Generations from present")+ylab("Z score")+
#    scale_y_continuous(limits=c(0,1))+
theme(legend.position="none")
  png(pname,res=480,width=8,height=6,units='in')
  print(p)
  dev.off()

clues_test <- clues[rsID=="rs431420",]
clues_test[,smooth_freq:=predict(loess(allele_pp_max~epoch,span=0.2),epoch),by=pop]
clues_test[,smooth_lci:=predict(loess(lowerCI~epoch),epoch,span=0.2),by=pop]
clues_test[,smooth_uci:=predict(loess(upperCI~epoch),epoch,span=0.2),by=pop]
  pname <- sprintf("%s/testsFigure/CLUES/trajectories/%s_%s_%s_loess.png",PAP_DIR,clues_test[,unique(gene_name)],"rs431420","eQTL")
    p <- ggplot(clues_test)+
    geom_vline(xintercept=c(770,970),linetype="dashed")+
    geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
    geom_point(aes(epoch,smooth_freq,color=pop),alpha=1,size=0.2)+
    scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
    ggtitle(sprintf("%s (%s, %s, %s)",clues_test[,unique(gene_name)],"rs431420","eQTL",clues_test[,unique(pop)]))+
    xlab("Generations from present")+ylab("Expected allele frequency")+
    scale_y_continuous(limits=c(0,1))+theme(legend.position="none")
  png(pname,res=480,width=8,height=6,units='in')
  print(p)
  dev.off()

n_epochs=clues[,nrow(.SD),by=.(rsID,type,pop)][,unique(V1)]

tic('Computing derivatives')
clues[order(type,rsID,-epoch),dotallele:=c(diff(allele_pp_max),0),by=.(rsID,type,pop)]
clues[,dotallele_norm:=dotallele/sqrt(allele_pp_max*(1-allele_pp_max))]
clues[,z:=dotallele_norm/sd(dotallele_norm),by=epoch]
toc()

#tic('Computing derivatives')
#dotallele=lapply(1:(n_epochs-1),function(i){
# clues_chs_reqtl[,(allele_pp_max[i+1]-allele_pp_max[i])/(epoch[i]-epoch[i+1]),by=.(rsID)][,.(epoch=i,dotallele=V1,rsID)]
#})%>%rbindlist()%>%.[order(rsID,epoch),]
#toc()
#
#dotallele <- rbind(dotallele,dotallele[,list(unique(rsID,type))][,.(epoch=0,dotallele=0,rsID=V1,type)])[order(rsID,epoch),]
#
#clues_chs$dotallele=dotallele[,setNames(dotallele,sprintf("%s_%s",rsID,epoch))][clues_chs[,sprintf("%s_%s",rsID,epoch)]]
#
#tic('Computing second derivatives')
#dotdotallele=lapply(1:(n_epochs-1),function(i){
#clues_chs[,(dotallele[i+1]-dotallele[i])/(epoch[i]-epoch[i+1]),by=rsID][,.(epoch=i,dotdotallele=V1,rsID)]
#})%>%rbindlist()%>%.[order(rsID,epoch),]
#toc()
#
#dotdotallele <- rbind(dotdotallele,dotdotallele[,list(unique(rsID))][,.(epoch=0,dotdotallele=0,rsID=V1)])[order(rsID,epoch),]
#
#clues_chs$dotdotallele=dotdotallele[,setNames(dotdotallele,sprintf("%s_%s",rsID,epoch))][clues_chs[,sprintf("%s_%s",rsID,epoch)]]
#
#clues_chs_sub=clues_chs[rsID%in%c("rs664910","rs786433"),]
#clues_chs_sub[,q25:=quantile(dotallele,0.25),by=rsID]
#clues_chs_sub[,q75:=quantile(dotallele,0.75),by=rsID]
#clues_chs_sub[,iqr:=q75-q25,by=rsID]
#
#clues_chs[,sum(dotallele)/sd(dotallele),by=rsID][order(-abs(V1)),head(.SD,10)][,rsID]
#
#for (snp in c("rs664910","rs786433")) {
#
#model=lm(epoch~dotallele,clues_chs[rsID==snp,])
#q25=clues_chs[rsID==snp,quantile(dotallele,0.25)]
#q75=clues_chs[rsID==snp,quantile(dotallele,0.75)]
#
#beta0=model$coefficients["(Intercept)"]
#beta1=model$coefficients["dotallele"]
#
#q25epoch=round(beta0+beta1*q25)
#q75epoch=round(beta0+beta1*q75)
#
#}

snps_to_plot=clues_chs[!is.nan(z),max(abs(z)),by=.(rsID,type)][order(-abs(V1)),head(rsID,20)]
pname=sprintf("%s/testsFigure/CLUES/%s_%s.%s",PAP_DIR,"top20","CHS","png")
p <- ggplot(clues_chs[rsID%in%snps_to_plot,])+
geom_point(aes(epoch,dotallele_norm,color=gene_name),alpha=0.2,size=0.05)+
geom_vline(xintercept=c(770,970),linetype="dashed")+
xlab("Generations from present")+ylab("Expected allele frequency")+
scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
png(pname,res=480,width=8,height=6,units='in')
print(p)
dev.off()


for (p in c("YRI","CHS","CEU")){
snps_to_plot=clues[!is.nan(z_smooth)&pop==p,max(abs(z_smooth)),by=.(rsID,type)][order(-abs(V1)),head(rsID,20)]
pname=sprintf("%s/testsFigure/CLUES/%s_%s_loess.%s",PAP_DIR,"top20",p,"png")
p <- ggplot(clues[pop==p&rsID%in%snps_to_plot,])+
geom_point(aes(epoch,dotallele_norm_smooth,color=gene_name),alpha=0.2,size=0.05)+
geom_vline(xintercept=c(770,970),linetype="dashed")+
#facet_grid(cols=vars(pop))+
xlab("Generations from present")+ylab("Expected allele frequency")+
scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
png(pname,res=480,width=8,height=6,units='in')
print(p)
dev.off()
}

################################################################################

eqtl_type='eQTL'
filter='Signif'
eqtl_file_select=fread(sprintf("%s/3_eQTL_mapping/%s_%s_both_PBS_nSL_long_v2.tsv.gz",DATA_DIR,eqtl_type,filter))
eqtl_file_select[,type:='eQTL']

eqtl_type='reQTL'
filter='Signif'
reqtl_file_select=fread(sprintf("%s/3_eQTL_mapping/%s_%s_both_PBS_nSL_long_v2.tsv.gz",DATA_DIR,eqtl_type,filter))
reqtl_file_select[,type:='reQTL']

all_selected=rbind(eqtl_file_select,reqtl_file_select)

snps_of_interest=clues[epoch!=0,max(z_smooth),by=rsID][order(-abs(V1)),head(rsID,20)]

all_selected[snps%in%snps_of_interest,]

clues$PBS=all_selected[METRIC=="PBS",setNames(VALUE,sprintf("%s_%s",snps,POP))][clues[,sprintf("%s_%s",rsID,pop)]]
clues$P_PBS=all_selected[METRIC=="PBS",setNames(P,sprintf("%s_%s",snps,POP))][clues[,sprintf("%s_%s",rsID,pop)]]
clues[,P_PBS_CAT:=case_when(P_PBS<0.05~"0.05",P_PBS<0.01~"0.01",T~"weak")]
#clues_chs$P_PBS=all_selected[POP=="CHS"&METRIC=="PBS",setNames(P,snps)][clues_chs$rsID]
#clues_chs[,P_PBS_CAT:=case_when(P_PBS<0.05~"0.05",P_PBS<0.01~"0.01",T~"weak")]

clues_chs$ABS_NSL=all_selected[POP=="CHS"&METRIC=="ABS_NSL",setNames(VALUE,snps)][clues_chs$rsID]
clues_chs$P_ABS_NSL=all_selected[POP=="CHS"&METRIC=="ABS_NSL",setNames(P,snps)][clues_chs$rsID]

SNP_info=getMap(annotate=T)

clues_chs$gene_name=eQTL_Signif_both[,setNames(gene_name,snps)][clues_chs$rsID]

CLUES_DIR=sprintf("%s/4_natural_selection/clues",DAT_DIR)
fwrite(clues_chs,sprintf("%s/clues_trajectory__info.tsv.gz",CLUES_DIR),sep="\t")

################################################################################

eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EIP)
CIS_DIST_TEXT='100kb'
RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_220409"
RUN_REQTL_LINEAGE="lineage_condition_logFC__logFC__CellPropLineage_SVs_220409"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_220409"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__logFC__CellPropLineage_SVs_220409"

snpSets=fread(sprintf('%s/%s/dist_%s/All_eQTL_snpsSets.txt.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT))
allowed_celltypes=paste(c(lineage_5,celltype_22),collapse='|')
regex=sprintf('^(r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\3',set)]
snpSets[,state:=gsub(regex,'\\4',set)]
snpSets[,specificity:=gsub(regex,'\\5',set)]
snpSets[,type:=ifelse(specificity!='','reQTL_breakdown',type)]

sets_of_interest=sprintf("%s_%s",rep(c("eQTL","reQTL"),each=2),rep(c("COV","COV_specific"),2))

################################################################################

selection_date_reqtlcov=clues[(type=="reQTL"&rsID%in%snpSets[set=='reQTL_COV'&snps%nin%snpSets[set=="reQTL_COV_specific",unique(snps)],unique(snps)] & !is.na(z_smooth)),.(start=max(c(-1,epoch[abs(z_smooth)>3])),end=min(c(2001,epoch[abs(z_smooth)>3])),max_abs_z=max(abs(z_smooth))),by=.(rsID,gene_name,type,pop)]
selection_date_reqtlcov[,selected:=ifelse(max_abs_z>3,T,F)]
selection_date_reqtlcov[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date_reqtlcov[,COV_specific:=F]


################################################################################
# are CHS selection events for COV-specific reQTLs more likely to have started
# in Enard's window?

myset="reQTL_COV_specific"

selection_date_reqtlcovsp=clues[(type=="reQTL"&rsID%in%snpSets[set==myset,unique(snps)] & !is.na(z_smooth)),.(start=max(c(-1,epoch[abs(z_smooth)>3])),end=min(c(2001,epoch[abs(z_smooth)>3])),max_abs_z=max(abs(z_smooth))),by=.(rsID,gene_name,type,pop)]
selection_date_reqtlcovsp[,selected:=ifelse(max_abs_z>3,T,F)]
selection_date_reqtlcovsp[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date_reqtlcovsp[,COV_specific:=T]

selection_date_reqtlcov=rbind(selection_date_reqtlcovsp,selection_date_reqtlcov)

fisher.test(selection_date_reqtlcovsp[,pop=="CHS"],selection_date_reqtlcovsp[,window=="start_window"])

#	Fisher's Exact Test for Count Data
#
#data:  selection_date_reqtlcovsp[, pop == "CHS"] and selection_date_reqtlcovsp[, window == "start_window"]
#p-value = 0.0007292
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.447778 4.843938
#sample estimates:
#odds ratio
#  2.635947

DT <- selection_date_reqtlcovsp
PctCHS=PctCEU=PctYRI=c()
for(b in 1:1000){
PctCHS[b]=DT[pop=='CHS',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctYRI[b]=DT[pop=='YRI',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctCEU[b]=DT[pop=='CEU',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
}

quantile(PctCHS,c(0.025,0.975))*100
quantile(PctYRI,c(0.025,0.975))*100
quantile(PctCEU,c(0.025,0.975))*100

#selection_date_reqtlcovsp=clues[(type=="reQTL"&abs(z_smooth)>3)|(type=="reQTL"&rsID%in%snpSets[set==myset,unique(snps)]),.(start=max(epoch),end=min(epoch),max_abs_z=max(z_smooth)),by=.(rsID,gene_name,type,pop)]
#selection_date_reqtlcovsp[,max_abs_z:=ifelse(is.nan(max_abs_z),0,max_abs_z)]
#selection_date_reqtlcovsp[,selected:=ifelse(max_abs_z>3,T,F)]
#selection_date_reqtlcovsp[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]

fisher.test(selection_date_reqtlcovsp[,pop=="CHS"],selection_date_reqtlcovsp[,window=="start_window"])

selection_date$PBS=clues[,setNames(PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$P_PBS=clues[,setNames(P_PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
#selection_date$PBS_idx <- selection_date[order(-PBS),setNames(1:nrow(selection_date),sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$z_max <- clues[!is.na(z_smooth),max(abs(z_smooth)),by=.(rsID,type,pop)][,setNames(V1,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]

################################################################################

selection_date=clues[abs(z_smooth)>3,.(start=max(epoch),end=min(epoch)),by=.(rsID,gene_name,type,pop)]
selection_date$PBS=clues[,setNames(PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$P_PBS=clues[,setNames(P_PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
#selection_date$PBS_idx <- selection_date[order(-PBS),setNames(1:nrow(selection_date),sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$z_max <- clues[!is.na(z_smooth),max(abs(z_smooth)),by=.(rsID,type,pop)][,setNames(V1,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]

selection_date[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]

pname=sprintf("%s/testsFigure/CLUES/%s_%s_%s.%s",PAP_DIR,"sel_date__scatter","CHS","all","png")
png(pname,res=480,width=8,height=6,units='in')
p <- ggplot(selection_date[pop=="CHS"])+
geom_point(aes(z_max,PBS,color=type),size=0.7)+
geom_text_repel(aes(z_max,PBS,color=type,label=gene_name))+
facet_grid(rows=vars(pop),scales="free_y")+
xlab("Max. abs. Z score")+ylab("PBS")+
ggtitle("eQTLs and reQTLs")+
#scale_color_viridis_c()
scale_color_viridis_d()
#scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
print(p)
dev.off()

myset="reQTL_COV_specific"
plot_data=selection_date[type=="reQTL"&rsID%in%snpSets[set==myset,unique(snps)],]
plot_data$idx <- plot_data[order(pop,-start),setNames(1:nrow(plot_data),sprintf("%s%s%s",rsID,type,pop))][plot_data[,sprintf("%s%s%s",rsID,type,pop)]]

sum_sig_z=lapply(2:4,function(thr){
  data=clues[,sum(abs(z)>thr),by=.(pop,epoch)]
  data[,threshold:=thr]
  return(data)
})%>%rbindlist()

pname=sprintf("%s/testsFigure/CLUES/%s.%s",PAP_DIR,"z_threshold","png")
png(pname,res=480,width=8,height=6,units='in')
ggplot(sum_sig_z[!is.na(V1),])+
geom_line(aes(epoch,V1,color=as.character(threshold)),size=0.5)+
geom_vline(xintercept=c(770,970),linetype="dashed")+
facet_grid(rows=vars(pop),scales="free_y")+
ylab("Number of Z-scores above threshold")+xlab("Selection time frame")+
#ggtitle(myset)+
ggtitle("No filter")+
theme(panel.spacing=unit(0,"mm"))+
scale_color_viridis_d()
#scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
dev.off()

pname=sprintf("%s/testsFigure/CLUES/%s_%s.%s",PAP_DIR,"sel_date__z_smooth",myset,"png")
png(pname,res=480,width=8,height=6,units='in')
ggplot(plot_data[order(PBS)])+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=PBS),size=0.5)+
geom_vline(xintercept=c(770,970),linetype="dashed")+
facet_grid(rows=vars(pop),scales="free_y")+
ylab("SNPs")+xlab("Selection time frame")+
ggtitle(myset)+
theme(panel.spacing=unit(0,"mm"))+
scale_color_viridis_c()
#scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
dev.off()

for (mypop in c("YRI","CHS","CEU")){
pname=sprintf("%s/testsFigure/CLUES/%s_%s_%s.%s",PAP_DIR,"sel_date__z_smooth",mypop,myset,"png")
png(pname,res=480,width=8,height=6,units='in')
p <- ggplot(plot_data[order(PBS),][pop==mypop,])+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=(start<=970&start>=770)),size=0.5)+
geom_vline(xintercept=c(770,970),linetype="dashed")+
facet_grid(rows=vars(pop),scales="free_y")+
ylab("SNPs")+xlab("Selection time frame")+
ggtitle(myset)+
#scale_color_viridis_c()
scale_color_viridis_d()
#scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
print(p)
dev.off()
}

#> fisher.test(plot_data[pop=="CHS",selected=="PBSSig"],plot_data[pop=="CHS",start>970])
#
#        Fisher's Exact Test for Count Data
#
#data:  plot_data[pop == "CHS", selected == "PBSSig"] and plot_data[pop == "CHS", start > 970]
#p-value = 0.01893
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.080956 3.520881
#sample estimates:
#odds ratio
#  1.940617

plot_data[,type:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]

for (mypop in c("YRI","CHS","CEU")){
scatterplot_data <- plot_data[pop==mypop,]
scatterplot_data[,type:=case_when(start<=970&start>=770~"start_window",start>970&end<770~"overlap_window",T~"not_window")]
pname=sprintf("%s/testsFigure/CLUES/%s_%s_%s.%s",PAP_DIR,"sel_date__scatter",mypop,myset,"png")
png(pname,res=480,width=8,height=6,units='in')
p <- ggplot(scatterplot_data)+
geom_point(aes(z_max,PBS,color=type),size=1)+
geom_text_repel(aes(z_max,PBS,color=type,label=gene_name))+
facet_grid(rows=vars(pop),scales="free_y")+
ylab("Max. abs. Z score")+xlab("PBS")+
ggtitle(myset)+
#scale_color_viridis_c()
scale_color_viridis_d()
#scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
print(p)
dev.off()
}

pname=sprintf("%s/testsFigure/CLUES/%s_%s_%s.%s",PAP_DIR,"sel_date__lilrb1",myset,pop,ext)
png(pname,res=480,width=8,height=6,units='in')
ggplot(plot_data)+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=gene_name=="LILRB1"),size=0.5)+
geom_vline(xintercept=c(770,970),linetype="dashed")+
ylab("SNPs")+xlab("Selection time frame")+
ggtitle(myset)#+
#scale_color_viridis_c()
#scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
dev.off()

myset="eQTL_COV"
expresp=gsub("_.+$","",myset)
plot_data=clues_chs[type==expresp&rsID%in%snpSets[set==myset,unique(snps)],]
plot_data$idx <- plot_data[order(-PBS),setNames(1:nrow(plot_data),sprintf("%s%s%s",rsID,type,pop))][plot_data[,sprintf("%s%s%s",rsID,type,pop)]]
plot_data[,z_scale:=scale(z,center=T,scale=T)]
plot_data[!is.na(z_scale),z_scale_alpha:=(z_scale+abs(min(z_scale)))/max(z_scale+abs(min(z_scale)))]

#mat=dcast(plot_data[!is.na(z)],rsID+type+pop~epoch,value.var="z",fill=0)
#matnames=sprintf("%s%s%s",mat$rsID,mat$type,mat$pop)
#mat=as.matrix(mat[,4:2001])
#rownames(mat)=matnames
#test=hclust(dist(mat))

plot_data$cluster_order=setNames(test$order,test$labels)[plot_data[,sprintf("%s%s%s",rsID,type,pop)]]

pname=sprintf("%s/testsFigure/CLUES/%s_%s_%s.%s",PAP_DIR,"sel_date_hmp",myset,pop,ext)
png(pname,res=480,width=8,height=6,units='in')
ggplot(plot_data[!is.na(z_scale)&type==expresp,])+
#geom_tile(aes(x=epoch,y=idx,color=z,alpha=z_scale_alpha))+
geom_tile(aes(x=epoch,y=idx,color=z))+
geom_vline(xintercept=c(770,970),linetype="dashed")+
ylab("SNPs")+xlab("Selection time frame")+
ggtitle(myset)+
scale_color_gradient2(low="blue",high="red",mid="white")
#scale_color_viridis_c()
#scale_color_discrete(guide=guide_legend(override.aes=list(size=1,alpha=1),title=""))
dev.off()




pop="CHS"
ext="png"

i <- 4

pname=sprintf("%s/testsFigure/CLUES/%s_%s.%s",PAP_DIR,sets_of_interest[i],pop,ext)
png(pname,res=480,width=8,height=6,units='in')
ggplot(clues_chs[order(PBS),][rsID%in%snpSets[set==sets_of_interest[i],unique(snps)],])+
#geom_point(aes(epoch,dotallele_norm,color=PBS),alpha=0.2,size=0.05)+
geom_point(aes(epoch,dotallele_norm,color=P_PBS_CAT),alpha=0.3,size=0.05)+
#scale_color_viridis_c()+
geom_vline(xintercept=c(800,1000,1200),linetype="dashed")
dev.off()

hiZscore_snps=clues_chs[!is.nan(z),max(abs(z)),by=.(rsID,type,gene_name)][abs(V1)>2,unique(rsID)]

#pname=sprintf("%s/testsFigure/CLUES/%s_%s.%s",PAP_DIR,"all_eQTLs",pop,ext)
pname=sprintf("%s/testsFigure/CLUES/%s_%s.%s",PAP_DIR,"hiZscore_zscore",pop,ext)
png(pname,res=480,width=8,height=6,units='in')
ggplot(clues_chs[rsID%in%hiZscore_snps,])+
geom_line(aes(epoch,z,color=PBS,group=rsID),alpha=0.2,size=0.4)+
scale_color_viridis_c()+
geom_vline(xintercept=c(800,1000,1200),linetype="dashed")
dev.off()

################################################################################

clues_lilrb1 <- lapply(c("CHS"),function(p){
  data <- lapply(c("reQTL"),function(expresp){
  data=fread(sprintf("%s/%s_CutOff2000/SNPs/%s/allele_freq_max__rs4806787.tsv.gz",CLUES_DIR,p,expresp))
  colnames(data)[1:5]=c("epoch","allele_pp_max","rsID","lowerCI","upperCI")
  data[,allele_pp_max_smooth:=predict(loess(allele_pp_max~epoch,span=0.1),epoch),by=rsID]
#  data[,allele_pp_max_smooth:=ksmooth(y=allele_pp_max,x=epoch,kernel="normal",bandwidth=0.1),by=rsID]
  data[order(rsID,-epoch),dotallele:=c(diff(allele_pp_max),0),by=.(rsID)]
  data[,dotallele_norm:=dotallele/sqrt(allele_pp_max*(1-allele_pp_max))]
  data[,z:=dotallele_norm/sd(dotallele_norm),by=epoch]
  data[order(rsID,-epoch),dotallele_smooth:=c(diff(allele_pp_max_smooth),0),by=.(rsID)]
  data[,dotallele_norm_smooth:=dotallele_smooth/sqrt(allele_pp_max_smooth*(1-allele_pp_max_smooth))]
  data[,z_smooth:=dotallele_norm_smooth/sd(dotallele_norm_smooth),by=epoch]
  data[,type:=expresp]
  data[,pop:=p]
  data$gene_name=get(sprintf("%s_Signif_both",expresp))[,setNames(gene_name,snps)][data$rsID]
  return(data)
  })%>%rbindlist()
})%>%rbindlist()

rsid="rs4806787"
DT <- clues_lilrb1

pname <- sprintf("%s/testsFigure/CLUES/trajectories_selbins/%s_%s_%s_smooth.png",PAP_DIR,DT[rsID==rsid,unique(gene_name)],rsid,expresp)
  p <- ggplot(DT[rsID==rsid&type==expresp,])+
  geom_vline(xintercept=c(770,970),linetype="dashed")+
  geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
  geom_point(aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.2)+
  scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
  ggtitle(sprintf("%s (%s, %s, %s)",DT[rsID==rsid,unique(gene_name)],rsid,expresp,DT[rsID==rsid,unique(pop)]))+
  xlab("Generations from present")+ylab("Expected allele frequency")+
  scale_y_continuous(limits=c(0,1))+theme(legend.position="none")
png(pname,res=480,width=8,height=6,units='in')
print(p)
dev.off()
