################################################################################
################################################################################
# File name: 4b8__clues_analysis.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Analyse CLUES results
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./4b__clues_analysis__lib.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))
source(sprintf("./load_eQTLs.R",MISC_DIR))
source(sprintf("./querySNPs.R",MISC_DIR))

################################################################################
# Calculate selection metric (Z-score)

# load (r)eQTL allele frequency trajectories calculated by CLUES per population
# and summarized in 4b7__summarize_clues.py
# smooth allele frequency trajectories using loess regression
# calculate first derivative (dotallele), adjust for frequency differences
# (dotallele_norm) and compute Z-score (z)

clues <- lapply(c("YRI","CEU","CHS"),function(p){
  data <- lapply(c("eQTL","reQTL"),function(expresp){
  data=fread(sprintf("%s/%s_CutOff2000/SNPs/%s/allele_freq_max__all_SNPs.tsv.gz",DAT_CLUES_DIR,p,expresp))
  colnames(data)[1:5]=c("epoch","allele_pp_max","rsID","lowerCI","upperCI")
  data[,allele_pp_max_smooth:=predict(loess(allele_pp_max~epoch,span=0.1),epoch),by=rsID]
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

################################################################################
# Add natural selection metrics

eqtl_type='eQTL'
filter='Signif'
eqtl_file_select=fread(sprintf("%s/3_eQTL_mapping/%s_%s_both_PBS_nSL_long_v2.tsv.gz",DATA_DIR,eqtl_type,filter))
eqtl_file_select[,type:='eQTL']

eqtl_type='reQTL'
filter='Signif'
reqtl_file_select=fread(sprintf("%s/3_eQTL_mapping/%s_%s_both_PBS_nSL_long_v2.tsv.gz",DATA_DIR,eqtl_type,filter))
reqtl_file_select[,type:='reQTL']

all_selected=rbind(eqtl_file_select,reqtl_file_select)

clues$PBS=all_selected[METRIC=="PBS",setNames(VALUE,sprintf("%s_%s",snps,POP))][clues[,sprintf("%s_%s",rsID,pop)]]
clues$P_PBS=all_selected[METRIC=="PBS",setNames(P,sprintf("%s_%s",snps,POP))][clues[,sprintf("%s_%s",rsID,pop)]]
clues[,P_PBS_CAT:=case_when(P_PBS<0.05~"0.05",P_PBS<0.01~"0.01",T~"weak")]

SNP_info=getMap(annotate=T)

clues$gene_name=eQTL_Signif_both[,setNames(gene_name,snps)][clues$rsID]

# write results
fwrite(clues,sprintf("%s/selection_zscore.tsv.gz",DAT_CLUES_DIR),sep="\t")

################################################################################
# Compute enrichment for SARS-CoV-2 (r)eQTLs

eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EIP)
CIS_DIST_TEXT='100kb'
RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_220409"
RUN_REQTL_LINEAGE="lineage_condition_logFC__logFC__CellPropLineage_SVs_220409"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_220409"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__logFC__CellPropLineage_SVs_220409"

# load resampling SNP sets
snpSets=fread(sprintf('%s/%s/dist_%s/All_eQTL_snpsSets.txt.gz',eQTL_DIR,RUN_EQTL_LINEAGE,CIS_DIST_TEXT))
allowed_celltypes=paste(c(lineage_5,celltype_22),collapse='|')
regex=sprintf('^(r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\3',set)]
snpSets[,state:=gsub(regex,'\\4',set)]
snpSets[,specificity:=gsub(regex,'\\5',set)]
snpSets[,type:=ifelse(specificity!='','reQTL_breakdown',type)]


selection_date_reqtlcov=clues[(type=="reQTL"&rsID%in%snpSets[set=='reQTL_COV'&snps%nin%snpSets[set=="reQTL_COV_specific",unique(snps)],unique(snps)] & !is.na(z_smooth)),.(start=max(c(-1,epoch[abs(z_smooth)>3])),end=min(c(2001,epoch[abs(z_smooth)>3])),max_abs_z=max(abs(z_smooth))),by=.(rsID,gene_name,type,pop)]
selection_date_reqtlcov[,selected:=ifelse(max_abs_z>3,T,F)]
selection_date_reqtlcov[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date_reqtlcov[,COV_specific:=F]

selection_date_reqtlcovsp=clues[(type=="reQTL"&rsID%in%snpSets[set=="reQTL_COV_specific",unique(snps)] & !is.na(z_smooth)),.(start=max(c(-1,epoch[abs(z_smooth)>3])),end=min(c(2001,epoch[abs(z_smooth)>3])),max_abs_z=max(abs(z_smooth))),by=.(rsID,gene_name,type,pop)]
selection_date_reqtlcovsp[,selected:=ifelse(max_abs_z>3,T,F)]
selection_date_reqtlcovsp[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date_reqtlcovsp[,COV_specific:=T]

selection_date_reqtlcov=rbind(selection_date_reqtlcovsp,selection_date_reqtlcov)

DT <- selection_date_reqtlcov[COV_specific==T,]
PctCHS=PctCEU=PctYRI=c()
for(b in 1:1000){
PctCHS[b]=DT[pop=='CHS',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctYRI[b]=DT[pop=='YRI',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctCEU[b]=DT[pop=='CEU',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
}

quantile(PctCHS,c(0.025,0.975))*100
quantile(PctYRI,c(0.025,0.975))*100
quantile(PctCEU,c(0.025,0.975))*100

selection_date$PBS=clues[,setNames(PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$P_PBS=clues[,setNames(P_PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
#selection_date$PBS_idx <- selection_date[order(-PBS),setNames(1:nrow(selection_date),sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$z_max <- clues[!is.na(z_smooth),max(abs(z_smooth)),by=.(rsID,type,pop)][,setNames(V1,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]

# write results
fwrite(selection_date,sprintf("%s/selection_date.tsv.gz",DAT_CLUES_DIR),sep="\t")

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


################################################################################


cluesFreq=fread(sprintf('%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/clues_trajectories.tsv.gz',EIP))

snpSets=fread(sprintf('%s/users/Javier/data/snp_sets/snpSets_June2022.txt',EIP))
COV_reQTL=snpSets[set=="reQTL_COV",unique(snps)]
COV_reQTL_GW=snpSets[set=="reQTL_GWsignif_COV",unique(snps)]





sum(COV_reQTL%chin% cluesFreq[,unique(rsID)])
# 1034

cluesFreq[,sum(COV_reQTL%chin% unique(rsID)),by=pop]
cluesFreq[,sum(COV_reQTL_GW%chin% unique(rsID)),by=pop]

myset='reQTL_COV'
myset='reQTL_COV_specific'
selection_date_reqtlcov=cluesFreq[(type=="reQTL"&rsID%in%snpSets[set==myset,unique(snps)] & !is.na(z_smooth)),.(start=max(c(-1,epoch[abs(z_smooth)>3])),end=min(c(2001,epoch[abs(z_smooth)>3])),max_abs_z=max(abs(z_smooth)),Selection_effect=sum(z_smooth[abs(z_smooth)>3])),by=.(rsID,gene_name,type,pop)]
selection_date_reqtlcov[,selected:=ifelse(max_abs_z>3,T,F)]
library(dplyr)
selection_date_reqtlcov[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
  selection_date_reqtlcov[,length(unique(rsID)),by=pop]






selection_date=cluesFreq[!is.na(z_smooth),.(start=max(c(-1,epoch[abs(z_smooth)>3])),
                                            end=min(c(2001,epoch[abs(z_smooth)>3])),
                                            max_abs_z=max(abs(z_smooth)),
                                            Selection_effect=sum(z_smooth[abs(z_smooth)>3]),
                                            Selection_effect_weak=mean(z_smooth),
                                            start_weak=max(c(-1,epoch[abs(z_smooth)>2])),
                                            end_weak=min(c(2001,epoch[abs(z_smooth)>2])))
                                            ,by=.(rsID,gene_name,type,pop)]
selection_date[,selected:=ifelse(max_abs_z>3,T,F)]
selection_date[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date[,length(unique(rsID)),by=.(pop,type,reQTL_COV_specific=(rsID%in%snpSets[set==myset,unique(snps)]),window)]

fwrite(selection_date,sprintf('%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/selection_date.tsv',EIP))

eQTL_Stats_celltype_mash[,beta_lower:=sign(beta)*(pmax(abs(beta)-2*se,0))]
reQTL_Stats_celltype_mash[,beta_lower:=sign(beta)*(pmax(abs(beta)-2*se,0))]
eQTL_best_celltype=eQTL_Stats_celltype_mash[order(snps,gene_name,-abs(beta_lower)),head(.SD[,.(best_celltype=celltype,best_state=state,beta_lower,pvalue_best=pvalue,type='eQTL')],1),by=.(snps,gene_name)]
reQTL_best_celltype=reQTL_Stats_celltype_mash[order(snps,gene_name,-abs(beta_lower)),head(.SD[,.(best_celltype=celltype,best_state=state,beta_lower,pvalue_best=pvalue,type='reQTL')],1),by=.(snps,gene_name)]

eQTL_Stats_celltype_mash[order(snps,gene_name,-abs(beta_lower)),head(.SD[,.(snps,gene_name,celltype,state,beta_lower,pvalue)],1),by=.(snps,gene_name)]

all_selected=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig5/test_PBS/geneList_all_selected_v2.tsv",EIP))
all_selected=merge(all_selected,selection_date[,.(snps=rsID,gene_name,POP=pop,start,end,max_abs_z,Selection_effect,Selection_effect_weak,start_weak,end_weak,selected)],by=c('snps','gene_name','POP'),all.x=TRUE)
all_selected=merge(all_selected,rbind(eQTL_best_celltype,reQTL_best_celltype),by=c('snps','gene_name','type'),all.x=TRUE)
selected_alleles=all_selected[METRIC=='PBS',.(snps,gene_name,type,POP,celltype,state,pvalue,best_celltype,best_state,beta_lower,
                                P_PBS=P,max_abs_z,Selection_effect_weak,start_weak,end_weak,
                                selected,start,end,Selection_effect,
                                GO_immune,ISG_effector,IEI,Targeted_viruses,
                                REF,ALT,ALT_DERANC,DAF_or_MAF_CEU, DAF_or_MAF_CHS, DAF_or_MAF_YRI,minP_perGene)]

selected_allele=selected_alleles[order(minP_perGene,-abs(Selection_effect_weak),pvalue),][!duplicated(paste(POP,snps,gene_name,type)),]

fwrite(selected_alleles,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V9/suppTables/TableS7/TableS7B_1_selection_dates_all_genes.tsv",EIP))
 dcast(selected_allele[!duplicated(snps,pop),],snps~POP,value.var='Selection_effect_weak',fill=0)

selected_allele[,sum(selected==TRUE,na.rm=T),keyby=.(POP,cut(end,seq(0,2200,by=200)))]

popDEGs=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V8/suppTables/TableS4/TableS4A_popDEG.tsv",EIP))
Signif_popDE=popDEGs[FDR_lineage_adj<.01 & abs(beta_lineage_adj)>.2,.(gene_name=Symbol,celltype,state)]

mediation=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V8/suppTables/TableS6/TableS6_mediation_analyses.tsv",EIP))
setnames(mediation,'mediator','rsID')
merge(mediation[mediator_type=="genetics",],selection_date[pop!='CHS',],by='rsID',all.y=TRUE)


fisher.test(as.matrix(dcast(tab,pop~frac_var)[2:1,-1])

tab=merge(mediation[mediator_type=="genetics" & !is.na(frac_var),],selection_date[pop!='CHS',],by='rsID')[,sum(abs(Selection_effect_weak)>1),by=.(pop,frac_var>0.5 & abs(betapop)>0.2)]
dcast(tab,frac_var~pop)
fisher.test(as.matrix(dcast(tab,pop~frac_var)[2:1,-1]))

#### comparison of the frequency of adaptive events between African and Europeans at popDEGs with mediation.
mediated_annot=merge(mediation[mediator_type=="genetics",],selection_date[pop!='CHS',],by='rsID')[frac_var>.5,][!duplicated(paste(rsID,pop)),]
mediated_annot[,mean(max_abs_z>3),by=pop]
# pop        V1
# 1: YRI 0.2091097
# 2: CEU 0.3405640
mediated_annot[,fisher.test(table(max_abs_z>3,pop))]
# Fisher's Exact Test for Count Data
# p-value = 7.765e-06
# odds ratio
#  0.5123159
# 95 percent confidence interval:
#  0.3779538 0.6923835

selection_date[rsID=='rs1142888',]
# rsID gene_name type pop start  end max_abs_z Selection_effect Selection_effect_weak start_weak end_weak selected
# 1: rs1142888      GBP7 eQTL YRI    -1 2001  1.320060           0.0000             0.5381256         -1     2001    FALSE
# 2: rs1142888      GBP7 eQTL CEU  1272  782  4.327977         638.5713             1.2034810       1315      753     TRUE
# 3: rs1142888      GBP7 eQTL CHS    -1 2001  2.904957           0.0000             1.1528575       1742      307    FALSE

mediated_annot[,mean(max_abs_z>4),by=pop]
# mediated_annot[,mean(max_abs_z>4),by=pop]
#   pop         V1
# 1: YRI 0.04140787
# 2: CEU 0.12147505

mediated_annot[,fisher.test(table(max_abs_z>4,pop))]


gene_list_strongPBS_coloc=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig5/test_PBS/gene_list_strongPBS_coloc.tsv",EIP),sep='\t')
merge(gene_list_strongPBS_coloc,selection_date,by.x=c('snps','POP'),by.y=c('rsID','pop'))
merge(gene_list_strongPBS_coloc,selection_date,by.x=c('snps'),by.y=c('rsID')
merge(gene_list_strongPBS_coloc[!duplicated(snps),],selection_date,by.x=c('snps'),by.y=c('rsID'))[abs(Selection_effect_weak)>1,.(snps,gene_name.x,POP,pop,Selection_effect_weak,max_abs_z,start,end)]
#          snps gene_name.x POP pop Selection_effect_weak max_abs_z start  end
# 1:  rs1559828         DR1 YRI YRI             -1.559947  2.535289    -1 2001
# 2:  rs1559828         DR1 YRI CEU              1.111553  2.717018    -1 2001
# 3:  rs1559828         DR1 YRI CHS              1.359449  3.086143  1825 1786
# 4:  rs2326569       RAB2A YRI CEU              1.159144  3.233121  1499 1440
# 5:  rs2625446       RAB2A YRI CHS              1.005959  3.103639   838  823
# 6:  rs4504381     C5orf56 CEU YRI             -1.808433  2.712711    -1 2001
# 7:   rs547915        LMNA YRI YRI              1.233701  2.542269    -1 2001
# 8:   rs547915        LMNA YRI CEU              1.630810  4.291111  1998 1774
# 9:   rs569414         DR1 YRI YRI             -1.401012  2.280138    -1 2001
# 10:   rs569414         DR1 YRI CEU              1.343146  2.978060    -1 2001
# 11:   rs569414         DR1 YRI CHS              1.262676  3.526592  1076  989
# 12: rs61542988      SNHG26 CEU CHS              1.669754  3.776341  1372  971
# 13:  rs7532549       TMED5 YRI CEU              1.755308  4.096023  1697 1559
# 14:  rs7532549       TMED5 YRI CHS              1.195595  2.906381    -1 2001


merge(selected_allele[GO_immune=='yes' & POP=='CEU' & abs(max_abs_z)>3,],Signif_popDE,by=c('gene_name','celltype','state'))

all_selected[METRIC=='PBS'& colocalized=='yes' & P<.01 & (covid_C2_pval<.01 | covid_B2_pval<.01 | covid_A2_pval<.01),][!duplicated(snps)]

# Selected_loci_window_reQTL_COV=merge(selection_date_reqtlcov[selected==TRUE & window=='start_window' & pop=='CHS',],unique(all_selected[POP=='CHS' & P<.05,.(P,METRIC,VALUE,DAF_or_MAF_CEU,DAF_or_MAF_CHS,DAF_or_MAF_YRI,snps,covid_A2_pval)]),by.x='rsID',by.y='snps')
# Selected_loci_window_reQTL_COV[,DeltaDAF:=DAF_or_MAF_CHS-(DAF_or_MAF_CEU+DAF_or_MAF_YRI)/2]

all_selected[selected==TRUE & POP=='CEU' & METRIC=='PBS',][order(P,pvalue),][!duplicated(paste(snps,gene_name,type)),][1:3,]
all_selected[selected==TRUE & POP=='CEU' & METRIC=='PBS',][P<.01,][order(-abs(beta_lower)),][!duplicated(paste(snps,gene_name,type)),][1:3,]

all_selected[METRIC=='PBS',][order(-abs(beta_lower)),][!duplicated(paste(POP,snps,gene_name,type)),][ISG_effector=='yes',.(snps,gene_name,celltype,best_celltype,state,best_state,beta_lower,pvalue,pvalue_best,max_abs_z,start,end,P,POP,Targeted_viruses,EffectOn,type,POP,Selection_effect_weak,start_weak,end_weak)][gene_name=='IFITM3',]
all_selected[METRIC=='PBS',][order(-abs(beta_lower)),][!duplicated(paste(POP,snps,gene_name,type)),][METRIC=='PBS' & P<.01,][order(-abs(Selection_effect_weak)),][1:10,.(snps,gene_name,celltype,best_celltype,state,best_state,beta_lower,pvalue,pvalue_best,max_abs_z,start,end,P,POP,Targeted_viruses,EffectOn,type,POP,Selection_effect_weak,start_weak,end_weak)][gene_name=='IFITM3',]
#all_selected[,beta_lower:=sign(beta)*(pmax(abs(beta)-2*se,0))]
all_selected[selected==TRUE & POP=='CEU' & METRIC=='PBS',][order(P),][1:3,]
dcast(selection_date[,length(unique(rsID)),by=.(pop,type,reQTL_COV_specific=(rsID%in%snpSets[set==myset,unique(snps)]),window)],pop+type+reQTL_COV_specific~window)
dcast(selection_date[,length(unique(rsID)),by=.(pop,type,reQTL_COV_specific=type,window)],pop+type+reQTL_COV_specific~window)

Selected_loci_window=merge(selection_date[selected==TRUE & window=='start_window' & pop=='CHS',],unique(all_selected[POP=='CHS' & P<.05,.(P,METRIC,VALUE,DAF_or_MAF_CEU,DAF_or_MAF_CHS,DAF_or_MAF_YRI,snps,covid_A2_pval,GO_immune,IEI,COV_VIPs)]),by.x='rsID',by.y='snps')

selection_date[rsID=='rs1142888',]
