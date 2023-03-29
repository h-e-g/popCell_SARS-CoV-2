################################################################################
################################################################################
# File name: Fig4.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Figure 4
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

eQTL_DIR="3__eQTL_mapping/data"
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

################################################################################
# Fig. 4a

# load PBS score resamplings
resamples_PBS=paste0('resampling_newq99_30_06_2022/',dir(sprintf('%s/users/Javier/results/resampling_newq99_30_06_2022',EVO_IMMUNO_POP_ZEUS)))

# extract relevant information
resamp_results=list()
for (i in resamples_PBS){
  cat(i,'\n')
  resamp_results[[i]]=fread(sprintf('%s/users/Javier/results/%s',EVO_IMMUNO_POP_ZEUS,i))
  resamp_results[[i]]=resamp_results[[i]][,.(RESAMP_NUM_SEL=mean(RESAMP_NUM_SEL),
      RESAMP_MEAN_SCORE=mean(RESAMP_MEAN_SCORE),
      FE_NUM_SEL=mean((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL)),
      lowerCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.025),
      upperCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.975),
      FE_MEAN_SCORE=mean(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,na.rm=T),
      lowerCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.025,na.rm=T),
      upperCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.975,na.rm=T)),
      by=.(OBS_NUM_SEL,NUM_SEL_PVAL,OBS_MEAN_SCORE,MEAN_SCORE_PVAL,NUM_TARGET_EQTLS,NUM_RESAMP)]
}

resamp_results=rbindlist(resamp_results,idcol='snp_set_resamp')
REGEX=sprintf(".*/(CEU|YRI|CHS)_(ABSnSL|PBS)_((r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?)(_NumTest[0-9]+)?_NumResamps10000.txt",allowed_celltypes)
resamp_results[,POP:=gsub(REGEX,'\\1',snp_set_resamp)]
resamp_results[,stat:=gsub(REGEX,'\\2',snp_set_resamp)]
resamp_results[,set:=gsub(REGEX,'\\3',snp_set_resamp)]
resamp_results[,num:=gsub('_NumTest','',gsub(REGEX,'\\9',snp_set_resamp))]

resamp_results[,FDR_MEAN_SCORE:=p.adjust(MEAN_SCORE_PVAL,'fdr'),by=stat]
resamp_results[,FDR_NUM_SEL:=p.adjust(NUM_SEL_PVAL,'fdr'),by=stat]

resamp_results=merge(resamp_results,unique(snpSets[,.(set, type, celltype, state, specificity)]),by='set')

fig4a_data=resamp_results[grepl('eQTL_(COV|IAV|NS|shared)',set) & specificity!='stronger' & stat=="PBS",.(point=FE_MEAN_SCORE,
    lowci=lowerCI_MEAN_SCORE,
    highci=upperCI_MEAN_SCORE,
    pval=MEAN_SCORE_PVAL),
    by=.(POP,type,celltype,state,specificity)]

fig4a_data[,FDR:=p.adjust(pval,'fdr')]
fig4a_data[state=='NS',state:='b']
fig4a_data[,state2:=paste(celltype,state,specificity)]
fig4a_data[,POP:=factor(POP,c('YRI','CEU','CHS'))]
fig4a_data[,state2:=case_when(state2==" COV specific"~"COV-specific",state2==" IAV specific"~"IAV-specific",state2=="  shared"~"Shared",T~state2)]

fig4a_plot <- ggplot(fig4a_data,aes(x = state2,y=point,color=POP,fill=POP)) +
  geom_errorbar(aes(ymin=lowci, ymax=highci),width=0.2,position=position_dodge(width = 0.6),size=0.2) +
  geom_point(aes(fill=POP),position=position_dodge(width = 0.6),alpha=0.2,pch=21) +
  geom_point(aes(color=POP),position=position_dodge(width = 0.6),fill=NA,pch=21) +
  facet_grid(cols=vars(type),scales="free",labeller=labeller(type=c("eQTL"="eQTL","reQTL"="reQTL","reQTL_breakdown"="reQTL breakdown"))) +
  scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg) +
  geom_hline(aes(yintercept=1), colour="black",linetype="dotted") +
  xlab("")+ylab("Fold enrichment")+
  labs(color="Population")+
  theme_yann() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),panel.spacing=unit(0,'pt'))

################################################################################
# Fig. 4b

CLUES_DIR="4__natural_selection/data/CLUES"
clues=fread(sprintf("%s/clues_trajectories.tsv.gz",CLUES_DIR))

selection_date=clues[abs(z_smooth)>3,.(start=max(epoch),end=min(epoch)),by=.(rsID,gene_name,type,pop)]
selection_date$PBS=clues[,setNames(PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$P_PBS=clues[,setNames(P_PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$z_max <- clues[!is.na(z_smooth),max(abs(z_smooth)),by=.(rsID,type,pop)][,setNames(V1,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date[,set:=ifelse(rsID%in%snpSets[set=="reQTL_COV",unique(snps)],'reQTLCOV','other')]
selection_date[,selected:=ifelse(P_PBS<0.05,"PBSSig",'PBSNSig')]

fig4b_data=selection_date[(set=="reQTLCOV")&type=="reQTL",]
fig4b_data$idx <- fig4b_data[order(pop,-start),setNames(1:nrow(fig4b_data),sprintf("%s%s%s",rsID,type,pop))][fig4b_data[,sprintf("%s%s%s",rsID,type,pop)]]

mypop="CHS"
fig4b1_plot <- ggplot(fig4b_data[order(PBS),][pop==mypop,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=200,ymax=800),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=selected),size=0.2)+
geom_segment(data=data.table(x=c(848,1203),y=c(328,320),yend=c(800,800)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(848,1203),y=c(328,320),yend=c(0,0)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(1384,1634),y=c(328,320),yend=c(800,800)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(1384,1634),y=c(328,320),yend=c(0,0)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_text_repel(data=fig4b_data[pop==mypop&type=="reQTL"&rsID%in%snps_to_plot],aes(x=end,y=idx,label=gene_name,color=selected),seed=1,xlim=0,direction="x",segment.size=0.1,segment.alpha=0.5,show.legend=F,size=1.5,segment.linetype='dashed')+
scale_color_manual(values=c("PBSSig"="#71458d","PBSNSig"="#AAAAAA"),breaks="PBSSig",labels="Top 5% PBS",guide=guide_legend(title=""))+
scale_y_continuous(breaks=c(280,524),labels=c("0","1"))+
coord_cartesian(ylim=c(280,524))+
ylab("SNPs")+xlab("Generations before present")

rsid="rs4806787"
gene=fig4b_data[pop==mypop & rsID==rsid,gene_name][1]
expresp="reQTL"
start_select=fig4b_data[pop==mypop & rsID==rsid,start]
end_select=fig4b_data[pop==mypop & rsID==rsid,end]
freq_start=clues[rsID==rsid&type==expresp & epoch==start_select & pop==mypop,allele_pp_max_smooth]
freq_end=clues[rsID==rsid&type==expresp & epoch==end_select & pop==mypop,allele_pp_max_smooth]

fig4b2_plot <- ggplot(clues[rsID==rsid&type==expresp,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-50,ymax=50),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
geom_line(aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.2)+
geom_line(data=clues[rsID==rsid&type==expresp&epoch>=end_select&epoch<=start_select&pop=="CHS",],mapping=aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.5)+
geom_segment(data=data.table(x=c(end_select,start_select),y=c(50,50),yend=c(freq_end,freq_start)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
scale_y_continuous(breaks=c(0,1),labels=c("0","1"))+
xlab("Generations from present")+ylab("Allele frequency")+
coord_cartesian(ylim=c(0,1))+
theme(legend.position="none")

rsid="rs1028396"
expresp="reQTL"
start_select=fig4b_data[pop==mypop & rsID==rsid,start]
end_select=fig4b_data[pop==mypop & rsID==rsid,start]
freq_start=clues[rsID==rsid&type==expresp & epoch==start_select & pop==mypop,allele_pp_max_smooth]
freq_end=clues[rsID==rsid&type==expresp & epoch==end_select & pop==mypop,allele_pp_max_smooth]

fig4b3_plot <- ggplot(clues[rsID==rsid&type==expresp,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-50,ymax=50),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
geom_line(aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.2)+
geom_line(data=clues[rsID==rsid&type==expresp&epoch>=end_select&epoch<=start_select&pop=="CHS",],mapping=aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.5)+
geom_segment(data=data.table(x=c(end_select,start_select),y=c(-50,-50),yend=c(freq_end,freq_start)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(end_select,start_select),y=c(50,50),yend=c(freq_end,freq_start)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
scale_y_continuous(breaks=c(0,1),labels=c("0","1"))+
xlab("Generations from present")+ylab("Allele frequency")+
coord_cartesian(ylim=c(0,1))+
theme(legend.position="none")

################################################################################
# Fig. 4c

fig4c_data=selection_date[(set=="reQTLCOV")&type=="reQTL",]

PctCHS=PctCEU=PctYRI=c()
for(b in 1:1000){
PctCHS[b]=fig4c_data[pop=='CHS',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctYRI[b]=fig4c_data[pop=='YRI',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctCEU[b]=fig4c_data[pop=='CEU',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
}

cibo=setNames(
  c(quantile(PctCHS,0.025),quantile(PctYRI,0.025),quantile(PctCEU,0.025)),
  c("CHS","YRI","CEU")
)
cito=setNames(
  c(quantile(PctCHS,0.975),quantile(PctYRI,0.975),quantile(PctCEU,0.975)),
  c("CHS","YRI","CEU")
)

fig4c_data=fig4c_data[,sum(window=="start_window")/nrow(.SD)*100,by=pop]
setnames(fig4c_data,"V1","prop")
fig4c_data$cibo=cibo[fig4c_data$pop]*100
fig4c_data$cito=cito[fig4c_data$pop]*100
fig4c_data[,pop:=factor(pop,c("YRI","CEU","CHS"))]

fig4c_plot=ggplot(fig4c_data,aes(pop,prop,color=pop,fill=pop))+
  geom_col(color=NA,alpha=0.2)+
  geom_col(fill=NA,size=0.2)+
  geom_segment(aes(y=cibo,yend=cito,x=pop,xend=pop,color=pop),size=0.2)+
  scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
  ylab("Percentage COV-specific reQTL")+xlab("")


################################################################################
# Fig. 4d

genoDT=get_eQTL('rs4806787','LILRB1')
genoDT$ID.y=NULL
setnames(genoDT,c("ID.x","Number_of_ALT_alelle"),c("ID","genotype"))

Expr=getExpr('LILRB1',resolution='celltype',metric='logCPM')
Expr=Expr[celltype=="pDC",]
Expr[,state:=factor(state,c('NS','COV','IAV'))]

fig4d_data=merge(Expr,genoDT,by="IID")


fig4d_data[,logCPM_adj:=AdjustOn(rankTransform(logCPM),Number_of_ALT_alelle,.SD[,.(Age,Gender,POP)]),by=.(celltype,state)]

fig4d_plot=ggplot(fig4de_data,aes(genotype,logCPM_adj,color=state,fill=state))+
  geom_violin(alpha=0.5,color=NA,scale="width")+
  geom_violin(fill=NA,size=0.1,scale="width")+
  geom_boxplot(fill="white",alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,size=0.1,notch=T)+
  facet_grid(cols=vars(state))+
  scale_color_manual(aesthetics=c("color","fill"),values=color_conditions)+
  theme_plot()+
  theme(legend.position="none")+
  xlab("Genotype at rs4806787")+ylab(expression(italic("LILRB1")*" log"*[2]*"CPM"))

################################################################################
# Fig. 4e

genoDT=get_eQTL('rs1028396','SIRPA')
genoDT$ID.y=NULL
setnames(genoDT,c("ID.x","Number_of_ALT_alelle"),c("ID","genotype"))

Expr=getExpr('SIRPA',resolution='celltype',metric='logCPM')
Expr=Expr[celltype=="MONO.CD14",]
Expr[,state:=factor(state,c('NS','COV','IAV'))]

fig4e_data=merge(Expr,genoDT,by="IID")


fig4e_data[,logCPM_adj:=AdjustOn(rankTransform(logCPM),Number_of_ALT_alelle,.SD[,.(Age,Gender,POP)]),by=.(celltype,state)]

fig4e_plot=ggplot(fig4de_data,aes(genotype,logCPM_adj,color=state,fill=state))+
  geom_violin(alpha=0.5,color=NA,scale="width")+
  geom_violin(fill=NA,size=0.1,scale="width")+
  geom_boxplot(fill="white",alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,size=0.1,notch=T)+
  facet_grid(cols=vars(state))+
  scale_color_manual(aesthetics=c("color","fill"),values=color_conditions)+
  theme_plot()+
  theme(legend.position="none")+
  xlab("Genotype at rs1028396")+ylab(expression(italic("SIRPA")*" log"*[2]*"CPM"))

################################################################################
