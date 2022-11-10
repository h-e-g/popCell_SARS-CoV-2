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
allowed_celltypes=paste(c(lineage_5,celltype_22),collapse='|')
regex=sprintf('^(r?eQTL(_GWsignif)?)_?(%s)?_*(NS|IAV|COV)?_?(shared|specific|stronger|same_strength)?',allowed_celltypes)
snpSets[,type:=gsub(regex,'\\1',set)]
snpSets[,celltype:=gsub(regex,'\\3',set)]
snpSets[,state:=gsub(regex,'\\4',set)]
snpSets[,specificity:=gsub(regex,'\\5',set)]
snpSets[,type:=ifelse(specificity!='','reQTL_breakdown',type)]

sets_of_interest=sprintf("%s_%s",rep(c("eQTL","reQTL"),each=2),rep(c("COV","COV_specific"),2))

################################################################################
# Fig. 5a

####### resample PBS scores
resamples_PBS=paste0('resampling_newq99_30_06_2022/',dir(sprintf('%s/users/Javier/results/resampling_newq99_30_06_2022',EVO_IMMUNO_POP_ZEUS)))
####### resample nSL scores
resamples_nSL=paste0('resampling/',dir(sprintf('%s/users/Javier/results/resampling',EVO_IMMUNO_POP_ZEUS),pattern='nSL'))

resamp_results=list()
for (i in c(resamples_PBS,resamples_nSL)){
  cat(i,'\n')
  resamp_results[[i]]=fread(sprintf('%s/users/Javier/results/%s',EVO_IMMUNO_POP_ZEUS,i))
  resamp_results[[i]]=resamp_results[[i]][,.(RESAMP_NUM_SEL=mean(RESAMP_NUM_SEL),
      RESAMP_MEAN_SCORE=mean(RESAMP_MEAN_SCORE),
      FE_NUM_SEL=mean((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL)),
      lowerCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.025),
      upperCI_NUM_SEL=quantile((1+OBS_NUM_SEL)/(1+RESAMP_NUM_SEL),0.975),
      FE_MEAN_SCORE=mean(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,na.rm=T),
      lowerCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.025,na.rm=T),
      upperCI_MEAN_SCORE=quantile(OBS_MEAN_SCORE/RESAMP_MEAN_SCORE,0.975,na.rm=T))
      ,by=.(OBS_NUM_SEL,NUM_SEL_PVAL,OBS_MEAN_SCORE,MEAN_SCORE_PVAL,NUM_TARGET_EQTLS,NUM_RESAMP)]
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

fig5a_data=resamp_results[grepl('eQTL_(COV|IAV|NS|shared)',set) & specificity!='stronger' & stat=="PBS",.(point=FE_MEAN_SCORE,
    lowci=lowerCI_MEAN_SCORE,
    highci=upperCI_MEAN_SCORE,
    pval=MEAN_SCORE_PVAL),
    by=.(POP,type,celltype,state,specificity)]

fig5a_data[,FDR:=p.adjust(pval,'fdr')]
fig5a_data[state=='NS',state:='b']
fig5a_data[,state2:=paste(celltype,state,specificity)]
fig5a_data[,POP:=factor(POP,c('YRI','CEU','CHS'))]
fig5a_data[,state2:=case_when(state2==" COV specific"~"COV-specific",state2==" IAV specific"~"IAV-specific",state2=="  shared"~"Shared",T~state2)]

fig5a_plot <- ggplot(fig5a_data,aes(x = state2,y=point,color=POP,fill=POP)) +
  geom_errorbar(aes(ymin=lowci, ymax=highci),width=0.2,position=position_dodge(width = 0.6),size=0.2) +
  geom_point(aes(fill=POP),position=position_dodge(width = 0.6),alpha=0.2,pch=21) +
  geom_point(aes(color=POP),position=position_dodge(width = 0.6),fill=NA,pch=21) +
  facet_grid(cols=vars(type),scales="free",labeller=labeller(type=c("eQTL"="eQTL","reQTL"="reQTL","reQTL_breakdown"="reQTL breakdown"))) +
  scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg) +
  geom_hline(aes(yintercept=1), colour="black",linetype="dotted") +
  xlab("") + ylab("Fold enrichment") + labs(color="Population") +
  theme_yann() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),panel.spacing=unit(0,'pt'))

################################################################################
# Fig. 5b
clues=fread(sprintf("%s/clues_trajectories.tsv.gz",CLUES_DIR))

selection_date=clues[abs(z_smooth)>3,.(start=max(epoch),end=min(epoch)),by=.(rsID,gene_name,type,pop)]
selection_date$PBS=clues[,setNames(PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$P_PBS=clues[,setNames(P_PBS,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date$z_max <- clues[!is.na(z_smooth),max(abs(z_smooth)),by=.(rsID,type,pop)][,setNames(V1,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]
selection_date[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
selection_date[,set:=ifelse(rsID%in%snpSets[set=="reQTL_COV",unique(snps)],'reQTLCOV','other')]
selection_date[,selected:=ifelse(P_PBS<0.05,"PBSSig",'PBSNSig')]

plot_data=selection_date[(set=="reQTLCOV")&type=="reQTL",]
plot_data$idx <- plot_data[order(pop,-start),setNames(1:nrow(plot_data),sprintf("%s%s%s",rsID,type,pop))][plot_data[,sprintf("%s%s%s",rsID,type,pop)]]

snps_to_plot=c("rs4806787","rs1028396","rs3827763","rs443099","rs2617157","rs281874806","rs7937334","rs11645448")
snps_to_plot=c("rs4806787","rs1028396")
mypop="CHS"
fig5b1_plot <- ggplot(plot_data[order(PBS),][pop==mypop,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=200,ymax=800),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
geom_segment(aes(x=start,xend=end,y=idx,yend=idx,color=selected),size=0.1)+
#geom_segment(data=data.table(x=c(721,1203),y=c(411,320),yend=c(800,800)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(721,1203),y=c(411,320),yend=c(0,0)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
#geom_segment(data=data.table(x=c(951,1634),y=c(411,320),yend=c(800,800)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
#geom_segment(data=data.table(x=c(951,1634),y=c(411,320),yend=c(800,800)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(951,1634),y=c(411,320),yend=c(0,0)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_text_repel(data=plot_data[pop==mypop&type=="reQTL"&rsID%in%snps_to_plot],aes(x=end,y=idx,label=gene_name,color=selected),seed=1,xlim=0,direction="x",segment.size=0.1,segment.alpha=0.5,show.legend=F,size=1.5,segment.linetype='dashed')+
scale_color_manual(values=c("PBSSig"="#71458d","PBSNSig"="#888888"),breaks="PBSSig",labels="Top 5% PBS",guide=guide_legend(title=""))+
scale_y_continuous(breaks=c(280,524),labels=c("0","1"))+
coord_cartesian(ylim=c(280,524))+
theme(legend.position=c(0.15,0.32),axis.ticks.y=element_line(color=NA),axis.text.y=element_text(color=NA))+
ylab("SNPs")+xlab("Generations before present")

pname=sprintf("%s/%s_%s.%s",FIG_DIR,"sel_date_",mypop,"png")
png(pname,res=480,width=8,height=6,units='in')
print(fig5b1_plot)
dev.off()

rsid="rs4806787"
expresp="reQTL"
fig5b2_plot <- ggplot(clues[rsID==rsid&type==expresp,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-50,ymax=50),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
geom_line(aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.2)+
geom_line(data=clues[rsID==rsid&type==expresp&epoch>=721&epoch<=951&pop=="CHS",],mapping=aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.5)+
#geom_segment(data=data.table(x=c(721,951),y=c(-50,-50),yend=c(0.4257766,0.1647577)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(721,951),y=c(50,50),yend=c(0.4257766,0.1647577)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
scale_y_continuous(breaks=c(0,1),labels=c("0","1"))+
xlab("Generations from present")+ylab("Allele frequency")+
coord_cartesian(ylim=c(0,1))+
theme(legend.position="none")

rsid="rs1028396"
expresp="reQTL"
fig5b3_plot <- ggplot(clues[rsID==rsid&type==expresp,])+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-50,ymax=50),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
geom_ribbon(aes(epoch,ymin=lowerCI,ymax=upperCI,fill=pop),alpha=0.25)+
geom_line(aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.2)+
geom_line(data=clues[rsID==rsid&type==expresp&epoch>=1203&epoch<=1634&pop=="CHS",],mapping=aes(epoch,allele_pp_max_smooth,color=pop),alpha=1,size=0.5)+
#geom_segment(data=data.table(x=c(1203,1634),y=c(-50,-50),yend=c(0.2540148,0.06785238)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(1203,1634),y=c(50,50),yend=c(0.2540148,0.06785238)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
scale_y_continuous(breaks=c(0,1),labels=c("0","1"))+
xlab("Generations from present")+ylab("Allele frequency")+
coord_cartesian(ylim=c(0,1))+
theme(legend.position="none")

################################################################################
# initially for lluis

selection_date_reqtlcovsp[,sum(window=="start_window")/nrow(.SD)*100,by=pop]

DT <- selection_date_reqtlcovsp
PctCHS=PctCEU=PctYRI=c()
for(b in 1:1000){
PctCHS[b]=DT[pop=='CHS',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctYRI[b]=DT[pop=='YRI',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
PctCEU[b]=DT[pop=='CEU',][sample(1:.N,.N,replace=TRUE),mean(window=="start_window" & max_abs_z>3 )]
}

cibo=setNames(
  c(quantile(PctCHS,0.025),quantile(PctYRI,0.025),quantile(PctCEU,0.025)),
  c("CHS","YRI","CEU")
)
cito=setNames(
  c(quantile(PctCHS,0.975),quantile(PctYRI,0.975),quantile(PctCEU,0.975)),
  c("CHS","YRI","CEU")
)

fig5c_data=selection_date_reqtlcovsp[,sum(window=="start_window")/nrow(.SD)*100,by=pop]
setnames(fig5c_data,"V1","prop")
fig5c_data$cibo=cibo[fig5c_data$pop]*100
fig5c_data$cito=cito[fig5c_data$pop]*100
fig5c_data[,pop:=factor(pop,c("YRI","CEU","CHS"))]

fig5c_plot=ggplot(fig5c_data,aes(pop,prop,color=pop,fill=pop))+
  geom_col(color=NA,alpha=0.2)+
  geom_col(fill=NA,size=0.2)+
  geom_segment(aes(y=cibo,yend=cito,x=pop,xend=pop,color=pop),size=0.2)+
  scale_color_manual(aesthetics=c("color","fill"),values=color_populations_1kg)+
  ylab("Percentage COV-specific reQTL")+xlab("")


pname=sprintf("%s/Fig5/Fig5cbis_lluis.pdf",FIG_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig5cbis_plot+theme(legend.position="none")),
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(2,2,2,1),
      c(2,2,2,2),
      c(2,2,2,2),
      c(2,2,2,2)
    ), heights=c(1.5,1.5,3,4), widths=c(2.25,2.25,2.5,3)
  )
dev.off()



################################################################################
# Fig. 5d
rankTransform = function(x){
    percentile=rank(x,ties.method='random',na.last = NA)/(length(x)+1)
    mean_level=mean(x,na.rm=TRUE)
    sd_level=sd(x,na.rm=TRUE)
    qnorm(percentile,mean_level,sd_level)
}

AdjustOn=function(y,x,DT){
  mod=lm(y~x+.,data=DT)
  results=residuals(mod)+mod$coef[1]+mod$coef[2]*x
  results
}

SIRPA=get_eQTL('rs1028396','SIRPA')
SIRPA$ID.y=NULL
setnames(SIPRA,c("ID.x","Number_of_ALT_alelle"),c("ID","genotype"))
SIRPA[,summary(lm(rankTransform(logCPM)~genotype+Age+Gender+POP))$coeff[2,4],by=.(celltype,state)]

SIRPA[,res_2:=AdjustOn(rankTransform(logCPM),Number_of_ALT_alelle,.SD[,.(Age,Gender,POP)]),by=.(celltype,state)]

pdf(sprintf('%s/testPlot_SIRPA.pdf',PAP_DIR))
ggplot(SIRPA[celltype=='MONO',],aes(y=res_2,x=as.factor(Number_of_ALT_alelle),fill=state))+geom_violin()+geom_boxplot(fill='white',notch=TRUE,alpha=0.5)+facet_grid(~state)+theme_yann()
print(p)
dev.off()

fig4de_data=fread(sprintf("%s/Fig4/fig54de_data.tsv.gz",FIG_DIR))
fig4de_data[,state:=factor(state,c("NS","COV","IAV"))]

fig4d_plot=ggplot(fig4de_data[celltype=="pDC"&Symbol=="LILRB1"],aes(genotype_s,logCPM_adj,color=state,fill=state))+
# option 1
  #geom_boxplot(alpha=0.5,color=NA,outlier.shape=NA,notch=F)+
  #geom_boxplot(fill=NA,outlier.shape=NA,size=0.1,notch=F)+
# option 2
  #geom_boxplot(alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  #geom_boxplot(fill=NA,outlier.shape=NA,size=0.1,notch=T)+
# option 3
  geom_violin(alpha=0.5,color=NA,scale="width")+
  geom_violin(fill=NA,size=0.1,scale="width")+
  geom_boxplot(fill="white",alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,size=0.1,notch=T)+
  facet_grid(cols=vars(state))+
  scale_color_manual(aesthetics=c("color","fill"),values=color_conditions)+
  theme_plot()+
  theme(legend.position="none")+
  xlab("Genotype at rs4806787")+ylab(expression(italic("LILRB1")*"'s logCPM (pDCs)"))

#pn="fig4d_lilrb1" # option 1
#pn="fig4d_lilrb1_notch" # option 2
pn="fig4d_lilrb1_violin" # option 3
pname=sprintf("%s/Fig4/%s.pdf",FIG_DIR,pn)
pdf(pname,width=2,height=2)
print(fig4d_plot)
dev.off()

fig4e_plot=ggplot(fig4de_data[celltype=="MONO.CD14"&Symbol=="SIRPA"],aes(genotype_s,logCPM_adj,color=state,fill=state))+
# option 1
  #geom_boxplot(alpha=0.5,color=NA,outlier.shape=NA,notch=F)+
  #geom_boxplot(fill=NA,outlier.shape=NA,size=0.1,notch=F)+
# option 2
  geom_boxplot(alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  geom_boxplot(fill=NA,outlier.shape=NA,size=0.1,notch=T)+
# option 3
  #geom_violin(alpha=0.5,color=NA,scale="width")+
  #geom_violin(fill=NA,size=0.1,scale="width")+
  #geom_boxplot(fill="white",alpha=0.5,color=NA,outlier.shape=NA,notch=T)+
  #geom_boxplot(color="black",fill=NA,outlier.shape=NA,size=0.1,notch=T)+

  facet_grid(cols=vars(state))+
  scale_color_manual(aesthetics=c("color","fill"),values=color_conditions)+
  theme_plot()+
  theme(legend.position="none")+
  xlab("Genotype at rs1028396")+ylab(expression(italic("SIRPA")*"'s logCPM (CD14+ monocytes)"))

#pn="fig4e_sirpa" # option 1
pn="fig4e_sirpa_notch" # option 2
#pn="fig4e_sirpa_violin" # option 3
pname=sprintf("%s/Fig4/%s.pdf",FIG_DIR,pn)
pdf(pname,width=2,height=2)
print(fig4e_plot)
dev.off()

################################################################################

pname=sprintf("%s/Fig5/Fig5.pdf",FIG_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig5a_plot+theme(legend.position="none")),
      ggplotGrob(fig5b1_plot+theme(legend.position="none")),
      ggplotGrob(fig5b2_plot+theme(legend.position="none",axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank())),
      ggplotGrob(fig5b3_plot+theme(legend.position="none",axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank())),
      ggplotGrob(fig5c_plot+theme(legend.position="none")),
      ggplotGrob(fig5d1_plot+theme(legend.position="none")),
      ggplotGrob(fig5d2_plot+theme(legend.position="none")),
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(1,1,1,8,2),
      c(5,6,7,8,4),
      c(5,6,7,8,4),
      c(5,6,7,8,3),
      c(8,8,8,8,3),
      c(8,8,8,8,8)
    ), heights=c(3,0.75,0.75,0.75,0.75,4), widths=c(1,2,2,1.2,2.8)
  )
dev.off()


################################################################################
