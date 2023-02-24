################################################################################
################################################################################
# File name: Fig5.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Figure 5
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
source(sprintf("%s/querySNPs.R",MISC_DIR))
source(sprintf("%s/load_eQTLs.R",MISC_DIR))

# declare useful functions and variables for plotting
source(sprintf("%s/set_colors.R",MISC_DIR))
source(sprintf("%s/misc_plots.R",MISC_DIR))

# read-in library ID
args <- commandArgs(TRUE)
LIB=args[1]


########################################
########### Final Figure 5a  ###########
########################################

plot_resamp=function(resamp_results,stat=('NUM'){
  # a function that takes a subset of resampling results as input, and plots them.
  # eg. resamp_results[grepl('GWsignif_(COV|IAV|NS)',set) & stat=='PBS',]
  df=resamp_results
  df_STAT=df[,.(point=FE_NUM_SEL,
        lowci=lowerCI_NUM_SEL,
        highci=upperCI_NUM_SEL,
        pval=NUM_SEL_PVAL_ENRICHMENT),
        by=.(POP,type,celltype,state,specificity)]
  df_STAT[,FDR:=p.adjust(pval,'fdr')]
  df_STAT[,state=='NS',state:='basal']
  df_STAT[,state2:=paste(celltype,state,specificity)]
  p <- ggplot(df_STAT,aes(x = state2,y=point,color=POP)) +
    geom_errorbar(aes(ymin=lowci, ymax=highci),width=0.2,position=position_dodge(width = 0.6)) +
    geom_point(aes(ymin=lowci, ymax=highci),position=position_dodge(width = 0.6)) +
    facet_grid(. ~ type, scales = "free") +
    scale_color_manual(values=c("#eba206","#71458d","#008000")) +
    geom_hline(aes(yintercept=1), colour="black",linetype="dotted") +
    xlab("") + ylab("Fold enrichment") + labs(color="Population") +
    theme_yann() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
    p
}

# eQTL per condition
pdf(sprintf('%s/Fig5/Fig5a_condition_eQTLs_aSNP.pdf',FIG_DIR),height=6.7*.4,width=7.2*.4)
fig5a <- plot_resamp(resamp_results[grepl('eQTL_(COV|IAV|NS)$',set),])
print(fig5a)
dev.off()


########################################
########### Final Figure 5b  ###########
########################################

################################################################
###### comparison of frequency between eQTL and non eQTL #######
################################################################

comparisons=dir(sprintf('%s/tagaSNP_freqs/',DAT_RESAMP_ASNP_DIR))

stats=list()
for (i in comparisons){
  stats[[i]]=fread(sprintf('%s/tagaSNP_freqs/%s',DAT_RESAMP_ASNP_DIR,i))
}
stats=rbindlist(stats,idcol='comparison')
stats=merge(stats,unique(snpSets[,.(eQTL_set=set, num, type, celltype, state, specificity)]),by='eQTL_set')


plot_FreqDiff=function(stats){
  tags_SNPs=list()
  for (i in stats[,comparison]){
    tags_SNPs[[i]]=fread(sprintf('%s/tagaSNP_freqs/%s.gz',EVO_IMMUNO_POP_ZEUS,gsub('stats_(.*)','\\1',i)))
  }
  tags_SNPs=rbindlist(tags_SNPs,idcol='comparison')
  tags_SNPs=merge(tags_SNPs,stats,by='comparison')
  tags_SNPs[specificity!='',type:=paste(type,'breakdown')]
  tags_SNPs[state=='NS',state:='b']
  tags_SNPs[,state2:=paste(celltype,state,specificity)]
  tags_SNPs[,eQTL2:=factor(ifelse(eQTL,'eQTL','no eQTL'),c('no eQTL','eQTL'))]
  p <- ggplot(tags_SNPs,aes(x = state2,y=pmin(ASNP_FREQ,0.5),color=eQTL2,fill=POP)) +
    geom_violin(scale='width',size=.5) + geom_boxplot(fill='white',alpha=0.5,notch=TRUE,position=position_dodge(width=0.85),size=.3,outlier.size=.3) +
    facet_grid( ~ POP, scales = "free_x") +
    scale_fill_manual(values=c("#eba206","#71458d","#008000")) +
    scale_color_manual(values=c('no eQTL'="#000000",'eQTL'="#e03342"))+
    xlab("") + ylab("archaic SNP frequency") + labs(color="Population") +
    theme_yann() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + ylim(c(0,.6))
    p

}

fig5b_data=stats[grepl('eQTL_(COV|IAV|NS)',eQTL_set) & specificity=='' & type=='eQTL',]

  pdf(sprintf('%s/Fig5/Fig5b_Condition_eQTLs_wide.pdf',FIG_DIR_FREQ),height=6.7*0.4,width=7.2*0.5)
  fig5b <- plot_FreqDiff(fig5b_data)
  print(fig5b)
  dev.off()

################################################################################
# Fig. 5C inset

SNP_info=getMap(annotate=T)

snp="rs58964929"
eQTL_data=get_eQTL(snp,'UBE2F',resolution='celltype')
eQTL_data[,state:=factor(state,c('NS','COV','IAV'))]
alleles_RSID=unlist(SNP_info[ID==snp,.(REF,ALT)])
genos_RSID=paste0(alleles_RSID[c(1,1,2)],alleles_RSID[c(1,2,2)])
geno_vector=genos_RSID[1+eQTL_data$Number_of_ALT_alelle]
eQTL_data[,geno:=factor(geno_vector,genos_RSID)]
eQTL_data[,logCPM:=rankTransform(logCPM),by=.(state,celltype)]

fig5c_eQTL_data <- eQTL_data

library(DescTools)
# p <- ggplot(figs10d[state=='IAV'&celltype%chin%c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16'),],aes(x=geno,y=logCPM,color=celltype,fill=celltype))
# p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=T)
# p <- p + theme_yann()
# p <- p + scale_fill_manual(values=color_cellTypes_24level[c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16')])
# p <- p + facet_grid(state~celltype)

fig5c_eQTL_plot <- ggplot(eQTL_data[celltype=="MONO.CD14",],aes(x=geno,y=logCPM,color=state,fill=state))+
# option 1
  geom_boxplot(alpha=0.5,color=NA,outlier.size=0.1,notch=F)+
  geom_boxplot(fill=NA,outlier.size=0.1,size=0.1,notch=F)+
# option 2
  #geom_boxplot(alpha=0.5,color=NA,outlier.size=0.1,notch=T)+
  #geom_boxplot(fill=NA,outlier.size=0.1,size=0.1,notch=T)+
# option 3
  #geom_violin(alpha=0.5,color=NA,scale="width")+
  #geom_violin(fill=NA,size=0.1,scale="width")+
  #geom_boxplot(fill="white",alpha=0.5,color=NA,outlier.size=0.1,notch=T)+
  #geom_boxplot(color="black",fill=NA,outlier.size=0.1,size=0.1,notch=T)+
#geom_boxplot(alpha=1,notch=F,outlier.size=0.1) +
theme_plot(rotate.x=90)+
scale_fill_manual(aesthetics=c("color","fill"),values=condition_color) +
facet_grid(cols=vars(state))+
theme(legend.position="none",strip.background=element_blank(),strip.text=element_blank(),axis.title=element_blank())

pn="fig5c_eQTL_ube2f" # option 1
#pn="fig5c_eQTL_ube2f_notch" # option 2
#pn="fig5c_eQTL_ube2f_violin" # option 3
pname=sprintf("%s/Fig5/%s.pdf",FIG_DIR,pn)
pdf(pname,width=1.5,height=1.5)
print(fig5c_eQTL_plot)
dev.off()

#####
