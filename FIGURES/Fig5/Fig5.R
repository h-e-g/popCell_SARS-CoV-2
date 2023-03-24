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
source(sprintf("./00_one_lib_to_rule_them_all.R",LIB_DIR))

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

################################################################################
# Figure 5a


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
    theme_plot() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
    p
}
#------------------------------------------------------------------------------#
fig5a_data=fread(sprintf('%s/fig5a_data.tsv.gz',SOURCE_DATA_DIR))
# eQTL per condition
pdf(sprintf('%s/Fig5/Fig5a_condition_eQTLs_aSNP.pdf',FIG_DIR),height=6.7*.4,width=7.2*.4)
fig5a <- plot_resamp(fig5a_data[grepl('eQTL_(COV|IAV|NS)$',set),])
print(fig5a)
dev.off()
#------------------------------------------------------------------------------#

################################################################################
# Fig 5b

#------------------------------------------------------------------------------#
# extract stats for Fig5b
comparisons=dir(sprintf('%s/tagaSNP_freqs/',DAT_RESAMP_ASNP_DIR))

stats=list()
for (i in comparisons){
  stats[[i]]=fread(sprintf('%s/tagaSNP_freqs/%s',DAT_RESAMP_ASNP_DIR,i))
}
stats=rbindlist(stats,idcol='comparison')
stats=merge(stats,unique(snpSets[,.(eQTL_set=set, num, type, celltype, state, specificity)]),by='eQTL_set')
fig5b_stats=stats[grepl('eQTL_(COV|IAV|NS)',eQTL_set) & specificity=='' & type=='eQTL',]

#------------------------------------------------------------------------------#
# extract SNP level data
tags_SNPs=list()
for (i in fig5b_stats[,comparison]){
    tags_SNPs[[i]]=fread(sprintf('%s/tagaSNP_freqs/%s.gz',DAT_RESAMP_ASNP_DIR,gsub('stats_(.*)','\\1',i)))
}
tags_SNPs=rbindlist(tags_SNPs,idcol='comparison')
tags_SNPs=merge(tags_SNPs,stats,by='comparison')
tags_SNPs[specificity!='',type:=paste(type,'breakdown')]
tags_SNPs[state=='NS',state:='b']
tags_SNPs[,state2:=paste(celltype,state,specificity)]
tags_SNPs[,eQTL2:=factor(ifelse(eQTL,'eQTL','no eQTL'),c('no eQTL','eQTL'))]

#------------------------------------------------------------------------------#
# generate source data 5b
fig5b_data=tags_SNPs
fwrite(fig5b_data,file=sprintf('%s/fig5b_data.tsv.gz',SOURCE_DATA_DIR),sep='\t')

#------------------------------------------------------------------------------#
# generate source fig 5b

fig5b<- ggplot(tags_SNPs,aes(x = state2,y=pmin(ASNP_FREQ,0.5),color=eQTL2,fill=POP)) +
    geom_violin(scale='width',size=.5) + geom_boxplot(fill='white',alpha=0.5,notch=TRUE,position=position_dodge(width=0.85),size=.3,outlier.size=.3) +
    facet_grid( ~ POP, scales = "free_x") +
    scale_fill_manual(values=c("#eba206","#71458d","#008000")) +
    scale_color_manual(values=c('no eQTL'="#000000",'eQTL'="#e03342"))+
    xlab("") + ylab("archaic SNP frequency") + labs(color="Population") +
    theme_plot() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + ylim(c(0,.6))

pdf(sprintf('%s/Fig5/Fig5b_Condition_eQTLs_wide.pdf',FIG_DIR),height=6.7*0.4,width=7.2*0.5)
print(fig5b)
dev.off()

################################################################################
# Fig. 5c upper panel

color_archaics=c('NEAND'='#E59866','SHARED'='#F9E79F','TARGET'='red')
fig5c_hi_data=fread(sprintf('%s/fig5c_hi_data.tsv.gz',SOURCE_DATA_DIR))


RS_ID='rs58964929'
POS_ID='2:237965031'
chr=2
POS=237965031
GENE='UBE2F'
buff=200000
xmin=POS-buff
xmax=POS+buff

p_asnp <- ggplot() +
		geom_linerange(data=fig5c_hi_data,aes(x=POS/1e6, ymax=ASNP_FREQ, ymin=0,color=Population),alpha=0.5) +
		geom_point(data=fig5c_hi_data, aes(x=POS/1e6, y=ASNP_FREQ,color=Population,fill=ORIGIN),alpha=0.8,shape=21) +
		geom_point(data=fig5c_hi_data[ORIGIN=='TARGET'], aes(x=POS/1e6, y=ASNP_FREQ,color=Population,fill=ORIGIN),alpha=0.8,shape=21) +
	  scale_color_manual("Population",values=setNames(color_populations[2:3],c('CEU','CHS'))) + ylim(c(0,0.55)) +
		scale_fill_manual("Archaic",values=color_archaics) +
		scale_fill_manual("Archaic",values=color_archaics) + facet_grid('Introgressed'~.) +
		geom_label_repel(data=fig5c_hi_data,aes(x=POS/1e6, y=ASNP_FREQ,label=label,color=Population),box.padding = 0.8, max.overlaps = Inf,size=2, segment.size = 0.2, show_guide = FALSE,fill="white") +
	  # geom_hline(data=asnps_summ,aes(yintercept = q99,color=Population),linetype="dashed",size=0.2) +
	  xlab(paste0("Chromosome ",chr," (bp)")) + ylab("aSNP frequency") +
	  xlim(c(xmin-buff,xmax+buff)/1e6) +
	  theme_plot() +
	  theme(panel.grid = element_blank(),
	        text = element_text(size=8),
	        legend.position="right",
	        legend.background = element_rect(fill = "white"),
	        panel.background = element_rect(fill = "white",colour = NA),
	        plot.background = element_rect(fill = "white",colour = NA),
	        strip.background=element_rect(fill="#012158"),
	        strip.text=element_text(color="white"),
	        axis.title.x=element_blank(),
	        axis.text.x=element_blank(),
	        axis.ticks.x=element_blank())

################################################################################
# Fig. 5c mid panel


RS_ID='rs58964929'
POS_ID='2:237965031'
chr=2
POS=237965031
GENE='UBE2F'
buff=200000
xmin=POS-buff
xmax=POS+buff
TOP_CELLTYPE='MONO.CD14'
TOP_STATE='IAV'
RESOLUTION_TOP='celltype'


fig5c_data=eQTL_stats[celltype=='MONO.CD14',]
fig5c_data[,effect:=ifelse(is.na(INTROGRESSED_ALLELE),'Non-archaic',ifelse(sign(beta)>0,'Increased','Decreased'))]

p_eqtls <- ggplot(data=fig5c_data, aes(x=POS, y=-log10(pvalue), col=ifelse(is.na(label),state,'target'), fill=state, label=label,
                                                                                                          shape=effect,
                                                                                                          size=r2))
p_eqtls <- p_eqtls + rasterize(geom_point(alpha=0.5),dpi=400) + geom_point(data=fig5c_data[!is.na(label),],stroke = 1)
p_eqtls <- p_eqtls +  scale_color_manual("Condition",values = c(color_conditions,'target'='black')) + scale_fill_manual("Condition",values =color_conditions) + scale_size("Condition",range=c(0,2),limits=c(0,1))
p_eqtls <- p_eqtls + scale_shape_manual("Lineage",values = c('Increased'=24,'Decreased'=25,'Non-archaic'=16))
p_eqtls <- p_eqtls + xlab(paste0("Chromosome ",chr," (bp)")) + ylab("-log10(P-value)") + xlim(c(xmin-buff,xmax+buff))
p_eqtls <- p_eqtls + theme_plot() + facet_grid(celltype ~ .) +
  theme(panel.grid = element_blank(),
        legend.position="right",
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        strip.background=element_rect(fill="#012158"),
        strip.text=element_text(color="white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

################################################################################
# Fig. 5c lower panel
library(GenomicRanges)
library(ggbio)
grl <- readRDS(sprintf('%s/fig5c_data_GRangeList_UBE2F_locusgenes.RDS',SOURCE_DATA_DIR))
p_genes <- autoplot(grl, aes(type = "model"), gap.geom = "chevron",size=0.1) +
  xlab(paste0("Chromosome ",chr," (Mb)")) +
  xlim(c(xmin-buff,xmax+buff))
p_genes <- p_genes@ggplot
p_genes$layers[[8]]$aes_params$size <- 2
p_genes <- p_genes + theme_plot()

################################################################################
# Fig. 5c assemble (no inset)

fig_5c_hi=p_asnp+theme(legend.position="none",text=element_text(size=10))
fig_5c_hi_legend=get_legend(p_asnp)
fig_5c_mid=p_eqtls+theme(legend.position="none",text=element_text(size=10))
fig_5c_mid_legend=get_legend(p_eqtls+ guides(color='none',alpha='none'))
fig_5c_lo=p_genes+theme(legend.position="none",text=element_text(size=10))

pname=sprintf("%s/Fig5/Fig5c.pdf",FIG_DIR)
pdf(pname,width=7.2,height=6.7)
grid.arrange(
    grobs=list(
      ggplotGrob(fig_5c_hi),
          fig_5c_up_legend,
      ggplotGrob(fig_5c_mid),
          fig_5c_mid_legend,
      ggplotGrob(fig_5c_lo),
                        grid.rect(gp=gpar(col="white")),
    #  grid.rect(gp=gpar(col="white")),
                        grid.rect(gp=gpar(col="white"))),
    layout_matrix=rbind(
      c(1,1,1,2,2),
      c(3,3,3,4,4),
      c(5,5,5,6,6),
      c(7,7,7,7,7)), heights=c(.25,.3,.2,.25), widths=c(.2,.3,.3,.1,.1))
dev.off()


################################################################################
# Fig. 5C (mid) inset

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

fig5c_eQTL_plot <- ggplot(fig5c_eQTL_data[celltype=="MONO.CD14",],aes(x=geno,y=logCPM,color=state,fill=state))+
  geom_boxplot(alpha=0.5,color=NA,outlier.size=0.1,notch=F)+
  geom_boxplot(fill=NA,outlier.size=0.1,size=0.1,notch=F)+
  theme_plot(rotate.x=90)+
  scale_fill_manual(aesthetics=c("color","fill"),values=condition_color) +
  facet_grid(cols=vars(state))+
  theme(legend.position="none",strip.background=element_blank(),strip.text=element_blank(),axis.title=element_blank())

pname=sprintf("%s/Fig5/fig5c_eQTL_ube2f.pdf",FIG_DIR)
pdf(pname,width=1.5,height=1.5)
print(fig5c_eQTL_plot)
dev.off()

#####
