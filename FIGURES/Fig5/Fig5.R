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
