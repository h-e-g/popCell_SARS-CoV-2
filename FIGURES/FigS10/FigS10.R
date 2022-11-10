################################################################################
################################################################################
# File name: FigS10.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Extended Data Figure 10
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
# Fig. S10D

SNP_info=getMap(annotate=T)

snp="rs11119346"
eQTL_data=get_eQTL(snp,'TRAF3IP3',resolution='celltype')
eQTL_data[,state:=factor(state,c('NS','COV','IAV'))]
alleles_RSID=unlist(SNP_info[ID==snp,.(REF,ALT)])
genos_RSID=paste0(alleles_RSID[c(1,1,2)],alleles_RSID[c(1,2,2)])
geno_vector=genos_RSID[1+eQTL_data$Number_of_ALT_alelle]
eQTL_data[,geno:=factor(geno_vector,genos_RSID)]
eQTL_data[,logCPM:=rankTransform(logCPM),by=.(state,celltype)]

figs10d_data <- eQTL_data

library(DescTools)
# p <- ggplot(figs10d[state=='IAV'&celltype%chin%c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16'),],aes(x=geno,y=logCPM,color=celltype,fill=celltype))
# p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=T)
# p <- p + theme_yann()
# p <- p + scale_fill_manual(values=color_cellTypes_24level[c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16')])
# p <- p + facet_grid(state~celltype)

figs10d_plot <- ggplot(eQTL_data[state=='IAV' & celltype%chin% c('MONO.CD14','MONO.CD14.INFECTED'),],aes(x=geno,y=logCPM,color=celltype,fill=celltype))+
# option 1
  #geom_boxplot(alpha=0.5,color=NA,outlier.size=0.1,notch=F)+
  #geom_boxplot(fill=NA,outlier.size=0.1,size=0.1,notch=F)+
# option 2
  #geom_boxplot(alpha=0.5,color=NA,outlier.size=0.1,notch=T)+
  #geom_boxplot(fill=NA,outlier.size=0.1,size=0.1,notch=T)+
# option 3
  geom_violin(alpha=0.5,color=NA,scale="width")+
  geom_violin(fill=NA,size=0.1,scale="width")+
  geom_boxplot(fill="white",alpha=0.5,color=NA,outlier.size=0.1,notch=T)+
  geom_boxplot(color="black",fill=NA,outlier.size=0.1,size=0.1,notch=T)+
#geom_boxplot(alpha=1,notch=F,outlier.size=0.1) +
theme_plot(rotate.x=90)+
xlab(sprintf("Genotype at %s",snp))+ylab("TRAF3IP3 (logCPM)")+
scale_fill_manual(aesthetics=c("color","fill"),values=celltype_color) +
facet_grid(cols=vars(celltype),labeller=labeller(celltype=c("MONO.CD14"="Bystanders","MONO.CD14.INFECTED"="Infected")))+theme(legend.position="none")

#pn="figs10d_traf3ip3" # option 1
#pn="figs10d_traf3ip3_notch" # option 2
pn="figs10d_traf3ip3_violin" # option 3
pname=sprintf("%s/FigS10/%s.pdf",FIG_DIR,pn)
pdf(pname,width=2,height=2)
print(figs10d_plot)
dev.off()

################################################################################
# Fig. S10F

snp="rs9520848"
eQTL_data=get_eQTL(snp,'TNFSF13B',resolution='celltype')
eQTL_data[,state:=factor(state,c('NS','COV','IAV'))]
alleles_RSID=unlist(SNP_info[ID==snp,.(REF,ALT)])
genos_RSID=paste0(alleles_RSID[c(1,1,2)],alleles_RSID[c(1,2,2)])
geno_vector=genos_RSID[1+eQTL_data$Number_of_ALT_alelle]
eQTL_data[,geno:=factor(geno_vector,genos_RSID)]
eQTL_data[,logCPM:=rankTransform(logCPM),by=.(state,celltype)]

figs10f_data <- eQTL_data

#eQTL_data=get_eQTL("rs9520848",'TNFSF13B',resolution='celltype')
#eQTL_data[,state:=factor(state,c('NS','COV','IAV'))]
#alleles_RSID=unlist(SNP_info[ID=="rs9520848",.(REF,ALT)])
#genos_RSID=paste0(alleles_RSID[c(1,1,2)],alleles_RSID[c(1,2,2)])
#geno_vector=genos_RSID[1+eQTL_data$Number_of_ALT_alelle]
#eQTL_data[,geno:=factor(geno_vector,genos_RSID)]
#eQTL_data[,logCPM:=rankTransform(logCPM),by=.(state,celltype)]

library(DescTools)
# p <- ggplot(eQTL_data[state=='IAV' & celltype%chin% c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16'),],aes(x=geno,y=logCPM,fill=celltype))
# p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=T)
# p <- p + theme_yann()
# p <- p + scale_fill_manual(values=color_cellTypes_24level[c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16')])
# p <- p + facet_grid(state~celltype)

figs10f_plot <- ggplot(eQTL_data[state=='NS'&celltype%chin%c('MAIT'),],aes(x=geno,y=logCPM,fill=celltype)) +
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
theme_plot(rotate.x=90) +
xlab(sprintf("Genotype at %s",snp))+ylab("TNFSF13B (logCPM)")+
scale_fill_manual(values=celltype_color['MAIT'])

pn="figs10f_tnfsf13b" # option 1
#pn="figs10f_tnfsf13b_notch" # option 2
#pn="figs10f_tnfsf13b_violin" # option 3
pname=sprintf("%s/FigS10/%s.pdf",FIG_DIR,pn)
pdf(pname,width=1,height=3)
print(figs10f_plot)
dev.off()

pdf(sprintf('%s/Fig6/eQTL_%s_%s_small.pdf',FIG_DIR,'TNFSF13B','rs9520848'),width=6.7*0.15,height=7.2*0.3)
print(p)
dev.off()
