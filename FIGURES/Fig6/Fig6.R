
# declare shortcuts
MISC_DIR="MISC"
source(sprintf("%s/shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("%s/misc_functions.R",MISC_DIR))

# declare useful functions and variables for plotting
source(sprintf("%s/misc_plots.R",MISC_DIR))
source(sprintf("%s/load_eQTLs.R",MISC_DIR))
source(sprintf("%s/querySNPs.R",MISC_DIR))

# declare outputDIR
FIGURE_DIR=sprintf('%s/Fig6/',FIG_DIR)

################################################################################
# Fig 6a. enrichments

########## fucntion to plot enrichments
 plot_enrichment=function(Enrich_table,facet=FALSE,angle=0,ymax=15,...){
   color_COVID=c(reported="#F1BABA", hospitalized="#E68686",critical="#DD5555")
   Enrich_table[,GWAS_type_full:=factor(GWAS_type_full,names(color_COVID))]
   p <- ggplot(Enrich_table)+geom_point(aes(x=eQTL_group,y=FE,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5))
   p <- p + geom_errorbar(aes(x=eQTL_group,y=FE,ymin=FE_q025,ymax=FE_q975,col=GWAS_type_full),width=0.3, position=position_dodge(width = .5)) + theme_plot(rotate.x=45)
   p <- p + geom_hline(yintercept=1,col='grey') + xlab('p-value threshold')+ ylab('Fold enrichment in GWAS loci') + coord_cartesian(ylim=c(0,ymax))
   p <- p + scale_color_manual(values=color_COVID)
   if(facet){
     p <- p + facet_grid(rows=vars(eQTLtype))
     p <- p + theme(axis.text.x=element_text(angle=angle,hjust=1,vjust=1))
   }
   print(p)
   return(p)
 }

fig6a_data=fread(sprintf('%s/fig6a_data.tsv.gz',SOURCE_DATA_DIR))
textSize=8;
Fig6A <- plot_enrichment(fig6a_data,facet=FALSE,angle=45,ymax=10)
Fig6A=Fig6A+guides(col=guide_legend(ncol=1,byrow = TRUE))+theme(text = element_text(size = textSize),axis.title.x=element_blank(),legend.key.height = unit(.01, 'mm'))+ xlab('')


pdf(file=sprintf('%s/Fig6A_enrichment_COVID_loci.pdf',FIGURE_DIR),width=7.2*.2,height=6.7*.4)
print(Fig6A)
dev.off()


################################################################################
# Fig. 6b & c : colocalization plot (i.e locusZoom for eQTL & GWAS)

# Figure 6b
myCT_S='NK.CD56dim__NS'
TRAIT='hospitalized'
GENE='IRF1'
SNP="rs10066378"
MAPPING="celltype_expr"
FIG_NB='fig6b'
# Figure 6c
myCT_S='T.CD4__NS'
TRAIT='hospitalized'
GENE='IFNAR2'
SNP="rs2834158"
MAPPING="lineage_expr"
FIG_NB='fig6c'

################################################################################
# Fig. 6b & c : upper panel

fig6bc_data=fread(sprintf('%s/%s_data.tsv.gz',SOURCE_DATA_DIR,FIG_NB))
locusZoomCols=c(target="#8A39B2", highLD="#C04625", midLD="#EFC54E", lowLD="#5BB9DD", noLD="#433693")

p1 <- ggplot(fig6bc_data,aes(x=POS,y=-log10(pval.covid),col=r2_target_group,fill=r2_target_group,shape=ifelse(beta.covid>0,'pos','neg'),size=ifelse(r2_target_group=='target',TRUE,FALSE)))
p1 <- p1 + rasterize(geom_point(alpha=.7),dpi=400) + facet_grid(rows=vars(GWAS_type_full)) + scale_fill_manual(values=locusZoomCols)
p1 <- p1 + scale_color_manual(values=locusZoomCols) + scale_shape_manual(values=c("pos"=24,"neg"=25))+ ylab('COVID GWAS p-values') + scale_size_manual(values=c("TRUE"=2,"FALSE"=.5))
p1 <- p1 + theme(axis.title.x = element_text(colour = "#00000000"),axis.ticks.x = element_blank(),axis.text.x = element_text(colour = "#00000000"))
p1 <- p1 + xlim(range(fig6bc_data$POS)) + scale_x_continuous(labels = unit_format(unit = " ", scale = 1e-6))
fig6bc_hi <- p1 + theme(text = element_text(size = 5),legend.position='none')+ theme_plot()#,axis.title.x=element_blank()) #+ guides(shape=guide_legend(title='none',ncol=1),col=guide_legend(title='none',ncol=2))

pdf(sprintf('%s/%s_hi_%s_%s_covid_%s.pdf',FIG_DIR,FIG_NB,GENE,SNP,TRAIT),width=7.2*.25, height=6.7*.3)
print(fig6bc_hi)
dev.off()


################################################################################
# Fig. 6b & c : middle panel
locusZoomCols=c(target="#8A39B2", highLD="#C04625", midLD="#EFC54E", lowLD="#5BB9DD", noLD="#433693")

p2 <- ggplot(fig6bc_data,aes(x=POS,y=-log10(pval.eQTL),col=r2_target_group,fill=r2_target_group,shape=ifelse(beta.eQTL>0,'pos','neg'),size=ifelse(r2_target_group=='target',TRUE,FALSE)))+facet_grid(rows=vars(paste(gsub('\\.',' ',celltype),state,sep=' - ')))
p2 <- p2 + rasterize(geom_point(alpha=.7),dpi=400)+ scale_fill_manual(values=locusZoomCols) + ylab('eQTL p-values')
p2 <- p2 + scale_color_manual(values=locusZoomCols) + scale_shape_manual(values=c('pos'=24,'neg'=25)) + scale_size_manual(values=c("TRUE"=2,"FALSE"=.5))
p2 <- p2 + xlab(paste('chromosome',chr,'position (Mb)')) + xlim(range(fig6bc_data$POS)) + scale_x_continuous(labels = unit_format(accuracy=0.1,unit = " ", scale = 1e-6))
fig6bc_mid <- p2 + theme(text = element_text(size = 5),legend.position='none') + theme_plot()

pdf(sprintf('%s/%s_mid_%s_%s_eQTL_%s.pdf',FIG_DIR,FIG_NB,GENE,SNP,myCT_S),width=7.2*.25, height=6.7*.3)
print(fig6bc_mid)
dev.off()


################################################################################
# Fig. 6b & c : lower panel
grl=readRDS(sprintf('%s/%s_data_GRangeList_%s_locusgenes.RDS',SOURCE_DATA_DIR,FIG_NB,GENE))

# deactivate to show all gene names
names(grl)[names(grl)!=GENE]=''
# end deactivate

p_genes <- autoplot(grl, gap.geom = "chevron",size=0.1) +
    theme_plot() +
    theme(legend.position = "none",legend.text = element_text(size=6),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.1,0.2,0.1,0.1), "cm"),
          text = element_text(size = 7)) +
    xlab(paste0("Chromosome ",chr," (Mb)"))
    p_genes <- p_genes@ggplot
    p_genes$layers[[length(p4$layers)]]$aes_params$size <- 1.8
  p_genes <- p_genes + scale_x_continuous(labels = unit_format(unit = " ", scale = 1e-6)) + ylim(c(.5,4.5))
  p_genes <- p_genes + theme(axis.ticks.y=element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(), axis.line.x=element_line(), panel.border=element_blank())

  pdf(sprintf('%s/%s_locusAnnot_genes_%s_%s.pdf',FIGURE_DIR,FIG_NB,GENE,SNP),height=1, width=3)
  print(p_genes)
  dev.off()
