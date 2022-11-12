# SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
# EQTL_SCRIPT_DIR="09_eQTLmapping"
# EQTL_SCRIPT_DIR=""
# sbatch --array=1-59 --qos=geh -p geh --parsable --mem=30G -o ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_eQTLstats_%A_%a.log -J eQTLstats ${SCRIPT_DIR}/00_Rscript.sh 14c_Fig_ARCHAIC_eQTLs.R --gene
	# sbatch --array=1-92 --qos=geh -p geh --parsable --mem=30G -o ${SCRIPT_DIR}/09_eQTLmapping/log/logfile_eQTLstats_%A_%a.log -J eQTLstats ${SCRIPT_DIR}/00_Rscript.sh 14c_Fig_ARCHAIC_eQTLs.R --snps --gene

options(stringsAsFactors=FALSE, max.print=9999, width=200, datatable.fread.input.cmd.message=FALSE)
EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
EIP=EVO_IMMUNO_POP_ZEUS

.libPaths(sprintf("%s/single_cell/resources/R_libs/4.1.0",EVO_IMMUNO_POP_ZEUS))

suppressMessages(library(data.table))
suppressMessages(library(tictoc))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrastr))
suppressMessages(library(ggrepel))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(units))
suppressMessages(library(scales))
suppressMessages(library(ggnewscale))

eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)
RES_DIR=sprintf("%s/single_cell/resources",EVO_IMMUNO_POP_ZEUS)

CIS_DIST_TEXT='100kb'
RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_220409"
RUN_REQTL_LINEAGE="lineage_condition_logFC__logFC__CellPropLineage_SVs_220409"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_220409"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__logFC__CellPropLineage_SVs_220409"


source(sprintf("%s/template_scripts/processing_pipeline/00_set_colors.R",RES_DIR))
source(sprintf("%s/template_scripts/querySNPs.R",RES_DIR))
source(sprintf("%s/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",RES_DIR))

SNP_info=getMap(annotate=TRUE)

all_ARCHAIC=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/all_ARCHAIC_genelist.tsv",EIP))

SNP=FALSE
buff=200000
RESOLUTION=''

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--gene' | cmd[i]=='-g' ){GENE_NUM = as.numeric(cmd[i+1])} # ID of the gene to test
  if (cmd[i]=='--snps' | cmd[i]=='-s' ){SNP = TRUE} # ID of the snps to test
  if (cmd[i]=='--buff' | cmd[i]=='-b' ){buff = cmd[i+1]} # ID of the snps to test
	if (cmd[i]=='--resolution' | cmd[i]=='-r' ){RESOLUTION = cmd[i+1]} # ID of the snps to test
}

if(SNP==FALSE){
	gene_list=all_ARCHAIC[FreqARCHAIC>q95 | colocalized=='yes',unique(gene_name)]
	GENE=gene_list[GENE_NUM]
	cat(GENE)

}else{
	gene_list=unique(all_ARCHAIC[FreqARCHAIC>q95 | colocalized=='yes',.(gene_name,snps)])
	GENE=gene_list[GENE_NUM,gene_name]
	RS_ID=gene_list[GENE_NUM,snps]
	POS_ID=SNP_info[match(RS_ID,ID),chromposID_hg38]
	cat(GENE, RS_ID)
	}

# GENE="OAS1"
print(eQTL_Signif_both[gene_name==GENE,][order(pvalue),])
print(reQTL_Signif_both[gene_name==GENE,][order(pvalue),])
GENE_ID=Feature_annot[gene_name==GENE,gene_id]
CHROM=Feature_annot[gene_name==GENE,seqnames]
POS=Feature_annot[gene_name==GENE,ifelse(strand=='+',start,end)]

# load("/pasteur/appa/homes/yaaquino/checkpoints/bog_2022/0025_plot_IL10RA_CRFprop_aSNPs_eQTLs_Genes_4Yann.RData")
eQTL_stat_file=sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/eQTL_stats_%s.tsv.gz',EIP,GENE)
cat('reading', eQTL_stat_file)

if(file.exists(eQTL_stat_file)){
  eQTL_stats=fread(eQTL_stat_file)
}else{
  eQTL_stats=list()
  for (CELLTYPE in celltype_22){
    RESOLUTION='celltype'
    eQTL_stats[[CELLTYPE]]=list()
    tic(paste('load eQTLs',CELLTYPE))
    for (STATE in c('NS','IAV','COV')){
      if(CELLTYPE!='MONO.CD14.INFECTED' | STATE=='IAV'){
        filename=sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/%s_condition___CellPropLineage_SVs_220409/%s__%s/eQTL_ALL_%s_assoc.txt.gz',EIP,RESOLUTION,CELLTYPE,STATE,CHROM)
        eQTL_stats[[CELLTYPE]][[STATE]]=fread(filename)[gene==GENE_ID,]
      }
    }
    toc()
    eQTL_stats[[CELLTYPE]]=rbindlist(eQTL_stats[[CELLTYPE]],idcol='state')
  }

  for (CELLTYPE in lineage_5){
    RESOLUTION='lineage'
    eQTL_stats[[CELLTYPE]]=list()
    tic(paste('load eQTLs',CELLTYPE))
    for (STATE in c('NS','IAV','COV')){
      filename=sprintf('%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/%s_condition___CellPropLineage_SVs_220409/%s__%s/eQTL_ALL_%s_assoc.txt.gz',EIP,RESOLUTION,CELLTYPE,STATE,CHROM)
      eQTL_stats[[CELLTYPE]][[STATE]]=fread(filename)[gene==GENE_ID,]
    }
    toc()
    eQTL_stats[[CELLTYPE]]=rbindlist(eQTL_stats[[CELLTYPE]],idcol='state')
  }

  eQTL_stats=rbindlist(eQTL_stats)
  eQTL_stats[,Resolution:=ifelse(celltype%in%lineage_5,'lineage','celltype'),by=celltype]

  eQTL_stats=merge(eQTL_stats,SNP_info[,.(snps=ID,POS,ID=chromposID_hg38,posID)],by='snps')
  eQTL_stats=eQTL_stats[order(celltype,state,POS),]
  fwrite(eQTL_stats,file=sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/eQTL_stats_%s.tsv.gz',EIP,GENE))
}

eQTL_stats[,beta_low:=sign(beta)*(abs(beta)-2*se)]
eQTL_stats[,beta_high:=sign(beta)*(abs(beta)+2*se)]
eQTL_stats[,beta_low0:=sign(beta)*pmax(0,abs(beta)-2*se)]
eQTL_stats[,beta_high0:=sign(beta)*pmax(0,abs(beta)+2*se)]

if(SNP==FALSE){
  q('no')
}



if(tolower(RESOLUTION)=='lineage'){
	eQTL_Signif_both=eQTL_Signif_both[celltype%in%lineage_5,]
}
if(tolower(RESOLUTION)=='celltype'){
	eQTL_Signif_both=eQTL_Signif_both[celltype%in%celltype_22,]
}

# RS_ID='rs10774671'
TOP_CELLTYPE=eQTL_Signif_both[gene_name==GENE & snps==RS_ID,][order(pvalue)[1],celltype]
TOP_STATE=eQTL_Signif_both[gene_name==GENE & snps==RS_ID,][order(pvalue)[1],state]
RESOLUTION_TOP=ifelse(TOP_CELLTYPE%in%celltype_22,'celltype','lineage')

lowPowerCelltypes=c("ILC","cDC","pDC","MAIT",'T.gd','Plasmablast','MONO.CD14.INFECTED','NK.CD56brt')

BOTTOM_CELLTYPE=eQTL_stats[snps==RS_ID & Resolution==RESOLUTION_TOP & !(celltype %chin% lowPowerCelltypes),.(minP=min(pvalue)),by=celltype][order(-minP)[1],celltype]
BOTTOM_STATE=eQTL_stats[snps==RS_ID & celltype==BOTTOM_CELLTYPE,][order(pvalue)[1],state]


## plot eQTL for chosen tissues
eQTL_data=rbind(get_eQTL(RS_ID,GENE,resolution='lineage'),
					get_eQTL(RS_ID,GENE,resolution='celltype'))



xmin <- POS-buff
xmax <- POS+buff
chr=gsub('chr','',CHROM)

d_POP=list()
for (POP in c('CEU','CHS')){
	for (ARCHAIC in c('Vindija33.19','Denisova')){
		d_POP[[paste(POP,ARCHAIC)]]=fread(sprintf('%s/users/Javier/results/CRF/meanProp_perSNP/%s_%s_%s_Post0.9_meanPropPerSNP.txt',EIP,POP,CHROM,ARCHAIC))
		d_POP[[paste(POP,ARCHAIC)]][,Population:=POP]
		d_POP[[paste(POP,ARCHAIC)]][,Archaic:=ARCHAIC]
	}
}
d_POP=rbindlist(d_POP)

## plot crf
p_crf <- ggplot(data=d_POP, aes(x=POS/1e6, y=MEAN,color=Population,alpha=Archaic)) +
  geom_step() +
  ylim(c(0,0.5)) +
#  geom_hline(data=d_summ,aes(yintercept = Q99,color=Population),linetype="dashed",size=0.2) +
  scale_color_manual("Population",values=setNames(color_populations[2:3],c('CEU','CHS'))) + #c("#66a61e","#7570b3"))
  xlab(paste0("Chromosome ",chr," (bp)")) + ylab("Neanderthal\nancestry") +
  xlim(c(xmin-buff,xmax+buff)/1e6) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=14),
        legend.position="right",
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        strip.background=element_rect(fill="#012158"),
        strip.text=element_text(color="white"),
        axis.title.x=element_blank())
				# ,
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank())
pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_crf.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
print(p_crf)
dev.off()

d_betas=eQTL_stats[snps==RS_ID,]
# p_betas <- ggplot(data=d_betas) +
#   geom_pointrange(aes(x=celltype, y=beta, ymin=beta_low, ymax=beta_high,color=factor(state,names(color_conditions))))+
#   facet_grid(~celltype) +
#   scale_color_manual("Condition",values=color_conditions) + #c("#66a61e","#7570b3"))
#   xlab('celltype') + ylab("eQTL effect size") + theme_yann() + theme(panel.spacing=unit(0,'pt'))
#
# pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_TopSNPEffectSizes.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
# print(p_betas)
# dev.off()
celltype_order=c("B","B.M.K", "B.M.L", "B.N.K", "B.N.L", "Plasmablast", "MONO","MONO.CD14",
"MONO.CD14.INFECTED", "MONO.CD16", "cDC", "pDC",  "T.CD4","T.CD4.N","T.CD4.E", "T.Reg", "T.CD8", "T.CD8.N", "T.CD8.CM.EM",
"T.CD8.EMRA","MAIT", "T.gd","ILC","NK", "NK.CD56dim", "NK.M.LIKE","NK.CD56brt")

BothCTLIN_colors=setNames(c(lineage_color,celltype_color)[celltype_order],celltype_order)

MAX_ABS_BETA=max(abs(c(d_betas$beta_low,d_betas$beta_high)))
RG_beta=c(-MAX_ABS_BETA,MAX_ABS_BETA)
p_betas <- ggplot()
p_betas <- p_betas + geom_hline(aes(yintercept=0),col='lightgrey',linetype=2 )
p_betas <- p_betas + geom_pointrange(data=d_betas,aes(x=factor(celltype,celltype_order),y=beta, ymin=beta_low, ymax=beta_high,color=factor(state,names(color_conditions))),position=position_dodge(width=.5),size=.5,fatten =.5)
p_betas <- p_betas + scale_color_manual("Condition",values=color_conditions) #c("#66a61e","#7570b3"))
p_betas <- p_betas + xlab('celltype') + ylab("eQTL effect size") + theme_yann(rotate.x=90)
p_betas <- p_betas + geom_rect(data=data.table(xmin=1:27-0.5,xmax=1:27+0.5,ymin=-2*(1+MAX_ABS_BETA),ymax=2*(1+MAX_ABS_BETA),celltype=celltype_order),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=celltype),col=NA,alpha=.2) + scale_fill_manual("celltype",values=BothCTLIN_colors)+
 coord_cartesian(xlim=c(1,27),ylim = RG_beta)
 p_betas <- p_betas + guides(fill='none')
#p_betas <- p_betas + geom_vline(data=data.table(x=0.5*1:26),aes(xintercept=x),col='lightgrey',linetype=2)

# DONE add rectangles for each condition
pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_TopSNPEffectSizes.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
print(p_betas)
dev.off()


celltype_order2=c("B.M.K", "B.M.L", "B.N.K", "B.N.L", "Plasmablast", "MONO.CD14",
"MONO.CD14.INFECTED", "MONO.CD16", "cDC", "pDC", "T.CD4.N","T.CD4.E", "T.Reg", "T.CD8.N", "T.CD8.CM.EM",
"T.CD8.EMRA","MAIT", "T.gd","ILC","NK.CD56dim", "NK.M.LIKE","NK.CD56brt")
d_betas2=d_betas[Resolution=='celltype',]

MAX_ABS_BETA=max(abs(c(d_betas2$beta_low,d_betas2$beta_high)))
RG_beta=c(-MAX_ABS_BETA,MAX_ABS_BETA)
p_betas <- ggplot()
p_betas <- p_betas + geom_hline(aes(yintercept=0),col='lightgrey',linetype=2 )
p_betas <- p_betas + geom_pointrange(data=d_betas2,aes(x=factor(celltype,celltype_order2),y=beta, ymin=beta_low, ymax=beta_high,color=factor(state,names(color_conditions))),position=position_dodge(width=.5),size=.5,fatten =.5)
p_betas <- p_betas + scale_color_manual("Condition",values=color_conditions) #c("#66a61e","#7570b3"))
p_betas <- p_betas + xlab('celltype') + ylab("eQTL effect size") + theme_yann(rotate.x=90)
p_betas <- p_betas + geom_rect(data=data.table(xmin=1:22-0.5,xmax=1:22+0.5,ymin=-2*(1+MAX_ABS_BETA),ymax=2*(1+MAX_ABS_BETA),celltype=celltype_order2),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=celltype),col=NA,alpha=.2) + scale_fill_manual("celltype",values=celltype_color[celltype_order2]) +coord_flip()
 p_betas <- p_betas + guides(fill='none')
#p_betas <- p_betas + geom_vline(data=data.table(x=0.5*1:26),aes(xintercept=x),col='lightgrey',linetype=2)

# DONE add rectangles for each condition
pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_TopSNPEffectSizes_noLineage.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.7*6.7,width=7.2*.3)
print(p_betas)
dev.off()



## plot asnps
introgressed_snps <- fread(paste0(EIP,"/users/Javier/results/asnps/Final_lenient_Adaptive_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt"))

asnps=fread('/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/Final_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt')

asnps$label <- ifelse(asnps$ID==POS_ID,RS_ID,NA) ## highlight only best eQTL SNP
asnps=asnps[,.(CHROM,POS,ID,ORIGIN,ASNP_FREQ,METHOD,Population,label)]
asnps[ORIGIN=='ARCHAIC',ORIGIN:='SHARED']
#eQTL SNP
if(!POS_ID%in%asnps$ID){
	target_snp=unique(all_ARCHAIC[snps==RS_ID,.(CHROM,POS,ID=paste(CHROM,POS,sep=':'),ORIGIN='TARGET',METHOD='INTROGRESSED',ASNP_FREQ=FreqARCHAIC,Population,label=snps)])
	asnps=rbind(asnps,target_snp)
}else{
asnps[ID==POS_ID,ORIGIN:='TARGET']
}

asnps[setdiff(which(!is.na(label)),order(is.na(label),-ASNP_FREQ)[1]),label:=NA]

library(ggnewscale)
p_asnp <- ggplot() +
  geom_point(data=asnps, aes(x=POS/1e6, y=ASNP_FREQ,color=Population,alpha=ifelse(METHOD=='CRF','low','high'))) +
  scale_color_manual("Population",values=setNames(color_populations[2:3],c('CEU','CHS'))) +
	scale_alpha_manual("Confidence",values=setNames(c(0.2,1),c('low','high'))) +
	new_scale("alpha")+
	geom_step(data=d_POP, aes(x=POS/1e6, y=MEAN, color=Population,alpha=Archaic)) +
	geom_label_repel(data=asnps,aes(x=POS/1e6, y=ASNP_FREQ,label=label,color=Population),box.padding = 0.8, max.overlaps = Inf,size=2, segment.size = 0.2,show_guide = FALSE,fill="white") +
  # geom_hline(data=asnps_summ,aes(yintercept = q99,color=Population),linetype="dashed",size=0.2) +
  xlab(paste0("Chromosome ",chr," (bp)")) + ylab("Neanderthal \naSNP frequency") +
  xlim(c(xmin-buff,xmax+buff)/1e6) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=14),
        legend.position="right",
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        strip.background=element_rect(fill="#012158"),
        strip.text=element_text(color="white"),
        axis.title.x=element_blank())
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank())

pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_asnpfreq.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
print(p_asnp)
dev.off()

color_archaics2=c(color_archaics[c('DENI','NEAND')],'SHARED'='#FEE08B','TARGET'='red')

color_archaics2=c('NEAND'='#E59866','DENI'='#73C6B6','SHARED'='#F9E79F','TARGET'='red')
 #
	p_asnp <- ggplot() +
		geom_linerange(data=asnps,aes(x=POS/1e6, ymax=ASNP_FREQ, ymin=0,color=Population),alpha=0.5) +
		geom_point(data=asnps, aes(x=POS/1e6, y=ASNP_FREQ,color=Population,fill=ORIGIN),alpha=0.8,shape=21) +
	  scale_color_manual("Population",values=setNames(color_populations[2:3],c('CEU','CHS'))) +
		scale_fill_manual("Archaic",values=color_archaics2) +
		geom_label_repel(data=asnps,aes(x=POS/1e6, y=ASNP_FREQ,label=label,color=Population),box.padding = 0.8, max.overlaps = Inf,size=2, segment.size = 0.2, show_guide = FALSE,fill="white") +
	  # geom_hline(data=asnps_summ,aes(yintercept = q99,color=Population),linetype="dashed",size=0.2) +
	  xlab(paste0("Chromosome ",chr," (bp)")) + ylab("Neanderthal \naSNP frequency") +
	  xlim(c(xmin-buff,xmax+buff)/1e6) +
	  theme_bw() +
	  theme(panel.grid = element_blank(),
	        text = element_text(size=14),
	        legend.position="right",
	        legend.background = element_rect(fill = "white"),
	        panel.background = element_rect(fill = "white",colour = NA),
	        plot.background = element_rect(fill = "white",colour = NA),
	        strip.background=element_rect(fill="#012158"),
	        strip.text=element_text(color="white"),
	        axis.title.x=element_blank())
	        # axis.text.x=element_blank(),
	        # axis.ticks.x=element_blank())

	pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_asnpfreq_step.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
	print(p_asnp)
	dev.off()
## plot eqtls
#color_conditions=c("NS"="#888888","IAV"="#6699CC","COV"="#DD5555")
eQTL_stats[,label:= ifelse(snps==RS_ID & ((celltype==TOP_CELLTYPE & state==TOP_STATE) | (celltype==BOTTOM_CELLTYPE & state==BOTTOM_STATE)), RS_ID, NA)] ## highlight best eQTL SNP

p_eqtls <- ggplot(data=eQTL_stats[celltype %in% c(TOP_CELLTYPE,BOTTOM_CELLTYPE),], aes(x=POS, y=-log10(pvalue),shape=ifelse(sign(beta)>0,'Increased','Decreased'),col=state,fill=state,label=label))
p_eqtls <- p_eqtls + rasterize(geom_point(alpha=0.5),dpi=400) + scale_color_manual("Condition",values = color_conditions) + scale_fill_manual("Condition",values = color_conditions)
  # ylim(c(0,10)) +
  # geom_label_repel(box.padding = 0.8, max.overlaps = Inf,size=2, segment.size = 0.2,show_guide = FALSE,fill="white") +

p_eqtls <- p_eqtls + scale_shape_manual("Lineage",values = c('Increased'=24,'Decreased'=25))
p_eqtls <- p_eqtls + xlab(paste0("Chromosome ",chr," (bp)")) + ylab("-log10(P-value)") + xlim(c(xmin-buff,xmax+buff))
p_eqtls <- p_eqtls + theme_yann() + facet_grid(celltype ~ .) +
  theme(panel.grid = element_blank(),
        text = element_text(size=14),
        legend.position="right",
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        strip.background=element_rect(fill="#012158"),
        strip.text=element_text(color="white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_eqtl.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
print(p_eqtls)
dev.off()


## plot genes using ggbio
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biovizBase)
library(ggbio)
library(org.Hs.eg.db)
library(AnnotationDbi)


dtop_gr <- GRanges(CHROM, IRanges(start = xmin-buff , end = xmax+buff))
dtop_gr <- range(dtop_gr, ignore.strand = TRUE)

## hg38 transcripts
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# retrieve transcript lengths
txlen <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
setDT(txlen)
txlen$len <- rowSums(as.matrix(txlen[, .(tx_len, utr5_len, utr3_len)]))
setkey(txlen, gene_id, len, tx_id)

# filter longest transcript by gene_id
ltx <- txlen[!is.na(gene_id)][, tail(.SD,1), by=gene_id]$tx_id

# filter txdb object
txb <- as.list(txdb)
txb$transcripts <- txb$transcripts[txb$transcripts$tx_id %in% ltx, ]
txb$splicings <- txb$splicings[txb$splicings$tx_id %in% ltx,]
txb$genes <- txb$genes[txb$genes$tx_id %in% ltx,]
txb <- do.call(makeTxDb, txb)

range <- GRanges(CHROM, IRanges(start = xmin-buff , end = xmax+buff))

gr.txdb <- crunch(txb, which = range)

colnames(values(gr.txdb))[4] <- "model"
grl <- split(gr.txdb, gr.txdb$gene_id)
symbols <- AnnotationDbi::select(org.Hs.eg.db, keys=names(grl), columns="SYMBOL", keytype="ENTREZID")
names(grl) <- symbols[match(symbols$ENTREZID, names(grl), nomatch=0),"SYMBOL"]

p4_chs <- autoplot(grl, aes(type = "model"), gap.geom = "chevron",size=0.1) +
  xlab(paste0("Chromosome ",chr," (Mb)")) +
  xlim(c(xmin-buff,xmax+buff))
p4_chs <- p4_chs@ggplot
p4_chs$layers[[8]]$aes_params$size <- 2
p4_chs <- p4_chs + theme_yann()
  scale_x_continuous(labels = unit_format(unit = " ", scale = 1e-6))
	pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_genes.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.2*6.7,width=7.2*.7)
	print(p4_chs)
	dev.off()
