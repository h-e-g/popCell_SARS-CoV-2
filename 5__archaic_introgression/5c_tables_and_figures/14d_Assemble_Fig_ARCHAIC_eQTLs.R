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
# p_crf <- ggplot(data=d_POP, aes(x=POS/1e6, y=MEAN,color=Population,alpha=Archaic)) +
#   geom_step() +
#   ylim(c(0,0.5)) +
# #  geom_hline(data=d_summ,aes(yintercept = Q99,color=Population),linetype="dashed",size=0.2) +
#   scale_color_manual("Population",values=setNames(color_populations[2:3],c('CEU','CHS'))) + #c("#66a61e","#7570b3"))
#   xlab(paste0("Chromosome ",chr," (bp)")) + ylab("Neanderthal\nancestry") +
#   xlim(c(xmin-buff,xmax+buff)/1e6) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         text = element_text(size=14),
#         legend.position="right",
#         legend.background = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white",colour = NA),
#         plot.background = element_rect(fill = "white",colour = NA),
#         strip.background=element_rect(fill="#012158"),
#         strip.text=element_text(color="white"),
#         axis.title.x=element_blank())
# 				# ,
#         # axis.text.x=element_blank(),
#         # axis.ticks.x=element_blank())
# pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_crf.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
# print(p_crf)
# dev.off()

celltype_order=c("B","B.M.K", "B.M.L", "B.N.K", "B.N.L", "Plasmablast", "MONO","MONO.CD14",
"MONO.CD14.INFECTED", "MONO.CD16", "cDC", "pDC",  "T.CD4","T.CD4.N","T.CD4.E", "T.Reg", "T.CD8", "T.CD8.N", "T.CD8.CM.EM",
"T.CD8.EMRA","MAIT", "T.gd","ILC","NK", "NK.CD56dim", "NK.M.LIKE","NK.CD56brt")

BothCTLIN_colors=setNames(c(lineage_color,celltype_color)[celltype_order],celltype_order)

d_betas=eQTL_stats[snps==RS_ID,]

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
asnps=asnps[order(factor(ORIGIN,c('SHARED','NEAND','DENI','TARGET'))),]

library(ggnewscale)

color_archaics2=c(color_archaics[c('DENI','NEAND')],'SHARED'='#FEE08B','TARGET'='red')
color_archaics2=c('NEAND'='#E59866','DENI'='#73C6B6','SHARED'='#F9E79F','TARGET'='red')
 #
	p_asnp <- ggplot() +
		geom_linerange(data=asnps,aes(x=POS/1e6, ymax=ASNP_FREQ, ymin=0,color=Population),alpha=0.5) +
		geom_point(data=asnps, aes(x=POS/1e6, y=ASNP_FREQ,color=Population,fill=ORIGIN),alpha=0.8,shape=21) +
		geom_point(data=asnps[ORIGIN=='TARGET'], aes(x=POS/1e6, y=ASNP_FREQ,color=Population,fill=ORIGIN),alpha=0.8,shape=21) +
	  scale_color_manual("Population",values=setNames(color_populations[2:3],c('CEU','CHS'))) + ylim(c(0,0.55)) +
		scale_fill_manual("Archaic",values=color_archaics2) +
		scale_fill_manual("Archaic",values=color_archaics2) + facet_grid('Introgressed'~.) +
		geom_label_repel(data=asnps,aes(x=POS/1e6, y=ASNP_FREQ,label=label,color=Population),box.padding = 0.8, max.overlaps = Inf,size=2, segment.size = 0.2, show_guide = FALSE,fill="white") +
	  # geom_hline(data=asnps_summ,aes(yintercept = q99,color=Population),linetype="dashed",size=0.2) +
	  xlab(paste0("Chromosome ",chr," (bp)")) + ylab("aSNP frequency") +
	  xlim(c(xmin-buff,xmax+buff)/1e6) +
	  theme_yann() +
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

	# pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_asnpfreq_step.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
	# print(p_asnp)
	# dev.off()
## plot eqtls
#color_conditions=c("NS"="#888888","IAV"="#6699CC","COV"="#DD5555")
eQTL_stats[,label:= ifelse(snps==RS_ID & ((celltype==TOP_CELLTYPE & state==TOP_STATE) | (celltype==BOTTOM_CELLTYPE & state==BOTTOM_STATE)), RS_ID, NA)] ## highlight best
eQTL_stats[,label:= ifelse(snps==RS_ID , RS_ID, NA)]
freqs=fread(sprintf("%s/users/Javier/results/freq/CHS_CEU_YRI_allCHR_popCellSNVs_nomono_hg38.frq.strat.gz",EIP))
# A1 = ALT alleles
# A2 = REF allele
# MAF = frequency of ALT allele

freqs_B=freqs[SNP%chin%asnps$ID,]
prov=merge(asnps,freqs_B[,.(ID=SNP,ALT=A1,REF=A2,ALT_FREQ=MAF,Population=CLST)],by=c('ID','Population')) #
asnps=merge(prov,freqs_B[CLST=='YRI',.(ID=SNP,ALT_FREQ_YRI=MAF)],by=c('ID'))
asnps[,INTROGRESSED_ALLELE:=ifelse(ALT_FREQ_YRI>.5,'REF','ALT')]

geno_region=queryRange(CHROM,start=eQTL_stats[,min(POS)],end=eQTL_stats[,max(POS)])
geno_region=melt(geno_region,id.vars=c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','CLASS','DR2','AF'))
MinCell_perCOND=500
keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
setnames(geno_region,'variable','IID')
geno_region=geno_region[as.character(IID)%chin%keptIID,]
geno_region=merge(geno_region,geno_region[ID==RS_ID,.(IID,target_SNP=value)],by='IID')
geno_region[,POP:=substr(IID,1,3)]
LD_snps=geno_region[,.(phase=sign(cor(value,target_SNP)),r2=cor(lm(value~POP)$res,lm(target_SNP~POP)$res)^2),by=.(snps=ID)]


eQTL_stats=merge(eQTL_stats,unique(asnps[,.(ID,ALT,REF,INTROGRESSED_ALLELE,ORIGIN)]),by='ID',all.x=TRUE)


eQTL_stats=merge(eQTL_stats,LD_snps,by='snps',all.x=TRUE)



p_eqtls <- ggplot(data=eQTL_stats[celltype %in% c(TOP_CELLTYPE),], aes(x=POS, y=-log10(pvalue), col=ifelse(is.na(label),state,'target'), fill=state, label=label,
                                                                                                          shape=ifelse(is.na(INTROGRESSED_ALLELE),'Non-archaic',ifelse(sign(beta)>0,'Increased','Decreased')),
                                                                                                          size=r2))

p_eqtls <- p_eqtls + rasterize(geom_point(alpha=0.5),dpi=400) + geom_point(data=eQTL_stats[celltype %in% c(TOP_CELLTYPE) & !is.na(label),],stroke = 1)
p_eqtls <- p_eqtls +  scale_color_manual("Condition",values = c(color_conditions,'target'='black')) + scale_fill_manual("Condition",values =color_conditions) + scale_size("Condition",range=c(0,2),limits=c(0,1))
  # ylim(c(0,10)) +
  # geom_label_repel(box.padding = 0.8, max.overlaps = Inf,size=2, segment.size = 0.2,show_guide = FALSE,fill="white") +

p_eqtls <- p_eqtls + scale_shape_manual("Lineage",values = c('Increased'=24,'Decreased'=25,'Non-archaic'=16))
p_eqtls <- p_eqtls + xlab(paste0("Chromosome ",chr," (bp)")) + ylab("-log10(P-value)") + xlim(c(xmin-buff,xmax+buff))
p_eqtls <- p_eqtls + theme_yann() + facet_grid(celltype ~ .) +
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







# pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_eqtl.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.4*6.7,width=7.2*.7)
# print(p_eqtls)
# dev.off()


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

#scale_x_continuous(labels = unit_format(unit = " ", scale = 1e-6))
# pdf(sprintf('%s/single_cell/project/pop_eQTL/data/5_archaics/plots/plot_%s_%s_%s_genes.pdf',EIP,GENE,RS_ID,TOP_CELLTYPE),height=0.2*6.7,width=7.2*.7)
# print(p4_chs)
# dev.off()

library(cowplot)

fig_6c=p_asnp+theme(legend.position="none",text=element_text(size=10))
fig6c_legend=get_legend(p_asnp)
fig_6d=p_eqtls+theme(legend.position="none",text=element_text(size=10))
fig6d_legend=get_legend(p_eqtls+ guides(color='none',alpha='none'))
fig_6e=p4_chs+theme(legend.position="none",text=element_text(size=10))

FIG_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/"
dir.create(paste0(FIG_DIR,'Fig6'))
pname=sprintf("%s/Fig6/Fig6_v2.pdf",FIG_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      ggplotGrob(fig_6c),
                        #grid.rect(gp=gpar(col="white")),
      fig6c_legend,
                        #grid.rect(gp=gpar(col="white")),
      ggplotGrob(fig_6d),
                        #grid.rect(gp=gpar(col="white")),
                        fig6d_legend,
                        #grid.rect(gp=gpar(col="white")),
                  ggplotGrob(fig_6e),
                        grid.rect(gp=gpar(col="white")),
    #  grid.rect(gp=gpar(col="white")),
                        grid.rect(gp=gpar(col="white"))),
    layout_matrix=rbind(
      c(1,1,1,2,2),
      c(3,3,3,4,4),
      c(5,5,5,6,6),
      c(7,7,7,7,7)), heights=c(.25,.3,.2,.25), widths=c(.2,.3,.3,.1,.1))
dev.off()

## plot eQTL for chosen tissues
rankTransform = function(x){
    percentile=rank(x,ties.method='random',na.last = NA)/(length(x)+1)
    mean_level=mean(x,na.rm=TRUE)
    sd_level=sd(x,na.rm=TRUE)
    qnorm(percentile,mean_level,sd_level)
    }


eQTL_data=rbind(get_eQTL(RS_ID,GENE,resolution='lineage'),
					get_eQTL(RS_ID,GENE,resolution='celltype'))
eQTL_data[,state:=factor(state,c('NS','COV','IAV'))]
alleles_RSID=unlist(SNP_info[ID==RS_ID,.(REF,ALT)])
genos_RSID=paste0(alleles_RSID[c(1,1,2)],alleles_RSID[c(1,2,2)])
geno_vector=genos_RSID[1+eQTL_data$Number_of_ALT_alelle]
eQTL_data[,geno:=factor(geno_vector,genos_RSID)]
eQTL_data[,logCPM:=rankTransform(logCPM),by=.(state,celltype)]



library(DescTools)
p <- ggplot(eQTL_data[celltype%chin% c('MONO'),],aes(x=geno,y=logCPM,fill=state))
p <- p + geom_boxplot(alpha=1,notch=F,outlier.size=0.1)
p <- p + theme_yann(rotate.x=90)
p <- p + scale_fill_manual(values=color_conditions)
p <- p + facet_grid(~state) +guides(fill='none') +theme(panel.spacing.x=unit(0,'mm'))

pdf(sprintf('%s/Fig6/eQTL_%s_%s_small.pdf',FIG_DIR,GENE,RS_ID),width=6.7*0.2,height=7.2*0.2)
print(p)
dev.off()

eQTL_data=get_eQTL("rs11119346",'TRAF3IP3',resolution='celltype')
eQTL_data[,state:=factor(state,c('NS','COV','IAV'))]
alleles_RSID=unlist(SNP_info[ID=="rs11119346",.(REF,ALT)])
genos_RSID=paste0(alleles_RSID[c(1,1,2)],alleles_RSID[c(1,2,2)])
geno_vector=genos_RSID[1+eQTL_data$Number_of_ALT_alelle]
eQTL_data[,geno:=factor(geno_vector,genos_RSID)]
eQTL_data[,logCPM:=rankTransform(logCPM),by=.(state,celltype)]

library(DescTools)
# p <- ggplot(eQTL_data[state=='IAV' & celltype%chin% c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16'),],aes(x=geno,y=logCPM,fill=celltype))
# p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=T)
# p <- p + theme_yann()
# p <- p + scale_fill_manual(values=color_cellTypes_24level[c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16')])
# p <- p + facet_grid(state~celltype)

p <- ggplot(eQTL_data[state=='IAV' & celltype%chin% c('MONO.CD14','MONO.CD14.INFECTED'),],aes(x=geno,y=logCPM,fill=celltype))
p <- p + geom_boxplot(alpha=1,notch=F,outlier.size=0.1)
p <- p + theme_yann(rotate.x=90)
p <- p + scale_fill_manual(values=color_conditions)
p <- p + facet_grid(~state) +guides(fill='none') +theme(panel.spacing.x=unit(0,'mm'))


pdf(sprintf('%s/Fig6/eQTL_%s_%s_small.pdf',FIG_DIR,'TRAF3IP3','rs11119346'),width=6.7*0.15,height=7.2*0.2)
print(p)
dev.off()



eQTL_data=get_eQTL("rs9520848",'TNFSF13B',resolution='celltype')
eQTL_data[,state:=factor(state,c('NS','COV','IAV'))]
alleles_RSID=unlist(SNP_info[ID=="rs9520848",.(REF,ALT)])
genos_RSID=paste0(alleles_RSID[c(1,1,2)],alleles_RSID[c(1,2,2)])
geno_vector=genos_RSID[1+eQTL_data$Number_of_ALT_alelle]
eQTL_data[,geno:=factor(geno_vector,genos_RSID)]
eQTL_data[,logCPM:=rankTransform(logCPM),by=.(state,celltype)]

library(DescTools)
# p <- ggplot(eQTL_data[state=='IAV' & celltype%chin% c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16'),],aes(x=geno,y=logCPM,fill=celltype))
# p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=T)
# p <- p + theme_yann()
# p <- p + scale_fill_manual(values=color_cellTypes_24level[c('MONO.CD14','MONO.CD14.INFECTED','MONO.CD16')])
# p <- p + facet_grid(state~celltype)

p <- ggplot(eQTL_data[state=='NS' & celltype%chin% c('MAIT'),],aes(x=geno,y=logCPM,fill=celltype))
p <- p + geom_boxplot(alpha=1,notch=F,outlier.size=0.1)
p <- p + theme_yann(rotate.x=90)
p <- p + scale_fill_manual(values=color_cellTypes_24level[c('MAIT')])
p <- p + facet_grid(~celltype) +guides(fill='none') +theme(panel.spacing.x=unit(0,'mm'))

pdf(sprintf('%s/Fig6/eQTL_%s_%s_small.pdf',FIG_DIR,'TNFSF13B','rs9520848'),width=6.7*0.15,height=7.2*0.3)
print(p)
dev.off()







###################


plot_geneList_ARCHAIC=function(gene_list,
                        addLineageState=TRUE,
                        addImmuneAnnot=TRUE,
                        addCOVID=TRUE,
                        addCOLOC=TRUE,
                        addARCHAIC=TRUE,
                        addFreqs=TRUE,
                        addAdaptive=TRUE,
                        mySize=2,legend=TRUE
                        ){
    eQTL_set=unique(gene_list[,paste(gene_name,snps,sep=' - ')])
    library(ggnewscale)

    p <- ggplot()
    p <- p + scale_y_continuous(breaks=match(eQTL_set,eQTL_set), labels=eQTL_set)
    p <- p + theme_yann()+theme(panel.border = element_blank(),axis.ticks.y = element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + xlab('') + ylab('')

    # add lineage
    x_pos_min=0
    x_pos_annot=c()
    x_pos_labels=c()

    #########@ extract list of conditions/celltypes where the eQTL is present
    # nominal eQTL
    eQTL_nominal=rbind(eQTL_Stats_celltype_mash[,.(snps,gene_name,celltype,state,beta,pvalue)],eQTL_Stats_lineage_mash[,.(snps,gene_name,celltype,state,beta,pvalue)])[pvalue<0.01,]
    eQTL_nominal[,type:='eQTL']
    reQTL_nominal=rbind(reQTL_Stats_celltype_mash[,.(snps,gene_name,celltype,state,beta,pvalue)],reQTL_Stats_lineage_mash[,.(snps,gene_name,celltype,state,beta,pvalue)])[pvalue<0.01,]
    reQTL_nominal[,type:='reQTL']
    eQTL_nominal=rbind(eQTL_nominal,reQTL_nominal)
    # GW signif eQTL
    eQTL_Signif_both[,type:='eQTL']
    reQTL_Signif_both[,type:='reQTL']
    eQTL_Signif_all=rbind(eQTL_Signif_both,reQTL_Signif_both)
    eQTL_nominal[,GWSignif:=paste(celltype,state,type,gene_name,snps)%chin%eQTL_Signif_all[,paste(celltype,state,type,gene_name,snps)]]
    # add lineage annotations
    eQTL_nominal[,lineage:=ifelse(celltype%in%lineage_5,celltype,celltype_2lineage[match(eQTL_nominal$celltype,celltype_2lineage$celltype),lineage])]
    # add DER/ANC annotation
    eQTL_nominal=merge(eQTL_nominal,unique(all_ARCHAIC[,.(snps,ALT_DERANC,INTROGRESSED_DERANC)]),by='snps',all.x=TRUE,allow.cartesian=TRUE)

    color_lineage=color_cellTypes_6level[lineage_5]
    color_celltype=c(color_lineage,color_cellTypes_24level[celltype_22])
    #color_archaics=c('AMH'=grey(0.7),'DENI'='#35978F','NEAND'='#F46D43','ARCHAIC'="#FEE08B")
    color_archaics=c('AMH'=grey(0.8),'DENI'='#80CDC1','NEAND'='#FDAE61','ARCHAIC'="#FEE08B",'UNKNOWN_ARCHAIC'='#FE8754')

    if( addLineageState ){
    #########@ assemble lineage x state panel
        DF_Lineage_state=eQTL_nominal[,.(snps,gene_name,beta,pvalue,lineage,state,type,GWSignif,celltype,ALT_DERANC,INTROGRESSED_DERANC)]
        DF_Lineage_state=DF_Lineage_state[order(snps,gene_name,lineage,-GWSignif,pvalue,type),][!duplicated(paste(snps,gene_name,lineage,state,type))]
        DF_Lineage_state=DF_Lineage_state[!duplicated(paste(snps,gene_name,lineage,state))]
        DF_Lineage_state[,x_pos:=match(paste(lineage,state,sep='-'),paste(rep(lineage_order,e=3),c('NS','COV','IAV'),sep='-'))]
        DF_Lineage_state[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_Lineage_state[,sign:=ifelse((beta*ifelse(ALT_DERANC=='DER',1,-1)*ifelse(INTROGRESSED_DERANC=='DER',1,-1))>0,'Increased','Decreased')]
        LINEAGE_STATE_GRID=data.table(expand.grid(x_pos=1:(3*length(lineage_order)),y_pos=1:length(eQTL_set)))

    #########@ add cell state panel to plot
        x_pos_min=pmax(max(x_pos_annot)+1,1)
        x_pos_annot=c(x_pos_annot,x_pos_min+1:15)
        x_pos_labels=c(x_pos_labels,paste(rep(lineage_order,e=3),c('NS','COV','IAV'),sep='-'))

        LINEAGE_STATE_GRID[,x_pos:=x_pos+x_pos_min]
        DF_Lineage_state[,x_pos:=x_pos+x_pos_min]
        p <- p + geom_tile(data=LINEAGE_STATE_GRID,aes(x_pos,y_pos),fill='white',col='black')
        p <- p + new_scale("fill")
        for (i in which(!is.na(DF_Lineage_state$y_pos))){
          p <- p + geom_tile(data=DF_Lineage_state[i,],aes(x_pos,y_pos,fill=lineage),alpha=.25,col='#000000')
         }
        p <- p + scale_fill_manual(values=color_lineage)
        p <- p + new_scale("fill") + geom_point(data=DF_Lineage_state[!is.na(y_pos),],aes(x_pos,y_pos,fill=state,alpha=log10(-log10(pvalue)),pch=sign),size=mySize)
        p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0.2,1))
        p <- p + scale_fill_manual(values=color_conditions)
				p <- p + scale_shape_manual(values=c('Increased'=24,'Decreased'=25))
				}


    if( addImmuneAnnot ){
    #########@ assemble immune panel
        DF_IMMUNE=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,][!duplicated(paste(gene_name,snps)),.(gene_name,snps,IEI,COV_VIPs,GO_immune)]
        DF_IMMUNE=melt(DF_IMMUNE,id.vars=c('gene_name','snps'))
        DF_IMMUNE[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_IMMUNE[,x_pos:=match(variable,c('GO_immune','COV_VIPs','IEI'))]
        # IMMUNE_GRID=data.table(expand.grid(x_pos=1:length(lineage_order),y_pos=1:length(eQTL_set)))

    #########@ add immune group panel to plot
        x_pos_min=pmax(max(x_pos_annot)+1,1)
        x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
        x_pos_labels=c(x_pos_labels,c('GO: immune response','COV VIPs','Inborn Errors Immunity'))

        DF_IMMUNE[,x_pos:=x_pos+x_pos_min]
        p <- p + geom_tile(data=DF_IMMUNE,aes(x_pos,y_pos),fill='white',col='black')
        for (i in which(DF_IMMUNE$value=='yes')){
          p <- p + geom_tile(data=DF_IMMUNE[i,],aes(x_pos,y_pos),fill='black')
        }
    }
    if( addARCHAIC){
      #########@ assemble Archaic panel
          DF_ARCHAIC=unique(gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,][,.(gene_name,snps,ALLELE_ORIGIN,INTROGRESSED_DERANC,IS_REINTROGRESSED=reintrogressed,Population,HAP_ORIGIN)])
          DF_ARCHAIC=dcast(DF_ARCHAIC,gene_name+snps+ALLELE_ORIGIN+INTROGRESSED_DERANC+IS_REINTROGRESSED~Population,fill='',value.var='HAP_ORIGIN')
          DF_ARCHAIC=melt(DF_ARCHAIC,id.vars=c('gene_name','snps','INTROGRESSED_DERANC'))
          DF_ARCHAIC[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
          DF_ARCHAIC[,x_pos:=match(variable,c('ALLELE_ORIGIN','CEU','CHS','IS_REINTROGRESSED'))]
          # IMMUNE_GRID=data.table(expand.grid(x_pos=1:length(lineage_order),y_pos=1:length(eQTL_set)))
          #DF_ARCHAIC[,sign:=ifelse(INTROGRESSED_DERANC=='DER','Increased','Decreased')]
      #########@ add Archaic panel to plot
          x_pos_min=pmax(max(x_pos_annot)+1,1)
          x_pos_annot=c(x_pos_annot,x_pos_min+1:4)
          x_pos_labels=c(x_pos_labels,c('Origin of eQTL SNP','Origin of introgressed haplotype (CEU)','Origin of introgressed haplotype (CHS)','Is the eQTL reintrogressed ?'))

          DF_ARCHAIC[,x_pos:=x_pos+x_pos_min]
          p <- p + new_scale("fill")
          p <- p + scale_fill_manual(values=setNames(c('white','black',color_archaics),c('','yes',names(color_archaics))))
          for (i in which(DF_ARCHAIC$value!='')){
            p <- p + geom_tile(data=DF_ARCHAIC[i,],aes(x_pos,y_pos,fill=value))
          }
					p <- p + geom_tile(data=DF_ARCHAIC,aes(x_pos,y_pos),fill='#FFFFFF00',col='black')
					  #p <- p + geom_point(data=DF_ARCHAIC[value!='',],aes(x_pos,y_pos,pch=sign),size=mySize,fill='black',alpha=.5)
            #p <- p + scale_shape_manual(values=c('Increased'=24,'Decreased'=25))
        }
    if( addFreqs ){
      #########@ assemble Frequency panel
      DF_FREQ=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,]
      DF_FREQ=unique(DF_FREQ[,.(gene_name,snps,CEU=FreqARCHAIC_CEU,CHS=FreqARCHAIC_CHS,YRI=FreqARCHAIC_YRI)])
      DF_FREQ=melt(DF_FREQ,id.vars=c('gene_name','snps'))
      DF_FREQ[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
      DF_FREQ[,x_pos:=match(variable,c('CEU','CHS','YRI'))]
      #########@ add immune group panel to plot
      x_pos_min=pmax(max(x_pos_annot)+1,1)
      x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
      x_pos_labels=c(x_pos_labels,c('Frequency of introgressed allele (CEU)','Frequency of introgressed allele (CHS)','Frequency of introgressed allele (YRI)'))

      DF_FREQ[,x_pos:=x_pos+x_pos_min]
      p <- p + geom_tile(data=DF_FREQ,aes(x_pos,y_pos),fill='white',col='black')
      p <- p + new_scale("fill")
      p <- p + scale_fill_manual(values=setNames(c('black',color_populations),c('IS_DERIVED','YRI','CEU','CHS')))
      p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0,1),limits=c(0,1))
      for (i in which(DF_FREQ$value!='')){
        p <- p + geom_tile(data=DF_FREQ[i,],aes(x_pos,y_pos,fill=variable,alpha=value))
      }
    }
    if( addAdaptive ){
      ADAPT_GRID=data.table(expand.grid(x_pos=1:2,y_pos=1:length(eQTL_set)))
      #########@ assemble Adaptive panel
      DF_ADAPT=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,]
      DF_ADAPT=unique(DF_ADAPT[,.(gene_name,snps,Population, Adaptive=(FreqARCHAIC>q95)+(FreqARCHAIC>q99))])
      DF_ADAPT=melt(DF_ADAPT,id.vars=c('gene_name','snps','Population'))
      DF_ADAPT[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
      DF_ADAPT[,x_pos:=match(Population,c('CEU','CHS'))]
      #########@ add Adaptive panel to plot
      x_pos_min=pmax(max(x_pos_annot)+1,1)
      x_pos_annot=c(x_pos_annot,x_pos_min+1:2)
      x_pos_labels=c(x_pos_labels,c('Adaptive in CEU','Adaptive in CHS'))

      DF_ADAPT[,x_pos:=x_pos+x_pos_min]
      ADAPT_GRID[,x_pos:=x_pos+x_pos_min]
      p <- p + geom_tile(data=ADAPT_GRID,aes(x_pos,y_pos),fill='white',col='grey')
      p <- p + geom_tile(data=DF_ADAPT,aes(x_pos,y_pos),fill='white',col='black')
      p <- p + new_scale("fill")
      p <- p + scale_fill_manual(values=setNames(color_populations[2:3],c('CEU','CHS')))
      p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0,1),limits=c(0,2))
      for (i in which(DF_ADAPT$value!='')){
        p <- p + geom_tile(data=DF_ADAPT[i,],aes(x_pos,y_pos,fill=Population,alpha=value))
      }
    }
    if( addCOVID ){
    #########@ assemble COVID panel
        DF_COVID=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,][!duplicated(paste(snps,gene_name,Population,type)),]
        DF_COVID=DF_COVID[,.(snps,gene_name,covid_A2_beta,covid_A2_pval,covid_B2_beta,covid_B2_pval,covid_C2_beta,covid_C2_pval,ALT_DERANC,INTROGRESSED_DERANC)]
        DF_COVID=melt(DF_COVID,id.vars=c('snps','gene_name','ALT_DERANC','INTROGRESSED_DERANC'))
        DF_COVID[,variable:=gsub('covid_','',variable)]
        DF_COVID[,stat:=gsub('([A-C]2)_(pval|beta)','\\2',variable)]
        DF_COVID[,trait:=gsub('([A-C]2)_(pval|beta)','\\1',variable)]
        DF_COVID=dcast(unique(DF_COVID),gene_name+snps+trait+ALT_DERANC+INTROGRESSED_DERANC~stat)
        DF_COVID[,x_pos:=match(trait,c('C2','B2','A2'))]
        DF_COVID[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_COVID[,sign:=ifelse((beta*ifelse(ALT_DERANC=='DER',1,-1)*ifelse(INTROGRESSED_DERANC=='DER',1,-1))>0,'Increased','Decreased')]

    #########@ add COVID panel to plot
        x_pos_min=pmax(max(x_pos_annot)+1,1)
        x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
        x_pos_labels=c(x_pos_labels,c('Reported','Hospitalized','Critical'))

        DF_COVID[,x_pos:=x_pos+x_pos_min]
        p <- p + geom_tile(data=DF_COVID,aes(x_pos,y_pos),fill='white',col='black')
        p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0.2,1))
        p <- p + geom_point(data=DF_COVID[pval<0.01,],aes(x_pos,y_pos,pch=sign,alpha=log10(-log10(pval))),fill='black',size=mySize)
      }

      if( addCOLOC ){
      #########@ assemble colocalization panel
          DF_COLOC=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,]
          DF_COLOC=DF_COLOC[,.(   PP.H4.abf_reported=max(c(PP.H4.abf_reported,0),na.rm=T),
                                  PP.H4.abf_hospitalized=max(c(PP.H4.abf_hospitalized,0),na.rm=T),
                                  PP.H4.abf_critical=max(c(PP.H4.abf_critical,0),na.rm=T)),by=.(gene_name,snps)]

          DF_COLOC=melt(DF_COLOC,id.vars=c('gene_name','snps'))
          DF_COLOC[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
          DF_COLOC[,x_pos:=match(variable,c('PP.H4.abf_reported','PP.H4.abf_hospitalized','PP.H4.abf_critical'))]
          # IMMUNE_GRID=data.table(expand.grid(x_pos=1:length(lineage_order),y_pos=1:length(eQTL_set)))

      #########@ add colocalization panel to plot
          x_pos_min=pmax(max(x_pos_annot)+1,1)
          x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
          x_pos_labels=c(x_pos_labels,c('max PPH4 (reported)','max PPH4 (hospitalized)','max PPH4 (critical)'))

          DF_COLOC[,x_pos:=x_pos+x_pos_min]
          p <- p + geom_tile(data=DF_COLOC,aes(x_pos,y_pos),fill='white',col='black')
           #for (i in which(     )){
            p <- p + new_scale("alpha") + scale_alpha_continuous(limits=c(0,1),range=c(0,1))
            p <- p + geom_tile(data=DF_COLOC,aes(x_pos,y_pos,alpha=value),fill='black')

           #}
      }
    p <- p + scale_x_continuous(breaks=x_pos_annot, labels=x_pos_labels)
    p <- p + coord_fixed()
    if(legend){
    p <- p + theme(legend.spacing.x= unit(1, 'mm'),
                    legend.spacing.y=unit(2, 'mm'),
                    legend.text=element_text(size=6),
                    legend.key = element_blank(),legend.box = "vertical",
                    legend.margin=margin(),legend.key.size = unit(2, 'mm'))
              }else{
    p <- p + theme(legend.position='none')
              }
    return(p)
}

gene_list_introgressed_coloc=all_ARCHAIC[(covid_A2_pval<0.01 | covid_C2_pval<0.01 | covid_B2_pval<0.01) & colocalized=='yes',]
p <- plot_geneList_ARCHAIC(gene_list_introgressed_coloc,
                        addLineageState=TRUE,
                        addImmuneAnnot=TRUE,
                        addCOVID=TRUE,
                        addCOLOC=TRUE,
                        addARCHAIC=TRUE,
                        addFreqs=TRUE,
                      addAdaptive=FALSE,mySize=1.5,legend=FALSE)
pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_introgressed_coloc_v8.pdf",EIP))
print(p)
dev.off()


gene_list_adaptive_introgressed_immune=all_ARCHAIC[ FreqARCHAIC>q95 & (IEI=='yes' | GO_immune=='yes' | COV_VIPs=='yes'),]

p <- plot_geneList_ARCHAIC(gene_list_adaptive_introgressed_immune,
                      addLineageState=TRUE,
                      addImmuneAnnot=TRUE,
                      addCOVID=TRUE,
                      addCOLOC=TRUE,
                      addARCHAIC=TRUE,
                      addFreqs=TRUE,
                    addAdaptive=TRUE,mySize=1.5,legend=FALSE)
pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_introgressed_adaptive_immune_v8.pdf",EIP))
print(p)
dev.off()

all_ARCHAIC=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/all_ARCHAIC_genelist.tsv",EIP))
gene_list_adaptive_introgressed=all_ARCHAIC[ FreqARCHAIC>q95,.(snps,type, gene_name,GO_immune, COV_VIPs,IEI,celltype,state, REF, ALT,which_INTROGRESSED,pvalue,beta=beta*ifelse(which_INTROGRESSED=='REF',-1,1),POP,rsID_tag_aSNP,HAP_ORIGIN,reintrogressed,Length_kb,FreqARCHAIC_CEU, FreqARCHAIC_CHS, FreqARCHAIC_YRI,colocalized,PP.H4.abf_critical, PP.H4.abf_hospitalized, PP.H4.abf_reported)]
gene_list_adaptive_introgressed=gene_list_adaptive_introgressed[order(snps,gene_name,pvalue),]
gene_list_adaptive_introgressed=gene_list_adaptive_introgressed[!duplicated(paste(POP,snps,gene_name)),]
fwrite(gene_list_adaptive_introgressed,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/SupptableS8e_Adaptively_introgressed_eQTLs.tsv",EIP),sep='\t')
# Fig6d
gene_list_adaptive_introgressed_immune=all_ARCHAIC[ FreqARCHAIC>q95 & (IEI=='yes' | GO_immune=='yes' | COV_VIPs=='yes'),]
# FigS8a
gene_list_introgressed_coloc=all_ARCHAIC[(covid_A2_pval<0.01 | covid_C2_pval<0.01 | covid_B2_pval<0.01) & colocalized=='yes',]
