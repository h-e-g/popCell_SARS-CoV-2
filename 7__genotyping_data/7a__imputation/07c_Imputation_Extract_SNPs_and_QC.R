
EVO_IMMUNO_POP_ZEUS = "/pasteur/zeus/projets/p02/evo_immuno_pop"
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

GENO_DIR=sprintf("%s/popCell_data/01_GenotypingData",EVO_IMMUNO_POP_ZEUS)

IMPUTE_DIR=sprintf("%s/popCell_data/02_ImputedGenotypes",EVO_IMMUNO_POP_ZEUS)
FIGURE_DIR=sprintf("%s/single_cell/project/pop_eQTL/figures",EVO_IMMUNO_POP_ZEUS)

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(data.table))
suppressMessages(library(cowplot))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggrastr))
suppressMessages(library(scales))

suppressMessages(library(VariantAnnotation))
suppressMessages(library(vcfR))
suppressMessages(library(rtracklayer))


#Define custom ggplot2 theme
theme_mary <- function(lpos="bottom") {
  font <- "Arial"
  theme_light() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=14),
          #axis.text.x=element_text(angle=45, hjust=1),
          legend.position=lpos,
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA))
}
# add 12 entry divergent colorblind_palette # (taken from rcartocolor package: safe palette)
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
color_conditions=c('#888888','#6699CC','#AA4499')
names(color_conditions)=c('NS','COV','IAV')

color_cellTypes=c('T'="#169FD8",'NK'="#B7DC2A",'DC'="#E7b315",'Mono'="#B7637E",'B'="#2D04E5",'unassigned'='lightgrey')

# function to add names to a vector on the fly.
add_names=function(vector,name){
    names(vector)=name;return(vector)
  }

# readIn the full Genotype data
GenoFile=sprintf("%s//PopCell_473x3723840_reheaded.sorted.vcf",GENO_DIR)

VCF_full=readVcf(GenoFile, "hg38")
VCF_full=VCF_full[!duplicated(rowRanges(VCF_full)),]
GenoNum=apply(geno(VCF_full)$GT,2,function(x){
    case_when( x == '0/0' ~ 0,
               x == '0/1' ~ 1,
               x == '1/1' ~ 2,
               x == './.' ~ -9,
                TRUE~-99)})
GenoNum[GenoNum<0]=NA
dimnames(GenoNum)=dimnames(geno(VCF_full)$GT)
GenoNum_full=GenoNum
colnames(GenoNum_full)=gsub('PopCell_|EvoImmunoPop_','',colnames(GenoNum_full))

INFO=as.data.table(info(VCF_full))
CHRPOS=rowRanges(VCF_full)
# check that there are no multi-allelic SNPs
ALT_unlist=unlist(rowRanges(VCF_full)$ALT)
No_MultiAllelic=( length(ALT_unlist)==nrow(VCF_full) )
# if there is no multi-allelic SNP, replace ALT by a DNAStringSet object and unlist AF info field, to speed up & facilitate handling
INFO$AF=apply(GenoNum_full,1,mean,na.rm=T)/2
if(No_MultiAllelic){
  CHRPOS$ALT=as.character(ALT_unlist)
  CHRPOS$REF=as.character(CHRPOS$REF)
  }
CHRPOS$rsid=rownames(GenoNum_full)
INFO=cbind(as.data.table(CHRPOS),INFO)

################################################################
######## aggregate Genotype data for all other SNP_Sets ########
################################################################

list_GenoNum=list()
list_INFO=list()
for (CHR in 1:22){
  for (set_i in 1:100){
    setID=paste('chr',CHR,'_set',set_i,sep='')
    cat(setID,'\n')
    VCF_set_i=sprintf("%s/Imputed/b38/reImputationQC/_chr%s_shapeit4_beagle5_nodup_snpSet%s.vcf.gz",IMPUTE_DIR,CHR,set_i)
    VCF <- readVcf(VCF_set_i, "hg38")
    VCF=VCF[!duplicated(rowRanges(VCF)),]
    GenoNum=matrix(unlist(geno(VCF)$DS),nrow(geno(VCF)$DS),ncol(geno(VCF)$DS))
    dimnames(GenoNum)=dimnames(geno(VCF)$GT)
    INFO=as.data.table(info(VCF))
    CHRPOS=rowRanges(VCF)
  # check that there are no multi-allelic SNPs
     ALT_unlist=unlist(rowRanges(VCF)$ALT)
     No_MultiAllelic=( length(ALT_unlist)==nrow(VCF) )
  # if there is no multi-allelic SNP, replace ALT by a DNAStringSet object and unlist AF info field, to speed up & facilitate handling
     if(No_MultiAllelic){
     CHRPOS$ALT=as.character(ALT_unlist)
     INFO$AF=unlist(INFO$AF)
     CHRPOS=as.data.table(CHRPOS)
     }
    list_GenoNum[[setID]]=GenoNum
    list_INFO[[setID]]=cbind(CHRPOS,INFO)
  }
}
GenoNum_Imputed=do.call(rbind,list_GenoNum)
INFO_Imputed=rbindlist(list_INFO,idcol='setID')
colnames(GenoNum_Imputed)=gsub('PopCell_|EvoImmunoPop_','',colnames(GenoNum_Imputed))
INFO_Imputed$rsid=rownames(GenoNum_Imputed)

#################################################################################
######## compare genotyped and Imputed DATA to assess Imputation quality ########
#################################################################################

ii=intersect(rownames(GenoNum_Imputed),rownames(GenoNum_full))
INFO_Imputed[,common:=rsid%chin%ii]
INFO[,common:=rsid%chin%ii]

GenoNum_full_matched=GenoNum_full[match(ii,rownames(GenoNum_full)),match(colnames(GenoNum_Imputed),colnames(GenoNum_full))]
GenoNum_Imputed_known=GenoNum_Imputed[match(ii,rownames(GenoNum_Imputed)),]

INFO_Imputed[match(ii,rownames(GenoNum_Imputed)),Pct_Match:=apply(GenoNum_full_matched==round(GenoNum_Imputed_known),1,mean,na.rm=T)]
ComputeR2_paired=function(x){
  n=length(x)/2
    r=cor(x[1:n],x[n+1:n],use='p')
    sign(r)*r^2
}
INFO_Imputed[match(ii,rownames(GenoNum_Imputed)),R2_Imputation:=apply(cbind(GenoNum_full_matched,GenoNum_Imputed_known),1,ComputeR2_paired)]

fwrite(INFO_Imputed,file=sprintf("%s/Imputed/b38/reImputationQC/00_reImputationQC_results_run1.tsv.gz",IMPUTE_DIR),sep='\t')

pdf(sprintf("%s/ImputationQC/Pct_Match_VS_R2_byAF.pdf",FIGURE_DIR))
p <- ggplot(INFO_Imputed[sample(which(AF>0),100000),],aes(x=R2_Imputation,y=Pct_Match)) + rasterize(geom_point(size=.5,alpha=.5),dpi=300)
p <- p + theme_mary() + facet_wrap(~cut(AF,c(0,0.01,0.05,0.1,1)))
print(p)
dev.off()

pdf(sprintf("%s/ImputationQC/DR2_VS_R2_byAF.pdf",FIGURE_DIR))
p <- ggplot(INFO_Imputed[sample(which(AF>0),100000),],aes(x=DR2,y=R2_Imputation)) + rasterize(geom_point(size=.5,alpha=.5),dpi=300)
p <- p + theme_mary() + facet_wrap(~cut(AF,c(0,0.01,0.05,0.1,1)))
print(p)
dev.off()


pdf(sprintf("%s/ImputationQC/DR2_VS_R2_byAF_hexbin.pdf",FIGURE_DIR))
p <- ggplot(INFO_Imputed[sample(which(AF>0),100000),],aes(x=DR2,y=R2_Imputation)) + stat_binhex(aes(fill=log2(..count..)))
p <- p + theme_mary() + facet_wrap(~cut(AF,c(0,0.01,0.05,0.1,1)))
print(p)
dev.off()


pdf(sprintf("%s/ImputationQC/Pct_Match_VS_R2_byAF_hexbin.pdf",FIGURE_DIR))
p <- ggplot(INFO_Imputed[sample(which(AF>0),100000),],aes(x=R2_Imputation,y=Pct_Match)) + stat_binhex(aes(fill=log2(..count..)))
p <- p + theme_mary() + facet_wrap(~cut(AF,c(0,0.01,0.05,0.1,1)))
print(p)
dev.off()


INFO_Imputed[,rank_R2:=100*rank(R2_Imputation)/.N,by=cut(AF,c(0,0.01,0.05,0.1,1))]
INFO_Imputed[,rank_Pct:=100*rank(Pct_Match)/.N,by=cut(AF,c(0,0.01,0.05,0.1,1))]

pdf(sprintf("%s/ImputationQC/Pct_Match_byAF.pdf",FIGURE_DIR),height=5,width=5)
p <- ggplot(INFO_Imputed[sample(which(AF>0),100000),],aes(x=rank_Pct,y=Pct_Match,col=cut(AF,c(0,0.01,0.05,0.1,1)))) + geom_line()
p <- p + theme_mary() + scale_color_manual(values=brewer.pal(n = 9, "Greens")[c(3,5,7,9)])
p <- p + geom_hline(yintercept=1,col='lightgrey') + geom_vline(xintercept=0,col='lightgrey')
p <- p + xlab('% of all SNPs') + ylab('% match between Genotyping & imputation')
print(p)
dev.off()

pdf(sprintf("%s/ImputationQC/R2_Imputation_byAF.pdf",FIGURE_DIR),height=5,width=5)
p <- ggplot(INFO_Imputed[sample(which(AF>0),100000),],aes(x=rank_R2,y=R2_Imputation,col=cut(AF,c(0,0.01,0.05,0.1,1)))) + geom_line()
p <- p + theme_mary() + scale_color_manual(values=brewer.pal(n = 9, "Blues")[c(3,5,7,9)])
p <- p + geom_hline(yintercept=1,col='lightgrey') + geom_vline(xintercept=0,col='lightgrey')
p <- p + xlab('% of all SNPs') + ylab('(signed) R2 between Genotyping & imputation')
print(p)
dev.off()
