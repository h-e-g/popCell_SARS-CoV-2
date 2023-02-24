
################################################################################
################################################################################
# File name: querySNPs.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: functions allowing fast loading of genome wide SNP annotations
# and rapid query of genotypes at specific SNPs, or in a target region
# also contains function for rapid loading of a specific gene across multiple contexts.
#
# Effector script
#################################################################################
# Requires bcftools to work.
# make sure that bcftools is available (and on the path) before opening R
################################################################################
################################################################################
MISC_DIR='./MISC'
source(sprintf("%s/shortcuts.R",MISC_DIR))

queryRange=function(chr,start,end=NULL,add.prefix=FALSE,phased=FALSE){
  require(tictoc)
  require(stringr)
  require(data.table)
  GENO_DIR=sprintf("%s/Imputed/b38",GENOTYPE_DIR)
  # query a single range
  ncols=482 # predefined for the dataset to query
  if(add.prefix==TRUE | any(!is.na(as.numeric(chr)))){chr[!is.na(as.numeric(chr))]=paste('chr',chr[!is.na(as.numeric(chr))],sep='')}
  if(is.null(end)){
    range=paste(chr,':',format(start,scientific=FALSE,trim=TRUE),sep='')
  }else{
    range=paste(chr,':',format(start,scientific=FALSE,trim=TRUE),'-',format(end,scientific=FALSE,trim=TRUE),sep='')
  }
  DT=data.table(chr,range)
  require(stringr)
  Genotypes=list()
  for (CHR in unique(DT$chr)){
    #cat(CHR,'')
    tic('load Genotype')
    cmd=sprintf('bcftools view -h %s/Geno_b38_473Ind_3723840snps_%s_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz| tail -n 1',GENO_DIR,CHR)
    NAMES=system(cmd,intern=TRUE)
    snpList=DT[chr==CHR,range]
    ncut=1+nchar(paste(snpList,collapse=','))%/%2^16
    genoSTRING=list()
    for (cut_num in 1:ncut){
      if(ncut>1){
        snpList_subset=snpList[cut(1:length(snpList),ncut,labels=FALSE)==cut_num]
      }else{
        snpList_subset=snpList
      }
      cmd=sprintf('bcftools view -H %s/Geno_b38_473Ind_3723840snps_%s_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz %s',GENO_DIR,CHR,paste(snpList_subset,collapse=','))
      genoSTRING[[cut_num]]=system(cmd,intern=TRUE)
    }
    genoSTRING=unlist(genoSTRING)
    toc()
    tic('format Genotype - part1')
    GENO=stringr:::str_split_fixed(genoSTRING,'\t',n=ncols)
    colnames(GENO)=str_split_fixed(NAMES,'\t',n=ncols)
    # cat(dim(GENO))
    colnames(GENO)=gsub('PopCell_|EvoImmunoPop_|#','',colnames(GENO))
    toc()
    tic('format Genotype - part2')
    GENO=data.table(GENO)
    GENO[,DR2:=as.numeric(gsub('DR2=(.*);AF=([0-9\\.]*);?(.*)','\\1',INFO))]
    GENO[,AF:=as.numeric(gsub('DR2=(.*);AF=([0-9\\.]*);?(.*)','\\2',INFO))]
    GENO[,CLASS:=gsub('DR2=(.*);AF=([0-9\\.]*);?(.*)','\\3',INFO)]
    GENO[,INFO:=NULL]
    GENO[,FORMAT:=NULL]
    GENO=melt(GENO,id.vars=c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','DR2','AF','CLASS'))
    GENO[,phased:=gsub('(.*):(.*)','\\1',value)]
    GENO[,variable:=substr(variable,1,6)] # remove terminal r from replivate samples (eg AFB035r & AFB198r))
    GENO[,dosage:=as.numeric(gsub('(.*):(.*)','\\2',value))]
    toc()
    #GENO[,FullID:=paste(ID,CHROM,POS,REF,ALT,sep='_')]
    tic('format Genotype - part3')
    if(nrow(GENO)>0){
      if(phased==TRUE){
        Genotypes[[CHR]]=dcast(GENO,CHROM+POS+ID+REF+ALT+QUAL+FILTER+DR2+AF+CLASS~variable,value.var='phased')
      }else{
        Genotypes[[CHR]]=dcast(GENO,CHROM+POS+ID+REF+ALT+QUAL+FILTER+DR2+AF+CLASS~variable,value.var='dosage')
      }
    }
    toc()
  }
  rbindlist(Genotypes)
}

#### test function
# start=c(133622477,133622491,133622519,133622696,133623308,133627179,133627288,133628049,133629266,25e6)
# end=c(133622477,133622491,133622519,133622696,133623308,133627179,133627288,133628049,133629266,25.1e6)
# chr=c(rep(10,9),22)
#
# x=queryRange(chr,start,end)
# x=queryRange(chr,start,end,phased=TRUE)
# x=queryRange(paste('chr',chr,sep=''),start,end,phased=TRUE)
# x[,1:10][pmin(AF,1-AF)>0.05,.N,by=CLASS]

getSNP=function(rsID_list,vector=TRUE,keep.info=FALSE,Map=NULL){
 require(data.table)
 options(datatable.fread.input.cmd.message=FALSE)
 add_names=function(vector,name){
     names(vector)=name;return(vector)
   }
 GENO_DIR=sprintf("%s/Imputed/b38",GENOTYPE_DIR)
 if(is.null(Map)){
   SNP_info=fread(sprintf('gunzip -c %s/SNP_info_basics.tsv.gz | grep -e "%s\\|CHROM"',GENO_DIR, paste(rsID_list,collapse='\\|')))
   SNP_info=SNP_info[which(rsID%chin%rsID_list),]
   CHROM_NAME='#CHROM'
 }else{
   SNP_info=Map[which(rsID%chin%rsID_list),]
   CHROM_NAME='CHROM'
 }
 SNP_geno=queryRange(SNP_info[,get(CHROM_NAME)],SNP_info[,get('POS')])
 INFOS=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","DR2","AF","CLASS")
 IID=colnames(SNP_geno)[grep('AFB|EUB|ASH',colnames(SNP_geno))]
if(length(rsID_list)>1 & vector==TRUE){
  warning("can't output vector format with >1 SNP, setting vector=FALSE")
  vector=FALSE
}
 if(vector==FALSE){
   if(keep.info==TRUE){
     return(SNP_geno[,mget(c(INFOS,IID))])
   }else{
       return(melt(SNP_geno[,mget(c('ID',IID))],id.vars='ID',variable.name='IID',value.name='Number_of_ALT_alelle'))
   }
 }else{
   SNP_geno=melt(SNP_geno[,mget(IID)],measure.vars=IID,variable.name='IID',value.name='SNP')
   SNP_geno=add_names(SNP_geno$SNP,SNP_geno$IID)
   return(SNP_geno)
   }
}

getSNP_info=function(rsID_list){
  require(data.table)
  SNP_info=fread(sprintf('gunzip -c %s/SNP_info_basics.tsv.gz | grep -e "%s\\|CHROM"',GENO_DIR, paste(rsID_list,collapse='\\|')))
  SNP_info[which(rsID%chin%rsID_list),]
  }

getMap=function(annotate=FALSE, COVID=annotate, POPGEN=annotate, ARCHAIC=annotate){
  require(data.table)
    SNP_info=fread(sprintf('%s/SNP_info_basics.tsv.gz',SNP_INFO_DIR))

    if(COVID){
      SNP_info_covid=fread(sprintf('%s/SNP_info_covid.tsv.gz',SNP_INFO_DIR))
      SNP_info=merge(SNP_info,SNP_info_covid,by=c('rsID','posID'))
    }
    if(POPGEN){
      SNP_info_popgen=fread(sprintf('%s/SNP_info_popgen.tsv.gz',SNP_INFO_DIR))
      SNP_info=merge(SNP_info,SNP_info_popgen,by=c('rsID','posID'))
    }
    if(ARCHAIC){
      SNP_info_archaics=fread(sprintf('%s/SNP_info_archaics.tsv.gz',SNP_INFO_DIR))
      SNP_info=merge(SNP_info,SNP_info_archaics,by=c('rsID','posID'))
    }
    SNP_info
    }

 getExpr=function(gene_name,resolution=c('lineage','celltype'),metric=c('logCPM','logFC'))
   resolution=match.arg(resolution)
   metric=match.arg(metric)
   Expression=fread(sprintf('gunzip -c %s/adjusted_pseudobulk_125libs__per_%s_condition_IID.tsv.gz | grep -e "%s\\|Symbol"',EXPR_DIR, resolution,gene_name))
   keptIID=fread(sprintf('%s/data/IID_individual_to_use.tsv',DAT_EQTL_DIR),header=F)[,V1]
   Expression=Expression[IID%chin%keptIID,]
   if(metric=='logFC'){
     Expression_NS=Expression[state=="NS",]
   	if(resolution=="celltype"){
   		Expression_NS_CD14.INFECTED=Expression_NS[celltype=='MONO.CD14',]
   		Expression_NS_CD14.INFECTED[,celltype:='MONO.CD14.INFECTED']
   		Expression_NS=rbind(Expression_NS,Expression_NS_CD14.INFECTED)
   		}
   	  Expression=Expression[state!="NS",]
   	# obtain logFC
     Expression=merge(Expression,Expression_NS,by=c('IID','celltype','ID','Symbol','POP'),suffix=c('','.NS'))
   	 Expression[,logFC:=logCPM-logCPM.NS]
     }
   Expression
   }

 get_eQTL=function(rs_id,gene_name,...){
   myGene=getExpr(gene_name,...)
   mySNP=getSNP(rs_id,vector=FALSE)
   DT=merge(myGene,mySNP,by='IID')
 }
