# faire module load samtools avant d'ouvrir R

queryRange=function(chr,start,end=NULL,add.prefix=FALSE,phased=FALSE){
  require(tictoc)
  require(stringr)
  require(data.table)
  GENO_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"
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


###################################################################################################
###########  Extract all SNPs to allow identifcation of position based on the snp.name ############
###################################################################################################
#
# EVO_IMMUNO_POP_ZEUS="/pasteur/zeus/projets/p02/evo_immuno_pop/"
# GENO_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"
# Map=list()
# require(stringr)
# require(data.table)
# require(VariantAnnotation)
# MinCell_perCOND=500
# keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]
#
# for (CHR in 22:1){
#   ImputedFILE_CHR=sprintf("Geno_b38_473Ind_3723840snps_chr%s_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct",CHR)
#   VCF=readVcf(sprintf("%s/%s.vcf.gz",GENO_DIR,ImputedFILE_CHR))
#   Freq=matrix(unlist(geno(VCF)$DS),nrow(geno(VCF)$DS),ncol(geno(VCF)$DS))
#   dimnames(Freq)=dimnames(geno(VCF)$GT)
#   rm(VCF);gc()
#   ### compute genotype frequencies
#   Freq=data.table(reshape2::melt(Freq,varnames=c('variant_id','IID')))
#   Freq[,IID:=gsub('(EvoImmunoPop|PopCell)_(ASH|EUB|AFB|PIL)([0-9]+)r?','\\2\\3',IID)]
#   Freq[,POP:=substr(IID,1,3)]
#   Freq_global=Freq[IID%chin%keptIID,.(AF_global=mean(value)/2),keyby=.(variant_id)]
#   Freq_byPop=Freq[IID%chin%keptIID,.(Freq=mean(value)/2),keyby=.(variant_id,POP)]
#   Freq_byPop=dcast(Freq_byPop,variant_id~POP)
#   Freq_byPop[,max_MAF:=pmax(pmin(AFB,1-AFB),pmin(EUB,1-EUB),pmin(ASH,1-ASH))]
#   Freq_byPop=merge(Freq_byPop,Freq_global,by='variant_id')
#   #### load further annotations
#   cmd=sprintf('bcftools view -h %s/Geno_b38_473Ind_3723840snps_chr%s_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz| tail -n 1 | cut -f 1-8', GENO_DIR, CHR)
#   NAMES=system(cmd,intern=TRUE)
#   cmd=sprintf('bcftools view -H %s/Geno_b38_473Ind_3723840snps_chr%s_shapeit4_beagle5_nodup_filtered.DR2_90pct.AF_1pct.vcf.gz | cut -f 1-8', GENO_DIR, CHR)
#   genoSTRING=system(cmd,intern=TRUE)
#   GENO=stringr:::str_split_fixed(genoSTRING,'\t',n=8)
#   colnames(GENO)=stringr:::str_split_fixed(NAMES,'\t',n=8)
#   GENO=data.table(GENO)
#   GENO[,DR2:=as.numeric(gsub('DR2=(.*);AF=([0-9\\.]*);?(.*)','\\1',INFO))]
#   GENO[,AF:=as.numeric(gsub('DR2=(.*);AF=([0-9\\.]*);?(.*)','\\2',INFO))]
#   GENO[,CLASS:=gsub('DR2=(.*);AF=([0-9\\.]*);?(.*)','\\3',INFO)]
#   GENO[,INFO:=NULL]
#   Map[[CHR]]=merge(GENO,Freq_byPop,by.x='ID',by.y='variant_id')
#   }
#  Map=rbindlist(Map)
 # fwrite(Map,file=sprintf('%s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct.tsv.gz',GENO_DIR),sep='\t')

#Map=fread(sprintf('%s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct_annotated.covid.balancing.iHS.FST.Age.tsv.gz',GENO_DIR))

getSNP=function(rsID,vector=TRUE,keep.info=FALSE,Map=NULL){
 require(data.table)
 options(datatable.fread.input.cmd.message=FALSE)
 add_names=function(vector,name){
     names(vector)=name;return(vector)
   }
 GENO_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"
 if(is.null(Map)){
   SNP_info=fread(sprintf('gunzip -c %s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct.tsv.gz| grep -e "%s\\|CHROM"',GENO_DIR, paste(rsID,collapse='\\|')))
   SNP_info=SNP_info[which(ID%chin%rsID),]
   CHROM_NAME='#CHROM'
 }else{
   SNP_info=Map[which(ID%chin%rsID),]
   CHROM_NAME='CHROM'
 }
 SNP_geno=queryRange(SNP_info[,get(CHROM_NAME)],SNP_info[,get('POS')])
 INFOS=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","DR2","AF","CLASS")
 IID=colnames(SNP_geno)[grep('AFB|EUB|ASH',colnames(SNP_geno))]
if(length(rsID)>1 & vector==TRUE){
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

getSNP_info=function(rsID){
  require(data.table)
  GENO_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"
  SNP_info=fread(sprintf('gunzip -c %s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct.tsv.gz| grep -e "%s\\|CHROM"',GENO_DIR, paste(rsID,collapse='\\|')))
  SNP_info[which(ID%chin%rsID),]
  }

getMap=function(annotate=FALSE){
  require(data.table)
  GENO_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/popCell_data/02_ImputedGenotypes/Imputed/b38"
    if(annotate){
      SNP_info=fread(sprintf('%s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct_annotated.covid.balancing.iHS.FST.Age.DistGene.ANCESTRAL.CHS.tsv.gz',GENO_DIR))
      #SNP_info=fread(sprintf('%s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct_annotated.covid.balancing.iHS.FST.Age.DistGene.tsv.gz',GENO_DIR))
    }else{
      SNP_info=fread(sprintf('%s/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct.tsv.gz',GENO_DIR))
    }
    SNP_info
    }

# genes_to_use=fread(sprintf('%s/single_cell/project/pop_eQTL/CAFEH/genes_to_use.tsv',EVO_IMMUNO_POP_ZEUS),header=F)$V1
# gene_id=genes_to_use[1]s


 getExpr=function(gene_name,resolution=c('lineage','celltype'),metric=c('logCPM','logFC'), EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"){
   DATA_DIR=sprintf('%s/single_cell/project/pop_eQTL/data/',EIP)
   resolution=match.arg(resolution)
   metric=match.arg(metric)
   gene_effect=fread(sprintf('gunzip -c %s/2_population_differences/BatchAdjusted_%s_125libs__per_%s_condition_annotated.tsv.gz | grep -e "%s\\|Symbol"',DATA_DIR,metric,resolution,gene_name))
   MinCell_perCOND=500
   keptIID=fread(sprintf('%s/1_dataset_description/keptIID_moreThan%scells.tsv',DATA_DIR,MinCell_perCOND),sep='\t',header=F)[,V1]
   gene_effect[IID%chin%keptIID,]
   }

 get_eQTL=function(rs_id,gene_name,...){
   myGene=getExpr(gene_name,...)
   mySNP=getSNP(rs_id,vector=FALSE)
   DT=merge(myGene,mySNP,by='IID')
 }
 # myGene=getExpr('MIR155HG',resolution='celltype',metric='logCPM')
 # mySNP=getSNP('rs114273142',vector=FALSE)



#
# getExpr=function(gene_id,broad=TRUE,nLibs=112){
#   celltype_category = ifelse(broad==TRUE,'cellstateBroad','cellstate')
#   celltype_category_lowercase = ifelse(broad==TRUE,'cellstate_broad','cellstate')
#   IID_effect=fread(file=sprintf('gunzip -c %s/IID_effect__pseudoBulk_full_%slibs__per_%s_clusterType_IID_LIB_annotated.tsv.gz | grep -e "%s\\|cellstate"', DATA_DIR, nLibs, celltype_category, gene_id))
#   setnames(IID_effect,'variable','IID')
#   setnames(IID_effect,celltype_category_lowercase,'cell_type')
#   IID_effect
#   }
#
#   myGene=getExpr(gene_id)
#
# sprintf('%s/IID_effect__pseudoBulk_full_%slibs__per_%s_clusterType_IID_LIB_annotated.tsv.gz',DATA_DIR,nLibs,celltype_category)

###########################################################################################
######### STEP 1: annotate data (Age, Sex, ncells, POP) and generate TPM_Matrix  ##########
###########################################################################################

# TODO make sure that paths/file names are still correct
# load raw and batch-adjusted CPM data
# Counts=fread(sprintf('%s/Counts_and_CPM__pseudoBulk_full_%slibs__per_%s_clusterType_IID_LIB.tsv.gz',DATA_DIR,nLibs,celltype_category))
