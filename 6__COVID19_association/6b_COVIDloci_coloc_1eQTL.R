#!/bin/bash
################################################################################
################################################################################
# File name: 6b_COVIDloci_coloc_launcher.R
# Author: J.MR., Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: extract sumstats and test colocalization for one eQTL
# effector script
################################################################################

# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare usefule functions
source(sprintf("%s/misc_plots.R",MISC_DIR))
source(sprintf("%s/querySNPs.R",MISC_DIR))
source(sprintf("%s/load_eQTLs.R",MISC_DIR))

###### DIRECTORY DEFINITIONS
# COLOC_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID',EVO_IMMUNO_POP_ZEUS)
# DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)
# eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
# OUT_DIR = eQTL_DIR
# SCRATCH="/pasteur/appa/scratch/mrotival/"
# RUN_NAME=""

#### read command line arguments
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
  if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test
  if (cmd[i]=='--snp' | cmd[i]=='-s' ){SNP = cmd[i+1]} # eQTL snps
  if (cmd[i]=='--gene' | cmd[i]=='-g' ){GENE = cmd[i+1]} # eQTL gene
  if (cmd[i]=='--name' | cmd[i]=='-n' ){SYMBOL = cmd[i+1]} # eQTL symbol
  if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  }


#### Intitialze parameters
if(!is.na(as.numeric(CHR))){
  CHR=paste0('chr',CHR)
}
cellstates=dir(sprintf('%s/%s',EQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV)')
feature_toUse=fread(sprintf('%s/genes_to_use.tsv',EQTL_DIR),header=F)$V1
allTRAITs=c('reported','hospitalized','critical')

WINDOW_SIZE=2*CIS_DIST
CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))

# create output DIRECTORY
dir.create(sprintf('%s/%s/',COLOC_DIR,RUN_NAME))
COLOC_OUTPUT_DIR=sprintf('%s/%s/%s_%s/',COLOC_DIR,RUN_NAME,SYMBOL,SNP)
dir.create(COLOC_OUTPUT_DIR)

# read SNP_info
SNP_info=getMap(annotate=TRUE)
myPOS=SNP_info[match(SNP,rsID),POS_B38]

# format covid summary statistics
COVID_stats=SNP_info[CHROM==CHR & POS_B38> myPOS-WINDOW_SIZE & POS_B38< myPOS+WINDOW_SIZE,.(rsID,CHROM,POS_B38, REF, ALT,covid_A2_pval, covid_B2_pval, covid_C2_pval, covid_A2_beta, covid_B2_beta, covid_C2_beta, covid_A2_sebeta, covid_B2_sebeta ,covid_C2_sebeta)]
COVID_stats=melt(COVID_stats,id.vars=c('rsID','CHROM','POS_B38','REF','ALT'))
COVID_stats[,stat:=gsub('covid_([ABC]2)_(beta|pval|sebeta)','\\2',variable)]
COVID_stats[,GWAS_type:=gsub('covid_([ABC]2)_(beta|pval|sebeta)','\\1',variable)]
COVID_stats=dcast(COVID_stats, rsID+CHROM+POS_B38+REF+ALT+GWAS_type~stat)
COVID_stats[,GWAS_type_full:=case_when(GWAS_type=='A2'~'critical',
                                        GWAS_type=='B2'~'hospitalized',
                                        GWAS_type=='C2'~'reported',
                                        TRUE~'NA')]

fwrite(COVID_stats,file=sprintf('%s/%s_%s_COVID.tsv.gz',COLOC_OUTPUT_DIR,SYMBOL,SNP),sep='\t')

# read an format eQTL summary statistics
eQTL_stats=list()
for (CELLTYPE__STATE in cellstates){
  cat(CELLTYPE__STATE,'')
  cmd=sprintf('gunzip -c %s/%s/%s/eQTL_ALL_%s_assoc.txt.gz | grep %s',EQTL_DIR,RUN_NAME,CELLTYPE__STATE,CHR,GENE)
  eQTL_stats[[CELLTYPE__STATE]]=fread(cmd)
}

eQTL_stats=rbindlist(eQTL_stats,idcol='cellstate')
setnames(eQTL_stats,c('V1','V2','V3','V4','V5','V6','V7','V8','V9'),c('snps','gene','statistic','pval','beta','celltype','R2','sebeta','CisDist'),skip_absent=TRUE)
eQTL_stats[,state:=gsub('(.*)__(.*)','\\2',cellstate)]
mm=match(eQTL_stats$snps,SNP_info$rsID)
eQTL_stats[,POS_B38:=SNP_info[mm,POS_B38]]
eQTL_stats=eQTL_stats[POS_B38>myPOS-WINDOW_SIZE & POS_B38<myPOS+WINDOW_SIZE,]

pname=sprintf("%s_%s_allCond_allCT",SYMBOL,SNP)
fwrite(eQTL_stats,file=sprintf('%s/%s_eQTL.tsv.gz',COLOC_OUTPUT_DIR,pname),sep='\t')

# extract genotypes in the region
keptIID=fread(sprintf('%s/data/IID_individual_to_use.tsv',DAT_EQTL_DIR),sep='\t',header=F)[,V1]

Genotypes=queryRange(CHR,myPOS-WINDOW_SIZE,myPOS+WINDOW_SIZE)
Genotypes=melt(Genotypes,id.vars=c("ID"),variable.name='IID')[IID%in%keptIID,]
setkey(Genotypes,ID,IID)
Genotypes[,value:=as.numeric(value)]

# perform colocalization analyses
coloc.summary=list()
coloc_eQTLstats=list()
coloc_COVIDstats=list()
for(TRAIT in allTRAITs){
  snp_set_COVID=COVID_stats[GWAS_type_full==TRAIT & !is.na(beta),rsID]
  for (myCT_S in cellstates){
    cat('\ntesting coloc',TRAIT,myCT_S,'\n')
    snp_set_eQTL=eQTL_stats[cellstate==myCT_S,snps]
    snp_set=intersect(snp_set_COVID,snp_set_eQTL)

    Geno_mat=dcast(Genotypes[ID%chin%snp_set,],IID~ID)
    Geno_res=apply(as.matrix(Geno_mat[,-'IID']),2,function(x){lm(x~substr(Geno_mat$IID,1,3))$res})
    LDMAP=cor(Geno_res)

    eQTL_stats_locus=eQTL_stats[cellstate==myCT_S & snps%chin%snp_set,][order(POS_B38),]
    COVID_stats_locus=COVID_stats[GWAS_type_full==TRAIT & %chin%snp_set,][order(POS_B38),]
    eQTL_stats_locus[is.na(sebeta),sebeta:=1e-7]

    compare_beta=merge(COVID_stats,eQTL_stats,by.x=c('rsID','POS_B38'),by.y=c('snps','POS_B38'),suffix=c('.covid','.eQTL'),allow.cartesian=TRUE)

    snp_set=eQTL_stats_locus[,snps]

    data_eQTL=list(beta=setNames(eQTL_stats_locus[,beta],eQTL_stats_locus[,snps]),
                varbeta=setNames(eQTL_stats_locus[,sebeta^2],eQTL_stats_locus[,snps]),
                snp=eQTL_stats_locus[,snps],
                position=eQTL_stats_locus[,POS_B38],
                LD=LDMAP[snp_set,snp_set],
                type='quant',
                sdY=1)

    data_COVID=list(beta=setNames(COVID_stats_locus[,beta],COVID_stats_locus[,]),
                varbeta=setNames(COVID_stats_locus[,sebeta^2],COVID_stats_locus[,]),
                snp=COVID_stats_locus[,],
                position=COVID_stats_locus[,POS_B38],
                LD=LDMAP[snp_set,snp_set],
                type='quant',
                sdY=1)

    coloc.res=coloc.signals(data_eQTL,data_COVID,p12=1e-5)
    coloc.res$summary$run=RUN_NAME

    coloc.summary_1cond=coloc.res$summary
    coloc.summary_1cond[,trait:=TRAIT]
    coloc.summary_1cond[,celltype:=gsub('(.*)__(.*)','\\1',myCT_S)]
    coloc.summary_1cond[,state:=gsub('(.*)__(.*)','\\2',myCT_S)]
    coloc.summary_1cond[,gene:=GENE]
    coloc.summary_1cond[,symbol:=SYMBOL]
    coloc.summary[[paste(TRAIT,myCT_S)]]=coloc.summary_1cond

    targetLine=summary_1cond
    RUN_NAME=targetLine[,run]
    CT=targetLine[,celltype]
    STATE=targetLine[,state]
    BEST1=targetLine[,best1] #
    BEST2=targetLine[,best2]
    BEST4=targetLine[,best4]
    HIT1=targetLine[,hit1]
    HIT2=targetLine[,hit2]

  coloc_EQTL_1_cond=eQTL_stats[celltype==CT & state==STATE,][match(c(SNP,BEST1,BEST2,BEST4),snps),]
  coloc_EQTL_1_cond[,type:=c('SNP_eQTL','BEST1_eQTL','BEST2_COVID','BEST4_COLOC')]
  coloc_EQTL_1_cond=coloc_EQTL_1_cond[,.(type,snps,celltype,state,trait=TRAIT,
                                        beta_eQTL=beta,
                                        beta_eQTL.lowerCI=sign(beta)*(abs(beta)-2*sebeta), beta_eQTL.upperCI=sign(beta)*(abs(beta)+2*sebeta),
                                        pvalue_eQTL=pval)]
  coloc_eQTLstats[[paste(TRAIT,myCT_S)]]=coloc_EQTL_1_cond

  # extarct
  coloc_COVID_1_cond=COVID_stats[GWAS_type_full==TRAIT,][match(c(SNP,BEST1,BEST2,BEST4),),]
  coloc_COVID_1_cond[,type:=c('SNP_eQTL','BEST1_eQTL','BEST2_COVID','BEST4_COLOC')]
  coloc_COVID_1_cond=coloc_COVID_1_cond[,.(type,snps,celltype=CT,state=STATE,
                                            trait=GWAS_type_full,
                                            beta_covid=beta,
                                            beta_covid.lowerCI=sign(beta)*(abs(beta)-2*sebeta), beta_covid.upperCI=sign(beta)*(abs(beta)+2*sebeta),
                                            pvalue_covid=pval)]

  coloc_COVIDstats[[paste(TRAIT,myCT_S)]]=coloc_COVID_1_cond
  }
}

coloc_COVIDstats=rbindlist(coloc_COVIDstats)
fwrite(coloc_eQTLstats,file=sprintf('%s/allTraits_allStates_coloc_signal_eQTLstats.tsv',COLOC_OUTPUT_DIR),sep='\t')
coloc_eQTLstats=rbindlist(coloc_eQTLstats)
fwrite(coloc_COVIDstats,file=sprintf('%s/allTraits_allStates_coloc_signal_Covidstats.tsv',COLOC_OUTPUT_DIR),sep='\t')
coloc.summary=rbindlist(coloc.summary)
fwrite(coloc.summary,file=sprintf('%s/allTraits_allStates_coloc_signal.tsv',COLOC_OUTPUT_DIR),sep='\t')

coloc.summary=merge(coloc.summary,dcast(coloc_COVIDstats,celltype+state+trait~type,value.var=c('beta_covid','beta_covid.lowerCI','beta_covid.upperCI','pvalue_covid')),by=c('trait','state','celltype'))
coloc.summary=merge(coloc.summary,dcast(coloc_eQTLstats,celltype+state+trait~type,value.var=c('beta_eQTL','beta_eQTL.lowerCI','beta_eQTL.upperCI','pvalue_eQTL')),by=c('trait','state','celltype'))
fwrite(coloc.summary,file=sprintf('%s/allTraits_allStates_coloc_signal_detail.tsv.gz',COLOC_OUTPUT_DIR),sep='\t')
