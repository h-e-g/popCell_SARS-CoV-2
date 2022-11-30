
options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(mashr))
suppressMessages(library(tictoc))
suppressMessages(library(readr))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggrastr))
suppressMessages(library(coloc))

EVO_IMMUNO_POP_ZEUS='/pasteur/zeus/projets/p02/evo_immuno_pop'
FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/colocalization_COVID',EVO_IMMUNO_POP_ZEUS)
DATA_DIR = sprintf("%s/single_cell/project/pop_eQTL/data",EVO_IMMUNO_POP_ZEUS)
eQTL_DIR = sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping",EVO_IMMUNO_POP_ZEUS)
OUT_DIR = eQTL_DIR
SCRATCH="/pasteur/appa/scratch/mrotival/"

theme_set(theme_bw())
theme_update(
  text=element_text(family="serif",size=12),
  panel.grid=element_blank(),legend.position="bottom",
  strip.background=element_rect(fill="#012158"),strip.text=element_text(color="white")
)
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EVO_IMMUNO_POP_ZEUS))
source(sprintf("%s/single_cell/resources/template_scripts/GOSeq.R",EVO_IMMUNO_POP_ZEUS))
source(sprintf("%s/single_cell/resources/template_scripts/querySNPs.R",EVO_IMMUNO_POP_ZEUS))

RUN_NAME="/lineage_condition___CellPropLineage_SVs_220409"
CIS_DIST=1e5

#dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

TRAIT=NULL
myCELLTYPE=NULL
mySTATE=NULL

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--run_name' | cmd[i]=='-r' ){RUN_NAME = cmd[i+1]} # ID of the run
  if (cmd[i]=='--dist' | cmd[i]=='-d' ){CIS_DIST = as.numeric(cmd[i+1])} # distance to consider eQTLs
  if (cmd[i]=='--chr' | cmd[i]=='-c' ){CHR = cmd[i+1]} # CHR to test
  if (cmd[i]=='--snp' | cmd[i]=='-s' ){SNP = cmd[i+1]} # eQTL snps
  if (cmd[i]=='--gene' | cmd[i]=='-g' ){GENE = cmd[i+1]} # eQTL gene
  if (cmd[i]=='--name' | cmd[i]=='-n' ){SYMBOL = cmd[i+1]} # eQTL symbol
  if (cmd[i]=='--trait' | cmd[i]=='-t' ){TRAIT = cmd[i+1]} # COVID trait
  #if (cmd[i]=='--pos' | cmd[i]=='--position' |cmd[i]=='-p' ){myPOS = cmd[i+1]} # SNP position
  if (cmd[i]=='--celltype' |cmd[i]=='-l' ){myCELLTYPE = cmd[i+1]} # cell type
  if (cmd[i]=='--state' |cmd[i]=='-a' ){mySTATE = cmd[i+1]} # state
  }
# celltype_condition___CellPropLineage_SVs_220409 --snp rs9281525 --name AIF1 --trait reported --celltype B.N.K --state NS

if(!is.na(as.numeric(CHR))){
  CHR=paste0('chr',CHR)
}
# POS=overlap_unique_gene[i,as.numeric(gsub('[0-9]+_([0-9]+)','\\1',eQTL_peakSNP))]

#myCELLTYPE__STATE=paste(myCELLTYPE,mySTATE,sep='__')

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
WINDOW_SIZE=2*CIS_DIST

SNP_info=getMap(annotate=TRUE)

if(is.null(TRAIT)){TRAIT='allTraits'}

if(TRAIT=='allTraits' | TRAIT=='critical'){
  COVID_A2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
  COVID_A2[,SNP:=gsub(':','_',SNP)]
  mm=match(SNP_info$posID,COVID_A2$SNP)
  SNP_info[,covid_A2_pval:=COVID_A2[mm,all_inv_var_meta_p]]
  SNP_info[,covid_A2_beta:=COVID_A2[mm,all_inv_var_meta_beta]]
  SNP_info[,covid_A2_sebeta:=COVID_A2[mm,all_inv_var_meta_sebeta]]
}

if(TRAIT=='allTraits' | TRAIT=='hospitalized'){
  COVID_B2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
  COVID_B2[,SNP:=gsub(':','_',SNP)]
  mm=match(SNP_info$posID,COVID_B2$SNP)
  SNP_info[,covid_B2_pval:=COVID_B2[mm,all_inv_var_meta_p]]
  SNP_info[,covid_B2_beta:=COVID_B2[mm,all_inv_var_meta_beta]]
  SNP_info[,covid_B2_sebeta:=COVID_B2[mm,all_inv_var_meta_sebeta]]
}

if(TRAIT=='allTraits' | TRAIT=='reported'){
  COVID_C2=fread(sprintf('%s/single_cell/resources/references/hgCOVID19/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz',EVO_IMMUNO_POP_ZEUS))
  COVID_C2[,SNP:=gsub(':','_',SNP)]
  mm=match(SNP_info$posID,COVID_C2$SNP)
  SNP_info[,covid_C2_pval:=COVID_C2[mm,all_inv_var_meta_p]]
  SNP_info[,covid_C2_beta:=COVID_C2[mm,all_inv_var_meta_beta]]
  SNP_info[,covid_C2_sebeta:=COVID_C2[mm,all_inv_var_meta_sebeta]]
}

myPOS=SNP_info[match(SNP,ID),POS]

COVID_stats=SNP_info[CHROM==CHR & POS> myPOS-WINDOW_SIZE & POS< myPOS+WINDOW_SIZE,.(ID,CHROM,POS, REF, ALT,covid_A2_pval, covid_B2_pval, covid_C2_pval, covid_A2_beta,  covid_B2_beta,covid_C2_beta,covid_A2_sebeta,  covid_B2_sebeta,covid_C2_sebeta)]
COVID_stats=melt(COVID_stats,id.vars=c('ID','CHROM','POS','REF','ALT'))
COVID_stats[,stat:=gsub('covid_([ABC]2)_(beta|pval|sebeta)','\\2',variable)]
COVID_stats[,GWAS_type:=gsub('covid_([ABC]2)_(beta|pval|sebeta)','\\1',variable)]
COVID_stats=dcast(COVID_stats, ID+CHROM+POS+REF+ALT+GWAS_type~stat)
COVID_stats[,GWAS_type_full:=case_when(GWAS_type=='A2'~'critical',
                                        GWAS_type=='B2'~'hospitalized',
                                        GWAS_type=='C2'~'reported',
                                        TRUE~'NA')]

dir.create(sprintf('%s/%s/',FIGURE_DIR,RUN_NAME))
dir.create(sprintf('%s/%s/Dist_%s',FIGURE_DIR,RUN_NAME,CIS_DIST_TEXT))

COLOC_DIR=sprintf('%s/%s/Dist_%s/%s_%s/',FIGURE_DIR,RUN_NAME,CIS_DIST_TEXT,SYMBOL,SNP)
dir.create(COLOC_DIR)


fwrite(COVID_stats,file=sprintf('%s/%s_%s_COVID.tsv.gz',COLOC_DIR,SYMBOL,SNP),sep='\t')

if(TRAIT=='allTraits'){
  p1 <- ggplot(COVID_stats,aes(x=POS,y=-log10(pval),col=GWAS_type_full,fill=GWAS_type_full,shape=ifelse(beta>0,'pos','neg')))
  p1 <- p1 + rasterize(geom_point(alpha=.7),dpi=400) + facet_grid(GWAS_type_full~1) + scale_fill_manual(values=color_COVID)
  p1 <- p1 + scale_color_manual(values=color_COVID) + scale_shape_manual(values=c("pos"=24,"neg"=25))
  p1 <- p1 + geom_vline(col='grey',xintercept=myPOS)
#COVID_stats=fread(sprintf('%s/%s_%s_COVID.tsv.gz',COLOC_DIR,SYMBOL,SNP))
  saveRDS(p1,file=sprintf('%s/%s_%s_COVID_plot.RDS',COLOC_DIR,SYMBOL,SNP))
  pdf(sprintf('%s/%s_%s_COVID.pdf',COLOC_DIR,SYMBOL,SNP))
  print(p1)
  dev.off()
}

cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')
if(!is.null(myCELLTYPE)){
  cellstates=cellstates[gsub('(.*)__(.*)','\\1',cellstates)==myCELLTYPE]
}
if(!is.null(mySTATE)){
  cellstates=cellstates[gsub('(.*)__(.*)','\\2',cellstates)==mySTATE]
}

eQTL_stats=list()
for (CELLTYPE__STATE in cellstates){
  cat(CELLTYPE__STATE,'')
  cmd=sprintf('gunzip -c %s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/%s/%s/eQTL_ALL_%s_assoc.txt.gz | grep %s',EVO_IMMUNO_POP_ZEUS,RUN_NAME,CELLTYPE__STATE,CHR,GENE)
  eQTL_stats[[CELLTYPE__STATE]]=fread(cmd)
}

eQTL_stats=rbindlist(eQTL_stats,idcol='cellstate')
setnames(eQTL_stats,c('V1','V2','V3','V4','V5','V6','V7','V8','V9'),c('snps','gene','statistic','pval','beta','celltype','R2','sebeta','CisDist'),skip_absent=TRUE)
eQTL_stats[,state:=gsub('(.*)__(.*)','\\2',cellstate)]
mm=match(eQTL_stats$snps,SNP_info$ID)
eQTL_stats[,POS:=SNP_info[mm,POS]]
eQTL_stats=eQTL_stats[POS>myPOS-WINDOW_SIZE & POS<myPOS+WINDOW_SIZE,]

if(length(cellstates) <= 15 & is.null(mySTATE) & is.null(myCELLTYPE)){
   p2 <- ggplot(eQTL_stats,aes(x=POS,y=-log10(pval),col=state,fill=state,shape=ifelse(beta>0,'pos','neg')))+facet_grid(celltype~1)
   p2 <- p2 + rasterize(geom_point(alpha=.7),dpi=400)+ scale_fill_manual(values=color_conditions)
   p2 <- p2 + scale_color_manual(values=color_conditions) + scale_shape_manual(values=c('pos'=24,'neg'=25))
   p2 <- p2 + geom_vline(col='grey',xintercept=myPOS)

   STATE_char='allCond'
   CELLTYPE_char='allCT'
   pname=sprintf("%s_%s_%s_%s",SYMBOL,SNP,STATE_char,CELLTYPE_char)
   saveRDS(p2,file=sprintf('%s/%s_eQTL_plot.RDS',COLOC_DIR,pname))
   pdf(sprintf('%s/%s_eQTL_plot.pdf',COLOC_DIR,pname))
       print(p2)
    dev.off()
}else{
   if(is.null(mySTATE)){
     STATE_char='allCond'
     }else{
       STATE_char=mySTATE
  }
  celltypes=unique(gsub('(.*)__(.*)','\\1',cellstates))
  for( myCELLTYPE in celltypes ){
     p2 <- ggplot(eQTL_stats[celltype==myCELLTYPE,],aes(x=POS,y=-log10(pval),col=state,fill=state,shape=ifelse(beta>0,'pos','neg')))+facet_grid(state~celltype)
     p2 <- p2 + rasterize(geom_point(alpha=.7),dpi=400)+ scale_fill_manual(values=color_conditions)
     p2 <- p2 + scale_color_manual(values=color_conditions) + scale_shape_manual(values=c('pos'=24,'neg'=25))
     p2 <- p2 + geom_vline(col='grey',xintercept=myPOS)
     CELLTYPE_char=myCELLTYPE
     pname=sprintf("%s_%s_%s_%s",SYMBOL,SNP,STATE_char,CELLTYPE_char)
     saveRDS(p2,file=sprintf('%s/%s_eQTL_plot.RDS',COLOC_DIR,pname))
     pdf(sprintf('%s/%s_eQTL_plot.pdf',COLOC_DIR,pname))
      print(p2)
     dev.off()
   }
}

pname=sprintf("%s_%s_allCond_allCT",SYMBOL,SNP)
fwrite(eQTL_stats,file=sprintf('%s/%s_eQTL.tsv.gz',COLOC_DIR,pname),sep='\t')
#eQTL_stats=fread(sprintf('%s/%s_eQTL.tsv.gz',COLOC_DIR,pname))
#COVID_stats=fread(sprintf('%s/%s_%s_COVID.tsv.gz',COLOC_DIR,SYMBOL,SNP))

MinCell_perCOND=500
keptIID=fread(sprintf('%s/single_cell/project/pop_eQTL/data/1_dataset_description/keptIID_moreThan%scells.tsv',EVO_IMMUNO_POP_ZEUS,MinCell_perCOND),sep='\t',header=F)[,V1]

allTRAITs=TRAIT
if(TRAIT=='allTraits'){
  allTRAITs=c('reported','hospitalized','critical')
  allTRAITs_char="allTraits"
}

Genotypes=queryRange(CHR,myPOS-WINDOW_SIZE,myPOS+WINDOW_SIZE)
Genotypes=melt(Genotypes,id.vars=c("ID"),variable.name='IID')[IID%in%keptIID,]
setkey(Genotypes,ID,IID)
Genotypes[,value:=as.numeric(value)]

coloc.summary=list()
coloc_eQTLstats=list()
coloc_COVIDstats=list()
for(TRAIT in allTRAITs){
  snp_set_COVID=COVID_stats[GWAS_type_full==TRAIT & !is.na(beta),ID]
  for (myCT_S in cellstates){
    cat('\ntesting coloc',TRAIT,myCT_S,'\n')
    snp_set_eQTL=eQTL_stats[cellstate==myCT_S,snps]
    snp_set=intersect(snp_set_COVID,snp_set_eQTL)

    Geno_mat=dcast(Genotypes[ID%chin%snp_set,],IID~ID)
    Geno_res=apply(as.matrix(Geno_mat[,-'IID']),2,function(x){lm(x~substr(Geno_mat$IID,1,3))$res})
    LDMAP=cor(Geno_res)

    eQTL_stats_locus=eQTL_stats[cellstate==myCT_S & snps%chin%snp_set,][order(POS),]
    COVID_stats_locus=COVID_stats[GWAS_type_full==TRAIT & ID%chin%snp_set,][order(POS),]
    eQTL_stats_locus[is.na(sebeta),sebeta:=1e-7]

    compare_beta=merge(COVID_stats,eQTL_stats,by.x=c('ID','POS'),by.y=c('snps','POS'),suffix=c('.covid','.eQTL'),allow.cartesian=TRUE)

    snp_set=eQTL_stats_locus[,snps]

    data_eQTL=list(beta=setNames(eQTL_stats_locus[,beta],eQTL_stats_locus[,snps]),
                varbeta=setNames(eQTL_stats_locus[,sebeta^2],eQTL_stats_locus[,snps]),
                snp=eQTL_stats_locus[,snps],
                position=eQTL_stats_locus[,POS],
                LD=LDMAP[snp_set,snp_set],
                type='quant',
                sdY=1)

    data_COVID=list(beta=setNames(COVID_stats_locus[,beta],COVID_stats_locus[,ID]),
                varbeta=setNames(COVID_stats_locus[,sebeta^2],COVID_stats_locus[,ID]),
                snp=COVID_stats_locus[,ID],
                position=COVID_stats_locus[,POS],
                LD=LDMAP[snp_set,snp_set],
                type='quant',
                sdY=1)

    coloc.res=coloc.signals(data_eQTL,data_COVID,p12=1e-5)
    coloc.res$summary$run=RUN_NAME
    coloc.summary[[paste(TRAIT,myCT_S)]]=coloc.res$summary
    coloc.summary[[paste(TRAIT,myCT_S)]][,trait:=TRAIT]
    coloc.summary[[paste(TRAIT,myCT_S)]][,celltype:=gsub('(.*)__(.*)','\\1',myCT_S)]
    coloc.summary[[paste(TRAIT,myCT_S)]][,state:=gsub('(.*)__(.*)','\\2',myCT_S)]
    coloc.summary[[paste(TRAIT,myCT_S)]][,gene:=GENE]
    coloc.summary[[paste(TRAIT,myCT_S)]][,symbol:=SYMBOL]
    fwrite(coloc.res$summary,file=sprintf('%s/%s_%s_%s_%s_coloc_signal.tsv',COLOC_DIR,SYMBOL,SNP,TRAIT,myCT_S),sep='\t')

    targetLine=coloc.summary[[paste(TRAIT,myCT_S)]]
    RUN_NAME=targetLine[,run]
    CT=targetLine[,celltype]
    STATE=targetLine[,state]
    BEST1=targetLine[,best1] #
    BEST2=targetLine[,best2]
    BEST4=targetLine[,best4]
    HIT1=targetLine[,hit1]
    HIT2=targetLine[,hit2]

  coloc_eQTLstats[[paste(TRAIT,myCT_S)]]=eQTL_stats[celltype==CT & state==STATE,][match(c(SNP,BEST1,BEST2,BEST4),snps),]
  coloc_eQTLstats[[paste(TRAIT,myCT_S)]][,type:=c('SNP_eQTL','BEST1_eQTL','BEST2_COVID','BEST4_COLOC')]
  coloc_eQTLstats[[paste(TRAIT,myCT_S)]]=coloc_eQTLstats[[paste(TRAIT,myCT_S)]][,.(type,snps,celltype,state,trait=TRAIT,beta_eQTL=beta,
                                                                                                      beta_eQTL.lowerCI=sign(beta)*(abs(beta)-2*sebeta), beta_eQTL.upperCI=sign(beta)*(abs(beta)+2*sebeta),
                                                                                                      pvalue_eQTL=pval)]
  fwrite(coloc_eQTLstats[[paste(TRAIT,myCT_S)]],file=sprintf('%s/%s_%s_%s_%s_coloc_signal_eQTLstats.tsv',COLOC_DIR,SYMBOL,SNP,TRAIT,myCT_S),sep='\t')

  coloc_COVIDstats[[paste(TRAIT,myCT_S)]]=COVID_stats[GWAS_type_full==TRAIT,][match(c(SNP,BEST1,BEST2,BEST4),ID),]
  coloc_COVIDstats[[paste(TRAIT,myCT_S)]][,type:=c('SNP_eQTL','BEST1_eQTL','BEST2_COVID','BEST4_COLOC')]
  coloc_COVIDstats[[paste(TRAIT,myCT_S)]]=coloc_COVIDstats[[paste(TRAIT,myCT_S)]][,.(type,ID,celltype=CT,state=STATE,trait=GWAS_type_full,beta_covid=beta,
                                                                                                              beta_covid.lowerCI=sign(beta)*(abs(beta)-2*sebeta), beta_covid.upperCI=sign(beta)*(abs(beta)+2*sebeta),
                                                                                                              pvalue_covid=pval)]
  fwrite(coloc_COVIDstats[[paste(TRAIT,myCT_S)]],file=sprintf('%s/%s_%s_%s_%s_coloc_signal_Covidstats.tsv',COLOC_DIR,SYMBOL,SNP,TRAIT,myCT_S),sep='\t')
  }
}
coloc_COVIDstats=rbindlist(coloc_COVIDstats)
fwrite(coloc_eQTLstats,file=sprintf('%s/allTraits_allStates_coloc_signal_eQTLstats.tsv',COLOC_DIR),sep='\t')
coloc_eQTLstats=rbindlist(coloc_eQTLstats)
fwrite(coloc_COVIDstats,file=sprintf('%s/allTraits_allStates_coloc_signal_Covidstats.tsv',COLOC_DIR),sep='\t')
coloc.summary=rbindlist(coloc.summary)
fwrite(coloc.summary,file=sprintf('%s/allTraits_allStates_coloc_signal.tsv',COLOC_DIR),sep='\t')

coloc.summary=merge(coloc.summary,dcast(coloc_COVIDstats,celltype+state+trait~type,value.var=c('beta_covid','beta_covid.lowerCI','beta_covid.upperCI','pvalue_covid')),by=c('trait','state','celltype'))
coloc.summary=merge(coloc.summary,dcast(coloc_eQTLstats,celltype+state+trait~type,value.var=c('beta_eQTL','beta_eQTL.lowerCI','beta_eQTL.upperCI','pvalue_eQTL')),by=c('trait','state','celltype'))
fwrite(coloc.summary,file=sprintf('%s/allTraits_allStates_coloc_signal_detail.tsv.gz',COLOC_DIR),sep='\t')

#dcast(coloc_COVIDstats,trait~type,value.var=c('beta_covid','beta_covid.lowerCI','beta_covid.upperCI','pvalue_covid'))
#dcast(coloc_eQTLstats,celltype+state~type,value.var=c('beta_eQTL','beta_eQTL.lowerCI','beta_eQTL.upperCI','pvalue_eQTL'))


  #### For a chosen eQTL/condition: show a direct comparison of betas, and a plot colored by R2 with the peaks eQTL SNP
  # S_eQTL=try(runsusie(data_eQTL,repeat_until_convergence=FALSE))
  # S_COVID=try(runsusie(data_eQTL,repeat_until_convergence=FALSE))
  # if(class(S_eQTL)=='try-error' | class(S_COVID)=='try-error'){
  #   coloc.res=list(summary=data.table(nsnps=length(snp_set), hit1=NA, hit2=NA, PP.H0.abf=1, PP.H1.abf=0, PP.H2.abf=0, PP.H3.abf=0,PP.H4.abf=0,idx1=NA, idx2=NA,converged=FALSE))
  #   }else{
  #     coloc.res=coloc.susie(S_eQTL,S_COVID)
  #     if(is.null(coloc.res$summary)){
  #       coloc.res=list(summary=data.table(nsnps=length(snp_set), hit1=NA, hit2=NA, PP.H0.abf=1, PP.H1.abf=0, PP.H2.abf=0, PP.H3.abf=0,PP.H4.abf=0,idx1=NA, idx2=NA,converged=TRUE))
  #       }else{
          # coloc.res$summary$converged=TRUE
  #       }
  #     }
  # fwrite(coloc.res$summary,file=sprintf('%s/%s_%s_%s/%s_%s_%s_%s_coloc.tsv',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP,myCT_S),sep='\t')
  #### run another colocalization technique
  #p3 <- ggplot(compare_beta[pval.covid<1e-3 | pval.eQTL<1e-3,],aes(x=beta.eQTL,y=beta.covid,col=state,fill=state))+facet_grid(celltype~GWAS_type_full)
  # p3 <- ggplot(compare_beta[GWAS_type_full==TRAIT,],aes(x=beta.eQTL/sebeta.eQTL,y=beta.covid/sebeta.covid,col=state,fill=state))+facet_grid(celltype~state)
  # p3 <- p3 + rasterize(geom_point(alpha=.7),dpi=400)+ scale_fill_manual(values=color_conditions)
  # p3 <- p3 + scale_color_manual(values=color_conditions) + scale_shape_manual(values=c('pos'=24,'neg'=25))
  # p3 <- p3 + geom_smooth(method='lm',formula=y~0+x,col='black',fill='black') + geom_hline(col='grey',yintercept=0)+ geom_vline(col='grey',xintercept=0)
  # p3 <- p3 + geom_hline(col='lightgrey',yintercept=3,linetype="dashed")+ geom_vline(col='grey',xintercept=3,linetype="dashed")
  # p3 <- p3 + geom_hline(col='lightgrey',yintercept=-3,linetype="dashed")+ geom_vline(col='grey',xintercept=-3,linetype="dashed")
  #   saveRDS(p3,file=sprintf('%s/%s_%s_%s/%s_%s_%s_COVID_eQTL_betacompare_plot.RDS',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
  #   pdf(sprintf('%s/%s_%s_%s/%s_%s_%s_COVID_eQTL_betacompare_plot.pdf',FIGURE_DIR,SYMBOL,TRAIT,SNP,SYMBOL,TRAIT,SNP))
  #   print(p3)
  #   dev.off()

############################# check results for a given SNP_Symbol_RUN
# SYMBOL='MUC20'
# for(TRAIT in allTRAITs){
#   for (myCT_S in cellstates){
#   COLOC_DIR=sprintf('%s/%s/Dist_%s/%s_%s/',FIGURE_DIR,RUN_NAME,CIS_DIST_TEXT,SYMBOL,SNP)
#   pname=sprintf("%s_%s_allCond_allCT",SYMBOL,SNP)
#   eQTL_stats=fread(sprintf('%s/%s_eQTL.tsv.gz',COLOC_DIR,pname))
#   COVID_stats=fread(sprintf('%s/%s_%s_COVID.tsv.gz',COLOC_DIR,SYMBOL,SNP))
#
# targetLine=all_results[hit1!=best1 & PP.H4.abf>.8 & symbol==SYMBOL,][1,]
#
# RUN_NAME=targetLine[,run]
# SNP=targetLine[,eQTL_snp]
# CT=targetLine[,celltype]
# STATE=targetLine[,state]
# TRAIT=targetLine[,trait]
# BEST1=targetLine[,best1]
# BEST2=targetLine[,best2]
# BEST4=targetLine[,best4]
# BEST4=targetLine[,best4]
# HIT1=targetLine[,hit1]
# HIT2=targetLine[,hit2]
#
# COLOC_DIR=sprintf('%s/%s/Dist_%s/%s_%s/',FIGURE_DIR,RUN_NAME,CIS_DIST_TEXT,SYMBOL,SNP)
#
# pname=sprintf("%s_%s_allCond_allCT",SYMBOL,SNP)
# eQTL_stats=fread(sprintf('%s/%s_eQTL.tsv.gz',COLOC_DIR,pname))
# COVID_stats=fread(sprintf('%s/%s_%s_COVID.tsv.gz',COLOC_DIR,SYMBOL,SNP))
#
# coloc.summary=fread(sprintf('%s/allTraits_allStates_coloc_signal.tsv',COLOC_DIR),sep='\t')
#
# eQTL_stats[celltype==CT & state==STATE,][match(c(SNP,HIT1,HIT2,BEST1,BEST2,BEST4),snps),][,type:=c('SNP_eQTL','HIT1_eQTL','HIT2_COVID','BEST1_eQTL','BEST2_COVID','BEST4_COLOC')][1:.N]
# COVID_stats[GWAS_type_full==TRAIT,][match(c(SNP,HIT1,HIT2,BEST1,BEST2,BEST4),ID),][,type:=c('SNP_eQTL','HIT1_eQTL','HIT2_COVID','BEST1_eQTL','BEST2_COVID','BEST4_COLOC')][1:.N]
