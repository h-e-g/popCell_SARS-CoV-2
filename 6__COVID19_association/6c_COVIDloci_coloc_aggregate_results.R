
################################################################################
################################################################################
# File name: 6c_COVIDloci_coloc_aggregate_results.R
# Author: J.MR, Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: aggregate colocalization results from 6b across all eQTLs/reQTLs
# Effector script
################################################################################
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

# declare parameters
CIS_DIST_TEXT='250kb'
RUNS=c('lineage_condition___CellPropLineage_SVs_RUN1',
      'lineage_condition_logFC__CellPropLineage_SVs_RUN1',
      'celltype_condition___CellPropLineage_SVs_RUN1',
      'celltype_condition_logFC__CellPropLineage_SVs_RUN1')

all_results=list()
failed=c()
for (RUN_NAME in RUNS){
    COLOC_tests=dir(sprintf('%s/%s/',COLOC_DIR,RUN_NAME))
    cat('\n',RUN_NAME,':',length(COLOC_tests),'tests\n')
    counter=0
    for (i in COLOC_tests){
      counter=counter+1
      ID=paste(RUN_NAME,i,sep='_')
      all_results[[ID]]=try(fread(sprintf('%s/%s/allTraits_allStates_coloc_signal_detail.tsv.gz',COLOC_DIR,RUN_NAME,i)))
      if('try-error'%in%class(all_results[[ID]])){
        all_results[[ID]]=NULL
        failed=c(ID,failed)
      }
      cat(counter,'')
  }
}
all_results=rbindlist(all_results,idcol='ID')
all_results[,Dist:=gsub('(.*)_RUN1_(.*)_(.*)','\\2',ID)]
all_results[,eQTL_snp:=gsub('(.*)_RUN1_(.*)_(.*)','\\4',ID)]

fwrite(all_results,file=sprintf('%s/colocalization_COVID_all_Traits_eQTL_and_cellStates.tsv.gz',COLOC_DIR),sep='\t')

  Signif_coloc=all_results[PP.H4.abf>.8 & Dist=='250kb'][order(-PP.H4.abf),]
  Signif_coloc[,colocID:=paste(symbol,eQTL_snp,trait,celltype, state)]
  Signif_coloc=Signif_coloc[!duplicated(colocID),]
  Signif_coloc=Signif_coloc[,.(trait,celltype,state,gene_name=symbol,eQTL_snp,
              best4, eQTL_effect_ALT=beta_eQTL_SNP_eQTL, effect_expressionIncreasingAllele_onCOVID=beta_covid_BEST4_COLOC/beta_eQTL_BEST4_COLOC,
              PP.H4.abf,run)]
  eQTL_effects=all_results[order(-abs(beta_eQTL.lowerCI_SNP_eQTL))]
  eQTL_effects=eQTL_effects[pvalue_eQTL_SNP_eQTL<.01,]
  eQTL_effects=eQTL_effects[,.(eQTL_effect_ALT_allTissues=paste(
                  unique( paste0(celltype,',',
                                      state,
                                      ifelse(beta_eQTL_SNP_eQTL>0,' (+)',' (-)')) )
                                                          ,collapse=' // ')),
                                    by=.(eQTL_snp,gene_name=symbol)]

  Trait_effects=all_results[PP.H4.abf>.8 & Dist=='250kb'][order(-PP.H4.abf),][order(pvalue_covid_BEST4_COLOC),][pvalue_covid_BEST4_COLOC<.01, .(covid_effect_increased_expression=paste(unique(paste0(trait,
                                      ifelse(sign(beta_eQTL_BEST4_COLOC)*sign(beta_covid_BEST4_COLOC)>0,
                                      ifelse(pvalue_covid_BEST4_COLOC<1e-5,' (++)',' (+)'),
                                      ifelse(pvalue_covid_BEST4_COLOC<1e-5,' (--)',' (-)')))),collapse=' // ')),by=.(eQTL_snp,gene_name=symbol,celltype, state)]

  Signif_coloc_eQTL=merge(Signif_coloc,eQTL_Signif_both[,.(eQTL_snp=snps,gene_name,celltype,state)],by=c('eQTL_snp','gene_name','celltype','state'))
  Signif_coloc_eQTL[,mapping:=paste0(ifelse(grepl('celltype',run),'celltype','lineage'),ifelse(grepl('logFC',run),'_logFC','_expr'))]
  Signif_coloc_eQTL[,run:=NULL]
  Signif_coloc_eQTL=merge(Signif_coloc_eQTL,Trait_effects,by=c('eQTL_snp','gene_name','celltype','state'))
  Signif_coloc_eQTL=merge(Signif_coloc_eQTL,eQTL_effects,by=c('eQTL_snp','gene_name'))

fwrite(Signif_coloc_eQTL[order(-PP.H4.abf),],file=sprintf('%s/TableS9A_coloc_results_full.tsv',COLOC_DIR),sep='\t')

fwrite(Signif_coloc_eQTL[order(-PP.H4.abf),][!duplicated(paste(gene_name,eQTL_snp))],file=sprintf('%s/TableS9A_coloc_results.tsv',COLOC_DIR),sep='\t')
