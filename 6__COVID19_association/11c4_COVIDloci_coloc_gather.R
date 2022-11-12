
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
EIP=EVO_IMMUNO_POP_ZEUS
FIGURE_DIR=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3/',EVO_IMMUNO_POP_ZEUS)
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
source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",EVO_IMMUNO_POP_ZEUS))

RUN_NAME="/lineage_condition_logFC___CellPropLineage_SVs_220409"
CIS_DIST=1e5

CIS_DIST_TEXT=ifelse(CIS_DIST<1e6, paste0(CIS_DIST/1000,'kb'),paste0(CIS_DIST/1e6,'Mb'))
#dir.create(sprintf('%s/%s/dist_%s',OUT_DIR,RUN_NAME,CIS_DIST_TEXT))
cellstates=dir(sprintf('%s/%s',eQTL_DIR,RUN_NAME),pattern='.*__(NS|COV|IAV|RESTING|ACTIVE)')

feature_toUse=fread(sprintf('%s/genes_to_use.tsv',DATA_DIR),header=F)$V1

TRAIT=NULL
myCELLTYPE=NULL
mySTATE=NULL



RUNS=c('celltype_condition___CellPropLineage_SVs_220409','lineage_condition___CellPropLineage_SVs_220409','lineage_condition_logFC__logFC__CellPropLineage_SVs_220409','celltype_condition_logFC__logFC__CellPropLineage_SVs_220409')

all_results=list()
failed=c()
for (RUN_NAME in RUNS){
  for (CIS_DIST_TEXT in c('100kb','250kb')){ #,'250kb')){
    COLOC_tests=dir(sprintf('%s/colocalization_COVID/%s/Dist_%s/',FIGURE_DIR,RUN_NAME,CIS_DIST_TEXT))
    cat('\n',RUN_NAME,':',length(COLOC_tests),'tests\n')
    counter=0
    for (i in COLOC_tests){
      counter=counter+1
      ID=paste(RUN_NAME,CIS_DIST_TEXT,i,sep='_')
      all_results[[ID]]=try(fread(sprintf('%s/colocalization_COVID/%s/Dist_%s/%s/allTraits_allStates_coloc_signal_detail.tsv.gz',FIGURE_DIR,RUN_NAME,CIS_DIST_TEXT,i)))
      if('try-error'%in%class(all_results[[ID]])){
        all_results[[ID]]=NULL
        failed=c(ID,failed)
      }
      cat(counter,'')
    }
  }
}
all_results=rbindlist(all_results,idcol='ID')
all_results[,Dist:=gsub('(.*)_220409_([0-9kM]+b)_(.*)_(.*)','\\2',ID)]
all_results[,eQTL_snp:=gsub('(.*)_220409_([0-9kM]+b)_(.*)_(.*)','\\4',ID)]

fwrite(all_results,file=sprintf('%s/colocalization_COVID/allTraits_allStates_all_eQTLs_near_1e5_covid_peak_coloc_signal_v2.tsv.gz',FIGURE_DIR),sep='\t')

all_results=fread(sprintf('%s/colocalization_COVID/allTraits_allStates_all_eQTLs_near_1e5_covid_peak_coloc_signal_v2.tsv.gz',FIGURE_DIR))

all_results[,.N,by=run]

  all_results[,.N,by=.(run,state)]
  all_results[,.(.N, mean(PP.H4.abf[abs(hit1.margz)>4 & abs(hit1.margz)>4]>.8)),by=.(run,state)]
  all_results[run=='lineage_condition_logFC___CellPropLineage_SVs_220409' & state=='COV' & PP.H4.abf>.8][order(-PP.H4.abf),][!duplicated(symbol),][,.(trait,celltype,state,symbol,best4,hit1.margz, hit2.margz,run)]
  all_results[run=='lineage_condition___CellPropLineage_SVs_220409' & state=='COV' & PP.H4.abf>.8][order(-PP.H4.abf),][!duplicated(symbol),][,.(trait,celltype,state,symbol,best4,hit1.margz, hit2.margz,run)]
  all_results[run=='lineage_condition___CellPropLineage_SVs_220409' & state=='COV' & PP.H4.abf>.8][order(-PP.H4.abf),][!duplicated(symbol),][,.(trait,celltype,state,symbol,best4,hit1.margz, hit2.margz,run)]

  Signif_coloc=all_results[PP.H4.abf>.8 & Dist=='250kb'][order(-PP.H4.abf),]
  Signif_coloc[,colocID:=paste(symbol,eQTL_snp,trait,celltype, state)]
  Signif_coloc=Signif_coloc[!duplicated(colocID),][,.(trait,celltype,state,gene_name=symbol,eQTL_snp,best4, eQTL_effect_ALT=beta_eQTL_SNP_eQTL, effect_expressionIncreasingAllele_onCOVID=beta_covid_BEST4_COLOC/beta_eQTL_BEST4_COLOC, PP.H4.abf,run)]

  eQTL_effects=all_results[order(-abs(beta_eQTL.lowerCI_SNP_eQTL))][pvalue_eQTL_SNP_eQTL<.01, .(eQTL_effect_ALT_allTissues=paste(unique(paste0(celltype,', ',state,ifelse(beta_eQTL_SNP_eQTL>0,' (+)',' (-)'))),collapse=' // ')),by=.(eQTL_snp,gene_name=symbol)]
  Trait_effects=all_results[PP.H4.abf>.8 & Dist=='250kb'][order(-PP.H4.abf),][order(pvalue_covid_BEST4_COLOC),][pvalue_covid_BEST4_COLOC<.01, .(covid_effect_increased_expression=paste(unique(paste0(trait,
                                      ifelse(sign(beta_eQTL_BEST4_COLOC)*sign(beta_covid_BEST4_COLOC)>0,
                                      ifelse(pvalue_covid_BEST4_COLOC<1e-5,' (++)',' (+)'),
                                      ifelse(pvalue_covid_BEST4_COLOC<1e-5,' (--)',' (-)')))),collapse=' // ')),by=.(eQTL_snp,gene_name=symbol,celltype, state)]

# trait (reported, hospitalized, or critcal covid-19 case) where the colocalization is detected
# celltype/state where the colocalization is detected
# eQTL index SNP,
# best4 : best consensus SNP:
# effect of eQTL index SNP on gene expression
# beta_covid_BEST4_COLOC/beta_eQTL_BEST4_COLOC effect of 2 fold increase in gene expression on Covid Risk
# celltype/state where the eQTL is active with effect of ALT allelel on gene expression, celltype/conditions are orderebd by decreasing lower bound of eQTL effect size
# trait associated with the locus and effect of increased expression on COVID risk

  Signif_coloc_eQTL=merge(Signif_coloc,eQTL_Signif_both[,.(eQTL_snp=snps,gene_name,celltype,state)],by=c('eQTL_snp','gene_name','celltype','state'))
  Signif_coloc_eQTL[,mapping:=paste0(ifelse(grepl('celltype',run),'celltype','lineage'),ifelse(grepl('logFC',run),'_logFC','_expr'))]
  Signif_coloc_eQTL[,run:=NULL]
  Signif_coloc_eQTL=merge(Signif_coloc_eQTL,Trait_effects,by=c('eQTL_snp','gene_name','celltype','state'))
  Signif_coloc_eQTL=merge(Signif_coloc_eQTL,eQTL_effects,by=c('eQTL_snp','gene_name'))

fwrite(Signif_coloc_eQTL[order(-PP.H4.abf),],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V11/TableS9/TableS9A_coloc_results_full.tsv',EVO_IMMUNO_POP_ZEUS),sep='\t')

fwrite(Signif_coloc_eQTL[order(-PP.H4.abf),][!duplicated(paste(gene_name,eQTL_snp))],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V11/TableS9/TableS9A_coloc_results.tsv',EVO_IMMUNO_POP_ZEUS),sep='\t')
#76 eQTL for 40 genes

  # COVID_info=melt(unique(SNP_info[,.(ID,covid_A2_beta,covid_B2_beta,covid_C2_beta,covid_A2_pval,covid_B2_pval,covid_C2_pval)]),id.vars='ID')
  # COVID_info[,GWAS_type:=gsub('covid_([ABC2]+)_(beta|pval)','\\1',variable)]
  # COVID_info[,stat:=gsub('covid_([ABC2]+)_(beta|pval)','\\2',variable)]
  # COVID_info[,trait:=case_when(GWAS_type=='A2'~'critical',
  #                                         GWAS_type=='B2'~'hospitalized',
  #                                         GWAS_type=='C2'~'reported',
  #                                         TRUE~'NA')]
  #
  # COVID_info=dcast(unique(COVID_info),trait+GWAS_type+ID~stat)
  # setnames(COVID_info,'ID','snps')
  #
  # Signif_coloc_eQTL_bestTrait=merge(Signif_coloc_eQTL_bestTrait,COVID_info,by=c('trait','snps'),suffix=c('.eQTL','.covid'))
  # Signif_coloc_eQTL_bestTrait=Signif_coloc_eQTL_bestTrait[,.(trait, gene, gene_name, mapping, celltype, state, PP.H4.abf, eQTL_SNP=snps, eQTL_beta_ALT=beta.eQTL, eQTL_pvalue=pvalue,covid_beta_ALT=beta.covid,covid_pvalue=pval, GWAS_peakSNP=best2, candidate_snp_coloc=best4)]
  # fwrite(Signif_coloc_eQTL_bestTrait[order(-PP.H4.abf),],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/suppTables//TableS6/TableS6E_coloc_results.tsv',EVO_IMMUNO_POP_ZEUS))
  # #fwrite(Signif_coloc_eQTL_bestTrait[order(-PP.H4.abf),][!duplicated(gene_name),],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/suppTables/Table_SXX_coloc_results_best_celltype_cond.tsv',EVO_IMMUNO_POP_ZEUS))

  #
  # eQTL_lineage=eQTL_Signif_both[celltype%chin%lineage_5,.(gene,gene_name,snps,CisDist,lineage=celltype,state,beta,R2,pvalue,FDR,pip_top_component)]
  # fwrite(eQTL_lineage[order(pvalue),],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/suppTables//TableS6/TableS6A_eQTL_lineage.tsv',EVO_IMMUNO_POP_ZEUS))
  # eQTL_celltype=eQTL_Signif_both[celltype%chin%celltype_22,.(gene,gene_name,snps,CisDist,lineage=celltype,state,beta,R2,pvalue,FDR,pip_top_component)]
  # fwrite(eQTL_celltype[order(pvalue),],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/suppTables//TableS6/TableS6B_eQTL_celltype.tsv',EVO_IMMUNO_POP_ZEUS))
  #
  #
  # reQTL_lineage=reQTL_Signif_both[celltype%chin%lineage_5,.(gene,gene_name,snps,CisDist,lineage=celltype,state,beta,R2,pvalue,FDR,pip_top_component)]
  # fwrite(reQTL_lineage[order(pvalue),],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/suppTables//TableS6/TableS6C_reQTL_lineage.tsv',EVO_IMMUNO_POP_ZEUS))
  # reQTL_celltype=reQTL_Signif_both[celltype%chin%celltype_22,.(gene,gene_name,snps,CisDist,lineage=celltype,state,beta,R2,pvalue,FDR,pip_top_component)]
  # fwrite(reQTL_celltype[order(pvalue),],file=sprintf('%s/single_cell/project/pop_eQTL/paper_draft/V3/suppTables//TableS6/TableS6D_reQTL_celltype.tsv',EVO_IMMUNO_POP_ZEUS))

  # all_results[PP.H4.abf>.8][order(-PP.H4.abf),][!duplicated(symbol),][,.(trait,celltype,state,symbol,best4,hit1.margz, hit2.margz,run,gene)][paste(celltype,state,gene)%chin%eQTL_Signif_both[,paste(celltype,state,gene)],]
