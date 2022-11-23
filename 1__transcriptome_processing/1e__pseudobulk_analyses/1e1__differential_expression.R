################################################################################
################################################################################
# File name: 1e1__differential expression.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Various differential expression/response tests
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./1a__quality_control__lib.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))

# define default parameters

NLIBS=125
CELLTYPE='celltype' # celltype variable to use. Will be used for naming of output files
STATE='condition'

# load batch-adjusted counts computed in 1c2__pseudobulk_batch_correction.R
Expr_ct=fread(sprintf("1__transcriptome_processing/data/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz",NLIBS,CELLTYPE,STATE))
CELLTYPE='lineage'
Expr_ln=fread(sprintf("1__transcriptome_processing/data/adjusted_pseudobulk_%slibs__per_%s_%s_IID.tsv.gz",NLIBS,CELLTYPE,STATE))

################################################################################
# Differential expression between conditions: lineage-level

# reshape batch-adjusted pseudobulk data
Expr_test=dcast(Expr_ln,IID+celltype+Symbol+ID~state,value.var="logCPM")

# difference in means, Wilcoxon's p-value by gene, lineage and condition
DE_COND=Expr_test[,.(logFC_COV=mean(COV)-mean(NS),
                  pval_COV=wilcox.test(COV-NS,data=.SD)$p.value,
                  logFC_IAV=mean(IAV)-mean(NS),
                  pval_IAV=wilcox.test(IAV-NS,data=.SD)$p.value)
                  ,by=.(ID,Symbol,celltype)]

DE_COND[,FDR_COV:=p.adjust(pval_COV,"fdr")]
DE_COND[,FDR_IAV:=p.adjust(pval_IAV,"fdr")]

# keep only significantly differentially expressed genes
logfcthr=0.5
DE_COND[grepl('ENSG', ID) & ((FDR_COV<0.01 & logFC_COV>logfcthr)|(FDR_IAV<0.01 & logFC_IAV>logfcthr)),length(unique(Symbol))]

fwrite(DE_COND,"1__transcriptome_processing/data/DE_condition_by_lineage.tsv.gz",sep="\t")

################################################################################
# Differential expression between conditions: cell type-level

# reshape batch-adjusted pseudobulk data
Expr_test=dcast(Expr_ct,IID+celltype+Symbol+ID~state,value.var="logCPM")

# difference in means, Wilcoxon's p-value by gene, cell type and condition
DE_COND=Expr_test[,.(logFC_COV=mean(COV)-mean(NS),
                  pval_COV=wilcox.test(COV-NS,data=.SD)$p.value,
                  logFC_IAV=mean(IAV)-mean(NS),
                  pval_IAV=wilcox.test(IAV-NS,data=.SD)$p.value)
                  ,by=.(ID,Symbol,celltype)]

DE_COND[,FDR_COV:=p.adjust(pval_COV,"fdr")]
DE_COND[,FDR_IAV:=p.adjust(pval_IAV,"fdr")]

# keep only significantly differentially expressed genes
logfcthr=0.5
DE_COND[grepl('ENSG', ID) & ((FDR_COV<0.01 & logFC_COV>logfcthr)|(FDR_IAV<0.01 & logFC_IAV>logfcthr)),length(unique(Symbol))]

fwrite(DE_COND,"1__transcriptome_processing/data/DE_condition_by_celltype.tsv.gz",sep="\t")

################################################################################
# Differential expression between memory and non-memory NK compartments

query="NK.M.LIKE"
de__nk_mem=lapply(c("NK.CD56dim","NK.CD56brt"),function(ref){
  NK=Expr_ct[grepl(sprintf("%s|%s",query,ref),celltype),]
  NK[,memory:=ifelse(celltype==query,1,0)]
  genes=Expr_ct[,unique(Symbol)]

  round3=lapply(c("NS","COV","IAV"),function(cn){
    round2=lapply(genes,function(g){
      data=NK[POP%in%c("AFB","EUB")&state==cn&Symbol==g]

      res_lm=lm(logCPM~memory+Age+POP+Mortality,data=data)
      ind_order=data[memory==1,IID]
      res_ttest=t.test(data[memory==1,logCPM],data[memory==0,setNames(logCPM,IID)][ind_order],paired=T)
      res_wilcox=wilcox.test(data[memory==1,logCPM],data[memory==0,setNames(logCPM,IID)][ind_order],paired=T)
      round1=data.table(
	Symbol=g,
	type="DE",
	beta.lm=summary(res_lm)$coefficients["memory","Estimate"],
	std_err.lm=summary(res_lm)$coefficients["memory","Std. Error"],
	statistic.lm=summary(res_lm)$coefficients["memory","t value"],
	pval.lm=summary(res_lm)$coefficients["memory","Pr(>|t|)"],
	beta.tt=res_ttest$estimate,
	statistic.tt=res_ttest$statistic,
	pval.tt=res_ttest$p.value,
	statistic.wx=res_wilcox$statistic,
	pval.wx=res_wilcox$p.value
      )
      return(round1)
    })%>%rbindlist()
    print(sprintf("Done: %s.",cn))
    round2[,FDR.lm:=p.adjust(pval.lm,method="hochberg")]
    round2[,FDR.tt:=p.adjust(pval.tt,method="hochberg")]
    round2[,FDR.wx:=p.adjust(pval.wx,method="hochberg")]
    round2[,state:=cn]
    round2[,state:=cn]
  })%>%rbindlist()
  round3[,background:=ref]
  return(round3)
})%>%rbindlist()

################################################################################
# Differential response between memory and non-memory NK compartments

# compute logFC
Expr_ct_NS=Expr_ct[state=="NS",]
if(CELLTYPE=="celltype"){
  Expr_ct_NS_CD14.INFECTED=Expr_ct_NS[celltype=='MONO.CD14',]
  Expr_ct_NS_CD14.INFECTED[,celltype:='MONO.CD14.INFECTED']
  Expr_ct_NS=rbind(Expr_ct_NS,Expr_ct_NS_CD14.INFECTED)
  }
Resp_ct=Expr_ct[state!="NS",]
# obtain logFC
Resp_ct=merge(Resp_ct,Expr_ct_NS,by=c('IID','celltype','ID','Symbol','POP'),suffix=c('','.NS'))
Resp_ct[,logFC:=logCPM-logCPM.NS]

query="NK.M.LIKE"
dr__nk_mem=lapply(c("NK.CD56dim","NK.CD56brt"),function(ref){
  NK=Resp_ct[grepl(sprintf("%s|%s",query,ref),celltype),]
  NK[,memory:=ifelse(celltype==query,1,0)]
  genes=Resp_ct[,unique(Symbol)]

  round3=lapply(c("COV","IAV"),function(cn){
    round2=lapply(genes,function(g){
      data=NK[POP%in%c("AFB","EUB")&state==cn&Symbol==g]

      res_lm=lm(logCPM~memory+Age+POP+Mortality,data=data)
      ind_order=data[memory==1,IID]
      res_ttest=t.test(data[memory==1,logCPM],data[memory==0,setNames(logCPM,IID)][ind_order],paired=T)
      res_wilcox=wilcox.test(data[memory==1,logCPM],data[memory==0,setNames(logCPM,IID)][ind_order],paired=T)
      round1=data.table(
	Symbol=g,
	type="DE",
	beta.lm=summary(res_lm)$coefficients["memory","Estimate"],
	std_err.lm=summary(res_lm)$coefficients["memory","Std. Error"],
	statistic.lm=summary(res_lm)$coefficients["memory","t value"],
	pval.lm=summary(res_lm)$coefficients["memory","Pr(>|t|)"],
	beta.tt=res_ttest$estimate,
	statistic.tt=res_ttest$statistic,
	pval.tt=res_ttest$p.value,
	statistic.wx=res_wilcox$statistic,
	pval.wx=res_wilcox$p.value
      )
      return(round1)
    })%>%rbindlist()
    print(sprintf("Done: %s.",cn))
    round2[,FDR.lm:=p.adjust(pval.lm,method="hochberg")]
    round2[,FDR.tt:=p.adjust(pval.tt,method="hochberg")]
    round2[,FDR.wx:=p.adjust(pval.wx,method="hochberg")]
    round2[,state:=cn]
    round2[,state:=cn]
  })%>%rbindlist()
  round3[,background:=ref]
  return(round3)
})%>%rbindlist()

de__nk_mem[,type:="DE"]
dr__nk_mem[,type:="DR"]
der__nk_mem <- rbind(de__nk_mem,dr__nk_mem)

fwrite(der__nk_mem,"1__transcriptome_processing/data/NK.M.LIKE_markers.tsv",sep="\t")
