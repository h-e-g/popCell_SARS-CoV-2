################################################################################
################################################################################
# File name: 3b1__mediation_analysis_celltype.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Use mediation analysis to quantify var. explained by cell proportions
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./3b__mediation_analyses__lib.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))
POP_DIR=DAT_POPDIFF_DIR
OUT_DIR=DAT_MEDIATION_DIR

# declare useful functions
source(sprintf("./misc_functions.R",MISC_DIR))

################################################################################
# Aggregate cell type mediation results

run_id="220617"

mediation_celltype=lapply(c("MONO","B","T.CD4","T.CD8","NK"),function(ln){
round4=lapply(c("DE","DR"),function(er){
   round3=lapply(c("EUB_AFB"),function(cp){
	   if(exp_resp=="DE"){ state_list=c("NS","COV","IAV") }
	   if(exp_resp=="DR"){ state_list=c("COV","IAV") }
	   round2=lapply(state_list,function(cn){
		   round1=lapply(c(T,F),function(pd){
         fread(sprintf("%s/celltype_frequencies/%s_%s__%s%s%s_%s.tsv",OUT_DIR,ln,cn,cp,ifelse(pd==T,"_popD","_"),gsub("^D","",er),run_id))
       })%>%rbindlist()
			 round1[,popdiff:=pd]
		 })%>%rbindlist()
		 round2[,condition:=cn]
	 })%>%rbindlist()
	 round3[,comparison:=cp]
 })%>%rbindlist()
 round4[,exp_resp:=er]
})%>%rbindlist()
mediation_celltype[,frac_var_std:=case_when(frac_var>1~1,frac_var<0~0,T~frac_var)]
mediation_celltype[,FDR:=p.adjust(pval,"fdr")]
setnames(mediation_celltype,"celltype","mediator")
setnames(mediation_celltype,"lineage","celltype")

vars_to_keep=c("Symbol","ID","comparison","exp_resp","popdiff","celltype","condition","mediator","frac_var","frac_var_std","pval","FDR")

fwrite(mediation_celltype[,mget(vars_to_keep)],sprintf("%s/celltype_frequencies/mediation_analysis_celltype.tsv",OUT_DIR),sep="\t")

################################################################################
# Aggregate genetics mediation results

mediation_genetics=lapply(c("MONO","B","T.CD4","T.CD8","NK"),function(ln){
round4=lapply(c("DE","DR"),function(er){
   round3=lapply(c("EUB_AFB"),function(cp){
	   if(exp_resp=="DE"){ state_list=c("NS","COV","IAV") }
	   if(exp_resp=="DR"){ state_list=c("COV","IAV") }
	   round2=lapply(state_list,function(cn){
		   round1=lapply(c(T,F),function(pd){
         fread(sprintf("%s/genetics/%s_%s__%s%s%s_%s.tsv",OUT_DIR,ln,cn,cp,ifelse(pd==T,"_popD","_"),gsub("^D","",er),run_id))
       })%>%rbindlist()
			 round1[,popdiff:=pd]
		 })%>%rbindlist()
		 round2[,condition:=cn]
	 })%>%rbindlist()
	 round3[,comparison:=cp]
 })%>%rbindlist()
 round4[,exp_resp:=er]
})%>%rbindlist()
mediation_genetics[,frac_var_std:=case_when(frac_var>1~1,frac_var<0~0,T~frac_var)]
mediation_genetics[,FDR:=p.adjust(pval,"fdr")]
setnames(mediation_celltype,"celltype","mediator")
setnames(mediation_celltype,"lineage","celltype")

fwrite(mediation_genetics[,mget(vars_to_keep)],sprintf("%s/genetics/mediation_analysis_genetics.tsv",OUT_DIR),sep="\t")

################################################################################
# Aggregate all mediation results

mediation=rbind(mediation_celltype,mediation_genetics)

fwrite(mediation,sprintf("%s/mediation_analysis_genetics.tsv",OUT_DIR),sep="\t")
