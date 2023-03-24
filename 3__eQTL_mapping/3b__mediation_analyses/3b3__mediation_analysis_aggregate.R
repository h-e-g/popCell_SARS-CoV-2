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

# load eQTLs
source(sprintf("./load_eQTLs.R",MISC_DIR))

################################################################################
# Aggregate cell type mediation results

run_id="RUN1"

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
 round4[,exp_resp:=ifelse(er=='DE','expression','response')]
})%>%rbindlist()
mediation_celltype[,frac_var_std:=case_when(frac_var>1~1,frac_var<0~0,T~frac_var)]
mediation_celltype[,FDR:=p.adjust(pval,"fdr")]
setnames(mediation_celltype,"celltype","mediator")
setnames(mediation_celltype,"lineage","celltype")
mediation_celltype[,mediator_type:='celltype']

vars_to_keep=c("Symbol","ID","comparison","exp_resp","popdiff","celltype","condition","mediator","frac_var","frac_var_std","pval","FDR","mediator_type")


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
 round4[,exp_resp:=ifelse(er=='DE','expression','response')]
})%>%rbindlist()
mediation_genetics[,frac_var_std:=case_when(frac_var>1~1,frac_var<0~0,T~frac_var)]
mediation_genetics[,FDR:=p.adjust(pval,"fdr")]
setnames(mediation_genetics,"snps","mediator")
mediation_genetics[,mediator_type:='genetics']

fwrite(mediation_genetics[,mget(vars_to_keep)],sprintf("%s/genetics/mediation_analysis_genetics.tsv",OUT_DIR),sep="\t")

################################################################################
# Aggregate all mediation results
mediation=rbind(mediation_celltype,mediation_genetics)
setnames(mediation_genetics,"expr_resp","type")

fwrite(mediation,sprintf("%s/mediation_analysis_genetics_and_celltype.tsv",OUT_DIR),sep="\t")

# focus on one celltype per lineage
celltype_mediators<-c("MONO.CD16","B.M.K","T.CD4.E","T.CD8.EMRA","NK.M.LIKE")
mediation=mediation[mediator%in%celltype_mediators | mediator_type=='genetics',]
mediation[,celltype:=factor(celltype,lineage_order)]
mediation[,state:=factor(state,cond_order)]

### load popDEGs
ANALYSE="popdiff___lineage_condition__noCellProps_noSVs___220409_perm0.tsv"
popdegs=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popdegs[,type:="expression"]
popdegs<-popdegs[FDR<0.01 & abs(beta)>0.2,]
bin_order=c("VL","L","M","H","VH")
breaks_DE=popdegs[,quantile(abs(beta),c(0,0.2,0.4,0.6,0.8,1))]
popdegs[,bin:=cut(abs(beta),breaks_DE,include.lowest=T,right=F,labels=bin_order)]

### load popDRGs
ANALYSE="popdiff__lineage_condition_logFC__noCellProps_noSVs___220409_perm0.tsv"
popdrgs=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popdrgs[,type:="response"]
popdrgs<-popdrgs[FDR<0.01 & abs(beta)>0.2,]
bin_order=c("VL","L","M","H","VH")
breaks_DR=popdrgs[,quantile(abs(beta),c(0,0.2,0.4,0.6,0.8,1))]
popdrgs[,bin:=cut(abs(beta),break_DRs,include.lowest=T,right=F,labels=bin_order)]

### concatenate them in a single object
pop_degs_drgs=rbind(pop_degs,popdrgs)

# add popDE/popDR to mediation table
mediation$betapop=pop_degs_drgs[,setNames(beta,sprintf("%s%s%s%s",celltype,state,Symbol,type))][mediation[,sprintf("%s%s%s%s",celltype,state,Symbol,type)]]
mediation$bin=pop_degs_drgs[,setNames(bin,sprintf("%s%s%s%s",celltype,state,Symbol,type))][mediation[,sprintf("%s%s%s%s",celltype,state,Symbol,type)]]


eQTL_file=fread(sprintf("%s/2_population_differences/best_eQTL_per_gene_and_comparison_MAF005.tsv",DAT_DIR))
eQTL_file[,type:='expression']
reQTL_file=fread(sprintf("%s/2_population_differences/best_reQTL_per_gene_and_comparison_MAF005.tsv",DAT_DIR))
eQTL_file[,type:='response']

mediation$betaeqtl=eQTL_file[comparison=="EUB_AFB",setNames(beta,sprintf("%s%s%s%s",celltype,state,snps,type))][mediation[,sprintf("%s%s%s",celltype,state,mediator,type)]]
mediation$peqtl=eQTL_file[comparison=="EUB_AFB",setNames(pvalue,sprintf("%s%s%s%s",celltype,state,snps,type))][mediation[,sprintf("%s%s%s",celltype,state,mediator,type)]]

eqtl_signif_info=merge(eQTL_Signif_both,SNP_info[,.(snps=ID,AFB,EUB,ASH)],by="snps",all.x=T)
eqtl_signif_info[,type:='DE']
reqtl_signif_info=merge(reQTL_Signif_both,SNP_info[,.(snps=ID,AFB,EUB,ASH)],by="snps",all.x=T)
reqtl_signif_info[,type:='DR']
eqtl_signif_info=rbind(eqtl_signif_info,reqtl_signif_info)
## compute MAF
eqtl_signif_info[,MAF_AFB:=pmin(AFB,1-AFB)]
eqtl_signif_info[,MAF_EUB:=pmin(EUB,1-EUB)]
eqtl_signif_info[,MAF_ASH:=pmin(ASH,1-ASH)]

mediation$eGene=eqtl_signif_info[celltype%in%lineage_order&(MAF_AFB>0.05|MAF_EUB>0.05),setNames(rep(T,nrow(eQTL_Signif_both)),sprintf("%s%s%s%s",celltype,state,gene_name,type))][mediation[,sprintf("%s%s%s%s",celltype,state,Symbol,type)]]
mediation[,eGene:=ifelse(is.na(eGene),F,eGene)]
mediation[,FDR_global:=p.adjust(pval,"fdr")]

fwrite(mediation[,.(comp,type,celltype,state,Symbol,ID,betapop,mediator_type,mediator,frac_var,pval,FDR=FDR_global,eGene,betaeqtl,peqtl)],file=sprintf("%s/mediation_w_betapop_tableS6.tsv",OUT_DIR))

# mediation_DE=fread(sprintf("%s/2_population_differences/mediation_analysis/mediation_w_betapop_popDE.tsv",DAT_DIR))
# mediation_DR=fread(sprintf("%s/2_population_differences/mediation_analysis/mediation_w_betapop_popDR.tsv",DAT_DIR))
