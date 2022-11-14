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

# define default values
ln="MONO"
cn="NS"
comparison="EUB_AFB" # which populations to compare
run_id="220409" # date at which the script is launched, used to identify the results unambiguously
popdiff=FALSE # focus only on population differences?
exp_resp="DE" # focus on effects on expression (DE) or response (DR)?

# update parameter values based on provided arguments
cmd=commandArgs()
#print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--lineage' | cmd[i]=='-t' ){ln = cmd[i+1]}
	if (cmd[i]=='--comparison' | cmd[i]=='-p' ){comparison = cmd[i+1]}
	if (cmd[i]=='--exp_resp' | cmd[i]=='-e' ){exp_resp = cmd[i+1]}
	if (cmd[i]=='--popdiff' | cmd[i]=='-f' ){popdiff = cmd[i+1]}
	if (cmd[i]=='--state' | cmd[i]=='-s' ){cn = cmd[i+1]}
	if (cmd[i]=='--runid' | cmd[i]=='-d' ){run_id = cmd[i+1]} # ID of the run (eg. date)
}

print(sprintf("Lineage: %s",ln))
print(sprintf("State: %s",cn))
print(sprintf("Comparison: %s",comparison))
print(sprintf("Only popDs?: %s",popdiff))
print(sprintf("DE or DR?: %s",exp_resp))

# load and format expression data and covariates
Mortality=fread("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/project/pop_eQTL/data/2_population_differences/Covariates/lineage_condition_/Covariates__B_COV.tsv.gz")
Mortality=Mortality[!is.na(MORTALITY)&!grepl("^ASH",IID),.(IID,Mortality=MORTALITY)]

if (exp_resp=="DE"){
  Expression=fread(sprintf("%s/data/BatchAdjusted_logCPM_125libs__per_lineage_condition_annotated.tsv.gz",DAT_POPDIFF_DIR))
  Expression=merge(Expression[POP!="ASH",],Mortality,by="IID")
} else if(exp_resp=="DR"){
  Expression=fread(sprintf("%s/data/BatchAdjusted_logFC_125libs__per_lineage_condition_annotated.tsv.gz",DAT_POPDIFF_DIR))
  Expression=merge(Expression[POP!="ASH",],Mortality,by="IID")
}

# consider only population differences in expression or reponse
if (popdiff==T){
  ANALYSE=ifelse(exp_resp=="DE",
    "lineage_condition_AFBEUB__lineage_condition__220409",
    "lineage_condition_AFBEUB__logFC_lineage_condition_logFC__220409")
  popde=fread(sprintf("%s/%s/perm0/popDiffResults_with_mash.tsv.gz",POP_DIR,ANALYSE))
  popde=popde[celltype==ln&state==cn&FDR<0.01&abs(beta)>0.2,]
  Expression=Expression[Symbol%chin%popde[,unique(Symbol)]]
}

# load and add cell type frequency information
celltype_frequencies=fread(sprintf("%s/data/Cellcounts/Pct_celltype_by_IID_condition.tsv.gz",DAT_POPDIFF_DIR))
celltype_frequencies_melt=melt(celltype_frequencies,measure.vars=3:24,variable.name="celltype",value.name="pct")
lineage_celltype=fread(sprintf("%s/data/lineage_celltype.tsv",AGGR_QC_DIR))
celltype_frequencies_melt$lineage=lineage_celltype[,setNames(lineage,celltype)][as.character(celltype_frequencies_melt$celltype)]
pct_sum_lineage=celltype_frequencies_melt[,sum(pct),by=.(IID,state,lineage)]
setnames(pct_sum_lineage,"V1","pct_sum_lineage")
celltype_frequencies_melt=merge(celltype_frequencies_melt,pct_sum_lineage,by=c("IID","state","lineage"))
celltype_frequencies_melt[,pct_std_lineage:=pct/pct_sum_lineage*100]

# extract feature information
feature_metadata=Expression[,.(ID,Symbol)]%>%unique

popq=str_split(comparison,"_",simplify=T)[,1]
popr=str_split(comparison,"_",simplify=T)[,2]
Expression_info=Expression[celltype==ln&(POP==popq|POP==popr),]

Expression_info=merge(
  Expression_info,
  dcast(celltype_frequencies_melt[lineage==ln,],IID+state~celltype,value.var="pct_std_lineage"),
  by=c("IID","state")
)

if(exp_resp=='DR'){
  setnames(Expression_info,'logFC','logCPM',skip_absent=TRUE) # change name to harmonize downstream code
}

ct_in_ln=celltype_frequencies_melt[lineage==ln,as.character(unique(celltype))]
if (ln=="MONO"&cn!="IAV") {
  ct_in_ln=ct_in_ln[-which(ct_in_ln=="MONO.CD14.INFECTED")]
}

################################################################################
# Run mediation analysis

# for each gene and condition, compute the fraction of variance in
# expression/response explained by the cell type with the most
# different frequencies between AFB and EUB in each lineage
result=lapply(ct_in_ln,function(ct){
    round1=lapply(feature_metadata[,unique(Symbol)],function(g){
      data=Expression_info[state==cn&Symbol==g,]
      if (data[,all(logCPM==0)]) {
        res=data.table(Symbol=g,ID=feature_metadata[Symbol==g,ID],frac_var=0)
      } else {
        if (!grepl("ASH",comparison)) {
          mod_1=data[,lm(logCPM~get(ct)+POP+Age+Mortality)]
          mod_2=data[,lm(get(ct)~POP+Age+Mortality)]
        } else if (grepl("ASH",comparison)) {
          mod_1=data[,lm(logCPM~get(ct)+POP+Age+GenderF+Mortality)]
          mod_2=data[,lm(get(ct)~POP+Age+GenderF+Mortality)]
        }
      med=mediate(mod_2,mod_1,treat="POP",mediator="get(ct)",control.value=popr,treat.value=popq)
      res=data.table(Symbol=g,ID=feature_metadata[Symbol==g,ID],frac_var=med$n0,pval=med$n0.p)
      }
      print(sprintf("Done: %s. %s, %s, %s",g,comparison,ct,cn))
      return(res)
    })%>%rbindlist()
    round1[,state:=cn]
    round1[,lineage:=ln]
    round1[,celltype:=ct]
    round1[,comp:=comparison]
    return(round1)
})%>%rbindlist()

fwrite(result,sprintf("%s/celltype_frequencies/%s_%s__%s%s%s_%s.tsv",OUT_DIR,ln,cn,comparison,ifelse(popdiff==T,"_popD","_"),gsub("^D","",exp_resp),run_id),sep="\t")
