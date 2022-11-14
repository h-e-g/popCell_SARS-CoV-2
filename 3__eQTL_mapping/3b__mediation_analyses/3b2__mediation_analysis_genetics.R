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
source(sprintf("./load_eQTLs.R",MISC_DIR))
source(sprintf("./querySNPs.R",MISC_DIR))

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

# load SNP info
SNP_info=getMap(annotate=TRUE)

RUN_EQTL_LINEAGE="lineage_condition___CellPropLineage_SVs_220409"
RUN_REQTL_LINEAGE="lineage_condition_logFC__logFC__CellPropLineage_SVs_220409"
RUN_EQTL_CELLTYPE="celltype_condition___CellPropLineage_SVs_220409"
RUN_REQTL_CELLTYPE="celltype_condition_logFC__logFC__CellPropLineage_SVs_220409"


# load eQTL information, keep common variants (MAF>0.05) in any of the populations compared
 Best_SNP_pergene_MAF005=lapply(c("EUB_AFB","ASH_AFB","ASH_EUB"),function(comp){
  popq=str_split(comp,"_",simplify=T)[,1]
  popr=str_split(comp,"_",simplify=T)[,2]
  round3=lapply(c("MONO","B","T.CD4","T.CD8","NK"),function(lineage){
  cond_to_test=c("NS","COV","IAV")
  if(exp_resp=='DR'){cond_to_test=setdiff(cond_to_test,'NS')}
    round2=lapply(cond_to_test,function(st){
      round1=lapply(1:22,function(ch){
         if(exp_resp=='DE'){RUN=RUN_EQTL_LINEAGE}else{RUN=RUN_REQTL_LINEAGE}
       	x=fread(sprintf("%s/3_eQTL_mapping/%s/%s__%s/FineMapping_100kb/eQTL_FineMapped_chr%s.txt.gz",DAT_DIR,RUN,lineage,st,ch))
       	x=merge(x,SNP_info[,.(snps=ID,AFB,EUB,ASH)],by="snps",all.x=T)
       	# compute MAF
       	x[,MAF_AFB:=pmin(AFB,1-AFB)]
       	x[,MAF_EUB:=pmin(EUB,1-EUB)]
       	x[,MAF_ASH:=pmin(ASH,1-ASH)]
       	x=x[get(sprintf("MAF_%s",popq))>0.05|get(sprintf("MAF_%s",popr))>0.05,]
       	print(sprintf("%s__%s, chr%s, %s",lineage,st,ch,comp))
       	x=x[order(pvalue),][!duplicated(gene),]
       	x[,chr:=ch]
       	x[,state:=gsub("^.+__","",cellstate)]
       	return(x)
      })%>%rbindlist()
      return(round1)
    })%>%rbindlist()
    return(round2)
  })%>%rbindlist()
 round3[,comparison:=comp]
 return(round3)
})%>%rbindlist()

eQTL_file=Best_SNP_pergene_MAF005
setnames(eQTL_file,c("snps","gene"),c("rsID","ID"))

# get genotypes at all these (r)eQTLs for all individuals
genotypes=getSNP(eQTL_file[,unique(rsID)],Map=SNP_info)
genotypes[,n_alt:=round(Number_of_ALT_alelle)]

comp=comparison
eQTL_file=eQTL_file[comparison==comp,]

# extract feature meta data and prepare data
feature_metadata=Expression[,.(ID,Symbol)]%>%unique

eQTL_file=eQTL_file[ID%chin%feature_metadata$ID & celltype==ln & state==cn,]
eQTL_file$Symbol=feature_metadata[,setNames(Symbol,ID)][eQTL_file$ID]
eQTL_file[,idx:=sprintf("%s_%s_%s",Symbol,celltype,state)]

genotypes$idx=unique(eQTL_file[,.(rsID,idx)])[,setNames(idx,rsID)][genotypes$ID]
genotypes[,Symbol:=str_split(idx,"_",simplify=T)[,1]]
genotypes[,celltype:=str_split(idx,"_",simplify=T)[,2]]
genotypes[,state:=str_split(idx,"_",simplify=T)[,3]]

popq=str_split(comparison,"_",simplify=T)[,1]
popr=str_split(comparison,"_",simplify=T)[,2]

Expression_info=Expression[celltype==ln & (POP==popq | POP==popr),]

if(exp_resp=='DR'){
  setnames(Expression_info,'logFC','logCPM',skip_absent=TRUE) # change name to harmonize downstream code
}
genotypes_info=genotypes[celltype==ln]

geno=genotypes_info[state==cn,]

################################################################################
# Run mediation analysis

# for each gene and condition, compute the fraction of variance in
# expression/response explained by the gene's best eQTL in each lineage
result=lapply(geno[,unique(Symbol)],function(g){
  data=Expression_info
  data$n_alt=geno[Symbol==g,setNames(n_alt,IID)][as.character(data$IID)]
  if (data[,length(table(n_alt))]>1) {
    if (!grepl("ASH",comparison)) {
      mod_1=data[state==cn & Symbol==g,lm(logCPM~n_alt+POP+Age+Mortality)]
      mod_2=data[state==cn & Symbol==g,lm(n_alt~POP+Age+Mortality)]
    } else if (grepl("ASH",comparison)) {
      mod_1=data[state==cn & Symbol==g,lm(logCPM~n_alt+POP+Age+GenderF+Mortality)]
      mod_2=data[state==cn & Symbol==g,lm(n_alt~POP+Age+GenderF+Mortality)]
    }
    rsID=geno[Symbol==g,unique(ID)]
    med=mediate(mod_2,mod_1,treat="POP",mediator="n_alt",control.value=popr,treat.value=popq)
    res=data.table(Symbol=g,ID=feature_metadata[Symbol==g,ID],frac_var=med$n0,pval=med$n0.p,rsID)
  } else {
    rsID=geno[Symbol==g,unique(ID)]
    res=data.table(Symbol=g,ID=feature_metadata[Symbol==g,ID],frac_var=0,pval=1,rsID)
  }
  print(sprintf("Done: %s. %s, %s",g,comparison,cn))
  return(res)
})%>%rbindlist()
result[,state:=cn]
result[,lineage:=ln]
result[,comp:=comparison]
toc()

fwrite(result,sprintf("%s/genetics/%s_%s__%s%s%s_%s.tsv",OUT_DIR,ln,cn,comparison,ifelse(popdiff==T,"_popD","_"),gsub("^D","",exp_resp),run_id),sep="\t")
