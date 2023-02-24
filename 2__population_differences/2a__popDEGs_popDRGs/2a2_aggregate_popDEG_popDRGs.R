################################################################################
################################################################################
# File name: 2a2__aggregate_popDEGs_popDRGs.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: combine popDE and popDR results (adjusted and raw into a single table)
# Effector script
################################################################################
################################################################################
### objective of the script
# Assess significance of effect size decrease after adjustments

# which population differences are actually lost (i.e., associated to signif.
# decrease in effect size) after adjusting for cellular composition?

# no adjustment > within-lineage (lineage) adjustment
# x := larger group (i.e., no)
# y := smaller group (i.e., lineage)

# estimate correlation between raw (beta_x) and adjusted (beta_y) under the null (permuted data)
# even if the data are permuted, the beta_x and beta_y are estimated from the same data and could hence be correlated
# this dependence needs to be accounted for when performing the t-test, lest we underestimate the number of significant changes


ANALYSE="popdiff___lineage_condition__noCellProps_noSVs___220409_perm0.tsv"
popDE_noCellProp=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDE_noCellProp[,type:="expression"]
ANALYSE="popdiff__lineage_condition_logFC__noCellProps_noSVs___220409_perm0.tsv"
popDR_noCellProp=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDR_noCellProp[,type:="response"]

# load adjusted population effect estimates on gene expression and response
ANALYSE="popdiff___lineage_condition__CellPropLineage_noSVs___220409_perm0.tsv"
popDE_lineageCellProp=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDE_lineageCellProp[,type:="expression"]
ANALYSE="popdiff___lineage_condition_logFC__CellPropLineage_noSVs___220409_perm0.tsv"
popDR_lineageCellProp=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDR_lineageCellProp[,type:="response"]


# load raw population effect estimates on gene expression and response
ANALYSE="popdiff___lineage_condition__noCellProps_noSVs___220409_perm1.tsv"
popDE_noCellProp_null=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDE_noCellProp_null[,type:="expression"]
ANALYSE="popdiff__lineage_condition_logFC__noCellProps_noSVs___220409_perm1.tsv"
popDR_noCellProp_null=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDR_noCellProp_null[,type:="response"]

# load adjusted population effect estimates on gene expression and response
ANALYSE="popdiff___lineage_condition__CellPropLineage_noSVs___220409_perm1.tsv"
popDE_lineageCellProp_null=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDE_lineageCellProp_null[,type:="expression"]
ANALYSE="popdiff___lineage_condition_logFC__CellPropLineage_noSVs___220409_perm1.tsv"
popDR_lineageCellProp_null=fread(sprintf("2__population_differences/data/%s",ANALYSE))
popDR_lineageCellProp_null[,type:="response"]

# annotate, merge and reshape
popDE_noCellProp_null[,cellprop_adjust:="no"]
popDE_lineageCellProp_null[,cellprop_adjust:="lineage"]
popDE_adjusted_null=rbind(popDE_noCellProp,popDE_lineageCellProp)
popDE_adjusted_null[,idx:=sprintf("%s_%s_%s",ID,celltype,state)]
popDE_adjusted_null_comp[,group:=ifelse(cellprop_adjust=="no","x","y")]

popDR_noCellProp_null[,cellprop_adjust:="no"]
popDR_lineageCellProp_null[,cellprop_adjust:="lineage"]
popDR_adjusted_null=rbind(popDR_noCellProp,popDR_lineageCellProp)
popDR_adjusted_null[,idx:=sprintf("%s_%s_%s",ID,celltype,state)]
popDR_adjusted_null_comp[,group:=ifelse(cellprop_adjust=="no","x","y")]

popDE_adjusted_null_wide=dcast(popDE_adjusted_null_comp,idx~group,value.var=c("beta","FDR","se"))

# correlation between raw and adjusted population effects on expression under null
corr_beta_no_lineage_de=popDE_adjusted_null_wide[,cor.test(beta_x,beta_y)]$estimate
# x := no, y := lineage
# > popDE_adjusted_null_wide[,cor.test(beta_x,beta_y)]
#
#         Pearson's product-moment correlation
#
# data:  beta_x and beta_y
# t = 818.66, df = 190063, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8816486 0.8836351
# sample estimates:
#       cor
# 0.8826458

# test for significant differences in raw and adjusted effect sizes, given
# correlation under the null
significant_popde_losses=lapply(c("DE","DR"),function(exp_resp){
  noCellProp=get(sprintf("pop%s_noCellProp",exp_resp))
  lineageCellProp=get(sprintf("pop%s_lineageCellProp",exp_resp))
  noCellProp[,cellprop_adjust:="no"]
  lineageCellProp[,cellprop_adjust:="lineage"]
  adjusted=rbind(noCellProp,lineageCellProp)
  adjusted[,idx:=sprintf("%s_%s_%s",ID,celltype,state)]

  adjusted_comp[,group:=ifelse(cellprop_adjust=="no","x","y")]

  signif_in_xy=dcast(
    adjusted_comp,idx~group,
    value.var=c("beta","FDR","se")
  )%>%.[(FDR_x<0.01&abs(beta_x)>0.2)&(FDR_y<0.01&abs(beta_y)>0.2),idx]
  signif_in_x=dcast(
    adjusted_comp,idx~group,
    value.var=c("beta","FDR","se")
  )%>%.[(FDR_x<0.01&abs(beta_x)>0.2)&!(FDR_y<0.01&abs(beta_y)>0.2),idx]
  signif_in_y=dcast(
    adjusted_comp,idx~group,
    value.var=c("beta","FDR","se")
  )%>%.[!(FDR_x<0.01&abs(beta_x)>0.2)&!(FDR_y<0.01&abs(beta_y)>0.2),idx]

  adjusted_comp$signif_in=case_when(
    adjusted_comp$idx%in%signif_in_xy~"xy",
    adjusted_comp$idx%in%signif_in_x~"x",
    adjusted_comp$idx%in%signif_in_y~"y"
  )

  adjusted_wide=dcast(adjusted_comp[!is.na(signif_in)],idx+signif_in~group,value.var=c("beta","FDR","se"))

  adjusted_wide=adjusted_wide[!grepl("^IAV|^SARS",idx),]


  betas_t=adjusted_wide[,.(t=(beta_x-beta_y)/(sqrt(se_x^2+se_y^2-2*corr_beta_no_lineage_de*se_x*se_y))),by=.(idx)]

  adjusted_wide$t=betas_t[,setNames(t,idx)][adjusted_wide$idx]

  adjusted_wide$different=ifelse(abs(adjusted_wide$t)>1.96,"Different","Same")

  adjusted_wide$Symbol=adjusted[,setNames(Symbol,idx)][adjusted_wide$idx]
  adjusted_wide$state=str_split(adjusted_wide$idx,"_",simplify=T)[,3]
  adjusted_wide$celltype=str_split(adjusted_wide$idx,"_",simplify=T)[,2]

  adjusted_labels=adjusted_wide[different=="Different"&signif_in!="all",][order(-abs(t)),][,head(.SD,3),by=.(signif_in,state,celltype)][,.(beta_x,beta_y,signif_in,celltype,state,Symbol)]

  adjusted_wide[,difference_type:=case_when(signif_in=="x"~"difference_lost",signif_in=="y"~"difference_gained",signif_in=="xy"~"difference_kept")]
  setnames(adjusted_wide,"different","difference_signif")
  result=adjusted_wide[,-"signif_in"]
	result[,type:=exp_resp]
  return(result)
})%>%rbindlist()

fwrite(significant_popde_losses,"2__population_differences/data/cell_composition_adjustment_beta_changes.tsv",sep="\t")
