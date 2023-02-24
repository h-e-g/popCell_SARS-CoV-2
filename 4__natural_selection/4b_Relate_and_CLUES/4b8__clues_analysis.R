################################################################################
################################################################################
# File name: 4b8__clues_analysis.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Analyse CLUES results
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("./00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("./misc_plots.R",MISC_DIR))
source(sprintf("./misc_functions.R",MISC_DIR))
source(sprintf("./load_eQTLs.R",MISC_DIR))
source(sprintf("./querySNPs.R",MISC_DIR))

################################################################################
# Calculate selection metric (Z-score)

# load (r)eQTL allele frequency trajectories calculated by CLUES per population
# and summarized in 4b7__summarize_clues.py
# smooth allele frequency trajectories using loess regression
# calculate first derivative (dotallele), adjust for frequency differences
# (dotallele_norm) and compute Z-score (z)

clues <- lapply(c("YRI","CEU","CHS"),function(p){
  data <- lapply(c("eQTL","reQTL"),function(expresp){
  data=fread(sprintf("%s/%s_CutOff2000/SNPs/%s/allele_freq_max__all_SNPs.tsv.gz",DAT_CLUES_DIR,p,expresp))
  colnames(data)[1:5]=c("epoch","allele_pp_max","rsID","lowerCI","upperCI")
  data[,allele_pp_max_smooth:=predict(loess(allele_pp_max~epoch,span=0.1),epoch),by=rsID]
  # data[order(rsID,-epoch),dotallele:=c(diff(allele_pp_max),0),by=.(rsID)]
  # data[,dotallele_norm:=dotallele/sqrt(allele_pp_max*(1-allele_pp_max))]
  # data[,z:=dotallele_norm/sd(dotallele_norm),by=epoch]
  data[order(rsID,-epoch),dotallele_smooth:=c(diff(allele_pp_max_smooth),0),by=.(rsID)]
  data[,dotallele_norm_smooth:=dotallele_smooth/sqrt(allele_pp_max_smooth*(1-allele_pp_max_smooth))]
  data[,z_smooth:=dotallele_norm_smooth/sd(dotallele_norm_smooth),by=epoch]
  data[,type:=expresp]
  data[,pop:=p]
  data$gene_name=get(sprintf("%s_Signif_both",expresp))[,setNames(gene_name,snps)][data$rsID]
  return(data)
  })%>%rbindlist()
})%>%rbindlist()

################################################################################
# Add natural selection metrics

SNP_info=getMap(annotate=TRUE)
clues[,P_PBS_CAT:=case_when(P_PBS<0.05~"0.05",P_PBS<0.01~"0.01",T~"weak")]
clues$gene_name=eQTL_Signif_both[,match(cluesrsID,snps),gene_name]

# write results
fwrite(clues,sprintf("%s/clues_trajectories.tsv.gz",DAT_CLUES_DIR),sep="\t")

################################################################################
# date of selection
selection_date=clues[!is.na(z_smooth),.(start=max(c(-1,epoch[abs(z_smooth)>3])),
                                            end=min(c(2001,epoch[abs(z_smooth)>3])),
                                            max_abs_z=max(abs(z_smooth)),
                                            Selection_effect=sum(z_smooth[abs(z_smooth)>3]),
                                            ,by=.(rsID,gene_name,type,pop,PBS,P_PBS)]
selection_date[,selected:=ifelse(max_abs_z>3,T,F)]
selection_date[,window:=case_when(start<=970&start>=770~"start_window",start>970&end<970~"overlap_window",T~"not_window")]
# selection_date[,length(unique(rsID)),by=.(pop,type,reQTL_COV_specific=(rsID%in%snpSets[set==myset,unique(snps)]),window)]

selection_date$z_max <- clues[!is.na(z_smooth),max(abs(z_smooth)),by=.(rsID,type,pop)][,setNames(V1,sprintf("%s%s%s",rsID,type,pop))][selection_date[,sprintf("%s%s%s",rsID,type,pop)]]

# write results
fwrite(selection_date,sprintf("%s/selection_date.tsv.gz",DAT_CLUES_DIR),sep="\t")

###############################################################################

COV_reQTL=snpSets[set=="reQTL_COV",unique(snps)]
sum(COV_reQTL%chin% cluesFreq[,unique(rsID)])

eQTL_Stats_celltype_mash[,beta_lower:=sign(beta)*(pmax(abs(beta)-2*se,0))]
reQTL_Stats_celltype_mash[,beta_lower:=sign(beta)*(pmax(abs(beta)-2*se,0))]
eQTL_best_celltype=eQTL_Stats_celltype_mash[order(snps,gene_name,-abs(beta_lower)),head(.SD[,.(best_celltype=celltype,best_state=state,beta_lower,pvalue_best=pvalue,type='eQTL')],1),by=.(snps,gene_name)]
reQTL_best_celltype=reQTL_Stats_celltype_mash[order(snps,gene_name,-abs(beta_lower)),head(.SD[,.(best_celltype=celltype,best_state=state,beta_lower,pvalue_best=pvalue,type='reQTL')],1),by=.(snps,gene_name)]

eQTL_Stats_celltype_mash[order(snps,gene_name,-abs(beta_lower)),head(.SD[,.(snps,gene_name,celltype,state,beta_lower,pvalue)],1),by=.(snps,gene_name)]

all_selected=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig5/test_PBS/geneList_all_selected_v2.tsv",EIP))
all_selected=merge(all_selected,selection_date[,.(snps=rsID,gene_name,POP=pop,start,end,max_abs_z,Selection_effect,Selection_effect_weak,start_weak,end_weak,selected)],by=c('snps','gene_name','POP'),all.x=TRUE)
all_selected=merge(all_selected,rbind(eQTL_best_celltype,reQTL_best_celltype),by=c('snps','gene_name','type'),all.x=TRUE)
selected_alleles=all_selected[METRIC=='PBS',.(snps,gene_name,type,POP,celltype,state,pvalue,best_celltype,best_state,beta_lower,
                                P_PBS=P,max_abs_z,Selection_effect_weak,start_weak,end_weak,
                                selected,start,end,Selection_effect,
                                GO_immune,ISG_effector,IEI,Targeted_viruses,
                                REF,ALT,ALT_DERANC,DAF_or_MAF_CEU, DAF_or_MAF_CHS, DAF_or_MAF_YRI,minP_perGene)]

selected_allele=selected_alleles[order(minP_perGene,-abs(Selection_effect_weak),pvalue),][!duplicated(paste(POP,snps,gene_name,type)),]

fwrite(selected_alleles,file=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V9/suppTables/TableS7/TableS7B_1_selection_dates_all_genes.tsv",EIP))


#### comparison of the frequency of adaptive events between African and Europeans at popDEGs with mediation.
mediation=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V8/suppTables/TableS6/TableS6_mediation_analyses.tsv",EIP))
setnames(mediation,'mediator','rsID')
mediated_annot=merge(mediation[mediator_type=="genetics",],selection_date[pop!='CHS',],by='rsID')[frac_var>.5,][!duplicated(paste(rsID,pop)),]
mediated_annot[,mean(max_abs_z>3),by=pop]
# pop        V1
# 1: YRI 0.2091097
# 2: CEU 0.3405640
mediated_annot[,fisher.test(table(max_abs_z>3,pop))]
# Fisher's Exact Test for Count Data
# p-value = 7.765e-06
# odds ratio
#  0.5123159
# 95 percent confidence interval:
#  0.3779538 0.6923835

selection_date[rsID=='rs1142888',]
# rsID gene_name type pop start  end max_abs_z Selection_effect Selection_effect_weak start_weak end_weak selected
# 1: rs1142888      GBP7 eQTL YRI    -1 2001  1.320060           0.0000             0.5381256         -1     2001    FALSE
# 2: rs1142888      GBP7 eQTL CEU  1272  782  4.327977         638.5713             1.2034810       1315      753     TRUE
# 3: rs1142888      GBP7 eQTL CHS    -1 2001  2.904957           0.0000             1.1528575       1742      307    FALSE
