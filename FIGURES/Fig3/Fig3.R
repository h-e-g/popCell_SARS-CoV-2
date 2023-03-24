################################################################################
################################################################################
# File name: Fig3.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: output panels for Figure 3
################################################################################
################################################################################

################################################################################
# Setup

# load required packages
LIB_DIR="LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="MISC"
source(sprintf("%s/shortcuts.R",MISC_DIR))

# declare useful functions
source(sprintf("%s/misc_functions.R",MISC_DIR))

# declare useful functions and variables for plotting
source(sprintf("%s/misc_plots.R",MISC_DIR))
source(sprintf("%s/load_eQTLs.R",MISC_DIR))
source(sprintf("%s/querySNPs.R",MISC_DIR))

# declare outputDIR
FIGURE_DIR=sprintf('%s/Fig3/',FIG_DIR)

################################################################################
# Fig. 3a : Number of eQTL_per Gene in each lineage

eQTL_perGene_lineage=eQTL_Signif_both[celltype%chin%lineage_order,.(N_eQTL_perGene=length(unique(snps))),by=.(celltype,gene)][,.N,keyby=.(celltype,N_eQTL_perGene)]
eQTL_perGene_lineage[,N_eQTL_perGene_maxed:=factor(ifelse(N_eQTL_perGene<4,as.character(N_eQTL_perGene),'4+'),rev(c('1','2','3','4+')))]
fwrite(eQTL_perGene_lineage,file=sprintf("%s/Fig3Adata_eQTL_perGene_lineage.tsv",SOURCE_DATA_DIR),sep='\t')

Fig3A_data=eQTL_perGene_lineage[,.(Ngene=sum(N)),keyby=.(celltype,N_eQTL_perGene_maxed)]

Fig3A <- ggplot(Fig3A_data,aes(x=N_eQTL_perGene_maxed, y=Ngene, fill=celltype))+geom_bar(stat='Identity')+facet_grid(row=vars(celltype))+ theme_plot()
Fig3A <- Fig3A + scale_fill_manual(values=lineage_color) + xlab('Number of eQTL per gene') + ylab('Number of genes') + guides(fill= "none") + coord_flip()
Fig3A

pdf(sprintf('%s/Fig3A_eQTL_perGene_lineage.pdf',FIGURE_DIR),width=7.2*.22,height=6.7*.4)
print(Fig3A)
dev.off()


################################################################################
# Fig. 3b : Number of eQTL per celltype

Fig3B_data=eQTL_Stats_celltype[celltype%chin%celltype_order,.(N_eQTL_GWsignif=length(unique(snps[!is.na(top_component)])),
                                                                          N_eQTL_GWsignif_CTspecific_praw1=length(unique(snps[!is.na(top_component) & Number_celltype_where_signif_praw1==1])),
                                                                          N_eQTL_signif_praw1=length(unique(snps[pvalue<0.01])))
                                                                            ,by=.(celltype,lineage)]

fwrite(Fig3B_data,file=sprintf('%s/Fig3Bdata_Nb_eQTL_by_celltype.tsv',SOURCE_DATA_DIR),sep='\t')
####------------------------------------------------------------------------------------####

Fig3B_data=Fig3B_data[order(lineage,-N_eQTL_GWsignif)]

Fig3B_data[,N_eQTL_GWsignif_shared_praw1:=N_eQTL_GWsignif-N_eQTL_GWsignif_CTspecific_praw1]
Fig3B_data[,N_eQTL_signif_praw1_notGW:=N_eQTL_signif_praw1-N_eQTL_GWsignif]
Fig3B_data[,sqrtN_eQTL_GWsignif_CTspecific_praw1:=sqrt(N_eQTL_GWsignif_CTspecific_praw1)]
Fig3B_data[,sqrtN_eQTL_GWsignif_shared_praw1:=sqrt(N_eQTL_GWsignif_shared_praw1)-sqrt(N_eQTL_GWsignif_CTspecific_praw1)]
Fig3B_data[,sqrtN_eQTL_signif_praw1_notGW:=sqrt(N_eQTL_signif_praw1)-sqrt(N_eQTL_GWsignif)]
Fig3B_data=Fig3B_data[,.(lineage,celltype,sqrtN_eQTL_GWsignif_CTspecific_praw1,sqrtN_eQTL_GWsignif_shared_praw1,sqrtN_eQTL_signif_praw1_notGW)]

Fig3B_data=melt(Fig3B_data,id.vars=c('celltype','lineage'))
Fig3B_data[,celltype:=factor(celltype,Nb_of_eQTL_per_celltype$celltype)]
Fig3B_data[,variable_clear:=case_when(variable=='sqrtN_eQTL_GWsignif_CTspecific_praw1'~"FDR<1%, cell-type-specific",
                                                        variable=='sqrtN_eQTL_GWsignif_shared_praw1'~"FDR<1%, shared",
                                                        variable=='sqrtN_eQTL_signif_praw1_notGW'~"p<0.01, shared",
                                                        TRUE~'NA')]
Fig3B_data[,variable_clear:=factor(variable_clear,c("p<0.01, shared","FDR<1%, shared","FDR<1%, cell-type-specific"))]


ticks_pos=c(10,100,250,500,1000,2500,5000)
p <- ggplot(Fig3B_data,aes(x=celltype, fill = celltype, y= value, pattern = variable_clear,alpha=variable_clear)) + scale_y_continuous(breaks=sqrt(ticks_pos),labels=ticks_pos)
p <- p + geom_bar_pattern(position='stack',stat='Identity', color='black', aes(pattern_angle=variable_clear,pattern_color=variable_clear), pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6,width=.7) #color = "black", pattern_fill = "black", pattern_angle = 45,
p <- p + scale_fill_manual(values=celltype_color) + scale_pattern_manual(values = c("FDR<1%, cell-type-specific" = "stripe", "FDR<1%, shared" = "none","p<0.01, shared"="stripe"))
p <- p + ylab('Number of eQTL') + guides(fill= "none") + theme_plot() + theme(legend.title = element_blank()) + theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
p <- p + scale_pattern_angle_manual(values = c("FDR<1%, cell-type-specific" = -45, "FDR<1%, shared" = 45, "p<0.01, shared"=45))
p <- p + scale_pattern_colour_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" = NA,"p<0.01, shared"="black"))
p <- p + scale_pattern_fill_manual(values = c("FDR<1%, cell-type-specific" = "white", "FDR<1%, shared" =NA,"p<0.01, shared"="black"))
p <- p + theme(axis.title.x=element_blank(),legend.position="top",text = element_text(size = 8))
Fig3B <- p +theme(legend.key.height = unit(.1,"cm"), legend.key.width = unit(.2,"cm"))

####------------------------------------------------------------------------------------####
pdf(sprintf('%s/Fig3B_Nb_eQTL_by_celltype.pdf',FIGURE_DIR),width=7.2*.4,height=6.7*.4)
print(Fig3B)
dev.off()
####------------------------------------------------------------------------------------####

################################################################################
# Fig. 3c : correlation of eQTL effect sizes

get_pairwise_sharing_raw=function (beta_mat, FDR_mat, factor = 0.5, FDR_thresh = 0.01,cor=FALSE,method='s'){
    R = ncol(beta_mat)
    S = matrix(NA, nrow = R, ncol = R)
    A = matrix(NA, nrow = R, ncol = R)
    for (i in 1:R) {
        for (j in i:R) {
            if(!all(is.na(beta_mat[, i]) | is.na(beta_mat[, j])) ){
              a = FDR_mat[,i]<FDR_thresh | FDR_mat[,j]<FDR_thresh
              if(cor){
              S[i, j] = cor(beta_mat[a, i],beta_mat[a, j],method=method)
              A[i, j]=sum(a)
               if(cor.test(beta_mat[a, i],beta_mat[a, j],method=method)$p.value>0.01){
                  S[i, j] = 0
                }
              }else{
              ratio = beta_mat[a, i]/beta_mat[a, j]
              S[i, j] = mean(ratio > factor & ratio < (1/factor))
              }
            }
        }
    }
    S[lower.tri(S, diag = FALSE)] = t(S)[lower.tri(S, diag = FALSE)]
    colnames(S) = row.names(S) = colnames(beta_mat)
    return(S)
}

celltype_2lineage_ordered=data.table(celltype=c("B.M.K","B.M.L","B.N.K","B.N.L","Plasmablast","MONO.CD14","MONO.CD14.INFECTED", "MONO.CD16", "cDC","pDC","NK.CD56brt","NK.CD56dim","NK.M.LIKE","T.CD8.EMRA","MAIT","T.CD8.CM.EM","T.CD8.N","T.CD4.E","T.CD4.N","T.Reg","T.gd","ILC"),
lineage=c("B","B","B","B","B","MONO","MONO", "MONO", "MONO","MONO","NK","NK","NK","T.CD8","T.CD8","T.CD8","T.CD8","T.CD4","T.CD4","T.CD4","T.CD8","T.CD8"))
celltype_2lineage_noINF=celltype_2lineage_ordered[celltype!='MONO.CD14.INFECTED',]

######### assess sharing of eQTL in NS condition
beta_mat=dcast(eQTL_Stats_celltype,snps+gene+gene_name+state~celltype,value.var='beta')
beta_mat_NS=as.matrix(beta_mat[state=='NS',-c('snps','gene','gene_name','state')])

pval_mat=dcast(eQTL_Stats_celltype,snps+gene+gene_name+state~celltype,value.var='pvalue')
pval_mat_NS=as.matrix(pval_mat[state=='NS',-c('snps','gene','gene_name','state')])

sharing_cor_NS=get_pairwise_sharing_raw(beta_mat_NS,pval_mat_NS,cor=TRUE,FDR_thresh=0.01)

######### assess sharing of reQTL in COV condition
beta_mat_reQTL=dcast(reQTL_Stats_celltype,snps+gene+gene_name+state~celltype,value.var='beta')
beta_mat_reQTL_COV=as.matrix(beta_mat_reQTL[state=='COV',-c('snps','gene','gene_name','state')])

pval_mat_reQTL=dcast(reQTL_Stats_celltype,snps+gene+gene_name+state~celltype,value.var='pvalue')
pval_mat_reQTL_COV=as.matrix(pval_mat_reQTL[state=='COV',-c('snps','gene','gene_name','state')])

sharing_cor_reQTL_COV=get_pairwise_sharing_raw(beta_mat_reQTL_COV,pval_mat_reQTL_COV,cor=TRUE,FDR_thresh=0.01)

######### create upper/lower diagonal matrix with NS eQTL /COV reQTL
sharing_cor_NSCOV=sharing_cor_NS[celltype_2lineage_noINF$celltype,celltype_2lineage_noINF$celltype]
for (i in 1:(nrow(celltype_2lineage_noINF)-1)){
  for (j in (i+1):nrow(celltype_2lineage_noINF)){
    sharing_cor_NSCOV[i,j]=sharing_cor_reQTL_COV[celltype_2lineage_noINF$celltype[i],celltype_2lineage_noINF$celltype[j]]
  }
}
####------------------------------------------------------------------------------------####
Fig3C_data=data.table(sharing_cor_NSCOV,celltype_2lineage_noINF,celltype_color=celltype_color[celltype_2lineage_noINF$celltype],lineage_color=lineage_color[celltype_2lineage_noINF$lineage])
fwrite(Fig3C_data,file=sprintf('%s/Fig3Cdata_eQTLsharing_raw_NSCOV.tsv',SOURCE_DATA_DIR),sep='\t')

####------------------------------------------------------------------------------------####
celltype_2lineage_noINF=Fig3C_data[,.(celltype,lineage)]
sharing=as.matrix(Fig3C_data[,mget(celltype_2lineage_noINF$celltype)])
rownames(sharing)=celltype_2lineage_noINF$celltype
ann_colors=list(celltype=Fig3C_data[,setNames(celltype_color,celltype)],lineage=Fig3C_data[!duplicated(lineage),setNames(lineage_color,lineage)])
pdf(sprintf('%s/Fig3C_eQTLsharing_raw_NSCOV.pdf',FIGURE_DIR),height=4,width=4)
Fig3C <- ComplexHeatmap::pheatmap(sharing,cluster_cols=FALSE,cluster_rows=FALSE,fontsize=8,
  annotation_col=celltype_2lineage_noINF,annotation_row=celltype_2lineage_noINF,annotation_colors=ann_colors,annotation_legend = FALSE,legend=TRUE)
dev.off()
####------------------------------------------------------------------------------------####

################################################################################
# Fig. 3c-related: testing for difference in effect size correlation (between/within lineage)

COR_comparion=as.data.table(reshape2::melt(sharing_cor_NSCOV,varnames=c('celltype_1','celltype_2'),value.name='correlation'))
COR_comparion=merge(COR_comparion,cbind(celltype_2lineage,N=1:nrow(celltype_2lineage)),by.x='celltype_1',by.y='celltype')
COR_comparion=merge(COR_comparion,cbind(celltype_2lineage,N=1:nrow(celltype_2lineage)),by.x='celltype_2',by.y='celltype',suffix=c('_1','_2'))
COR_comparion=COR_comparion[celltype_2!=celltype_1,]
COR_comparion[,COND:=ifelse(N_1<N_2,'COV','NS')]
COR_comparion[,mean(correlation,na.rm=T),by=.(lineage_1==lineage_2,COND)]
# compare correlation within/between lineage
COR_comparion[,wilcox.test(correlation[COND=='NS' & lineage_1==lineage_2],correlation[COND=='NS' & lineage_1!=lineage_2])$p.value]
COR_comparion[,wilcox.test(correlation[COND=='COV' & lineage_1==lineage_2],correlation[COND=='COV' & lineage_1!=lineage_2])$p.value]
# compare correlation eQTL VS reQTL
COR_comparion[,wilcox.test(correlation[COND=='COV'],correlation[COND=='NS'])$p.value]


################################################################################
# Fig. 3d : comparison of reQTL effect sizes

Fig3D_data=Signif_reQTL_compare[celltype%in%lineage_5,][sample(1:.N),]
Fig3D_data=Fig3D_data[!duplicated(paste(snps,gene,lineage)),.(snps,gene_name,lineage,beta_COV,pvalue_COV,beta_IAV,pvalue_IAV,p_diff_COV_IAV)]
fwrite(Fig3D_data,file=sprintf('%s/Fig3Ddata_reQTL_Signif_lineage_compared.tsv',SOURCE_DATA_DIR),sep='\t')

####------------------------------------------------------------------------------------####
p <- ggplot(Fig3D_data,aes(x=beta_COV,y=beta_IAV,col=lineage))+xlab(expression(beta[COV]))+ylab(expression(beta[IAV]))
p <- p + rasterize(geom_point(alpha=.5),dpi=400)+xlim(c(-1.5,1.5))+ylim(c(-1.5,1.5))#+theme_plot()
p <- p + scale_color_manual(values=lineage_color)+ theme_plot()
p <- p + geom_hline(yintercept=0,col='lightgrey',linetype=2)+ geom_vline(xintercept=0,col='lightgrey',linetype=2)
Fig3D <- p + theme(text = element_text(size = textSize))+guides(color=guide_legend(nrow=2,overide.aes=list(alpha=1), title.theme =element_blank()))
pdf(sprintf('%s/Fig3D_reQTL_Signif_lineage_compared.pdf',FIGURE_DIR),width=7.2*.28,height=6.7*.4)
print(Fig3D)
dev.off()

################################################################################
# Fig. 3e : number of  eQTL with different effect sizes

Fig3E_data=Fig3D_data[,.(number_different_praw1=sum(p_diff_COV_IAV<0.01)),keyby=.(lineage,stronger_in=ifelse(abs(beta_COV)>abs(beta_IAV),'COV','IAV'))]
fwrite(Fig3E_data,file=sprintf('%s/Fig3Edata_number_virus_specific_reQTLs_bylineage.tsv',SOURCE_DATA_DIR),sep='\t')

####------------------------------------------------------------------------------------####
 p <- ggplot(Fig3E_data,aes(x=lineage,y=number_different_praw1,alpha=gsub('COV','Sars-CoV-2',stronger_in),fill=lineage))
 p <- p + geom_bar(position='stack',stat='Identity')+ theme_plot() + scale_alpha_manual(values=c('Sars-CoV-2'=1,IAV=.5))
 p <- p + scale_fill_manual(values=lineage_color) + ylab('Number of virus-dependent reQTL')
 p <- p + guides(fill= "none") + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1))
 pdf(sprintf('%s/Fig3E_number_virus_specific_reQTLs_bylineage.pdf',FIGURE_DIR),height=5,width=3)
 print(p)
 dev.off()


 ################################################################################
 # Fig. 3f : MMP1 reQTL

Fig3F_data=get_eQTL('rs534191','MMP1',resolution='lineage',metric='logCPM')
Fig3F_data=Fig3F_data[Symbol=='MMP1',]
Fig3F_data[,state:=factor(state,c('NS','COV','IAV'))]
nIQR=3
Fig3F_data[,logCPM_IQR:=pmax(pmin(logCPM,median(logCPM)+nIQR*IQR(logCPM)),median(logCPM)-nIQR*IQR(logCPM)),by=.(state,celltype)]
Fig3F_data[,Number_of_ALT_alelle:=as.factor(round(Number_of_ALT_alelle,0))]
Fig3F_data=Fig3F_data[,.(donorID=IID,POP,Symbol,celltype,state,logCPM,logCPM_IQR,Number_of_ALT_alelle)]
fwrite(Fig3F_data,file=sprintf('%s/Fig3F_data_eQTL_MMP1_MONO.tsv',SOURCE_DATA_DIR),sep='\t')

####------------------------------------------------------------------------------------####
p <- ggplot(eQTL_MMP1[Symbol=='MMP1' & celltype%in%c('MONO'),],aes(x=Number_of_ALT_alelle,y=logCPM_IQR,fill=state))
p <- p + geom_violin(scale='width') + geom_boxplot(fill='white',alpha=.5,notch=TRUE)
p <- p + theme_plot() + scale_fill_manual(values=color_conditions)# + scale_color_manual(values=color_populations)
p <- p + facet_grid(state~celltype) +ylab('logCPM')
Fig3F=p+theme(legend.position='none',text = element_text(size = textSize))
pdf(sprintf('%s/Fig3F_eQTL_MMP1_MONO.pdf',FIGURE_DIR),width=7.2*.2,height=6.7*.4)
print(Fig3F)
dev.off()

################################################################################
# Fig. 3G

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
# load data on beta_change after adjustments
popdiff_losses=fread(sprintf("2_population_differences/data/cell_composition_adjustment_beta_changes.tsv"))

# fraction of all genes with an eQTL
eqtls<-eQTL_Signif_both[celltype%in%lineage_order,]
setnames(eqtls,"gene_name","Symbol")
eqtls_genome_wide=eqtls[,length(unique(Symbol))/12672*100,by=celltype]
setnames(eqtls_genome_wide,"V1","frac_eqtl")
eqtls_genome_wide$adj="gw"
eqtls_genome_wide$type="DE"

# fraction of all genes with an reQTL
reqtls<-reQTL_Signif_both[celltype%in%lineage_order,]
setnames(reqtls,"gene_name","Symbol")
reqtls_genome_wide=reqtls[,length(unique(Symbol))/12672*100,by=celltype]
setnames(reqtls_genome_wide,"V1","frac_eqtl")
reqtls_genome_wide$adj="gw"
reqtls_genome_wide$type="DR"

# fraction of popDEGs (raw) with an eQTL
fig3g_data=lapply(c("DE","DR"),function(z){
  lapply(c("no","lineage"),function(y){
    popdegs=get(sprintf("pop%s_%sCellProp",z,y))
    lapply(lineage_order,function(x){
      if (y=="lineage"){
      pdegs_in_lineage=popdegs[FDR<0.01&abs(beta)>0.2&celltype==x,unique(Symbol)]
      signif_losses=popdiff_losses[type==z&celltype==x&adjustment_phase=="no_lineage"&difference_signif=="Different"&difference_type=="difference_lost",unique(Symbol)]
      pdegs_in_lineage=pdegs_in_lineage[pdegs_in_lineage%nin%signif_losses]
      } else {
      pdegs_in_lineage=popdegs[FDR<0.01&abs(beta)>0.2&celltype==x,unique(Symbol)]
      }
      pdegs_in_lineage_w_eqtl=eqtls[celltype==x&Symbol%in%pdegs_in_lineage,length(unique(Symbol))]
      n=pdegs_in_lineage_w_eqtl/length(pdegs_in_lineage)*100
      data.table(celltype=x,frac_eqtl=n,adj=y,type=z)
    })%>%rbindlist()
  })%>%rbindlist()
})%>%rbindlist()
fig3g_data=rbind(fig3g_data,eqtls_genome_wide,reqtls_genome_wide)
fig3g_data$celltype=factor(fig3g_data$celltype,rev(lineage_order))
fig3g_data$adj=factor(fig3g_data$adj,rev(c("gw","no","lineage")))

fwrite(fig3g_data,sprintf("%s/fig3g_data.tsv",SOURCE_DATA_DIR),sep="\t")

fig3g_plot=ggplot(fig3g_data[type=="DE"],aes(frac_eqtl,celltype,fill=celltype,alpha=adj,group=adj))+
  geom_vline(data=fig3g_data[type=="DE",mean(frac_eqtl),by=adj],mapping=aes(xintercept=V1,alpha=adj),linetype="dashed",size=0.2)+
  geom_col(position="dodge")+
  scale_x_continuous(limits=c(0,60),breaks=c(0,15,30,45,60),position="bottom")+
  scale_fill_manual(values=lineage_color,guide="none")+
  scale_y_discrete(labels=rev(c("MONO","B","CD4+ T","CD8+ T","NK")),position="left")+
  scale_alpha_manual(values=c("no"=0.6,"lineage"=1,"gw"=0.2),breaks=c("gw","no","lineage"),labels=c("Genome-wide","Raw","Intra-lineage"),
		     guide=guide_legend(override.aes=list(shape=21,color="black")))+
  xlab("Genes with an eQTL (%)")+
  theme_plot()+
  theme(text=element_text(size=5),axis.title.y=element_text(color="white"))+
  theme(panel.spacing=unit(0,"pt"),legend.position="bottom")

p=ggplot(fig3g_data,aes(celltype,frac_eqtl,fill=celltype,alpha=adj,group=adj))+
  geom_point()+
  scale_fill_discrete(guide="none")+
  scale_alpha_manual(values=c("no"=0.6,"lineage"=1,"gw"=0.2),breaks=c("gw","no","lineage"),
		     labels=c("Genome-wide        ","Raw             ","Adjusted"),guide=guide_legend(override.aes=list(shape=21,color=NA,fill="black",size=2),nrow=1,title.position="top",title.hjust=1),name="Cell proportion adjustment")+
  theme_plot()+
  theme(legend.position="bottom",text=element_text(size=5),legend.title=element_text())
fig3g_legend=get_legend(p+theme(text=element_text(size=5),legend.spacing.x=unit(-1,'mm')))


################################################################################
# Fig. 3H

mediation=fread(sprintf('%s/mediation_w_betapop_tableS6.tsv',DAT_MEDIATION_DIR),sep='\t')

celltype_mediators<-c("MONO.CD16","B.M.K","T.CD4.E","T.CD8.EMRA","NK.M.LIKE")
color_mediator=lineage_color[lineage_celltype[,setNames(lineage,celltype)][celltype_mediators]]
attr(color_mediator,"names")=celltype_mediators
color_mediator=c(color_mediator,'Genetics'=grey(.3))


fig3h_data_All<-mediation[,.(N=mean(FDR_global<.01), Npraw1=mean(pval<.01),
			     frac_var_mean=mean(frac_var_std),frac_var_meanpraw1=mean(frac_var_std[pval<.01])),by=.(celltype,state,mediator_type,type,mediator=ifelse(grepl('rs|ss|esv',mediator),'Genetics',mediator))]
fig3h_data_All[,state:=factor(state,c('NS','COV','IAV'))]
fig3h_data_All[,celltype:=factor(celltype,lineage_order)]
#fig3h_data_All[,set:=ifelse(type=='expression','popDE','popDR')]
fig3h_data_All[,set:='popDE/DR']


fig3h_data_eGene <- mediation[eGene==TRUE,.(N=mean(FDR_global<.01), Npraw1=mean(pval<.01),frac_var_mean=mean(frac_var_std),frac_var_meanpraw1=mean(frac_var_std[pval<.01])),by=.(celltype,state,mediator_type,type,mediator=ifelse(grepl('rs|ss|esv',mediator),'Genetics',mediator))]
fig3h_data_eGene[,state:=factor(state,c('NS','COV','IAV'))]
fig3h_data_eGene[,celltype:=factor(celltype,lineage_order)]
#fig3h_data_eGene[,set:=ifelse(type=='DE','eQTL-popDE','reQTL-popDR')]
fig3h_data_eGene[,set:='eQTL-popDE/DR']

fig3h_data=rbind(fig3h_data_All,fig3h_data_eGene)
fig3h_data[,set:=factor(set,c('popDE/DR','eQTL-popDE/DR'))]

fwrite(fig3h_data,sprintf("%s/Fig3h_data.tsv",SOURCE_DATA_DIR),sep="\t")

fig3h_data[,celltype:=factor(celltype,rev(lineage_order))]

fig3h_plot<-ggplot(fig3h_data[type=="expression",],aes(frac_var_mean*100,celltype, size=N*100))+
  geom_point(color="black",pch=21,position=position_dodge(width=0))+
  geom_point(aes(fill=mediator),pch=21,position=position_dodge(width=0),alpha=0.6)+
  scale_fill_manual(values=color_mediator,guide='none')+
  scale_size_continuous(name="Significantly mediated popDEGs (%)",breaks = (c(10, 50, 100, 200, 500, 1000)/10),labels=paste0(c(10, 50, 100, 200, 500, 1000)/10,'%'),
		  guide=guide_legend(nrow=1,title.position="top",title.hjust=0.5))+
  facet_grid(cols=vars(state),rows=vars(set),labeller=labeller(set=c("popDE/DR"="All popDEGs","eQTL-popDE/DR"="eQTL popDEGs")))+
  scale_x_continuous(limits=(c(0,1)*100),breaks=(c(0,0.2,0.4,0.6,0.8,1)*100))+
  scale_y_discrete(labels=rev(c("MONO","B","CD4+ T","CD8+ T","NK")),breaks=rev(c("MONO","B","T.CD4","T.CD8","NK")),position="left")+
  xlab("Average popDE explained (%)")+
  theme_plot()+
  theme(text=element_text(size=7),panel.spacing=unit(0,"mm"),
				axis.text.x=element_text(angle=0,hjust=0.5),legend.title=element_text(size=5),axis.title.y=element_blank())#+
#  theme(legend.spacing.x= unit(1, 'mm'),legend.spacing.y=unit(2, 'mm'),legend.text=element_text(size=6),
#	legend.box = "vertical", legend.margin=margin(),legend.key.size = unit(1, 'mm'),axis.title.y=element_blank())

fig3h_legend<-cowplot::get_legend(fig3h_plot)
