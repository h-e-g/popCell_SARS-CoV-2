which_INTROGRESSED.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(dynamicTreeCut))
suppressMessages(library(BiocNeighbors))
suppressMessages(library(CelliD))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(ggrastr))
suppressMessages(library(scales))
suppressMessages(library(data.table))
suppressMessages(library(DropletUtils))
suppressMessages(library(SoupX))
suppressMessages(library(batchelor))
suppressMessages(library(kBET))
suppressMessages(library(HGC))
suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

EIP="/pasteur/zeus/projets/p02/evo_immuno_pop"
DAT_DIR=sprintf("%s/single_cell/project/pop_eQTL/data",EIP)
PAP_DIR=sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V2",EIP)
FIG_DIR=sprintf("%s/figureMaterial",PAP_DIR)
RES_DIR=sprintf("%s/single_cell/resources",EIP)

source(sprintf("%s/template_scripts/processing_pipeline/00_set_colors.R",RES_DIR))
source(sprintf("%s/template_scripts/querySNPs.R",RES_DIR))
source(sprintf("%s/template_scripts/processing_pipeline/09_eQTLmapping/ZZ_load_eQTLs.R",RES_DIR))

################################################################################

#eqtl_type="eQTL"
eqtl_type="reQTL"
#filter="signif"
filter="mash"

#snp_metadata=fread(sprintf("%s/popCell_data/02_ImputedGenotypes/Imputed/b38/Map_allChr_b38_imputed_filtered.DR2_90pct.AF_1pct.tsv.gz",EIP))
snp_metadata=getMap(annotate=TRUE)
#setnames(snp_metadata,"#CHROM","CHROM")
setnames(snp_metadata,"ID","rsID")
setnames(snp_metadata,"chromposID_hg38","ID")


# SNPs of archaic origin (high confidence)
aSNP_good=fread('/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/Final_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt')
CHS_LD=fread("/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/CHS_LDblocks_GOOD_aSNPs_popCellSNVs_nomono_hg38_biSNPs_REF2AA_phased.txt")
CEU_LD=fread("/pasteur/zeus/projets/p02/evo_immuno_pop/users/Javier/results/asnps/CEU_LDblocks_GOOD_aSNPs_popCellSNVs_nomono_hg38_biSNPs_REF2AA_phased.txt")
# introgressed SNPs (high confidence)
aSNP_lenient=fread(sprintf('%s/users/Javier/results/asnps/Final_lenient_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt',EIP))

# TODO : for 740 SNPs, FILTER_ORIGIN is AMH BUT the SNP is a high confidence aSNP , why ? (740 out of >288k introgressed SNPs with FILTER_ORIGIN==AMH),
# aSNP_good[ID%chin% aSNP_lenient[FILTER_ORIGIN=='AMH' & ID%chin%aSNP_good$ID,ID]]

if (tolower(filter)=="signif"){
  if (eqtl_type=="eQTL"){
    eqtl_file=eQTL_Signif_both
  } else if (eqtl_type=="reQTL"){
    eqtl_file=reQTL_Signif_both
  }
} else if (filter=="mash"){
  if (eqtl_type=="eQTL"){
    eqtl_file_lineage=eQTL_Stats_lineage_mash
    eqtl_file_celltype=eQTL_Stats_celltype_mash
    colnames(eqtl_file_lineage)%<>%{gsub("lineage","celltype",.)}
    eqtl_file=rbind(eqtl_file_celltype,eqtl_file_lineage[,mget(colnames(eqtl_file_celltype))])
  } else if (eqtl_type=="reQTL"){
    eqtl_file_lineage=reQTL_Stats_lineage_mash
    eqtl_file_celltype=reQTL_Stats_celltype_mash
    colnames(eqtl_file_lineage)%<>%{gsub("lineage","celltype",.)}
    eqtl_file=rbind(eqtl_file_celltype,eqtl_file_lineage[,mget(colnames(eqtl_file_celltype))])
  }
}
eqtl_file=merge(eqtl_file,snp_metadata[,.(snps=rsID,REF,ALT,ANCESTRAL,AF,AFB,EUB,ASH,AF_global,ID,covid_A2_pval,covid_A2_beta,covid_B2_pval,covid_B2_beta,covid_C2_pval,covid_C2_beta)],by=c("snps"),all.x=T)

eqtl_file=merge(eqtl_file,aSNP_lenient,by=c('ID','REF','ALT'))
# DONE: for each eQTL, identify which allele is present on the archaic haplotype. Files have been created by Javier In
# ./users/Javier/results/archaics_eqtls/Vindija33.19_gtformat.txt
# ./users/Javier/results/archaics_eqtls/Denisova_gtformat.txt
# check that REF ALT match when you’re quering these files


fwrite(eqtl_file,sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/%s_%s_both_ARCHAIC_long.tsv.gz",EIP,eqtl_type,filter),sep="\t")


AllelesV=fread(sprintf("%s/users/Javier/results/archaics_eqtls/Vindija33.19_gtformat.txt",EIP))
colnames(AllelesV)=c('CHROM','POS','REF','ALT','V1','V2')
AllelesV[,Nb_ALT_Vindija:=ifelse(ALT!='.',grepl(unique(V1),ALT) + grepl(unique(V2),ALT) ,0),by=.(V1,V2)]
AllelesV[,ID:=paste(CHROM,POS,sep=':')]

AllelesD=fread(sprintf("%s/users/Javier/results/archaics_eqtls/Denisova_gtformat.txt",EIP))
colnames(AllelesD)=c('CHROM','POS','REF','ALT','D1','D2')
AllelesD[,Nb_ALT_Denisova:=ifelse(ALT!='.',grepl(unique(D1),ALT) + grepl(unique(D2),ALT) ,0),by=.(D1,D2)]
AllelesD[,ID:=paste(CHROM,POS,sep=':')]

# Q: how to assess
# 1. which is the nearest aSNP in LD
# 2. which allele is on the archaic haplotype
AllelesVD=merge(AllelesV,AllelesD,by=c('ID','CHROM','POS','REF'),suffix=c('.V','.D'))
AllelesVD=AllelesVD[order(CHROM,POS)]
AllelesVD[,ALT.VD:=ifelse(ALT.V=='.' & ALT.D=='.','.',paste(setdiff(c(str_split(ALT.V,',',simplify=T),str_split(ALT.D,',',simplify=T)),'.'),collapse=',')),by=.(ALT.V,ALT.D)]
#AllelesVD[,.N,by=.(Nb_ALT_Vindija,Nb_ALT_Denisova)][order(-N)]
rm(AllelesV,AllelesD)
## read aSNPs
# asnps <- fread("results/asnps/Final_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased_GOOD_R2haps_10Kb.txt")
# asnps[asnps$ID=="20:63644491",]


##############################################################################################################
################             strong candidates  natural  selection                 ###########################
##############################################################################################################

eqtl_type='eQTL'
filter='Signif'
eqtl_file_archaic=fread(sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/%s_%s_both_ARCHAIC_long.tsv.gz",EIP,eqtl_type,filter))
eqtl_file_archaic[,type:='eQTL']

eqtl_type='reQTL'
filter='Signif'
reqtl_file_archaic=fread(sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/%s_%s_both_ARCHAIC_long.tsv.gz",EIP,eqtl_type,filter))
reqtl_file_archaic[,type:='reQTL']


q99=fread(sprintf("%s/users/Javier/results/asnps/Percentiles_Final_aSNPs_SPrimeorCRF_popCellSNVs_biSNPs_nodups_nomono_hg38_REF2AA_phased.txt",EIP))[METHOD=="MATCH_CRF/S'" & ORIGIN=='ANY']
####################################
all_ARCHAIC=rbind(eqtl_file_archaic,reqtl_file_archaic)
all_ARCHAIC=merge(all_ARCHAIC,AllelesVD[,.(ID,REF,Nb_ALT_Vindija,Nb_ALT_Denisova,ALT.VD)],by=c('ID','REF'),suffix=c('','.archaic'),all.x=T)

library(tictoc)
hap_origin=list()
for(pop in c('CEU','CHS')){
  hap_origin[[pop]]=list()
 for(chr in 1:22){
   tic(paste(pop,chr))
    ## read all snps
    df <- fread(paste0(EIP,"/users/Javier/data/POS/CHS_CEU_YRI_chr",chr,"_popCellSNVs_nomono_hg38.pos.txt"))
    colnames(df) <- c("CHROM","POS","REF","ALT")
    df$ID <- paste0(df$CHROM,":",df$POS)

    ## read LD files
    ld <- fread(paste0(EIP,"/users/Javier/results/ld/",pop,"_chr",chr,"_popCellSNVs_nomono_hg38.ld.gz"),h=T,stringsAsFactors = F, data.table = F)
    ld_b <- ld[,c("CHR_B","BP_B","SNP_B","CHR_A","BP_A","SNP_A","R2")] ## to make symmetric
    colnames(ld_b) <- colnames(ld)
    ld <- rbind(ld,ld_b)
    ld <- ld[ld$SNP_A%chin%all_ARCHAIC$ID,]
    rm(ld_b)

    ## add intro info
    ld$ID <- ld$SNP_B
    ld <- merge(ld,aSNP_good[,c("ID","ORIGIN")],by="ID")
    ld$ID <- NULL
    ld$ORIGIN <- ifelse(is.na(ld$ORIGIN),"AMH",ld$ORIGIN)
    hap_origin[[pop]][[chr]]=ld
    toc()
  }
  hap_origin[[pop]]=rbindlist(hap_origin[[pop]])
}
hap_origin=rbindlist(hap_origin,idcol='pop')
setnames(hap_origin,'ORIGIN','ORIGIN_B')

freqs=fread(sprintf("%s/users/Javier/results/freq/CHS_CEU_YRI_allCHR_popCellSNVs_nomono_hg38.frq.strat.gz",EIP))
# A1 = ALT alleles
# A2 = REF allele
# MAF = frequency of ALT allele

freqs_B=freqs[SNP%chin%hap_origin$SNP_B,]
prov=merge(hap_origin,freqs_B[,.(SNP_B=SNP,ALT=A1,REF=A2,ALT_FREQ=MAF,pop=CLST)],by=c('SNP_B','pop')) #
prov2=merge(prov,freqs_B[CLST=='YRI',.(SNP_B=SNP,ALT_FREQ_YRI=MAF)],by=c('SNP_B'))
prov3=merge(prov2,AllelesVD[,.(SNP_B=ID,REF,Nb_ALT_Vindija,Nb_ALT_Denisova,ALT.VD)],all.x=TRUE,by=c('SNP_B','REF')) #
prov3[,INTROGRESSED_ALLELE_B:=ifelse(ALT_FREQ_YRI>.5,REF,ALT)]
prov3[,INTROGRESSED_ALLELE_FREQ:=ifelse(INTROGRESSED_ALLELE_B==ALT,ALT_FREQ,1-ALT_FREQ)]
prov3[,INTROGRESSED_ALLELE_TAG_ALTREF:=ifelse(INTROGRESSED_ALLELE_B==ALT,"ALT","REF")]

introgressed_annot=prov3[order(SNP_A,!SNP_B%chin%snp_metadata$ID,-R2,-INTROGRESSED_ALLELE_FREQ),][!duplicated(paste(SNP_A,pop)),]
introgressed_annot=introgressed_annot[,.(ID=SNP_A, Population=pop,
                                      tag_aSNP=SNP_B,
                                      r2_tag_aSNP=R2,
                                      INTROGRESSED_ALLELE_TAG=INTROGRESSED_ALLELE_B,
                                      INTROGRESSED_ALLELE_TAG_ALTREF,
                                      ALT.VD_tag=ALT.VD,
                                      Nb_ALT_V_tag=Nb_ALT_Vindija,
                                      Nb_ALT_D_tag=Nb_ALT_Denisova,
                                      INTROGRESSED_ALLELE_FREQ)]

all_ARCHAIC=merge(all_ARCHAIC,introgressed_annot,by=c('ID','Population'),all.x=T)

FSTs=fread(sprintf("%s/users/Javier/results/pbs/FST_CEU_CHS_YRI_allCHR_popCellSNVs_nomono_hg38.pbs.txt.gz",EIP))
setnames(FSTs,'SNP','ID')
all_ARCHAIC=merge(all_ARCHAIC,FSTs,by='ID',all.x=TRUE)

#### annotate which genes are COV VIPs
COV_VIPs=fread(sprintf("%s/single_cell/resources/references/Souilmi_SARS-CoV-2-VIPs/SARS_CoV_2_VIPs_332.txt",EIP))
all_ARCHAIC[,COV_VIPs:=ifelse(gene_name%chin%COV_VIPs$`HGNC symbol`,'yes','')]

#### annotate which genes are immune genes (inborn errors)
IEI_table=fread(sprintf('%s/single_cell/resources/references/Inborn_Error_Immunity_2022/IUIS_IEI_website_April_2022.txt',EIP))
all_ARCHAIC[,IEI:=ifelse(gene_name%chin%IEI_table$`Genetic defect`,'yes','')]
all_ARCHAIC[,IEI_type:=ifelse(gene_name%chin%IEI_table$`Genetic defect`,gsub("Table [0-9]+","",IEI_table$`Major category`)[match(gene_name,IEI_table$`Genetic defect`)],'')]
all_ARCHAIC[,IEI_Disease:=ifelse(gene_name%chin%IEI_table$`Genetic defect`,gsub("Table [0-9]+","",IEI_table$`Disease`)[match(gene_name,IEI_table$`Genetic defect`)],'')]
all_ARCHAIC[,IEI_Features:=ifelse(gene_name%chin%IEI_table$`Genetic defect`,gsub("Table [0-9]+","",IEI_table$`Associated features`)[match(gene_name,IEI_table$`Genetic defect`)],'')]

allGOterms=fread(sprintf('%s/single_cell/resources/references/PATHWAYS/allGOterms_12393_expressed_genes.tsv',EIP))
all_ARCHAIC[,GO_immune:=ifelse(gene_name%chin%allGOterms[go=='GO:0006955',unique(gene)],'yes','')]

#### annotate which genes are immune genes (Schoggins Cur Op Vir)
ISG_effector_table=fread(sprintf('%s/single_cell/resources/references/PATHWAYS/SchogginsCurrOpVir2011.tsv',EIP))
all_ARCHAIC[,ISG_effector:=ifelse(gene_name%chin%ISG_effector_table$`Gene symbol`,'yes','')]
all_ARCHAIC[,Targeted_viruses:=ifelse(gene_name%chin%ISG_effector_table$`Gene symbol`,ISG_effector_table$`Targeted viruses`[match(gene_name,ISG_effector_table$`Gene symbol`)],'')]
all_ARCHAIC[,EffectOn:=ifelse(gene_name%chin%ISG_effector_table$`Gene symbol`,ISG_effector_table$`Viral life cycle`[match(gene_name,ISG_effector_table$`Gene symbol`)],'')]

# all_ARCHAIC[order(P,pvalue),][P<.05][!duplicated(paste(gene_name,POP)),][METRIC=='PBS' & (IEI=='yes' & grepl('vir|Vir|influenza|COV',IEI_Features)==TRUE | ISG_effector=='yes'),][,.(POP,gene_name,ISG_effector,IEI,P,type,snps,celltype,state,ASH,AFB,EUB,pvalue,covid_A2_pval,covid_B2_pval,covid_C2_pval)][order(POP,P),]
#
# all_ARCHAIC[order(P,pvalue),][P<.05][!duplicated(paste(gene_name,POP)),][METRIC=='ABS_NSL' & (IEI=='yes' & grepl('vir|Vir|influenza|COV',IEI_Features)==TRUE | ISG_effector=='yes'),][,.(POP,gene_name,ISG_effector,IEI,P,type,snps,celltype,state,ASH,AFB,EUB,pvalue,covid_A2_pval,covid_C2_pval)][order(POP,P),]
#
##### annotate allele Age in the considered population
# all_ARCHAIC[,Age_lowerCI:=case_when(Population=='CHS'~Age_lowerCI_CHS,
#                                     Population=='CEU'~Age_lowerCI_CEU,
#                                     Population=='YRI'~Age_lowerCI_YRI)]
# all_ARCHAIC[,Age_upperCI:=case_when(Population=='CHS'~Age_upperCI_CHS,
#                                      Population=='CEU'~Age_upperCI_CEU,
#                                      Population=='YRI'~Age_upperCI_YRI)]
# all_ARCHAIC[,Age_log10pvalue:=case_when(Population=='CHS'~Age_log10pvalue_CHS,
#                                     Population=='CEU'~Age_log10pvalue_CEU,
#                                     Population=='YRI'~Age_log10pvalue_YRI)]

#### Define genes that are up regulated upon stimulation ####

DR=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/suppTables/TableS2/TableS2D_logFC_perlineage.tsv.gz",EVO_IMMUNO_POP_ZEUS))
DR= melt(DR) %>% mutate(sumstat=gsub('(.*)_logFC_(.*)_(.*)','\\1',variable)) %>%
            mutate(celltype=gsub('(.*)_logFC_(.*)_(.*)','\\2',variable)) %>%
            mutate(state=gsub('(.*)_logFC_(.*)_(.*)','\\3',variable)) %>%
            dcast(ID+Symbol+celltype+state~sumstat,value.var='value')
 setnames(DR,'mean','logFC')
DR[,FDR:=p.adjust(P,'fdr')]


DR_celltype=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/suppTables/TableS2/TableS2C_logFC_perCelltype.tsv.gz",EVO_IMMUNO_POP_ZEUS))
DR_celltype= melt(DR_celltype) %>% mutate(sumstat=gsub('(.*)_logFC_(.*)_(.*)','\\1',variable)) %>%
            mutate(celltype=gsub('(.*)_logFC_(.*)_(.*)','\\2',variable)) %>%
            mutate(state=gsub('(.*)_logFC_(.*)_(.*)','\\3',variable)) %>%
            dcast(ID+Symbol+celltype+state~sumstat,value.var='value')
 setnames(DR_celltype,'mean','logFC')
DR_celltype[,FDR:=p.adjust(P,'fdr')]

DR=rbind(DR,DR_celltype)
DR=dcast(DR,ID+Symbol+celltype~state,value.var=list('logFC','FDR'))
setnames(DR,c('ID','Symbol'),c('gene','gene_name'),skip_absent=TRUE)
DR[,response:=case_when(logFC_COV>.5 & FDR_COV<.01 & logFC_IAV>.5 & FDR_IAV<.01~'up both',
                                  logFC_COV>.5 & FDR_COV<.01 & (abs(logFC_IAV)<.5 | FDR_IAV<.01)~'up COV',
                                  logFC_IAV>.5 & FDR_IAV<.01 & (abs(logFC_COV)<.5 | FDR_COV<.01)~'up IAV',
                                  logFC_COV< -.5 & FDR_COV<.01 & logFC_IAV>.5 & FDR_IAV<.01~'divergent IAV',
                                  logFC_COV> .5 & FDR_COV<.01 & logFC_IAV< -.5 & FDR_IAV<.01~'divergent COV',
                                  logFC_COV< -.5 & FDR_COV<.01 & logFC_IAV< -.5 & FDR_IAV<.01~'down both',
                                  logFC_COV< -.5 & FDR_COV<.01 & (abs(logFC_IAV)<.5 | FDR_IAV<.01)~'down COV',
                                  logFC_IAV< -.5 & FDR_IAV<.01 & (abs(logFC_COV)<.5 | FDR_COV<.01)~'down IAV',
                                  TRUE~'')]

all_ARCHAIC=merge(all_ARCHAIC,DR,by=c('gene','gene_name','celltype'),all.x=TRUE)
#### Define the effect of selected alleles at eQTL/ reQTL on response to IAV/COV

all_ARCHAIC[,effect_of_ALT_on_response_COV:=ifelse(FDR_COV<.01,sign(beta)*sign(logFC_COV),0)]
all_ARCHAIC[,effect_of_ALT_on_response_IAV:=ifelse(FDR_IAV<.01,sign(beta)*sign(logFC_IAV),0)]

all_ARCHAIC[,unique(snps)]
all_ARCHAIC[,unique(tag_aSNP)%in%snp_metadata$ID]
all_ARCHAIC[,rsID_tag_aSNP:=snp_metadata[match(tag_aSNP,ID),rsID],by=tag_aSNP]
Geno_tag_aSNP=getSNP(all_ARCHAIC[,na.omit(unique(rsID_tag_aSNP))],Map=snp_metadata[,.(ID=rsID,CHROM,POS)])
fwrite(Geno_tag_aSNP,file=sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/All_eQTL_tag_aSNP_genotypes.tsv.gz",EIP),sep='\t')

Geno_tag_aSNP=fread(sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/All_eQTL_tag_aSNP_genotypes.tsv.gz",EIP))

Geno_eQTL=fread(sprintf("%s/single_cell/project/pop_eQTL/data/3_eQTL_mapping/All_eQTL_and_reQTL_genotypes.tsv.gz",EIP))

all_ARCHAIC[,which_INTROGRESSED:='NA']
SNP_pairs=all_ARCHAIC[!duplicated(paste(snps,rsID_tag_aSNP,Population)),.(snps,rsID_tag_aSNP,Population,INTROGRESSED_ALLELE_TAG_ALTREF)]
for( i in 1:SNP_pairs[,.N]){
  eQTL_rsID=SNP_pairs[i,snps]
  tag_aSNP_rsID=SNP_pairs[i,rsID_tag_aSNP]
  POP_target=SNP_pairs[i,Population]
  cat(i,eQTL_rsID,tag_aSNP_rsID,POP_target,':')
  POP_EIP=c('CEU'='EUB','YRI'='AFB','CHS'='ASH')
  Geno_test=merge(Geno_eQTL[ID==eQTL_rsID & grepl(POP_EIP[POP_target],IID),],Geno_tag_aSNP[ID==tag_aSNP_rsID & grepl(POP_EIP[POP_target],IID),],by=c('IID'),suffix=c('.eQTL','.aSNP'))
  r_aSNP_eQTL=Geno_test[,cor(Number_of_ALT_alelle.eQTL,Number_of_ALT_alelle.aSNP)]
  wINTROGRESSED=ifelse(r_aSNP_eQTL>0,SNP_pairs[i,INTROGRESSED_ALLELE_TAG_ALTREF],ifelse(SNP_pairs[i,INTROGRESSED_ALLELE_TAG_ALTREF]=='REF','ALT','REF'))
    cat(r_aSNP_eQTL,wINTROGRESSED,'\n')
  all_ARCHAIC[snps==eQTL_rsID & rsID_tag_aSNP==tag_aSNP_rsID & Population==POP_target, which_INTROGRESSED:=wINTROGRESSED]
}
all_ARCHAIC[which_INTROGRESSED=='NA' & snps=='rs117961446',which_INTROGRESSED:='ALT'] # manual annotation bease on ALT allele frequency of the eQTL/tag aSNP (~5% each in EAS)
all_ARCHAIC[which_INTROGRESSED=='NA' & snps=='rs75801422',which_INTROGRESSED:='ALT'] # manual annotation bease on allele freqiency of the eQTL/tag aSNP (~5% each in EAS)




# all_ARCHAIC[,which_INTROGRESSED:=case_when(POP=='CHS' & P<.05~ifelse(ASH>EUB & ASH>AFB,'ALT','REF'),
#                                         POP=='CEU' & P<.05~ifelse(EUB>ASH & EUB>AFB,'ALT','REF'),
#                                         POP=='YRI' & P<.05~ifelse(AFB>ASH & AFB>EUB,'ALT','REF'),
#                                         TRUE~'NA')]
all_ARCHAIC[,effect_of_INTROGRESSED_on_response_COV:=ifelse(which_INTROGRESSED=='ALT',effect_of_ALT_on_response_COV,-effect_of_ALT_on_response_COV)]
all_ARCHAIC[,effect_of_SEL_on_response_IAV:=ifelse(which_INTROGRESSED=='ALT',effect_of_ALT_on_response_IAV,-effect_of_ALT_on_response_IAV)]

all_ARCHAIC[,effect_of_INTROGRESSED_on_A2:=ifelse(covid_A2_pval>0.01,0,ifelse(which_INTROGRESSED=='ALT',sign(covid_A2_beta),-sign(covid_A2_beta)))]
all_ARCHAIC[,effect_of_SEL_on_C2:=ifelse(covid_C2_pval>0.01,0,ifelse(which_INTROGRESSED=='ALT',sign(covid_C2_beta),-sign(covid_C2_beta)))]

# check results in CHS (considering only the cell state with strongest effect fopr each gene)
# all_ARCHAIC[ METRIC=='PBS' & P<.01 & type=='reQTL' & state=='COV' & POP=='CHS',][order(gene_name,pvalue)][!duplicated(gene_name),.(.N,mean(effect_of_SEL_on_response_COV>0))),keyby=.(response)]
# all_ARCHAIC[ METRIC=='PBS' & P<.05 & type=='reQTL' & state=='COV' & POP=='CHS',][order(gene_name,pvalue)][!duplicated(gene_name),.(.N,mean(effect_of_SEL_on_response_COV>0))),keyby=.(response)]
# all_ARCHAIC[ METRIC=='PBS' & P<.01 & type=='reQTL' & POP=='CHS',.(.N,mean(effect_of_SEL_on_response_COV>0),mean(effect_of_SEL_on_response_COV<0)),keyby=.(response)]

freqs_wide=merge(freqs[SNP%chin%all_ARCHAIC$ID,],all_ARCHAIC[,.(ID,REF,ALT,ANCESTRAL,which_INTROGRESSED,Population)],by.x='SNP',by.y='ID',all.x=TRUE,allow.cartesian=TRUE)
freqs_wide[,DAF_or_MAF:=ifelse(A2==toupper(ANCESTRAL),MAF,ifelse(A1==toupper(ANCESTRAL),1-MAF,MAF))]
freqs_wide[,ALT_DERANC:=ifelse(A2==toupper(ANCESTRAL),'DER',ifelse(A1==toupper(ANCESTRAL),"ANC","DER"))]
freqs_wide[,INTROGRESSED_DERANC:=case_when(which_INTROGRESSED=='ALT' & A2==toupper(ANCESTRAL)~'DER',
                                       which_INTROGRESSED=='ALT' & A1==toupper(ANCESTRAL)~'ANC',
                                       which_INTROGRESSED=='ALT'~'DER',
                                       which_INTROGRESSED=='REF' & A2==toupper(ANCESTRAL)~'ANC',
                                       which_INTROGRESSED=='REF' & A1==toupper(ANCESTRAL)~'DER',
                                       which_INTROGRESSED=='REF'~'ANC')]
freqs_wide=dcast(unique(freqs_wide),CHR+SNP+A1+A2+Population+which_INTROGRESSED+ANCESTRAL+INTROGRESSED_DERANC+ALT_DERANC~CLST,value.var=c('MAF','DAF_or_MAF','MAC','NCHROBS'))

all_ARCHAIC=merge(all_ARCHAIC,freqs_wide[,.(SNP,A2,Population,INTROGRESSED_DERANC,ALT_DERANC,DAF_or_MAF_CEU,DAF_or_MAF_CHS,DAF_or_MAF_YRI)],by.x=c('ID','REF','Population'),by.y=c('SNP','A2','Population'),all.x=TRUE)


############################ (weak evidence P select <0.01 , no local enrichment) & colocalize with COVID

Coloc_table_full=fread(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V3/figureMaterial/Fig3//colocalization_COVID/allTraits_allStates_all_eQTLs_near_1e5_covid_peak_coloc_signal_v2.tsv.gz",EIP))
Coloc_table_full=Coloc_table_full[Dist=="250kb",]
Coloc_table_full[,type:=ifelse(grepl('logFC',run),'reQTL','eQTL')]
Coloc_table_full[,snps:=gsub('.*220409_.*_(.*)','\\1',ID)]
setnames(Coloc_table_full,'symbol','gene_name')
Coloc_table_full_wide=dcast(Coloc_table_full[,.(snps, gene, gene_name, type, celltype, state,trait,PP.H3.abf,PP.H4.abf)],snps+gene+gene_name+type+celltype+state~trait,value.var=list('PP.H3.abf','PP.H4.abf'))
all_ARCHAIC=merge(all_ARCHAIC,Coloc_table_full_wide,by=c('gene','gene_name','celltype','state','snps','type'),all.x=TRUE)
all_ARCHAIC[,tested_coloc:=ifelse(!is.na(PP.H4.abf_critical) | !is.na(PP.H4.abf_hospitalized) | !is.na(PP.H4.abf_reported),'yes','')]
all_ARCHAIC[,colocalized:=ifelse(tested_coloc=='','',ifelse(PP.H4.abf_critical>.8 | PP.H4.abf_hospitalized>.8 | PP.H4.abf_reported>.8,'yes',''))]

############################ (weak evidence P select <0.01 , P local enrichment < 0.01 )
all_ARCHAIC[,FreqARCHAIC:=case_when(INTROGRESSED_DERANC=='DER' & Population=='CEU'~DAF_or_MAF_CEU,
                                    INTROGRESSED_DERANC=='DER' & Population=='CHS'~DAF_or_MAF_CHS,
                                    INTROGRESSED_DERANC=='ANC' & Population=='CEU'~1-DAF_or_MAF_CEU,
                                    INTROGRESSED_DERANC=='ANC' & Population=='CHS'~1-DAF_or_MAF_CHS)]
all_ARCHAIC[,maxFreq_perGene:=max(FreqARCHAIC),by=gene_name]
prov_CHS=merge(all_ARCHAIC[Population=='CHS',],CHS_LD[,.(ID,Length_kb=LEN/1000,START,END,NUM_aSNPS_IN_HAP)],by.x='tag_aSNP',by.y='ID',all.x=TRUE)
prov_CEU=merge(all_ARCHAIC[Population=='CEU',],CEU_LD[,.(ID,Length_kb=LEN/1000,START,END,NUM_aSNPS_IN_HAP)],by.x='tag_aSNP',by.y='ID',all.x=TRUE)
all_ARCHAIC=rbind(prov_CHS,prov_CEU)
all_ARCHAIC=all_ARCHAIC[order(-maxFreq_perGene,pvalue),]


all_ARCHAIC[,aSNP_perKb:=NUM_aSNPS_IN_HAP/Length_kb]
all_ARCHAIC[,FreqARCHAIC_CEU:=ifelse(INTROGRESSED_DERANC=='DER',DAF_or_MAF_CEU,1-DAF_or_MAF_CEU)]
all_ARCHAIC[,FreqARCHAIC_CHS:=ifelse(INTROGRESSED_DERANC=='DER',DAF_or_MAF_CHS,1-DAF_or_MAF_CHS)]
all_ARCHAIC[,FreqARCHAIC_YRI:=ifelse(INTROGRESSED_DERANC=='DER',DAF_or_MAF_YRI,1-DAF_or_MAF_YRI)]
all_ARCHAIC[,ALLELE_ORIGIN:=MATCH]
all_ARCHAIC[FreqARCHAIC_YRI<.01 & ALLELE_ORIGIN=='AMH',ALLELE_ORIGIN:='UNKNOWN_ARCHAIC']
all_ARCHAIC[,reintrogressed:=ifelse(ALLELE_ORIGIN=='AMH' | INTROGRESSED_DERANC=='ANC','yes','')]

all_ARCHAIC=merge(all_ARCHAIC,q99[,.(Population,q95,q99)],all.x=TRUE)
fwrite(all_ARCHAIC,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/all_ARCHAIC_genelist_v2.tsv",EIP),sep='\t')


rbind(all_ARCHAIC[colocalized=='yes',.(snps,gene_name)],gene_list_strongPBS_coloc[,.(snps,gene_name)])[!duplicated(paste(snps,gene_name)),.(snps,gene_name)]
#################################################################################################
########################## define sets of SNPs to analyse #######################################
#################################################################################################
dir.create(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/",EIP))
# gene_list_immune_introgressed=all_ARCHAIC[(IEI=='yes' | ISG_effector=='yes'| COV_VIPs=='yes'),]
# fwrite(gene_list_immune_introgressed,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/gene_list_immune_introgressed.tsv",EIP),sep='\t')
#
# gene_list_immuneGO_introgressed=all_ARCHAIC[(IEI=='yes' | ISG_effector=='yes' | GO_immune=='yes' | COV_VIPs=='yes'),]
# fwrite(gene_list_immuneGO_introgressed,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/geneList_ImmuneGO_introgressed.tsv",EIP),sep='\t')
#
# gene_list_immuneGO_introgressed_COVID=gene_list_immuneGO_introgressed[(covid_A2_pval<0.01 | covid_C2_pval<0.01 | covid_B2_pval<0.01),]
# fwrite(gene_list_immuneGO_introgressed_COVID,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/gene_list_immuneGO_introgressed_COVID.tsv",EIP),sep='\t')
#
gene_list_introgressed_coloc=all_ARCHAIC[(covid_A2_pval<0.01 | covid_C2_pval<0.01 | covid_B2_pval<0.01) & colocalized=='yes',]
fwrite(gene_list_introgressed_coloc,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/gene_list_introgressed_coloc.tsv",EIP),sep='\t')
#
# gene_list_archaic_immune=all_ARCHAIC[(IEI=='yes' | ISG_effector=='yes' | GO_immune=='yes' | COV_VIPs=='yes') & ALLELE_ORIGIN!='AMH',]
# fwrite(gene_list_archaic_immune,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/gene_list_archaic_immune.tsv",EIP),sep='\t')
#
# gene_list_adaptive_introgressed=all_ARCHAIC[FreqARCHAIC>q95,]
# fwrite(gene_list_adaptive_introgressed,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/gene_list_adaptive_introgressed.tsv",EIP),sep='\t')
#
gene_list_adaptive_introgressed_immune=all_ARCHAIC[ FreqARCHAIC>q95 & (IEI=='yes' | ISG_effector=='yes' | GO_immune=='yes' | COV_VIPs=='yes'),]
fwrite(gene_list_adaptive_introgressed_immune,sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/gene_list_adaptive_introgressed_immune.tsv",EIP),sep='\t')

# GWSignif_th_table=eQTL_Signif_all[,.(TH=max(pvalue)),by=.(celltype,state,type)]
# eQTL_nominal=merge(eQTL_nominal,GWSignif_th_table,by=c('celltype','state','type'))
# eQTL_nominal[,GWSignif_th:=pvalue <= TH]
# eQTL_nominal[,GWSignif_th_unique:= (pvalue <= median(GWSignif_th_table[,TH]))]



# define criteria for the list of genes to plot
eQTL_set=unique(gene_list_introgressed_coloc[,paste(gene_name,snps,sep=' - ')])

######################################################################################################
########################  Implement plotting function to visualize eQTLs  ############################
######################################################################################################

plot_geneList_ARCHAIC=function(gene_list,
                        addLineage=FALSE,
                        addLineageState=TRUE,
                        addLineageStateSummary=FALSE,
                        addImmuneAnnot=TRUE,
                        addLOGFC=TRUE,
                        addCOVID=TRUE,
                        addCOLOC=TRUE,
                        addARCHAIC=TRUE,
                        addLength=TRUE,
                        addFreqs=TRUE,
                        addAdaptive=TRUE,
                        mySize=2,legend=TRUE
                        ){
  eQTL_set=unique(gene_list[,paste(gene_name,snps,sep=' - ')])
  eQTL_set_allele=unique(gene_list[,paste0(gene_name,' (',snps,'-',ifelse(which_INTROGRESSED=='ALT',ALT,REF),')')])
    library(ggnewscale)

    p <- ggplot()
    p <- p + scale_y_continuous(breaks=match(eQTL_set,eQTL_set), labels=eQTL_set_allele)
    p <- p + theme_yann()+theme(panel.border = element_blank(),axis.ticks.y = element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + xlab('') + ylab('')

    # add lineage
    x_pos_min=0
    x_pos_annot=c()
    x_pos_labels=c()

    #########@ extract list of conditions/celltypes where the eQTL is present
    # nominal eQTL
    eQTL_nominal=rbind(eQTL_Stats_celltype_mash[,.(snps,gene_name,celltype,state,beta,pvalue)],eQTL_Stats_lineage_mash[,.(snps,gene_name,celltype,state,beta,pvalue)])[pvalue<0.01,]
    eQTL_nominal[,type:='eQTL']
    reQTL_nominal=rbind(reQTL_Stats_celltype_mash[,.(snps,gene_name,celltype,state,beta,pvalue)],reQTL_Stats_lineage_mash[,.(snps,gene_name,celltype,state,beta,pvalue)])[pvalue<0.01,]
    reQTL_nominal[,type:='reQTL']
    eQTL_nominal=rbind(eQTL_nominal,reQTL_nominal)
    # GW signif eQTL
    eQTL_Signif_both[,type:='eQTL']
    reQTL_Signif_both[,type:='reQTL']
    eQTL_Signif_all=rbind(eQTL_Signif_both,reQTL_Signif_both)
    eQTL_nominal[,GWSignif:=paste(celltype,state,type,gene_name,snps)%chin%eQTL_Signif_all[,paste(celltype,state,type,gene_name,snps)]]
    # add lineage annotations
    eQTL_nominal[,lineage:=ifelse(celltype%in%lineage_5,celltype,celltype_2lineage[match(eQTL_nominal$celltype,celltype_2lineage$celltype),lineage])]
    # add DER/ANC annotation
    eQTL_nominal=merge(eQTL_nominal,unique(all_ARCHAIC[,.(snps,which_INTROGRESSED)]),by='snps',all.x=TRUE,allow.cartesian=TRUE)

    color_lineage=color_cellTypes_6level[lineage_5]
    color_celltype=c(color_lineage,color_cellTypes_24level[celltype_22])
    #color_archaics=c('AMH'=grey(0.7),'DENI'='#35978F','NEAND'='#F46D43','ARCHAIC'="#FEE08B")
    color_archaics=c('AMH'=grey(0.8),'DENI'='#80CDC1','NEAND'='#FDAE61','ARCHAIC'="#FEE08B",'UNKNOWN_ARCHAIC'=grey(0.5))
    if( addLineage ){
      #########@ assemble lineage panel
        DF_Lineage=eQTL_nominal[,.(snps,gene_name,beta,pvalue,lineage,state,type,GWSignif,celltype,which_INTROGRESSED)]
        DF_Lineage=DF_Lineage[order(snps,gene_name,lineage,-GWSignif,pvalue,type),][!duplicated(paste(snps,gene_name,lineage,state,type))]
        DF_Lineage_state=DF_Lineage[!duplicated(paste(snps,gene_name,lineage,state))]
        DF_Lineage=DF_Lineage[!duplicated(paste(snps,gene_name,lineage))]
        DF_Lineage[,x_pos:=match(lineage,lineage_order)]
        DF_Lineage[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_Lineage[,sign:=ifelse((beta*ifelse(which_INTROGRESSED=='ALT',1,-1))>0,'Increased','Decreased')]
        LINEAGE_GRID=data.table(expand.grid(x_pos=1:length(lineage_order),y_pos=1:length(eQTL_set)))

    #########@ add lineage panel to plot
        x_pos_annot=c(x_pos_annot,x_pos_min+1:length(lineage_order))
        x_pos_labels=c(x_pos_labels,lineage_order)
        p <- p + geom_tile(data=LINEAGE_GRID,aes(x_pos,y_pos),fill='white',col='black')
        #p <- p + geom_point(data=DF_Lineage[!is.na(y_pos),],aes(x_pos,y_pos,fill=lineage,alpha=(1+GWSignif)/2,pch=sign),size=2)
        p <- p + geom_point(data=DF_Lineage[!is.na(y_pos),],aes(x_pos,y_pos,fill=lineage,alpha=log10(-log10(pvalue)),pch=sign),size=mySize)
        p <- p + scale_fill_manual(values=color_lineage) + scale_shape_manual(values=c('Increased'=24,'Decreased'=25))
    }

    if( addLineageState ){
    #########@ assemble lineage x state panel
        DF_Lineage_state=eQTL_nominal[,.(snps,gene_name,beta,pvalue,lineage,state,type,GWSignif,celltype,which_INTROGRESSED)]
        DF_Lineage_state=DF_Lineage_state[order(snps,gene_name,lineage,-GWSignif,pvalue,type),][!duplicated(paste(snps,gene_name,lineage,state,type))]
        DF_Lineage_state=DF_Lineage_state[!duplicated(paste(snps,gene_name,lineage,state))]
        DF_Lineage_state[,x_pos:=match(paste(lineage,state,sep='-'),paste(rep(lineage_order,e=3),c('NS','COV','IAV'),sep='-'))]
        DF_Lineage_state[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_Lineage_state[,sign:=ifelse((beta*ifelse(which_INTROGRESSED=='ALT',1,-1))>0,'Increased','Decreased')]
        LINEAGE_STATE_GRID=data.table(expand.grid(x_pos=1:(3*length(lineage_order)),y_pos=1:length(eQTL_set)))

    #########@ add cell state panel to plot
        x_pos_min=pmax(max(x_pos_annot)+1,1)
        x_pos_annot=c(x_pos_annot,x_pos_min+1:15)
        x_pos_labels=c(x_pos_labels,paste(rep(lineage_order,e=3),c('NS','COV','IAV'),sep='-'))

        LINEAGE_STATE_GRID[,x_pos:=x_pos+x_pos_min]
        DF_Lineage_state[,x_pos:=x_pos+x_pos_min]
        p <- p + geom_tile(data=LINEAGE_STATE_GRID,aes(x_pos,y_pos),fill='white',col='black')
        p <- p + new_scale("fill")
        for (i in which(!is.na(DF_Lineage_state$y_pos))){
          p <- p + geom_tile(data=DF_Lineage_state[i,],aes(x_pos,y_pos,fill=lineage),alpha=.25,col='#000000')
         }
        p <- p + scale_fill_manual(values=color_lineage)
        p <- p + new_scale("fill") + scale_fill_manual(values=color_conditions)
        p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0.2,1))
        p <- p + geom_point(data=DF_Lineage_state[!is.na(y_pos),],aes(x_pos,y_pos,fill=state,alpha=log10(-log10(pvalue)),pch=sign),size=mySize)
        p <- p + scale_shape_manual(values=c('Increased'=24,'Decreased'=25))
    }

    if( addLineageStateSummary ){
    #########@ assemble lineage x state panel
        DF_Lineage_state_sum=eQTL_nominal[,.(snps,gene_name,beta,pvalue,celltype,state,type,GWSignif,lineage,which_INTROGRESSED)]
        DF_Lineage_state_sum=DF_Lineage_state_sum[order(snps,gene_name,type,pvalue,celltype,state),][!duplicated(paste(snps,gene_name,type))]
        DF_Lineage_state_sum=DF_Lineage_state_sum[!duplicated(paste(snps,gene_name,type))]
        DF_Lineage_state_sum[,x_pos:=1]
        DF_Lineage_state_sum[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_Lineage_state_sum[,sign:=ifelse((beta*ifelse(which_INTROGRESSED=='ALT',1,-1))>0,'Increased','Decreased')]
        LINEAGE_STATE_GRID=data.table(expand.grid(x_pos=1,y_pos=1:length(eQTL_set)))

    #########@ add cell state panel to plot
        x_pos_min=pmax(max(x_pos_annot)+1,1)
        x_pos_annot=c(x_pos_annot,x_pos_min+1)
        x_pos_labels=c(x_pos_labels,'Lineage, celltype & state\nwhere the eQTL is strongest')

        LINEAGE_STATE_GRID[,x_pos:=x_pos+x_pos_min]
        DF_Lineage_state_sum[,x_pos:=x_pos+x_pos_min]
        p <- p + geom_tile(data=LINEAGE_STATE_GRID,aes(x_pos,y_pos),fill='white',col='black')
        p <- p + new_scale("fill") + geom_tile(data=DF_Lineage_state_sum[!is.na(y_pos),],aes(x_pos,y_pos,fill=celltype),alpha=.25,col='#000000')
        p <- p + scale_fill_manual(values=color_celltype)
        p <- p + new_scale("fill") + geom_point(data=DF_Lineage_state_sum[!is.na(y_pos),],aes(x_pos,y_pos,fill=state,pch=sign),size=mySize)
        p <- p + scale_fill_manual(values=color_conditions)+ scale_shape_manual(values=c('Increased'=24,'Decreased'=25))
    }
    if( addImmuneAnnot ){
    #########@ assemble immune panel
        DF_IMMUNE=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,][!duplicated(paste(gene_name,snps)),.(gene_name,snps,IEI,COV_VIPs,GO_immune)]
        DF_IMMUNE=melt(DF_IMMUNE,id.vars=c('gene_name','snps'))
        DF_IMMUNE[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_IMMUNE[,x_pos:=match(variable,c('IEI','COV_VIPs','GO_immune'))]
        # IMMUNE_GRID=data.table(expand.grid(x_pos=1:length(lineage_order),y_pos=1:length(eQTL_set)))

    #########@ add immune group panel to plot
        x_pos_min=pmax(max(x_pos_annot)+1,1)
        x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
        x_pos_labels=c(x_pos_labels,c('Inborn Errors Immunity','COV VIPs','GO: immune response'))

        DF_IMMUNE[,x_pos:=x_pos+x_pos_min]
        p <- p + geom_tile(data=DF_IMMUNE,aes(x_pos,y_pos),fill='white',col='black')
        for (i in which(DF_IMMUNE$value=='yes')){
          p <- p + geom_tile(data=DF_IMMUNE[i,],aes(x_pos,y_pos),fill='black')
        }
    }
    if( addLOGFC ){

    }
    if( addARCHAIC){
      #########@ assemble Archaic panel
          DF_ARCHAIC=unique(gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,][,.(gene_name,snps,ALLELE_ORIGIN,INTROGRESSED_DERANC,Population,HAP_ORIGIN)])
          DF_ARCHAIC=dcast(DF_ARCHAIC,gene_name+snps+ALLELE_ORIGIN+INTROGRESSED_DERANC~Population,fill='',value.var='HAP_ORIGIN')
          DF_ARCHAIC=melt(DF_ARCHAIC,id.vars=c('gene_name','snps','INTROGRESSED_DERANC'))
          DF_ARCHAIC[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
          DF_ARCHAIC[,x_pos:=match(variable,c('ALLELE_ORIGIN','CEU','CHS'))]
          # IMMUNE_GRID=data.table(expand.grid(x_pos=1:length(lineage_order),y_pos=1:length(eQTL_set)))
          DF_ARCHAIC[,sign:=ifelse(INTROGRESSED_DERANC=='DER','Increased','Decreased')]
      #########@ add Archaic panel to plot
          x_pos_min=pmax(max(x_pos_annot)+1,1)
          x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
          x_pos_labels=c(x_pos_labels,c('Origin of Derived Allele','Origin of introgressed haplotype (CEU)','Origin of introgressed haplotype (CHS)'))

          DF_ARCHAIC[,x_pos:=x_pos+x_pos_min]
          p <- p + geom_tile(data=DF_ARCHAIC,aes(x_pos,y_pos),fill='white',col='black')
          p <- p + new_scale("fill")
          p <- p + scale_fill_manual(values=color_archaics)
          for (i in which(DF_ARCHAIC$value!='')){
            p <- p + geom_tile(data=DF_ARCHAIC[i,],aes(x_pos,y_pos,fill=value))
          }
            p <- p + geom_point(data=DF_ARCHAIC[value!='',],aes(x_pos,y_pos,pch=sign),size=mySize,fill='black',alpha=.5)
            p <- p + scale_shape_manual(values=c('Increased'=24,'Decreased'=25))
        }
    if( addFreqs ){
      #########@ assemble Frequency panel
      DF_FREQ=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,]
      DF_FREQ=unique(DF_FREQ[,.(gene_name,snps,IS_DERIVED=ifelse(INTROGRESSED_DERANC=='DER',1,0),CEU=FreqARCHAIC_CEU,CHS=FreqARCHAIC_CHS,YRI=FreqARCHAIC_YRI)])
      DF_FREQ=melt(DF_FREQ,id.vars=c('gene_name','snps'))
      DF_FREQ[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
      DF_FREQ[,x_pos:=match(variable,c('IS_DERIVED','CEU','CHS','YRI'))]
      #########@ add immune group panel to plot
      x_pos_min=pmax(max(x_pos_annot)+1,1)
      x_pos_annot=c(x_pos_annot,x_pos_min+1:4)
      x_pos_labels=c(x_pos_labels,c('Introgressed Allele = Derived','Frequency of introgressed allele (CEU)','Frequency of introgressed allele (CHS)','Frequency of introgressed allele (YRI)'))

      DF_FREQ[,x_pos:=x_pos+x_pos_min]
      p <- p + geom_tile(data=DF_FREQ,aes(x_pos,y_pos),fill='white',col='black')
      p <- p + new_scale("fill")
      p <- p + scale_fill_manual(values=setNames(c('black',color_populations),c('IS_DERIVED','YRI','CEU','CHS')))
      p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0,1),limits=c(0,1))
      for (i in which(DF_FREQ$value!='')){
        p <- p + geom_tile(data=DF_FREQ[i,],aes(x_pos,y_pos,fill=variable,alpha=value))
      }
    }
    if( addAdaptive ){
      ADAPT_GRID=data.table(expand.grid(x_pos=1:2,y_pos=1:length(eQTL_set)))
      #########@ assemble Adaptive panel
      DF_ADAPT=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,]
      DF_ADAPT=unique(DF_ADAPT[,.(gene_name,snps,Population, Adaptive=(FreqARCHAIC>q95)+(FreqARCHAIC>q99))])
      DF_ADAPT=melt(DF_ADAPT,id.vars=c('gene_name','snps','Population'))
      DF_ADAPT[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
      DF_ADAPT[,x_pos:=match(Population,c('CEU','CHS'))]
      #########@ add Adaptive panel to plot
      x_pos_min=pmax(max(x_pos_annot)+1,1)
      x_pos_annot=c(x_pos_annot,x_pos_min+1:2)
      x_pos_labels=c(x_pos_labels,c('Adaptive in CEU','Adaptive in CHS'))

      DF_ADAPT[,x_pos:=x_pos+x_pos_min]
      ADAPT_GRID[,x_pos:=x_pos+x_pos_min]
      p <- p + geom_tile(data=ADAPT_GRID,aes(x_pos,y_pos),fill='white',col='grey')
      p <- p + geom_tile(data=DF_ADAPT,aes(x_pos,y_pos),fill='white',col='black')
      p <- p + new_scale("fill")
      p <- p + scale_fill_manual(values=setNames(color_populations[2:3],c('CEU','CHS')))
      p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0,1),limits=c(0,2))
      for (i in which(DF_ADAPT$value!='')){
        p <- p + geom_tile(data=DF_ADAPT[i,],aes(x_pos,y_pos,fill=Population,alpha=value))
      }
    }
    if( addLength ){
      NUM_GRID=data.table(expand.grid(x_pos=c(1:3,5:7,9:10),y_pos=1:length(eQTL_set)))
      #########@ assemble haplotype length panel
      DF_NUM=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,]
      DF_NUM=unique(DF_NUM[,.(gene_name,snps,Population,NUM_aSNPs_NEAND, NUM_aSNPs_DENI, NUM_aSNPs_ARCHAIC,Length_kb)])
      #DF_NUM=dcast(DF_NUM,gene_name+snps~Population,fill=0,value.var=c('NUM_aSNPs_NEAND','NUM_aSNPs_DENI','NUM_aSNPs_ARCHAIC','Length_kb'))
      DF_NUM=melt(DF_NUM,id.vars=c('gene_name','snps','Population'))
      DF_NUM[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]

      DF_LENGTH=DF_NUM[variable=='Length_kb',]
      DF_NUM=DF_NUM[variable!='Length_kb',]
      DF_NUM[,x_pos:=match(variable,c('NUM_aSNPs_NEAND','NUM_aSNPs_DENI','NUM_aSNPs_ARCHAIC'))+ifelse(Population=='CHS',4,0)]
      DF_LENGTH[,x_pos:=match(Population,c('CEU','CHS'))+8]
      #########@ add haplotype length panel to plot
      x_pos_min=pmax(max(x_pos_annot)+1,1)
      x_pos_annot=c(x_pos_annot,x_pos_min+c(1:3,5:7,9:10))
      x_pos_labels=c(x_pos_labels,'Nb of linked Vindija aSNPs (CEU)','Nb of linked Denisova aSNPs (CEU)','Nb of linked ancient aSNPs (CEU)','Nb of linked Vindija aSNPs (CHS)','Nb of linked Denisova aSNPs (CHS)','Nb of linked ancient aSNPs (CHS)','length of archaic haplotype (CEU)','length of archaic haplotype (CHS)')

      DF_NUM[,x_pos:=x_pos+x_pos_min]
      DF_LENGTH[,x_pos:=x_pos+x_pos_min]
      NUM_GRID[,x_pos:=x_pos+x_pos_min]
      p <- p + geom_tile(data=NUM_GRID[x_pos<x_pos_min+4,],aes(x_pos,y_pos),fill='white',col='grey')
      p <- p + geom_tile(data=NUM_GRID[x_pos>x_pos_min+4 & x_pos<x_pos_min+8,],aes(x_pos,y_pos),fill='white',col='grey')
      p <- p + geom_tile(data=NUM_GRID[x_pos>x_pos_min+8,],aes(x_pos,y_pos),fill='white',col='grey')

      p <- p + geom_tile(data=DF_NUM,aes(x_pos,y_pos),fill='white',col='black')
      p <- p + geom_tile(data=DF_LENGTH,aes(x_pos,y_pos),fill='white',col='black')
      p <- p + new_scale("fill")
      p <- p + scale_fill_manual(values=setNames(c(color_archaics[2:4]),c('NUM_aSNPs_DENI','NUM_aSNPs_NEAND','NUM_aSNPs_ARCHAIC')))
      p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0,1),limits=c(0,2.6))
      for (i in which(DF_NUM$value!='')){
        p <- p + geom_tile(data=DF_NUM[i,],aes(x_pos,y_pos,fill=variable,alpha=log10(1+value)))
      }
      p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0,1),limits=c(0,4))
      for (i in which(DF_LENGTH$value!='')){
        p <- p + geom_tile(data=DF_LENGTH[i,],aes(x_pos,y_pos,alpha=log10(1+value)),fill='black')
      }
    }
    if( addCOVID ){
    #########@ assemble COVID panel
        DF_COVID=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,][!duplicated(paste(snps,gene_name,Population,type)),]
        DF_COVID=DF_COVID[,.(snps,gene_name,covid_A2_beta,covid_A2_pval,covid_B2_beta,covid_B2_pval,covid_C2_beta,covid_C2_pval,which_INTROGRESSED)]
        DF_COVID=melt(DF_COVID,id.vars=c('snps','gene_name','which_INTROGRESSED'))
        DF_COVID[,variable:=gsub('covid_','',variable)]
        DF_COVID[,stat:=gsub('([A-C]2)_(pval|beta)','\\2',variable)]
        DF_COVID[,trait:=gsub('([A-C]2)_(pval|beta)','\\1',variable)]
        DF_COVID=dcast(unique(DF_COVID),gene_name+snps+trait+which_INTROGRESSED~stat)
        DF_COVID[,x_pos:=match(trait,c('C2','B2','A2'))]
        DF_COVID[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
        DF_COVID[,sign:=ifelse((beta*ifelse(which_INTROGRESSED=='ALT',1,-1))>0,'Increased','Decreased')]

    #########@ add COVID panel to plot
        x_pos_min=pmax(max(x_pos_annot)+1,1)
        x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
        x_pos_labels=c(x_pos_labels,c('Reported','Hospitalized','Critical'))

        DF_COVID[,x_pos:=x_pos+x_pos_min]
        p <- p + geom_tile(data=DF_COVID,aes(x_pos,y_pos),fill='white',col='black')
        p <- p + new_scale("alpha")+ scale_alpha_continuous(range=c(0.2,1))
        p <- p + geom_point(data=DF_COVID[pval<0.01,],aes(x_pos,y_pos,pch=sign,alpha=log10(-log10(pval))),fill='black',size=mySize)
      }

      if( addCOLOC ){
      #########@ assemble colocalization panel
          DF_COLOC=gene_list[paste(gene_name,snps,sep=' - ')%in%eQTL_set,]
          DF_COLOC=DF_COLOC[,.(   tested_coloc=ifelse(any(tested_coloc=='yes'),1,0),
                                  colocalized=ifelse(any(colocalized=='yes'),1,0),
                                  PP.H4.abf_reported=max(c(PP.H4.abf_reported,0),na.rm=T),
                                  PP.H4.abf_hospitalized=max(c(PP.H4.abf_hospitalized,0),na.rm=T),
                                  PP.H4.abf_critical=max(c(PP.H4.abf_critical,0),na.rm=T)),by=.(gene_name,snps)]

          DF_COLOC=melt(DF_COLOC,id.vars=c('gene_name','snps'))
          DF_COLOC[,y_pos:=match(paste(gene_name,snps,sep=' - '),eQTL_set)]
          #DF_COLOC[,x_pos:=match(variable,c('tested_coloc','colocalized','PP.H4.abf_reported','PP.H4.abf_hospitalized','PP.H4.abf_critical'))]
          DF_COLOC[,x_pos:=match(variable,c('PP.H4.abf_reported','PP.H4.abf_hospitalized','PP.H4.abf_critical'))]
          # IMMUNE_GRID=data.table(expand.grid(x_pos=1:length(lineage_order),y_pos=1:length(eQTL_set)))

      #########@ add colocalization panel to plot
          x_pos_min=pmax(max(x_pos_annot)+1,1)
          # x_pos_annot=c(x_pos_annot,x_pos_min+1:5)
          # x_pos_labels=c(x_pos_labels,c('colocalization tested','signif. colocalization','max PPH4 (reported)','max PPH4 (hospitalized)','max PPH4 (critical)'))
          x_pos_annot=c(x_pos_annot,x_pos_min+1:3)
          x_pos_labels=c(x_pos_labels,c('max PPH4 (reported)','max PPH4 (hospitalized)','max PPH4 (critical)'))

          DF_COLOC[,x_pos:=x_pos+x_pos_min]
          p <- p + geom_tile(data=DF_COLOC,aes(x_pos,y_pos),fill='white',col='black')
           #for (i in which(     )){
            p <- p + new_scale("alpha") + scale_alpha_continuous(limits=c(0,1),range=c(0,1))
            p <- p + geom_tile(data=DF_COLOC,aes(x_pos,y_pos,alpha=value),fill='black')

           #}
      }
    p <- p + scale_x_continuous(breaks=x_pos_annot, labels=x_pos_labels)
    p <- p + coord_fixed()
    if(legend){
    p <- p + theme(legend.spacing.x= unit(1, 'mm'),
                    legend.spacing.y=unit(2, 'mm'),
                    legend.text=element_text(size=6),
                    legend.key = element_blank(),legend.box = "vertical",
                    legend.margin=margin(),legend.key.size = unit(2, 'mm'))
              }else{
    p <- p + theme(legend.position='none')
              }
    return(p)
}


p <- plot_geneList_ARCHAIC(gene_list_adaptive_introgressed_immune,
                        addLineage=FALSE,
                        addLineageState=TRUE,
                        addLineageStateSummary=FALSE,
                        addImmuneAnnot=TRUE,
                        addLOGFC=FALSE,
                        addCOVID=TRUE,
                        addCOLOC=TRUE,
                        addARCHAIC=FALSE,
                        addFreqs=TRUE,
                        addLength=FALSE,
                        addAdaptive=TRUE,mySize=1.5,legend=FALSE)
pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_adaptive_introgressed_immune_v3.pdf",EIP))
print(p)
dev.off()



  p <- plot_geneList_ARCHAIC(gene_list_introgressed_coloc,
                          addLineage=FALSE,
                          addLineageState=TRUE,
                          addLineageStateSummary=FALSE,
                          addImmuneAnnot=TRUE,
                          addLOGFC=TRUE,
                          addCOVID=TRUE,
                          addCOLOC=TRUE,
                          addARCHAIC=TRUE,
                          addFreqs=TRUE,
                          addLength=TRUE,
                        addAdaptive=TRUE,mySize=1.5,legend=FALSE)
  pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_introgressed_coloc_v8_effectIntrogressed.pdf",EIP))
  print(p)
  dev.off()

p <- plot_geneList_ARCHAIC(gene_list_immuneGO_introgressed,
                        addLineage=FALSE,
                        addLineageState=TRUE,
                        addLineageStateSummary=FALSE,
                        addImmuneAnnot=TRUE,
                        addLOGFC=TRUE,
                        addCOVID=TRUE,
                        addCOLOC=TRUE,
                        addARCHAIC=TRUE,
                        addFreqs=TRUE,
                        addLength=TRUE,
                        addAdaptive=TRUE,mySize=1.5,legend=FALSE)
pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_immuneGO_introgressed_v2.pdf",EIP))
print(p)
dev.off()


  p <- plot_geneList_ARCHAIC(gene_list_adaptive_introgressed_immune,
                          addLineage=FALSE,
                          addLineageState=TRUE,
                          addLineageStateSummary=FALSE,
                          addImmuneAnnot=TRUE,
                          addLOGFC=TRUE,
                          addCOVID=TRUE,
                          addCOLOC=TRUE,
                          addARCHAIC=TRUE,
                          addFreqs=TRUE,
                          addLength=TRUE,
                          addAdaptive=TRUE,mySize=1.5,legend=FALSE)
  pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_adaptive_introgressed_immune_v2.pdf",EIP))
  print(p)
  dev.off()


p <- plot_geneList_ARCHAIC(gene_list_adaptive_introgressed,
                        addLineage=FALSE,
                        addLineageState=TRUE,
                        addLineageStateSummary=FALSE,
                        addImmuneAnnot=TRUE,
                        addLOGFC=TRUE,
                        addCOVID=TRUE,
                        addCOLOC=TRUE,
                        addARCHAIC=TRUE,
                        addFreqs=TRUE,
                        addLength=TRUE,
                        addAdaptive=TRUE,mySize=1.5,legend=FALSE)
pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_adaptive_introgressed_v2.pdf",EIP))
print(p)
dev.off()


p <- plot_geneList_ARCHAIC(gene_list_archaic_immune,
                        addLineage=FALSE,
                        addLineageState=TRUE,
                        addLineageStateSummary=FALSE,
                        addImmuneAnnot=TRUE,
                        addLOGFC=TRUE,
                        addCOVID=TRUE,
                        addCOLOC=TRUE,
                        addARCHAIC=TRUE,
                        addFreqs=TRUE,
                        addLength=TRUE,
                        addAdaptive=TRUE,mySize=1.5,legend=FALSE)
pdf(sprintf("%s/single_cell/project/pop_eQTL/paper_draft/V7/testsFigure/Fig6/geneLists/00_gene_list_archaic_immune_v2.pdf",EIP))
print(p)
dev.off()
