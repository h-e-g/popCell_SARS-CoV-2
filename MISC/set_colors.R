suppressMessages(library(tictoc))

#Define custom ggplot2 theme
theme_mary <- function(lpos="bottom") {
  font <- "Arial"
  theme_light() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=14),
          #axis.text.x=element_text(angle=45, hjust=1),
          legend.position=lpos,
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA))
}


theme_yann=function(lpos="bottom",rotate.x=0) {
  ptheme <- theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=8,family="Helvetica"),
          legend.position=lpos,
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background=element_rect(fill="#012158"),
          strip.text=element_text(color="white"))
          if(rotate.x>0){
            ptheme <- ptheme+ theme(axis.text.x=element_text(angle=rotate.x, hjust=1,vjust=0.5))
          }
    ptheme
}

mergeCols=function(col1,col2,pct2=0.5){colorRampPalette(c(col1,col2))(100)[round(pct2*100)]}

# function to add names to a vector on the fly.
add_names=function(vector,name){
    names(vector)=name;return(vector)
  }

# add 12 entry divergent colorblind_palette # (taken from rcartocolor package: safe palette)
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

# Maxime"s condition colors
color_conditions=c("NS"="#888888","IAV"="#6699CC","COV"="#DD5555")
cond_color=color_conditions
# define activation colors
color_activation <- c("NS"="#888888","COV"="#DD5555","IAV"="#6699CC","ACTIVE"="#B7520D","SHARED"="lightgrey","RESTING"="gray")

# define activation colors
# color_activation <- c("NS"="black","COV"="#c2a016","IAV"="#ac0505","ACTIVE"="#B7520D","SHARED"="lightgrey","RESTING"="gray")
# cond_colors <- c("COV"="#c2a016","IAV"="#ac0505","NS"="black")
# color_condition <- cond_colors

# define broad cell types colors
#color_cellTypes_5level <- c("Tcells"="#169FD8","NK"="#B7DC2A","DC"="#E7b315","Mono"="#B7637E","Bcells"="#2D04E5","unassigned"="lightgrey")
color_cellTypes_8level <- c("T.CD4"="#005274", "T.CD8"="#048c41", "NK"= "#7D9E00", "DC"="#E7b315", "MONO"="#B7637E", "B"="#9281DE", "Plasmablast"="#6801A4", "Progenitors"="#8E500A")
color_cellTypes_6level <- color_cellTypes_8level[c("MONO","DC","B","T.CD4","T.CD8","NK")]
celltype.6level_colors=color_cellTypes_6level

#### define celltype annotations
color_cellTypes_15level <- c("T.CD4.N"="#169FD8", "T.CD4.E"="#005274", "T.CD8.N"="#00BB54", "T.CD8.CM.EM"="#048c41", "T.CD8.TE"="#19cf00","T.Reg"="#03C2C0","T.gd"="#ffce73", "MAIT"="#004A21", "NK.CD56dim"= "#B7DC2A", "NK.CD56brt"= "#7D9E00","B.N"="#9281DE", "B.M"="#2D04E5", "Plasmablast"="#6801A4","Myeloid.CD14"="#B7637E", "Myeloid.CD16"="#8F568A", "pDC"="#8E6B00", "cDC"="#E7b315", "mDC"="#e78815","Progenitors"="#8E500A")# "TDB"="black", "Low"="#7e7d7d")
color_cellTypes_19level <- color_cellTypes_15level
color_cellTypes_20level <- c(color_cellTypes_19level,"INF.MONO"="#B9364B")

color_cellid=c("B intermediate"="#5E41E1", "B memory"="#2D04E5", "B naive"="#9281DE", "CD14 Mono"="#B7637E",
"CD16 Mono"="#8F568A",  "CD4 Naive"="#169FD8", "CD4 TCM"="#0B78A6", "CD4 TEM"="#005274", "Treg"="#03C2C0", "CD4 CTL"="#016E5A","dnT"="#83D5C9","CD8 Naive"="#00BB54",
"CD8 TCM"="#01A34A", "CD8 TEM"="#048c41", "ILC"="#3F4F00", "MAIT"="#004A21", "gdT"="#ffce73", "NK"="#B7DC2A",   "cDC1"="#E79D15", "cDC2"="#E7b315",
 "pDC"="#8E6B00","ASDC"="#644B00", "unassigned"="#7e7d7d")


color_populations=c("AFB"="#008000","EUB"="#eba206","ASH"="#71458d")

cond_tp_colors <- c("T0_T0"="#888888","NS_T6"="#888888","NS_T24"="#888888","IAV_T6"="#6699CC","IAV_T24"="#254d74","COV_T6"="#DD5555","COV_T24"="#8f0a0a")
cond_tp_order <- c("T0_T0","NS_T6","NS_T24","COV_T6","COV_T24","IAV_T6","IAV_T24")

color_cellTypes_morelevels <- c(color_cellTypes_19level,"T.CD8.TEMRA"="#19cf00","NK.M.LIKE"=mergeCols("#B7DC2A","#19cf00"), "Myeloid.CD14.INFECTED"="#B9364B","B.INFECTED"=mergeCols("#B9364B","#9281DE"),"NK.INFECTED"=mergeCols("#B7DC2A","#B9364B"),"T.N.INFECTED"=mergeCols(mergeCols("#048c41","#005274"),"#B9364B"),
"Myeloid.CD14.other"=mergeCols("#B7637E","#EEEEEE",0.4),"NK.other"=mergeCols("#7D9E00","#EEEEEE",0.4),"B.other"=mergeCols("#9281DE","#EEEEEE",0.4),"T.other"=mergeCols(mergeCols("#048c41","#005274"),"#EEEEEE",0.4),"MIX_NK_CD8TE"=mergeCols(mergeCols("#7D9E00","#19cf00"),"#EEEEEE",0.4),"ILC?"="#3F4F00","LowQ"="#7e7d7d")
color_cellTypes_YA=color_cellTypes_morelevels
color_cellTypes_YA = setNames(color_cellTypes_YA,gsub("T.CD8.TEMRA","T.CD8.EMRA",names(color_cellTypes_YA)))

whiteningFactor=0.6

color_cellTypes_YA=c("T.CD4.N"="#169FD8", "T.CD4.E"="#005274", "T.Reg"="#03C2C0","T.gd"="#ffce73", "MAIT"="#004A21",'ILC?'="#3F4F00", "T.other"=mergeCols(mergeCols('#048c41','#005274'),'#EEEEEE',0.4),
"MIX_T.CD4_CD8"=mergeCols(mergeCols("#169FD8","#00BB54"),"#EEEEEE",whiteningFactor),
"T.CD8.N"="#00BB54", "T.CD8.CM.EM"="#048c41", "T.CD8.EMRA"="#0EAD20","MIX_NK_CD8TE"=mergeCols(mergeCols('#B7DC2A',"#0EAD20",.75),'#EEEEEE',whiteningFactor),
'NK.M.LIKE'=mergeCols('#B7DC2A',"#0EAD20"), "NK.CD56dim"= "#B7DC2A", "NK.CD56brt"= "#7D9E00",'NK.other'=mergeCols('#7D9E00','#EEEEEE',0.4),
"B.N"="#9281DE", "B.M"="#2D04E5", "Plasmablast"="#6801A4",'B.other'=mergeCols('#9281DE','#EEEEEE',0.4),
"Myeloid.CD14"="#B7637E", "Myeloid.CD16"="#8F568A", "pDC"="#8E6B00", "cDC"="#E7b315", "mDC"="#e78815",'Myeloid.CD14.other'=mergeCols('#B7637E','#EEEEEE',0.4),
"Progenitors"="#8E500A",'LowQC'='#7e7d7d')

color_cellTypes_AB <- c("T.CD4.N"="#169FD8", "T.CD4.E"="#005274", "T.Reg"="#03C2C0","T.gd"="#ffce73", "MAIT"="#004A21",'ILC'="#3F4F00",
"MIX_T.CD4_CD8"=mergeCols(mergeCols("#169FD8","#00BB54"),"#EEEEEE",whiteningFactor),
"MIX_T.CD4.N_T.CD4.E"=mergeCols(mergeCols("#169FD8","#005274"),"#EEEEEE",whiteningFactor),
"MIX_T.CD8.CM.EM_MAIT"=mergeCols(mergeCols("#048c41","#004A21"),"#EEEEEE",whiteningFactor),
"T.CD8.N"="#00BB54", "T.CD8.CM.EM"="#048c41", "T.CD8.EMRA"="#0EAD20",
"MIX_CD8.CM.EM_NK.M.LIKE"=mergeCols(mergeCols("#048c41",mergeCols("#B7DC2A","#0EAD20")),"#EEEEEE",whiteningFactor),
"MIX_CD8.EMRA_NK.M.LIKE"=mergeCols(mergeCols("#0EAD20",mergeCols("#B7DC2A","#0EAD20")),"#EEEEEE",whiteningFactor),
"MIX_T.CD8.CM.EM_NK"=mergeCols(mergeCols("#048c41","#B7DC2A"),"#EEEEEE",whiteningFactor),
'NK.M.LIKE'=mergeCols('#B7DC2A',"#0EAD20"), "NK.CD56dim"= "#B7DC2A", "NK.CD56brt"= "#7D9E00",'NK.other'=mergeCols('#7D9E00','#EEEEEE',0.4),
"B.N"="#9281DE", "B.M"="#2D04E5", "Plasmablast"="#6801A4",'B.other'=mergeCols('#9281DE','#EEEEEE',0.4),
"Myeloid.CD14"="#B7637E", "Myeloid.CD16"="#8F568A", "pDC"="#8E6B00", "cDC"="#E7b315", "mDC"="#e78815",'Myeloid.CD14.other'=mergeCols('#B7637E','#EEEEEE',0.4),
"Progenitors"="#8E500A",'LowQC'='#7e7d7d')


color_cellTypes_21level <- c("T.CD4.N"="#169FD8", "T.CD4.E"="#005274", "T.Reg"="#03C2C0",
"T.gd"="#ffce73","T.CD8.N"="#00BB54", "T.CD8.CM.EM"="#0EAD20", "T.CD8.EMRA"="#048c41","MAIT"="#004A21",
'NK.M.LIKE'=mergeCols('#B7DC2A',"#0EAD20"), "NK.CD56dim"= "#B7DC2A", "NK.CD56brt"= "#7D9E00",'ILC'="#3F4F00",
"B.N"="#9281DE", "B.M"="#2D04E5", "Plasmablast"="#6801A4",
"MONO.CD14"="#B7637E", "MONO.CD16"="#8F568A", "pDC"="#8E6B00", "cDC"="#E7b315", "mDC"="#e78815",
"Progenitors"="#8E500A","MONO.CD14.INFECTED"="#B9364B")

color_cellTypes_24level <- c("T.CD4.N"="#169FD8", "T.CD4.E"="#005274", "T.Reg"="#03C2C0",
"T.gd"="#ffce73","T.CD8.N"="#00BB54", "T.CD8.CM.EM"="#0EAD20", "T.CD8.EMRA"="#048c41","MAIT"="#004A21",
'NK.M.LIKE'=mergeCols('#B7DC2A',"#0EAD20"), "NK.CD56dim"= "#B7DC2A", "NK.CD56brt"= "#7D9E00",'ILC'="#3F4F00",
"B.N.K"="#9281DE", "B.N.L"="#a681de", "B.M.K"="#6045d9", "B.M.L"="#8045d9", "B.M.CCL22"="#2D04E5", "Plasmablast"="#6801A4",
"MONO.CD14"="#B7637E", "MONO.CD16"="#8F568A", "pDC"="#8E6B00", "cDC"="#E7b315", "mDC"="#e78815",
"Progenitor"="#8E500A","MONO.CD14.INFECTED"="#B9364B")

celltype_color=c(
"T.CD4.N"="#169FD8","T.CD4.E"="#005274","T.Reg"="#03C2C0","T.gd"="#ffce73",
"T.CD8.N"="#00BB54","T.CD8.CM.EM"="#0EAD20","T.CD8.EMRA"="#048c41","MAIT"="#004A21",
"NK.M.LIKE"="#B7DC2A","NK.CD56dim"="#7D9E00","NK.CD56brt"="#3F4F00","ILC"="#63C425",
"B.N.K"="#9281DE","B.N.L"="#a681de","B.M.K"="#6045d9","B.M.L"="#8045d9",
"Plasmablast"="#6801A4","MONO.CD14"="#B7637E","MONO.CD16"="#8F568A","pDC"="#8E6B00",
"cDC"="#E7b315","MONO.CD14.INFECTED"="#B9364B")

celltype_order=c("MONO.CD14","MONO.CD16","MONO.CD14.INFECTED","cDC","pDC","B.N.K","B.N.L","B.M.K","B.M.L","Plasmablast","T.CD4.N","T.CD4.E","T.Reg","T.gd","T.CD8.N","T.CD8.CM.EM","T.CD8.EMRA","ILC","MAIT","NK.CD56dim","NK.CD56brt","NK.M.LIKE")
celltype_label=c("MONO CD14+","MONO CD16+","MONO IAV+","cDC","pDC","B N k","B N l","B M k","B M l","Plasmablast","T CD4+ N","T CD4+ E","T Reg","T gd","T CD8+ N","T CD8+ CM/EM","T CD8+ EMRA","ILC","MAIT","NK CD56dim","NK CD56brt","NK mem")

lineage_order=c("MONO","B","T.CD4","T.CD8","NK")
lineage_color <- c("B"="#9281DE","MONO"="#B7637E", "T.CD4"="#005274", "T.CD8"="#048c41", "NK"= "#7D9E00")
lineage_label=c("MONO","B","CD4+ T","CD8+ T", "NK")
attr(lineage_label,'names')=c("MONO","B","T.CD4","T.CD8","NK")

color_cellTypes_6level <- color_cellTypes_8level[c("MONO","DC","B","T.CD4","T.CD8","NK")]

cond_order=c("NS","COV","IAV")

color_COVID=c(critical="#DD5555", hospitalized="#E68686", reported="#F1BABA")
color_archaics=c('AMH'=grey(0.8),'DENI'='#80CDC1','NEAND'='#FDAE61','ARCHAIC'="#FEE08B",'UNKNOWN_ARCHAIC'=grey(0.5))

# meta_raw[cluster_seurat==24,updated_celltype:='MAIT']
# meta_raw[cluster_seurat%in%c(40,60),updated_celltype:='Myeloid.CD14.INFECTED']
# meta_raw[cluster_seurat%in%c(83,88),updated_celltype:='B.INFECTED']
# meta_raw[cluster_seurat==4,updated_celltype:='T.CD8.TEMRA']
# meta_raw[cluster_seurat==15,updated_celltype:='T.CD8.TEMRA']
# meta_raw[cluster_seurat==9,updated_celltype:='NK.M.LIKE'] # KLRC2+ KLRC1- NK
# meta_raw[cluster_seurat==26,updated_celltype:='NK.M.LIKE']
# meta_raw[cluster_seurat==69,updated_celltype:='NK.INFECTED']
# meta_raw[cluster_seurat==64,updated_celltype:='T.N.INFECTED']
# meta_raw[cluster_seurat%in%c(59,48,57),updated_celltype:='MIX_NK_CD8TE']
# meta_raw[cluster_seurat%in%c(90,52),updated_celltype:='NK.other']
# meta_raw[cluster_seurat%in%c(42,45,47,48,80,51,49,68,67),updated_celltype:='T.other']
# meta_raw[cluster_seurat%in%c(72,75,86),updated_celltype:='B.other']
# meta_raw[cluster_seurat%in%c(74),updated_celltype:='Myeloid.CD14.other']
# meta_raw[cluster_seurat%in%c(46),updated_celltype:='T.CD4.E']
# meta_raw[cluster_seurat%in%c(82),updated_celltype:='ILC?']
# meta_raw[cluster_seurat%in%c(47,53,50,61,65,70,71,77,87),updated_celltype:='LowQ']
