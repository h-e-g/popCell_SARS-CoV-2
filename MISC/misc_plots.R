################################################################################
# Useful functions

# custom plot theme
theme_plot=function(lpos="bottom",rotate.x=0) {
  ptheme <- theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=7,family="sans"),
          legend.position=lpos,
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill="#012158"),
          panel.spacing = unit(0,"pt"),
          strip.text = element_text(color="white"))
        if(rotate.x>0){
          ptheme <- ptheme+ theme(axis.text.x=element_text(angle=rotate.x, hjust=1,vjust=0.5))
        }
  ptheme
}

# function to add names to a vector on the fly.
add_names=function(vector,name){
    names(vector)=name;return(vector)
  }

# function to shift legend to empty panel

shift_legend <- function(p){

  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }

  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }

  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")

  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")

  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")

  return(gp)
}

################################################################################
# Color definitions

# condition colors, order
condition_color=c("NS"="#888888","IAV"="#6699CC","COV"="#DD5555")
condition_order=c("NS","COV","IAV")

# condition and time point colors, order (time course)
cond_tp_color <- c("T0_T0"="#888888","NS_T6"="#888888","NS_T24"="#888888","IAV_T6"="#6699CC","IAV_T24"="#254d74","COV_T6"="#DD5555","COV_T24"="#8f0a0a")
cond_tp_order <- c("T0_T0","NS_T6","NS_T24","COV_T6","COV_T24","IAV_T6","IAV_T24")

# lineage colors, order
lineage_order=c("MONO","B","T.CD4","T.CD8","NK")
lineage_color <- c("B"="#9281DE","MONO"="#B7637E", "T.CD4"="#005274", "T.CD8"="#048c41", "NK"= "#7D9E00")

# celltype colors, order
celltype_color=c(
  "MONO.CD14"="#B7637E","MONO.CD16"="#8F568A","pDC"="#8E6B00","cDC"="#E7b315","MONO.CD14.INFECTED"="#B9364B",
  "B.N.K"="#9281DE","B.N.L"="#a681de","B.M.K"="#6045d9","B.M.L"="#8045d9","Plasmablast"="#6801A4",
  "T.CD4.N"="#169FD8","T.CD4.E"="#005274","T.Reg"="#03C2C0","T.gd"="#ffce73",
  "T.CD8.N"="#00BB54","T.CD8.CM.EM"="#0EAD20","T.CD8.EMRA"="#048c41","MAIT"="#004A21",
  "NK.M.LIKE"="#B7DC2A","NK.CD56dim"="#7D9E00","NK.CD56brt"="#3F4F00","ILC"="#63C425"
)
celltype_order=c(
  "MONO.CD14","MONO.CD16","MONO.CD14.INFECTED","cDC","pDC",
  "B.N.K","B.N.L","B.M.K","B.M.L","Plasmablast",
  "T.CD4.N","T.CD4.E","T.Reg","T.gd",
  "T.CD8.N","T.CD8.CM.EM","T.CD8.EMRA","ILC","MAIT",
  "NK.CD56dim","NK.CD56brt","NK.M.LIKE"
)

# population colors, order
population_color=c("AFB"="#008000","EUB"="#eba206","ASH"="#71458d")
population_order=c("AFB","EUB","ASH")

population_1kg_color=c("YRI"="#008000","CEU"="#eba206","CHS"="#71458d")
population_1kg_order=c("YRI","CEU","CHS")

# COVID-19 susceptibility and severity colors, order
COVID_color=c("critical"="#DD5555","hospitalized"="#E68686","reported"="#F1BABA")
COVID_order=c("reported","hospitalized","critical")

# archaic colors, order
archaic_color=c('AMH'=grey(0.8),'DENI'='#80CDC1','NEAND'='#FDAE61','ARCHAIC'="#FEE08B",'UNKNOWN_ARCHAIC'=grey(0.5))

# add 12 entry divergent colorblind_palette # (taken from the rcartocolor package: safe palette)
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
