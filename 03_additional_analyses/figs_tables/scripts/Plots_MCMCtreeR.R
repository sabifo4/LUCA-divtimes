#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#--------------#
# LOAD PACKAGE #
#--------------#
library( MCMCtreeR )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
main_wd            <- gsub( pattern = "figs_tables/scripts/", replacement = "", x = wd )
wd_conc_cb_GBM     <- gsub( pattern = "figs_tables/scripts/", replacement = "../01_PAML/02_MCMCtree/sum_analyses/01_posterior/mcmc_files_GBM", x = wd )
wd_conc_cb_ILN     <- gsub( pattern = "figs_tables/scripts/", replacement = "../01_PAML/02_MCMCtree/sum_analysessum_analyses/01_posterior/mcmc_files_ILN", x = wd )
wd_indATP_GBM      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_ATP_GBM", x = wd )
wd_indATP_ILN      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_ATP_GBM", x = wd )
wd_indEF_GBM       <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_EF_GBM", x = wd )
wd_indEF_ILN       <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_EF_GBM", x = wd )
wd_indLeu_GBM      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_Leu_GBM", x = wd )
wd_indLeu_ILN      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_Leu_GBM", x = wd )
wd_indSRP_GBM      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_SRP_GBM", x = wd )
wd_indSRP_ILN      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_SRP_GBM", x = wd )
wd_indTyr_GBM      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_Tyr_GBM", x = wd )
wd_indTyr_ILN      <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_analyses/01_posterior/mcmc_files_Tyr_GBM", x = wd )
wd_noATP_GBM       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noatp/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_GBM", x = wd )
wd_noATP_ILN       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noatp/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_ILN", x = wd )
wd_noEF_GBM        <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noef/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_GBM", x = wd )
wd_noEF_ILN        <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noef/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_ILN", x = wd )
wd_noLeu_GBM       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noleu/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_GBM", x = wd )
wd_noLeu_ILN       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noleu/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_ILN", x = wd )
wd_noSRP_GBM       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_nosrp/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_GBM", x = wd )
wd_noSRP_ILN       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_nosrp/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_ILN", x = wd )
wd_noTyr_GBM       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_notyr/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_GBM", x = wd )
wd_noTyr_ILN       <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_notyr/01_PAML/01_MCMCtree/sum_analyses/01_posterior/mcmc_files_ILN", x = wd )
wd_bsinbv_GBM      <- gsub( pattern = "figs_tables/scripts/", replacement = "02_blunch_analysis/02_PAML/MCMCtree/sum_analyses/01_posterior/mcmc_files_GBM", x = wd )
wd_bsinbv_ILN      <- gsub( pattern = "figs_tables/scripts/", replacement = "02_blunch_analysis/02_PAML/MCMCtree/sum_analyses/01_posterior/mcmc_files_ILN", x = wd )
all_wds <- c( main_wd,
              wd_conc_cb_GBM, wd_conc_cb_ILN,
              wd_indATP_GBM, wd_indATP_ILN, wd_indEF_GBM, wd_indEF_ILN,
              wd_indLeu_GBM, wd_indLeu_ILN, wd_indSRP_GBM, wd_indSRP_ILN,
              wd_indTyr_GBM, wd_indTyr_ILN,
              wd_noATP_GBM, wd_noATP_ILN, wd_noEF_GBM, wd_noEF_ILN,
              wd_noLeu_GBM, wd_noLeu_ILN, wd_noSRP_GBM, wd_noSRP_ILN,
              wd_noTyr_GBM, wd_noTyr_ILN, 
              wd_bsinbv_GBM, wd_bsinbv_ILN )
)
name_dirs <- c( "main_wd",
                "conc_cb_GBM", "conc_cb_ILN", "indATP_GBM", "indATP_ILN",
                "indEF_GBM", "indEF_ILN", "indLeu_GBM", "indLeu_ILN",
                "indSRP_GBM", "indSRP_ILN", "indTyr_GBM", "indTyr_ILN",
                "noATP_GBM", "noATP_ILN", "noEF_GBM", "noEF_ILN", "noLeu_GBM",
                "noLeu_ILN", "noSRP_GBM", "noSRP_ILN", "noTyr_GBM", "noTyr_ILN",
                "bsinbv_GBM", "bsinbv_ILN" )
source( "../../../src/Functions_plots_MCMCtreeR.R" )
source( "../../../src/Functions_plots.R" )

#### LOAD DATA ----
#-------------#
# START TASKS #
#-------------#
# Create objects for each dataset that will be required for plotting
plot_obj <- create_plot_obj( abs_path = all_wds, name_entries = name_dirs,
                             tree_pattern = "95CI\\.tree" )
# Create dir to save RData object
if( ! dir.exists( paste( main_wd, "figs_tables/out_RData", sep = "" ) ) ){
  dir.create( paste( main_wd, "figs_tables/out_RData", sep = "" ) )
}
# Save object, then comment the next line and just load the RData file!
# save( plot_obj, file = paste( main_wd,
#                               "figs_tables/out_RData/plot_obj.RData",
#                               sep = "" ) )
load( file = paste( main_wd, "figs_tables/out_RData/plot_obj.RData",
                    sep = "" ) )

#### PLOTS ----
#----------------#
# START PLOTTING #
#----------------#

## [[ PLOT MAIN FIGURE (i.e., cross-bracing A)]] ----

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 0. PARTITIONED DATASETS, GBM vs ILN, cross-bracing #
#----------------------------------------------------#
# Saved in PDF: plots/LUCAdup_GBMnILNmodels_mainfig.pdf
nodes_2plot <- c( 248,368, # LUCA
                  269,457, # LECA
                  249,438, # LACA
                  304,369, # LBCA
                  312,338,402,468, # ARCHAEPLASTIDA
                  329,391          # TG-EUKARYA-MITO
)
##> GBM + PART + CB
ind_mat <- which( names( plot_obj$node_ages$GBM_part_cb ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$GBM_part_cb, 
                                      #xlim.scale = c(-800,5000),  # wide
                                      #xlim.scale = c(-500,6000),  # less wide
                                      xlim.scale = c(-5000,7000), # narrow
                                      #xlim.scale = c(-6000,7000),  # very narrow
                                      node.ages = plot_obj$node_ages$GBM_part_cb[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.2,
                                      time.correction = 1000, 
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col,
                                      density.border.col = "blue"
)

# Add now the other distribution with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange",
               "pink", "lightblue", "cyan", "brown", "purple", "lightgreen" )
##> ILN + PART + CB
transp.col2  <- adjustcolor( col = col.plot[ 11 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_part_cb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_part_cb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 11 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "GBM, partitioned",
                                "ILN, partitioned" ),
        col = c( "blue", col.plot[c(11)]),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( rep( "LUCA", 2 ), rep( "LECA", 2 ),
                           rep( "LACA", 2 ), rep( "LBCA", 2 ),
                           rep( "ARCHAEPLASTIDA", 4),
                           rep( "TG-EUKARYA-MITO", 2 ) )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x),
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 2,
        names( sort_nodes2plot )[i] )    
}


## [[ PLOT SUPP FIGURE - without cross-bracing ]] ----

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 1. Part VS Conc, GBM vs ILN, cross-bracing A #
#----------------------------------------------#
# Saved in PDF: plots/LUCAdup_cb_suppfig.pdf
nodes_2plot <- c( 248,368, # LUCA
                  269,457, # LECA
                  249,438, # LACA
                  304,369, # LBCA
                  312,338,402,468, # ARCHAEPLASTIDA
                  329,391          # TG-EUKARYA-MITO
)
##> GBM + PART + NOT-CB
ind_mat <- which( names( plot_obj$node_ages$GBM_part_cb ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$GBM_part_cb, 
                                      #xlim.scale = c(-800,5000),  # wide
                                      #xlim.scale = c(-500,6000),  # less wide
                                      xlim.scale = c(-5000,7000), # narrow
                                      #xlim.scale = c(-6000,7000),  # very narrow
                                      node.ages = plot_obj$node_ages$GBM_part_cb[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.2,
                                      time.correction = 1000, 
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col,
                                      density.border.col = "blue"
)

# Add now the other distribution with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange",
               "pink", "lightblue", "cyan", "brown", "purple", "lightgreen" )
##> ILN + PART + CB
transp.col2  <- adjustcolor( col = col.plot[ 11 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_part_cb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_part_cb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 11 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )
##> GBM + CONC + CB
coords_2plot <- add.extra.dists( phy = plot_obj$phy$GBM_conc_cb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$GBM_conc_cb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = "white",
                                 density.border.col = col.plot[ 3 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )
##> ILN + CONC + CB
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_conc_cb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_conc_cb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = "white",
                                 density.border.col = col.plot[ 10 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "GBM, part, cross-bracing A",
                                "ILN, part, cross-bracing A",
                                "GBM, conc, cross-bracing A",
                                "ILN, conc, cross-bracing A" ),
        col = c( "blue", col.plot[c(11, 3, 10)]),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( rep( "LUCA", 2 ), rep( "LECA", 2 ),
                           rep( "LACA", 2 ), rep( "LBCA", 2 ),
                           rep( "ARCHAEPLASTIDA", 4),
                           rep( "TG-EUKARYA-MITO", 2 ) )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x),
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 2,
        names( sort_nodes2plot )[i] )    
}


## [[ PLOT SUPP FIGURE - without cross-bracing ]] ----

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 2. Part VS Conc, GBM vs ILN, no cross-bracing #
#-----------------------------------------------#
# Saved in PDF: plots/LUCAdup_notcb_suppfig.pdf
nodes_2plot <- c( 248,368, # LUCA
                  269,457, # LECA
                  249,438, # LACA
                  304,369, # LBCA
                  312,338,402,468, # ARCHAEPLASTIDA
                  329,391          # TG-EUKARYA-MITO
)
##> GBM + PART + NOT-CB
ind_mat <- which( names( plot_obj$node_ages$GBM_part_notcb ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$GBM_part_notcb, 
                                      #xlim.scale = c(-800,5000),  # wide
                                      #xlim.scale = c(-500,6000),  # less wide
                                      xlim.scale = c(-5000,7000), # narrow
                                      #xlim.scale = c(-6000,7000),  # very narrow
                                      node.ages = plot_obj$node_ages$GBM_part_notcb[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.2,
                                      time.correction = 1000, 
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col,
                                      density.border.col = "blue"
)

# Add now the other distribution with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange",
               "pink", "lightblue", "cyan", "brown", "purple", "lightgreen" )
##> ILN + PART + NOT-CB
transp.col2  <- adjustcolor( col = col.plot[ 11 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_part_notcb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_part_notcb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 11 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )
##> GBM + CONC + NOT-CB
coords_2plot <- add.extra.dists( phy = plot_obj$phy$GBM_conc_notcb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$GBM_conc_notcb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = "white",
                                 density.border.col = col.plot[ 3 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )
##> ILN + CONC + NOT-CB
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_conc_notcb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_conc_notcb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = "white",
                                 density.border.col = col.plot[ 10 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "GBM, part, no cross-bracing",
                                "ILN, part, no cross-bracing",
                                "GBM, conc, no cross-bracing",
                                "ILN, conc, no cross-bracing" ),
        col = c( "blue", col.plot[c(11, 3, 10)]),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( rep( "LUCA", 2 ), rep( "LECA", 2 ),
                           rep( "LACA", 2 ), rep( "LBCA", 2 ),
                           rep( "ARCHAEPLASTIDA", 4),
                           rep( "TG-EUKARYA-MITO", 2 ) )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x),
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 2,
        names( sort_nodes2plot )[i] )    
}


## [[ PLOT SUPP FIGURE - cross-bracing B (only nodes with fossil cals) ]] ----

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 3. Part VS Conc, GBM vs ILN, fpssil cross-bracing #
#---------------------------------------------------#
# Saved in PDF: plots/LUCAdup_fosscb_suppfig.pdf
nodes_2plot <- c( 248,368, # LUCA
                  269,457, # LECA
                  249,438, # LACA
                  304,369, # LBCA
                  312,338,402,468, # ARCHAEPLASTIDA
                  329,391          # TG-EUKARYA-MITO
)
##> GBM + PART + FOSSCB
ind_mat <- which( names( plot_obj$node_ages$GBM_part_fosscb ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$GBM_part_fosscb, 
                                      #xlim.scale = c(-800,5000),  # wide
                                      #xlim.scale = c(-500,6000),  # less wide
                                      xlim.scale = c(-5000,7000), # narrow
                                      #xlim.scale = c(-6000,7000),  # very narrow
                                      node.ages = plot_obj$node_ages$GBM_part_fosscb[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.2,
                                      time.correction = 1000, 
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col,
                                      density.border.col = "blue"
)

# Add now the other distribution with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange",
               "pink", "lightblue", "cyan", "brown", "purple", "lightgreen" )
##> ILN + PART + FOSSCB
transp.col2  <- adjustcolor( col = col.plot[ 11 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_part_fosscb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_part_fosscb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 11 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )
##> GBM + CONC + FOSSCB
coords_2plot <- add.extra.dists( phy = plot_obj$phy$GBM_conc_fosscb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$GBM_conc_fosscb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = "white",
                                 density.border.col = col.plot[ 3 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )
##> ILN + CONC + FOSSCB
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_conc_fosscb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_conc_fosscb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = "white",
                                 density.border.col = col.plot[ 10 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "GBM, part, cross-bracing B",
                                "ILN, part, cross-bracing B",
                                "GBM, conc, cross-bracing B",
                                "ILN, conc, cross-bracing B" ),
        col = c( "blue", col.plot[c(11, 3, 10)]),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( rep( "LUCA", 2 ), rep( "LECA", 2 ),
                           rep( "LACA", 2 ), rep( "LBCA", 2 ),
                           rep( "ARCHAEPLASTIDA", 4),
                           rep( "TG-EUKARYA-MITO", 2 ) )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x),
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 2,
        names( sort_nodes2plot )[i] )    
}
## -- ADDITIONAL PLOTS FOR BENCHMARKING ----
## [[ PLOT CONCATENATED vs PARTITIONED DATASETS | GBM ]] ----

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# GBM | cb VS not-cb & conc VS part #
#-----------------------------------#
# Saved in PDF: plots/LUCAdup_GBMmodels.pdf
nodes_2plot <- c( 248,368, # LUCA
                  269,457, # LECA
                  249,438, # LACA
                  304,369,  # LBCA
                  312,338,402,468, # ARCHAEPLASTIDA
                  329,391          # TG-EUKARYA-MITO
)
ind_mat <- which( names( plot_obj$node_ages$GBM_conc_cb ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$GBM_conc_cb, 
                                      #xlim.scale = c(-800,5000),  # wide
                                      #xlim.scale = c(-500,6000),  # less wide
                                      xlim.scale = c(-5000,7000), # narrow
                                      #xlim.scale = c(-6000,7000),  # very narrow
                                      node.ages = plot_obj$node_ages$GBM_conc_cb[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.2,
                                      time.correction = 1000, 
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col, #density.col = "white",
                                      density.border.col = "blue"
)

# Add now the other distributions with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange",
               "pink", "lightblue", "cyan", "brown", "purple", "lightgreen" )
transp.col2  <- adjustcolor( col = col.plot[ 3 ], alpha.f = 0.2 )
# CONC + NO CROSS-BRACING
coords_2plot <- add.extra.dists( phy = plot_obj$phy$GBM_conc_notcb, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$GBM_conc_notcb[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 1000,
                                 density.col = transp.col2, #density.col = "white",
                                 density.border.col = col.plot[ 3 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )
# CONC + FOSSCB
add.extra.dists( phy = plot_obj$phy$GBM_conc_fosscb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$GBM_conc_fosscb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white",
                 density.border.col = col.plot[6],
                 distribution.height = 0.8, transparency = transparency )
# PART + CROSS-BRACED
add.extra.dists( phy = plot_obj$phy$GBM_part_cb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$GBM_part_cb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white", density.border.col = col.plot[11],
                 distribution.height = 0.8, transparency = transparency )
# PART + NO CROSS-BRACING
add.extra.dists( phy = plot_obj$phy$GBM_part_notcb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$GBM_part_notcb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white", density.border.col = col.plot[12],
                 distribution.height = 0.8, transparency = transparency )
# PART + FOSSCB
add.extra.dists( phy = plot_obj$phy$GBM_part_fosscb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$GBM_part_fosscb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white",
                 density.border.col = col.plot[4],
                 distribution.height = 0.8, transparency = transparency )
# Add legend
legend( locator(1), legend = c( "GBM, conc, cross-bracing A",
                                "GBM, conc",
                                "GBM, conc, cross-bracing B",
                                "GBM, part, cross-bracing A",
                                "GBM, part",
                                "GBM, part, cross-bracing B"
                                ),
        col = c( "blue", col.plot[c(3,6,11,12,10)]),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( rep( "LUCA", 2 ), rep( "LECA", 2 ),
                           rep( "LACA", 2 ), rep( "LBCA", 2 ),
                           rep( "ARCHAEPLASTIDA", 4),
                           rep( "TG-EUKARYA-MITO", 2 ) )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x),
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 2,
        names( sort_nodes2plot )[i] )    
}

## [[ PLOT CONCATENATED vs PARTITIONED DATASETS | ILN ]] ----

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# ILN | cb VS not-cb & conc VS part #
#-----------------------------------#
# Saved in PDF: plots/LUCAdup_ILNmodels.pdf
nodes_2plot <- c( 248,368, # LUCA
                  269,457, # LECA
                  249,438, # LACA
                  304,369, # LBCA
                  312,338,402,468, # ARCHAEPLASTIDA
                  329,391          # TG-EUKARYA-MITO
)
ind_mat <- which( names( plot_obj$node_ages$ILN_conc_cb ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$ILN_conc_cb, 
                                      #xlim.scale = c(-800,5000),  # wide
                                      #xlim.scale = c(-500,6000),  # less wide
                                      xlim.scale = c(-5000,7000), # narrow
                                      #xlim.scale = c(-6000,7000),  # very narrow
                                      node.ages = plot_obj$node_ages$ILN_conc_cb[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.2,
                                      time.correction = 1000, 
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col, #density.col = "white",
                                      density.border.col = "blue"
)

# Add now the other distributions with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange", "pink",
               "lightblue", "cyan", "brown", "purple", "lightgreen" )
# CONC + NOT CROSS-BRACED
transp.col2 <- adjustcolor( col = col.plot[ 3 ], alpha.f = 0.2 )
coords_2plot_ILN <- add.extra.dists( phy = plot_obj$phy$ILN_conc_notcb, num.models = 1,
                                     last.plot = last.plot,
                                     node.ages = plot_obj$node_ages$ILN_conc_notcb[ind_mat],
                                     plot.type = "distributions",
                                     time.correction = 1000,
                                     density.col = transp.col2, #density.col = "white",
                                     density.border.col = col.plot[ 3 ],
                                     distribution.height = 0.8, transparency = transparency,
                                     return_coords = TRUE )
# CONC + FOSSCB
add.extra.dists( phy = plot_obj$phy$ILN_conc_fosscb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$ILN_conc_fosscb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white",
                 density.border.col = col.plot[6],
                 distribution.height = 0.8, transparency = transparency )
# PART + CROSS-BRACING
add.extra.dists( phy = plot_obj$phy$ILN_part_cb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$ILN_part_cb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white", density.border.col = col.plot[11],
                 distribution.height = 0.8, transparency = transparency )
# PART + NOT CROSS-BRACED
add.extra.dists( phy = plot_obj$phy$ILN_part_notcb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$ILN_part_notcb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white", density.border.col = col.plot[12],
                 distribution.height = 0.8, transparency = transparency )
# PART + FOSSCB
add.extra.dists( phy = plot_obj$phy$ILN_part_fosscb, num.models = 1,
                 last.plot = last.plot,
                 node.ages = plot_obj$node_ages$ILN_part_fosscb[ind_mat],
                 plot.type = "distributions",
                 time.correction = 1000,
                 density.col = "white",
                 density.border.col = col.plot[4],
                 distribution.height = 0.8, transparency = transparency )

# Add legend
legend( locator(1),
        legend = c( "ILN, conc, cross-bracing A",
                    "ILN, conc",
                    "ILN, conc, cross-bracing B",
                    "ILN, part, cross-bracing A",
                    "ILN, part",
                    "ILN, part, cross-bracing B" ),
        col = c( "blue", col.plot[c(3,6,11,12,10)]),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( rep( "LUCA", 2 ), rep( "LECA", 2 ),
                           rep( "LACA", 2 ), rep( "LBCA", 2 ),
                           rep( "ARCHAEPLASTIDA", 4),
                           rep( "TG-EUKARYA-MITO", 2 ) )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x),
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 2,
        names( sort_nodes2plot )[i] )    
}
