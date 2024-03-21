#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
main_wd    <- gsub( pattern = "figs_tables/scripts/", replacement = "", x = wd )
wd_conc_cb <- gsub( pattern = "figs_tables/scripts/", replacement = "../01_PAML/02_MCMCtree/sum_files_post/", x = wd )
wd_ind     <- gsub( pattern = "figs_tables/scripts/", replacement = "00_indgenes_analyses/01_PAML/00_MCMCtree/sum_files_post/", x = wd )
wd_noATP   <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noatp/01_PAML/01_MCMCtree/sum_files_post/", x = wd )
wd_noEF    <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noef/01_PAML/01_MCMCtree/sum_files_post/", x = wd )
wd_noLeu   <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_noleu/01_PAML/01_MCMCtree/sum_files_post/", x = wd )
wd_noSRP   <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_nosrp/01_PAML/01_MCMCtree/sum_files_post/", x = wd )
wd_noTyr   <- gsub( pattern = "figs_tables/scripts/", replacement = "01_lou_analyses/LUCA-divtimes_notrptyr/01_PAML/01_MCMCtree/sum_files_post/", x = wd )
wd_bsinbv  <- gsub( pattern = "figs_tables/scripts/", replacement = "02_blunc_analysis/02_PAML/MCMCtree/sum_files_post/", x = wd )
all_wds <- c( wd_conc_cb,
              rep( wd_ind, 5),
              wd_noATP, wd_noEF, wd_noLeu, wd_noSRP, wd_noTyr, 
              wd_bsinbv )
source( "../../../src/Functions_plots.R" )

#### LOAD DATA ----
#-------------#
# START TASKS #
#-------------#
# Load all output csv tables previously generated
dir_names <- c( "conc_cb",
                "indATP", "indEF", "indLeu", "indSRP", "indTyr",
                "noATP", "noEF", "noLeu", "noSRP", "noTyr",
                "bsinBV" )
patt_csv <- c( "",
               "ATP", "EF", "Leu", "SRP", "Tyr",
               rep( "", 5),
               "" )
all_divt  <- divt_csv( all_paths = all_wds, name_dirs = dir_names,
                       pattern_csv = patt_csv )
# Select LUCA nodes that are to be plotted:
nodes_2plot <- paste( "t_n", c( 248, 368, # conc
                                rep( c(248, 368), 5), # indgenes
                                243, 359, # ATP
                                238, 351, # EF
                                245, 364, # Leu
                                248, 368, # SRP
                                246, 365, # Tyr
                                248, 368  # bsinbv
), sep = "" )
plots_per_fig <- rep( 2, 12 )
names( nodes_2plot ) <- c( paste( c( "LUCA", "LUCA-dup" ), "_main-conc",
                                  sep = "" ),
                           paste( paste( c( "LUCA", "LUCA-dup" ), "_",
                                         sep = "" ),
                                  c( "ATP", "ATP", "EF", "EF", "Leu", "Leu",
                                     "SRP", "SRP", "Tyr", "Tyr" ),
                                  sep = ""  ),
                           paste( paste( c( "LUCA", "LUCA-dup" ), "_",
                                         sep = "" ),
                                  c( "noATP", "noATP", "noEF", "noEF", 
                                     "noLeu", "noLeu", "noSRP", "noSRP", 
                                     "noTyr","noTyr" ),
                                  sep = ""  ),
                           paste( c( "LUCA", "LUCA-dup" ), "_main-bsinbv",
                                  sep = "" ) )

#### TABLES ----
#------------------#
# COMPARISON TABLE #
#------------------#
# Get the total number of columns, subtract one as we are not going to keep
# the fourth column with the calibrations (only a check when running 
# MCMC diagnostics)
short_GBM <- short_ILN <- matrix( 0, nrow = length( nodes_2plot ),
                                  ncol = 3 )
rownames( short_GBM ) <- rownames( short_ILN ) <- 
  paste( nodes_2plot, "|", names( nodes_2plot ), sep = "" )
# Get colnames ready
colnames( short_GBM ) <- colnames( short_ILN ) <-
  c( "mean_t", "2.5%q", "97.5%q" )
# Start saving the output results in the empty vectors for which memory has
# already been allocated in the corresponding cells
start <- end <- count <- start2 <- end2 <- 0
conc_names <- vector( mode = "character", length = length(dir_names)*3 )
for( i in 1:length( dir_names ) ){
  if( i == 1 ){
    start <- 1 
    end   <- start + 1
    start2 <- 1
    end2   <- start2 + 3
  }else{
    start <- end + 1
    end   <- start + 1
    start2 <- end2 + 1
    end2   <- start2 + 3
  }
  cat( "round ", i, "start = ", start, " | end = ", end, "\n" )
  for( j in nodes_2plot[c(start:end)] ){
    # Get indexes
    count <- count + 1
    only_nums <- gsub( x = rownames( all_divt$GBM[[ i ]] ),
                       pattern = "_[A-Z]..*", replacement = "" )
    tmp_ind <- which( only_nums %in% j )
    short_GBM[count,1:3] <- as.numeric( all_divt$GBM[[ i ]][tmp_ind,1:3] )
    short_ILN[count,1:3] <- as.numeric( all_divt$ILN[[ i ]][tmp_ind,1:3] )
  }
  conc_names[start2:end2] <- rep( dir_names[i], 4 )
}
# Output the summary tsv files
if( ! dir.exists( paste( main_wd, "figs_tables/tables", sep = "" ) ) ){
  dir.create( paste( main_wd, "figs_tables/tables", sep = "" ) ) 
}
write.table( x = short_GBM, file = paste( main_wd, "figs_tables/tables/",
                                          "compare_divtimes_GBM.tsv",
                                          sep = "" ),
             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE )
write.table( x = short_ILN, file = paste( main_wd, "figs_tables/tables/",
                                          "compare_divtimes_ILN.tsv",
                                          sep = "" ),
             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE )

# Do the same but for all the nodes
num_nodes_long <- length( grep( pattern = "[0-9]_[A-Z]*",
                                x = rownames( all_divt$GBM[[ 1 ]] ) ) )
long_GBM <- long_ILN <- matrix( 0, nrow = num_nodes_long,
                                ncol = 4*(length( all_divt$GBM )) )
colnames( long_GBM ) <- colnames( long_ILN ) <- 
  paste( rep( c( "label", "mean_t", "2.5%q", "97.5%q" ), length( dir_names ) ),
         conc_names, sep = "-" )
start <- end <- 0
for( i in 1:length( all_divt$GBM ) ){
  ind_GBM <- grep( pattern = "[0-9]_[A-Z]*",
                   x = rownames( all_divt$GBM[[ i ]] ) )
  ind_ILN <- grep( pattern = "[0-9]_[A-Z]*",
                   x = rownames( all_divt$ILN[[ i ]] ) )
  start   <- end + 1
  end     <- start+ 2 + 1
  end_row <- length( rownames( all_divt$GBM[[ i ]] )[ind_GBM] )
  long_GBM[1:end_row,start]   <- rownames( all_divt$GBM[[ i ]] )[ind_GBM]
  long_GBM[1:end_row,start+1] <- as.numeric( all_divt$GBM[[ i ]][ind_GBM,1] )
  long_GBM[1:end_row,start+2] <- as.numeric( all_divt$GBM[[ i ]][ind_GBM,2] )
  long_GBM[1:end_row,end]     <- as.numeric( all_divt$GBM[[ i ]][ind_GBM,3] )
  long_ILN[1:end_row,start]   <- rownames( all_divt$ILN[[ i ]] )[ind_ILN]
  long_ILN[1:end_row,start+1] <- as.numeric( all_divt$ILN[[ i ]][ind_ILN,1] )
  long_ILN[1:end_row,start+2] <- as.numeric( all_divt$ILN[[ i ]][ind_ILN,2] )
  long_ILN[1:end_row,end]     <- as.numeric( all_divt$ILN[[ i ]][ind_ILN,3] )
}
# Replace `0` (no cals) with label
ind_0s_GBM <- which( long_GBM == 0, arr.ind = TRUE )
ind_0s_ILN <- which( long_ILN == 0 )
long_GBM[ind_0s_GBM] <- "nocal"
long_ILN[ind_0s_ILN] <- "nocal"
write.table( x = long_GBM, file = paste( main_wd, "figs_tables/tables/",
                                          "compare_divtimes_allcalnodes_GBM.tsv",
                                          sep = "" ),
             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
write.table( x = long_ILN, file = paste( main_wd, "figs_tables/tables/",
                                          "compare_divtimes_allcalnodes_ILN.tsv",
                                          sep = "" ),
             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )

#### PLOTS ----
# Create `plots` directory and PDF file
if( ! dir.exists( paste( wd, "../plots", sep = "" ) ) ){
  dir.create( paste( wd, "../plots", sep = "" ) )
}
pdf( paste( wd, "../plots/divtimes_LUCA_compare_all.pdf", sep = "" ),
     width=15, height=10 )
# Find max and min CIs
start <- end <- 0
divtmeanD_2plot <- divtlowqD_2plot <- divtupqD_2plot <- 
  divtmeanM_2plot <- divtlowqM_2plot <- divtupqM_2plot <-
  rep( 0 , length( nodes_2plot ))
for( j in 1:length(all_divt$GBM) ){
  start <- end + 1
  end   <- start + 1
  lab_nodes <- gsub( x = rownames( all_divt$GBM[[ j ]] ),
                     pattern = "_[A-Z]..*", replacement = "" )
  tmp_ind_lab <- which( lab_nodes %in% nodes_2plot[start:end] )
  names( tmp_ind_lab ) <- names( nodes_2plot[start:end] )
  print( tmp_ind_lab )
  #Mean divt
  divtmeanD_2plot[start] <- all_divt$GBM[[ j ]][tmp_ind_lab[1],1]
  divtmeanD_2plot[end]   <- all_divt$ILN[[ j ]][tmp_ind_lab[1],1]
  divtmeanM_2plot[start] <- all_divt$GBM[[ j ]][tmp_ind_lab[2],1]
  divtmeanM_2plot[end]   <- all_divt$ILN[[ j ]][tmp_ind_lab[2],1]
  names( divtmeanD_2plot )[start:end] <- 
    paste( nodes_2plot[start], c( "_GBM", "_ILN" ), sep = "" )
  names( divtmeanM_2plot )[start:end] <- 
    paste( nodes_2plot[end], c( "_GBM", "_ILN" ), sep = "" )
  #Low quantiles
  divtlowqD_2plot[start] <- all_divt$GBM[[ j ]][tmp_ind_lab[1],2]
  divtlowqD_2plot[end]   <- all_divt$ILN[[ j ]][tmp_ind_lab[1],2]
  divtlowqM_2plot[start] <- all_divt$GBM[[ j ]][tmp_ind_lab[2],2]
  divtlowqM_2plot[end]   <- all_divt$ILN[[ j ]][tmp_ind_lab[2],2]
  names( divtlowqD_2plot )[start:end] <- 
    paste( nodes_2plot[start], c( "_GBM", "_ILN" ), sep = "" )
  names( divtlowqM_2plot )[start:end] <- 
    paste( nodes_2plot[end], c( "_GBM", "_ILN" ), sep = "" )
  #Upper quantiles
  divtupqD_2plot[start]  <- all_divt$GBM[[ j ]][tmp_ind_lab[1],3]
  divtupqD_2plot[end]    <- all_divt$ILN[[ j ]][tmp_ind_lab[1],3]
  divtupqM_2plot[start]  <- all_divt$GBM[[ j ]][tmp_ind_lab[2],3]
  divtupqM_2plot[end]    <- all_divt$ILN[[ j ]][tmp_ind_lab[2],3]
  names( divtupqD_2plot )[start:end] <- 
    paste( nodes_2plot[start], c( "_GBM", "_ILN" ), sep = "" )
  names( divtupqM_2plot )[start:end] <- 
    paste( nodes_2plot[end], c( "_GBM", "_ILN" ), sep = "" )
  if( j == 1 ){
    tmpy_CImax_GBM <- all_divt$GBM[[ j ]][tmp_ind_lab,3]
    tmpy_CImin_GBM <- all_divt$GBM[[ j ]][tmp_ind_lab,2]
    tmpy_CImax_ILN <- all_divt$ILN[[ j ]][tmp_ind_lab,3]
    tmpy_CImin_ILN <- all_divt$ILN[[ j ]][tmp_ind_lab,2]
  }else{
    tmpy_CImax_GBM <- max( tmpy_CImax_GBM, all_divt$GBM[[ j ]][tmp_ind_lab,3] )
    tmpy_CImin_GBM <- min( tmpy_CImin_GBM, all_divt$GBM[[ j ]][tmp_ind_lab,2] )
    tmpy_CImax_ILN <- max( tmpy_CImax_ILN, all_divt$ILN[[ j ]][tmp_ind_lab,3] )
    tmpy_CImin_ILN <- min( tmpy_CImin_ILN, all_divt$ILN[[ j ]][tmp_ind_lab,2] )
  }
}
# Now, plot the CIs
par( mfrow = c( 2, 1 ) )
min_y  <- min( tmpy_CImin_GBM, tmpy_CImin_ILN )
max_y  <- max( tmpy_CImax_GBM, tmpy_CImax_ILN )
x_vals <- 1:length( nodes_2plot )
cols_plot <- c( "#E2ED71", "#E2ED71",
                "#D81B60", "#D81B60",
                "#1E88E5", "#1E88E5",
                "#7B38FD", "#7B38FD",
                "#FFC107", "#FFC107",
                "#004D40", "#004D40",
                "#FCA095", "#FCA095",
                "#281ACD", "#281ACD",
                "#6A5F7E", "#6A5F7E",
                "#964067", "#964067", 
                "#31DE93", "#31DE93",
                "#BD7254", "#BD7254" )
xlabs <- gsub( x = gsub( pattern = "LUCA_", replacement =  "GBM_",
                         x = names(nodes_2plot) ),
               pattern = "LUCA-dup", replacement = "ILN" )
## PLOT DRIVER NODES
plot( x = x_vals, 
      y = divtmeanD_2plot,
      pch = 16, ylim = c( min_y, max_y ),
      xaxt = "n", xlab = "",
      ylab = "Estimated mean divergence times",
      col = cols_plot )
title( main = "Comparing LUCA time estimates under various strategies - driver node" )
axis( side = 1, at = 1:length( nodes_2plot ),
      labels = FALSE, las = 2 )
text( x = x_vals, # offset calculated by eye
      labels = xlabs,
      ## Rotate the labels by 25 degrees.
      par("usr")[3]-3, 
      srt = 25, adj = 1, xpd = TRUE )
# For 97.5%CI, use triangle facing down
points( x = x_vals, y = divtlowqD_2plot, pch = 25,
        cex = 0.7,
        col = cols_plot )
# For 2.75%CI, use triangle facing up
points( x = x_vals, y = divtupqD_2plot, pch = 2,
        cex = 0.7,
        col = cols_plot )
# Add lines!
start <- end <- 0
for( i in 1:(length(nodes_2plot)/2) ){
  start <- end + 1
  end   <- start + 1
  lines( x = x_vals[start:end], y = divtmeanD_2plot[start:end], pch = 25,
         lty = 1, col = "black" )
  lines( x = x_vals[start:end], y = divtlowqD_2plot[start:end], pch = 25,
         lty = 2, col = "black" )
  lines( x = x_vals[start:end], y = divtupqD_2plot[start:end], pch = 25,
         lty = 2, col = "black" )
}
## PLOT MIRROR NODES
plot( x = x_vals, 
      y = divtmeanM_2plot,
      pch = 16, ylim = c( min_y, max_y ),
      xaxt = "n", xlab = "",
      ylab = "Estimated mean divergence times",
      col = cols_plot )
title( main = "Comparing LUCA time estimates under various strategies - mirrored node" )
axis( side = 1, at = 1:length( nodes_2plot ),
      labels = FALSE, las = 2 )
text( x = x_vals, # offset calculated by eye
      labels = xlabs,
      ## Rotate the labels by 25 degrees.
      par("usr")[3]-3, 
      srt = 25, adj = 1, xpd = TRUE )
# For 97.5%CI, use triangle facing down
points( x = x_vals, y = divtlowqM_2plot, pch = 25,
        cex = 0.7,
        col = cols_plot )
# For 2.75%CI, use triangle facing up
points( x = x_vals, y = divtupqM_2plot, pch = 2,
        cex = 0.7,
        col = cols_plot )
# Add lines!
start <- end <- 0
for( i in 1:(length(nodes_2plot)/2) ){
  start <- end + 1
  end   <- start + 1
  lines( x = x_vals[start:end], y = divtmeanM_2plot[start:end], pch = 25,
         lty = 1, col = "black" )
  lines( x = x_vals[start:end], y = divtlowqM_2plot[start:end], pch = 25,
         lty = 2, col = "black" )
  lines( x = x_vals[start:end], y = divtupqM_2plot[start:end], pch = 25,
         lty = 2, col = "black" )
}
# Close graphics
dev.off()

