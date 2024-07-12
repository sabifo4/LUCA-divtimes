#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# This package lets you find automatically the path to a specific location
# in your file structure
# If you have not installed this package, you will need to install it. 
# You can uncomment the following line to do this:
#install.packages( "rstudioapi" )
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
# This script and the `Functions.R` script are saved under
# a directory called `scripts`. We will remove this part from 
# the path to find our home directory
home_dir <- gsub( pattern = "scripts/", replacement = "", x = wd )
# Load main script with all functions required below
source( file = "../../../src/Functions.R" )

#-------------------------------------------------------------#
# DEFINE GLOBAL VARIABLES -- modify according to your dataset #
#-------------------------------------------------------------#
# First, we will define global variables that we will keep using throughout this
# script.

# 1. Label the file with calibrations. If you have tested different calibrations
# and have more than one file with the corresponding calibrations, give as 
# many labels as files you have.
dat <- c( "LUCAdup_part_notcb" )

# 2. Number of divergence times that have been estimated. One trick to find
# this out quickly is to subtract 1 to the number of species. In this case,
# there are 246 taxa (246), so the number of internal nodes
# is `n_taxa-=246-1=245`.
# Another way to verify this is by opening the `mcmc.txt` file and check the
# header. The first element after `Gen` will have the format of `t_nX`, where
# X will be an integer (i.e., 247). Subtract two to this number 
# (i.e., 247-2=245) and this will be your number of divergence times that are 
# parameters of the MCMC. Please modify the number below so it fits to the 
# dataset you are using. 
num_divt <- 245

# 3. Total number of samples that you collected after generating the
# final `mcmc.txt` files with those from the chains that passed the filters. 
# You can check these numbers in scripts `MCMC_diagnostics_posterior.R` and
# `MCMC_diagnostics_prior.R`. E.g., `sum_post_QC$<name_dataset>$total_samples`
# or `sum_prior_QC$<name_dataset>$total_samples`
#
# CLK: The number of lines is 60003, and so you need to specify one less
# GBM: The number of lines is 45649, and so you need to specify one less
# ILN: The number of lines is 48216, and so you need to specify one less
#
# NOTE: If you had more than one dataset, you would add another vector of three
# values with the samples for CLK, GBM, and ILN to create `def_samples`
# E.g. two datasts: c( c( 120005, 120005, 120005), c( 120005, 120005, 120005) )
def_samples <- c( 60002, 45648, 48215 )

# 4. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles (i.e., 95%CI). If you want to
# change this, however, just modify the value.
perc <- 0.975

# 5. Number of columns in the `mcmc.txt` that are to be deleted as they do not 
# correspond to sample values for divergence times (i.e., the entries are not 
# names following the format `t_nX`). To figure out this number quickly, you 
# can open the `mcmc.txt` file, read the header, and count the number of `mu*`
# and `sigma2*` elements. Do not count the `lnL` value when looking at 
# `mcmc.txt` files generated when sampling from the posterior -- this is 
# automatically accounted for in the in-house R functions that you will 
# subsequently use. E.g., assuming an MCMC ran under a relaxed-clock model with  
# no partitions, we would see `mu` and `sigma2` columns. Therefore, the variable  
# would be set to `delcol_post <- 2`. Please modify the value/s below 
# (depending on having one or more datasets) according to the `mcmc.txt` file
# generated when sampling from the posterior (`delcol_obj`). When running
# from the prior and `clock = 1`, you will only see `mu*` columns but, if you
# ran it with options `clock = 2` or `clock = 3`, you shall also see `sigma2*`
# columns.
##> NOTE: If you ran `MCMCtree` with `clock = 2` or `clock = 3` when
##> sampling from the prior, you will also need to count the `sigma2*`
##> columns! We ran `clock = 1` so that the analyses ran quicker, and thus
##> we only have `mu*` columns.
delcol_obj <- c( 5, 10, 10 ) # prior (1), posterior (2), posterior (2)

# 6. Path to the directory where the concatenated `mcmc.txt` file has been 
# generated. Note that, if you have run more than one chain in `MCMCtree` for
# each hypothesis tested, you are expected to have generated a concatenated 
# `mcmc.txt` file with the bash script `Combine_MCMC_prior.sh` or any similar 
# approaches.
num_dirs     <- 3
paths_dat    <- c( paste( home_dir,
                          "sum_analyses/00_prior/mcmc_files_part_notcb_CLK",
                          sep = "" ),
                   paste( home_dir,
                          "sum_analyses/01_posterior/mcmc_files_part_notcb_GBM",
                          sep = "" ),
                   paste( home_dir,
                          "sum_analyses/01_posterior/mcmc_files_part_notcb_ILN",
                          sep = "" )
)

# 7. Load a semicolon-separated file with info about calibrated nodes. Note that
# this file is output by script `Merge_node_labels.R`. A summary of its content
# in case you are to generate your own input files:
#
# Each column needs to be separated with semicolons and an extra blank line
# after the last row with calibration information needs to be added. If the
# extra blank is not added, R will complain and will not load the file!
# If you add a header, please make sure you name the column elements as 
# `Calib;node;Prior`. If not, the R function below will deal with the header,
# but make sure you set `head_avail = FALSE` when running `read_calib_f` 
# function below. An example of the content of this file is given below:
#
# ```
# Calib;node;Prior
# ex_n5;5;ST(5.8300,0.0590,0.1120,109.1240)
# ex_n7;7;B(4.1200,4.5200,0.0250,0.0250)
#
# ```
#
# The first column will have the name of the calibration/s that can help you
# identify which node belongs to which calibration. The second column is the
# number given to this node by`MCMCtree` (this information is automatically
# found when you run the script `Merge_node_labels.R`, otherwise you will need
# to keep checking the output file `node_trees.tree` to figure out which node
# is which). The third column is the calibration used for that node in
# `MCMCtree` format.
# 
# [[ NOTES ABOUT ALLOWED CALIBRATION FORMATS]]
#
# Soft-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the default tail probabilities would have the following equivalent 
#         formats:
#         >> B(0.6,0.8) | B(0.6,0.8,0.025,0.025)
#  E.g.2: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the pL=0.001 and pU=0.025 would have the following format. Note that, 
#         whenever you want to modify either pL or pU, you need to write down 
#         the four  parameters in the format of "B(min,max,pL,pU)":
#         >> B(0.6,0.8,0.001,0.025)
#
# Lower-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and the default parameters for
#         p = 0.1, c = 1, pL = 0.025:
#         >> L(0.6) | L(0.6,0.1,1,0.025)
#  E.g.2: A calibration with a hard minimum at 0.6, and so pL = 1e-300. 
#         Note that, whenever you want to modify either pL or pU, you need to  
#         write down the four parameters in the format of "L(min,p,c,pL)":
#         >> L(0.6,0.1,1,1e-300)
#
# Upper-bound calibrations: 
#  E.g.1: A calibration with a maximum of 0.8 and the default parameters for
#         pU = 0.025:
#         >> U(0.8) | U(0.8,0.025)
#  E.g.2: A calibration with a hard maximum at 0.8, and so pU = 1e-300. 
#         Note that, if you want to modify pU, you need to write down the two
#         parameters in the format of "U(max,pU)":
#         >> U(0.8,1e-300)
#
# ST distributions: 
#  The format accepted has four parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape), nu (df). Accepted format: 
#  >> ST(5.8300,0.0590,0.1120,109.1240)
#
# SN distributions: 
#  The format accepted has three parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape). Accepted format: 
#  >> SN(5.8300,0.0590,0.1120)  
#
#
# The next command executes the `read_calib_f` in-house function, which reads
# your input files (semicolon-separated files). The path to this directory is 
# what the argument `main_dir` needs. The argument `f_names` requires the name 
# of the file/s that you have used. Argument `dat` requires the same global 
# object that you have created at the beginning of the script.
calib_nodes <- read_calib_f( main_dir = paste( home_dir, "calib_files/",
                                               sep = "" ),
                             f_names = "Calibnodes.csv",
                             dat = dat, head_avail = TRUE )

#-----------#
# LOAD DATA #
#-----------#
# Load mcmc files for all datasets
mcmc_obj <- vector( "list", num_dirs )
names( mcmc_obj ) <- c( "CLK", "GBM", "ILN" )
prior <- c( TRUE, FALSE, FALSE )
count <- 0
for( i in c( paths_dat ) ){
  count <- count + 1
  cat( "[[ Parsing file for dataset", names( mcmc_obj )[count], " ]]\n" )
  mcmc_obj[[count]] <- load_dat( mcmc = paste( i, "/mcmc.txt", sep = "" ),
                                 delcol = delcol_obj[count], perc = perc,
                                 def_samples = def_samples[count],
                                 prior = prior[count] )
}

#---------------------------#
# PLOTS: prior VS posterior #
#---------------------------#
# Prepare object with all the cross-braced nodes and the corresponding labels
dup_dat          <- vector( mode = "list", 13 )
names( dup_dat ) <- c( "LUCA", "TG-EUKARYA-ARCH", "LECA", "FUNGI", "METAZOA",
                       "EUMETAZOA", "ARCHAEPLASTIDA", "EMBRYOPHYTA",
                       "EUDICOT-MONOCOT", "CG-FORAMINIFERA",
                       "TG-OXYPHOTOBACTERIA", "CG-OXYPHOTOBACTERIA",
                       "TG-EUKARYA-MITO" )
dup_dat[[1]] <- paste( "t_n", c(248,368), sep = "" )
dup_dat[[2]] <- paste( "t_n", c(268,456), sep = "" )
dup_dat[[3]] <- paste( "t_n", c(269,457), sep = "" )
dup_dat[[4]] <- paste( "t_n", c(272,332,395,460), sep = "" )
dup_dat[[5]] <- paste( "t_n", c(274,333,397,462), sep = "" )
dup_dat[[6]] <- paste( "t_n", c(275,334,398,463), sep = "" )
dup_dat[[7]] <- paste( "t_n", c(280,312,338,402,426,468), sep = "" )
dup_dat[[8]] <- paste( "t_n", c(283,314,340,404,430,471), sep = "" )
dup_dat[[9]] <- paste( "t_n", c(284,315,341,405,472), sep = "" )
dup_dat[[10]] <- paste( "t_n", c(292,413,480), sep = "" )
dup_dat[[11]] <- paste( "t_n", c(308,422), sep = "" )
dup_dat[[12]] <- paste( "t_n", c(309,423), sep = "" )
dup_dat[[13]] <-  paste( "t_n", c(329,391), sep = "" )
# Get objects with the labels in the csv file and the corresponding calib
# names
name_labs <- colnames( mcmc_obj$CLK$divt )
labs_in_csv <- paste( "t_n", calib_nodes$LUCAdup_part[,2], sep = "" )
# Start plotting!
for( i in 1:length(dup_dat) ){
  pdf( file = paste( home_dir, "plots/priorVSpost_", names(dup_dat)[i],
                     "_plot.pdf",
                     sep = "" ), 
       paper = "a4r", width = 0, height = 0 )
  # Find total num of duplications and indexes
  nod_ind   <- which( name_labs %in% dup_dat[[ i ]] )
  tot_nodes <- length( nod_ind )
  if( tot_nodes == 2 ){
    par( mfrow = c(1, 2 ) )
  }else if( tot_nodes == 3 |  tot_nodes == 4 ){
    par( mfrow = c(2,2) )
  }else if( tot_nodes == 5 | tot_nodes == 6 ){
    par( mfrow = c(3,2) )
  }
  # Get calibration density in correct formatting (4 nums)
  tn_ind   <- which( labs_in_csv %in% dup_dat[[ i ]] )
  tmp_allB <- as.numeric( stringr::str_split( string = gsub( pattern = "B\\(|\\)", 
                                                             replacement = "",
                                                             x = calib_nodes$LUCAdup_part[tn_ind,3] ),
                                              pattern = "," )[[1]] )
  # Get max/min x and y
  ##> NOTE: Only the first number is used, because the other number/s in 
  ##> `dup_dat[[i]` have the same fossil info!
  tn_matrix  <- which( colnames(mcmc_obj$GBM$divt) %in% dup_dat[[ i ]] )
  max_y <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$y ),
                   max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$y ),
                   max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$y ) ) )
  max_x <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                   max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                   max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                   tmp_allB[2]) )
  min_x <- min( c( min( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                   min( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                   min( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                   tmp_allB[1] ) )
  # Get name of calibration, which corresponds to the first calibration
  # in the list
  tmp_ind_name_calib <- which( labs_in_csv %in% dup_dat[[ i ]][1] )
  tmp_name_calib     <- calib_nodes$LUCAdup_part[tmp_ind_name_calib,1]
  # Start plotting calib density, marginal density, and posterior densities for
  # each cross-braced node -- easy for comparisons!
  for( j in dup_dat[[ i ]] ){
    # Get temporary matching column in `mcmc_obj` to `j`
    tmp_tn         <- which( colnames(mcmc_obj$GBM$divt) %in% j )
    plot( density( mcmc_obj$GBM$divt[[ tmp_tn ]], adj = 1 ),
          main = paste( "Comparing densities (col",
                        tmp_tn, ") | ", tmp_name_calib, sep = "" ),
          col = "blue", ylim = c(0,max_y), xlim = c(min_x, max_x) )
    lines( density( mcmc_obj$CLK$divt[[ tmp_tn ]], adj = 1 ), col = "brown" )
    lines( density( mcmc_obj$ILN$divt[[ tmp_tn ]], adj = 1 ),
           col = "darkolivegreen3" )
    curve( mcmc3r::dB( x, tL = tmp_allB[1], tU = tmp_allB[2],
                       pL = tmp_allB[3], pU = tmp_allB[4] ),
           from = tmp_allB[1]-0.5, to = tmp_allB[2]+0.5,
           n = 1e5, add = TRUE, col = "black" )
    legend( "topleft", legend = c( "Calibration density", "Marginal density",
                                   "Post-GBM", "Post-ILN"),
            lwd = 1, bty = "n", cex = 1,
            col = c( "black", "brown", "blue", "darkolivegreen3" ) )
    cols_vec <- c( "blue", "red", "purple", "darkgreen", "brown" )
  }
  # Close pdf
  dev.off()
}

