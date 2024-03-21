#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#---------------------------------------------------#
# SET WORKING DIRECTORY AND LOAD IN-HOUSE FUNCTIONS #
#---------------------------------------------------#
# Load colour-friendly package
# If you have not installed this package, you will need to install it. 
# You can uncomment the following line to do this:
#install.packages( "ColorBlindness" )
library( colorBlindness )
# This package lets you find automatically the path
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
source( file = "../../../../../src/Functions.R" )

#--------------------------------#
# DEFINE USER'S GLOBAL VARIABLES #
#--------------------------------#
# First, we will define global variables that we will keep using throughout this
# script.

# 1. Label the file with calibrations. If you have tested different calibrations
# and have more than one file with the corresponding calibrations, give as 
# many labels as files you have.
dat <- c( "ATP", "EF", "Leu", "SRP", "Tyr" )

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

# 3. Number of samples that you specified in the `MCMCtree` control file to 
# collect. NOTE: you may have not collect them all, but do not worry!
# In this case, you are using the concatenated files. The number of
# lines is 120006, and so you need to specify one less here:
def_samples <- 120005

# 4. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles. If you want to change this,
# however, just modify the value.
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
# would be set to `delcol = 2`. Please modify the values below according to your  
# the `mcmc.txt` file generated when sampling from the prior (`delcol_prior`) 
# dataset for and when sampling from the posterior (`delcol_posterior`).
delcol_prior   <- 1 # There is one column: mu

# Path to the directory where the concatenated `mcmc.txt` file has been 
# generated. Note that, if you have run more than one chain in `MCMCtree` for
# each hypothesis tested, you are expected to have generated a concatenated 
# `mcmc.txt` file with the bash script `Combine_MCMC_prior.sh` or any similar 
# approaches.
num_dirs    <- 5
path_prior <- vector( mode = "character", length = length( dat ) )
for( i in 1:length( dat ) ){
  path_prior[i]  <- paste( home_dir, "sum_analyses/00_prior/mcmc_files_",
                           dat[i], "_CLK/", sep = "" )
}

# Load semicolon-separated file/s with info about calibrated nodes. Note that
# each column needs to be separated with semicolons and an extra blank line
# after the last row with calibration information needs to be added (i.e., files
# need to have an extra blank line so R does not complain when reading them). 
# If you add a header, please make sure you name the column elements as 
# `Calib;node;Prior`. If not, the R function below will deal with the header. 
# An example of the format you need to follow to summarise the calibration info
# for each node is the following:
#
# ```
# Calib;node;Prior
# ex_n5;5;ST(5.8300,0.0590,0.1120,109.1240)
# ex_n7;7;B(4.1200,4.5200,0.0250,0.0250)
#
# ```
#
# The first column should have the name of the calibration (e.g., Afrotheria, 
# Laurasiatheria, etc.) as it will help you identify which plot belongs to which
# calibration. The second column is the node used in MCMCtree. The third column
# is the calibration used for that node in MCMCtree format.
# 
# [[ NOTES ABOUT ALLOWED CALIBRATION FORMATS]]
#
# Soft-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the default tail probabilities would have the following equivalent 
#         formats:
##        >> B(0.6,0.8) | B(0.6,0.8,0.025,0.025)
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
##        pU = 0.025:
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
# your input files (semicolon-separated files). Please save all your calibration 
# files (if more than one) in the same directory. The path to this directory is 
# what the argument `main_dir` needs. The argument `f_names` requires the name 
# of the file/s that you have used. Argument `dat` requires the same global 
# object that you have created at the beginning of the script. If your input  
# files have a header, please keep `head_avail = TRUE`. Otherwise, change this  
# to FALSE.
dat_ff <- list.files( path = "../calib_files/", pattern = "_effVSuser.csv",
                      full.names = FALSE )
calib_nodes <- read_calib_f( main_dir = paste( home_dir, "calib_files/",
                                               sep = "" ),
                             f_names = dat_ff,
                             dat = c("LUCAdup_indgenes"), head_avail = TRUE )

#-----------#
# LOAD DATA #
#-----------#
# Load mcmc files for all datasets
mcmc_priors <- vector( "list", num_dirs )
names( mcmc_priors ) <- paste( "CLK_", dat, sep = "" )
count <- 0
for( i in c( path_prior) ){
  count <- count + 1
  cat( "[[ Parsing file for dataset", names( mcmc_priors )[count], " ]]\n" )
  mcmc_priors[[count]] <- load_dat( mcmc = paste( i, "mcmc.txt", sep = "" ),
                                    delcol = delcol_prior, perc = perc,
                                    def_samples = def_samples, prior = TRUE )
}

#------------------------------------------------#
# PLOTS: calibration density VS marginal density #
#------------------------------------------------#
# Plot calibration density VS marginal density
if( ! dir.exists( paste( home_dir, "plots", sep = "" ) ) ){
  dir.create( paste( home_dir, "plots", sep = "" ) )
}
if( ! dir.exists( paste( home_dir, "plots/margVScalib", sep = "" ) ) ){
  dir.create( paste( home_dir, "plots/margVScalib", sep = "" ) )
}

# Please write down how many rows and how many columns
# the plot with the summary of all calibrated nodes should have.
num_rows <- 4
num_cols <- 4
# Run function so the calibration density VS marginal density plots are 
# generated for each dataset
for( i in 1:num_dirs ){
  cat( "\n[[ Generating plots for dataset ",
       names( mcmc_priors )[i], " ]]\n" )
  cat( "\n[[ Output plots in PDF format]]\n" )
  plot_check_calibnodes( calibs = calib_nodes[[ 1 ]], # same calibs for all!
                         divt_list = mcmc_priors[[ i ]], 
                         dat = dat[i], out = names( mcmc_priors )[i],
                         clock = "CLK", main_wd = home_dir, ind = TRUE,
                         n_row = num_rows, n_col = num_cols+1,
                         out_format = "pdf" )
  cat( "\n[[ Output plots in JPG format]]\n" )
  plot_check_calibnodes( calibs = calib_nodes[[ 1 ]], # same calibs for all!
                         divt_list = mcmc_priors[[ i ]], 
                         dat = dat[i], out = names( mcmc_priors )[i],
                         clock = "CLK", main_wd = home_dir, ind = TRUE,
                         n_row = num_rows, n_col = num_cols+1,
                         out_format = "jpg" )
  cat( "\n[[ Output plots in TIFF format]]\n" )
  plot_check_calibnodes( calibs = calib_nodes[[ 1 ]], # same calibs for all!
                         divt_list = mcmc_priors[[ i ]], 
                         dat = dat[i], out = names( mcmc_priors )[i],
                         clock = "CLK", main_wd = home_dir, ind = TRUE,
                         n_row = num_rows, n_col = num_cols+1,
                         out_format = "tiff" )
  
}

# Plot nodes which age is constrained by a cross-calibration
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
name_labs <- colnames( mcmc_priors[[ k ]]$divt ) # same name labs for all datasets
count <- 0
for( k in 1:length(dat) ){
  for( i in 1:length(dup_dat) ){
    pdf( file = paste( home_dir, "plots/dupnodes_", dat[k], "_",
                       names(dup_dat)[i], "_plot.pdf", sep = "" ), 
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
    max_x <- max( density( mcmc_priors[[ k ]]$divt[[ nod_ind[1] ]], adj = 1 )$x ) + 0.2
    min_x <- min( density( mcmc_priors[[ k ]]$divt[[ nod_ind[1] ]], adj = 1 )$x ) - 0.2
    plot( density( mcmc_priors[[ k ]]$divt[[ nod_ind[1] ]], adj = 1 ),
          xlim = c(min_x,max_x),
          main = paste( "Compare ages for ", names(dup_dat)[i], " - ",
                        dup_dat[[ i ]][1],
                        sep = "" ) )
    cols_vec <- c( "blue", "red", "purple", "darkgreen", "brown" )
    count <- 0 
    for( j in 2:tot_nodes ){
      count <- count + 1
      plot( density( mcmc_priors[[ k ]]$divt[[ nod_ind[1] ]], adj = 1 ),
            xlim = c(min_x,max_x),
            main = paste( "Compare ages for ", names(dup_dat)[i], " - ",
                          dup_dat[[ i ]][j],
                          sep = "" ))
    }
    dev.off()
  }
}

#> NOTE: It seems that the everything is fine! :) Braced nodes have the same
#> distributions!

