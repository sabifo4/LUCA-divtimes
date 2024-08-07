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
scripts_dir   <- gsub( pattern = "scripts..*", replacement = "scripts/",
                       x = getActiveDocumentContext()$path )
setwd( scripts_dir )
# Load the file with all the functions used throughout this script
source( file = "../../../../../src/Functions.R" )
# Run in-house function to set home directory and output directory for ESS
# and convergence tests
# NOTE: This function will create a directory called `plots` and another called
# `ESS_and_chains_convergence` inside the `analyses` directory if you have
# not created them yet
home_dir      <- set_homedir()$home
outchecks_dir <- set_homedir()$ESS
# By now, set the working directory to `home_dir`
setwd( home_dir )

#-------------------------------------------------------------#
# DEFINE GLOBAL VARIABLES -- modify according to your dataset #
#-------------------------------------------------------------#
# First, we will define the global variables that match the settings in our 
# analysis.

# 1. Number of chains
num_chains <- 6

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
def_samples <- 20000

# 4. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles (i.e., 95%CI). If you want to
# change this, however, just modify the value.
perc <- 0.975

# 5. Load a semicolon-separated file with info about calibrated nodes. Note that
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
ind_dat <- c( "LUCA_indgenes" ) # same calib file for all!
dat_ff  <- list.files( path = "calib_files", pattern = "Calibnodes.csv",
                       full.names = FALSE )
calib_nodes <- read_calib_f( main_dir = paste( home_dir, "calib_files/",
                                               sep = "" ),
                             f_names = dat_ff,
                             dat = ind_dat, head_avail = TRUE )

# 6. Number of columns in the `mcmc.txt` that are to be deleted as they do not 
# correspond to sample values for divergence times (i.e., the entries are not 
# names following the format `t_nX`). To figure out this number quickly, you 
# can open the `mcmc.txt` file, read the header, and count the number of `mu*`
# elements. Do not count the `lnL` value when looking at 
# `mcmc.txt` files generated when sampling from the posterior -- this is 
# automatically accounted for in the in-house R functions that you will 
# subsequently use. E.g., you expect to see as many `mu[0-9]` as alignment
# blocks you have in your sequence file! E.g., if you had two alignment blocks,
# you would speciy `delcol_prior <- 2`. Please modify the value/s below 
# (depending on having one or more datasets) according to the `mcmc.txt` file
# generated when sampling from the prior (`delcol_prior`)
##> NOTE: If you ran `MCMCtree` with `clock = 2` or `clock = 3` when
##> sampling from the prior, you will also need to count the `sigma2*`
##> columns! We ran `clock = 1` so that the analyses ran quicker, and thus
##> we only have `mu*` columns.
delcol_prior <- 1 # There is one column: mu

# 7. Path to the directory where the subdirectories where each chain ran for 
# each dataset are saved In this case, we will have the path to the directory 
# where the analyses when sampling from the prior took place (i.e., directory
# that contains the subdirectories from `1` to `n`, where `n` is the number
# of chains we ran).
dat <- c( "ATP", "EF", "Leu", "SRP", "Tyr" )
path_prior <- vector( mode = "character", length = length( dat ) )
for( i in 1:length( dat ) ){
  path_prior[i]  <- paste( home_dir, "sum_analyses/00_prior/CLK/", i, "/",
                           sep = "" )
}

#--------------#
# ANALYSE DATA #
#--------------#
# Define object names
num_dirs      <- num_chains
delcol        <- delcol_prior
path          <- path_prior
num_divt      <- num_divt
node_calib    <- calib_nodes[[ 1 ]] # same for all ind genes!
dataset       <- dat
perc          <- perc
def_samples   <- def_samples
prior         <- TRUE
out_dat       <- path_prior
time_unit     <- 1000
out_file_pdf  <- paste( "Convergence_plot_prior_", dat, sep = "" )
out_title_pdf <- paste( "CLK_", dat, sep = "" )
th            <- 0.1
# Generate a list to keep all the returned objects
sum_prior_QC          <- vector( "list", length( path_prior ) )
names( sum_prior_QC ) <- dataset
# Now, run the `QC_conv` function in a loop
for( i in 1:length( dataset ) ){
  # NOTE: Only vectors with more than one element rely on `i`
  sum_prior_QC[[ i ]] <- QC_conv( num_dirs = num_dirs, delcol = delcol, 
                                  path = path[i], num_divt = num_divt, 
                                  node_calib = node_calib,
                                  dataset = dataset[i],
                                  perc = perc, def_samples = def_samples,
                                  prior = prior, out_dat = out_dat[i],
                                  time_unit = time_unit, 
                                  out_file_pdf = out_file_pdf[i],
                                  out_title_pdf = out_title_pdf[i], 
                                  th = th,
                                  outchecks_dir = outchecks_dir )
}

# Save object!
if( ! dir.exists( paste( home_dir, "out_RData", sep = "" ) ) ){
  dir.create( paste( home_dir, "out_RData", sep = "" ) ) 
}
save( file = paste( home_dir, "out_RData/sum_prior_QC.Rdata", sep = "" ),
      sum_prior_QC )

# Print out sum stats in the script as a log

## [[ CLK - ATP ]] ####
sum_prior_QC$ATP$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59282    59681
# Min.    57378    57785
# Max.    60490    61641
# 
# $minRhat
# [1] 0.9999518
# 
# $maxRhat
# [1] 1.000157
# 
# $total_samples
# [1] 120006    245
length( sum_prior_QC$ATP$not_conv_nodes ) # 0
sum_prior_QC$ATP$num_chains # 6

# [[ CLK - EF ]] ####
sum_prior_QC$EF$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59352    59680
# Min.    56084    56557
# Max.    60702    60666
# 
# $minRhat
# [1] 0.9999583
# 
# $maxRhat
# [1] 1.000157
# 
# $total_samples
# [1] 120006    245
length( sum_prior_QC$EF$not_conv_nodes ) # 0
sum_prior_QC$EF$num_chains # 6

# [[ CLK - Leu ]] ####
sum_prior_QC$Leu$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59456    59886
# Min.    56286    56447
# Max.    60735    61475
# 
# $minRhat
# [1] 0.999945
# 
# $maxRhat
# [1] 1.000144
# 
# $total_samples
# [1] 120006    245
length( sum_prior_QC$Leu$not_conv_nodes ) # 0
sum_prior_QC$Leu$num_chains # 6


# [[ CLK - SRP ]] ####
sum_prior_QC$SRP$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59400    59801
# Min.    56285    57948
# Max.    60524    61090
# 
# $minRhat
# [1] 0.9999416
# 
# $maxRhat
# [1] 1.00021
# 
# $total_samples
# [1] 120006    245
length( sum_prior_QC$SRP$not_conv_nodes ) # 0
sum_prior_QC$SRP$num_chains # 6

# [[ CLK - Tyr ]] ####
sum_prior_QC$Tyr$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59261    59864
# Min.    55382    56935
# Max.    60613    61001
# 
# $minRhat
# [1] 0.9999512
# 
# $maxRhat
# [1] 1.000205
# 
# $total_samples
# [1] 120006    245
length( sum_prior_QC$Tyr$not_conv_nodes ) # 0
sum_prior_QC$Tyr$num_chains # 6
