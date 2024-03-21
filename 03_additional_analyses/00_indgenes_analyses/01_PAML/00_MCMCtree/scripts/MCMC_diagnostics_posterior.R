#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# Load package needed to automatically find the path to the "scripts"
# directory
library( rstudioapi )
scripts_dir   <- gsub( pattern = "scripts..*", replacement = "scripts/",
                       x = getActiveDocumentContext()$path )
setwd( scripts_dir )
# Load the file with all the functions used throughout this script
source( file = "../../../../../src/Functions.R" )
# Run in-house function to set home directory and output directory for ESS
# and convergence tests
# NOTE: It will create a directory called `plots` and another called
# `ESS_and_chains_convergence` inside the `analyses` directory if you have
# not created them yet
home_dir      <- set_homedir()$home
outchecks_dir <- set_homedir()$ESS
# By now, set the working directory to `home_dir`
setwd( home_dir )


#--------------------------------#
# DEFINE USER'S GLOBAL VARIABLES #
#--------------------------------#
# First, we will define the global variables that match the settings in our 
# analysis.

# 1. Number of chains
num_chains <- 16

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
# set to 0.975 so the 97.5% and 2.5% quantiles. If you want to change this,
# however, just modify the value.
perc <- 0.975

# 5. Load a semicolon-separated file with info about calibrated nodes. Note that
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
# and `sigma2*` elements. Do not count the `lnL` value when looking at 
# `mcmc.txt` files generated when sampling from the posterior -- this is 
# automatically accounted for in the in-house R functions that you will 
# subsequently use. E.g., assuming an MCMC ran under a relaxed-clock model with  
# no partitions, we would see `mu` and `sigma2` columns. Therefore, the variable  
# would be set to `delcol = 2`. Please modify the values below according to your  
# the `mcmc.txt` file generated when sampling from the prior (`delcol_prior`) 
# dataset for and when sampling from the posterior (`delcol_posterior`).
delcol_post <- 2 # There are two column: mu, sigma2

# Path to the directory of each data alignment.
# In this case, we will have the path to the directory where
# the analyses when sampling from the prior took place (i.e., this is the path
# to the `CLK` directory that contains the subdirectories from `1` to `5` as
# we ran 5 chains in this example).
## POSTERIOR ##
dat <- c( "ATP", "EF", "Leu", "SRP", "Tyr" )
all_paths <- vector( mode = "character", length = length( dat )*2 )
count <- 0
for( i in 1:length( dat ) ){
  count <- count + 1
  all_paths[count]  <- paste( home_dir, "sum_analyses/01_posterior/GBM/", i, "/",
                           sep = "" )
  count <- count + 1
  all_paths[count]  <- paste( home_dir, "sum_analyses/01_posterior/ILN/", i, "/",
                          sep = "" )
}
#--------------#
# ANALYSE DATA #
#--------------#

## [[ POSTERIOR DATASET ]]
# Define object names
num_dirs      <- num_chains
delcol        <- delcol_post
path          <- all_paths
num_divt      <- num_divt
dataset       <- sort( c( paste( dat, "_GBM", sep = "" ),
                          paste( dat, "_ILN", sep = "" ) ) )
node_calib    <- calib_nodes[[ 1 ]] # same for all ind genes!
perc          <- perc
def_samples   <- def_samples
prior         <- FALSE
out_dat       <- all_paths
time_unit     <- 1000
out_file_pdf  <- sort( c( paste( "Convergence_plot_post_", dat, "_GBM",
                                 sep = "" ),
                          paste( "Convergence_plot_post_", dat, "_ILN",
                                 sep = "" ) ) )
out_title_pdf <- rep( c( "GBM","ILN" ), 5 )
th            <- 0.1
# Generate a list to keep all the returned objects
sum_post_QC          <- vector( "list", length( all_paths ) )
names( sum_post_QC ) <- dataset
# Now, run the `QC_conv` function in a loop
for( i in 1:length( dataset ) ){
  # NOTE: Only vectors with more than one element rely on `i`
  sum_post_QC[[ i ]] <- QC_conv( num_dirs = num_dirs, delcol = delcol, 
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
save( file = paste( home_dir, "out_RData/sum_post_QC.RData", sep = "" ),
      sum_post_QC )

# Print out sum stats in the script as a log

### ATP ------------
## [GBM -- ATP] 
sum_post_QC$ATP_GBM$ESS_results
# $median
# [1] 18787
# 
# $min
# [1] 17096
# 
# $max
# [1] 19205
# 
# $num_samples_for_ESS
# [1] 273536
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   134229   132235
# Min.   124994   105145
# Max.   136899   137777
# 
# $minRhat
# [1] 0.9999713
# 
# $maxRhat
# [1] 1.000117
# 
# $total_samples
# [1] 295945    245
length( sum_post_QC$ATP_GBM$not_conv_nodes ) # 0
sum_post_QC$ATP_GBM$num_chains # 16

## [ILN -- ATP]
sum_post_QC$ATP_ILN$ESS_results
# $median
# [1] 19530.5
# 
# $min
# [1] 14471
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 231536
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   115022   115100
# Min.   112440   107422
# Max.   116964   117104
# 
# $minRhat
# [1] 0.9999621
# 
# $maxRhat
# [1] 1.000135
# 
# $total_samples
# [1] 304994    245
length( sum_post_QC$ATP_ILN$not_conv_nodes ) # 0
sum_post_QC$ATP_ILN$num_chains # 16


### EF ------------
## [GBM -- EF]
sum_post_QC$EF_GBM$ESS_results
# $median
# [1] 18165
# 
# $min
# [1] 14613
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 233808
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   115006   114558
# Min.   103659    84408
# Max.   117753   118424
# 
# $minRhat
# [1] 0.9999657
# 
# $maxRhat
# [1] 1.000111
# 
# $total_samples
# [1] 289369    245
length( sum_post_QC$EF_GBM$not_conv_nodes ) # 0
sum_post_QC$EF_GBM$num_chains # 16

## [ILN -- EF]
sum_post_QC$EF_ILN$ESS_results
# $median
# [1] 19570
# 
# $min
# [1] 17857
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 285712
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   141687   141848
# Min.   136029   137343
# Max.   143126   143936
# 
# $minRhat
# [1] 0.9999551
# 
# $maxRhat
# [1] 1.000071
# 
# $total_samples
# [1] 306880    245
length( sum_post_QC$EF_ILN$not_conv_nodes ) # 0
sum_post_QC$EF_ILN$num_chains # 16

### Leu ------------
## [GBM -- Leu]
sum_post_QC$Leu_GBM$ESS_results
# $median
# [1] 19020
# 
# $min
# [1] 18017
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 288272
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   141492   139906
# Min.   132970   112723
# Max.   145057   145142
# 
# $minRhat
# [1] 0.9999604
# 
# $maxRhat
# [1] 1.000072
# 
# $total_samples
# [1] 305029    245
length( sum_post_QC$Leu_GBM$not_conv_nodes ) # 0
sum_post_QC$Leu_GBM$num_chains # 16

## [ILN -- Leu]
sum_post_QC$Leu_ILN$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 18116
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 289856
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   143352   142933
# Min.   139677   135362
# Max.   145765   145568
# 
# $minRhat
# [1] 0.999958
# 
# $maxRhat
# [1] 1.000085
# 
# $total_samples
# [1] 315278    245
length( sum_post_QC$Leu_ILN$not_conv_nodes ) # 0
sum_post_QC$Leu_ILN$num_chains # 16


### SRP ------------
## [GBM -- SRP]
sum_post_QC$SRP_GBM$ESS_results
# $median
# [1] 18455
# 
# $min
# [1] 17066
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 273056
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   134504   133824
# Min.   123178   115796
# Max.   136972   138054
# 
# $minRhat
# [1] 0.9999502
# 
# $maxRhat
# [1] 1.000112
# 
# $total_samples
# [1] 296763    245
length( sum_post_QC$SRP_GBM$not_conv_nodes ) # 0
sum_post_QC$SRP_GBM$num_chains # 16

## [ILN -- SRP]
sum_post_QC$SRP_ILN$ESS_results
# $median
# [1] 19034
# 
# $min
# [1] 15379
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 246064
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   122246   122466
# Min.   118732   119189
# Max.   123855   124071
# 
# $minRhat
# [1] 0.9999514
# 
# $maxRhat
# [1] 1.000115
# 
# $total_samples
# [1] 301147    245
length( sum_post_QC$SRP_ILN$not_conv_nodes ) # 0
sum_post_QC$SRP_ILN$num_chains # 16


### Tyr ------------
## [GBM -- Tyr]
sum_post_QC$Tyr_GBM$ESS_results
# $median
# [1] 19203.5
# 
# $min
# [1] 17159
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 274544
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   134867   133943
# Min.   119639   101071
# Max.   137726   138259
# 
# $minRhat
# [1] 0.9999669
# 
# $maxRhat
# [1] 1.00009
# 
# $total_samples
# [1] 301955    245
length( sum_post_QC$Tyr_GBM$not_conv_nodes ) # 0
sum_post_QC$Tyr_GBM$num_chains # 16

## [ILN -- Tyr]
sum_post_QC$Tyr_ILN$ESS_results
# $median
# [1] 19288.5
# 
# $min
# [1] 17369
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 277904
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.   138063   138028
# Min.   133835   131304
# Max.   139601   140689
# 
# $minRhat
# [1] 0.9999649
# 
# $maxRhat
# [1] 1.000076
# 
# $total_samples
# [1] 306000    245
length( sum_post_QC$Tyr_ILN$not_conv_nodes ) # 0
sum_post_QC$Tyr_ILN$num_chains # 16
