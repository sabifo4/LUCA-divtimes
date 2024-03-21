#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#--------------#
# LOAD OBJECTS #
#--------------#
# Load tree 
raw_tt <- ape::read.tree( file = "../../00_data_formatting/00_raw_data/trees/00_IQTREE//LUCAdup_topo_bl_rooted.tree" )

#-------#
# TASKS #
#-------#
# 1. Find tree height. You can use the function `phytools::nodeHeights` to calculate 
#    all the heights of the tree. Then, we can extract the maximum height calculated,
#    which will correspond to the length from the root to the highest tip.
tree_height <- max( phytools::nodeHeights( raw_tt ) ) #3.573201

# 2. Get an estimate of the calibration set for the root to have the 
#    time for the speciation event at the root, what we will use 
#    to estimate the mean evolutionary rate later. As we have an upper-bound
#    calibration (maximum age of 4.520), we will use such age as an approximate
#    values for the root age in 1,000Ma = 1 Ga.
root_age <- 4.520 # Ga

# 3. Estimate mean rate based on the two different time units
#
#    tree_height = mean_rate x root_age --> mean_rate = tree_height / root_age
mean_rate <- tree_height / root_age # 0.7905311 subst/site per time unit

# If we want to know the mean rate in subst/site/year, we apply the time unit. We
# We should get the same estimate regardless the time unit used:
#
# Time unit 1Ga (1e+09y): 0.7905311 subst/site/1e+09 = 7.91e-10 subst/site/year

# 4. Now, we can build the gamma distribution given that we now have an estimate
#    for the mean rate. We will also use `alpha = 2` as we will start with a 
#    vague distribution. Nevertheless, if you were very sure about the mean 
#    rate, you could build a more constraint prior.
#
#    mean_Gamma = mean_rate = alpha / beta --> beta = alpha / mean_rate
alpha <- 2
beta  <- alpha/mean_rate # 2.529945 ~ 2.5

# We can plot this distribution
curve( dgamma( x, shape = 2, rate = beta ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,2.5) " ), 
        col = "black", lty = 1, box.lty = 2 )

# 5. Plot the gamma distributions
if ( ! dir.exists( "out_RData" ) ){
  dir.create( "out_RData" )
}
pdf( file = "out_RData/gamma_dists.pdf", paper = "a4" )
curve( dgamma( x, alpha, rate = beta ), from = 0, to = 2, col = "black" )
legend( "topright", legend = c( "G(2,2.5) " ), 
        col = "black", lty = 1, box.lty = 2 )
dev.off()

