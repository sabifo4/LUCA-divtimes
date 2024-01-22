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
main_wd       <- gsub( pattern = "figs_tables/scripts/", replacement = "", x = wd )
wd_conc_cb    <- gsub( pattern = "figs_tables/scripts/", replacement = "02_MCMCtree/sum_files_post/", x = wd )
wd_conc_notcb <- gsub( pattern = "figs_tables/scripts/", replacement = "03_MCMCtree_nonbraced/sum_files_post/", x = wd )
wd_conc_fosscb <- gsub( pattern = "figs_tables/scripts/", replacement = "04_MCMCtree_fossbraced/sum_files_post/", x = wd )
wd_part_cb    <- gsub( pattern = "figs_tables/scripts/", replacement = "05_MCMCtree_part/sum_files_post/", x = wd )
wd_part_notcb <- gsub( pattern = "figs_tables/scripts/", replacement = "06_MCMCtree_part_nonbraced/sum_files_post/", x = wd )
wd_part_fosscb <- gsub( pattern = "figs_tables/scripts/", replacement = "07_MCMCtree_part_fossbraced/sum_files_post/", x = wd )
all_wds <- c( wd_conc_cb, 
              wd_conc_notcb,
              wd_conc_fosscb,
              wd_part_cb,
              wd_part_notcb,
              wd_part_fosscb )
source( "../../../src/Functions_plots.R" )

#### LOAD DATA ----
#-------------#
# START TASKS #
#-------------#
# Load all output csv tables previously generated
dir_names <- c( "conc_cb", "conc_notcb", "conc_fosscb",
                "part_cb", "part_notcb", "part_fosscb" )
all_divt  <- divt_csv( all_paths = all_wds, name_dirs = dir_names )
# Select nodes that are to be plotted:
nodes_2plot <- paste( "t_n", c( 247,
                                248, 368,
                                249, 438,
                                267, 455,
                                268, 456,
                                269, 457,
                                272, 332, 395, 460,
                                274, 333, 397, 462,
                                275, 334, 398, 463,
                                280, 312, 338, 402, 426, 468,
                                283, 314, 340, 404, 430, 471,
                                284, 315, 341, 405, 472,
                                287, 318, 344, 408, 431, 475,
                                292, 413, 480,
                                294, 482,
                                296, 484,
                                304, 369,
                                308, 422,
                                309, 423,
                                328, 390,
                                329, 391 ), sep = "" )
plots_per_fig <- c( 1, 2, 2, 2, 2, 2, 4, 4, 4, 6,
                    6,5, 6, 3, 2, 2, 2, 2, 2, 2, 2 )
names( nodes_2plot ) <- c( "ROOT",
                           c( "LUCA", "LUCA-dup" ),
                           c( "LACA", "LACA-dup" ),
                           c( "LAsCA", "LAsCA-dup" ),
                           c( "TG-EUKARYA-ARCH", "TG-EUKARYA-ARCH-DUP" ),
                           c( "LECA", "LECA-DUP" ),
                           c( "FUNGI", rep( "FUNGI-DUP", 3) ),
                           c( "METAZOA", rep( "METAZOA-DUP", 3 ) ),
                           c( "EUMETAZOA", rep( "EUMETAZOA-DUP", 3 ) ),
                           c( "ARCHAEPLASTIDA", rep( "ARCHAEPLASTIDA-DUP", 5 ) ),
                           c( "EMBRYOPHYTA", rep( "EMBRYOPHYTA-DUP", 5 ) ),
                           c( "EUDICOT-MONOCOT", rep( "EUDICOT-MONOCOT-DUP", 4 ) ),
                           c( "RHODOPHYTA", rep( "RHODOPHYTA-DUP", 5 ) ),
                           c( "FORAMINIFERA", rep( "FORAMINIFERA-DUP", 2 ) ),
                           c( "HEIMDAL", "HEIMDAL-DUP" ),
                           c( "ASGARD", "ASGARD-DUP" ),
                           c( "LBCA", "LBCA-DUP" ),
                           c( "TG-OXYPHOTOBACTERIA", "TG-OXYPHOTOBACTERIA-DUP" ),
                           c( "CG-OXYPHOTOBACTERIA", "CG-OXYPHOTOBACTERIA-DUP" ),
                           c( "LECA", "LECA-DUP" ),
                           c( "TG-EUKARYA-MITO", "TG-EUKARYA-MITO-DUP" ) )
# Get indexes
only_nums <- gsub( x = rownames( all_divt$GBM$conc_cb ),
                   pattern = "_[A-Z]..*", replacement = "" )


#### TABLES ---
#------------------#
# COMPARISON TABLE #
#------------------#
# Get the total number of columns, subtract one as we are not going to keep
# the fourth column with the calibrations (only a check when running 
# MCMC diagnostics)
cols_per_dataset <- length( colnames( all_divt$GBM$conc_cb ) )-1
short_GBM <- short_ILN <- matrix( 0, nrow = length( nodes_2plot ),
                                  ncol = cols_per_dataset*length(dir_names) )
rownames( short_GBM ) <- rownames( short_ILN ) <- 
  paste( nodes_2plot, "|", names( nodes_2plot ), sep = "" )
# Prepare a vector to create colnames based on number of datasets
conc_names <- vector( mode = "character", length = length(dir_names)*3 )
start <- end <- 0
for( i in 1:length( dir_names ) ){
  if( i == 1 ){
    start <- 1 
    end   <- start + 2
  }else{
    start <- end + 1
    end <- start + 2
  }
  conc_names[start:end] <- rep( dir_names[i], 3 )
}
# Get colnames ready
colnames( short_GBM ) <- colnames( short_ILN ) <- 
  paste( rep( c( "mean_t", "2.75%CI", "97.5%CI" ), length( dir_names ) ),
         conc_names, sep = "-" )
# Start saving the output results in the empty vectors for which memory has
# already been allocated in the corresponding cells
start <- end <- 0
for( i in 1:length( dir_names ) ){
  if( i == 1 ){
    start <- 1 
    end   <- cols_per_dataset
  }else{
    start <- end + 1
    end   <- start + cols_per_dataset - 1
  }
  cat( "round ", i, "start = ", start, " | end = ", end, "\n" )
  count <- 0
  for( j in nodes_2plot ){
    count <- count + 1
    tmp_ind <- which( only_nums %in% j )
    short_GBM[count,start:end] <- as.numeric( all_divt$GBM[[ i ]][tmp_ind,1:3] )
    short_ILN[count,start:end] <- as.numeric( all_divt$ILN[[ i ]][tmp_ind,1:3] )
  }
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

#### PLOTS ----
#----------------#
# START PLOTTING #
#----------------#
# Generate one plot per node indicated above -- each PDF file will have as
# many graphs as mirrored nodes!
comparison_plots( sum_obj = all_divt,
                  # Directory `plots` is created inside this directory
                  out_dir = paste( main_wd, "figs_tables/", sep = "" ),
                  plots_per_doc = plots_per_fig,
                  all_nodes = nodes_2plot, lab_nodes = only_nums,
                  data_perclock = 6,
                  # Same order as datasets inside object `all_divt`
                  x_labs = c( "GBM, conc, cb", "GBM, conc", "GBM, conc, fosscb",
                              "GBM, part, cb", "GBM, part", "GBM, part, fosscb",
                              "ILN, conc, cb", "ILN, conc", "ILN, conc, fosscb",
                              "ILN, part, cb", "ILN, part", "ILN, part, fosscb" ),
                  # Colours for those nodes that are cross-braced
                  # so it is easier to identify that the estimates are the same
                  # As the datasets in positions 1, 3, 5, 7 are cross-braced,
                  # these positions have different colours from black
                  points_col =  c( "red", "black", "orange",
                                   "pink", "black", "lightgreen",
                                   "blue", "black", "cyan",
                                   "purple", "black", "darkgrey" ),
                  sep_space = rep( c( 0.2, 0.1, 0.2, 0.1, 0.2, 0.1 ), 2 ),
                  suffix = "" )


