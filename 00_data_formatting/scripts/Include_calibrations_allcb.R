#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd. Note that the wd will be the same directory where this 
# script is saved.
setwd( wd )

#----------------------------#
# DEFINE GLOBAL VARS BY USER #
#----------------------------#
# Name of output calibrated tree file ready to be used by `MCMCtree`.
# Note that the file name will have the following format
# "<your_selected_out_name>_calib_MCMCtree.tree".
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files. This file needs to be
# already in PHYLIP format. Please follow the same format as used in the 
# example tree file provided.
out_name <- c( "LUCAdup_246sp_allcb" )

# Path to your input tree with calibrations. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file. You need to include the
# flags within square brackets (e.g., [Mammalia]) and write them on the node
# that is to be calibrated.
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
path_tree <- c( "../00_raw_data/trees/01_calibrations/LUCAdup_allcb_calibnames.tree" )

# Path to your input text file that allows you to match the flags you have 
# used to label the nodes that are to be calibrated with the calibration you 
# want to use in `MCMCtree` format. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file.
# The format of this file meets the following criteria:
#
#   - No header.
#   - One row per calibration.
#   - No spaces at all.
#   - The pipe character, "|", is used to separate the flag you are using on
#     the node you want to calibrate (left) from the `MCMCtree` calibration
#     that you will use to calibrate the node (right).
#   - The `MCMCtree` calibrations are written in `MCMCtree` format and with 
#     single quotation marks. No spaces.
#   - No spaces when writing the name of your flag either.
# 
# E.g.: row in this text file to calibrate node "Mammalia". The flag used in the
# tree to locate the node "Mammalia" is "Mammalia" and the `MCMCtree`
# calibration is a soft-bound calibration:
#
# ```
# Mammalia|'B(1.649,2.51254)'
# ```
#`
# If in doubt, please follow the same format as used in the example text 
# file provided.
path_textconv <- c( "Calib_converter_allcb.txt" )

#---------------------------------#
# READ TREE AND CALIBRATIONS FILE #
#---------------------------------#
# Read tree and get phylip header
# NOTE: Make always sure that there is at least one blank line at the 
# end of the tree file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
tt_name       <- path_tree
tt            <- readLines( tt_name )
phylip.header <- tt[1]
tt            <- tt2 <- tt[2]

# Read file converter that will be used to rename calibrations
calibrations <- read.table( file = path_textconv,
                            stringsAsFactors = F, sep = "|", blank.lines.skip = T )
colnames( calibrations ) <- c( "name", "MCMCtree calib" )

#--------------------------------#
# REPLACE TAGS WITH CALIBRATIONS #
#--------------------------------#
# Replace calibration names with corresponding calibration
for( j in 1:dim( calibrations )[1] ){
  # Conditional is used so that the single quotation marks are only kept 
  # in the upper-bound calibration for the root. Inequality calibrations
  # do not require single quotation marks
  tmp_calib <- gsub( x = calibrations[j,2], pattern = "\\(..*",
                     replacement = "" )
  tmp_calib2 <- gsub( x = calibrations[j,2], pattern = "[0-9]*",
                     replacement = "" )
  if( tmp_calib == "U" ){
    tt <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                x = tt,
                replacement = paste( "'", calibrations[j,2], "'", sep = "" ) )
  }else{
    tt <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                x = tt,
                replacement = paste( calibrations[j,2], sep = "" ) )
  }
  # Copy to visualise in FigTree
  reps <- gsub( x = gsub( x = gsub( x = gsub( x = gsub( x = calibrations[j,2],
                                                        pattern = "\\{",
                                                        replacement = "(" ),
                                              pattern = "\\}",
                                              replacement = ")" ), 
                                    pattern = "\\[|\\]", replacement = "" ),
                          pattern = "\\#", replacement = "flag" ),
                pattern = " ", replacement = "-" )
  if( tmp_calib2 == '#' ){
    reps <- gsub( x = gsub( x = reps, pattern = "\\#", replacement = "flag" ),
                  pattern = "\\]", replacement = "" )
    tt2 <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                 x = tt2,
                 replacement = paste0( "'", reps, "-", calibrations[j,1], "'", 
                 collapse = "" ) )
  }else{
    tt2 <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                 x = tt2,
                 replacement = paste0( "'", reps, "-", calibrations[j,1], "'",
                                       collapse = "" ) )
  }
}


#-------------------------------#
# WRITE CALIBRATED TREE IN FILE #
#-------------------------------#
out_dir <- "../01_inp_data/"
if( ! dir.exists( "../01_inp_data/" ) ){
  dir.create( "../01_inp_data/" )
}
write( x = phylip.header, file = paste( out_dir, out_name, "_calib_MCMCtree.tree", sep = "" ) )
write( x = tt, file = paste( out_dir, out_name, "_calib_MCMCtree.tree", sep = "" ),
       append = TRUE )
write( x = phylip.header, file = paste( "../00_raw_data/trees/01_calibrations/",
                                        out_name,
                                        "_outR_fordisplay_calib_MCMCtree.tree",
                                        sep = "" ) )
write( x = tt2, file = paste( "../00_raw_data/trees/01_calibrations/", out_name,
                              "_outR_fordisplay_calib_MCMCtree.tree", sep = "" ),
       append = TRUE )
  




