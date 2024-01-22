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

#-------#
# TASKS #
#-------#
# Scan taxa
all_taxa <- readLines( con = "../00_raw_data/alignment/all_taxa.txt")
# Scan alignments
all_f <- list.files( path = "../00_raw_data/alignment/partitioned/", 
                     pattern = "fas$" )
for( j in all_f ){
  cat( "[[ Analysing dataset ", j, " ]]\n")
  conc <- gsub( x = readLines( con = paste( "../00_raw_data/alignment/partitioned/",
                                            j, sep = "" ) ),
                pattern = ">", replacement = "" )
  # Find which taxa are present in alignment
  avail_taxa <- which( conc %in% all_taxa )
  ind_taxa   <- c( avail_taxa )
  miss_taxa  <- which( ! all_taxa %in% conc )
  len_seq    <- stringr::str_count( conc[2], "" )
  gap_seq    <- paste( rep( "-", len_seq ), collapse = "" )
  # Write alignment
  tmp_name <- gsub( x = j, pattern = "_..*", replacement = "_filt" )
  count <- 0
  for( i in ind_taxa ){
    count <- count + 1
    if( count == 1 ){
      write( x = paste( ">", conc[i], sep = "" ),
             file = paste( "../00_raw_data/alignment/partitioned/", tmp_name,
                           ".fasta", sep = "" ), )
      write( x = conc[i+1],
             file = paste( "../00_raw_data/alignment/partitioned/", tmp_name,
                           ".fasta", sep = ""  ), append = TRUE )
    }else{
      write( x = paste( ">", conc[i], sep = "" ),
             file = paste( "../00_raw_data/alignment/partitioned/", tmp_name,
                           ".fasta", sep = ""  ),
             append = TRUE )
      write( x = conc[i+1],
             file = paste( "../00_raw_data/alignment/partitioned/", tmp_name,
                           ".fasta", sep = ""  ),
             append = TRUE )
    }
  }
  # Add now missing taxa
  for( i in miss_taxa ){
    write( x = paste( ">", all_taxa[i], sep = "" ),
           file = paste( "../00_raw_data/alignment/partitioned/", tmp_name,
                         ".fasta", sep = ""  ),
           append = TRUE )
    write( x = gap_seq,
           file = paste( "../00_raw_data/alignment/partitioned/", tmp_name,
                         ".fasta", sep = ""  ),
           append = TRUE )
    
  }
}






