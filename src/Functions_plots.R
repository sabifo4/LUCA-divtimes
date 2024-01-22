# Function to create objects to later be plotted using adapted functions 
# from the `MCMCtreeR` R package
#
# Arguments
#
# abs_path      Characater, length `n`. This object will have as many entries as 
#               datasets analysed. Absolute path required or, if not, relative
#               path from the working directory specified at the beginning of
#               the script.
# name_entries  Character, length `n`. This object will have as many entries as 
#               datasets analysed. E.g., `GBM_conc_cb`
# tree_pattern  Character. Pattern to match the tree with mean divergence times
#               that is to be used for plotting. E.g., "95CI", "95HPD", etc.
create_plot_obj <- function( abs_path, name_entries, tree_pattern )
{
  # Stop if condition is not met
  if( length( abs_path ) != length( name_entries ) ){
    stop( "The length of \"abs_path\", and \"name_entries\" are not equal\n" )
  }
  # Create TMP objects for all datasts (steps 1-2)
  phy_all <- edge_all <- mcmc_all <- vector( mode = "list",
                                             length = length( name_entries ) )
  names( phy_all ) <- name_entries
  # Create `node_ages_all` object (steps 3-8)
  node_ages_all <- vector( mode = "list", length = length( name_entries ) )
  names( node_ages_all ) <- name_entries
  # Run all steps!
  for( i in 1:length( abs_path ) ){
    #>> Run steps 1-2 to create `phy`, `mcmc`, and `edge` objects  
    ## 1. Load trees and mcmc files for all datasets
    cat( "Parsing files for dataset ", i, " ...\n" )
    cat( "[[ PATH: ", abs_path[i], " ]]\n" )
    tt   <- list.files( abs_path[i], pattern = tree_pattern )
    cat( "----> Parsing tree \"", tt, "\"... \n" )
    phy  <- MCMCtreeR::readMCMCtree( paste( abs_path[i], "/", tt,
                                            sep = "" ),
                                     from.file = TRUE )
    cat( "----> Parsing MCMC file ...\n" )
    mcmc <- read.table( paste( abs_path[i], "/mcmc.txt", sep = "" ),
                        sep = "\t", header = T, stringsAsFactors = F )
    ## 2. Create objects in APE format
    phy_all[[ i ]]  <- phy$apePhy
    edge_all[[ i ]] <- phy$apePhy$edge
    mcmc_all[[ i ]] <- mcmc
    ## 3. Extract ages with node age posteriors from column 2
    cat( "----> Extracting node information ...\n" )
    mcmc_node_ages <- mcmc_all[[ i ]][, 2:Ntip( phy_all[[ i ]] )] 
    all_nodes      <- as.numeric( gsub( "t_n", "",
                                        colnames( mcmc_node_ages ) ) )
    ## 4. Create a vector of names for each list element as internal nodes
    ##    Use APE tree, using phy$edge object  
    node_ages_names <- c( Ntip( phy_all[[ i ]] ) + 1,
                                 edge_all[[ i ]][ which( edge_all[[ i ]][,2] > Ntip( phy_all[[ i ]] ) ), 2 ] )
    ## 5. Find where each posterior node age appears in APE edge object
    cat( "----> Matching nodes with their corresponding posterior distributions ...\n" )
    match_nodes <- match( all_nodes, as.numeric( node_ages_names ) )
    ## 6. Create a list extracting the info from the mcmc
    ##    chain in APE node order
    cat( "----> Create a list with MCMC info and corresponding nodes ...\n" )
    node_ages <- lapply( match_nodes, function( uu ) mcmc_node_ages[, uu ] )
    ## 7. Name each element in list
    cat( "----> Generate output file for this dataset ...\n\n" )
    names( node_ages ) <- node_ages_names
    ## 8. Save object
    node_ages_all[[ i ]] <- node_ages
  }
  
  # Return objects
  cat( "Tasks finished!\n Returning objects\n\n" )
  return( list( phy = phy_all, node_ages = node_ages_all ) )
          
}


# Function to read all csv output files with sum stats about divtimes
#
# Arguments
# 
# all_paths Character, length `n`. As many entries as datasets.
# name_dirs Character, length `n`. As many entries as datasets. You can use
#           abbreviations to identify each dataset. E.g., "GBM_notcb".
divt_csv <- function( all_paths, name_dirs )
{
  # Create objects for each dataset that will be required for plotting
  GBM_all <- ILN_all <- vector( mode = "list", length = length( all_paths ) )
  names( GBM_all ) <- names( ILN_all ) <- name_dirs
  for( i in 1:length( all_paths ) ){
    cat( "Parsing dataset ", i, "...\n" )
    cat( "[[ PATH: ", all_paths[i], "]]\n" )
    all_csv <- list.files( all_paths[i], pattern = "all_mean_est.tsv" )
    ind_GBM <- grep( pattern = "GBM", x = all_csv )
    GBM_csv <- all_csv[ind_GBM]
    ind_ILN <- grep( pattern = "ILN", x = all_csv )
    ILN_csv <- all_csv[ind_ILN]
    ind_GBM_filt <- grep( pattern = "FILT", x = GBM_csv )
    ind_ILN_filt <- grep( pattern = "FILT", x = ILN_csv )
    # Read GBM file
    if( length( ind_GBM_filt ) > 0 ){
      cat( "   ---> This is a filtered dataset!\n" )
      GBM_f <- read.table( file = paste( all_paths[i], GBM_csv[ind_GBM_filt],
                                         sep = "" ), 
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t" )
    }else{
      GBM_f <- read.table( file = paste( all_paths[i], GBM_csv, sep = "" ), 
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t" )
    }
    # Read ILN file
    if( length( ind_ILN_filt ) > 0 ){
      ILN_f <- read.table( file = paste( all_paths[i], ILN_csv[ind_ILN_filt],
                                         sep = "" ), 
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t" )
    }else{
      ILN_f <- read.table( file = paste( all_paths[i], ILN_csv, sep = "" ), 
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t" )
    }
    # Save objects
    GBM_all[[ i ]] <- GBM_f
    ILN_all[[ i ]] <- ILN_f
    
  }
  
  # Return objects
  return( list( GBM = GBM_all, ILN = ILN_all ) )
}


# Function to plot divtimes and CIs for all datasets, one after the other.
# Useful plots for comparison.
# Those cross-braced nodes have the same colours!
#
# Arguments
#
# sum_obj        List. Object with all the summary stats generated prior to
#                running this function. There will be as many entries as
#                analysis have been run. E.g., see object `all_divt` in main
#                scripts for reference.
# out_dir        Character. Path to the output directory where a new directory
#                called `plots` will be generated. All the output graphs will
#                be stored there. The output file names will start with
#                `Divtimes_` and will be followed by the name given to the node
#                being plotted. E.g., see `names(nodes_2plot)` in main script
#                for reference.
# plots_per_doc  Integer. Number of plots for each specific divergence time that
#                is to be plotted. Useful when there are mirrored nodes and
#                they are to be plotted in the same page to ease comparison.
#                E.g., see object `plots_per_fig` in main script for reference.
# all_nodes      Character. Object generated prior to running this funcion.
#                This object should consist of the node labels in `MCMCtree`
#                format and, in addition, names should have been given to
#                each entry. E.g.: see object `nodes_2plot` in main script
#                for reference.
# lab_nodes      Character. Same as `all_nodes` but without names. This object
#                is generated differently. See object `only_nums` in main
#                script for reference
# data_perclock  Integer. Number of analyses carried out under each clock.
#                In this study, it was four: conc+cb, conc, part+cb, part.
# x_labs         Character. Abbreviations used to identify the datasets. Same
#                order is expected as the datasets found in `sum_obj`.
# points_col     Character. Colours chosen for the points plotted for each
#                dataset. It is suggested that, if cross-bracing used, you
#                use different colours from black to identify them in the
#                plots generated for each mirrored node.
comparison_plots <- function( sum_obj, out_dir, plots_per_doc, all_nodes, lab_nodes,
                              data_perclock = 4, x_labs,
                              points_col =  c( "red", "black", "pink", "black",
                                               "blue", "black", "purple", "black" ),
                              sep_space = rep( c( 0.2, 0.1, 0.2, 0.1 ), 2 ),
                              suffix = "" )
{
  # Plot mean age and CIs
  start <- end <- 0
  if( ! dir.exists( paste( out_dir, "/plots", sep = "" ) ) ){
    dir.create( paste( out_dir, "/plots", sep = "" ) ) 
  }
  for( k in 1:length( plots_per_doc ) ){
    # Go through `all_nodes` and plot them
    if( k == 1 ){
      start <- 1 
      end <- plots_per_doc[k]
    }else{
      start <- end + 1
      end   <- start + plots_per_doc[k] - 1
    }
    ##> START: Useful for debugging, do not delete
    #cat( "Round ", k, "start = ", start, "end = ", end, "\n" )
    ##> END
    plot_ns <- all_nodes[start:end]
    tmp_ind_lab <- which( lab_nodes %in% plot_ns )
    cat( "[[ Plotting results for node", names(plot_ns)[1], " ]]\n")
    cat( "---> Number of nodes = ", plots_per_doc[k], "\n" )
    cat( "---> Output directory:\n",
         paste( out_dir, "/plots/divtimes_", names(plot_ns)[1], suffix,
                ".pdf",
                sep = "" ), "\n\n" )
    # Determine num of plots per fig
    if( plots_per_doc[k] == 1 ){
      #par( mar = c(8, 8, 4, 2) )
      pdf( paste( out_dir, "/plots/divtimes_", names(plot_ns)[1], suffix,
                  ".pdf", sep = "" ), width=15, height=10 )
    }else if( plots_per_doc[k] == 2 ){
      #par( mar = c(8, 8, 4, 2), mfrow = c( 1, 2 ) )
      pdf( paste( out_dir, "/plots/divtimes_", names(plot_ns)[1], suffix,
                  ".pdf", sep = "" ), width=15, height=10 )
      par( mfrow = c( 1, 2 ) )
    }else if( plots_per_doc[k] == 3 ){
      #par( mar = c(8, 8, 4, 2), mfrow = c( 2, 2 ) )
      pdf( paste( out_dir, "/plots/divtimes_", names(plot_ns)[1], suffix,
                  ".pdf", sep = "" ), width=22.5, height=10 )
      par( mfrow = c( 1, 3 ) )
    }else if( plots_per_doc[k] == 4 ){
      #par( mar = c(8, 8, 4, 2), mfrow = c( 2, 2 ) )
      pdf( paste( out_dir, "/plots/divtimes_", names(plot_ns)[1], suffix,
                  ".pdf", sep = "" ), width=15, height=10 )
      par( mfrow = c( 2, 2 ) )
    }else if( plots_per_doc[k] == 5 || plots_per_doc[k] == 6 ){
      #par( mar = c(8, 8, 4, 2), mfrow = c( 2, 3 ) )
      pdf( paste( out_dir, "/plots/divtimes_", names(plot_ns)[1], suffix,
                  ".pdf", sep = "" ), width=22.5, height=10 )
      par( mfrow = c( 2, 3 ) )
    }
    count <- 0
    # Create empty vectors according to how many datasets per clock
    # are to be compared
    y_CImax_GBM <- y_CImin_GBM <- y_meant_GBM <- 
      y_CImax_ILN <- y_CImin_ILN <- y_meant_ILN <- 
      vector( mode = "numeric", length = data_perclock )
    tmpy_CImax_GBM <- tmpy_CImin_GBM <- tmpy_CImax_ILN <- tmpy_CImin_ILN <-
      vector( mode = "numeric", length = data_perclock*length(tmp_ind_lab) )
    st <- nd <- 0
    for( j in 1:data_perclock ){
      if( j == 1 ){
        st <- 1
        nd <- length( tmp_ind_lab )
        tmpy_CImax_GBM[st:nd] <- sum_obj$GBM[[ j ]][tmp_ind_lab,3]
        tmpy_CImin_GBM[st:nd] <- sum_obj$GBM[[ j ]][tmp_ind_lab,2]
        tmpy_CImax_ILN[st:nd] <- sum_obj$ILN[[ j ]][tmp_ind_lab,3]
        tmpy_CImin_ILN[st:nd] <- sum_obj$ILN[[ j ]][tmp_ind_lab,2]
      }else{
        st <- nd + 1
        nd <- st + length( tmp_ind_lab ) - 1
        tmpy_CImax_GBM[st:nd] <- sum_obj$GBM[[ j ]][tmp_ind_lab,3]
        tmpy_CImin_GBM[st:nd] <- sum_obj$GBM[[ j ]][tmp_ind_lab,2]
        tmpy_CImax_ILN[st:nd] <- sum_obj$ILN[[ j ]][tmp_ind_lab,3]
        tmpy_CImin_ILN[st:nd] <- sum_obj$ILN[[ j ]][tmp_ind_lab,2]
      }
      
    }
    # Get min and max y values
    min_y <- min( c( tmpy_CImin_GBM, tmpy_CImin_ILN ) ) - 10
    max_y <- max( c( tmpy_CImax_GBM, tmpy_CImax_ILN ) ) + 10
    for( i in tmp_ind_lab ){
      count <- count + 1
      # Define values to be plotted
      for( j in 1:data_perclock ){
        y_meant_GBM[j] <- sum_obj$GBM[[ j ]][i,1]
        y_meant_ILN[j] <- sum_obj$ILN[[ j ]][i,1]
        y_CImax_GBM[j] <- sum_obj$GBM[[ j ]][i,3]
        y_CImax_ILN[j] <- sum_obj$ILN[[ j ]][i,3]
        y_CImin_GBM[j] <- sum_obj$GBM[[ j ]][i,2]
        y_CImin_ILN[j] <- sum_obj$ILN[[ j ]][i,2]
      }
      # Get x values, leaving one empty slot between the datasets analysed
      # under each clock model
      x_vals <- c( 1:data_perclock, c(data_perclock+2):c(data_perclock*2+1) )
      # Plot results for node "i" for all GBM analyses
      plot( x = x_vals, 
            y = c( y_meant_GBM, y_meant_ILN ),
            pch = 16, ylim = c( min_y, max_y ),
            xaxt = "n", xlab = "",
            ylab = "Estimated mean divergence times",
            col = points_col )
      title( main = paste( "Node ", plot_ns[count], 
                           " | ", names( plot_ns )[count], sep = "" ) )
      axis( side = 1, at = x_vals,
            labels = FALSE, las = 2 )
      text( x = x_vals+c( sep_space ), # offset calculated by eye
            labels = x_labs,
            ## Rotate the labels by 25 degrees.
            par("usr")[3]-3, 
            srt = 25, adj = 1, xpd = TRUE )
      # Now, plot the CIs with dashed lines!
      lines( x = x_vals[1:c(length(x_vals)/2)], y = y_meant_GBM, lty = 1 )
      lines( x = x_vals[1:c(length(x_vals)/2)], y = y_CImax_GBM, lty = 2 )
      lines( x = x_vals[1:c(length(x_vals)/2)], y = y_CImin_GBM, lty = 2 )
      lines( x = x_vals[c(data_perclock+1):length(x_vals)], y = y_meant_ILN, lty = 1 )
      lines( x = x_vals[c(data_perclock+1):length(x_vals)], y = y_CImax_ILN, lty = 2 )
      lines( x = x_vals[c(data_perclock+1):length(x_vals)], y = y_CImin_ILN, lty = 2 )
      # For 97.5%CI, use triangle facing down
      points( x = x_vals, y = c( y_CImax_GBM, y_CImax_ILN ), pch = 25,
              cex = 0.7,
              col = points_col )
      # For 2.75%CI, use triangle facing up
      points( x = x_vals, y = c( y_CImin_GBM, y_CImin_ILN ), pch = 2,
              cex = 0.7,
              col = points_col )
    }
    # Close graphics
    dev.off()
    
  }
}


