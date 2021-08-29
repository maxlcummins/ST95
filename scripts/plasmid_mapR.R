#' A function for mapping plasmid blast hits for datavis and calculating discontiguous coverage.
#
#' Read the documents on the github page https://github.com/maxlcummins/plasmid_mapR for more info
#' @param path_to_abricate Filename of your abricate input
#' @param plasmid_reference_name Reference plasmid name. Determines name of your output files
#' @param output_directory Path to output directory. Directory should already exist.
#' @param min_hit_id Minimum nucleotide ID (also as a percentage, i.e. 90 = 90%)
#' @param min_hit_length Minimum hit length (as a percentage, i.e. 0.5 = 0.5%)
#' @param writecsv TRUE/FALSE: Write/Dont write to csv. Default = FALSE
#' @param path_to_tree Optional - Leave blank if you dont want to filter. Path to a tree for filtering the dataset based on names in the tree.
#' @param tree_reference_name Optional - Leave blank if you dont want to filter. If you are filtering your dataframe based on a tree and this tree is reference based, you can provide a reference name.



plasmid_mapR <-
  function(path_to_abricate,
           plasmid_reference_name,
           output_directory,
           min_hit_id = 90,
           min_hit_length = 0.5,
           writecsv = FALSE,
           path_to_tree = FALSE,
           tree_reference_name = "Placeholder") {
    
    require(pheatmap)
    require(ggplot2)
    require(tidyverse)
    require(magrittr)
    require(reshape2)
    require(tidytree)
    require(ggtree)
    
    #### Read in, process and subset abricate data ####
    #Read in the abricate genotype data sheet
    #(small number of rows for colname reassignment)
    #This is to reduce memory requirements
    abricate_hits <-
      read_delim(
        path_to_abricate,
        "\t",
        escape_double = FALSE,
        trim_ws = TRUE,
        n_max = 10
      )
    
    #Colname reassignment
    colnames(abricate_hits)[c(1, 10:11)] <-
      c("name", "perc_coverage", "perc_identity")
    
    #Extract column names for later reassignment
    abricate_hits_colnames <- colnames(abricate_hits)
    
    #Re-read in PAI abricate genotype data sheet
    abricate_hits <-
      read_delim(
        path_to_abricate,
        "\t",
        escape_double = FALSE,
        trim_ws = TRUE,
        col_names = FALSE,
        skip = 1
      )
    
    #Remove cases where there are multiple headers from concatenation of abricate reports
    abricate_hits <- abricate_hits %>% filter(X2 != "SEQUENCE")
    
    #Colname reassignment
    colnames(abricate_hits) <- abricate_hits_colnames
    
    #Convert percent coverage and identity to numeric type to allow filtering
    abricate_hits$perc_coverage <-
      as.numeric(abricate_hits$perc_coverage)
    abricate_hits$perc_identity <-
      as.numeric(abricate_hits$perc_identity)
    
    #Filter to perc_identity > 95%
    #abricate_hits <-
    abricate_hits <-
      abricate_hits %>% filter(perc_identity > min_hit_id)
    abricate_hits <-
      abricate_hits %>% filter(perc_coverage > min_hit_length)
    
    #Trim excess characters the assembly names and reassign this to rownames
    abricate_hits$name <- gsub("\\..*", "", abricate_hits$name)
    
    if (path_to_tree != FALSE) {
      #Read in the tree file
      tree <-
        read.tree(file = path_to_tree)
      
      #trim the names of the assemblies in the tree tip labels
      tree$tip.label <- gsub("\\..*", "", tree$tip.label)
      
      #Replace Reference with refname
      if(tree_reference_name != "Placeholder"){
        tree$tip.label <- gsub("Reference", tree_reference_name, tree$tip.label)
      }
      
      
      #Subset the hits to strains within the tree (this saves memory)
      abricate_hits <-
        abricate_hits %>% filter(name %in% tree$tip.label)
    }
    
    #Extract from coverage column the coverage range and the length of the reference
    abricate_hits$coverage_range <-
      gsub("\\/.*", "", abricate_hits$COVERAGE)
    abricate_hits$ref_length <-
      gsub(".*\\/", "", abricate_hits$COVERAGE)
    
    #Save length of plasmid reference as a variable for later use
    ref_length <- as.numeric(unique(abricate_hits$ref_length))
    
    #Replace the '-' with a ':' in the coverage
    abricate_hits$coverage_range <-
      gsub("-", ":", abricate_hits$coverage_range)
    
    #Create a column for start and end coordinates of hits
    abricate_hits$end <-
      gsub("[0-9]+:", "", abricate_hits$coverage_range)
    abricate_hits$start <-
      gsub(":[0-9]+", "", abricate_hits$coverage_range)
    
    #Select columns of interest
    abricate_hits <- abricate_hits %>%
      select(name,
             gene = GENE,
             ref_length,
             start,
             end,
             percentage = perc_coverage)
    
    #Convert start and end coordinates to numeric
    abricate_hits$start <- as.numeric(abricate_hits$start)
    abricate_hits$end <- as.numeric(abricate_hits$end)
    
    #Create an empty matrix equal to length of ref plasmid
    empty_plasrow <-
      rep(0, times = unique(abricate_hits$ref_length))
    
    #Create an empty matrix with n rows (n = sample size) with ncol == length(ref plasmid)
    empty_plasmatrix <- matrix(rep(empty_plasrow,
                                   times = length(unique(abricate_hits$name))),
                               nrow = length(unique(abricate_hits$name)))
    
    #Create a list of levels for sample names, a list of start coords and a list of end coords
    #and bind these in a list of lists
    start_ends <-
      list(as.list(as.integer(as.factor(abricate_hits$name))),
           as.list(as.integer(abricate_hits$start)),
           as.list(as.integer(abricate_hits$end)))
    
    #Create a counter
    counter <- 0
    
    #Create a progress bar
    pb1 = txtProgressBar(min = 0, max = nrow(abricate_hits), initial = 0) 
    message(paste("Mapping all abricate hits to coordinates of the reference plasmid..."))
    
    #Map the BLAST hits to our matrix of bp coordinates
    for (i in 1:nrow(abricate_hits)) {
      sample <- start_ends[[1]][[i]]
      start_coord <- start_ends[[2]][[i]]
      end_coord <- start_ends[[3]][[i]]
      empty_plasmatrix[sample, start_coord:end_coord] <- 1
      
      #Keep track of our progress with a counter and a completion bar
      #Blank message is for spacing and formatting of output
      counter <- counter + 1
      setTxtProgressBar(pb1,counter)
    }
    
    #Blank message is for spacing and formatting of output
    message("")
    
    #Rename matrix
    base_matrix <- empty_plasmatrix
    
    #Remove old matrix
    rm(empty_plasmatrix)
    
    # Convert matrix to a dataframe
    base_matrix <-
      as.data.frame(base_matrix, stringsAsFactors = FALSE)
    
    # Assign sample names to rows
    rownames(base_matrix) <- unique(abricate_hits$name)
    
    # Get the length of the reference sequence
    df_length <- length(abricate_hits)
    
    # Generate indices that cover blocks of 100 columns in base_matrix
    bin_ranges <- c(seq(from = 1, to = ref_length, by = 100))
    
    bin_ranges2 <-
      c(seq(from = 100, to = ref_length, by = 100), ref_length)
    
    
    # Split the ranges into two lists of  [1] start and [2] end indices
    bin_splits <- list(bin_ranges, bin_ranges2)
    
    # Initialise empty vector for loop below
    binned_hits <- vector()
    
    #Create a progress bar
    pb2 = txtProgressBar(min = 0, max = length(bin_ranges2), initial = 0) 
    message("Binning hits into chunks of 100bp...")
    
    #Create a counter
    counter <- 0
    
    # Binning loop
    for (i in 1:length(bin_ranges2)) {
      # If the last bin is only 1 base long then the as.matrix line won't work,
      #so we have to include the if statement below:
      if (i == length(bin_ranges) &
          bin_ranges[length(bin_ranges)] == bin_ranges2[length(bin_ranges2)]) {
        row_sum <- 1
      } else{
        # Generate row sums (i.e. Number of matching bases) for 100 column chunks of the base_matrix
        row_sum <-
          as.matrix(rowSums(base_matrix[, bin_splits[[1]][i]:bin_splits[[2]][i]]))
      }
      # Bind them together in a new vector
      binned_hits <- cbind(binned_hits, row_sum)
      
      #Keep track of our progress with a counter and a completion bar
      
      counter <- counter + 1
      setTxtProgressBar(pb2,counter)
    }
    
    #Blank message is for spacing and formatting of output
    message("")
    
    #Save rownames
    nems <- rownames(binned_hits)
    
    #Convert binned_hits to a data frame
    binned_hits <- as.data.frame(binned_hits)
    
    df6 <- as.data.frame(rowSums(binned_hits))
    
    df6$working_names <- rownames(df6)
    
    colnames(df6) <- c(plasmid_reference_name, "working_name")
    
    df6$Percent_hit <- round((df6[,1] / ref_length) * 100)
    
    if(writecsv != FALSE){
      # Write the DF
      write.csv(
        x = binned_hits,
        file = paste0(
          output_directory,
          "/",
          plasmid_reference_name,
          "_plasmid_coverage.csv"
        ),
        row.names = TRUE
      )
      
      
      message(
        paste0(
          "Writing to file ",
          output_directory,
          "/",
          plasmid_reference_name,
          "_plasmid_coverage.csv"
        )
      )
      
      write.csv(
        x = df6,
        file = paste0(
          output_directory,
          "/",
          plasmid_reference_name,
          "_plasmid_coverage_percentage.csv"
        ),
        row.names = TRUE
      )
      
      message(
        paste0(
          "Writing to file ",
          output_directory,
          "/",
          plasmid_reference_name,
          "_plasmid_coverage_percentage.csv"
        )
      )
      
    }else{
      message("Write to file left as FALSE - files unwritten")
    }
    
    #Assign the binned hits dataframe to the local environment
    assign(paste(plasmid_reference_name, "binned_hits", sep= "_"), binned_hits, envir=globalenv())
    
    if (path_to_tree != FALSE){
      if("Placeholder" %in% tree$tip.label){
        warning("It appears you have enabled the filter by tree option by providing a path to a tree file,
              however it appears that you should have provided the name of your reference strain for this tree and have not.
              Set the variable tree_reference_name to be equal to the name of the reference of your strain and run this tool again.")
      }
      
    }
  }
