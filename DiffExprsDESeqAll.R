

#####
### Function to calculate differential expression contrasts using DESeq2
#####

DiffExprsDESeq <- function(DESeq_object, 
                           formula = design(DESeq_object), 
                           contrast_groups,
                           stratification = NULL){
  
  ##arguments:
  #1 the DESeq objct generated from DESeqDataSetFromMatrix()
  #2 an R formula or a string of the R formula for GLM fitting
  #3 a list or vector containing the groups between which the contrasts need to
  #     be made
  #4 an optional vector which can be used to fine tune analyses by stratifying
  #     the provided factors further by their respective values ; e.g., here it 
  #     has been used to create tissue---.cellTypeT and tissue---.cellTypeNT 
  #     groups for the formula ~0+tissue:cellType which were otherwise being 
  #     grouped into a single object and would have become confounded between T 
  #     and NT
  
  #checking and loading libraries
  require(DESeq2)
  require(PTXQC)
  require(gtools)
  
  #formatting 'formula' as formula in case it is provided as string
  formula = formula(formula)
  
  #assigning formula to the design matrix
  design(DESeq_object) <- formula
  
  #performing stock DESeq2 pipeline on DESeq object
  DESeq_object <- DESeq(object = DESeq_object, 
                        test = "Wald", 
                        betaPrior = F, 
                        minReplicatesForReplace = Inf, 
                        parallel = T)
  
  #available result names to draw contrasts
  available_contrasts <- resultsNames(DESeq_object)
  
  #dynamically generating a list of sets of possible contrast names that can be 
  #used to draw multiple combinations of interesting contrasts between the given
  #contrast groups
  available_contrasts_list <- list()
  
  #patterns available to search for from the factor names
  available_patterns <- colnames(colData(DESeq_object))
  
  #removing 'sizeFactor' and 'replaceable' column names as these are not 
  #annotations we use for generating contrast patterns
  available_patterns <- available_patterns[
    !(available_patterns %in% c("sizeFactor", "replaceable"))
    ]
  
  #making sub-patterns for further stratified factors of interest if needed
  if(!is.null(stratification) & !all(is.na(stratification))){
    for (i in seq(1, length(stratification))){

      new_params <-
        as.vector(unique(colData(DESeq_object)[
          which(colnames(colData(DESeq_object)) %in% stratification[i])
          ])[,1])

      available_patterns <- c(available_patterns,
                              as.vector(outer(X = stratification[i],
                                              Y = new_params,
                                              FUN = paste0))
                              )
    }
  }
  
  writeLines(paste0("Patterns available for generating contrasts are: ", 
             paste(available_patterns, collapse = ", "))
             )
  
  #permutations of negation of a factor name to allow a pool of all possible 
  #contrasts
  permutations <- gtools::permutations(n = 2, 
                                       r = length(available_patterns), 
                                       v = c(TRUE, FALSE), 
                                       repeats.allowed = T)
  
  #using grep recursively to generate lists of patterns which permute though the
  #inclusion and exclusion of different patterns using inversion vector
  #permutation[1,] is always all FALSE
  for (i in seq(2, nrow(permutations))) {
    
    intersection <- RecurGrep(available_contrasts = available_contrasts,
                              available_patterns = available_patterns,
                              inversion_vector = which(!permutations[i,]),
                              iter = 1)
    
    #available contrasts are added only if an intersection result is found
    if(length(intersection)>0){
      available_contrasts_list[[length(available_contrasts_list)+1]] <- 
      intersection
      
      #constructing the name/label for the comparison
      name_table <- do.call("rbind", strsplit(x = intersection, split = "\\."))
      label <- PTXQC::LCSn(strings = name_table[,1])
      
      if(ncol(name_table)>1){
        for(k in seq(2, ncol(name_table))){
          label <- paste0(label, ".", PTXQC::LCSn(strings = name_table[,k]))
        }
      }
      
      #assigning constructed name
      names(available_contrasts_list)[length(available_contrasts_list)] <- 
        paste0("Comparison on ", label)
    
    }
  }
  
  #initialise a list that will store the contrast results
  results_list <- list()
  
  #looping through the list of available contrasts and generating the results
  #list
  for (i in seq(1, length(available_contrasts_list))){
    
    #move to the next iteration if only a single group is found for the contrast
    if(length(available_contrasts_list[[i]])<2) next
    
    writeLines(paste0("Drawing contrasts for the permutation: ", 
                      tail(x = unlist(strsplit(
                        x = names(available_contrasts_list[i]), 
                        split = " ")), 
                        n = 1)
                      )
               )
    
    #otherwise call the function to generate the contrasts on the groups
    results_list[[length(results_list)+1]] <- list()
    results_list[[length(results_list)]] <- 
      GenerateContrasts(contrast_groups = contrast_groups, 
                        available_contrasts = available_contrasts_list[[i]],
                        results_list = results_list[[length(results_list)]],
                        DESeq_object = DESeq_object,
                        #prefix = names(available_contrasts_list[i])
                        )
    
    #if no result was calculated, remove the last object and go to the next 
    #iteration
    if(length(results_list[[length(results_list)]])<1){
      results_list[[length(results_list)]] <- NULL
      next
    }
    
    #otherwise, write the label of the new result
    names(results_list)[length(results_list)] <- 
      names(available_contrasts_list[i])
    
  }

  #returning the contrast results list
  return(results_list)
  
}


###
## Companion function to reduce code repetition further
###
GenerateContrasts <- function(contrast_groups,
                              available_contrasts,
                              results_list,
                              DESeq_object,
                              prefix = ""){
  
  #check if contrast_groups is a list, otherwise convert it to a list
  if(typeof(contrast_groups) != "list"){
    contrast_groups <- as.list(contrast_groups)
  }
  
  #generate contrasts iteratively between provided contrast groups
  #iterate i from 1 to n-1
  for(i in seq(1, length(contrast_groups)-1)){
    
    print(paste0("Analysing for: ", 
                 paste(contrast_groups[[i]], collapse = "_")))
    
    #extracting the contrast name(s) for group i
    contrast_i <- available_contrasts[
      grep(pattern = paste(contrast_groups[[i]], collapse = "|"), 
           x = available_contrasts, 
           ignore.case = T)
      ]
    
    #moving onto next iteration if no contrasts matched the pattern
    if(length(contrast_i)<1){
      writeLines(paste0("Could not find any available contrasts matching the ", 
                        "pattern: ", paste(
                          contrast_groups[[i]], collapse = " or ")
                        )
                 )
      next
    }
    
    #extracting the title of the contrast group i
    name_i <- paste(contrast_groups[[i]], collapse = "_")
    
    #iterate j from i+1 to n
    for(j in seq(i+1, length(contrast_groups))){
      #extracting the contrast_name(s) for group j
      contrast_j <- available_contrasts[
        grep(pattern = paste(contrast_groups[[j]], collapse = "|"), 
             x = available_contrasts, 
             ignore.case = T)
        ]
      
      #moving onto next iteration if no contrasts matched the pattern
      if(length(contrast_j)<1){
        writeLines(paste0("Could not find any available contrasts matching ", 
                          "the pattern: '", paste(
                            contrast_groups[[j]], collapse = " or "), 
                          "'"))
        next
      }
      
      #extracting the title of the contrast group j
      name_j <- paste(contrast_groups[[j]], collapse = "_")
      
      #check and skip to next iteration if both contrast_i and contrast_j are 
      #the same
      if(length(intersect(x = contrast_i, y = contrast_j))>0){
        writeLines(paste0("Overlapping term(s): ", 
                          paste(intersect(x = contrast_i, y = contrast_j), 
                                collapse = ", "), 
                          " found. Terms cannot overlap for drawing ", 
                          "contrasts. This contrast will be skipped."))
        next
      } 
      
      #performing the contrast
      results_list[[length(results_list)+1]] <-
        results(object = DESeq_object, test = "Wald",
                contrast = list(contrast_i, contrast_j),
                listValues = c(1/length(contrast_i), -1/length(contrast_j)),
                cooksCutoff = F,
                independentFiltering = F,
                parallel = T)
      
      #the below code should be commented; it performs the contrast and applies
      # log fold change shrinkage
      # results_list[[length(results_list)]] <-
      #   lfcShrink(dds = DESeq_object,
      #             contrast = list(contrast_i, contrast_j),
      #             res = results_list[[length(results_list)]],
      #             type = "ashr",
      #             lfcThreshold = 0,
      #             parallel = T)
      
      #creating name of the contrast result
      names(results_list)[length(results_list)] <- 
        paste0(prefix, name_i, "__Vs__", name_j)
    }
    
  }
  
  return(results_list)
  
}


###
## Recursive grep function to extract all possible 
###
RecurGrep <- function(available_contrasts = NULL, 
                      available_patterns = c(""), 
                      inversion_vector = c(0), 
                      iter = 1){
  
  #check if deseq object is missing
  if(is.null(available_contrasts)) stop("No available contrast names provided")
  
  #check if list of available patterns to search are less than two
  #the recursive function is equipped to use the last two available patterns 
  #together for the intersection
  if(length(available_patterns)<2){
    if(available_patterns == "" | is.null(available_patterns)){
      stop("Vector of available patterns is empty")
    }
    else stop("You need at least two names in the vector of available 
                patterns")
  }
  
  #check if the x component of the intersection needs to have its results 
  #inverted
  if(iter %in% inversion_vector) invert_i <- TRUE
  else invert_i <- FALSE
  
  #check if there are still more than two patterns available so that the
  #function can continue calling itself recursively
  if(length(available_patterns)>2){
    grep_result <- RecurGrep(available_contrasts = available_contrasts,
                             available_patterns = available_patterns[-1], 
                             inversion_vector = inversion_vector,
                             iter = iter+1)
    
    #once the recursion comes back, the intersection is performed
    return(available_contrasts[intersect(
      x = grep(pattern = available_patterns[1], 
               x = available_contrasts, 
               ignore.case = T, 
               invert = invert_i), 
      y = which(available_contrasts %in% grep_result))
      ])
    
  }
  
  #the following condition will only be triggered when the number of patterns 
  #in available patterns is exactly two - that is the last two patterns
  else{
    
    #check if the results of the last pattern would need to be inverted
    if ((iter+1) %in% inversion_vector) invert_j <- TRUE
    else invert_j <- FALSE
    
    #perform intersection
    return(available_contrasts[
      intersect(x = grep(pattern = available_patterns[1], 
                         x = available_contrasts, 
                         ignore.case = T, 
                         invert = invert_i), 
                y = grep(pattern = available_patterns[2], 
                         x = available_contrasts, 
                         ignore.case = T, 
                         invert = invert_j)
      )
      ]
    )
    
  }
  
}


