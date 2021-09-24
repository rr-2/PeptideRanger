#' Peptide Detectability Predictor
#' @description Predicts detectability scores for given peptide(s) using provided RF model.
#' @param peptides list of peptide amino acid sequences of interest.
#' @param prediction_model RF model for used for detectability predictions.
#' @param missed_cleavages list of peptide missed cleavages in same order as peptides. If missing, missed cleavages will be calculated from peptide sequence for model predictions.
#' @return dataframe of RF detectability scores with peptide sequences.
#' @export

peptide_predictions <- function(peptides, prediction_model, missed_cleavages){




  # ----- Missed Cleavages -----


  if(missing(missed_cleavages)){
  # if MCs not provided for peptides, calculates them from sequence

    missed_cleavages <- c()
    # creates a df for MC values

    for(i in 1:length(peptides)){
      # loops through all peptide sequences

      curr_pep <- peptides[i]
      # extracts current peptide sequence
      kr_locs <- stringr::str_locate_all(curr_pep, "K|R")
      # locates positions of all K and R residues in current peptide

      if(is.na(unlist(kr_locs)[1]) == TRUE){
      # if no K/R in sequence

        mc_count <- 0
        # set missed cleavage count to 0 for peptide

      }else{
      # if K/R in sequence

        mc_count <- as.numeric(sum(unlist(kr_locs)[1:dim(kr_locs[[1]])[1]] < nchar(curr_pep)))
        # counts num of non-C-term K/R residues, uses as num of missed cleavages

      }

      missed_cleavages <- c(missed_cleavages, mc_count)
      # adds MC of current peptide to missed_cleavages list

    }

  }



  if(length(missed_cleavages) != length(peptides)){
  # if the n of missed cleavage values is not equal to the n of peptides
    stop("Length of missed cleavages list must match the length of peptide list")

  }



  # ----- Phyiochemical Descriptor Calculation -----


  pc_descriptors <- data.frame(peptides, missed_cleavages)
  # turns peptide list and missed cleavages into a df for calculating physiochem properties
  colnames(pc_descriptors)[1] <- "sequence"

  pc_descriptors <- calculate_pc_descriptors(pc_descriptors = pc_descriptors)
  # calculates physiochemical descriptors needed for model predictions


  # ----- Detectability Predictions -----

  # adds RF scores to list of peptides for export
  pc_descriptors$RF_score <- stats::predict(prediction_model,pc_descriptors)



  # ----- Output Prep -----

  # filters for relevant columns for output
  output_list <- dplyr::select(pc_descriptors, sequence, RF_score)

  return(output_list)


}
