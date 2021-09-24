#' Peptide Prioritizer
#' @description Selects the top n peptides for a list of protein uniprot IDs given a set of priorities
#' @param uniprot_list a list of uniprot IDs for proteins of interest
#' @param max_n max number of peptides selected for each protein
#' @param peptidome list of all possible peptides for relevant proteome, their uniprot IDs, and the values for their non-RF priorities
#' @param prediction_model RF model used to predict scores for peptide detectability
#' @param priorities list of priority column names in peptidome
#' @param priority_thresholds list of minimum values that priorities must reach to be considered
#' @return list of n prioritized peptides for proteins of interest
#' @import randomForest
#' @import dplyr
#' @import Peptides
#' @export



prioritize_peptides <- function(uniprot_list, max_n, peptidome, prediction_model, priorities, priority_thresholds){
  # peptidome: list of all possible peptides (minimum cols: "sequences", "uniprot")
  # uniprot_list: list of uniprot IDs for proteins of interest
  # max_n: max number of peptides prioritized per protein in final output
  # prediction_model: RF model used to predict/rank peptide observability
  # priorities: list of column names in peptidome to consider ahead of RF scores during prioritization
  # priority_thresholds: list of value cutoffs or booleans for priorities


  # ----- Errors and Defaults -----

  if(missing(peptidome)){

    # default peptidome if none specified
    peptidome <- peptidome_SwissProt2018

    print('default peptidome: 2018 SwissProt Trypsin digest, 0-2 MCs, 6-25aa, unique to 1 protein')
  }



  if(missing(max_n)){

    # default max num of peptides if none specified
    max_n <- 5
    print(paste("default: max of", max_n, "peptides prioritized per protein"))
  }



  if(missing(prediction_model) & sum(grepl("RF_score", priorities)) != 0){

    # need model if "RF_score" is in priorities
    stop("Must set a prediction_model if RF_score used as a priority")

  }



  if(missing(prediction_model) & sum(grepl("RF_score", priorities)) == 0){

    # no model used if not specified and "RF_score" not in priorities
    prediction_model <- NA

  }



  if(sum(c("size", "missed_cleavages") %in% priorities) > 0){

    stop("size and missed_cleavages are not valid priorities")

  }



  if(length(priorities) != length(priority_thresholds)){

    # num of thresholds must match num of priorities
    stop("number of priorities does not match number of thresholds")

  }



  # ----- Subset Relevant Peptidome -----

  # removes duplicates from uniprot ID list
  uniprot_list <- unique(as.character(uniprot_list))

  # creates a subset of peptidome to only include peptides from proteins in the uniprot_list
  peps_of_interest <- peptidome[peptidome$uniprot %in% uniprot_list,]



  # ----- Random Forest Predictions -----

  if(isTRUE("RF_score" %in% priorities)){

    # creates a new df of the subset, only including peptide sequences and missed cleavages
    pc_descriptors <- dplyr::select(peps_of_interest, sequence, missed_cleavages)

    # calculates relevant physiochemical descriptors for model predictions
    pc_descriptors <- calculate_pc_descriptors(pc_descriptors = pc_descriptors)

    # RF uses peptide descriptors to predict detectability
    # adds predictions to column in list of peptides and corresponding parent protein IDs
    peps_of_interest$RF_score <- stats::predict(prediction_model,pc_descriptors)

  }



  # ----- Peptide Prioritization -----

  # creates empty df for final prioritized peptides
  peptide_predictions <- data.frame()

  # loop goes through each prot ID in uniprot_list
  for(i in 1:length(uniprot_list)){

    # filters for all peptides from current prot
    curr_peps <- dplyr::filter(peps_of_interest, uniprot == uniprot_list[i])
    # creates a df to add prioritized peptides to from current protein - meant for addition to peptide_predictions
    peps_toAdd <- data.frame()


    # loop through priorities
    for(k in 1:length(priorities)){

      # sets current priority
      curr_priority <- as.character(priorities[k])
      # sets corresponding threshold
      curr_threshold <- priority_thresholds[k]

      # if threshold = boolean
      if(isTRUE(curr_threshold) == TRUE | isFALSE(curr_threshold) == TRUE){

        # filters for peptides from curr_peps that have value = threshold in boolean curr_priority column
        priority_peps <- dplyr::filter(curr_peps, curr_peps[[curr_priority]] == curr_threshold)

        # adds max_n peptides from filtered list to df of peptides to add to final prioritized list
        peps_toAdd <- BiocGenerics::rbind(peps_toAdd, priority_peps[1:(max_n-nrow(peps_toAdd)),])
        # removes empty rows: appear if < max_n peptides in filtered list
        peps_toAdd <- na.omit(peps_toAdd)

      # if threshold != boolean (should be numeric)
      }else{

        # filters for peptides from curr_peps that have value >= threshold in numeric curr_priority column
        priority_peps <- dplyr::filter(curr_peps, curr_peps[[curr_priority]] >= curr_threshold)
        # arranges peptides in descending order by value in the curr_priority column
        priority_peps <- dplyr::arrange(priority_peps, desc(priority_peps[[curr_priority]]))

        # adds peptides in the top max_n rows (minus nrow peptides_toAdd) in the ordered list
        peps_toAdd <- BiocGenerics::rbind(peps_toAdd, priority_peps[1:(max_n-nrow(peps_toAdd)),])
        # removes empty rows: appear if < max_n peptides in filtered list
        peps_toAdd <- na.omit(peps_toAdd)


      }

      # after selecting peptides for each priority, removes selected peptides from curr_peps (the current protein's pool)
      curr_peps <- curr_peps[(curr_peps$sequence %in% peps_toAdd$sequence) == FALSE,]

      # if list of pepitdes for current protein is already = max_n, break priority loop and go to next protein
      if(nrow(peps_toAdd) == max_n){

        break

      }

    } # end of priority loop

    # adds list of up to max_n prioritized peptides to peptide_predictions
    peptide_predictions = BiocGenerics::rbind(peptide_predictions, peps_toAdd)

  } # end of protein loop


  # returns peptide_predictions
  return(peptide_predictions)

}











