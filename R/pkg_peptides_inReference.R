#' Peptide Reference Merger
#' @description Matches peptide reference information with peptidome sequences
#' @param peptidome list of all possible peptides for proteins of interest and their uniprot IDs
#' @param pep_reference dataframe with a "sequence" column of peptides and a column of n observations in a reference
#' @param reference_name the name of the reference, which is used to label the output columns
#' @param exp_counts_col the name of the reference column containing the n of peptide observations
#' @param detection_freq If TRUE, calculates the each peptide's frequency of appearance in the reference
#' @param detection_ratio If TRUE, calculates ratio of peptide observations to parent protein observations
#' @return Dataframe of peptidome with reference information
#' @import dplyr
#' @export


peptides_inReference <- function(peptidome, pep_reference, reference_name, exp_counts_col, detection_freq, detection_ratio){
  # peptidome: list of all possible peptides (minimum cols: "sequences")
  # pep_reference: list of reference peptides (minimum cols: "sequences", "n_obs_pep" (if non boolean))
  # reference_name: name of reference (used for resulting column names)
  # exp_counts_col: name of n_obs_pep column
  # detection_freq: If TRUE, adds column of num experiments peptide was observed in divided by total num of experiments in reference
  # detection_ratio: If TRUE adds column of num experiments peptide was observed in divded by num experiments parent protein was observed in (also adds n_obs_prot column)


  # ----- Defaults -----

  if(missing(exp_counts_col)){

    # default: no column for n_obs_pep if none specified
    exp_counts_col <- NA

  }

  if(missing(detection_ratio)){

    # default: no comparison if none specified
    detection_ratio <- FALSE

  }

  if(missing(detection_freq)){

    detection_freq <- FALSE

  }



  # ----- Check if Peptide in Reference -----

  # if n_obs_pep column not present, all peptidome peptides found in reference labeled TRUE
  if(is.na(exp_counts_col)==TRUE){

    # creates a col in peptidome with the reference name and labels all peptides FALSE
    peptidome[[reference_name]] <- FALSE
    # if peptide is found in reference, col = TRUE
    peptidome[[reference_name]][peptidome$sequence %in% pep_reference$sequence] <- TRUE

  }else{

    # name of n_obs_pep col with reference name added
    obs_pep_colname <- paste0(reference_name, "_n_obs_pep")

    # if n_obs_pep col already present in peptidome for current reference, skips addition
    if(obs_pep_colname %in% colnames(peptidome) == FALSE){

      # updates n_obs_pep col name in reference df before joining
      colnames(pep_reference)[BiocGenerics::grep(exp_counts_col, colnames(pep_reference))] <- obs_pep_colname
      # adds num of peptide observations to peptidome
      peptidome <- dplyr::left_join(peptidome, pep_reference, by = "sequence")
      # replaces NA values in num of observations with 0
      peptidome[[obs_pep_colname]][is.na(peptidome[[obs_pep_colname]])] <- 0

    }



    # ----- Calculate Frequency of Appearance -----

    if(isTRUE(detection_freq)){

      # largest n of experiments in reference that any peptide has been observed in
      max_n_obs <- max(peptidome[[obs_pep_colname]])
      # name of freq col for peptidome
      freq_colname <- paste0(reference_name, "_freq")
      # calculates proportion of experiments in reference that peptidome peptides were observed in
      peptidome[,freq_colname] <- peptidome[,obs_pep_colname]/max_n_obs

    }



    # ----- Calculate Peptide Detectability Ratio -----

    if(isTRUE(detection_ratio)){

      # determines largest num of experiments in reference that each uniprot ID/protein was observed in
      prot_nums <- dplyr::group_by(peptidome, uniprot) %>% dplyr::summarise(n_obs_prot = max(!!rlang::sym(obs_pep_colname)))
      # name of n_obs_prot col for peptidome
      obs_prot_colname <- paste0(reference_name, "_n_obs_prot")
      colnames(prot_nums)[2] <- obs_prot_colname

      # matches num of experiments in reference that proteins were observed in to uniprot IDs
      peptidome <- dplyr::left_join(peptidome, prot_nums, by = 'uniprot')
      # name of ratio col for peptidome
      ratio_colname <- paste0(reference_name, "_ratio")
      # calculates the ratio of experiments in reference that a peptide was observed in to the experiments that its parent protein was observed in
      peptidome[[ratio_colname]] <- peptidome[[obs_pep_colname]]/peptidome[[obs_prot_colname]]
      # converts NA's to 0 (result from n_obs_prot = 0)
      peptidome[[ratio_colname]][is.na(peptidome[[ratio_colname]])] <- 0


    }



  }

  return(peptidome)

}





