#' RF Model Training
#' @description Trains a random forest model using a peptidome with reference information
#' @param peptidome List of all possible peptide sequences for relevant proteome (total_trainingSet = FALSE) OR a list of peptide sequences in a custom training set (total_trainingSet = TRUE). Must also include peptide n of observations (n_obs_pep) and parent protein n of observations (n_obs_prot) in a desired reference, and peptide missed cleavages.
#' @param reference_name Name of reference (same as column names)
#' @param total_trainingSet If TRUE, model will be trained on all peptides in a custom input peptidome. If FALSE, peptides will be selected from whole input peptidome to create a training set of 10k peptides (same method as models provided in package).
#' @return RF model for predicting peptide detectability
#' @import randomForest
#' @export


train_RFmodel <- function(peptidome, reference_name, total_trainingSet){


  # ----- Errors and Defaults -----

  if(missing(total_trainingSet)){
  # if total_trainingSet not stated, set to FALSE
    total_trainingSet <- FALSE

  }



  # ----- Training Set Filtering -----


  ratio_name <- paste0(reference_name, "_ratio")
  # save name of ratio column for desired reference

  if(total_trainingSet == FALSE){
  # if total peptidome input/not a custom training set

    n_obs_prot_name <- paste0(reference_name, "_n_obs_prot")
    # save name of n_obs_prot column for desired reference

    top10 <- quantile(peptidome[[n_obs_prot_name]], prob = 0.9)
    # calculate value for 0.9 quantile of n_obs_prot

    training_peptides <- peptidome[peptidome[[n_obs_prot_name]] >= top10,]
    # filter peptides reaching n_obs_prot value = 0.9 quantile (peptides from frequently observed proteins)

    detected_peptides <- training_peptides[training_peptides[[ratio_name]] >= quantile(training_peptides[[ratio_name]], prob = 0.9),]
    # filter peptides reaching 0.9 quantile of reference ratio values (peptides with good detectability)
    detected_peptides <- sample_n(detected_peptides, 5000)
    # randomly sample 5k peptides from pool of peptides with good detectability

    undetected_peptides <- training_peptides[training_peptides[[ratio_name]] == 0,]
    # filter for unobserved peptides (ratio = 0)
    undetected_peptides <- dplyr::sample_n(training_peptides, 5000)
    # randomly sample 5k peptides from pool of unobserved peptides

    training_peptides <- rbind(detected_peptides, undetected_peptides)
    # create training 10k training set of peptides with good detectability and no observations for frequently detected proteins

  }else{

    training_peptides <- peptidome
    # If total_trainingSet = TRUE, training peptides = entire input peptidome

  }


  # ----- Descriptor Preparation -----

  pc_descriptors <- training_peptides[,c(ratio_name, "sequence", "missed_cleavages")]
  # create df of peptide detectability ratios, sequences, and missed cleavages
  colnames(pc_descriptors)[1] <- "ratio"
  # generalize ratio column name

  pc_descriptors <- calculate_pc_descriptors(pc_descriptors = pc_descriptors)
  # calculate pysio chemical descriptors for training peptides
  pc_descriptors$sequence <- NULL
  # remove sequences from df


  # ----- Model Training -----

  start.time<-proc.time()
  # record start time
  RFmodel <- randomForest::randomForest(ratio~., data = pc_descriptors, ntree = 500, importance = TRUE )
  # train RF model to predict ratios
  stop.time<-proc.time()
  # record stop time
  run.time<-stop.time-start.time
  # calculate run time

  print(paste0("Model training time = ", round(run.time[3]/60,2), " minutes"))

  return(RFmodel)


}


