#' @noRd

calculate_pc_descriptors <- function(pc_descriptors){
# calculates relevant physiochemical descriptors for peptide prediction models
# pc_descriptors: df with a column of peptide sequences with column name "sequences"

  # --- peptide physiochem descriptors ---

  # naive descriptors
  pc_descriptors$mw=Peptides::mw(pc_descriptors$sequence)
  pc_descriptors$charge_ph7=Peptides::charge(pc_descriptors$sequence)
  pc_descriptors$hmoment=Peptides::hmoment(pc_descriptors$sequence)

  # saves col num before group calculations
  pre_cols <- ncol(pc_descriptors)

  # statistically derived descriptor groups
  pc_descriptors$mswhimScores=Peptides::mswhimScores(pc_descriptors$sequence)
  pc_descriptors$blosumIndices=Peptides::blosumIndices(pc_descriptors$sequence)
  pc_descriptors$fasVectors=Peptides::fasgaiVectors(pc_descriptors$sequence)
  pc_descriptors$VHSE_Scales=Peptides::vhseScales(pc_descriptors$sequence)
  pc_descriptors$z_scales=Peptides::zScales(pc_descriptors$sequence)
  pc_descriptors$st_scales=Peptides::stScales(pc_descriptors$sequence)
  pc_descriptors$prot_FP=Peptides::protFP(pc_descriptors$sequence)

  # saves col num after group calculations
  post_cols <- ncol(pc_descriptors)

  pc_descriptors$blsm5=unlist(lapply(pc_descriptors$blosumIndices,function(X) X[[5]]))
  pc_descriptors$blsm10=unlist(lapply(pc_descriptors$blosumIndices,function(X) X[[10]]))
  pc_descriptors$fasV3=unlist(lapply(pc_descriptors$fasVectors,function(X) X[[3]]))
  pc_descriptors$fasV5=unlist(lapply(pc_descriptors$fasVectors,function(X) X[[5]]))
  pc_descriptors$FP5=unlist(lapply(pc_descriptors$prot_FP,function(X) X[[5]]))
  pc_descriptors$FP8=unlist(lapply(pc_descriptors$prot_FP,function(X) X[[8]]))
  pc_descriptors$MSWHIM2=unlist(lapply(pc_descriptors$mswhimScores,function(X) X[[2]]))
  pc_descriptors$MSWHIM3=unlist(lapply(pc_descriptors$mswhimScores,function(X) X[[3]]))
  pc_descriptors$ST7=unlist(lapply(pc_descriptors$st_scales,function(X) X[[7]]))
  pc_descriptors$VHSE8=unlist(lapply(pc_descriptors$VHSE_Scales,function(X) X[[8]]))
  pc_descriptors$zSc1=unlist(lapply(pc_descriptors$z_scales,function(X) X[[1]]))


  # removes sequence column and unused group calculation columns
  pc_descriptors <- pc_descriptors[-c((pre_cols+1):post_cols)]

  return(pc_descriptors)


}
