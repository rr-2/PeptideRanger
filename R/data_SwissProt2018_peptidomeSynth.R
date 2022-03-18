#' Uniprot/Swiss-Prot Synthetic-Friendly Peptidome
#'
#' Filtered version of SwissProt2018_peptidome that does not include peptides that are challening to synthesize. This currently excludes peptides with N-term Q and E residues and with C-term P residues.
#' @format A data frame with 1694311 rows and 7 variables:
#' \describe{
#'   \item{symbol}{gene symbol of parent protein}
#'   \item{uniprot}{uniprot ID of parent protein}
#'   \item{missed_cleavages}{number of MCs that peptide resulted from}
#'   \item{sequence}{amino acid sequence of peptide}
#'   \item{start}{parent protein position of first residue in peptide}
#'   \item{end}{parent protein position of last residue in peptide}
#'   \item{size}{number of amino acid residues in peptide}

#'
#' }
#'
"SwissProt2018_peptidome_synth"
