#' List of Transmembrane Peptides
#'
#' List of peptides that overlap with transmembrane domain sequence annotations downloaded from SwissProt-UniProt on 2021/01/14.
#'
#' @format A data frame with 14121 rows and 6 variables:
#' \describe{
#'   \item{symbol}{gene symbol of parent protein}
#'   \item{uniprot}{uniprot ID of parent protein}
#'   \item{sequence}{amino acid sequence of peptide}
#'   \item{size}{number of amino acid residues in peptide}
#'   \item{start}{parent protein position of first residue in peptide}
#'   \item{end}{parent protein position of last residue in peptide}
#'
#' }
#'
"tm_peptides"
