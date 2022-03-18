#' Peptidome Generator
#' @description Generates a tryptic peptidome of unique peptides using a provided proteome. Requires bioconductor pacakge "cleaver".
#' @param proteome_dir directory to proteome fasta file
#' @param missed_cleavages range of missed cleavage events considered during in-silico digestion
#' @param synth_peps If TRUE, filters out peptides that are problematic for peptide synthesis. Currnetly includes peptides with: N-term Q/E residues and C-term P residues
#' @param aa_range range of peptide size by amino acid count included
#' @return Peptidome of unique peptides with parent protein uniprot IDs, sequence locations in parent proteins, and n missed cleavages
#' @examples
#' \dontrun{
#' create_peptidome(proteome_dir = "./proteome.fasta",
#'                  missed_cleavages = c(0,2),
#'                  synth_peps = FALSE,
#'                  aa_range = c(7-25))
#'}
#'
#' @export


create_peptidome <- function(proteome_dir, missed_cleavages, synth_peps, aa_range){
  # proteome_dir: directory to fasta file of all proteins in desired proteome
  # missed_cleavages: range of missed cleavages considered in in-silico digestion
  # synth_peps: if TRUE, then filters out peptides starting with Q and E
  # aa_range: range of number of amino acids in peptides that are filtered for



  # ----- Defaults -----

  size <- symbol <- uniprot <- start <- end <- NULL


  if(missing(missed_cleavages)){

    # default missed cleavages
    missed_cleavages <- c(0,2)

    print('default: 0-2 missed cleavages included in digestion')
  }


  if(missing(aa_range)){

    # default peptide size range
    aa_range <- c(6,25)

    print('default: 6-25 amino acid range')
  }


  if(missing(synth_peps)){

    # default synth_peps status
    synth_peps <- FALSE

    #print('default: synthetic peptide constraints not considered')
  }



  # ----- Digestions -----

  # directly accesses fasta file using directory and loads in proteome as large AAStringSet
  proteome <- Biostrings::readAAStringSet(file=proteome_dir)

  # creates empty peptidome df
  peptidome <- data.frame()

  ### CUSTOM DIGESTION CODE:
  # sets cut sites (Trypsin + loss of N-term M)
  # cut_sites <- c("^M|([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))", "((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))")
  # digestion with custom cut sites
  # curr_cleavage <- as.data.frame(cleaver::cleave(proteome, custom = cut_sites, missedCleavages = i, unique = FALSE))
  # curr_cleavage <- BiocGenerics::as.data.frame(cleaver::cleave(proteome, enzym = 'trypsin' , missedCleavages = i, unique = FALSE))

  cut_sites <- "^M|K|R"

  # loops through missed cleavage range
  for(i in missed_cleavages[1]:missed_cleavages[2]){

    # performs in silico trypsin digestion of proteome yielding peptides with i missed cleavages
    curr_cleavage <- BiocGenerics::as.data.frame(cleaver::cleave(proteome, custom = cut_sites, missedCleavages = i, unique = FALSE))

    # provides locations of peptides within parent proteins
    curr_positions <- BiocGenerics::as.data.frame(cleaver::cleavageRanges(proteome, custom = cut_sites, missedCleavages = i))

    # adds peptide start and end locations
    curr_cleavage$start <- curr_positions$start
    curr_cleavage$end <- curr_positions$end
    # adds the number of amino acids in peptides
    curr_cleavage$size <- curr_positions$width

    # adds the number of missed cleavages that resulted in the peptide
    curr_cleavage$missed_cleavages <- i
    # adds peptides from i missed cleavages and info to peptidome
    peptidome <- rbind(peptidome, curr_cleavage)

  }



  # ----- Column Renaming -----

  # extracts gene symbols from fasta headers
  peptidome$symbol <- stringr::str_match(peptidome$group_name,pattern = 'GN=(.+) PE')[,2]
  # extracts uniprot IDs from fasta headers
  peptidome$uniprot <- stringr::str_match(peptidome$group_name,pattern = 'sp\\|(.+)\\|')[,2]
  colnames(peptidome)[3] <- 'sequence'



  # ----- Filtering by Sequence Constraints -----

  # filters for peptides in the amino acid range specified
  peptidome = dplyr::filter(peptidome, size >= aa_range[1], size <= aa_range[2])

  # removes peptides with sequences containing "U"
  peptidome <- dplyr::filter(peptidome, grepl("U", peptidome$sequence) != TRUE)

  # if synthetic peptide parameter = TRUE
  if(synth_peps == TRUE){

    # remove peptides starting with Q or E
    peptidome = peptidome[!grepl('^Q|^E', peptidome$sequence),]
    # removes peptides ending with P
    peptidome = peptidome[!grepl("P$", peptidome$sequence),]

  }



  # ----- Removal of Non-Unique Peptides -----

  # removes duplicates, leaving only peptides unique to one protein
  peptidome = peptidome[!(BiocGenerics::duplicated(peptidome$sequence) | BiocGenerics::duplicated(peptidome$sequence, fromLast = TRUE)),]



  # ----- Final Output -----

  # selects a subset of columns
  peptidome <- dplyr::select(peptidome, symbol, uniprot, missed_cleavages, sequence, start, end, size)


  return(peptidome)

}
