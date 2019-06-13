#' Assign NTE wrapper
#'
#' This function is multicore wrapper for gene_NTE_assigner().
#' @param genome BSgenome name, compatible with gtf
#' @param gtf_table gtf data-frame annotation
#' @param start_codons vector of triplets to be searched as start codons.
#' @param stop_codons vector of triplets to be searched as stop codons
#' @param max_nte_length maximum NTE length in aminoacids.
#' @param genetic_code GENETIC_CODE in biostrings genetic code format.
#' @param context_span two element vector, corresponding to start and end of the sequence relative to first position of start codon 
#' @param truncations logical, if TRUE truncations will be searched for.
#' @param truncation_start_codons triplets to be searched for as truncations start codons, "ATG" by default.
#' @param max_truncation_length maximum truncation length in aminoacids.
#' @param cores number of cores to use. Passes argument to mcmapply() mc.cores argument. 
#' @return MitoCrypt_table of all isoforms called all genes in gtf table.
#' @keywords NTE extension truncation
#' @seealso gene_NTE_assigner
#' @family MitoCrypt NTE_assign
#' @export
#' @examples
#'


  assign_NTE_wrapper <- function(genome, gtf_table, start_codons = c("ACG", "ATA", "ATC", "ATG", "ATT", "CTG","GTG", "TTG"), stop_codons = (c("TAA", "TAG", "TGA")), max_nte_length = 300, genetic_code = GENETIC_CODE, context_span = c(-6,4), truncations = TRUE, truncation_start_codons = "ATG",max_truncation_length = 50, cores = 1){
    gtf_table <- gtf_table[which(gtf_table$seqid %in% seqnames(genome)),]
    gtf_table <- gtf_table[which(gtf_table$feature == "CDS"),]
    gtf_table <- gtf_table[which(gtf_table$gene_biotype == "protein_coding"),]
    
    NTE_table <- data.frame(Gene.Name = character(),
                            Standard.Name = character(),
                            Chromosome = character(),
                            Start.Coordinate = integer(),
                            End.Coordinate = integer(),
                            Strand = character(),
                            Start.Codon = character(),
                            Start.Codon.Context = character(),
                            NTE.Sequence = character(),
                            NTE.Length = integer(),
                            Isoform.Sequence = character(),
                            blockCount = integer(),
                            blockSizes = character(),
                            blockStarts = character(),
                            stringsAsFactors = FALSE)
    
    list_of_dataframes <- mclapply(levels(as.factor(gtf_table$transcript_id)), function(x) gene_NTE_assigner(genome = genome, gtf_table = gtf_table, start_codons = start_codons, stop_codons = stop_codons, max_nte_length = max_nte_length, genetic_code = genetic_code, context_span = context_span, truncations = truncations, truncation_start_codons = truncation_start_codons,max_truncation_length = max_truncation_length, x),mc.cores = cores)
    NTE_table <- do.call("rbind", list_of_dataframes) 
    return(NTE_table)
  }
