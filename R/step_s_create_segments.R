#' Create fixed-size genomic segments (SAF format)
#'
#' Splits autosomes (chr1-chr22) into fixed-size segments and writes a SAF file
#' compatible with `Rsubread::featureCounts`.
#'
#' This function uses built-in chromosome sizes for hg19 or hg38.
#' It strictly filters the genome to include **only** chromosomes named "chr1" through "chr22".
#' Sex chromosomes (chrX, chrY) and mitochondrial DNA (chrM) are excluded.
#'
#' @param file_saf_out Output path for the SAF file.
#' @param genome Genome build to use: "hg19" (default) or "hg38".
#' @param segment_size Character like "5kb" or numeric bp (e.g. 5000).
#'   Default is "5kb".
#'
#' @return Invisibly returns the segment data frame. Writes \code{file_saf_out} to disk.
#'
#' @examples
#' \dontrun{
#' # Create 5kb segments for hg19 (default)
#' cacti_s_create_segments("inst/extdata/test_results/test_segments.saf")
#'
#' # Create 1kb segments for hg38
#' cacti_s_create_segments("inst/extdata/test_results/test_segments.saf", genome = "hg38", segment_size = "1kb")
#' }
#' @export
cacti_s_create_segments <- function(
    file_saf_out,
    genome = c("hg19", "hg38"),
    segment_size = "5kb"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr required.")

  # --- 1. Parse Segment Size ---
  if (is.character(segment_size) && stringr::str_detect(segment_size, "^\\d+[.]?\\d*[kK][bB]$")) {
    # Convert "5kb" -> 5000
    seg_bp <- as.numeric(stringr::str_extract(segment_size, "^\\d+[.]?\\d*")) * 1000
  } else if (is.numeric(segment_size)) {
    seg_bp <- segment_size
  } else {
    stop("segment_size must be character (e.g. '5kb') or numeric bp.")
  }

  # --- 2. Select Built-in Chromosome Sizes ---
  genome <- match.arg(genome)
  message("Using built-in chromosome sizes for: ", genome)

  if (genome == "hg38") {
    # Standard UCSC hg38 sizes
    sizes <- c(
      chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555, chr5=181538258,
      chr6=170805979, chr7=159345973, chr8=145138636, chr9=138394717, chr10=133797422,
      chr11=135086622, chr12=133275309, chr13=114364328, chr14=107043718, chr15=101991189,
      chr16=90338345, chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
      chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415, chrM=16569
    )
  } else {
    # Standard UCSC hg19 sizes
    sizes <- c(
      chr1=249250621, chr2=243199373, chr3=198022430, chr4=191154276, chr5=180915260,
      chr6=171115067, chr7=159138663, chr8=146364022, chr9=141213431, chr10=135534747,
      chr11=135006516, chr12=133851895, chr13=115169878, chr14=107349540, chr15=102531392,
      chr16=90354753, chr17=81195210, chr18=78077248, chr19=59128983, chr20=63025520,
      chr21=48129895, chr22=51304566, chrX=155270560, chrY=59373566, chrM=16571
    )
  }

  chr_df <- data.frame(chr = names(sizes), size = as.numeric(sizes), stringsAsFactors = FALSE)

  # --- 3. Filter for Autosomes (chr1-chr22) Only ---
  target_chrs <- paste0("chr", 1:22)
  chr_df <- chr_df |> dplyr::filter(chr %in% target_chrs)

  # Sort numerically (chr1, chr2, ..., chr22)
  chr_df <- chr_df[match(target_chrs, chr_df$chr), ]
  chr_df <- na.omit(chr_df)

  message("Generating ", segment_size, " segments for ", nrow(chr_df), " autosomes (chr1-chr22)...")

  # --- 4. Generate Segments ---
  res_list <- lapply(seq_len(nrow(chr_df)), function(i) {
    cc <- chr_df$chr[i]
    sz <- chr_df$size[i]

    starts <- seq(1, sz, by = seg_bp)
    ends   <- pmin(starts + seg_bp - 1, sz)

    # Construct SAF columns
    # GeneID: "chr_start_end"
    ids <- paste0(cc, "_", starts, "_", ends)

    data.frame(
      GeneID = ids,
      Chr    = cc,
      Start  = starts,
      End    = ends,
      Strand = ".",
      stringsAsFactors = FALSE
    )
  })

  saf_df <- dplyr::bind_rows(res_list)

  # --- 5. Write Output ---
  saf_df$Start <- as.integer(saf_df$Start)
  saf_df$End   <- as.integer(saf_df$End)

  data.table::fwrite(saf_df, file_saf_out, sep = "\t", scipen = 999)

  message("Success! Created ", nrow(saf_df), " segments. Saved to: ", file_saf_out)
  invisible(saf_df)
}
