#!/usr/bin/env Rscript
#
# BSF R script to refine the GATK Callable Loci analyis.
#
# The coverage assessment loads (Ensembl) exon information, collates (or projects)
# all (overlapping) exons into transcribed regions on the genome, overlaps those with
# the target regions to get transcribed target regions and finally overlaps those with
# non-callable loci to get the minimal set of problematic regions.
# Each problematic region is annotated with target name, as well as exon, transcript and
# gene information. A summary data frame of metrics collected along the procedure is also
# written to disk.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility of
# the Research Center for Molecular Medicine (CeMM) of the Austrian Academy of
# Sciences and the Medical University of Vienna (MUW).
#
#
# This file is part of BSF R.
#
# BSF R is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BSF R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BSF R.  If not, see <http://www.gnu.org/licenses/>.

suppressPackageStartupMessages(expr = library(package = "optparse"))

argument_list <- parse_args(object = OptionParser(
  option_list = list(
    make_option(
      opt_str = c("--verbose", "-v"),
      action = "store_true",
      default = TRUE,
      help = "Print extra output [default]",
      type = "logical"
    ),
    make_option(
      opt_str = c("--quiet", "-q"),
      action = "store_false",
      default = FALSE,
      dest = "verbose",
      help = "Print little output",
      type = "logical"
    ),
    make_option(
      opt_str = c("--exons"),
      dest = "exon_path",
      help = "File path to the gene, transcript and exon annotation GTF",
      type = "character"
    ),
    make_option(
      opt_str = c("--exon-flanks"),
      default = 0L,
      dest = "exon_flanks",
      help = "Exon flanking regions [0]",
      type = "integer"
    ),
    make_option(
      opt_str = c("--no-filter"),
      action = "store_true",
      default = FALSE,
      dest = "no_filter",
      help = "Do not filter GTF annotation by 'tag \"basic\";' [FALSE]",
      type = "logical"
    ),
    make_option(
      opt_str = c("--callable-loci"),
      dest = "callable_loci_path",
      help = "File path to the GATK Callable Loci BED",
      type = "character"
    ),
    make_option(
      opt_str = c("--targets"),
      dest = "target_path",
      help = "File path to the enrichment targets BED",
      type = "character"
    ),
    make_option(
      opt_str = c("--plot-width"),
      default = 7.0,
      dest = "plot_width",
      help = "Plot width in inches [7.0]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--plot-height"),
      default = 7.0,
      dest = "plot_height",
      help = "Plot height in inches [7.0]",
      type = "numeric"
    )
  )
))

if (file.size(argument_list$callable_loci_path) > 1e+09) {
  message("The file specified by the --callable-loci option is too large to process")
  print(x = sessionInfo())
  quit()
}

suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

# Keep overall statistics in a summary data frame.
summary_frame <- data.frame(stringsAsFactors = FALSE)
i <- 1L

# Import Ensembl annotation -----------------------------------------------


# Import Ensembl gene, transcript and exon annotation as GRanges vector object,
# where only exon components are relevant for this analysis.
summary_frame[i, "exon_path"] <- argument_list$exon_path
message(paste0("Importing exon ranges: ", summary_frame[i, "exon_path"]))
exon_ranges <-
  import(con = summary_frame[i, "exon_path"],
         genome = "hs37d5",
         feature.type = "exon")
# Ensembl now annotates a "tag" in GFF files with value "basic" indicating standard (basic) transcript models.
# Can the exon ranges be subset by such a tag?
if ("tag" %in% names(x = mcols(x = exon_ranges)) &
    !argument_list$no_filter) {
  message("Filtering by GTF 'tag = \"basic\"' annotation")
  summary_frame[i, "exon_number_raw"] <- length(x = exon_ranges)
  message(paste0("Number of raw exon ranges: ", summary_frame[i, "exon_number_raw"]))
  summary_frame[i, "exon_width_raw"] <- sum(width(x = exon_ranges))
  message(paste0("Cumulative width of raw exon ranges: ", summary_frame[i, "exon_width_raw"]))
  # Use the %in% operator for character matching as it sets NA values to FALSE, automatically.
  basic_exon_ranges <-
    exon_ranges[mcols(x = exon_ranges)$tag %in% "basic", ]
  if (length(x = basic_exon_ranges) > 0L) {
    # Only replace the original Exon ranges if there were any "basic" matches.
    exon_ranges <- basic_exon_ranges
  }
  rm(basic_exon_ranges)
}
summary_frame[i, "exon_number"] <- length(x = exon_ranges)
message(paste0("Number of exon ranges: ", summary_frame[i, "exon_number"]))
summary_frame[i, "exon_width"] <- sum(width(x = exon_ranges))
message(paste0("Cumulative width of exon ranges: ", summary_frame[i, "exon_width"]))

# Apply flanking regions --------------------------------------------------


# Apply flanking regions, by default 0L, to the exon ranges.
exon_ranges <-
  resize(
    x = exon_ranges,
    width = width(x = exon_ranges) + argument_list$exon_flanks,
    fix = "end"
  )
exon_ranges <-
  resize(
    x = exon_ranges,
    width = width(x = exon_ranges) + argument_list$exon_flanks,
    fix = "start"
  )
summary_frame[i, "exon_flank_width"] <- sum(width(x = exon_ranges))
message(paste0("Cumulative width of exon ranges with flanks: ", summary_frame[i, "exon_flank_width"]))

# Reduce non-redundant Ensembl exons --------------------------------------


# Reduce the non-redundant Ensembl exons to their footprint on the genome to get
# transcribed regions. The revmap column contains the mapping to original
# exon_ranges components.
message("Reducing exon ranges to transcribed ranges.")
transcribed_ranges <-
  reduce(
    x = exon_ranges,
    drop.empty.ranges = TRUE,
    with.revmap = TRUE,
    ignore.strand = TRUE
  )
summary_frame[i, "transcribed_number"] <-
  length(x = transcribed_ranges)
message(paste0("Number of transcribed ranges: ", summary_frame[i, "transcribed_number"]))
summary_frame[i, "transcribed_width"] <-
  sum(width(x = transcribed_ranges))
message(paste0("Cumulative width of transcribed ranges: ", summary_frame[i, "transcribed_width"]))

# Read target regions -----------------------------------------------------


# Read the file of targeted regions, if available. Although the GATK Callable
# Loci analysis is generally only run on these target regions, this file
# provides the target (probe) names of the enrichment design.
target_ranges <- NULL
constrained_ranges <- NULL
if (!is.null(x = argument_list$target_path)) {
  summary_frame[i, "target_path"] <- argument_list$target_path
  message(paste0("Importing target range annotation: ", summary_frame[i, "target_path"]))
  # The rtrackayer::import() function reads the genome version from the "db" attribute of the BED "track" line.
  target_ranges <- import(con = summary_frame[i, "target_path"])
  summary_frame[i, "target_number_raw"] <-
    length(x = target_ranges)
  message(paste0("Number of target ranges: ", summary_frame[i, "target_number_raw"]))
  summary_frame[i, "target_width_raw"] <-
    sum(width(x = target_ranges))
  message(paste0("Cumulative width of target ranges: ", summary_frame[i, "target_width_raw"]))
  
  message("Overlapping target and transcribed ranges.")
  overlap_frame <-
    mergeByOverlaps(query = target_ranges, subject = transcribed_ranges)
  overlap_ranges <- overlap_frame[, "target_ranges"]
  # Adjust start and end to the minimally overlapping regions.
  constrained_ranges <-
    GRanges(
      seqnames = seqnames(x = overlap_frame$target_ranges),
      ranges = IRanges(
        start = pmax(
          start(x = overlap_frame$target_ranges),
          start(x = overlap_frame$transcribed_ranges)
        ),
        end = pmin(
          end(x = overlap_frame$target_ranges),
          end(x = overlap_frame$transcribed_ranges)
        )
      )
    )
  rm(overlap_frame)
} else {
  # If target regions are not avaiable, all transcribed GRanges count.
  summary_frame[i, "target_path"] <- NA
  message("Not importing target range annotation.")
  target_ranges <- transcribed_ranges
  summary_frame[i, "target_width_raw"] <-
    summary_frame[i, "transcribed_width"]
  summary_frame[i, "target_number_raw"] <- 0L
  constrained_ranges <- transcribed_ranges
}
summary_frame[i, "target_number_constrained"] <-
  length(x = constrained_ranges)
message(paste0("Number of transcribed target ranges: ", summary_frame[i, "target_number_constrained"]))
summary_frame[i, "target_width_constrained"] <-
  sum(width(x = constrained_ranges))
message(paste0("Cumulative width of transcribed target ranges: ", summary_frame[i, "target_width_constrained"]))

# Read callable loci ------------------------------------------------------


# Import the callable loci BED file produced by the GATK CallableLoci analysis.
summary_frame[i, "callable_loci_path"] <-
  argument_list$callable_loci_path
# Store the sample name in the summary frame.
summary_frame[i, "sample_name"] <-
  gsub(pattern = "_callable_loci.bed",
       replacement = "",
       x = summary_frame[i, "callable_loci_path"])
message(paste0("Processing sample name: ", summary_frame[i, "sample_name"]))
callable_ranges <-
  import(con = summary_frame[i, "callable_loci_path"])
non_callable_ranges <-
  callable_ranges[callable_ranges$name != "CALLABLE", ]
non_callable_ranges$name <-
  as.factor(x = non_callable_ranges$name)
summary_frame[i, "non_callable_number_raw"] <-
  length(x = non_callable_ranges)
message(paste0("Number of non-callable raw ranges: ", summary_frame[i, "non_callable_number_raw"]))
summary_frame[i, "non_callable_width_raw"] <-
  sum(width(x = non_callable_ranges))
message(paste0("Cumulative width of non-callable raw ranges: ", summary_frame[i, "non_callable_width_raw"]))

# Overlap constrained GRanges ---------------------------------------------


# To get accurate non-callable statistics with regards to the target regions
# that are transcribed, non-callable GRanges need overlapping with the
# constrained GRanges.
overlap_frame <-
  mergeByOverlaps(query = non_callable_ranges, constrained_ranges)
# Constrain the GRanges to the minimum overlap.
overlap_ranges <- GRanges(
  seqnames = seqnames(x = overlap_frame$non_callable_ranges),
  ranges = IRanges(start = pmax(
    start(x = overlap_frame$non_callable_ranges),
    start(x = overlap_frame$constrained_ranges)
  ),
  end = pmin(
    end(x = overlap_frame$non_callable_ranges),
    end(x = overlap_frame$constrained_ranges)
  )),
  mapping_status = overlap_frame$name
)
summary_frame[i, "non_callable_number_constrained.TOTAL"] <-
  length(x = overlap_ranges)
message(paste0("Number of non-callable constrained ranges: ",
               summary_frame[i, "non_callable_number_constrained.TOTAL"]))
summary_frame[i, "non_callable_width_constrained.TOTAL"] <-
  sum(width(x = overlap_ranges))
message(paste0(
  "Cumulative width of non-callable constrained ranges: ",
  summary_frame[i, "non_callable_width_constrained.TOTAL"]
))
# summary_frame[i, "non_callable_constrained_fraction.TOTAL"] <-
#   summary_frame[i, "non_callable_width_constrained.TOTAL"] / summary_frame[i, "target_width_constrained"]

# Summarise by mapping status ---------------------------------------------


# Summarise also separately by mapping status.
# Populate the summary frame with columns for each mapping status level,
# regardless of whether it is associated with data or not.
# Use fixed mapping status levels emitted by GATK CallableLoci.
for (level in c(
  "REF_N",
  "NO_COVERAGE",
  "LOW_COVERAGE",
  "EXCESSIVE_COVERAGE",
  "POOR_MAPPING_QUALITY"
)) {
  summary_frame[i, paste("non_callable_number_constrained", level, sep = ".")] <-
    0L
  summary_frame[i, paste("non_callable_width_constrained", level, sep = ".")] <-
    0L
}
rm(level)
if (length(x = overlap_ranges) > 0L) {
  # Count the number of entries for each mapping status level.
  aggregate_frame <-
    as.data.frame(x = table(mcols(x = overlap_ranges)$mapping_status))
  # Assign the result levels (rows) as summary frame columns.
  for (j in seq_len(length.out = nrow(x = aggregate_frame))) {
    summary_frame[i, paste("non_callable_number_constrained", aggregate_frame[j, 1L], sep = ".")] <-
      aggregate_frame[j, 2L]
  }
  # Sum the widths of entries for each mapping_status level.
  aggregate_frame <-
    aggregate.data.frame(
      x = data.frame(width = width(x = overlap_ranges)),
      by = list(mapping_status = mcols(x = overlap_ranges)$mapping_status),
      FUN = "sum"
    )
  # Assign the result levels (rows) as summary frame columns.
  for (j in seq_len(length.out = nrow(x = aggregate_frame))) {
    summary_frame[i, paste("non_callable_width_constrained", aggregate_frame[j, 1L], sep = ".")] <-
      aggregate_frame[j, 2L]
  }
  rm(j, aggregate_frame)
}
rm(overlap_ranges, overlap_frame)

# Annotate non-callable GRanges with target region names ------------------


# Annotate the table of non-callable GRanges with target region names if available.
diagnose_ranges <- NULL
if (!is.null(x = summary_frame[i, "target_path"])) {
  # If the target GRanges are available, merge by overlap with the non-callable GRanges into a new DataFrame.
  overlap_frame <-
    mergeByOverlaps(query = target_ranges, subject = non_callable_ranges)
  # Extract the "non_callable_ranges".
  diagnose_ranges <- overlap_frame[, c("non_callable_ranges")]
  # Annotate with the target_name from the target GRanges object.
  mcols(x = diagnose_ranges)$target_name <-
    overlap_frame[, "name"]
  # Rename the "name" column of the mcols() DataFrame into "mapping_status".
  colnames(x = mcols(x = diagnose_ranges))[colnames(x = mcols(x = diagnose_ranges)) == "name"] <-
    "mapping_status"
  rm(overlap_frame)
} else {
  # Diagnose GRanges without target annotation.
  diagnose_ranges <- non_callable_ranges
}
# To annotate non-callable regions, merge by overlap with the exon GRanges.
overlap_frame <-
  mergeByOverlaps(query = diagnose_ranges, subject = exon_ranges)
# The mergeByOverlaps() function returns a DataFrame of query (diagnose GRanges) and
# subject (exon GRanges) variables, as well as all mcols() variables that were present
# in either GRanges object. To remove these redundant mcols() variables, keep only the
# GRanges objects themselves.
write.table(
  x = overlap_frame[, c("diagnose_ranges", "exon_ranges")],
  file = paste(
    summary_frame[i, "sample_name"],
    "non_callable_loci.tsv",
    sep = "_"
  ),
  quote = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# Group GRanges by non-callable regions -----------------------------------


# Since the mergeByOverlaps() function provides the cartesian product of diagnosis and exon GRanges,
# the table can be rather long and unwieldy. Annotate the diagnosis GRanges with exon GRanges meta
# information before grouping by diagnosis GRanges so that there is a single observation for each problematic region.
# Gene and transcript identifiers and names are turned into comma-separated lists.
message("Annotate the diagnostic ranges.")
overlap_diagnose_ranges <- overlap_frame$diagnose_ranges
overlap_diagnose_frame <- mcols(x = overlap_diagnose_ranges)
overlap_diagnose_frame$gene_id <-
  mcols(x = overlap_frame$exon_ranges)$gene_id
overlap_diagnose_frame$gene_name <-
  mcols(x = overlap_frame$exon_ranges)$gene_name
overlap_diagnose_frame$transcript_id <-
  mcols(x = overlap_frame$exon_ranges)$transcript_id
overlap_diagnose_frame$transcript_name <-
  mcols(x = overlap_frame$exon_ranges)$transcript_name
overlap_diagnose_frame$exon_id <-
  mcols(x = overlap_frame$exon_ranges)$exon_id
mcols(x = overlap_diagnose_ranges) <- overlap_diagnose_frame
# This returns a Grouping object (CompressedManyToOneGrouping) of the IRanges package,
# specifying which groups contain which indices to the original object.
message("Group annotated diagnose ranges by region.")
overlap_diagnose_grouping <-
  as(object = overlap_diagnose_ranges, "Grouping")

grouped_ranges <- unlist(x = GRangesList(lapply(
  X = overlap_diagnose_grouping,
  FUN = function(x) {
    sub_ranges <- overlap_diagnose_ranges[x]
    sub_mcols <- mcols(x = sub_ranges)
    selected_range <- sub_ranges[1L]
    
    mcols(x = selected_range) <- DataFrame(
      "mapping_status" = sub_mcols[1L, "mapping_status"],
      
      "gene_ids" = paste(unique(x = sort(x = sub_mcols$gene_id)), collapse = ","),
      
      "gene_names" = paste(unique(x = sort(x = sub_mcols$gene_name)), collapse = ","),
      
      "transcript_ids" = paste(unique(x = sort(x = sub_mcols$transcript_id)), collapse = ","),
      
      "transcript_names" = paste(unique(x = sort(
        x = sub_mcols$transcript_name
      )), collapse = ","),
      
      "exon_ids" = paste(unique(x = sort(x = sub_mcols$exon_id)), collapse = ",")
    )
    return(selected_range)
  }
)))
write.table(
  x = grouped_ranges,
  file = paste(
    summary_frame[i, "sample_name"],
    "non_callable_regions.tsv",
    sep = "_"
  ),
  quote = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

rm(
  grouped_ranges,
  overlap_diagnose_grouping,
  overlap_diagnose_frame,
  overlap_diagnose_ranges,
  overlap_frame,
  diagnose_ranges,
  callable_ranges,
  non_callable_ranges
)

write.table(
  x = summary_frame,
  file = paste(
    summary_frame[i, "sample_name"],
    "non_callable_summary.tsv",
    sep = "_"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

rm(
  i,
  exon_ranges,
  transcribed_ranges,
  target_ranges,
  constrained_ranges,
  summary_frame,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())
