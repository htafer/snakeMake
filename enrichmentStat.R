#!/usr/bin/env Rscript
# ===========================================================================
# enrichmentStat.R
# ===========================================================================
#
# Description:
#   Performs gene set enrichment analysis for either GO terms or custom
#   annotations (PFAM, KEGG, etc.). Identifies significantly enriched
#   terms in a gene list compared to a background set.
#
# Usage:
#   Rscript enrichmentStat.R -b background.txt -d genes.txt -t <GO|PFAM|KEGG> [options]
#
# Required Arguments:
#   -b, --background  Background data file with gene annotations
#                     Format: tab-delimited with headers 'IDs' and 'ANN'
#   -d, --geneFile    File containing genes to test for enrichment
#                     Format: tab-delimited with gene IDs
#   -t, --type       Analysis type: 'GO' or other annotation (PFAM/KEGG)
#
# Optional Arguments:
#   -p, --pValue     P-value threshold for significance (default: 0.1)
#   -D, --dictionary File mapping IDs to descriptions (optional)
#
# Output:
#   Creates a CSV file with enrichment results: {input}.{type}.csv
#
# Dependencies:
#   R packages: optparse, AnnotationForge, GOstats, GSEABase, data.table
#
# Author: htafer
# Last Updated: 2025-07-28
# ===========================================================================

# Load required packages with error handling
required_packages <- c("optparse", "AnnotationForge", "GOstats", "GSEABase")

# Function to check and load required packages
load_required_packages <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("Package", pkg, "is required but not installed."))
        }
    }
    suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
}

# Load base packages
load_required_packages("optparse")

# Define command line arguments
option_list <- list(
    make_option(
        c("-b", "--background"),
        help = paste(
            "Background data file containing gene annotations.",
            "Must have headers 'IDs' and 'ANN'."
        ),
        type = "character"
    ),
    make_option(
        c("-d", "--geneFile"),
        help = "File containing genes for enrichment analysis",
        type = "character"
    ),
    make_option(
        c("-t", "--type"),
        help = "Analysis type: 'GO' or other annotation type",
        type = "character"
    ),
    make_option(
        c("-p", "--pValue"),
        default = 0.1,
        help = "P-value threshold for significance (default: 0.1)",
        type = "double"
    ),
    make_option(
        c("-D", "--dictionary"),
        default = "NA",
        help = "Optional: dictionary file for ID to description mapping",
        type = "character"
    )
)

# Parse and validate command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$background) || is.null(opt$geneFile) || is.null(opt$type)) {
    stop("Error: --background, --geneFile, and --type arguments are required")
}

# Validate file existence
for (file in c(opt$background, opt$geneFile)) {
    if (!file.exists(file)) {
        stop(paste("Error: File does not exist:", file))
    }
}

# Validate p-value
if (opt$pValue <= 0 || opt$pValue > 1) {
    stop("Error: P-value must be between 0 and 1")
}

# Load additional required packages based on analysis type
if (opt$type == "GO") {
    load_required_packages(c("AnnotationForge", "GOstats", "GSEABase", "xtable"))
} else {
    load_required_packages(c("data.table", "AnnotationForge", "GOstats", "GSEABase"))
}

# Helper Functions
# ===========================================================================

#' Get genes corresponding to enriched terms
#'
#' @param term_id The ID of the enriched term to find genes for
#' @return Comma-separated string of gene IDs
#' @note Reads background and gene files each time - could be optimized
get_corresponding_genes <- function(term_id) {
    tryCatch({
        # Read annotation data
        annotations <- read.table(
            opt$background,
            header = TRUE,
            stringsAsFactors = FALSE
        )
        
        # Read gene list
        gene_data <- read.table(
            opt$geneFile,
            header = TRUE,
            row.names = 1,
            stringsAsFactors = FALSE
        )
        
        # Find genes annotated with the term
        matching_genes <- row.names(gene_data)[
            row.names(gene_data) %in% annotations$IDs[annotations$ANN %in% term_id]
        ]
        
        return(toString(matching_genes))
    }, error = function(e) {
        warning(paste("Error getting genes for term", term_id, ":", e$message))
        return(NA)
    })
}

#' Perform GO enrichment analysis
#'
#' @return Data frame with enrichment results
perform_go_enrichment <- function() {
    message("Performing GO enrichment analysis...")
    
    # Read and prepare GO data
    go_data <- read.table(opt$background, header = TRUE, stringsAsFactors = FALSE)
    
    # Create GO frame
    go_frame <- GOFrame(go_data, organism = "Exophiala dermatitidis")
    go_all_frame <- GOAllFrame(go_frame)
    gene_set_collection <- GeneSetCollection(go_all_frame, setType = GOCollection())
    
    # Get universe of genes
    universe <- getGOFrameData(go_all_frame)
    universe <- unique(universe$gene_id)
    
    # Read and prepare gene list
    diff_data <- read.table(opt$geneFile, header = TRUE, row.names = 1)
    gene_list <- row.names(diff_data)
    gene_list <- intersect(universe, gene_list)
    
    if (length(gene_list) == 0) {
        stop("No genes in input list match the universe")
    }
    
    # Set up and perform enrichment test
    params <- GSEAGOHyperGParams(
        name = "GO Enrichment Analysis",
        geneSetCollection = gene_set_collection,
        geneIds = gene_list,
        universeGeneIds = universe,
        ontology = c("BP", "CC", "MF"),
        pvalueCutoff = 1,
        conditional = FALSE,
        testDirection = "over"
    )
    
    # Run analysis
    results <- hyperGTest(params)
    results_summary <- summary(results)
    
    # Add FDR correction
    results_summary$fdr <- p.adjust(results_summary$Pvalue, method = "fdr")
    
    # Filter by significance
    significant_results <- results_summary[results_summary$fdr < opt$pValue, ]
    
    if (nrow(significant_results) == 0) {
        warning("No significantly enriched GO terms found")
    }
    
    return(significant_results)
}

#' Perform enrichment analysis for other annotation types (PFAM, KEGG)
#'
#' @return Data frame with enrichment results
perform_other_enrichment <- function() {
    message(paste("Performing", opt$type, "enrichment analysis..."))
    
    tryCatch({
        # Read and prepare annotation data
        annotations <- read.table(
            opt$background,
            header = TRUE,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        
        # Create gene sets
        gene_sets <- split(annotations$IDs, annotations$ANN)
        gene_set_collection <- GeneSetCollection(Map(function(pid, gids) {
            GeneSet(
                gids,
                setName = pid,
                collectionType = PfamCollection(pid)
            )
        }, names(gene_sets), gene_sets))
        
        # Get universe of genes
        universe <- unique(annotations$IDs)
        
        # Read and prepare gene list
        diff_data <- read.table(
            opt$geneFile,
            header = TRUE,
            row.names = 1,
            stringsAsFactors = FALSE
        )
        gene_list <- row.names(diff_data)
        gene_list <- intersect(universe, gene_list)
        
        if (length(gene_list) == 0) {
            stop("No genes in input list match the universe")
        }
        
        # Set up and perform enrichment test
        params <- GSEAKEGGHyperGParams(
            name = paste(opt$type, "Enrichment Analysis"),
            geneSetCollection = gene_set_collection,
            geneIds = gene_list,
            universeGeneIds = universe,
            testDirection = "over",
            pvalueCutoff = 1
        )
        
        # Run analysis
        results <- hyperGTest(params)
        results_summary <- summary(results)
        
        # Add FDR correction
        results_summary$fdr <- p.adjust(results_summary$Pvalue, method = "fdr")
        
        # Filter by significance
        significant_results <- results_summary[results_summary$fdr < opt$pValue, ]
        
        if (nrow(significant_results) == 0) {
            warning("No significantly enriched terms found")
            return(NULL)
        }
        
        # Add descriptions if dictionary provided
        if (opt$dictionary != "NA" && file.exists(opt$dictionary)) {
            descriptions <- fread(opt$dictionary, header = FALSE)
            significant_results$desc <- descriptions$V2[
                match(significant_results$KEGGID, descriptions$V1)
            ]
        }
        
        # Add corresponding genes
        significant_results$genes <- sapply(
            significant_results$KEGGID,
            get_corresponding_genes
        )
        
        return(significant_results)
    }, error = function(e) {
        stop(paste("Error in enrichment analysis:", e$message))
    })
}

# Main Execution
# ===========================================================================

# Perform analysis based on type
results <- if (opt$type == "GO") {
    perform_go_enrichment()
} else {
    perform_other_enrichment()
}

# Save results if any found
if (!is.null(results) && nrow(results) > 0) {
    output_file <- paste(opt$geneFile, opt$type, "csv", sep = ".")
    write.table(
        results,
        file = output_file,
        sep = ",",
        row.names = TRUE,
        col.names = TRUE,
        quote = TRUE
    )
    message(paste("Results written to:", output_file))
} else {
    warning("No significant results found. No output file created.")
}

