prepare_promoter_annotation <- function(gtf, sqlite) {
  # Check if GTF file exists
  if (!file.exists(gtf)) {
    stop(paste("GTF file not found:", gtf))
  }
  
  # Create or load TxDb object
  if (!file.exists(sqlite)) {
    txdb <- txdbmaker::makeTxDbFromGFF(
      gtf, 
      format = "gtf", 
      dataSource = "gencode", 
      organism = "Homo sapiens"
    )
    saveDb(txdb, file = sqlite)
  } else {
    txdb <- loadDb(sqlite)
  }
  
  # Prepare promoter annotation
  promoterAnnotation <- proActiv::preparePromoterAnnotation(
    txdb = txdb, 
    species = "Homo_sapiens"
  )
  
  return(promoterAnnotation)
}

process_salmon_data <- function(input_dir, promoterAnnotation, sample_info, 
                               condition_col = "condition", 
                               sample_col = "sample_id",
                               apply_filter = TRUE, 
                               min_count = 10, min_prop = 0.7) {

  # Find quant.sf files
  quant.files <- list.files(input_dir, pattern = "quant.sf", 
                            recursive = TRUE, full.names = TRUE)
  
  if (length(quant.files) == 0) {
    stop("No quant.sf files found in ", input_dir)
  }
  
  quant.dirs <- unique(dirname(quant.files))
  message("Found ", length(quant.dirs), " Salmon quantification directories")
  
  # Import using catchSalmon
  catch <- edgeR::catchSalmon(paths = quant.dirs)
  
  # Prepare count matrix
  if (!"Overdispersion" %in% colnames(catch$annotation)) {
    count_matrix <- catch$counts
  } else {
    count_matrix <- catch$counts / catch$annotation$Overdispersion
  }
  
  # Clean column names
  colnames(count_matrix) <- gsub(".*/", "", quant.dirs)
  
  # Create data frame for merging
  count_df <- data.frame(
    transcriptName = rownames(count_matrix),
    count_matrix,
    stringsAsFactors = FALSE
  )
  
  # Get promoter mapping
  promoterIdMap <- promoterIdMapping(promoterAnnotation)
  
  # Merge with promoter annotation
  transcripts_matched <- sum(count_df$transcriptName %in% promoterIdMap$transcriptName)
  message("Direct transcript ID matches: ", transcripts_matched, "/", nrow(count_df))
  
  merged_data <- merge(count_df, promoterIdMap, by = "transcriptName", all.x = FALSE)
  
  # Identify sample columns
  sample_cols <- intersect(colnames(merged_data), sample_info[[sample_col]])
  
  if (length(sample_cols) == 0) {
    sample_cols <- colnames(count_matrix)
  }
  
  message("Using ", length(sample_cols), " samples for analysis")
  
  # Group by promoter
  merged_dt <- as.data.table(merged_data)
  promoter_data <- merged_dt[,
    lapply(.SD, sum, na.rm = TRUE),
    by = .(promoterId, geneId),
    .SDcols = sample_cols
  ]
  promoter_data <- as.data.frame(promoter_data)
  rownames(promoter_data) <- promoter_data$promoterId
  
  # Create DGEList
  count_matrix_final <- as.matrix(promoter_data[, sample_cols])
  dge_list <- edgeR::DGEList(counts = count_matrix_final)
  
  # Apply filtering if requested
  if (apply_filter) {
    message("Applying edgeR filtering")
    
    if (condition_col %in% colnames(sample_info)) {
      design <- model.matrix(as.formula(paste("~", condition_col)), data = sample_info)
    } else {
      design <- model.matrix(~ 1, data = sample_info)
    }
    
    keep <- edgeR::filterByExpr(dge_list, design = design, 
                               min.count = min_count, min.prop = min_prop)
    message("Filtering: kept ", sum(keep), "/", length(keep), " promoters")
    
    dge_list <- dge_list[keep, , keep.lib.sizes = FALSE]
  }
  
  promoterCounts <- round(dge_list$counts)
  
  # Quality control summary
  cat("\n=== FINAL SUMMARY ===\n")
  cat("Original transcripts:", nrow(count_df), "\n")
  cat("Transcripts with promoter annotation:", nrow(merged_data), "\n")
  cat("Unique promoters after aggregation:", nrow(promoter_data), "\n")
  cat("Final promoters:", nrow(promoterCounts), "\n")
  cat("Samples analyzed:", ncol(promoterCounts), "\n")
  
  cat("\n=== Quality Control Summary ===\n")
  cat("Promoters with zero counts:", sum(rowSums(promoterCounts) == 0), "\n")
  cat("Promoters with counts > 10:", sum(apply(promoterCounts, 1, function(x) any(x > 10))), "\n")
  cat("Median counts per promoter:", median(rowSums(promoterCounts)), "\n")
  cat("Median counts per sample:", median(colSums(promoterCounts)), "\n")
  
  return(promoterCounts)
}

buildSummarizedExperiment_edited <- function(promoterCounts, promoterAnnotation, 
                                     sample_info, 
                                     sample_col = "sample_id", 
                                     condition_col = "condition") {

    # Extract file labels and conditions from sample_info
    fileLabels <- sample_info[[sample_col]]
    condition <- if (condition_col %in% colnames(sample_info)) sample_info[[condition_col]] else NULL
    
    # Normalize promoter counts
    normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
    
    # Get absolute promoter activity
    absolutePromoterActivity <- getAbsolutePromoterActivity(
        normalizedPromoterCounts, 
        promoterAnnotation
    )
    
    # Get gene expression
    geneExpression <- getGeneExpression(absolutePromoterActivity)
    rownames(geneExpression) <- rownames(promoterCounts)
    
    # Get relative promoter activity
    relativePromoterActivity <- getRelativePromoterActivity(
        absolutePromoterActivity, 
        geneExpression
    )
    
    # Create SummarizedExperiment
    result <- SummarizedExperiment(assays = list(
        promoterCounts = promoterCounts,
        normalizedPromoterCounts = normalizedPromoterCounts,
        absolutePromoterActivity = absolutePromoterActivity[, fileLabels, drop = FALSE],
        relativePromoterActivity = relativePromoterActivity[, fileLabels, drop = FALSE],
        geneExpression = geneExpression[, fileLabels, drop = FALSE]
    ))
    
    message('Calculating positions of promoters...')
    
    # Get promoter coordinates and mapping
    promoterCoords <- promoterCoordinates(promoterAnnotation)
    promoterIdMap <- promoterIdMapping(promoterAnnotation)
    
    # Add geneId to coordinates
    promoterCoords$geneId <- promoterIdMap$geneId[match(
        promoterCoords$promoterId, promoterIdMap$promoterId
    )]
    
    # Calculate promoter positions
    promoterCoords <- as.data.table(promoterCoords)
    promoterPosition <- geneId <- strand <- NULL
    promoterCoords[, promoterPosition := ifelse(strand == '+', seq_len(.N), 
                                               rev(seq_len(.N))), by = geneId]
    
    # Filter promoter coordinates to match filtered promoters
    filtered_promoter_ids <- absolutePromoterActivity$promoterId
    promoterCoords_filtered <- promoterCoords[promoterCoords$promoterId %in% filtered_promoter_ids, ]
    
    # Build row data
    rowData(result) <- data.frame(
        absolutePromoterActivity[, c('promoterId', 'geneId')], 
        promoterCoords_filtered[match(absolutePromoterActivity$promoterId, 
                                     promoterCoords_filtered$promoterId), 
                               c("seqnames", "start", "strand", 
                                 "internalPromoter", "promoterPosition")]
    )
    
    # Add transcript information
    transcriptByPromoter <- split(promoterIdMap$transcriptName, 
                                 promoterIdMap$promoterId)
    rowData(result)$txId <- transcriptByPromoter[match(
        rowData(result)$promoterId, 
        names(transcriptByPromoter)
    )]
    
    # Summarize across conditions if provided
    if (!is.null(condition)) {
        if (length(condition) != length(fileLabels)) {
            warning('Condition length does not match sample length. 
                    Returning results not summarized across conditions.')
            return(result)
        } else {
            result <- summarizeAcrossCondition(result, condition) 
        }
    }
    
    return(result)
}

# Helper function to summarize results across condition
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom SummarizedExperiment rowData colData 'colData<-' assays
summarizeAcrossCondition <- function(result, condition) {
        if (any(make.names(condition) != condition)) {
            warning("Condition is modified to be syntactically valid")
            condition <- make.names(condition)
        }
        colData(result) <- DataFrame(sampleName = colnames(result), 
                                    condition=condition)
        colnames(result) <- colData(result)$sampleName
        message('Summarising expression and activity across conditions...')
        for (group in unique(condition)) {
            rowData(result)[,paste0(group, '.mean')] <- 
                rowMeans(assays(result[,colData(result)$condition==group])$abs) 
            rowData(result)[,paste0(group, '.gene.mean')] <- 
                rowMeans(assays(result[,colData(result)$condition==group])$gene)
        }
        rowData(result) <- categorizePromoters(rowData(result), condition)
    return(result)
}

# Helper function to categorize promoters
#' @importFrom rlang .data
#' @importFrom dplyr as_tibble '%>%' slice_max
categorizePromoters <- function(rdata, condition) {
    rdata <- as_tibble(rdata)
    for (group in unique(condition)) {
        message(paste0('Categorizing ', group, ' promoters...'))
        mean <- paste0(group, '.mean')
        class <- paste0(group, '.class')
        
        max_rows <- rdata %>%
            group_by(.data$geneId) %>%
            slice_max(!!as.name(mean), with_ties = FALSE)
        rdata[[class]] <- ifelse(rdata[[mean]] < 0.25, 'Inactive', 'Minor')
        rdata[[class]][match(max_rows$promoterId, rdata$promoterId)] <- "Major"
        rdata[[class]][which(rdata[[mean]] < 0.25)] <- "Inactive"
        rdata[[class]][which(rdata$internalPromoter)] <- NA
        rdata[[class]][which(is.na(rdata$internalPromoter))] <- NA
    }
    return(rdata)
}
