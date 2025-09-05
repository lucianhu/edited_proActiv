buildSummarizedExperiment_edited <- function(promoterCounts, promoterAnnotation, 
                                     sample_info, 
                                     sample_col = "sample_id", 
                                     condition_col = "condition") {

    # Extract file labels and conditions from sample_info
    fileLabels <- sample_info[[sample_col]]
    condition <- if (condition_col %in% colnames(sample_info)) sample_info[[condition_col]] else NULL
    
   normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
   absolutePromoterActivity <- getAbsolutePromoterActivity(
                                                    normalizedPromoterCounts, 
                                                    promoterAnnotation)
   geneExpression <- getGeneExpression(absolutePromoterActivity)
   rownames(geneExpression) <- rownames(promoterCounts)
   relativePromoterActivity <- getRelativePromoterActivity(
                                                    absolutePromoterActivity, 
                                                    geneExpression)
  result <- SummarizedExperiment(assays = list(
          promoterCounts = promoterCounts,
          normalizedPromoterCounts = normalizedPromoterCounts,
          absolutePromoterActivity = absolutePromoterActivity[, 
                                                  fileLabels, drop = FALSE],
          relativePromoterActivity = relativePromoterActivity[, 
                                                  fileLabels, drop = FALSE],
          geneExpression = geneExpression[, fileLabels, drop = FALSE]))
    
    message('Calculating positions of promoters...')
    promoterCoordinates <- promoterCoordinates(promoterAnnotation)
    promoterIdMapping <- promoterIdMapping(promoterAnnotation)
    promoterCoordinates$geneId <- promoterIdMapping$geneId[match(
                promoterCoordinates$promoterId, promoterIdMapping$promoterId)]
    promoterCoordinates <- as.data.table(promoterCoordinates)
    promoterPosition <- geneId <- strand <- NULL
    promoterCoordinates[, promoterPosition := ifelse(strand == '+', seq_len(.N), 
                                                rev(seq_len(.N))), by=geneId]
  
    # Filter promoter coordinates to match filtered promoters
    filtered_promoter_ids <- absolutePromoterActivity$promoterId
    promoterCoords_filtered <- promoterCoordinates[promoterCoordinates$promoterId %in% filtered_promoter_ids, ]
    
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
