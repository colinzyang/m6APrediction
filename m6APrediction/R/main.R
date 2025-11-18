#' DNA Sequence Encoding
#'
#' @description
#' Converts a vector of DNA sequences into a data frame of factors suitable for machine learning models.
#' Each position in the sequence becomes a column with factor levels A, T, C, G.
#'
#' @param dna_strings A character vector of DNA sequences. All sequences must be of the same length.
#'
#' @return A data frame where each column represents a nucleotide position.
#' @export
#'
#' @examples
#' dna_encoding(c("ATCGG", "TTAGC"))
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Multiple Sample Prediction for m6A Sites
#'
#' @description
#' Predicts the probability and status of m6A modification for multiple sites provided in a data frame.
#'
#' @param ml_fit A trained Random Forest model object.
#' @param feature_df A data frame containing the genomic features. Must contain columns:
#' "gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
#' "evolutionary_conservation", "DNA_5mer".
#' @param positive_threshold A numeric value between 0 and 1. Cutoff for "Positive" classification. Default is 0.5.
#'
#' @return The input data frame with two additional columns:
#' \item{predicted_m6A_prob}{The probability of the site being m6A positive.}
#' \item{predicted_m6A_status}{The classification ("Positive" or "Negative").}
#' @import randomForest
#' @export
#'
#' @examples
#' \dontrun{
#' # Load internal model and example data
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' input_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' 
#' # Run prediction
#' results <- prediction_multiple(rf_model, input_df)
#' head(results)
#' }
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  
  # Set factor levels for categorical variables to match training data
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  
  # Encode DNA sequences
  dna_encoded <- dna_encoding(feature_df$DNA_5mer)
  feature_df <- cbind(feature_df, dna_encoded)
  
  # Make predictions using the model
  predicted_probs <- predict(ml_fit, newdata = feature_df, type = "prob")
  positive_probs <- predicted_probs[, "Positive"]
  
  # Apply threshold to determine status
  predicted_status <- ifelse(positive_probs > positive_threshold, "Positive", "Negative")
  
  # Add prediction results to the original dataframe
  feature_df$predicted_m6A_prob <- positive_probs
  feature_df$predicted_m6A_status <- predicted_status
  
  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}

#' Single Sample Prediction for m6A Sites
#'
#' @description
#' Predicts the probability and status of m6A modification for a single site based on individual features.
#'
#' @param ml_fit A trained Random Forest model object.
#' @param gc_content Numeric. GC content of the site.
#' @param RNA_type Character. Type of RNA (e.g., "mRNA", "lincRNA").
#' @param RNA_region Character. Region of RNA (e.g., "CDS", "3'UTR").
#' @param exon_length Numeric. Length of the exon.
#' @param distance_to_junction Numeric. Distance to the nearest splice junction.
#' @param evolutionary_conservation Numeric. Conservation score.
#' @param DNA_5mer Character. 5-nucleotide DNA sequence (e.g., "GGACA").
#' @param positive_threshold Numeric. Cutoff for "Positive" classification. Default is 0.5.
#'
#' @return A named vector containing:
#' \item{predicted_m6A_prob}{The probability score.}
#' \item{predicted_m6A_status}{The classification status.}
#' @export
#'
#' @examples
#' \dontrun{
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(rf_model, 0.5, "mRNA", "CDS", 10, 8, 0.5, "GGACA")
#' }
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  
  # Convert inputs to proper factors with correct levels
  RNA_type <- factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  RNA_region <- factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  
  # Create a single-row data frame with the input values
  single_row_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  
  # Use prediction_multiple function to get the prediction
  result_df <- prediction_multiple(ml_fit, single_row_df, positive_threshold)
  
  # Extract the prediction results
  predicted_prob <- result_df$predicted_m6A_prob[1]
  predicted_status <- result_df$predicted_m6A_status[1]
  
  # Create named vector to return
  returned_vector <- c(predicted_prob, predicted_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}