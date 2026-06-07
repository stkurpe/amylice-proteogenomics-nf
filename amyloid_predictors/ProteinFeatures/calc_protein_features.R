#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript calc_protein_features.R <input.fasta> <output.csv> [window] [pH]", call. = FALSE)
}

input_fasta <- args[1]
output_csv <- args[2]
window <- ifelse(length(args) >= 3, as.integer(args[3]), 7L)
pH <- ifelse(length(args) >= 4, as.numeric(args[4]), 7.0)

required_packages <- c("Peptides", "protr", "Biostrings", "zoo")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    sprintf(
      "Missing R package(s): %s. Install with install.packages(c(\"Peptides\", \"protr\", \"zoo\")) and BiocManager::install(\"Biostrings\").",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(Peptides))
suppressPackageStartupMessages(library(protr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(zoo))

data(AAindex, package = "protr")

alpha_index <- "CHOP780201"
beta_index <- "CHOP780202"

clean_seq <- function(seq) {
  seq <- toupper(as.character(seq))
  gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
}

get_aaindex_scale <- function(index_id) {
  aa_order <- c(
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )
  if (is.data.frame(AAindex)) {
    hit <- AAindex[AAindex$AccNo == index_id, aa_order, drop = FALSE]
    if (nrow(hit) != 1) stop(sprintf("AAindex ID not found: %s", index_id), call. = FALSE)
    scale <- as.numeric(hit[1, aa_order])
    names(scale) <- aa_order
    return(scale)
  }

  scale <- AAindex[[index_id]]$I[aa_order]
  if (is.null(scale)) stop(sprintf("AAindex ID not found or unsupported structure: %s", index_id), call. = FALSE)
  scale
}

sliding_mean <- function(values, window = 7L) {
  if (length(values) < window) return(rep(NA_real_, length(values)))
  zoo::rollapply(
    values,
    width = window,
    FUN = mean,
    align = "center",
    fill = NA_real_
  )
}

mean_or_na <- function(values) {
  if (!length(values) || all(is.na(values))) return(NA_real_)
  mean(values, na.rm = TRUE)
}

max_or_na <- function(values) {
  if (!length(values) || all(is.na(values))) return(NA_real_)
  max(values, na.rm = TRUE)
}

sd_or_na <- function(values) {
  if (!length(values) || sum(!is.na(values)) < 2) return(NA_real_)
  sd(values, na.rm = TRUE)
}

calc_features <- function(seq, window = 7L, pH = 7.0) {
  seq <- clean_seq(seq)
  aa <- strsplit(seq, "", fixed = TRUE)[[1]]
  n <- length(aa)

  if (n == 0) {
    return(data.frame(
      protein_length = NA_integer_,
      KD_mean = NA_real_,
      KD_max = NA_real_,
      KD_sd = NA_real_,
      charge = NA_real_,
      charge_density = NA_real_,
      aromaticity = NA_real_,
      alpha_propensity = NA_real_,
      beta_propensity = NA_real_,
      chameleon_score = NA_real_,
      pI = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  kd_values <- vapply(
    aa,
    function(x) Peptides::hydrophobicity(x, scale = "KyteDoolittle"),
    numeric(1)
  )
  kd_window <- sliding_mean(kd_values, window)

  alpha_scale <- get_aaindex_scale(alpha_index)
  beta_scale <- get_aaindex_scale(beta_index)

  alpha_values <- alpha_scale[aa]
  beta_values <- beta_scale[aa]

  alpha_window <- sliding_mean(alpha_values, window)
  beta_window <- sliding_mean(beta_values, window)

  denom <- alpha_window + beta_window
  chameleon_window <- ifelse(
    is.na(denom) | denom == 0,
    NA_real_,
    1 - abs(alpha_window - beta_window) / denom
  )

  net_charge <- Peptides::charge(seq, pH = pH)

  data.frame(
    protein_length = n,
    KD_mean = mean_or_na(kd_window),
    KD_max = max_or_na(kd_window),
    KD_sd = sd_or_na(kd_window),
    charge = net_charge,
    charge_density = net_charge / n,
    aromaticity = Peptides::aIndex(seq),
    alpha_propensity = mean_or_na(alpha_window),
    beta_propensity = mean_or_na(beta_window),
    chameleon_score = mean_or_na(chameleon_window),
    pI = Peptides::pI(seq),
    stringsAsFactors = FALSE
  )
}

seqs <- Biostrings::readAAStringSet(input_fasta)
features <- do.call(rbind, lapply(as.character(seqs), calc_features, window = window, pH = pH))
features$protein_id <- names(seqs)
features <- features[, c(
  "protein_id",
  "protein_length",
  "KD_mean",
  "KD_max",
  "KD_sd",
  "charge",
  "charge_density",
  "aromaticity",
  "alpha_propensity",
  "beta_propensity",
  "chameleon_score",
  "pI"
)]

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(features, output_csv, row.names = FALSE, na = "")
