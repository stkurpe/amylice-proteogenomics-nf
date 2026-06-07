#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_file <- ifelse(length(args) >= 1, args[1], "tests/fixtures/amylogram_reference.json")

suppressPackageStartupMessages(library(AmyloGram))
suppressPackageStartupMessages(library(biogram))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(jsonlite))

model <- AmyloGram_model

feature_matrix <- function(seq_chars) {
  seqs_m <- tolower(matrix(seq_chars, nrow = 1))
  gl <- do.call(rbind, lapply(1L:nrow(seqs_m), function(i) {
    res <- do.call(
      rbind,
      strsplit(decode_ngrams(seq2ngrams(seqs_m[i, ][!is.na(seqs_m[i, ])], 6, seqinr::a()[-1])), "")
    )
    cbind(res, id = paste0("P", rep(i, nrow(res))))
  }))

  bitrigrams <- as.matrix(count_multigrams(
    ns = c(1, rep(2, 4), rep(3, 3)),
    ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
    seq = degenerate(gl[, -7], model[["enc"]]),
    u = as.character(1L:length(model[["enc"]]))
  ))
  test_ngrams <- bitrigrams > 0
  storage.mode(test_ngrams) <- "integer"
  test_ngrams
}

examples <- list(
  VQIVYK = strsplit("VQIVYK", "")[[1]],
  ACDEFGHIK = strsplit("ACDEFGHIK", "")[[1]],
  GGGGGG = strsplit("GGGGGG", "")[[1]],
  KLVFFA = strsplit("KLVFFA", "")[[1]],
  NQSTDE = strsplit("NQSTDE", "")[[1]]
)

example_rows <- lapply(names(examples), function(example_id) {
  seq_chars <- examples[[example_id]]
  fm <- feature_matrix(seq_chars)
  pred <- predict(model, setNames(list(seq_chars), example_id))
  list(
    id = example_id,
    sequence = paste(seq_chars, collapse = ""),
    probability = pred$Probability[[1]],
    all_feature_names = colnames(fm),
    all_feature_values = as.integer(fm[1, ]),
    selected_feature_values = as.integer(fm[1, model[["imp_features"]]])
  )
})

payload <- list(
  source = "R AmyloGram 1.1",
  enc = model[["enc"]],
  imp_features = model[["imp_features"]],
  independent_variable_names = model[["rf"]][["forest"]][["independent.variable.names"]],
  forest = list(
    num_trees = model[["rf"]][["forest"]][["num.trees"]],
    child_node_ids = model[["rf"]][["forest"]][["child.nodeIDs"]],
    split_var_ids = model[["rf"]][["forest"]][["split.varIDs"]],
    split_values = model[["rf"]][["forest"]][["split.values"]],
    terminal_class_counts = model[["rf"]][["forest"]][["terminal.class.counts"]],
    levels = model[["rf"]][["forest"]][["levels"]]
  ),
  examples = example_rows
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write(toJSON(payload, auto_unbox = TRUE, pretty = TRUE), output_file)
cat("Wrote", output_file, "\n")
