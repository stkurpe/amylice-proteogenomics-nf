#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript predict_amylogram_fast.R <input.fasta> <output.csv> [chunk_dir] [summary.tsv]", call. = FALSE)
}

input_file <- args[1]
output_file <- args[2]
chunk_dir <- ifelse(length(args) >= 3, args[3], dirname(output_file))
summary_file <- ifelse(length(args) >= 4, args[4], file.path(dirname(output_file), "amylogram_fast_summary.tsv"))

suppressPackageStartupMessages(library(AmyloGram))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(parallel))

timestamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
log_msg <- function(...) cat(sprintf("[%s] %s\n", timestamp(), sprintf(...)))

get_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1) return(default)
  parsed
}

get_env_bool <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(default)
  value %in% c("1", "true", "yes", "y")
}

safe_id <- function(x) {
  x <- gsub("[\t\r\n]", " ", x)
  trimws(x)
}

predict_batch <- function(seqs) {
  input_list <- strsplit(seqs, "", fixed = TRUE)
  names(input_list) <- sprintf("seq_%05d", seq_along(seqs))
  pred <- predict(AmyloGram_model, input_list)
  probabilities <- as.numeric(pred[["Probability"]])
  data.frame(
    AmyloGram_Prob = probabilities,
    AmyloGram_Pred = ifelse(probabilities > 0.5, "AMYLOID", "Non-Amyloid"),
    Error_Message = "",
    stringsAsFactors = FALSE
  )
}

predict_resilient <- function(seqs, max_depth, depth = 0) {
  result <- tryCatch(
    predict_batch(seqs),
    error = function(e) e
  )

  if (!inherits(result, "error")) {
    return(result)
  }

  if (length(seqs) > 1 && depth < max_depth) {
    split_at <- ceiling(length(seqs) / 2)
    left <- predict_resilient(seqs[seq_len(split_at)], max_depth, depth + 1)
    right <- predict_resilient(seqs[(split_at + 1):length(seqs)], max_depth, depth + 1)
    return(rbind(left, right))
  }

  data.frame(
    AmyloGram_Prob = NA_real_,
    AmyloGram_Pred = "ERROR",
    Error_Message = conditionMessage(result),
    stringsAsFactors = FALSE
  )
}

required_chunk_columns <- c("Unique_Index", "AmyloGram_Prob", "AmyloGram_Pred", "Error_Message")

empty_chunk_stats <- function() {
  data.frame(
    chunk_id = integer(),
    status = character(),
    seconds = numeric(),
    rows = integer(),
    errors = integer(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

chunk_stat <- function(chunk_id, status, seconds, rows, errors, message = "") {
  data.frame(
    chunk_id = chunk_id,
    status = status,
    seconds = seconds,
    rows = rows,
    errors = errors,
    message = message,
    stringsAsFactors = FALSE
  )
}

read_chunk <- function(path) {
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE),
    error = function(e) e
  )
}

is_valid_chunk <- function(path, expected_rows) {
  if (!file.exists(path) || file.info(path)$size == 0) return(FALSE)
  df <- read_chunk(path)
  if (inherits(df, "error")) return(FALSE)
  if (!all(required_chunk_columns %in% names(df))) return(FALSE)
  if (nrow(df) != expected_rows) return(FALSE)
  if (!identical(as.integer(df$Unique_Index), seq_len(expected_rows))) return(FALSE)
  TRUE
}

predict_one_chunk <- function(chunk_id, seqs, out_path, max_depth) {
  chunk_start <- Sys.time()
  if (is_valid_chunk(out_path, length(seqs))) {
    existing <- read_chunk(out_path)
    return(chunk_stat(
      chunk_id = chunk_id,
      status = "SKIPPED",
      seconds = 0,
      rows = nrow(existing),
      errors = sum(existing$AmyloGram_Pred == "ERROR", na.rm = TRUE)
    ))
  }
  if (file.exists(out_path)) {
    unlink(out_path)
  }

  prediction <- predict_resilient(seqs, max_depth = max_depth)
  if (nrow(prediction) != length(seqs)) {
    stop(sprintf("Chunk %d row mismatch: got %d expected %d", chunk_id, nrow(prediction), length(seqs)), call. = FALSE)
  }

  result <- data.frame(
    Unique_Index = seq_along(seqs),
    prediction,
    stringsAsFactors = FALSE
  )

  tmp_path <- sprintf("%s.tmp.%s.%d", out_path, Sys.getpid(), chunk_id)
  on.exit(unlink(tmp_path), add = TRUE)
  write.csv(result, tmp_path, row.names = FALSE, quote = TRUE)
  if (!file.rename(tmp_path, out_path)) {
    file.copy(tmp_path, out_path, overwrite = TRUE)
    unlink(tmp_path)
  }
  if (!is_valid_chunk(out_path, length(seqs))) {
    stop(sprintf("Chunk %d output failed validation after write: %s", chunk_id, out_path), call. = FALSE)
  }
  elapsed <- as.numeric(difftime(Sys.time(), chunk_start, units = "secs"))
  chunk_stat(
    chunk_id = chunk_id,
    status = "DONE",
    seconds = elapsed,
    rows = nrow(result),
    errors = sum(result$AmyloGram_Pred == "ERROR", na.rm = TRUE)
  )
}

predict_one_chunk_safe <- function(chunk_id, seqs, out_path, max_depth) {
  tryCatch(
    predict_one_chunk(chunk_id, seqs, out_path, max_depth),
    error = function(e) {
      chunk_stat(
        chunk_id = chunk_id,
        status = "ERROR",
        seconds = 0,
        rows = 0,
        errors = length(seqs),
        message = conditionMessage(e)
      )
    }
  )
}

run_chunk_ids <- function(ids_to_run, run_cores) {
  if (!length(ids_to_run)) return(empty_chunk_stats())
  worker <- function(i) {
    idx <- which(chunk_ids == i)
    predict_one_chunk_safe(i, unique_seqs[idx], chunk_paths[i], max_depth = max_depth)
  }
  results <- if (run_cores > 1 && length(ids_to_run) > 1) {
    mclapply(ids_to_run, worker, mc.cores = min(run_cores, length(ids_to_run)), mc.preschedule = FALSE)
  } else {
    lapply(ids_to_run, worker)
  }
  rows <- lapply(seq_along(results), function(i) {
    result <- results[[i]]
    if (inherits(result, "try-error") || inherits(result, "error")) {
      return(chunk_stat(
        chunk_id = ids_to_run[i],
        status = "ERROR",
        seconds = 0,
        rows = 0,
        errors = 0,
        message = as.character(result)
      ))
    }
    result
  })
  do.call(rbind, rows)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

detected_cores <- max(1, detectCores())
cores <- get_env_int("AMYLOGRAM_CORES", max(1, detected_cores - 1))
cores <- min(cores, detected_cores)
chunk_size <- get_env_int("AMYLOGRAM_CHUNK_SIZE", 200)
max_depth <- get_env_int("AMYLOGRAM_RETRY_SPLIT_DEPTH", 10)
sequential_retries <- get_env_int("AMYLOGRAM_SEQUENTIAL_RETRIES", 3)
allow_errors <- get_env_bool("AMYLOGRAM_ALLOW_ERRORS", FALSE)

log_msg("FAST AmyloGram mode")
log_msg("input=%s", input_file)
log_msg("output=%s", output_file)
log_msg("chunk_dir=%s", chunk_dir)
log_msg(
  "cores=%d chunk_size=%d retry_split_depth=%d sequential_retries=%d allow_errors=%s",
  cores,
  chunk_size,
  max_depth,
  sequential_retries,
  allow_errors
)

read_start <- Sys.time()
seqs_list <- read.fasta(file = input_file, seqtype = "AA", as.string = TRUE)
ids <- safe_id(names(seqs_list))
raw_seqs <- toupper(unlist(seqs_list, use.names = FALSE))
clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", raw_seqs)
valid_mask <- nchar(clean_seqs) >= 6

valid_ids <- ids[valid_mask]
valid_seqs <- clean_seqs[valid_mask]
invalid_count <- length(raw_seqs) - length(valid_seqs)

if (!length(valid_seqs)) stop("No valid sequences found after cleaning/filtering", call. = FALSE)

unique_seqs <- unique(valid_seqs)
unique_map <- match(valid_seqs, unique_seqs)

chunk_ids <- ceiling(seq_along(unique_seqs) / chunk_size)
chunk_count <- max(chunk_ids)
chunk_paths <- file.path(chunk_dir, sprintf("amylogram_chunk_%04d.csv", seq_len(chunk_count)))

log_msg(
  "records=%d valid=%d invalid_or_short=%d unique=%d duplicates_saved=%d chunks=%d",
  length(raw_seqs),
  length(valid_seqs),
  invalid_count,
  length(unique_seqs),
  length(valid_seqs) - length(unique_seqs),
  chunk_count
)
log_msg("read_and_prepare_seconds=%.2f", as.numeric(difftime(Sys.time(), read_start, units = "secs")))

run_start <- Sys.time()
expected_chunk_rows <- as.integer(tabulate(chunk_ids, nbins = chunk_count))

chunk_stats <- run_chunk_ids(seq_len(chunk_count), cores)

invalid_chunks <- function() {
  which(!vapply(seq_len(chunk_count), function(i) {
    is_valid_chunk(chunk_paths[i], expected_chunk_rows[i])
  }, logical(1)))
}

missing_or_invalid <- invalid_chunks()
if (length(missing_or_invalid)) {
  log_msg(
    "chunk_retry_initial missing_or_invalid=%d chunks=%s",
    length(missing_or_invalid),
    paste(missing_or_invalid, collapse = ",")
  )
}

retry_round <- 1L
while (length(missing_or_invalid) && retry_round <= sequential_retries) {
  log_msg(
    "chunk_retry round=%d mode=sequential count=%d chunks=%s",
    retry_round,
    length(missing_or_invalid),
    paste(missing_or_invalid, collapse = ",")
  )
  retry_stats <- run_chunk_ids(missing_or_invalid, 1L)
  chunk_stats <- rbind(chunk_stats, retry_stats)
  missing_or_invalid <- invalid_chunks()
  retry_round <- retry_round + 1L
}

if (length(missing_or_invalid)) {
  log_msg(
    "chunk_retry_failed missing_or_invalid=%d chunks=%s",
    length(missing_or_invalid),
    paste(missing_or_invalid, collapse = ",")
  )
}

log_msg("chunk_status done=%d skipped=%d errors=%d",
        sum(chunk_stats$status == "DONE"),
        sum(chunk_stats$status == "SKIPPED"),
        sum(chunk_stats$errors))

chunk_results <- lapply(seq_len(chunk_count), function(i) {
  path <- chunk_paths[i]
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(sprintf("Missing chunk output: %s", path), call. = FALSE)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$Global_Unique_Index <- df$Unique_Index + ((i - 1) * chunk_size)
  df
})

unique_df <- do.call(rbind, chunk_results)
unique_df <- unique_df[order(unique_df$Global_Unique_Index), ]

if (nrow(unique_df) != length(unique_seqs)) {
  stop(sprintf("Unique prediction count mismatch: got %d expected %d", nrow(unique_df), length(unique_seqs)), call. = FALSE)
}

expanded <- data.frame(
  Sequence_ID = valid_ids,
  AmyloGram_Prob = unique_df$AmyloGram_Prob[unique_map],
  AmyloGram_Pred = unique_df$AmyloGram_Pred[unique_map],
  Unique_Index = unique_map,
  Error_Message = unique_df$Error_Message[unique_map],
  stringsAsFactors = FALSE
)

write.csv(expanded, output_file, row.names = FALSE, quote = TRUE)
elapsed_mins <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))
error_count <- sum(expanded$AmyloGram_Pred == "ERROR", na.rm = TRUE)

summary <- data.frame(
  metric = c(
    "input_records",
    "valid_records",
    "invalid_or_short_records",
    "unique_sequences",
    "duplicates_saved",
    "chunk_size",
    "chunk_count",
    "cores",
    "amyloid_predictions",
    "non_amyloid_predictions",
    "error_predictions",
    "runtime_minutes"
  ),
  value = c(
    length(raw_seqs),
    length(valid_seqs),
    invalid_count,
    length(unique_seqs),
    length(valid_seqs) - length(unique_seqs),
    chunk_size,
    chunk_count,
    cores,
    sum(expanded$AmyloGram_Pred == "AMYLOID", na.rm = TRUE),
    sum(expanded$AmyloGram_Pred == "Non-Amyloid", na.rm = TRUE),
    error_count,
    sprintf("%.3f", elapsed_mins)
  )
)
write.table(summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(chunk_stats, file.path(chunk_dir, "chunk_stats.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

log_msg("DONE runtime_minutes=%.2f errors=%d", elapsed_mins, error_count)
log_msg("results=%s", output_file)
log_msg("summary=%s", summary_file)

if (error_count > 0 && !allow_errors) {
  stop(sprintf("AmyloGram completed with %d sequence-level errors", error_count), call. = FALSE)
}
