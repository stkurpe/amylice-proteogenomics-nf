args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Usage: Rscript benchmark.R <input.fasta>", call.=FALSE)

input_file <- args[1]

suppressPackageStartupMessages(library(AmyloGram))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(parallel))

cat("--- BENCHMARK MODE ---\n")
cat("Reading file...\n")
seqs_list <- read.fasta(file = input_file, seqtype = "AA", as.string = TRUE)
raw_seqs <- toupper(unlist(seqs_list))

# Очистка
clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", raw_seqs)
valid_mask <- nchar(clean_seqs) >= 6
final_seqs <- clean_seqs[valid_mask]
final_ids <- names(seqs_list)[valid_mask]

total_seqs <- length(final_seqs)
test_size <- 500 # Тестируем на 500 штук

if (total_seqs < test_size) test_size <- total_seqs

# Берем подвыборку
test_seqs <- final_seqs[1:test_size]
test_ids <- final_ids[1:test_size]

# Настройка ядер
num_cores <- detectCores()
if(num_cores > 1) num_cores <- num_cores - 1 

cat(sprintf("Running test on %d cores using %d sequences...\n", num_cores, test_size))

process_one_seq <- function(idx) {
    seq_str <- test_seqs[idx]
    tryCatch({
        seq_vector <- strsplit(seq_str, "")[[1]]
        input_data <- list(seq_vector)
        res <- predict(AmyloGram_model, input_data)
        return(TRUE)
    }, error = function(e) return(FALSE))
}

# --- ЗАМЕР ВРЕМЕНИ ---
start_time <- Sys.time()

# Запускаем расчет
dummy_res <- mclapply(1:length(test_seqs), process_one_seq, mc.cores = num_cores)

end_time <- Sys.time()
# ---------------------

time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
avg_time_per_seq <- time_taken / test_size
predicted_total_time <- avg_time_per_seq * total_seqs

cat("\n--- RESULTS ---\n")
cat(sprintf("Time for %d seqs: %.2f seconds\n", test_size, time_taken))
cat(sprintf("Average per seq:  %.4f seconds\n", avg_time_per_seq))
cat("------------------------------------------------\n")
cat(sprintf("PREDICTED TIME for %d seqs:\n", total_seqs))
cat(sprintf(">> %.1f minutes (or %.2f hours) <<\n", predicted_total_time / 60, predicted_total_time / 3600))
cat("------------------------------------------------\n")
