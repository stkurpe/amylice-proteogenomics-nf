#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: Rscript predict_amylogram_parallel.R <input.fasta> <output.csv>", call.=FALSE)

input_file <- args[1]
output_file <- args[2]

# --- ПОДКЛЮЧЕНИЕ БИБЛИОТЕК ---
# Используем suppress... чтобы не засорять логи
suppressPackageStartupMessages(library(AmyloGram))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(parallel)) # Стандартная библиотека (вместо pbmcapply)

cat("Reading:", input_file, "\n")

tryCatch({
    seqs_list <- read.fasta(file = input_file, seqtype = "AA", as.string = TRUE)
}, error = function(e) stop("Fatal Error reading FASTA."))

ids <- names(seqs_list)
raw_seqs <- toupper(unlist(seqs_list))

# Очистка и фильтрация
clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", raw_seqs)
valid_mask <- nchar(clean_seqs) >= 6
final_seqs <- clean_seqs[valid_mask]
final_ids <- ids[valid_mask]

# Определяем ядра
num_cores <- detectCores()
if(num_cores > 1) num_cores <- num_cores - 1 

cat(sprintf("Starting PARALLEL analysis (Silent Mode) on %d cores for %d sequences...\n", num_cores, length(final_seqs)))

process_one_seq <- function(idx) {
    seq_str <- final_seqs[idx]
    seq_id <- final_ids[idx]
    
    tryCatch({
        seq_vector <- strsplit(seq_str, "")[[1]]
        input_data <- list(seq_vector)
        res <- predict(AmyloGram_model, input_data)
        
        return(list(
            Sequence_ID = seq_id,
            AmyloGram_Prob = res[["Probability"]],
            AmyloGram_Pred = ifelse(res[["Probability"]] > 0.5, "AMYLOID", "Non-Amyloid")
        ))
    }, error = function(e) {
        return(list(Sequence_ID = seq_id, AmyloGram_Prob = NA, AmyloGram_Pred = "ERROR"))
    })
}

# --- ЗАПУСК ЧЕРЕЗ mclapply (Быстро, но молча) ---
# mclapply входит в пакет parallel, он работает так же быстро, но без полоски загрузки
results_list <- mclapply(1:length(final_seqs), process_one_seq, mc.cores = num_cores)

cat("Aggregating results...\n")
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

if (any(results_df$AmyloGram_Pred == "ERROR")) {
    results_df <- results_df[results_df$AmyloGram_Pred != "ERROR", ]
}

write.csv(results_df, output_file, row.names = FALSE, quote = FALSE)
cat("✅ DONE! Results saved to:", output_file, "\n")
