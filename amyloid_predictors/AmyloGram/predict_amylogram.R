#!/usr/bin/env Rscript

# --- НАСТРОЙКИ ---
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: Rscript predict_amylogram.R <input.fasta> <output.csv>", call.=FALSE)

input_file <- args[1]
output_file <- args[2]

suppressPackageStartupMessages(library(AmyloGram))
suppressPackageStartupMessages(library(seqinr))

cat("Reading:", input_file, "\n")

# 1. Читаем файл
tryCatch({
    seqs_list <- read.fasta(file = input_file, seqtype = "AA", as.string = TRUE)
}, error = function(e) stop("Fatal Error reading FASTA."))

ids <- names(seqs_list)
raw_seqs <- toupper(unlist(seqs_list))

# 2. Очистка
clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", raw_seqs)

# 3. Фильтр длины
valid_mask <- nchar(clean_seqs) >= 6
final_seqs <- clean_seqs[valid_mask]
final_ids <- ids[valid_mask]

cat(sprintf("Starting analysis of %d sequences...\n", length(final_seqs)))

# Векторы для хранения результатов
probs <- numeric(length(final_seqs))
preds <- character(length(final_seqs))
has_error <- logical(length(final_seqs))

# Прогресс-бар
pb_step <- max(1, floor(length(final_seqs) / 20))

# --- ЦИКЛ ОБРАБОТКИ ---
for (i in 1:length(final_seqs)) {
    
    if (i %% pb_step == 0) cat(sprintf("  Progress: %d%%\n", round(i/length(final_seqs)*100)))
    
    tryCatch({
        # ВАЖНОЕ ИЗМЕНЕНИЕ: Разбиваем строку на вектор букв (как в тесте)
        # "MKA..." -> c("M", "K", "A"...)
        seq_vector <- strsplit(final_seqs[i], "")[[1]]
        
        # Оборачиваем в список, так как AmyloGram любит списки
        input_data <- list(seq_vector)
        
        # Предсказываем
        res <- predict(AmyloGram_model, input_data)
        
        probs[i] <- res[["Probability"]]
        preds[i] <- ifelse(res[["Probability"]] > 0.5, "AMYLOID", "Non-Amyloid")
        has_error[i] <- FALSE
        
    }, error = function(e) {
        # Если упало - просто пропускаем
        probs[i] <- NA
        preds[i] <- "ERROR"
        has_error[i] <- TRUE
    })
}

# Сборка результатов
results_df <- data.frame(
    Sequence_ID = final_ids,
    AmyloGram_Prob = probs,
    AmyloGram_Pred = preds
)

# Удаляем ошибочные строки
if (any(has_error)) {
    cat(sprintf("Warning: %d sequences caused errors and were removed.\n", sum(has_error)))
    results_df <- results_df[!has_error, ]
}

write.csv(results_df, output_file, row.names = FALSE, quote = FALSE)
cat("✅ DONE! Results saved to:", output_file, "\n")
