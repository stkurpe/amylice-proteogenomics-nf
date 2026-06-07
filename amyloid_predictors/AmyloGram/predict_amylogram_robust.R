#!/usr/bin/env Rscript

# --- НАСТРОЙКИ ---
BATCH_MODE <- FALSE # Если TRUE - быстро, но может упасть. FALSE - надежно, по одному.

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
# Преобразуем в вектор и сразу в верхний регистр
raw_seqs <- toupper(unlist(seqs_list))

# 2. Жесткая чистка регулярными выражениями
# Удаляем любые символы, кроме 20 стандартных аминокислот
clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", raw_seqs)

# 3. Фильтр длины (минимум 6 аминокислот)
valid_mask <- nchar(clean_seqs) >= 6
final_seqs <- clean_seqs[valid_mask]
final_ids <- ids[valid_mask]

cat(sprintf("Processing %d sequences (Skipped %d too short/empty)...\n", 
            length(final_seqs), length(raw_seqs) - length(final_seqs)))

# 4. ПОШТУЧНАЯ ОБРАБОТКА (Fault Tolerant Loop)
# Мы создадим пустые векторы для результатов
probs <- numeric(length(final_seqs))
preds <- character(length(final_seqs))
has_error <- logical(length(final_seqs))

# Прогресс-бар (текстовый)
pb_step <- max(1, floor(length(final_seqs) / 10))

cat("Starting robust prediction loop. Please wait...\n")

for (i in 1:length(final_seqs)) {
    # Отображение прогресса
    if (i %% pb_step == 0) cat(sprintf("  ... processed %d / %d\n", i, length(final_seqs)))
    
    # Попытка предсказать ОДНУ последовательность
    tryCatch({
        # AmyloGram требует вектор, подаем вектор длины 1
        # Важно: убираем имена (unname), иногда пакеты R сходят с ума от метаданных
        single_seq <- unname(final_seqs[i])
        
        # Предсказание
        res <- predict(AmyloGram_model, single_seq)
        
        probs[i] <- res[["Probability"]]
        preds[i] <- ifelse(res[["Probability"]] > 0.5, "AMYLOID", "Non-Amyloid")
        has_error[i] <- FALSE
        
    }, error = function(e) {
        # Если упало - логируем и идем дальше
        probs[i] <- NA
        preds[i] <- "ERROR_SKIPPED"
        has_error[i] <- TRUE
        
        # Вывод информации о "битом" белке (для диагностики)
        cat(sprintf("\n[WARNING] Crashed on Seq Index %d (ID: %s)\n", i, final_ids[i]))
        cat(sprintf("Seq Content: %s\nContinuing...\n", final_seqs[i]))
    })
}

# 5. Собираем результаты
results_df <- data.frame(
    Sequence_ID = final_ids,
    AmyloGram_Prob = probs,
    AmyloGram_Pred = preds
)

# Удаляем строки, где были ошибки (если они есть)
if (any(has_error)) {
    cat(sprintf("\nRemoved %d sequences that caused internal errors.\n", sum(has_error)))
    results_df <- results_df[!has_error, ]
}

# 6. Сохраняем
write.csv(results_df, output_file, row.names = FALSE, quote = FALSE)
cat("Success! Results saved to:", output_file, "\n")
