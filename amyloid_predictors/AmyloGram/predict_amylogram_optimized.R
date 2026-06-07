args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: Rscript predict_amylogram_optimized.R <input.fasta> <output.csv>", call.=FALSE)

input_file <- args[1]
output_file <- args[2]

# --- БИБЛИОТЕКИ ---
suppressPackageStartupMessages(library(AmyloGram))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(parallel))

cat("--- OPTIMIZED PARALLEL MODE ---\n")
cat("Reading:", input_file, "\n")

# Чтение
tryCatch({
    seqs_list <- read.fasta(file = input_file, seqtype = "AA", as.string = TRUE)
}, error = function(e) stop("Fatal Error reading FASTA."))

# Подготовка данных
ids <- names(seqs_list)
raw_seqs <- toupper(unlist(seqs_list))

# Очистка и фильтрация
clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", raw_seqs)
valid_mask <- nchar(clean_seqs) >= 6

final_seqs <- clean_seqs[valid_mask]
final_ids <- ids[valid_mask]
total_seqs <- length(final_seqs)

if (total_seqs == 0) stop("No valid sequences found!")

# Определяем ядра
num_cores <- detectCores()
if(num_cores > 1) num_cores <- num_cores - 1 

cat(sprintf("Processing %d sequences on %d cores using CHUNKING strategy...\n", total_seqs, num_cores))

# --- ГЛАВНАЯ ОПТИМИЗАЦИЯ: Разбиваем данные на куски (Chunks) ---
# Создаем индексы для разбиения на группы (по количеству ядер)
chunk_indices <- cut(seq_along(final_seqs), breaks = num_cores, labels = FALSE)
# Разбиваем сами последовательности и ID на списки
chunks_seqs <- split(final_seqs, chunk_indices)
chunks_ids <- split(final_ids, chunk_indices)

# Функция обработки ЦЕЛОГО КУСКА за один раз
process_chunk <- function(idx) {
    chunk_s <- chunks_seqs[[idx]]
    chunk_i <- chunks_ids[[idx]]
    
    # 1. Подготовка всего куска сразу (векторизация strsplit)
    # AmyloGram требует список, где каждый элемент - вектор букв
    input_list <- strsplit(chunk_s, "")
    
    # 2. ОДИН вызов predict для тысяч последовательностей
    # Это убирает overhead вызова функции
    tryCatch({
        res <- predict(AmyloGram_model, input_list)
        
        # 3. Сборка результата
        # AmyloGram возвращает объект, из которого можно достать векторы
        return(data.frame(
            Sequence_ID = chunk_i,
            AmyloGram_Prob = res[["Probability"]],
            AmyloGram_Pred = ifelse(res[["Probability"]] > 0.5, "AMYLOID", "Non-Amyloid"),
            stringsAsFactors = FALSE
        ))
    }, error = function(e) {
        # Если упало на всем чанке (редко), возвращаем ошибку для чанка
        return(data.frame(
            Sequence_ID = chunk_i,
            AmyloGram_Prob = NA,
            AmyloGram_Pred = "CHUNK_ERROR",
            stringsAsFactors = FALSE
        ))
    })
}

# --- ЗАПУСК ---
start_time <- Sys.time()

# Теперь mclapply запускается всего num_cores раз (например, 7 раз),
# но каждый процесс делает огромную работу
results_list <- mclapply(seq_along(chunks_seqs), process_chunk, mc.cores = num_cores)

# Сборка
cat("Aggregating results...\n")
results_df <- do.call(rbind, results_list)

# Удаление ошибок (если были)
if (any(results_df$AmyloGram_Pred == "CHUNK_ERROR")) {
    cat("Warning: Some chunks failed processing.\n")
    results_df <- results_df[results_df$AmyloGram_Pred != "CHUNK_ERROR", ]
}

end_time <- Sys.time()
run_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat(sprintf("✅ DONE in %.2f minutes!\n", run_time))
write.csv(results_df, output_file, row.names = FALSE, quote = TRUE)
cat("Results saved to:", output_file, "\n")
