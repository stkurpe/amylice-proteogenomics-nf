# AmyloGram Rust Port Analysis

## Current AmyloGram Core

The R `predict.ag_model` implementation does the following:

1. Validate amino acid alphabet.
2. Convert each input protein into overlapping 6-mer windows.
3. Encode windows through the AmyloGram degenerate alphabet.
4. Count selected mono/bi/tri-gram features.
5. Select 262 important features.
6. Run a `ranger` probability random forest with 500 trees.
7. Return the maximum amyloid probability across windows for each protein.

The current optimized runner already improves throughput by:

- deduplicating identical proteins;
- chunking unique sequences;
- dynamically parallelizing chunks;
- checkpointing every chunk.

## Why Rust Could Help

Rust can speed up:

- FASTA parsing;
- amino-acid validation;
- 6-mer window generation;
- feature encoding/counting;
- random forest inference;
- CSV writing;
- memory layout and parallel scheduling.

Expected gain after a faithful port is roughly **3-10x** for large proteomes, with much lower RAM use than R.

## Main Risk

The model is an R `ranger` object embedded in `AmyloGram_model`. A Rust port must reproduce its inference exactly enough for scientific use.

The blocker is not Rust itself; it is exporting and validating the random forest:

- tree split variables;
- split thresholds or categorical rules;
- terminal class probabilities;
- tree aggregation behavior;
- feature order and selected feature names.

## Recommended Migration Plan

### Phase 1: Extract Golden Test Corpus

Generate a small benchmark set:

- 100 short peptides;
- 100 medium proteins;
- 100 long proteins;
- known SRR32060234 sequences;
- edge cases with invalid/short sequences.

Save R outputs as golden fixtures.

### Phase 2: Export Model

Create an R exporter that serializes:

- `AmyloGram_model$enc`;
- `AmyloGram_model$imp_features`;
- `AmyloGram_model$rf$forest`;
- tree metadata required for inference.

Target format:

```text
amylogram_model.json
```

or a compact binary format after JSON validation.

### Phase 3: Rust Feature Engine

Implement:

- streaming FASTA reader;
- protein cleaning;
- 6-mer generation;
- degenerate alphabet mapping;
- selected feature counts.

Validate against R feature matrices before porting inference.

### Phase 4: Rust Random Forest Inference

Implement `ranger` probability forest traversal.

Validation target:

```text
abs(r_probability - rust_probability) <= 1e-9
```

or a documented tolerance if floating-point differences appear.

### Phase 5: CLI and Docker

Proposed CLI:

```bash
amylogram-rs predict \
  --model amylogram_model.json \
  --input combine_proteome.fasta \
  --output amylogram_prediction.csv \
  --threads 8
```

### Phase 6: Replace R Runner

Keep the R implementation as reference and fallback until the Rust implementation passes all contract tests.

## Recommendation

Do not rewrite the full package blindly. Build a Rust-compatible model exporter first, then port feature extraction and inference behind a contract-test suite. This gives speed without losing scientific traceability.

