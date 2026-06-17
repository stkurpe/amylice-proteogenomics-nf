# Pipeline Block Scheme

This diagram describes the pipeline as a set of biological analysis blocks, in the spirit of a no-code workflow such as `blockr`: inputs enter as data objects, biologically meaningful transformations happen in blocks, and only interpretable outputs continue to downstream interpretation.

## Biological Flow

```mermaid
flowchart TD
    A["RNA-seq sample"] --> B{"Available input type"}

    B -->|"Raw public data"| C1["Read quality control<br/>and preprocessing"]
    C1 --> C2["Expression quantification<br/>which transcripts are active?"]
    C1 --> C3["Variant discovery<br/>which loci are altered?"]

    B -->|"Prepared inputs"| D0["Expression table<br/>+ filtered variant calls<br/>+ reference annotation"]
    C2 --> D0
    C3 --> D0

    D0 --> E1["Expressed coding transcripts<br/>filter by TPM and CDS length"]
    E1 --> E2["Sample transcript whitelist<br/>biologically active CDS set"]

    E2 --> F1["SNP and nonsense branch"]
    D0 --> F1
    F1 --> F2["Haplotype-specific genomes<br/>maternal/paternal sequence context"]
    F2 --> F3["Translated haplotype proteins"]
    F3 --> F4["SNP-derived and<br/>nonsense-truncated proteins"]

    E2 --> G1["Frameshift branch"]
    D0 --> G1
    G1 --> G2["Coding indels mapped<br/>to expressed transcripts"]
    G2 --> G3["Novel frameshift proteins<br/>new C-terminal tails"]

    F4 --> H1["Sample-specific<br/>mutant proteome"]
    G3 --> H1
    H1 --> H2["Proteome quality gate<br/>complete and biologically usable FASTA"]

    H2 --> I1["Amyloidogenicity prediction<br/>aggregation-risk scoring"]
    H2 --> I2["Protein feature profiling<br/>length, charge, hydrophobicity,<br/>secondary-structure tendencies"]

    I1 --> J1["Amyloid risk table<br/>per protein sequence"]
    I2 --> J2["Interpretable protein<br/>feature table"]
    H2 --> J3["Mutant protein sequence set<br/>for downstream reuse"]

    J1 --> K["Integrated sample report"]
    J2 --> K
    J3 --> K

    B -->|"Proteome-only input"| H2

    classDef input fill:#e8f3ff,stroke:#9bb9d9,color:#1c2530
    classDef bio fill:#eaf7ee,stroke:#9cc5a8,color:#1c2530
    classDef variant fill:#fff2d8,stroke:#d8b56b,color:#1c2530
    classDef risk fill:#fdecef,stroke:#d7a6b0,color:#1c2530
    classDef output fill:#eef0ff,stroke:#aeb7df,color:#1c2530
    class A,B,C1,C2,C3,D0 input
    class E1,E2,H1,H2 bio
    class F1,F2,F3,F4,G1,G2,G3 variant
    class I1,I2 risk
    class J1,J2,J3,K output
```

## What Each Block Means

| Block | Biological purpose | Main input | Main output |
|---|---|---|---|
| Read quality control | Confirms that RNA-seq reads are suitable for expression and variant analysis. | FASTQ reads | Cleaned reads and QC evidence |
| Expression quantification | Identifies transcripts that are active in the sample. | RNA-seq reads or transcript abundance | Expression table |
| Variant discovery | Finds sequence changes that may alter proteins. | RNA-seq alignment and reference genome | Filtered variant calls |
| Expressed coding transcripts | Keeps biologically relevant CDS records and removes weak or unsafe candidates. | Expression table and GTF | Active CDS set |
| SNP and nonsense branch | Builds proteins affected by SNPs and early stop codons. | SNP variants and active CDS set | SNP-derived and nonsense-truncated proteins |
| Frameshift branch | Captures coding indels that can create novel protein tails. | Indel variants and active CDS set | Frameshift-derived proteins |
| Mutant proteome assembly | Combines all valid altered proteins into one sample-level proteome. | SNP/nonsense proteins and frameshift proteins | Mutant protein sequence set |
| Proteome quality gate | Ensures the proteome is complete enough for interpretation. | Mutant proteome outputs | Verified proteome |
| Amyloidogenicity prediction | Scores proteins for aggregation or amyloid-like risk. | Verified mutant proteome | Amyloid risk table |
| Protein feature profiling | Adds interpretable biochemical descriptors. | Verified mutant proteome | Protein feature table |
| Integrated sample report | Combines sequence, risk, and feature outputs for downstream review. | Proteome, predictions, features | Final sample-level outputs |

## Biological Conditions

| Condition | Meaning |
|---|---|
| Expression threshold | Only expressed transcripts are used to build the sample-specific proteome. |
| Coding sequence length threshold | Very short CDS records are excluded because they are unlikely to produce interpretable protein sequences. |
| Mitochondrial transcript exclusion | Mitochondrial coding rules differ from the standard nuclear genetic code, so these records are removed from this proteome model. |
| SNP/nonsense separation from frameshifts | SNPs and indels affect proteins differently, so they are modeled in separate branches before being recombined. |
| Frameshift protein length threshold | Very short frameshift products are not treated as meaningful downstream protein candidates. |
| Proteome quality gate | Amyloid prediction and feature profiling run only after the mutant proteome is complete and internally consistent. |

## Suggested `blockr` Representation

| blockr element | Suggested contents |
|---|---|
| Input blocks | RNA-seq reads, expression table, variant calls, reference annotation, optional prebuilt mutant proteome |
| Transformation blocks | read QC, expression quantification, variant discovery, expressed CDS filtering, haplotype translation, frameshift translation, proteome assembly |
| Analysis blocks | amyloidogenicity prediction, protein feature profiling |
| Decision blocks | input type, expression threshold, CDS length threshold, frameshift length threshold, proteome quality gate |
| Output blocks | mutant protein sequence set, amyloid risk table, protein feature table, integrated sample report |

## Visual Assessment

The simplified diagram renders more clearly than the first technical version because it removes low-level failure nodes, file existence checks, implementation names, and repeated output filenames. The remaining blocks are larger biological concepts, so the visual flow now reads as:

```text
RNA-seq -> expression + variants -> active coding transcripts
       -> SNP/nonsense proteins + frameshift proteins
       -> mutant protein sequence set
       -> amyloid risk + protein features
       -> integrated report
```

Visual check result: block sizes are balanced on a desktop-width preview, text wraps cleanly, and the main biological story is easy to scan. The only remaining limitation is that Mermaid may draw some cross-branch arrows tightly depending on renderer width; for slides or a blockr mockup, a four-lane layout with Inputs, Biological filtering, Protein generation, and Interpretation is the clearest version.
