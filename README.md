# Code Execution Instructions

## Overview
I developed two versions of the code:
1. **Normal Version**: Loads the entire genome into memory at once.
   - File: `CS249_Assignment0_Normal_AbdelRahman_Alsabbagh.py`
2. **Optimized Version**: Reads the genome line by line to optimize memory usage.
   - File: `CS249_Assignment0_Optimized_AbdelRahman_Alsabbagh.py`
   - Once a chromosome is fully processed, pattern matching is performed, and the chromosome is deleted from memory.

Additionally, to ensure correctness, I created a set of small unit tests on sample texts instead of running the code directly on full genomes. These unit tests are included in:
- `CS249_Assignment0_Unit_AbdelRahman_Alsabbagh.py`

## Installation
To install all dependencies, run:
```bash
pip install -r requirements.txt
```

## Running the Code
### Exact Match on GRCh38 (hg38) Genome
```bash
python CS249_Assignment0_Optimized_AbdelRahman_Alsabbagh.py --genome <hg38_Fasta> \
--pattern <Pattern_Fasta> --output results_exact_match_hg38.txt
```

### Approximate Match on GRCh38 (hg38) Genome
```bash
python CS249_Assignment0_Optimized_AbdelRahman_Alsabbagh.py --genome <hg38_Fasta> \
--pattern <Pattern_Fasta> --output results_approximate_match_hg38.txt --approx
```

### Exact Match on T2T CHM13v2.0 Genome
```bash
python CS249_Assignment0_Optimized_AbdelRahman_Alsabbagh.py --genome <t2t_Fasta> \
--pattern <Pattern_Fasta> --output results_exact_match_t2t.txt
```

### Approximate Match on T2T CHM13v2.0 Genome
```bash
python CS249_Assignment0_Optimized_AbdelRahman_Alsabbagh.py --genome <t2t_Fasta> \
--pattern <Pattern_Fasta> --output results_approximate_match_t2t.txt --approx
```

### Running Unit Tests
```bash
python CS249_Assignment0_Unit_AbdelRahman_Alsabbagh.py
```

## Disclaimer
I utilized Generative AI models to enhance the aesthetics of my code, particularly for writing documentation and discussing function arguments and return values, as well as improving this report. However, all other aspects of the work were completed independently.

