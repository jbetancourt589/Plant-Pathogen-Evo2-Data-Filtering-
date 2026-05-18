# Plant Pathogen Evo2 Data Processing

This project processes plant-pathogen genome-related datasets and compares them with Evo2/OpenGenome data.

Evo2 is a genomic model trained on biological sequence data. Here, I am checking which names from the plant-pathogen dataset show up in `evo2_eukaryotic_dataset.txt`, the eukaryotic organism list from the Evo2 dataset, and in `evo2_full_training_dataset.txt`, the full Evo2 training-organism list.

The program writes two result files:

```text
# Y = found in evo2_eukaryotic_dataset.txt; N = not found.
# Species_Name	Y/N	Assembly_IDs


Agrobacterium tumefaciens	Y	GCA_017744915.1;GCF_000834635.1
Alternaria alternata	N
```

`Y` means the plant-pathogen name appears to match an organism in the source file named in the output header.
`N` means the plant-pathogen name was not found in that source file.
Each row contains the plant-pathogen species name, the comparison-specific `Y/N`, and the matching Evo2 assembly ID or IDs on the right. `N` rows leave the assembly-ID column blank.

## What It Does

The script builds one plant-pathogen dataset from two sources: the UC IPM plant disease page and `combined_plant_pathogen_list.txt`. It then normalizes names so small formatting differences do not block a match, including capitalization, punctuation, parenthetical notes, genus-level names such as `Pseudomonas spp.`, and strain/pathovar text after the genus and species.

It runs two comparisons:

1. The UC IPM names plus `combined_plant_pathogen_list.txt` are checked against `evo2_eukaryotic_dataset.txt`.
2. The same UC IPM names plus `combined_plant_pathogen_list.txt` are checked against `evo2_full_training_dataset.txt`.

Each output row contains a plant-pathogen species name from the website plus `combined_plant_pathogen_list.txt`, `Y` or `N`, and the matching Evo2 assembly ID or IDs if found.

## Files

- `Scripts/compare_plant_pathogens_to_evo2.py` - compares plant-pathogen names against Evo2 organism datasets
- `Scripts/compare_evo2_outputs.py` - compares John and OpenGenome/Evo2 FASTA filtering outputs
- `Datasets/Plant Pathogen Preprocessing Datasets/combined_plant_pathogen_list.txt` - local plant pathogen names added to the UC IPM website names
- `Datasets/Plant Pathogen Preprocessing Datasets/evo2_eukaryotic_dataset.txt` - eukaryotic organism dataset from Evo2
- `Datasets/Plant Pathogen Preprocessing Datasets/evo2_full_training_dataset.txt` - full Evo2 training dataset used for comparison
- `Results/Plant Pathogen Preprocessing Results/plant_pathogens_vs._eukaryotes_evo2` - UC IPM plus local pathogen names checked against `evo2_eukaryotic_dataset.txt`
- `Results/Plant Pathogen Preprocessing Results/plant_pathogen_vs._entire_evo2` - UC IPM plus local pathogen names checked against `evo2_full_training_dataset.txt`

## How To Run

From this project folder, run:

```powershell
uv run python Scripts/compare_plant_pathogens_to_evo2.py
```

If Python is already installed and on your `PATH`, `python Scripts/compare_plant_pathogens_to_evo2.py` works too.

This creates:

```text
Results/Plant Pathogen Preprocessing Results/plant_pathogens_vs._eukaryotes_evo2
Results/Plant Pathogen Preprocessing Results/plant_pathogen_vs._entire_evo2
```

The script will:

1. Download the UC IPM plant disease list.
2. Add the names from `combined_plant_pathogen_list.txt`.
3. Read `evo2_eukaryotic_dataset.txt` as the Evo2 eukaryotic comparison dataset.
4. Compare the website plus combined-list names to that Evo2 eukaryotic dataset.
5. Read `evo2_full_training_dataset.txt` as the full Evo2 comparison dataset.
6. Compare the website plus combined-list names to the full Evo2 dataset.
7. Write both result files with matching Evo2 assembly IDs, plant-pathogen species names, and `Y/N`.

## Data Source

Plant pathogen names are downloaded from:

https://ipm.ucanr.edu/PMG/diseases/diseaseslist.html

## Evo2/OpenGenome Filtering Reproducibility Test

This project also includes a reproducibility check for one organism/assembly where three FASTA files are available in `Datasets/`:

- `Datasets/orginial_from_ncbi_GCA_000359685.2.fasta` - the original NCBI input FASTA
- `Datasets/john_filtered_GCA_000359685.2.fasta` - John's filtered output FASTA
- `Datasets/open_genome2_filtered_GCA_000359685.2.fasta` - the OpenGenome/Evo2 filtered output FASTA

Run the comparison script from the project root:

```powershell
uv run python Scripts/compare_evo2_outputs.py
```

The script uses the three `Datasets/` FASTA files above by default and writes results to `Results/Evo2 Data Reproduction/`. You can still pass `--original`, `--john`, `--opengenome`, or `--outdir` if you want to compare different files.

The script writes:

- `Results/Evo2 Data Reproduction/summary.txt` - human-readable summary of record counts, base counts, N/non-ACGT content, split-record validation, and reproduction results
- `Results/Evo2 Data Reproduction/comparison_report.csv` - record-level comparisons for John vs OpenGenome/Evo2 and recreated vs OpenGenome/Evo2
- `Results/Evo2 Data Reproduction/opengenome_chunk_validation.csv` - validation of OpenGenome-style `contig:start-end` records against the original NCBI contigs
- `Results/Evo2 Data Reproduction/recreated_opengenome_like.fasta` - recreated FASTA made by splitting original contigs into continuous A/C/G/T-only chunks of at least 10,000 bp

The analysis asks whether OpenGenome/Evo2 differs from John's output because it splits original contigs around `N` or other non-ACGT bases, keeps only long A/C/G/T-only chunks, and names those chunks with OpenGenome-style 0-based half-open coordinates.

