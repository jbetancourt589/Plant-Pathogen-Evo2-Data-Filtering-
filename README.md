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
- `Scripts/Data Reproductions/compare_evo2_outputs_#1.py` - compares John and OpenGenome/Evo2 FASTA filtering outputs for reproduction set 1
- `Scripts/Data Reproductions/compare_evo2_outputs_gtf.py` - compares OpenGenome/Evo2 FASTA output for the GTF/GFF reproduction dataset
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

This project also includes reproducibility checks for Evo2/OpenGenome-style FASTA filtering.

The comparison scripts recreate an OpenGenome/Evo2-like FASTA from the original NCBI FASTA, then compare that recreated file to the official OpenGenome/Evo2 output. The recreated FASTA is strict: only `A`, `C`, `G`, and `T` are kept. `N` and other ambiguity characters such as `R`, `Y`, `M`, `K`, `W`, `S`, `B`, `D`, `H`, and `V` are treated as break points and are not copied into the recreated output.

Run the first comparison script from the project root:

```powershell
uv run python "Scripts/Data Reproductions/compare_evo2_outputs_#1.py"
```

The scripts accept `--original`, `--john`, `--opengenome`, `--annotation`, and `--outdir`. The `--john` file is optional for datasets that do not have a John-filtered FASTA.

### Optional GTF/GFF Annotation Precheck

Before normal strict OpenGenome-style filtering, the comparison scripts now check whether a GTF/GFF annotation file is available.

If no GTF/GFF file is provided or the path is missing, annotation filtering is skipped and the script continues with the old behavior:

1. Read the original NCBI FASTA.
2. Split each contig at every non-ACGT character.
3. Keep only continuous A/C/G/T chunks at least 10,000 bp long.
4. Name each kept chunk with 0-based half-open coordinates, such as `contig_id:2470-24374`.
5. Compare the recreated FASTA to the official OpenGenome/Evo2 FASTA.

If a GTF/GFF file is found, the script parses it first and removes only true centromere-region annotations before the normal strict filtering steps. A region is removed if the feature type is exactly `centromere` or if the attributes clearly mark the region itself as a centromere region. The script is conservative: it does not remove genes, transcripts, exons, CDS records, or proteins just because their names mention `centromere protein`, `centrosomal protein`, or `kinetochore protein`.

When annotation filtering runs, the script also writes:

- `annotation_filter_report.csv` - every parsed annotation row, whether it was selected for removal, the reason, and the number of base pairs removed

The script writes:

- `Results/Evo2 Data Reproduction/summary.txt` - human-readable summary of record counts, base counts, N/non-ACGT content, split-record validation, and reproduction results
- `Results/Evo2 Data Reproduction/comparison_report.csv` - record-level comparisons for John vs OpenGenome/Evo2 and recreated vs OpenGenome/Evo2
- `Results/Evo2 Data Reproduction/opengenome_chunk_validation.csv` - validation of OpenGenome-style `contig:start-end` records against the original NCBI contigs
- `Results/Evo2 Data Reproduction/recreated_opengenome_like.fasta` - recreated FASTA made by splitting original contigs into continuous A/C/G/T-only chunks of at least 10,000 bp

The analysis asks whether OpenGenome/Evo2 differs from John's output because it splits original contigs around `N` or other non-ACGT bases, keeps only long A/C/G/T-only chunks, and names those chunks with OpenGenome-style 0-based half-open coordinates.

### GTF/GFF Dataset Script

The GTF/GFF reproduction dataset can be run directly with:

```powershell
uv run python "Scripts/Data Reproductions/compare_evo2_outputs_gtf.py"
```

By default, this script uses:

- `Datasets/Evo2 Data Reproduction with GTF File/original_yeast_from_NCBI_GCF_000313485.1.fasta`
- `Datasets/Evo2 Data Reproduction with GTF File/open_genome2_filtered_GCF_000313485.1.fasta`
- `Datasets/Evo2 Data Reproduction with GTF File/genomic.gtf`

It writes results to:

```text
Results/Evo2 Data Reproduction with GTF/
```

