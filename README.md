# Plant Pathogen Genome Filtering

This project checks a plant-pathogen dataset against Evo2 organism datasets.

Evo2 is a genomic model trained on biological sequence data. Here, I am checking which names from the plant-pathogen dataset show up in `evo2_eukaryotic_dataset.txt`, the eukaryotic organism list from the Evo2 dataset, and in `evo2_full_training_dataset.txt`, the full Evo2 training-organism list.

The program writes two result files:

```text
# Y = found in evo2_eukaryotic_dataset.txt; N = not found.
# Pathogen	Y/N


Agrobacterium tumefaciens	Y
Alternaria alternata	N
```

`Y` means the pathogen appears to match an organism in the source file named in the output header.
`N` means the pathogen was not found.

## What It Does

The script builds one plant-pathogen dataset from two sources: the UC IPM plant disease page and `combined_plant_pathogen_list.txt`. It then normalizes names so small formatting differences do not block a match, including capitalization, punctuation, parenthetical notes, genus-level names such as `Pseudomonas spp.`, and strain/pathovar text after the genus and species.

It runs two comparisons:

1. The UC IPM names plus `combined_plant_pathogen_list.txt` are checked against `evo2_eukaryotic_dataset.txt`.
2. The same UC IPM names plus `combined_plant_pathogen_list.txt` are checked against `evo2_full_training_dataset.txt`.

Each output row contains the pathogen name and `Y` or `N`. The output files do not include assembly IDs; they only report whether the pathogen matched the selected Evo2 source file.

## Files

- `main.py` - the Python script
- `combined_plant_pathogen_list.txt` - local plant pathogen names added to the UC IPM website names
- `evo2_eukaryotic_dataset.txt` - eukaryotic organism dataset from Evo2
- `evo2_full_training_dataset.txt` - full Evo2 training dataset used for comparison
- `plant_pathogens_vs._eukaryotes_evo2` - UC IPM plus local pathogen names checked against `evo2_eukaryotic_dataset.txt`
- `plant_pathogen_vs._entire_evo2` - UC IPM plus local pathogen names checked against `evo2_full_training_dataset.txt`

## How To Run

From this project folder, run:

```powershell
uv run python main.py
```

If Python is already installed and on your `PATH`, `python main.py` works too.

This creates:

```text
plant_pathogens_vs._eukaryotes_evo2
plant_pathogen_vs._entire_evo2
```

The script will:

1. Download the UC IPM plant disease list.
2. Add the names from `combined_plant_pathogen_list.txt`.
3. Read `evo2_eukaryotic_dataset.txt` as the Evo2 eukaryotic comparison dataset.
4. Compare the website plus combined-list names to that Evo2 eukaryotic dataset.
5. Read `evo2_full_training_dataset.txt` as the full Evo2 comparison dataset.
6. Compare the website plus combined-list names to the full Evo2 dataset.
7. Write both result files with `Y/N`.

## Data Source

Plant pathogen names are downloaded from:

https://ipm.ucanr.edu/PMG/diseases/diseaseslist.html

