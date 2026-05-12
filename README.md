# Plant Pathogen Genome Filtering

This project checks plant pathogen names against `evo2_full_training_dataset.txt`, the full Evo2 training-organism file.

Evo2 is a genomic model trained on biological sequence data. Here, I am using the Evo2 organism list to see which plant pathogens show up in that training data.

The program writes two result files:

```text
# Y = found in the Evo2 organism file; N = not found.
# Multiple assembly IDs mean multiple Evo2 assemblies matched the same pathogen.
# Pathogen	Y/N	Assembly_IDs


Agrobacterium tumefaciens	Y	GCA_017744915.1;GCF_000834635.1;GCF_001541315.1;GCF_005221325.1;GCF_005221385.1;GCF_009649785.1;GCF_013337285.1;GCF_017726655.1
Alternaria alternata	N
```

`Y` means the pathogen appears to match an organism in the Evo2 organism file.
`N` means the pathogen was not found.
The third field contains the matching Evo2 assembly ID or IDs when a match is found. `N` rows leave that field blank.

For the eukaryotic output, the original `Input_Assembly_ID` and `Species_Name` columns are kept, then `Y/N` and `Evo2_Assembly_IDs` are added.

## Files

- `main.py` - the Python script
- `evo2_full_training_dataset.txt` - the Evo2 list used for comparison
- `additional_bacteria_plant_pathogens.txt` - extra bacterial plant pathogen names added to the UC IPM website names
- `eukaryotic_plant_pathogens.txt` - eukaryotic plant pathogen names
- `bacteria_plant_pathogen_results` - bacterial output file
- `eukaryotic_plant_pathogen_results` - eukaryotic output file

## How To Run

From this project folder, run:

```powershell
python main.py
```

This creates:

```text
bacteria_plant_pathogen_results
eukaryotic_plant_pathogen_results
```

To run only one part:

```powershell
python main.py --organism-type bacterial
python main.py --organism-type eukaryotic
```

The script will:

1. Read `evo2_full_training_dataset.txt`.
2. Make the bacteria list from the UC IPM website plus `additional_bacteria_plant_pathogens.txt`.
3. Read the eukaryote list from `eukaryotic_plant_pathogens.txt`.
4. Compare both lists to the Evo2 list.
5. Write the two result files with matching assembly IDs.

## Data Source

Plant pathogen names are downloaded from:

https://ipm.ucanr.edu/PMG/diseases/diseaseslist.html

