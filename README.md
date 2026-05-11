# Plant Pathogen Genome Filtering

This project checks whether plant pathogens from the UC IPM disease list and `additional_pathogens.txt` appear in an Evo2 training-organism file.

The program creates a text file with one pathogen per line:

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

## Files

- `main.py` - the Python script
- `additional_pathogens.txt` - extra pathogen names to include with the UC IPM website names
- `all_bacterial_species_in_evo2.txt` - the Evo2 assembly ID lookup file
- `evo2_plant_pathogen_matches.txt` - the output file created by the script

## How To Run

From this project folder, run:

```powershell
python main.py
```

The script will:

1. Download the UC IPM plant disease list.
2. Pull out the pathogen scientific names.
3. Read extra pathogen names from `additional_pathogens.txt`.
4. Read the Evo2 assembly ID lookup file.
5. Compare all pathogen names to the Evo2 names.
6. Write the final `Y` / `N` results file with matching assembly IDs.

## Data Source

Plant pathogen names are downloaded from:

https://ipm.ucanr.edu/PMG/diseases/diseaseslist.html

