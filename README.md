# Plant Pathogen Genome Filtering

This project checks whether plant pathogens from the UC IPM disease list appear in an Evo2 training-organism file.

The program creates a text file with one pathogen per line:

```text
Agrobacterium tumefaciens	Y
Alternaria alternata	N
```

`Y` means the pathogen appears to match an organism in the Evo2 training file.
`N` means the pathogen was not found.

## Files

- `main.py` - the Python script
- `all_bacterial_species_in_evo2.txt` - the Evo2 training-organism file
- `evo2_plant_pathogen_matches.txt` - the output file created by the script

## How To Run

From this project folder, run:

```powershell
python main.py --evo2-file all_bacterial_species_in_evo2.txt --output evo2_plant_pathogen_matches.txt
```

The script will:

1. Read the Evo2 organism file.
2. Download the UC IPM plant disease list.
3. Pull out the pathogen scientific names.
4. Compare those names to the Evo2 names.
5. Write the final `Y` / `N` results file.

## Data Source

Plant pathogen names are downloaded from:

https://ipm.ucanr.edu/PMG/diseases/diseaseslist.html

