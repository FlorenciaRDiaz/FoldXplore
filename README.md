# AlphaFoldXplore

[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AngieLCerrutti/AlphaFoldXplore/blob/main/example/AlphaFoldXplore_example.ipynb)


## Usage

Install alphafoldxplore.py on a folder without anything else.

## API
```python
import alphafoldxplore
```
### Functions

```python
set_up()
#Downloads and installs AlphaFold.

predict(dir_string)
#Does not work without running set_up() beforehand.
#Predicts the terciary structure of a single protein (or a group of them) by reading a FASTA file.
#Simplified to default parameters to speed up the process. Creates two folders 'json_files' and 'pdb_files' with the results inside.

extract_zips(dir_string)
#Unneeded to use normally. 
#Creates two folders 'json_files' and 'pdb_files' and stores the .json and .pdb files from the folders inside.

get_pae_files(dir_string)
#Reads all .json files from a folder (by default, 'json_files'). Assumes the .json files are those from the predictions.
#Returns a dictionary.

get_plddt_files(dir_string)
#Reads all .pdb files from a folder (by default, 'pdb_files') and extracts the pLDDT values from its CA atoms.
#Returns a dictionary.

pae_results(protein_1,protein_2)
#Compares the PAE values of two proteins by reading .json files and creating heatmaps. Protein_2 is optional.
#Admits both strings and entries from the get_pae_files() dictionary.
#Prints heatmaps.

plddt_results(protein_1,protein_2)
#Compares the pLDDT values of two proteins by reading .pdb files and plotting values of all CA atoms. Protein_2 is optional.
#Admits both strings and entries from the get_plddt_files() dictionary.
#Prints a plot.

watch_proteins(protein_1, protein_2)
#Doesn't work on Colab.
#Creates a widget and allows the visualization of protein structures by reading .pdb files. Protein_2 is optional.
#Admits strings.
#Returns a NGLView view.

superimpose_proteins(protein_1,protein_2)
#Rotates and translates the molecular structure of protein_2 to superimpose (match) it with protein_1.
#Both proteins must be of the same length for it to work properly.
#Creates a .pdb file on the root folder named "superimposed_{filename}.pdb".
#Admits strings.
#Prints the mean RMSD.

calc_individual_rmsd(protein_1,protein_2,start,end)
#Calculates the individual RMSD (Root mean-square deviation) between CA atoms of both proteins.
#protein_1 and protein_2 must be strings. start and end must be positive int numbers and are optional.
#Plots the result and prints the mean RMSD.
#Returns a list with the values per CA pair.
```

