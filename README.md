### Metagenomics-lab
# **STRIM STRatification In Metagenomics** 
![logo2prova](https://repository-images.githubusercontent.com/366718180/d5a34c80-b7d8-11eb-8ca3-0c163cefd5d6)
## Using STRIM

### Cite Us
This code was developed as part of a project carried out during the Microbial Metagenomics course (Molecular Biology Master Degree) at the University of Padova. The project was supervised by Prof. Stefano Campanaro and Dr. Arianna Basile.

## Description
STRIM performs two analyses: (1) it calculates the abundance of each KEGG ortholog or module taking into account the relative abundance of species having that function and the number of genes that are part of the module (e.g. all the genes encoding a function according to the KEGG database, times the abundance of the species having that function). (2) It stratifies the taxonomies for each function; by doing this, STRIM calculates the contribution of each taxon to the KEGG orthologs/modules.

## Input files required
Mandatory input files:
* the relative abundance of each microbial species in a community of a particular sample, file format must correspond to that obtained using checkM (tested with version 1.0.12);
* the metabolic functions of each species of the community, annotated independently with eggNOG mapper (tested with version 2.0.1-1);
* the taxonomic assignment of each species according to the NCBI classification.
For more information on the files format see the example files provided.
The script uses these information and it provides as output the abundance of the KEGG modules and the stratification of the taxonomies contributing to each metabolic function.

Note: all the input files must be placed in the same directory.

## Requirements
Software requirements:
* pandas (version 1.2.4)
* numpy (version 1.20.3 )
* matplotlib.pyplot (version 3.4.2 )
* os module (part of the Standard Library of Python 3)
* glob module (part of the Standard Library of Python 3)
* fnmatch module (part of the Standard Library of Python 3)
* collections module (part of the Standard Library of Python 3)
* itertools module (part of the Standard Library of Python 3)
* argparse module (part of the Standard Library of Python 3)

If you do not have these libraries installed, please follow these procedures:
- pandas: https://pypi.org/project/pandas/
- numpy: https://pypi.org/project/numpy/
- matplotlib.pyplot: https://pypi.org/project/matplotlib/

## How to use 
In order to use the code, the user has to follow the subsequent steps:
* download the file "strim.py"
* from the directory where the file strim.py is saved run 'python3 strim.py <input directory> <output directory>'
* the program will ask you a series of parameters as input:
  - selection of the KEGG features you would like to analyze (KEGG orthologs or KEGG modules)
  - the taxonomic level for the stratification step (species, genus, etc.)
* the script will start by calculating KEGG orthologs/modules abundances for each sample;
* when the abundances calculation step is ended, the stratification step begins.
The computational time requested for a typical workflow with 100-200 genomes and 10-20 different conditions will be 5-20 minutes depending on the computer used.

## Output
In the output folder, the following files will be saved:
* two tabular files:
  - one defining the taxonomic assignment
  - one defining the abundances weighted fot the occurences for each metabolic function;
* image files, each corresponding to the weighted abundances for each sample;
* the stratification analysis provides the tabular files for each KEGG code analyzed.
