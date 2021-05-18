# Metagenomics-lab
# **STRIM STRatification In Metagenomics** 
![logo2prova](https://repository-images.githubusercontent.com/366718180/d5a34c80-b7d8-11eb-8ca3-0c163cefd5d6)
## Using STRIM

### Cite Us
This code was developed as part of a project carried out during the Microbial Metagenomics course (Molecular Biology Master Degree) at the University of Padova. The project was supervised by Prof. Stefano Campanaro and Dr. Arianna Basile.

## **Description**
STRIM performes two analyses: (1) it calculates the abundance of each KEGG ortholog or module taking into account the relative abundance of species having that function and the number of genes that are part of the module (e.g. all the genes encoding a function according to the KEGG database, times the abundance of the species having that function). (2) It stratifies the taxonomies for each function; by doing this, STRIM calculates the contribution of each taxon to the KEGG orthologs/modules.
Mandatory input files:
* relative abundance of each microbial species in a community of a particular sample, file format must corrspond to that obtained using checkM (tested with version 1.0.12);
* metabolic functions of each species of the community, annotated independently with eggNOG mapper (tested with version 2.0.1-1);
* the taxonomic assignment of each species according to NCBI classification.
For more information on the files format see the example files provided.
The script uses these information and it provides as output the abundance of the KEGG modules and the stratification of the taxonomies contributing to each metabolic function.

## Requirements
Software requirements:
* pandas: it is a library which allows you to perform data analysis and manipulation, particularly it is useful for manipulating numerical tables and time series.
* numpy: it is a library useful to manipulate large arrays and matrices and it offers high-level mathematical functions to work on these arrays.
* matplotlib.pyplot: it is a plotting library which enables to create static, animated or interactive plots and figures.
* os module: it enables the interaction between the user and the operating system.
* glob module: it is used to recover files or pathnames matching a specified pattern.
* fnmatch: supplies function to match files having a specific pattern or filter files having a specific pattern.
* collections: gives alternatives to create container data types such as lists, tuples and dictionaries.
* itertools: supplies functions to generate complex iterators. An iterator is an object that contains a countable number of values that you can pass through, meaning that the vaues can be iterated upon. 

## Generation of the inputs and inputs required
Mandatory input files:
* abundances of all species calculated with checkM;
* the eggNOG annotations results;
* the taxonomic assignments according to NCBI taxonomy.

## How to use 
In order to use the code, the user has to follow the subsequent steps:
* download the file "strim.py"
* as default the output files are saved in "", to modify this option, open the file "strim.py" with an editor and insert the path in where the output files should be saved (at line number...)
* change the columns names of the taxonomic assignments with the ones you have in your file containing the taxonomic assignments (e.g. column_list = ['average_Re1','average_Re2','average_Re3', 'average_Re4'] put the correct names between '')
* from the directory where the file strim.py is saved run 'python3 strim.py'
* the program will ask you a series of parameters as input:
  - selection of the KEGG features you would like to analyze (KEGG orthologs or KEGG modules)
  - the taxonomic level for the stratification step (species, genus, etc.)
* the script will start by calculating KEGG orthologs/modules abundances for each sample;
* when the abundances calculation step is ended, the stratification step begins.
The computation time requested for a typical workflow with 100-200 genomes and 10-20 different conditions will be within minutes depending on the computer used.

## Output
In the output folder, the following files will be saved:
* two tabular files:
  - one defining the taxonomic assignment
  - one defining the abundances weighted fot the occurences for each metabolic function;
* some images, each corresponding to the weighted abundances for each sample;
* the stratification gives you back tabular files for each KEGG code analyzed.
