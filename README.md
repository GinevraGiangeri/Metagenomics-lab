# Metagenomics-lab
# **STRIM STRatification In Metagenomics** 
![logo2prova](https://repository-images.githubusercontent.com/366718180/d5a34c80-b7d8-11eb-8ca3-0c163cefd5d6)
## Using STRIM

### Cite Us
This code was developed as part of a project carried out during the course of Microbial Metagenomics (Molecular Biology Master Degree) at the University of Padova. The project was supervised by Prof. Stefano Campanaro and Dr. Arianna Basile.

## **Description**
The STRIM code is useful to calculate the abundances weighted for the occurrencies for each metabolic function (e.g. a function identified in the KEGG database) and it allows to define the stratification of the taxonomies for every different function. In other words, the stratification of the taxonomies means that the code enables the user to determine which is the contribution of each taxon to each specific function.
Mandatory input files:
* relative abundance of each microbial species in a community of a particular sample, file format must corrspond to that obtained using checkM (vn.n.n) 
* metabolic functions of each species of the community, annotated independently with eggNOG mapper (vn.n.n)
* the taxonomic assignment of each species according to NCBI classification
The script uses these information and it provides as output the stratification of the taxonomies contributing to each metabolic function.

## Requirements
In order to use the software, you will need to import:
* pandas: it is a library which allows you to perform data analysis and manipulation, particularly it is useful for manipulating numerical tables and time series.
* numpy: it is a library useful to manipulate large arrays and matrices and it offers high-level mathematical functions to work on these arrays.
* matplotlib.pyplot: it is a plotting library which enables to create static, animated or interactive plots and figures.
* os module: it enables the interaction between the user and the operating system.
* glob module: it is used to recover files or pathnames matching a specified pattern.
* fnmatch: supplies function to match files having a specific pattern or filter files having a specific pattern.
* collections: gives alternatives to create container data types such as list, tuple and dictionary.
* itertools: supplies functions that operate on iterators to produce complex iterators. An iterator is an object that contains a countable number of values that you can pass through, meaning that the vaues can be iterated upon. 

## Generation of the inputs and inputs required
The inputs needed will be:
* the abundances calculated with checkM
* the eggNOG annotations
* the taxonomic assignments with the NCBI classification

## How to use 
In order to use the code, the user will need to follow the subsequent steps:
* download the file strim.py
* open the file strim.py with a text editor and modify the paths in which the files you need are saved (every time you will have to open a file or save an output that is going to be generated e.g.'examplefilename.to_csv', you will have to specify the path in which the file needs to be saved)
* check if you have all the files required 
* change the columns names of the taxonomic assignments with the ones you have in your file containing the taxonomic assignments (e.g. column_list = ['average_Re1','average_Re2','average_Re3', 'average_Re4'] put the correct names between '')
* open the terminal and enter in the directory where the file strim.py is saved
* run 'python3 strim.py'
* the program will ask you a series of parameters as input:
  - the choice between KEGG orthologs and KEGG modules
  - the taxonomic level desired for the stratification step
* the program will elaborate the data and will return a message on the screen every time it will end the abundances calculation of a defined sample
* when the abundances calculation step is ended, the stratification step begins
* at the end of the stratification step, the program will return a message on the screen, inviting the user to check the results obtained.

## Output
In the folder in which you decided to save the outputs, the following output files should be found:
* two tabular files:
  - one defining the taxonomic assignment
  - one defining the abundances weighted fot the occurences for each metabolic function
* a series of pictures, each corresponding to the weighted abundances for each sample
* the stratification gives you back tabular files for each KEGG code analyzed!
